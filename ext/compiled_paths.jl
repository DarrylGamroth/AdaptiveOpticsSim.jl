# Derived from SPIDERSSim.jl b479871e017e8f1df2e56880bba6f47f80d1273a
# src/optics/physical15.jl. Preserves the lowered spectre Fresnel and spectral
# arithmetic; instrument-specific selection and export remain cold data.
using FixedSizeArrays: FixedSizeVectorDefault, FixedSizeMatrixDefault
import KernelAbstractions
using KernelAbstractions: @index, @localmem, @synchronize

# The compiler's centred (y,x) arrays are the numerical contract. Only PROPER
# propagation crosses to its corner-origin storage, through prepared scratch.
@enum PathOperationKind begin
    PathMultiply
    PathNormalize
    PathTransferFunction
    PathFourierTransform
    PathSpectralMultiplier
    PathResample
    PathLiveScreen
end

struct PathOperation{T,A}
    kind::PathOperationKind
    array::A
    auxiliary::A
    post::A
    dx::T
    dz::T
    wavelength_m::T
    gouy::Complex{T}
    from_waist::Bool
    forward::Bool
    source::UInt8 # 1 = DM OPD; 2 = entrance OPD and log-amplitude attenuation
end

struct OpticalPathPlan{A,O}
    path_index::Int
    anchor_index::Int
    entrance::A
    prefix::O
    closed::O
    only::O
end

struct PathTransformPlan{A,V}
    source_index::Int
    shift_transfer::A
    phase_multiplier::V
    defocus_transfer::V
    ghost::V
end

struct PathPrefixPlan
    cache_index::Int
    first_path::Int
    count::Int
end

struct CompiledPathsPlan{T,D,P,F,A}
    target::D
    grid::Int
    image_size::Int
    crop_start::Int
    n_env::Int
    n_anchor::Int
    wavelengths::FixedSizeVectorDefault{T}
    paths::P
    shared_prefixes::FixedSizeVectorDefault{PathPrefixPlan}
    fast_paths::F
    resample_kernel::A
    main_phase::A
    pinhole_gain::FixedSizeVectorDefault{Complex{T}}
    pixel_mtf::A
    phantom_pad::Int
    pixel_pad::Int
    apply_pixel_mtf::Bool
end

struct PathOperationWorkspace{A,R,E}
    input::A
    temporary::A
    mapped::A
    opd::R
    electric_field::E
end

struct PathFFTWorkspace{A,F,I}
    field::A
    scratch::A
    forward::F
    inverse::I
end

struct PathNormalizationWorkspace{P,S}
    partials::P
    power::S
end

struct CompiledPathsWorkspace{A,W,C,O,F,V,N,M,L}
    field::A
    prefix::A
    shared_prefixes::V
    wavefront::W
    context::C
    operations::O
    physical_fft::F
    image_fft::F
    phantom_fft::F
    crop::A
    resample_temporary::A
    closed_envelopes::V
    only_envelopes::V
    open_envelopes::V
    closed_sources::V
    only_sources::V
    normalization::N
    pixel_mtf::M
    mul_coefficients::L
end

"""Canonical `(y,x)` compiled multi-path relative-intensity products, dimensionless intensity."""
struct CompiledPathsProducts{R}
    fringed::R
    unfringed::R
    reference::R
end

"""Prepared compiled optical path execution with exact input buffer bindings and products."""
struct PreparedCompiledPaths{P,W,I,O,C,S}
    plan::P
    workspace::W
    inputs::I
    products::O
    execution_context::C
    synthesis::S
end

# Prepared normalization launches use 256 lanes, cld(length(field), 256)
# partials and one power element. The field guard and fixed reduction tree
# bound every access, including padded lanes and local_i + stride.
# Two fixed reduction stages retain every result on the device. Each workgroup
# reduces in one defined tree order; no atomics or device-to-host scalar queries.
KernelAbstractions.@kernel function _compiled_paths_power_partials_kernel!(partials, field, count)
    i = @index(Global, Linear)
    local_i = @index(Local, Linear)
    group = @index(Group, Linear)
    shared = @localmem(eltype(partials), (256,))
    @inbounds shared[local_i] = i <= count ? abs2(field[i]) : zero(eltype(partials))
    @synchronize
    for stride in (128, 64, 32, 16, 8, 4, 2, 1)
        if local_i <= stride
            @inbounds shared[local_i] += shared[local_i + stride]
        end
        @synchronize
    end
    if local_i == 1
        @inbounds partials[group] = shared[1]
    end
end

KernelAbstractions.@kernel function _compiled_paths_power_finalize_kernel!(power, partials, count)
    local_i = @index(Local, Linear)
    shared = @localmem(eltype(partials), (256,))
    value = zero(eltype(partials))
    # count = cld(length(field), 256), so the final i + 256 is at most
    # count + 256 and remains representable as an Int.
    i = local_i
    while i <= count
        @inbounds value += partials[i]
        i += 256
    end
    @inbounds shared[local_i] = value
    @synchronize
    for stride in (128, 64, 32, 16, 8, 4, 2, 1)
        if local_i <= stride
            @inbounds shared[local_i] += shared[local_i + stride]
        end
        @synchronize
    end
    if local_i == 1
        @inbounds power[1] = shared[1]
    end
end

KernelAbstractions.@kernel function _compiled_paths_normalize_kernel!(field, power, count)
    i = @index(Global, Linear)
    if i <= count
        @inbounds field[i] = iszero(power[1]) ? zero(eltype(field)) : field[i]/sqrt(power[1])
    end
end

function _compiled_paths_normalization_workspace(::Matrix, target, ::Type{T}) where {T}
    return PathNormalizationWorkspace(
        FixedSizeVectorDefault{T}(undef, 0), FixedSizeVectorDefault{T}(undef, 0))
end

function _compiled_paths_normalization_workspace(field, target, ::Type{T}) where {T}
    return PathNormalizationWorkspace(
        _compiled_paths_array(target, T, cld(length(field), 256)),
        _compiled_paths_array(target, T, 1))
end

function _compiled_paths_normalize!(field::Matrix, workspace)
    power=sum(abs2,field)
    if iszero(power)
        fill!(field,zero(eltype(field)))
    else
        field ./= sqrt(power)
    end
    return field
end

function _compiled_paths_normalize!(field, workspace)
    style = Backends.execution_style(field)
    Backends.launch_kernel_async!(style, _compiled_paths_power_partials_kernel!,
        workspace.partials, field, length(field); ndrange=256length(workspace.partials))
    Backends.launch_kernel_async!(style, _compiled_paths_power_finalize_kernel!,
        workspace.power, workspace.partials, length(workspace.partials); ndrange=256)
    Backends.launch_kernel_async!(style, _compiled_paths_normalize_kernel!,
        field, workspace.power, length(field); ndrange=length(field))
    return field
end

# KA masks lanes outside the output axes. Prepared shifts have equal positive
# input/output sizes and 0 <= shift <= axis length. Each shifted coordinate
# lies between 1 and twice the axis length; one subtraction gives a valid index.
KernelAbstractions.@kernel function _compiled_paths_shift_kernel!(out, input, sy, sx)
    y, x = @index(Global, NTuple)
    ny, nx = size(input)
    source_y = y + sy
    source_x = x + sx
    source_y > ny && (source_y -= ny)
    source_x > nx && (source_x -= nx)
    @inbounds out[y, x] = input[source_y, source_x]
end

# Cold preparation bounds crop_start by input_size - output_size; phantom
# crops use nonnegative padding with input_size = output_size + 2padding.
KernelAbstractions.@kernel function _compiled_paths_crop_kernel!(out, input, start)
    y, x = @index(Global, NTuple)
    @inbounds out[y, x] = input[y + start, x + start]
end

# The inside guard proves both source coordinates are within the input axes.
KernelAbstractions.@kernel function _compiled_paths_pad_kernel!(out, input, pad)
    y, x = @index(Global, NTuple)
    inside = pad < y <= size(input, 1) + pad && pad < x <= size(input, 2) + pad
    @inbounds out[y, x] = inside ? input[y - pad, x - pad] : zero(eltype(out))
end

function _compiled_paths_shift!(out, input, inverse::Bool=false)
    n, m = size(input)
    sy = inverse ? fld(n, 2) : cld(n, 2)
    sx = inverse ? fld(m, 2) : cld(m, 2)
    Backends.launch_kernel_async!(Backends.execution_style(out),
        _compiled_paths_shift_kernel!, out, input, sy, sx; ndrange=size(out))
    return out
end

function _compiled_paths_crop!(out, input, start)
    Backends.launch_kernel_async!(Backends.execution_style(out),
        _compiled_paths_crop_kernel!, out, input, start; ndrange=size(out))
    return out
end

function _compiled_paths_pad!(out, input, pad)
    Backends.launch_kernel_async!(Backends.execution_style(out),
        _compiled_paths_pad_kernel!, out, input, pad; ndrange=size(out))
    return out
end


function _compiled_paths_shift!(out::Matrix, input::Matrix, inverse::Bool=false)
    n, m = size(input)
    size(out) == (n, m) || throw(DimensionMismatch("compiled optical paths shift shapes differ"))
    sy = inverse ? fld(n, 2) : cld(n, 2)
    sx = inverse ? fld(m, 2) : cld(m, 2)
    # Both arrays are one-based dense matrices of the checked equal dimensions;
    # each single wrap below maps to a valid input index.
    @inbounds for x in 1:m
        source_x = x + sx
        source_x > m && (source_x -= m)
        for y in 1:n
            source_y = y + sy
            source_y > n && (source_y -= n)
            out[y, x] = input[source_y, source_x]
        end
    end
    return out
end

function _compiled_paths_crop!(out::Matrix, input::Matrix, start)
    for x in axes(out, 2), y in axes(out, 1)
        out[y, x] = input[y + start, x + start]
    end
    return out
end

function _compiled_paths_pad!(out::Matrix, input::Matrix, pad)
    fill!(out, zero(eltype(out)))
    for x in axes(input, 2), y in axes(input, 1)
        out[y + pad, x + pad] = input[y, x]
    end
    return out
end

function _compiled_paths_fft!(out, input, fft_workspace, inverse=false)
    _compiled_paths_shift!(fft_workspace.field, input, true)
    if inverse
        Backends.enqueue_fft_plan!(fft_workspace.field, fft_workspace.inverse)
    else
        Backends.enqueue_fft_plan!(fft_workspace.field, fft_workspace.forward)
    end
    _compiled_paths_shift!(out, fft_workspace.field)
    return out
end

function _compiled_paths_filter!(out, input, transfer, fft_workspace)
    _compiled_paths_fft!(fft_workspace.scratch, input, fft_workspace)
    fft_workspace.scratch .*= transfer
    _compiled_paths_fft!(out, fft_workspace.scratch, fft_workspace, true)
    return out
end

_compiled_paths_seal(values::Vector{T}) where {T} = FixedSizeVectorDefault{T}(values)

# Array identity proves the frozen coefficients came from the same cold model
# binding; scalar equality includes wavelength, screen source and propagation
# direction. Distinct bindings are conservatively executed independently.
function _compiled_paths_same_operation(a, b)
    return a.kind == b.kind && a.array === b.array &&
        a.auxiliary === b.auxiliary && a.post === b.post &&
        isequal(a.dx, b.dx) && isequal(a.dz, b.dz) &&
        isequal(a.wavelength_m, b.wavelength_m) && isequal(a.gouy, b.gouy) &&
        a.from_waist == b.from_waist && a.forward == b.forward && a.source == b.source
end

function _compiled_paths_shared_prefixes(paths)
    prefixes = fill(PathPrefixPlan(0, 0, 0), length(paths))
    visited = falses(length(paths))
    cache_count = 0
    for first_path in eachindex(paths)
        visited[first_path] && continue
        reference = paths[first_path]
        group = [i for i in first_path:length(paths) if
            paths[i].anchor_index == reference.anchor_index &&
            paths[i].entrance === reference.entrance]
        visited[group] .= true
        length(group) > 1 || continue
        count = 0
        for index in 1:minimum(i -> length(paths[i].prefix), group)
            all(i -> _compiled_paths_same_operation(reference.prefix[index],
                paths[i].prefix[index]), group) || break
            count = index
        end
        count > 0 || continue
        cache_count += 1
        prefixes[group] .= Ref(PathPrefixPlan(cache_count, first_path, count))
    end
    return _compiled_paths_seal(prefixes), cache_count
end

function _compiled_paths_array(target, ::Type{T}, dims...) where {T}
    out = Backends.allocate_device_array(target, T, dims...)
    fill!(out, zero(T))
    return out
end

function _compiled_paths_convert(values, ::Type{T}, role) where {T}
    all(isfinite, values) || throw(ArgumentError("$role contains nonfinite source values"))
    converted = T.(values)
    all(isfinite, converted) || throw(ArgumentError("$role overflows or becomes nonfinite after conversion to $T"))
    return converted
end

function _compiled_paths_resident(target, values, ::Type{T}) where {T}
    converted = _compiled_paths_convert(values, T, "compiled optical paths array")
    out = _compiled_paths_array(target, T, size(values)...)
    copyto!(out, converted)
    return out
end

function _compiled_paths_fft_workspace(target, ::Type{T}, n) where {T}
    field = _compiled_paths_array(target, Complex{T}, n, n)
    scratch = similar(field)
    fill!(scratch, zero(Complex{T}))
    forward = Backends.plan_fft_backend!(field)
    inverse = Backends.plan_ifft_backend!(field)
    return PathFFTWorkspace(field, scratch, forward, inverse)
end

function _compiled_paths_electric_field(field, dx, wavelength)
    optics = AdaptiveOpticsSim.Optics
    metadata = optics.OpticalPlaneMetadata(optics.IntermediatePlane(), field;
        coordinate_domain=optics.MetricCoordinates(), sampling=(dx, dx),
        origin=(-fld(size(field, 1), 2) * dx, -fld(size(field, 2), 2) * dx),
        orientation=optics.PlaneAxisOrientation((:y, :x)),
        spectral=optics.MonochromaticChannel(wavelength),
        normalization=optics.DimensionlessNormalization())
    return optics.ElectricField(metadata, field)
end

function _compiled_paths_operation(table, loadarray, empty, ::Type{T}) where {T}
    kind = table["kind"]
    fields = if kind=="mult" || kind=="spectral_mult"
        (("kind","array"),("kind","array"))
    elseif kind=="normalize"
        (("kind",),("kind",))
    elseif kind=="tf"
        (("kind","dx","dz","lam","edge"),("kind","dx","dz","lam"))
    elseif kind=="ft"
        (("kind","dx","dz","lam","edge","from_waist","forward","gouy_real","gouy_imag"),
            ("kind","dx","dz","lam","from_waist","forward","gouy_real","gouy_imag"))
    elseif kind=="resample"
        (("kind","kernel","roll","post"),("kind","kernel","roll","post"))
    elseif kind=="live_screen"
        (("kind","input","pitch","dx","lam","size","kernel"),("kind","input","dx","lam","kernel"))
    else
        throw(ArgumentError("unsupported compiled optical operation: $kind"))
    end
    all(k->k in fields[1],keys(table)) && all(k->haskey(table,k),fields[2]) ||
        throw(ArgumentError("unknown or missing compiled optical operation fields for $kind"))
    kinds = (mult=PathMultiply, normalize=PathNormalize,
        tf=PathTransferFunction, ft=PathFourierTransform, spectral_mult=PathSpectralMultiplier,
        resample=PathResample, live_screen=PathLiveScreen)
    key = Symbol(kind)
    hasproperty(kinds, key) || throw(ArgumentError("unsupported compiled_paths operation: $kind"))
    array = haskey(table, "array") ? loadarray(table["array"]) :
            haskey(table, "kernel") ? loadarray(table["kernel"]) : empty
    auxiliary = haskey(table, "edge") ? loadarray(table["edge"]) :
                haskey(table, "roll") ? loadarray(table["roll"]) : empty
    post = haskey(table, "post") ? loadarray(table["post"]) : empty
    source = UInt8(0)
    if kind == "live_screen"
        input = table["input"]
        input in ("dm_opd", "pupil_opd") ||
            throw(ArgumentError("unsupported compiled optical live input: $input"))
        source = input == "dm_opd" ? UInt8(1) : UInt8(2)
    end
    if kind == "ft"
        table["forward"] == (table["dz"] >= 0) ||
            throw(ArgumentError("compiled optical paths FT direction disagrees with signed distance"))
    end
    for key in ("dx","lam")
        value=T(get(table,key,1))
        isfinite(value) && value>0 || throw(ArgumentError("compiled optical $key must be finite and positive"))
    end
    isfinite(T(get(table,"dz",0))) || throw(ArgumentError("compiled optical distance must be finite"))
    if kind=="live_screen" && haskey(table,"size")
        table["size"]==size(array,2) || throw(DimensionMismatch("compiled live input size differs from its interpolation kernel"))
    end
    return PathOperation(getproperty(kinds, key), array, auxiliary, post,
        T(get(table, "dx", 1)), T(get(table, "dz", 0)), T(get(table, "lam", 1)),
        complex(T(get(table, "gouy_real", 1)), T(get(table, "gouy_imag", 0))),
        get(table, "from_waist", false), get(table, "forward", true), source)
end

function _compiled_paths_operation_workspace(op, field, target, ::Type{T}) where {T}
    n = size(field, 1)
    m = op.kind == PathLiveScreen ? size(op.array, 2) : 0
    temporary = _compiled_paths_array(target, Complex{T},
        op.kind in (PathLiveScreen, PathResample) ? n : 0,
        op.kind == PathLiveScreen ? m : op.kind == PathResample ? n : 0)
    input = _compiled_paths_array(target, Complex{T}, m, m)
    mapped = _compiled_paths_array(target, Complex{T}, m == 0 ? 0 : n, m == 0 ? 0 : n)
    opd = _compiled_paths_array(target, T, m == 0 ? 0 : n, m == 0 ? 0 : n)
    electric_field = _compiled_paths_electric_field(field, op.dx, op.wavelength_m)
    return PathOperationWorkspace(input, temporary, mapped, opd, electric_field)
end

function _compiled_paths_bind(target, ::Type{T}, value, n, role) where {T}
    value === nothing && return _compiled_paths_array(target, T, n, n)
    _validate_compiled_paths_buffer(target, value, T, (n, n), role)
    return value
end

"""
    prepare_compiled_optical_paths(metadata, read_array; T=Float64, target=HostComputeDevice(),
        dm_opd=nothing, entrance_opd=nothing, entrance_attenuation=nothing)

Prepare the compiler's frozen numerical contract. `read_array(reference)` must
return a canonical `(y,x)` host array. Bind metre-valued OPD buffers once.
`entrance_attenuation` is the fitted model's log-amplitude attenuation `a`, applied
as `exp(-a)`; its zero-filled default is unit transmission. Live inputs add to
frozen fitted screens. Python and JAX are absent from frame execution.
"""
function AdaptiveOpticsSim.Optics.prepare_compiled_optical_paths(metadata, read_array;
    target=Backends.HostComputeDevice(), kwargs...)
    execution_context = Backends._prepare_device_execution_context(target)
    return Backends._with_prepared_device_execution_context(execution_context) do
        owner = _prepare_compiled_paths(metadata, read_array, execution_context; target, kwargs...)
        Backends._synchronize_prepared_device_execution_context!(execution_context)
        owner
    end
end

function _prepare_compiled_paths(metadata, read_array, execution_context;
    T::Type{<:Union{Float32,Float64}}=Float64,
    target=Backends.HostComputeDevice(), dm_opd=nothing, entrance_opd=nothing,
    entrance_attenuation=nothing, summation::Symbol=:path_pairs)
    metadata["schema_version"] == 1 || throw(ArgumentError("unsupported compiled_paths schema"))
    metadata["array_orientation"] == "yx" || throw(ArgumentError("compiled optical paths requires yx arrays"))
    n = Int(metadata["grid"])
    image_size = Int(metadata["image_size"])
    crop_start = Int(metadata["crop_start"])
    n > 0 && image_size > 0 && 0 <= crop_start <= n - image_size ||
        throw(ArgumentError("invalid compiled_paths crop/grid"))
    empty = _compiled_paths_array(target, Complex{T}, 0, 0)
    A = typeof(empty)
    cache = Dict{String,A}()
    function loadarray(reference)
        get!(cache, String(reference)) do
            _compiled_paths_resident(target, read_array(reference), Complex{T})
        end
    end
    operation_type = typeof(_compiled_paths_operation(Dict("kind" => "normalize"), loadarray, empty, T))
    paths = OpticalPathPlan{A,FixedSizeVectorDefault{operation_type}}[]
    for p in metadata["propagations"]
        sequences = map(("prefix", "closed", "only")) do name
            ops = operation_type[_compiled_paths_operation(op, loadarray, empty, T) for op in p[name]]
            _compiled_paths_seal(ops)
        end
        push!(paths, OpticalPathPlan(Int(p["path_index"]) + 1, Int(p["anchor_index"]) + 1,
            loadarray(p["entrance"]), sequences...))
    end
    na = length(metadata["wavelengths"])
    nc = length(metadata["channels"])
    ne = Int(metadata["n_env"])
    length(paths) == ne * na || throw(ArgumentError("compiled optical paths path/anchor coverage incomplete"))
    slots = sort([(p.path_index, p.anchor_index) for p in paths])
    slots == sort([(p, a) for p in 1:ne for a in 1:na]) ||
        throw(ArgumentError("compiled optical paths path/anchor slots are duplicated or missing"))
    fast = metadata["fast"]
    phantom_pad = Int(fast["phantom_pad"])
    pixel_pad = Int(fast["pixel_pad"])
    max_pad = (typemax(Int) - image_size) ÷ 2
    0 <= phantom_pad <= max_pad && 0 <= pixel_pad <= max_pad ||
        throw(ArgumentError("compiled optical paths padding must be nonnegative and its padded size must fit Int"))
    fastpath_type = PathTransformPlan{A,FixedSizeVectorDefault{A}}
    fast_paths = fastpath_type[]
    for p in metadata["fast_paths"]
        phases = _compiled_paths_seal(A[loadarray(ref) for ref in p["phase_multiplier"]])
        defocus = _compiled_paths_seal(A[loadarray(ref) for ref in get(p, "defocus_transfer", String[])])
        ghosts = _compiled_paths_seal(A[loadarray(ref) for ref in get(p, "ghost_field_anchor", String[])])
        length(phases) == na || throw(ArgumentError("compiled optical paths phase anchor count differs"))
        push!(fast_paths, PathTransformPlan(Int(p["source_index"]) + 1,
            loadarray(p["shift_transfer"]), phases, defocus, ghosts))
    end
    nk = length(fast_paths)
    basis_source = read_array(fast["basis"])
    all(isfinite, basis_source) || throw(ArgumentError("compiled optical paths spectral basis contains nonfinite source values"))
    basis = FixedSizeMatrixDefault{T}(_compiled_paths_convert(real.(basis_source), T, "compiled optical paths spectral basis"))
    coefficients = _compiled_paths_convert(read_array(fast["coefficient"]), Complex{T}, "compiled optical paths spectral coefficient")
    size(coefficients) == (nc, nk) || throw(DimensionMismatch("compiled optical paths coefficient shape differs"))
    coefficient_rows = FixedSizeMatrixDefault{Complex{T}}(coefficients)
    visibility = _compiled_paths_seal([FixedSizeMatrixDefault{Complex{T}}(
        _compiled_paths_convert(read_array(ref), Complex{T}, "compiled optical paths spectral visibility"))
        for ref in fast["visibility"]])
    pg = _compiled_paths_convert(complex.(fast["pinhole_gain_real"], fast["pinhole_gain_imag"]),
        Complex{T}, "compiled optical paths pinhole gain")
    size(basis) == (nc, na) || throw(DimensionMismatch("compiled optical paths basis shape differs"))
    length(visibility) == nc || throw(DimensionMismatch("compiled optical paths visibility channel count differs"))
    all(v -> size(v) == (nk, nk), visibility) ||
        throw(DimensionMismatch("compiled optical paths visibility path count differs"))
    length(pg) == na || throw(DimensionMismatch("compiled optical paths pinhole gain anchor count differs"))
    for path in paths
        size(path.entrance) == (n, n) || throw(DimensionMismatch("compiled optical paths entrance grid differs"))
        for sequence in (path.prefix, path.closed, path.only), op in sequence
            if op.kind in (PathMultiply, PathSpectralMultiplier, PathResample)
                size(op.array) == (n, n) || throw(DimensionMismatch("compiled optical paths operation grid differs"))
            elseif op.kind == PathLiveScreen
                size(op.array, 1) == n || throw(DimensionMismatch("compiled optical paths live kernel output differs"))
            end
            isempty(op.auxiliary) || size(op.auxiliary) == (n, n) ||
                throw(DimensionMismatch("compiled optical paths edge/roll grid differs"))
            isempty(op.post) || size(op.post) == (n, n) ||
                throw(DimensionMismatch("compiled optical paths resampling post grid differs"))
        end
    end
    for path in fast_paths
        0 <= path.source_index <= ne || throw(ArgumentError("compiled optical paths fast source index differs"))
        size(path.shift_transfer) == (image_size, image_size) ||
            throw(DimensionMismatch("compiled optical paths shift grid differs"))
        all(a -> size(a) == (image_size, image_size), path.phase_multiplier) ||
            throw(DimensionMismatch("compiled optical paths fast phase grid differs"))
        if path.source_index == 0
            length(path.ghost) == na || throw(DimensionMismatch("compiled optical paths ghost anchor count differs"))
            all(a -> size(a) == (image_size, image_size), path.ghost) ||
                throw(DimensionMismatch("compiled optical paths ghost grid differs"))
        end
        if !isempty(path.defocus_transfer)
            length(path.defocus_transfer) == na ||
                throw(DimensionMismatch("compiled optical paths defocus anchor count differs"))
            padded = image_size + 2phantom_pad
            all(a -> size(a) == (padded, padded), path.defocus_transfer) ||
                throw(DimensionMismatch("compiled optical paths defocus grid differs"))
        end
    end
    size(loadarray(fast["resample_kernel"])) == (image_size, image_size) ||
        throw(DimensionMismatch("compiled optical paths detector resampling grid differs"))
    size(loadarray(fast["main_phase"])) == (image_size, image_size) ||
        throw(DimensionMismatch("compiled optical paths main phase grid differs"))
    pixel_padded = image_size + 2pixel_pad
    size(loadarray(fast["pixel_mtf"])) == (pixel_padded, pixel_padded) ||
        throw(DimensionMismatch("compiled optical paths pixel MTF grid differs"))
    shared_prefixes, prefix_count = _compiled_paths_shared_prefixes(paths)
    plan = CompiledPathsPlan(target, n, image_size, crop_start, ne, na,
        FixedSizeVectorDefault{T}(T.(metadata["wavelengths"])),
        _compiled_paths_seal(paths), shared_prefixes, _compiled_paths_seal(fast_paths),
        loadarray(fast["resample_kernel"]), loadarray(fast["main_phase"]),
        FixedSizeVectorDefault{Complex{T}}(pg),
        loadarray(fast["pixel_mtf"]), phantom_pad, pixel_pad,
        get(fast, "apply_pixel_mtf", true))
    field = _compiled_paths_array(target, Complex{T}, n, n)
    prefix = similar(field)
    fill!(prefix, zero(Complex{T}))
    internal = similar(field)
    fill!(internal, zero(Complex{T}))
    wf = Proper.WaveFront(internal, T(first(metadata["wavelengths"])), one(T), zero(T), one(T))
    ctx = Proper.RunContext(wf)
    operations = _compiled_paths_seal([
        _compiled_paths_operation_workspace(op, field, target, T)
        for p in paths for seq in (p.prefix, p.closed, p.only) for op in seq])
    dm_size = 0
    entrance_size = 0
    for p in paths, seq in (p.prefix, p.closed, p.only), op in seq
        if op.kind == PathLiveScreen
            if op.source == 1
                dm_size == 0 || dm_size == size(op.array, 2) ||
                    throw(DimensionMismatch("compiled optical paths DM source sizes differ across paths"))
                dm_size = size(op.array, 2)
            else
                entrance_size == 0 || entrance_size == size(op.array, 2) ||
                    throw(DimensionMismatch("compiled optical paths entrance source sizes differ across paths"))
                entrance_size = size(op.array, 2)
            end
        end
    end
    inputs = (dm_opd=_compiled_paths_bind(target, T, dm_opd, dm_size, "compiled optical paths DM OPD"),
        entrance_opd=_compiled_paths_bind(target, T, entrance_opd, entrance_size, "compiled optical paths entrance OPD"),
        entrance_attenuation=_compiled_paths_bind(target, T, entrance_attenuation, entrance_size, "compiled optical paths entrance attenuation"))
    function fields(count)
        _compiled_paths_seal(A[_compiled_paths_array(target, Complex{T}, image_size, image_size) for _ in 1:count])
    end
    image_fft = _compiled_paths_fft_workspace(target, T, image_size)
    physical_fft = _compiled_paths_fft_workspace(target, T, n)
    phantom_fft = _compiled_paths_fft_workspace(target, T, image_size + 2plan.phantom_pad)
    products = CompiledPathsProducts(
        _compiled_paths_array(target, T, image_size, image_size),
        _compiled_paths_array(target, T, image_size, image_size),
        _compiled_paths_array(target, T, image_size, image_size))
    pixel_mtf = _compiled_paths_prepare_pixel_mtf(plan, products)
    prefix_fields = _compiled_paths_seal(A[
        _compiled_paths_array(target, Complex{T}, n, n) for _ in 1:prefix_count])
    workspace = CompiledPathsWorkspace(field, prefix, prefix_fields, wf, ctx, operations,
        physical_fft, image_fft, phantom_fft,
        _compiled_paths_array(target, Complex{T}, image_size, image_size),
        _compiled_paths_array(target, Complex{T}, image_size, image_size),
        fields(nk * na), fields(nk * na), fields(nk * na), fields(ne * na), fields(ne * na),
        _compiled_paths_normalization_workspace(field, target, T), pixel_mtf,
        Backends._prepare_mul_coefficients(field))
    channels = AdaptiveOpticsSim.Optics.SpectralBundle(T.(metadata["channels"]),
        ones(T, nc); T)
    reference_visibility = [Complex{T}[1;;] for _ in 1:nc]
    reference_inputs = _compiled_paths_seal(A[workspace.only_envelopes[a] for a in 1:na])
    synthesis = (
        AdaptiveOpticsSim.Optics.prepare_spectral_synthesis(workspace.closed_envelopes,
            basis, coefficient_rows, visibility, channels; output=products.unfringed, summation),
        AdaptiveOpticsSim.Optics.prepare_spectral_synthesis(workspace.open_envelopes,
            basis, coefficient_rows, visibility, channels; output=products.fringed, summation),
        AdaptiveOpticsSim.Optics.prepare_spectral_synthesis(reference_inputs,
            basis, coefficient_rows[:, 1:1], reference_visibility, channels;
            output=products.reference, summation),
    )
    owner = PreparedCompiledPaths(plan, workspace, inputs, products, execution_context, synthesis)
    # Prime PROPER's shape-dependent storage and plans outside execution.
    for p in paths, seq in (p.prefix, p.closed, p.only), op in seq
        if op.kind in (PathTransferFunction, PathFourierTransform)
            _compiled_paths_proper!(workspace, op)
        end
    end
    fft_workspace = Proper.fft_workspace(ctx)
    fft_workspace.plans_valid || throw(ArgumentError("compiled optical paths PROPER FFT plans were not prepared"))
    wf.field === fft_workspace.scratch || throw(ArgumentError(
        "compiled optical paths PROPER field must retain its prepared FFT scratch binding"))
    return owner
end

function _compiled_paths_proper!(workspace, op)
    wf, ctx = workspace.wavefront, workspace.context
    if !isempty(op.auxiliary)
        workspace.field .*= op.auxiliary
    end
    _compiled_paths_shift!(wf.field, workspace.field, true)
    wf.wavelength_m = op.wavelength_m
    wf.sampling_m = op.dx
    wf.z_m = zero(op.dx)
    if op.kind == PathTransferFunction
        wf.reference_surface = Proper.PLANAR
        Proper.prop_ptp(wf, op.dz, ctx)
    elseif op.from_waist
        wf.reference_surface = Proper.PLANAR
        Proper.prop_wts(wf, op.dz, ctx)
    else
        wf.reference_surface = Proper.SPHERICAL
        Proper.prop_stw(wf, op.dz, ctx)
    end
    _compiled_paths_shift!(workspace.field, wf.field)
    if op.kind == PathFourierTransform
        workspace.field .*= op.gouy
    end
    return nothing
end

function _compiled_paths_interpolate!(work, kernel, source, coefficients)
    work.input .= source
    Backends._enqueue_mul!(work.temporary, kernel, work.input, coefficients)
    Backends._enqueue_mul!(work.mapped, work.temporary, transpose(kernel), coefficients)
    work.opd .= real.(work.mapped)
    return work.opd
end

function _compiled_paths_execute!(owner, sequence, offset, first=1, last=length(sequence))
    w = owner.workspace
    for index in first:last
        op = sequence[index]
        work = w.operations[offset + index]
        if op.kind == PathMultiply
            w.field .*= op.array
        elseif op.kind == PathNormalize
            _compiled_paths_normalize!(w.field, w.normalization)
        elseif op.kind in (PathTransferFunction, PathFourierTransform)
            _compiled_paths_proper!(w, op)
        elseif op.kind == PathSpectralMultiplier
            _compiled_paths_filter!(w.field, w.field, op.array, w.physical_fft)
        elseif op.kind == PathResample
            _compiled_paths_fft!(w.physical_fft.scratch, w.field, w.physical_fft)
            if !isempty(op.auxiliary)
                w.physical_fft.scratch .*= op.auxiliary
            end
            Backends._enqueue_mul!(work.temporary, op.array, w.physical_fft.scratch,
                w.mul_coefficients)
            Backends._enqueue_mul!(w.field, work.temporary, transpose(op.array),
                w.mul_coefficients)
            w.field ./= size(op.array, 2)^2
            if !isempty(op.post)
                w.field .*= op.post
            end
        elseif op.kind == PathLiveScreen
            input = op.source == 1 ? owner.inputs.dm_opd : owner.inputs.entrance_opd
            _compiled_paths_interpolate!(work, op.array, input, w.mul_coefficients)
            AdaptiveOpticsSim.Optics.apply_phase_async!(work.electric_field, work.opd)
            if op.source == 2
                _compiled_paths_interpolate!(work, op.array, owner.inputs.entrance_attenuation,
                    w.mul_coefficients)
                @. w.field *= exp(-work.opd)
            end
        end
    end
    return offset + length(sequence)
end

@inline _compiled_paths_slot(path, anchor, n_anchor) = (path - 1) * n_anchor + anchor

function _compiled_paths_slow!(owner, share_prefix::Bool=true)
    plan, w = owner.plan, owner.workspace
    offset = 0
    for (index, path) in enumerate(plan.paths)
        shared = plan.shared_prefixes[index]
        if share_prefix && shared.count > 0
            cached = w.shared_prefixes[shared.cache_index]
            if index == shared.first_path
                Backends.copyto_backend_async!(w.field, path.entrance)
                _compiled_paths_execute!(owner, path.prefix, offset, 1, shared.count)
                Backends.copyto_backend_async!(cached, w.field)
            else
                Backends.copyto_backend_async!(w.field, cached)
            end
            # Each operation owns scratch at its original offset, including
            # skipped live-screen workspaces. Prefix fields are refreshed at
            # the first member on every invocation, also during GPU capture.
            # PROPER operations rebind their wavelength, sampling, position
            # and reference surface in _compiled_paths_proper!; no trajectory
            # geometry or FFT scratch binding is restored from this field.
            offset = _compiled_paths_execute!(owner, path.prefix, offset, shared.count + 1)
        else
            Backends.copyto_backend_async!(w.field, path.entrance)
            offset = _compiled_paths_execute!(owner, path.prefix, offset)
        end
        Backends.copyto_backend_async!(w.prefix, w.field)
        offset = _compiled_paths_execute!(owner, path.closed, offset)
        slot = _compiled_paths_slot(path.path_index, path.anchor_index, plan.n_anchor)
        _compiled_paths_crop!(w.crop, w.field, plan.crop_start)
        _compiled_paths_resample_image!(w.closed_sources[slot], w.crop, owner)
        Backends.copyto_backend_async!(w.field, w.prefix)
        offset = _compiled_paths_execute!(owner, path.only, offset)
        _compiled_paths_crop!(w.crop, w.field, plan.crop_start)
        _compiled_paths_resample_image!(w.only_sources[slot], w.crop, owner)
    end
    for a in 1:plan.n_anchor
        w.closed_sources[a] .*= plan.main_phase
        w.only_sources[a] .*= plan.main_phase
    end
    return nothing
end

function _compiled_paths_resample_image!(out, input, owner)
    plan, w = owner.plan, owner.workspace
    _compiled_paths_fft!(w.image_fft.scratch, input, w.image_fft)
    Backends._enqueue_mul!(w.resample_temporary, plan.resample_kernel, w.image_fft.scratch,
        w.mul_coefficients)
    Backends._enqueue_mul!(out, w.resample_temporary, transpose(plan.resample_kernel),
        w.mul_coefficients)
    out ./= plan.image_size^2
    return out
end

function _compiled_paths_fast!(envelopes, sources, owner, free_ghost)
    plan, w = owner.plan, owner.workspace
    for (k, path) in enumerate(plan.fast_paths), a in 1:plan.n_anchor
        out = envelopes[_compiled_paths_slot(k, a, plan.n_anchor)]
        if path.source_index == 0
            if free_ghost
                Backends.copyto_backend_async!(out, path.ghost[a])
            else
                fill!(out, zero(eltype(out)))
            end
        else
            source = sources[_compiled_paths_slot(path.source_index, a, plan.n_anchor)]
            if isempty(path.defocus_transfer)
                Backends.copyto_backend_async!(out, source)
            else
                _compiled_paths_pad!(w.phantom_fft.scratch, source, plan.phantom_pad)
                _compiled_paths_filter!(w.phantom_fft.scratch, w.phantom_fft.scratch,
                    path.defocus_transfer[a], w.phantom_fft)
                _compiled_paths_crop!(out, w.phantom_fft.scratch, plan.phantom_pad)
            end
        end
        _compiled_paths_filter!(out, out, path.shift_transfer, w.image_fft)
        out .*= path.phase_multiplier[a]
    end
    return nothing
end

function _compiled_paths_prepare_pixel_mtf(plan, products)
    plan.apply_pixel_mtf || return nothing
    centered = Array(plan.pixel_mtf)
    all(iszero, imag.(centered)) || throw(ArgumentError(
        "compiled optical paths pixel MTF transfer must have real values"))
    # The artifact stores fftshift(fftfreq) ordering; AOS accepts corner-origin
    # FFT ordering. Conversion belongs to cold calibration preparation.
    corner = circshift(real.(centered),
        (-fld(size(centered, 1), 2), -fld(size(centered, 2), 2)))
    return map((products.fringed, products.unfringed, products.reference)) do input
        output = _compiled_paths_array(plan.target, eltype(input), size(input)...)
        AdaptiveOpticsSim.Detectors.prepare_pixel_mtf(input, corner;
            padding=plan.pixel_pad, output)
    end
end

_compiled_paths_apply_pixel_mtf!(::Nothing) = nothing

function _compiled_paths_apply_pixel_mtf!(owners::Tuple)
    for owner in owners
        output = AdaptiveOpticsSim.Detectors.apply_pixel_mtf!(owner)
        Backends.copyto_backend_async!(owner.input, output)
    end
    return nothing
end

function _enqueue_compiled_paths!(owner::PreparedCompiledPaths, share_prefix::Bool=true)
    plan, w, products = owner.plan, owner.workspace, owner.products
    _compiled_paths_slow!(owner, share_prefix)
    _compiled_paths_fast!(w.closed_envelopes, w.closed_sources, owner, true)
    _compiled_paths_fast!(w.only_envelopes, w.only_sources, owner, false)
    for k in 1:length(plan.fast_paths), a in 1:plan.n_anchor
        slot = _compiled_paths_slot(k, a, plan.n_anchor)
        gain = plan.pinhole_gain[a]
        open = w.open_envelopes[slot]
        closed = w.closed_envelopes[slot]
        only = w.only_envelopes[slot]
        @. open = closed + gain * only
    end
    for synthesis in owner.synthesis
        AdaptiveOpticsSim.Optics.synthesize_spectrum!(synthesis)
    end
    _compiled_paths_apply_pixel_mtf!(w.pixel_mtf)
    return products
end

"""Update the prepared canonical `(y,x)` fringed, unfringed and reference images."""
function AdaptiveOpticsSim.Optics.propagate_compiled_optical_paths!(owner::PreparedCompiledPaths)
    owner.execution_context === nothing && throw(ArgumentError(
        "graph-owned compiled_paths propagation requires its caller's prepared execution context"))
    return Backends._with_prepared_device_execution_context(owner.execution_context) do
        products = _enqueue_compiled_paths!(owner)
        Backends._synchronize_prepared_device_execution_context!(owner.execution_context)
        products
    end
end


function _validate_compiled_paths_buffer(target,values,::Type{T},shape,role) where {T}
    eltype(values)===T || throw(ArgumentError("$role element type must be $T"))
    size(values)==shape || throw(DimensionMismatch("$role shape must be $shape"))
    Backends.compute_device(values)==target || throw(ArgumentError("$role occupies another compute device"))
    return values
end
