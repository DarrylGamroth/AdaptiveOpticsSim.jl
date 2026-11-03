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

struct PathPrefixCheckpoint
    index::Int
    cache_index::Int
end

struct PathPrefixPlan{C}
    cache_index::Int
    count::Int
    checkpoints::C
end

struct PathSharingPlan{S}
    prefix::S
    closed::S
    only::S
end

struct CompiledPathsPlan{T,D,P,F,A,S}
    target::D
    grid::Int
    image_size::Int
    crop_start::Int
    n_env::Int
    n_anchor::Int
    wavelengths::FixedSizeVectorDefault{T}
    paths::P
    shared_prefixes::S
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
    identity_roll::Bool
end

struct PathFFTWorkspace{A,F,I}
    field::A
    scratch::A
    forward::F
    inverse::I
end

# A plane is one fast path at one spectral anchor. Transforms cover only y,x;
# the third axis is an independent batch, never an optical coordinate.
struct PathFilterBatch{A,F,I}
    field::A
    transfer::A
    forward::F
    inverse::I
end

struct PathFastWorkspace{B,D,S}
    image::B
    defocus::D
    defocus_slots::S
end

struct PathNormalizationWorkspace{P,S}
    partials::P
    power::S
end

struct PathResampleBatch{N,A,F,O,C}
    batch_size::Val{N}
    field::A
    temporary::A
    kernel::A
    forward::F
    outputs::O
    coefficients::C
end

struct PathBranchSegment
    last::Int
    closed_cache::Int
    only_cache::Int
end

struct PathBranchPair{S}
    first::Int
    segments::S
end

struct PathBranchBatch{A,V,F,B,P}
    field::A
    planes::V
    forward::F
    backward::B
    pairs::P
end

struct CompiledPathsWorkspace{A,W,C,O,F,V,N,M,L,B,P,E,K,R,Q}
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
    fast::B
    propagation::P
    edge_multiplier::E
    prepared_mul_coefficients::K
    image_resample::R
    branch_batch::Q
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

# Pack a centred field into one corner-order FFT plane. Preparation establishes
# size(field,1:2) == size(input) .+ 2pad; the inside test bounds input accesses.
KernelAbstractions.@kernel function _compiled_paths_pack_filter_kernel!(field, input, slot, pad, enabled)
    y, x = @index(Global, NTuple)
    n = size(field, 1)
    sy = y + fld(n, 2)
    sx = x + fld(n, 2)
    sy > n && (sy -= n)
    sx > n && (sx -= n)
    inside = pad < sy <= size(input, 1) + pad && pad < sx <= size(input, 2) + pad
    @inbounds field[y, x, slot] = enabled && inside ? input[sy - pad, sx - pad] : zero(eltype(field))
end

# The defocus result remains corner ordered. Read its centred crop while
# packing the next image FFT; this combines crop, fftshift and ifftshift.
KernelAbstractions.@kernel function _compiled_paths_pack_defocus_kernel!(field, defocus, slot, defocus_slot, pad)
    y, x = @index(Global, NTuple)
    n, p = size(field, 1), size(defocus, 1)
    sy = y + fld(n, 2)
    sx = x + fld(n, 2)
    sy > n && (sy -= n)
    sx > n && (sx -= n)
    sy += pad + cld(p, 2)
    sx += pad + cld(p, 2)
    sy > p && (sy -= p)
    sx > p && (sx -= p)
    @inbounds field[y, x, slot] = defocus[sy, sx, defocus_slot]
end

# Equal positive square dimensions were checked during preparation. One wrap
# maps each output pixel to the centred inverse-transform result.
KernelAbstractions.@kernel function _compiled_paths_unpack_filter_kernel!(out, field, phase, slot)
    y, x = @index(Global, NTuple)
    n = size(out, 1)
    sy = y + cld(n, 2)
    sx = x + cld(n, 2)
    sy > n && (sy -= n)
    sx > n && (sx -= n)
    @inbounds out[y, x] = field[sy, sx, slot] * phase[y, x]
end

function _compiled_paths_pack_filter!(field, input, slot, pad, enabled=true)
    Backends.launch_kernel_async!(Backends.execution_style(field),
        _compiled_paths_pack_filter_kernel!, field, input, slot, pad, enabled;
        ndrange=(size(field, 1), size(field, 2)))
    return nothing
end

function _compiled_paths_pack_defocus!(field, defocus, slot, defocus_slot, pad)
    Backends.launch_kernel_async!(Backends.execution_style(field),
        _compiled_paths_pack_defocus_kernel!, field, defocus, slot, defocus_slot, pad;
        ndrange=(size(field, 1), size(field, 2)))
    return nothing
end

function _compiled_paths_unpack_filter!(out, field, phase, slot)
    Backends.launch_kernel_async!(Backends.execution_style(out),
        _compiled_paths_unpack_filter_kernel!, out, field, phase, slot; ndrange=size(out))
    return nothing
end

function _compiled_paths_pack_filter!(field::Array, input::Matrix, slot, pad, enabled=true)
    n = size(field, 1)
    for x in 1:n, y in 1:n
        sy, sx = mod1(y + fld(n, 2), n), mod1(x + fld(n, 2), n)
        inside = pad < sy <= size(input, 1) + pad && pad < sx <= size(input, 2) + pad
        field[y, x, slot] = enabled && inside ? input[sy - pad, sx - pad] : zero(eltype(field))
    end
    return nothing
end

function _compiled_paths_pack_defocus!(field::Array, defocus::Array, slot, defocus_slot, pad)
    n, p = size(field, 1), size(defocus, 1)
    for x in 1:n, y in 1:n
        sy = mod1(mod1(y + fld(n, 2), n) + pad + cld(p, 2), p)
        sx = mod1(mod1(x + fld(n, 2), n) + pad + cld(p, 2), p)
        field[y, x, slot] = defocus[sy, sx, defocus_slot]
    end
    return nothing
end

function _compiled_paths_unpack_filter!(out::Matrix, field::Array, phase::Matrix, slot)
    n = size(out, 1)
    for x in 1:n, y in 1:n
        sy, sx = mod1(y + cld(n, 2), n), mod1(x + cld(n, 2), n)
        out[y, x] = field[sy, sx, slot] * phase[y, x]
    end
    return nothing
end

function _compiled_paths_filter_batch!(batch::PathFilterBatch)
    Backends.enqueue_fft_plan!(batch.field, batch.forward)
    batch.field .*= batch.transfer
    Backends.enqueue_fft_plan!(batch.field, batch.inverse)
    return nothing
end
_compiled_paths_filter_batch!(::Nothing) = nothing

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

# Group independent numerical paths by exact root binding. Output slots retain
# their declared path/anchor indices; only the serial evaluation order changes.
# Completing each root before the next permits checkpoint storage to be reused.
function _compiled_paths_group_paths(paths)
    grouped=eltype(paths)[]
    visited=falses(length(paths))
    for first in eachindex(paths)
        visited[first] && continue
        reference=paths[first]
        for index in first:length(paths)
            path=paths[index]
            if !visited[index] && path.anchor_index==reference.anchor_index &&
                    path.entrance===reference.entrance
                push!(grouped,path)
                visited[index]=true
            end
        end
    end
    return grouped
end

function _compiled_paths_shared_prefixes(paths)
    empty = PathPrefixPlan(0,0,_compiled_paths_seal(PathPrefixCheckpoint[]))
    P = typeof(PathSharingPlan(empty,empty,empty))
    isempty(paths) && return _compiled_paths_seal(P[]),0
    O = eltype(first(paths).prefix)
    # Cold trie nodes identify exact numerical operations and their ancestry.
    # Roots additionally bind entrance identity and spectral anchor. No dictionary,
    # abstract operation container or mutable trie enters the prepared executor.
    operations = Union{Nothing,O}[]
    children = Vector{Int}[]
    origins, positions = Int[],Int[]
    roots, root_paths = Int[],Int[]
    chains = Vector{Int}[]
    function node!(op,stage,index)
        push!(operations,op);push!(children,Int[])
        push!(origins,stage);push!(positions,index)
        return length(operations)
    end
    function chain!(sequence,parent,stage)
        chain=Int[]
        for (index,op) in enumerate(sequence)
            found=findfirst(child->_compiled_paths_same_operation(operations[child]::O,op),children[parent])
            child=if found===nothing
                created=node!(op,stage,index)
                push!(children[parent],created)
                created
            else
                children[parent][found]
            end
            push!(chain,child)
            parent=child
        end
        push!(chains,chain)
        return parent
    end
    for (index,path) in enumerate(paths)
        found=findfirst(i->paths[root_paths[i]].anchor_index==path.anchor_index &&
            paths[root_paths[i]].entrance===path.entrance,eachindex(roots))
        root=if found===nothing
            created=node!(nothing,0,0)
            push!(roots,created);push!(root_paths,index)
            created
        else
            roots[found]
        end
        base=chain!(path.prefix,root,3index-2)
        chain!(path.closed,base,3index-1)
        chain!(path.only,base,3index)
    end
    counts=zeros(Int,length(chains))
    targets=Set{Int}()
    last_use=zeros(Int,length(operations))
    for (stage,chain) in enumerate(chains)
        count=0
        for (index,node) in enumerate(chain)
            origins[node]<stage || break
            count=index
        end
        counts[stage]=count
        if count>0
            node=chain[count]
            push!(targets,node)
            last_use[node]=stage
        end
    end
    cache_nodes=sort!(collect(targets))
    cache_indices=zeros(Int,length(operations))
    # Recycle storage only after its last consuming stage. Roots must finish
    # in this order; sibling batching must preserve initial cache reads and
    # disjoint checkpoint stores within these conservative stage intervals.
    slot_ends=Int[]
    for node in cache_nodes
        slot=findfirst(last->last<origins[node],slot_ends)
        if slot===nothing
            push!(slot_ends,last_use[node])
            slot=length(slot_ends)
        else
            slot_ends[slot]=last_use[node]
        end
        cache_indices[node]=slot
    end
    sequences=typeof(empty)[]
    for (stage,chain) in enumerate(chains)
        checkpoints=PathPrefixCheckpoint[
            PathPrefixCheckpoint(positions[node],cache_indices[node])
            for node in cache_nodes if origins[node]==stage]
        sort!(checkpoints;by=point->point.index)
        count=counts[stage]
        cache=count==0 ? 0 : cache_indices[chain[count]]
        push!(sequences,PathPrefixPlan(cache,count,_compiled_paths_seal(checkpoints)))
    end
    return _compiled_paths_seal(P[PathSharingPlan(sequences[3i-2],sequences[3i-1],sequences[3i])
        for i in eachindex(paths)]),length(slot_ends)
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

function _compiled_paths_filter_batch(target, ::Type{T}, n, transfers) where {T}
    isempty(transfers) && return nothing
    field = _compiled_paths_array(target, Complex{T}, n, n, length(transfers))
    transfer = similar(field)
    # ifftshift(fftshift(F) .* H) == F .* ifftshift(H). Keeping H in
    # corner order removes the cancelling shifts between each transform pair.
    host = Array{Complex{T}}(undef, size(field))
    for (slot, values) in enumerate(transfers)
        host[:, :, slot] .= ifftshift(Array(values))
    end
    copyto!(transfer, host)
    return PathFilterBatch(field, transfer,
        Backends.plan_fft_backend!(field, (1, 2)),
        Backends.plan_ifft_backend!(field, (1, 2)))
end

function _compiled_paths_fast_workspace(plan, ::Type{T}) where {T}
    image_transfers = [path.shift_transfer for path in plan.fast_paths for _ in 1:plan.n_anchor]
    defocus_transfers = [path.defocus_transfer[a] for path in plan.fast_paths
        for a in 1:plan.n_anchor if path.source_index != 0 && !isempty(path.defocus_transfer)]
    slots = zeros(Int, length(image_transfers))
    defocus_slot = 0
    for (k, path) in enumerate(plan.fast_paths), a in 1:plan.n_anchor
        if path.source_index != 0 && !isempty(path.defocus_transfer)
            defocus_slot += 1
            slots[_compiled_paths_slot(k, a, plan.n_anchor)] = defocus_slot
        end
    end
    return PathFastWorkspace(
        _compiled_paths_filter_batch(plan.target, T, plan.image_size, image_transfers),
        _compiled_paths_filter_batch(plan.target, T,
            plan.image_size + 2plan.phantom_pad, defocus_transfers),
        _compiled_paths_seal(slots))
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
    identity_roll = op.kind == PathResample && !isempty(op.auxiliary) &&
        all(value -> isequal(value, one(value)), Array(op.auxiliary))
    return PathOperationWorkspace(input, temporary, mapped, opd, electric_field, identity_roll)
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
    paths = _compiled_paths_group_paths(paths)
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
    coefficients = Backends._prepare_mul_coefficients(field)
    closed_sources, only_sources = fields(ne * na), fields(ne * na)
    batch_size = AdaptiveOpticsSim.Optics._compiled_paths_resample_batch_size(
        target, field, image_size, 2length(paths))
    image_resample = batch_size == 1 ? nothing : _compiled_paths_prepare_image_resample(
        plan, closed_sources, only_sources, coefficients, Val(batch_size))
    propagation = _compiled_paths_prepare_propagation(paths, wf, ctx)
    branch_batch = _compiled_paths_prepare_branch_batch(plan, field, propagation)
    workspace = CompiledPathsWorkspace(field, prefix, prefix_fields, wf, ctx, operations,
        physical_fft, image_fft, phantom_fft,
        _compiled_paths_array(target, Complex{T}, image_size, image_size),
        _compiled_paths_array(target, Complex{T}, image_size, image_size),
        fields(nk * na), fields(nk * na), fields(nk * na), closed_sources, only_sources,
        _compiled_paths_normalization_workspace(field, target, T), pixel_mtf,
        coefficients, _compiled_paths_fast_workspace(plan, T),
        propagation,
        _compiled_paths_prepare_edge_multiplier(paths, target, field),
        AdaptiveOpticsSim.Optics._prepare_compiled_paths_mul_coefficients(field, coefficients), image_resample, branch_batch)
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

# Fixed geometry is shared across optical branches. The cold dictionary never
# enters execution: the sealed vector has an explicit concrete union element.
# Enable this lowering on the CUDA backend qualified by the captured gates.
_compiled_paths_prepare_propagation(paths, wf, ctx) =
    _compiled_paths_prepare_propagation(Backends.array_backend_selector(typeof(wf.field)), paths, wf, ctx)
function _compiled_paths_prepare_propagation(::Backends.AbstractArrayBackend, paths, wf, ctx)
    count = sum(length(seq) for path in paths for seq in (path.prefix, path.closed, path.only))
    return _compiled_paths_seal(fill(nothing, count))
end
function _compiled_paths_prepare_propagation(::Backends.CUDABackend, paths, wf::Proper.WaveFront{T}, ctx) where {T}
    function prepare(op)
        wf.wavelength_m = op.wavelength_m
        wf.sampling_m = op.dx
        kind = op.kind == PathTransferFunction ? :ptp : op.from_waist ? :wts : :stw
        return Proper.prepare_propagation(wf, op.dz, ctx; kind)
    end
    first_prepared = Proper.prepare_propagation(wf, zero(T), ctx; kind=:ptp)
    key(op) = (op.kind, op.wavelength_m, op.dx, op.dz, op.from_waist)
    P = typeof(first_prepared)
    cache = Dict{Tuple{PathOperationKind,T,T,T,Bool},P}()
    values = Union{Nothing,P}[]
    for path in paths, seq in (path.prefix, path.closed, path.only), op in seq
        value = if op.kind in (PathTransferFunction, PathFourierTransform)
            get!(cache, key(op)) do
                prepare(op)
            end
        else
            nothing
        end
        push!(values, value)
    end
    return _compiled_paths_seal(values)
end

# Preserve the original complex coefficients for the retained numerical executor.
# A real, even multiplier can instead be read through a quarter-sized indexed
# view. Reconstruction is exact, including signed zeros; asymmetric or multiple
# distinct edges retain their original storage.
function _compiled_paths_prepare_edge_multiplier(paths, target, field)
    AdaptiveOpticsSim.Optics._compiled_paths_edge_cache_enabled(target, field) || return nothing
    edge = nothing
    for path in paths, sequence in (path.prefix, path.closed, path.only), op in sequence
        op.kind in (PathTransferFunction, PathFourierTransform) || continue
        isempty(op.auxiliary) && continue
        edge === nothing && (edge = op.auxiliary)
        edge === op.auxiliary || return nothing
    end
    edge === nothing && return nothing
    T = real(eltype(field))
    values = Array(edge)
    all(value -> isequal(imag(value), zero(T)), values) || return nothing
    ny, nx = size(values)
    iseven(ny) && iseven(nx) || return nothing
    iy = Int32[min(i, ny + 2 - i) for i in 1:ny]
    ix = Int32[min(i, nx + 2 - i) for i in 1:nx]
    quarter = real.(values[1:ny÷2+1, 1:nx÷2+1])
    isequal(complex.(quarter[iy, ix]), values) || return nothing
    resident = _compiled_paths_resident(target, quarter, T)
    rows = _compiled_paths_resident(target, iy, Int32)
    columns = ny == nx ? rows : _compiled_paths_resident(target, ix, Int32)
    return (edge, view(resident, rows, columns))
end
@inline _compiled_paths_edge(::Nothing, op) = isempty(op.auxiliary) ? nothing : op.auxiliary
@inline _compiled_paths_edge(cache::Tuple, op) =
    op.auxiliary === first(cache) ? last(cache) : _compiled_paths_edge(nothing, op)
@inline _compiled_paths_mul_coefficients(w, prepared::Bool) =
    prepared ? w.prepared_mul_coefficients : w.mul_coefficients

_compiled_paths_proper!(workspace, op, ::Nothing) = _compiled_paths_proper!(workspace, op)

function _compiled_paths_proper!(workspace, op, prepared::Proper.PreparedPropagation)
    edge = _compiled_paths_edge(workspace.edge_multiplier, op)
    Proper.propagate!(workspace.field, prepared, workspace.field; centered=true,
        input_multiplier=edge,
        output_multiplier=op.kind == PathFourierTransform ? op.gouy : one(op.gouy))
    return nothing
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

function _compiled_paths_proper_pair!(workspace, op, next, first::Proper.PreparedPropagation, second::Proper.PreparedPropagation)
    edge(value) = _compiled_paths_edge(workspace.edge_multiplier, value)
    gain(value) = value.kind == PathFourierTransform ? value.gouy : one(value.gouy)
    Proper.propagate!(workspace.field, (first,second), workspace.field; centered=true,
        input_multipliers=(edge(op),edge(next)), output_multipliers=(gain(op),gain(next)))
    return nothing
end
function _compiled_paths_proper_pair!(workspace, op, next, ::Nothing, ::Nothing)
    _compiled_paths_proper!(workspace, op)
    _compiled_paths_proper!(workspace, next)
    return nothing
end

function _compiled_paths_proper_group!(workspace, operations, masks,
    prepared::NTuple{N,P},input=workspace.field,initial_mask=nothing) where {N,P<:Proper.PreparedPropagation}
    edge(op)=_compiled_paths_edge(workspace.edge_multiplier, op)
    gain(op)=op.kind==PathFourierTransform ? op.gouy : one(op.gouy)
    inputs=ntuple(Val(N)) do i
        i==1 ? (initial_mask,edge(operations[i])) :
            isodd(i) ? (masks[fld(i-1,2)].array,edge(operations[i])) : edge(operations[i])
    end
    Proper.propagate!(workspace.field,prepared,input;centered=true,
        input_multipliers=inputs,output_multipliers=map(gain,operations))
    return nothing
end
function _compiled_paths_proper_group!(workspace,operations,masks,::NTuple{N,Nothing},
    input=workspace.field,initial_mask=nothing) where {N}
    input === workspace.field || Backends.copyto_backend_async!(workspace.field,input)
    initial_mask === nothing || (workspace.field .*=initial_mask)
    for i in 1:N
        if i>1 && isodd(i)
            workspace.field .*=masks[fld(i-1,2)].array
        end
        _compiled_paths_proper!(workspace,operations[i])
    end
    return nothing
end

function _compiled_paths_group_matches(sequence,index,last,::Val{N}) where {N}
    stop=index+3(N÷2)-2
    stop<=last || return false
    # Leave a final PTP and its Fourier resampling to the transform pair
    # lowering below. A checkpoint beyond `last` must still receive a field.
    stop<last && sequence[stop].kind==PathTransferFunction &&
        sequence[stop+1].kind==PathResample && return false
    for pair in 0:N÷2-1
        first=index+3pair
        sequence[first].kind in (PathTransferFunction,PathFourierTransform) &&
            sequence[first+1].kind in (PathTransferFunction,PathFourierTransform) || return false
        pair==N÷2-1 || sequence[first+2].kind==PathMultiply || return false
    end
    return true
end

function _compiled_paths_resample_frequency!(w,op,work,prepared::Bool=true)
    if !isempty(op.auxiliary) && !(prepared && work.identity_roll)
        w.physical_fft.scratch .*=op.auxiliary
    end
    Backends._enqueue_mul!(work.temporary,op.array,w.physical_fft.scratch,_compiled_paths_mul_coefficients(w,prepared))
    Backends._enqueue_mul!(w.field,work.temporary,transpose(op.array),_compiled_paths_mul_coefficients(w,prepared))
    normalization = size(op.array,2)^2
    if prepared && !isempty(op.post)
        @. w.field = (w.field / normalization) * op.post
    else
        w.field ./=normalization
        isempty(op.post) || (w.field .*=op.post)
    end
    return nothing
end
function _compiled_paths_resample!(w,op,work,prepared::Bool=true)
    _compiled_paths_fft!(w.physical_fft.scratch,w.field,w.physical_fft)
    return _compiled_paths_resample_frequency!(w,op,work,prepared)
end

function _compiled_paths_proper_resample!(w,op,tf,resample,work,
    first::Proper.PreparedPropagation,last::Proper.PreparedPropagation)
    if last.direction!=2
        _compiled_paths_proper_pair!(w,op,tf,first,last)
        return _compiled_paths_resample!(w,resample,work)
    end
    # PTP ends with BFFT followed by scale; resampling immediately applies FFT.
    # FFT(BFFT(G)*s)=length(G)*s*G. Preserve phases, edge/gain positions and the
    # centered Fourier coordinates; only the redundant FFT pair is removed.
    edge(value)=_compiled_paths_edge(w.edge_multiplier, value)
    gain=op.kind==PathFourierTransform ? op.gouy : one(op.gouy)
    backend=KernelAbstractions.get_backend(first.field)
    shift=fld(size(first.field,1),2)
    Proper._propagation_pack_kernel!(backend,(16,16))(
        first.field,w.field,edge(op),first.before,shift;ndrange=size(w.field))
    Proper._propagation_fft!(first)
    Proper._propagation_join_kernel!(backend,(16,16))(
        first.field,first.after,first.scale,gain,edge(tf),last.before,shift;ndrange=size(w.field))
    last.forward*last.field
    Proper._propagation_transfer_kernel!(backend,(16,16))(
        last.field,last.transfer;ndrange=size(w.field))
    Proper._propagation_unpack_kernel!(backend,(16,16))(
        w.physical_fft.scratch,last.field,last.after,length(last.field)*last.scale,
        one(tf.gouy),cld(size(last.field,1),2);ndrange=size(w.field))
    return _compiled_paths_resample_frequency!(w,resample,work)
end
function _compiled_paths_proper_resample!(w,op,tf,resample,work,::Nothing,::Nothing)
    _compiled_paths_proper_pair!(w,op,tf,nothing,nothing)
    return _compiled_paths_resample!(w,resample,work)
end
function _compiled_paths_execute_group!(owner,sequence,offset,index,::Val{N},
    input=owner.workspace.field,initial_mask=nothing) where {N}
    return _compiled_paths_execute_group!(owner,sequence,offset,index,Val(N),input,initial_mask,
        owner.workspace.propagation[offset+index])
end
function _compiled_paths_execute_group!(owner,sequence,offset,index,::Val{N},
    input,initial_mask,::P) where {N,P}
    positions=ntuple(i->index+3fld(i-1,2)+mod(i-1,2),Val(N))
    masks=ntuple(i->sequence[index+3i-1],Val(N÷2-1))
    # Preparation seals one propagation type (or all Nothing). Narrow it before
    # constructing the tuple, rather than boxing Union-valued tuple elements.
    prepared=ntuple(i->owner.workspace.propagation[offset+positions[i]]::P,Val(N))
    _compiled_paths_proper_group!(owner.workspace,map(i->sequence[i],positions),masks,
        prepared,input,initial_mask)
    return index+3(N÷2)-2
end

# Consume an initial sampled multiplier while packing the first transform.
# The source may be the current field or a retained prefix checkpoint.
function _compiled_paths_multiplied_group!(owner,sequence,offset,first,last,input)
    first+2<=last && sequence[first].kind==PathMultiply &&
        sequence[first+1].kind in (PathTransferFunction,PathFourierTransform) &&
        owner.workspace.propagation[offset+first+1] !== nothing || return first-1
    return _compiled_paths_multiplied_group!(owner,sequence,offset,first,last,input,
        (Val(10),Val(8),Val(6),Val(4),Val(2)))
end
@inline _compiled_paths_multiplied_group!(owner,sequence,offset,first,last,input,::Tuple{}) = first-1
@inline function _compiled_paths_multiplied_group!(owner,sequence,offset,first,last,input,groups::Tuple)
    group=Base.first(groups)
    if _compiled_paths_group_matches(sequence,first+1,last,group)
        return _compiled_paths_execute_group!(owner,sequence,offset,first+1,group,input,sequence[first].array)
    end
    return _compiled_paths_multiplied_group!(owner,sequence,offset,first,last,input,Base.tail(groups))
end

function _compiled_paths_interpolate!(work, kernel, source, coefficients)
    work.input .= source
    Backends._enqueue_mul!(work.temporary, kernel, work.input, coefficients)
    Backends._enqueue_mul!(work.mapped, work.temporary, transpose(kernel), coefficients)
    work.opd .= real.(work.mapped)
    return work.opd
end

function _compiled_paths_execute!(owner, sequence, offset, first=1, last=length(sequence), prepared::Bool=true)
    w = owner.workspace
    index = first
    while index <= last
        op = sequence[index]
        work = w.operations[offset + index]
        if op.kind == PathMultiply
            consumed = prepared ? _compiled_paths_multiplied_group!(owner,sequence,offset,index,last,w.field) : index-1
            if consumed < index
                w.field .*= op.array
            else
                index = consumed
            end
        elseif op.kind == PathNormalize
            _compiled_paths_normalize!(w.field, w.normalization)
        elseif op.kind in (PathTransferFunction, PathFourierTransform)
            if prepared && index+2<=last && sequence[index+1].kind==PathTransferFunction &&
                    sequence[index+2].kind==PathResample
                _compiled_paths_proper_resample!(w,op,sequence[index+1],sequence[index+2],
                    w.operations[offset+index+2],w.propagation[offset+index],w.propagation[offset+index+1])
                index+=2
            elseif prepared && _compiled_paths_group_matches(sequence,index,last,Val(10))
                index=_compiled_paths_execute_group!(owner,sequence,offset,index,Val(10))
            elseif prepared && _compiled_paths_group_matches(sequence,index,last,Val(8))
                index=_compiled_paths_execute_group!(owner,sequence,offset,index,Val(8))
            elseif prepared && _compiled_paths_group_matches(sequence,index,last,Val(6))
                index=_compiled_paths_execute_group!(owner,sequence,offset,index,Val(6))
            elseif prepared && _compiled_paths_group_matches(sequence,index,last,Val(4))
                index=_compiled_paths_execute_group!(owner,sequence,offset,index,Val(4))
            elseif prepared && index < last && sequence[index+1].kind in (PathTransferFunction, PathFourierTransform)
                _compiled_paths_proper_pair!(w, op, sequence[index+1],
                    w.propagation[offset+index], w.propagation[offset+index+1])
                index += 1
            elseif prepared
                _compiled_paths_proper!(w, op, w.propagation[offset + index])
            else
                _compiled_paths_proper!(w, op)
            end
        elseif op.kind == PathSpectralMultiplier
            _compiled_paths_filter!(w.field, w.field, op.array, w.physical_fft)
        elseif op.kind == PathResample
            _compiled_paths_resample!(w,op,work,prepared)
        elseif op.kind == PathLiveScreen
            input = op.source == 1 ? owner.inputs.dm_opd : owner.inputs.entrance_opd
            _compiled_paths_interpolate!(work, op.array, input, _compiled_paths_mul_coefficients(w,prepared))
            AdaptiveOpticsSim.Optics.apply_phase_async!(work.electric_field, work.opd)
            if op.source == 2
                _compiled_paths_interpolate!(work, op.array, owner.inputs.entrance_attenuation,
                    _compiled_paths_mul_coefficients(w,prepared))
                @. w.field *= exp(-work.opd)
            end
        end
        index += 1
    end
    return offset + length(sequence)
end

@inline _compiled_paths_slot(path, anchor, n_anchor) = (path - 1) * n_anchor + anchor

function _compiled_paths_sequence!(owner,sequence,offset,sharing,share_prefix,prepared)
    if !share_prefix
        return _compiled_paths_execute!(owner,sequence,offset,1,length(sequence),prepared)
    end
    w=owner.workspace
    first=sharing.count+1
    if sharing.count>0
        input=w.shared_prefixes[sharing.cache_index]
        last=isempty(sharing.checkpoints) ? length(sequence) : sharing.checkpoints[1].index
        consumed=prepared ? _compiled_paths_multiplied_group!(owner,sequence,offset,first,last,input) : first-1
        if consumed < first
            Backends.copyto_backend_async!(w.field,input)
        else
            first=consumed+1
        end
    end
    for checkpoint in sharing.checkpoints
        _compiled_paths_execute!(owner,sequence,offset,first,checkpoint.index,prepared)
        Backends.copyto_backend_async!(w.shared_prefixes[checkpoint.cache_index],w.field)
        first=checkpoint.index+1
    end
    return _compiled_paths_execute!(owner,sequence,offset,first,length(sequence),prepared)
end

# Inputs remain independent through the resampling FFT and both matrix
# products. Batch outputs are consumed only after all slow paths complete.
function _compiled_paths_prepare_image_resample(plan, closed, only, coefficients, ::Val{N}) where {N}
    N > 0 || throw(ArgumentError("compiled image resampling batch must be positive"))
    n = plan.image_size
    field = _compiled_paths_array(plan.target, eltype(plan.resample_kernel), n, n, N)
    temporary = similar(field)
    # fftshift(F) = P*F*transpose(P). Copy K*P by permuting K's columns,
    # so both products accept raw FFT order. Coefficients are not reevaluated;
    # the changed dot-product summation order can change roundoff.
    columns = mod1.((1:n) .+ fld(n,2), n)
    values = reshape(Array(plan.resample_kernel)[:,columns], n, n, 1)
    kernel = _compiled_paths_resident(plan.target, values, eltype(field))
    outputs = eltype(closed)[]
    for path in plan.paths
        slot = _compiled_paths_slot(path.path_index, path.anchor_index, plan.n_anchor)
        push!(outputs, closed[slot], only[slot])
    end
    return PathResampleBatch(Val(N), field, temporary, kernel,
        Backends.plan_fft_backend!(field, (1,2)), _compiled_paths_seal(outputs), coefficients)
end

# Preparation checks crop bounds. The wrapped coordinates stay in 1:n.
# A final partial batch clears unused planes to prevent stale values growing
# on replay; transforms and products never combine separate batch planes.
KernelAbstractions.@kernel function _compiled_paths_pack_resample_kernel!(out, input, start, slot, valid)
    y,x = @index(Global,NTuple)
    n = size(out,1)
    if y <= n && x <= n
        sy,sx = y + fld(n,2), x + fld(n,2)
        sy > n && (sy -= n)
        sx > n && (sx -= n)
        @inbounds out[y,x,slot] = input[sy+start,sx+start]
        for plane in valid+1:size(out,3)
            @inbounds out[y,x,plane] = zero(eltype(out))
        end
    end
end
KernelAbstractions.@kernel function _compiled_paths_unpack_resample_kernel!(outputs, field, valid, normalization)
    y,x,plane = @index(Global,NTuple)
    if y <= size(field,1) && x <= size(field,2) && plane <= valid
        @inbounds outputs[plane][y,x] = field[y,x,plane] / normalization
    end
end

function _compiled_paths_stage_image!(::Nothing, owner, number, out, prepared)
    w,plan = owner.workspace, owner.plan
    _compiled_paths_crop!(w.crop,w.field,plan.crop_start)
    _compiled_paths_resample_image!(out,w.crop,owner,prepared)
    return nothing
end
function _compiled_paths_stage_image!(batch::PathResampleBatch{N}, owner, number, out, prepared) where {N}
    prepared || return _compiled_paths_stage_image!(nothing,owner,number,out,prepared)
    n = owner.plan.image_size
    slot = mod1(number,N)
    valid = number == length(batch.outputs) ? slot : N
    backend = KernelAbstractions.get_backend(batch.field)
    _compiled_paths_pack_resample_kernel!(backend,(16,16))(
        batch.field, owner.workspace.field, owner.plan.crop_start, slot, valid; ndrange=(n,n))
    if slot == N || number == length(batch.outputs)
        Backends.enqueue_fft_plan!(batch.field,batch.forward)
        AdaptiveOpticsSim.Optics._compiled_paths_resample_mul!(batch.field,batch.temporary,
            batch.kernel,batch.coefficients)
        first = number - slot
        outputs = ntuple(i -> batch.outputs[min(first+i,length(batch.outputs))], Val(N))
        _compiled_paths_unpack_resample_kernel!(backend,(16,16,1))(
            outputs, batch.field, slot, n^2; ndrange=size(batch.field))
    end
    return nothing
end

# Match independent suffixes without changing the existing prefix trie or its
# cache allocation. Restrict cuts to common multiplier/transform-pair regions.
# Other sequences retain their ordinary executor.
function _compiled_paths_branch_pair(path, sharing, propagation, offset)
    a,b = path.closed,path.only
    length(a) == length(b) && !isempty(a) && sharing.closed.count == sharing.only.count || return nothing
    a[end].kind == b[end].kind == PathResample || return nothing
    _compiled_paths_same_operation(a[end],b[end]) || return nothing
    left,right = sharing.closed.checkpoints,sharing.only.checkpoints
    length(left) == length(right) || return nothing
    all(i -> left[i].index == right[i].index,eachindex(left)) || return nothing
    # Only may depend on a checkpoint produced by this closed suffix. It must
    # then wait for that producer in the original serial order.
    sharing.only.count > 0 && any(p -> p.cache_index == sharing.only.cache_index,left) && return nothing
    any(x -> any(y -> x.cache_index == y.cache_index,right),left) && return nothing
    segments = PathBranchSegment[PathBranchSegment(p.index,p.cache_index,right[i].cache_index)
        for (i,p) in enumerate(left)]
    push!(segments,PathBranchSegment(length(a)-1,0,0))
    first = sharing.closed.count+1
    for segment in segments
        last = segment.last
        last >= first && (last-first+1)%3 == 0 || return nothing
        for i in first:3:last
            a[i].kind == b[i].kind == PathMultiply || return nothing
            for j in i+1:i+2
                a[j].kind in (PathTransferFunction,PathFourierTransform) &&
                    _compiled_paths_same_operation(a[j],b[j]) || return nothing
                p,q = propagation[offset+j],propagation[offset+length(a)+j]
                p !== nothing && p === q || return nothing
            end
        end
        first = last+1
    end
    a[end-1].kind == PathTransferFunction || return nothing
    propagation[offset+length(a)-1].direction == 2 || return nothing
    return PathBranchPair(sharing.closed.count+1,_compiled_paths_seal(segments))
end

function _compiled_paths_branch_pairs(plan, propagation)
    S = typeof(_compiled_paths_seal(PathBranchSegment[]))
    pairs = Union{Nothing,PathBranchPair{S}}[]
    offset = 0
    for (i,path) in enumerate(plan.paths)
        offset += length(path.prefix)
        push!(pairs,_compiled_paths_branch_pair(path,plan.shared_prefixes[i],propagation,offset))
        offset += length(path.closed)+length(path.only)
    end
    return _compiled_paths_seal(pairs)
end
function _compiled_paths_branch_batch(field, pairs)
    n = size(field,1)
    planes = similar(field,n,n,2)
    fill!(planes,zero(eltype(planes)))
    return PathBranchBatch(planes,(view(planes,:,:,1),view(planes,:,:,2)),
        Backends.plan_fft_backend!(planes,(1,2)),
        Backends.plan_bfft_backend!(planes,(1,2)),pairs)
end
function _compiled_paths_prepare_branch_batch(plan,field,propagation)
    AdaptiveOpticsSim.Optics._compiled_paths_branch_batch_enabled(plan.target,field) || return nothing
    pairs = _compiled_paths_branch_pairs(plan,propagation)
    all(isnothing,pairs) && return nothing
    return _compiled_paths_branch_batch(field,pairs)
end
@inline _compiled_paths_branch_pair(::Nothing,index) = nothing
@inline _compiled_paths_branch_pair(batch::PathBranchBatch,index) = batch.pairs[index]
@inline _compiled_paths_branch_gain(op) = op.kind == PathFourierTransform ? op.gouy : one(op.gouy)

function _compiled_paths_branch_fft!(batch, propagation, terminal::Bool=false)
    if propagation.direction == 2
        Backends.enqueue_fft_plan!(batch.field,batch.forward)
        n = size(batch.field,1)
        Proper._propagation_batch_transfer_kernel!(KernelAbstractions.get_backend(batch.field),(16,16))(
            batch.field,propagation.transfer;ndrange=(n,n))
        terminal || Backends.enqueue_fft_plan!(batch.field,batch.backward)
    elseif propagation.direction == 1
        Backends.enqueue_fft_plan!(batch.field,batch.forward)
    elseif propagation.direction == -1
        Backends.enqueue_fft_plan!(batch.field,batch.backward)
    end
    return nothing
end

function _compiled_paths_branches!(owner,path,sharing,offset,index,pair)
    w = owner.workspace
    batch = w.branch_batch
    n = size(batch.field,1)
    backend = KernelAbstractions.get_backend(batch.field)
    a,b = path.closed,path.only
    other = offset+length(a)
    # Both initial cache reads precede either suffix's new stores. Stage-level
    # cache lifetimes therefore remain conservative; paired cache outputs are
    # disjoint and no current closed producer feeds the paired only input.
    inputs = (sharing.closed.count > 0 ? w.shared_prefixes[sharing.closed.cache_index] : w.field,
        sharing.only.count > 0 ? w.shared_prefixes[sharing.only.cache_index] : w.field)
    first = pair.first
    for segment in pair.segments
        initial = first+1
        previous = w.propagation[offset+initial]
        previous_op = a[initial]
        Proper._propagation_batch_pack_kernel!(backend,(16,16))(
            batch.field,inputs,(a[first].array,b[first].array),
            _compiled_paths_edge(w.edge_multiplier,previous_op),previous.before,fld(n,2);ndrange=(n,n))
        _compiled_paths_branch_fft!(batch,previous)
        for j in initial+1:segment.last
            a[j].kind == PathMultiply && continue
            next = w.propagation[offset+j]
            masks = a[j-1].kind == PathMultiply ? (a[j-1].array,b[j-1].array) : (nothing,nothing)
            Proper._propagation_batch_join_kernel!(backend,(16,16))(
                batch.field,previous.after,previous.scale,_compiled_paths_branch_gain(previous_op),
                masks,_compiled_paths_edge(w.edge_multiplier,a[j]),next.before,fld(n,2);ndrange=(n,n))
            _compiled_paths_branch_fft!(batch,next,j == length(a)-1)
            previous,previous_op = next,a[j]
        end
        if segment.closed_cache == 0
            # Preserve FFT(BFFT(G)*s)=length(G)*s*G at the terminal PTP.
            # Consume one plane at a time through the existing scalar mapping
            # scratch; neither resampling GEMM writes the paired raw field.
            slot = _compiled_paths_slot(path.path_index,path.anchor_index,owner.plan.n_anchor)
            for plane in 1:2
                Proper._propagation_unpack_kernel!(backend,(16,16))(
                    w.physical_fft.scratch,batch.planes[plane],previous.after,
                    n^2*previous.scale,one(previous_op.gouy),cld(n,2);ndrange=(n,n))
                operation = plane == 1 ? a[end] : b[end]
                work = w.operations[(plane == 1 ? offset : other)+length(a)]
                _compiled_paths_resample_frequency!(w,operation,work)
                output = plane == 1 ? w.closed_sources[slot] : w.only_sources[slot]
                _compiled_paths_stage_image!(w.image_resample,owner,2index-2+plane,output,true)
            end
        else
            outputs = (w.shared_prefixes[segment.closed_cache],w.shared_prefixes[segment.only_cache])
            Proper._propagation_batch_unpack_kernel!(backend,(16,16))(
                outputs,batch.field,previous.after,previous.scale,
                _compiled_paths_branch_gain(previous_op),cld(n,2);ndrange=(n,n))
            inputs = outputs
        end
        first = segment.last+1
    end
    return other+length(b)
end

function _compiled_paths_slow!(owner,share_prefix::Bool=true,prepared::Bool=true)
    plan,w=owner.plan,owner.workspace
    offset=0
    for (index,path) in enumerate(plan.paths)
        sharing=plan.shared_prefixes[index]
        # Each reused prefix is produced earlier in this invocation, including
        # during capture. Scratch offsets remain those of the original program;
        # cached fields do not restore wavefront geometry or FFT bindings.
        if !share_prefix || sharing.prefix.count==0
            Backends.copyto_backend_async!(w.field,path.entrance)
        end
        offset=_compiled_paths_sequence!(owner,path.prefix,offset,sharing.prefix,share_prefix,prepared)
        pair = prepared && share_prefix ? _compiled_paths_branch_pair(w.branch_batch,index) : nothing
        if pair !== nothing
            offset = _compiled_paths_branches!(owner,path,sharing,offset,index,pair)
            continue
        end
        if !share_prefix || sharing.only.count==0
            Backends.copyto_backend_async!(w.prefix,w.field)
        end
        offset=_compiled_paths_sequence!(owner,path.closed,offset,sharing.closed,share_prefix,prepared)
        slot=_compiled_paths_slot(path.path_index,path.anchor_index,plan.n_anchor)
        _compiled_paths_stage_image!(w.image_resample,owner,2index-1,w.closed_sources[slot],prepared)
        if !share_prefix || sharing.only.count==0
            Backends.copyto_backend_async!(w.field,w.prefix)
        end
        offset=_compiled_paths_sequence!(owner,path.only,offset,sharing.only,share_prefix,prepared)
        _compiled_paths_stage_image!(w.image_resample,owner,2index,w.only_sources[slot],prepared)
    end
    for a in 1:plan.n_anchor
        w.closed_sources[a] .*=plan.main_phase
        w.only_sources[a] .*=plan.main_phase
    end
    return nothing
end

function _compiled_paths_resample_image!(out, input, owner, prepared::Bool=true)
    plan, w = owner.plan, owner.workspace
    _compiled_paths_fft!(w.image_fft.scratch, input, w.image_fft)
    Backends._enqueue_mul!(w.resample_temporary, plan.resample_kernel, w.image_fft.scratch,
        _compiled_paths_mul_coefficients(w,prepared))
    Backends._enqueue_mul!(out, w.resample_temporary, transpose(plan.resample_kernel),
        _compiled_paths_mul_coefficients(w,prepared))
    out ./= plan.image_size^2
    return out
end

function _compiled_paths_fast_scalar!(envelopes, sources, owner, free_ghost)
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

function _compiled_paths_fast!(envelopes, sources, owner, free_ghost)
    plan, w = owner.plan, owner.workspace
    fast = w.fast
    # All bindings and batch slots are frozen at preparation. Closed and only
    # branches reuse the same scratch in stream order, including capture replay.
    for (k, path) in enumerate(plan.fast_paths), a in 1:plan.n_anchor
        slot = _compiled_paths_slot(k, a, plan.n_anchor)
        defocus_slot = fast.defocus_slots[slot]
        if defocus_slot != 0
            source = sources[_compiled_paths_slot(path.source_index, a, plan.n_anchor)]
            _compiled_paths_pack_filter!(fast.defocus.field, source, defocus_slot, plan.phantom_pad)
        end
    end
    _compiled_paths_filter_batch!(fast.defocus)
    for (k, path) in enumerate(plan.fast_paths), a in 1:plan.n_anchor
        slot = _compiled_paths_slot(k, a, plan.n_anchor)
        defocus_slot = fast.defocus_slots[slot]
        if defocus_slot != 0
            _compiled_paths_pack_defocus!(fast.image.field, fast.defocus.field,
                slot, defocus_slot, plan.phantom_pad)
        elseif path.source_index == 0
            _compiled_paths_pack_filter!(fast.image.field, path.ghost[a], slot, 0, free_ghost)
        else
            source = sources[_compiled_paths_slot(path.source_index, a, plan.n_anchor)]
            _compiled_paths_pack_filter!(fast.image.field, source, slot, 0)
        end
    end
    _compiled_paths_filter_batch!(fast.image)
    for (k, path) in enumerate(plan.fast_paths), a in 1:plan.n_anchor
        slot = _compiled_paths_slot(k, a, plan.n_anchor)
        _compiled_paths_unpack_filter!(envelopes[slot], fast.image.field, path.phase_multiplier[a], slot)
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

function _enqueue_compiled_paths!(owner::PreparedCompiledPaths, share_prefix::Bool=true, batch_fast::Bool=true, prepared::Bool=true)
    plan, w, products = owner.plan, owner.workspace, owner.products
    _compiled_paths_slow!(owner, share_prefix, prepared)
    if batch_fast
        _compiled_paths_fast!(w.closed_envelopes, w.closed_sources, owner, true)
        _compiled_paths_fast!(w.only_envelopes, w.only_sources, owner, false)
    else
        _compiled_paths_fast_scalar!(w.closed_envelopes, w.closed_sources, owner, true)
        _compiled_paths_fast_scalar!(w.only_envelopes, w.only_sources, owner, false)
    end
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
