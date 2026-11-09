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
    execution = _compiled_paths_prepare_execution(plan, propagation)
    owner = PreparedCompiledPaths(plan, workspace, inputs, products, execution_context, synthesis, execution)
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
