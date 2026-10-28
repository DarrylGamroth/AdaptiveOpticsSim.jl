# Polychromatic propagation reuses complete scalar trains. Each scalar
# intensity carries its cell-integrated power convention onto a common grid.
struct PolychromaticOpticalTrainNode{T} end

struct PolychromaticOpticalTrainPlan{P,C,W,R,S}
    reference::P
    channels::C
    weights::W
    transmission::R
    output_sampling_m::S
end

function _polychromatic_optical_train_node(::Type{T}, name, reference, config, bindings) where {T}
    channels = if hasproperty(config, :channels)
        AOG._proper_channel_configs(config)
    else
        ((wavelength_um=config.wavelength_um, weight=one(T), steps=config.steps),)
    end
    output_sampling_m = if hasproperty(config, :output_sampling_m)
        value = T(AOG._file_real(config.output_sampling_m, "output_sampling_m"))
        isfinite(value) && value > 0 || throw(AOG.AlgorithmGraphError(
            "output_sampling_m must be finite, positive and representable"))
        value
    else
        nothing
    end
    wavelengths = map(channels) do channel
        value = T(AOG._file_real(channel.wavelength_um, "channel wavelength_um"))
        isfinite(value) && value > 0 && isfinite(value*T(1e-6)) && value*T(1e-6) > 0 ||
            throw(AOG.AlgorithmGraphError("channel wavelength must be finite, positive and representable in metres"))
        value
    end
    weights = map(channels) do channel
        value = Float64(AOG._file_real(channel.weight, "channel weight"))
        isfinite(value) && value >= 0 ||
            throw(AOG.AlgorithmGraphError("channel weight must be finite and nonnegative"))
        value
    end
    total = sum(weights)
    isfinite(total) && total > 0 ||
        throw(AOG.AlgorithmGraphError("channel weights must have a finite positive sum"))
    transmission = map(channels) do channel
        value = hasproperty(channel, :transmission) ?
            AOG._file_real(channel.transmission, "channel transmission") : one(T)
        0 <= value <= 1 ||
            throw(AOG.AlgorithmGraphError("channel transmission must be in [0,1]"))
        T(value)
    end
    # The common AOS spectrum defines photon-number fractions, not energy
    # fractions. Convert the declared micrometres at this physical API boundary.
    samples = [AdaptiveOpticsSim.Optics.SpectralSample{T}(wavelengths[i]*T(1e-6),
        T(weights[i]/total)) for i in eachindex(channels)]
    spectrum = AdaptiveOpticsSim.Optics.SpectralBundle(samples)
    plans = map(channels, wavelengths) do channel, wavelength
        steps = map(AOG._file_tuple(channel.steps, "channel steps")) do step
            step isa NamedTuple || throw(AOG.AlgorithmGraphError("step must be a mapping"))
            _train_step(step, T, config, bindings)
        end
        OpticalTrainPlan(steps, reference.resolution, reference.pupil_resolution,
            reference.diameter_m, wavelength, reference.schemas, nothing, reference.output_quantity)
    end
    normalized_weights = Tuple(sample.weight for sample in spectrum)
    plan = PolychromaticOpticalTrainPlan(reference, plans, normalized_weights, transmission, output_sampling_m)
    return AOG.algorithm_node(name, PolychromaticOpticalTrainNode{T}, plan; props=NamedTuple())
end

function AOG.graph_node_ports(::Type{PolychromaticOpticalTrainNode{T}},
        plan::PolychromaticOpticalTrainPlan) where {T}
    return AOG.graph_node_ports(OpticalTrainNode{T}, plan.reference)
end

struct CellOverlap{T}
    source::Int
    destination::Int
    fraction::T
end

# Destination-grouped sparse overlaps; device kernels gather each output cell
# under a single writer. The sweep has O(n) entries, never a dense n×n map.
struct CellOverlapGather{I,F}
    offsets::I
    sources::I
    fractions::F
end

_prepare_cell_overlaps(::Nothing, n, target) = nothing
function _prepare_cell_overlaps(overlaps::AbstractVector{CellOverlap{T}}, n, target) where {T}
    target == Backends.HostComputeDevice() && return overlaps
    offsets = ones(Int, n + 1)
    for overlap in overlaps
        offsets[overlap.destination + 1] += 1
    end
    for i in 2:n+1
        offsets[i] += offsets[i - 1] - 1
    end
    cursor = copy(offsets)
    sources = Vector{Int}(undef, length(overlaps))
    fractions = Vector{T}(undef, length(overlaps))
    for overlap in overlaps
        i = cursor[overlap.destination]
        sources[i] = overlap.source
        fractions[i] = overlap.fraction
        cursor[overlap.destination] += 1
    end
    device_offsets = Backends.allocate_device_array(target, Int, length(offsets))
    device_sources = Backends.allocate_device_array(target, Int, length(sources))
    device_fractions = Backends.allocate_device_array(target, eltype(fractions), length(fractions))
    copyto!(device_offsets, offsets)
    copyto!(device_sources, sources)
    copyto!(device_fractions, fractions)
    return CellOverlapGather(device_offsets, device_sources, device_fractions)
end

Backends.KernelAbstractions.@kernel function _cell_rows!(scratch, intensity, offsets, sources, fractions)
    row, col = Backends.KernelAbstractions.@index(Global, NTuple)
    if row <= size(scratch, 1) && col <= size(scratch, 2)
        value = zero(eltype(scratch))
        for i in offsets[row]:offsets[row + 1]-1
            value += fractions[i] * intensity[sources[i], col]
        end
        scratch[row, col] = value
    end
end

Backends.KernelAbstractions.@kernel function _cell_columns!(output, scratch, offsets, sources, fractions, weight)
    row, col = Backends.KernelAbstractions.@index(Global, NTuple)
    if row <= size(output, 1) && col <= size(output, 2)
        value = output[row, col]
        for i in offsets[col]:offsets[col + 1]-1
            value += weight * fractions[i] * scratch[row, sources[i]]
        end
        output[row, col] = value
    end
end

# A source cell is uniform over [(i-center-1/2)*pitch,(i-center+1/2)*pitch].
# Fractions are overlap/source width; their tensor product conserves integrated
# source-cell power wherever the reference grid covers the source support.
function _cell_overlaps(n, source_pitch::T, destination_pitch::T) where {T}
    all(x -> isfinite(x) && x > 0, (source_pitch, destination_pitch)) ||
        throw(AOG.AlgorithmGraphError("polychromatic output sampling must be finite and positive"))
    source_pitch == destination_pitch && return nothing
    overlaps = CellOverlap{T}[]
    center = fld(n, 2)+1
    source, destination = 1, 1
    # Cold Float64 geometry reduces overlap roundoff for Float32 execution.
    sp, dp = Float64(source_pitch), Float64(destination_pitch)
    while source <= n && destination <= n
        slo, shi = (source-center-0.5)*sp, (source-center+0.5)*sp
        dlo, dhi = (destination-center-0.5)*dp, (destination-center+0.5)*dp
        overlap = min(shi, dhi)-max(slo, dlo)
        if overlap > 0
            push!(overlaps, CellOverlap(source, destination, T(overlap/sp)))
        end
        if shi <= dhi
            source += 1
        else
            destination += 1
        end
    end
    return FixedSizeVectorDefault{CellOverlap{T}}(overlaps)
end

struct PolychromaticOpticalTrainOwner{P,C,M,A,O,T}
    plan::P
    channels::C
    overlaps::M
    scratch::A
    output::O
    sampling_m::T
end

function AOG.prepare_graph_node(::Type{PolychromaticOpticalTrainNode{T}},
        plan::PolychromaticOpticalTrainPlan, ::NamedTuple{()},
        inputs::NamedTuple{(:pupil_opd,:pupil_amplitude)}, outputs::NamedTuple{(:output,)},
        ::NamedTuple{()}, target) where {T}
    n = plan.reference.resolution
    eltype(outputs.output) === T && size(outputs.output) == (n,n) &&
        Backends.compute_device(outputs.output) == target ||
        throw(AOG.AlgorithmGraphError("Proper optical train buffer contract mismatch"))
    reference_sampling = _prepare_masks(plan.reference, target).sampling_m
    output_sampling = something(plan.output_sampling_m, reference_sampling)
    prepared = map(plan.channels) do channel
        output = Backends.allocate_device_array(target, T, n, n)
        _prepare_optical_train(channel, inputs, (; output), target)
    end
    channels = map(p -> p.owner, prepared)
    overlaps = map(prepared) do channel
        _prepare_cell_overlaps(_cell_overlaps(n, channel.sampling_m, output_sampling), n, target)
    end
    scratch = Backends.allocate_device_array(target, T, n, n)
    fill!(scratch, 0)
    return PolychromaticOpticalTrainOwner(plan, channels, overlaps, scratch,
        outputs.output, output_sampling)
end

function _add_channel_intensity!(output, intensity, scratch, ::Nothing, weight)
    @. output += weight*intensity
    return nothing
end

function _add_channel_intensity!(output, intensity, scratch, overlaps, weight)
    fill!(scratch, zero(eltype(scratch)))
    for col in axes(intensity,2), overlap in overlaps
        scratch[overlap.destination,col] += overlap.fraction*intensity[overlap.source,col]
    end
    for overlap in overlaps, row in axes(output,1)
        output[row,overlap.destination] += weight*overlap.fraction*scratch[row,overlap.source]
    end
    return nothing
end

function _add_channel_intensity!(output, intensity, scratch, overlaps::CellOverlapGather, weight)
    backend = Backends.KernelAbstractions.get_backend(output)
    _cell_rows!(backend)(scratch, intensity, overlaps.offsets, overlaps.sources,
        overlaps.fractions; ndrange=size(output))
    _cell_columns!(backend)(output, scratch, overlaps.offsets, overlaps.sources,
        overlaps.fractions, weight; ndrange=size(output))
    return nothing
end

function AOG.step_graph_node!(owner::PolychromaticOpticalTrainOwner)
    fill!(owner.output, zero(eltype(owner.output)))
    foreach(owner.channels, owner.overlaps, owner.plan.weights,
            owner.plan.transmission) do channel, overlaps, weight, transmission
        AOG.step_graph_node!(channel)
        _add_channel_intensity!(owner.output, channel.output, owner.scratch, overlaps,
            weight*transmission)
    end
    return nothing
end

AOG.optical_sampling_m(owner::PolychromaticOpticalTrainOwner) = owner.sampling_m
function AOG.graph_node_capture_capability(owner::PolychromaticOpticalTrainOwner)
    return all(channel -> AOG.graph_node_capture_capability(channel) isa AOG.GraphNodeCaptureSafe,
        owner.channels) ? AOG.GraphNodeCaptureSafe() : AOG.GraphNodeCaptureUnsupported()
end

function AOG.reset_graph_node!(owner::PolychromaticOpticalTrainOwner)
    foreach(AOG.reset_graph_node!, owner.channels)
    fill!(owner.scratch, 0)
    fill!(owner.output, 0)
    return nothing
end
