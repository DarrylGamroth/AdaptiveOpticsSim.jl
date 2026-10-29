"""
    external_pupil_contract(; pupil_resolution, diameter_m, plane, curvature,
        registration="provisional", spectral_model="grey_pupil_amplitude",
        wavelength_um=nothing)

Return the transport-neutral `schema`, string `properties`, and `(2,p,p)`
`shape` of a Float32 external pupil packet. Logical axes are `(quantity,y,x)`;
quantities are OPD in metres and relative field amplitude. The origin is at
`fld(p,2)+1`, and sample spacing is `diameter_m/p` metres. Native AOS pupil
arrays use `[x,y]`: publishers copy `opd[x,y]` to `packet[1,y,x]` and amplitude
to `packet[2,y,x]`. Storage layout and transport belong to the caller.

The grey OPD/amplitude pair is shared by every receiving wavelength. The
collimated pupil excludes nominal convergence, which the receiving optical
train introduces. Registration is descriptive metadata and does not change
the existing schema hash. This contract does not infer a relay registration
or flux calibration. `wavelength_um` preserves the original monochromatic
packet contract for existing callers; omit it for a grey pupil.
"""
function external_pupil_contract(; pupil_resolution, diameter_m, plane, curvature,
        registration="provisional", spectral_model="grey_pupil_amplitude", wavelength_um=nothing)
    pupil_resolution isa Integer && !(pupil_resolution isa Bool) &&
        0 < pupil_resolution <= typemax(Int) ||
        throw(ArgumentError("external pupil_resolution must be a positive Int"))
    diameter_m isa Real && !(diameter_m isa Bool) && isfinite(diameter_m) && diameter_m > 0 ||
        throw(ArgumentError("external pupil diameter_m must be positive and finite"))
    all(value -> value isa AbstractString && !isempty(strip(value)) && !occursin('\n', value),
        (plane, curvature, registration)) ||
        throw(ArgumentError("external pupil plane, curvature and registration require nonempty single-line text"))
    curvature == "collimated" || throw(ArgumentError("unsupported external pupil curvature"))
    spectral_model == "grey_pupil_amplitude" || throw(ArgumentError("unsupported external pupil spectral_model"))
    if wavelength_um !== nothing
        wavelength_um isa Real && !(wavelength_um isa Bool) && isfinite(wavelength_um) && wavelength_um > 0 ||
            throw(ArgumentError("external pupil wavelength_um must be positive and finite"))
    end
    p = Int(pupil_resolution)
    spectral_contract = wavelength_um === nothing ? spectral_model : wavelength_um
    # Keep field order and numeric spelling identical to the original CAL2
    # publisher/receiver contract, independently of YAML mapping order.
    contract = join(("org.adaptiveopticssim.pupil-field/1", p, diameter_m, spectral_contract,
        plane, curvature, "quantity,y,x", "opd_m,amplitude",
        "relative_field_amplitude", "input_owned", "simulation_nanoseconds", "piecewise_constant"), "\n")
    digest = bytes2hex(sha256(contract))
    schema = "org.adaptiveopticssim.pupil-field/1;sha256=" * digest
    properties = Dict("optics.contract-sha256" => digest,
        "optics.sampling-m" => string(diameter_m/p), "optics.plane" => String(plane),
        "optics.axis-order" => "quantity,y,x", "optics.quantities" => "opd_m,amplitude",
        "optics.amplitude-units" => "relative_field_amplitude", "optics.pts-units" => "simulation_nanoseconds",
        "optics.registration" => String(registration))
    if wavelength_um === nothing
        properties["optics.spectral-model"] = String(spectral_model)
    else
        properties["optics.wavelength-m"] = string(wavelength_um*1e-6)
    end
    return (; schema, properties, shape=(2,p,p))
end

"""
    IntensityExposure

Single-writer storage for duration-weighted relative intensity over exactly
covered model intervals. `image` is borrowed output storage; copy an emitted
image before retaining it across subsequent integrations.
"""
mutable struct IntensityExposure{A<:AbstractMatrix{Float64}}
    duration_ns::Int64
    start_ns::Int64
    covered_ns::Int64
    next_ns::Int64
    integrated::A
    image::A
    sequence::UInt64
    incomplete::UInt64
end

"""
    IntensityExposure(shape, duration::ModelDuration; target=HostComputeDevice())

Create an exposure on the declared compute device with positive matrix dimensions
and a positive integer-nanosecond duration. Accumulation uses Float64; exposure
boundaries align to model time zero.
"""
function IntensityExposure(shape::Tuple{Int,Int}, duration::ModelDuration;
        target::AbstractComputeDevice=HostComputeDevice())
    ns=model_nanoseconds(duration)
    ns>0 && all(>(0),shape) || throw(AlgorithmGraphError("exposure duration and shape must be positive"))
    integrated=allocate_device_array(target,Float64,shape...)
    image=allocate_device_array(target,Float64,shape...)
    Backends._with_compute_device(target) do
        fill!(integrated,0);fill!(image,0)
    end
    return IntensityExposure(ns,0,0,0,integrated,image,UInt64(0),UInt64(0))
end

"""
    reset_intensity_exposure!(exposure, time::ModelTimestamp)

Discard partial coverage and align to the exposure containing `time`. Preserve
the emitted sequence and last borrowed image; integration resumes at `time`.
"""
function reset_intensity_exposure!(exposure::IntensityExposure,time::ModelTimestamp)
    exposure.covered_ns>0 && (exposure.incomplete+=1)
    exposure.start_ns=fld(model_nanoseconds(time),exposure.duration_ns)*exposure.duration_ns
    exposure.covered_ns=0
    exposure.next_ns=model_nanoseconds(time)
    Backends._with_compute_device(compute_device(exposure.integrated)) do
        fill!(exposure.integrated,0)
    end
    return exposure
end

function _validate_intensity(intensity)
    all(intensity) do value
        value isa Real && isfinite(value) && isfinite(Float64(value)) && value >= 0
    end || throw(AlgorithmGraphError("intensity must be real, finite, nonnegative and representable as Float64"))
    return nothing
end

_scaled_intensity_valid(::Type{T},intensity,scale) where {T} =
    all(value -> isfinite(T(value*scale)),intensity)

"""
    integrate_intensity!(emit, exposure, intensity, timestamp, duration)

Integrate piecewise-constant relative intensity over exact model-time overlaps.
`intensity` must be a one-based matrix matching the exposure shape, with real,
finite, nonnegative values representable as finite Float64 values. It must not
alias `exposure.integrated` or `exposure.image`, including through a view.
Validate these input contracts before changing exposure state.

Coverage must be positive and contiguous; call `reset_intensity_exposure!` after
a gap. For each complete exposure, synchronously call
`emit(image, timestamp, duration, sequence)` with its duration-weighted mean.
The emitted image is borrowed and may change again within this same call; copy
it before retaining it. Return `exposure`.
"""
function integrate_intensity!(emit,exposure::IntensityExposure,intensity::AbstractMatrix,
        timestamp::ModelTimestamp,duration::ModelDuration)
    size(intensity)==size(exposure.image) || throw(DimensionMismatch("intensity exposure shape"))
    Base.require_one_based_indexing(intensity)
    compute_device(intensity)==compute_device(exposure.image) ||
        throw(AlgorithmGraphError("intensity and exposure must use the same compute device"))
    (Base.mightalias(intensity,exposure.integrated) || Base.mightalias(intensity,exposure.image)) &&
        throw(AlgorithmGraphError("intensity input must not alias exposure storage"))
    return Backends._with_compute_device(compute_device(exposure.image)) do
        _integrate_intensity!(emit,exposure,intensity,timestamp,duration)
    end
end

function _integrate_intensity!(emit,exposure,intensity,timestamp,duration)
    _validate_intensity(intensity)
    start=model_nanoseconds(timestamp);span=model_nanoseconds(duration)
    span>0 || throw(AlgorithmGraphError("input coverage must be positive"))
    start==exposure.next_ns || throw(AlgorithmGraphError("intensity coverage must be contiguous; reset after a gap"))
    finish=Base.checked_add(start,span);cursor=start
    while cursor<finish
        boundary=Base.checked_add(exposure.start_ns,exposure.duration_ns)
        overlap=min(finish,boundary)-cursor
        overlap>0 || throw(AlgorithmGraphError("nonmonotonic exposure input"))
        exposure.integrated .+= Float64.(intensity) .* overlap
        exposure.covered_ns+=overlap;cursor+=overlap
        exposure.next_ns=cursor
        if cursor==boundary
            if exposure.covered_ns==exposure.duration_ns
                exposure.image .= exposure.integrated ./ exposure.duration_ns
                # A consumer may copy the borrowed image to the host or submit
                # a detector graph on another stream. Publish completed pixels.
                Backends.synchronize_backend!(execution_style(exposure.image))
                exposure.sequence+=1
                emit(exposure.image,ModelTimestamp(exposure.start_ns),ModelDuration(exposure.duration_ns),exposure.sequence)
            else
                exposure.incomplete+=1
            end
            exposure.start_ns=boundary;exposure.covered_ns=0
            fill!(exposure.integrated,0)
        end
    end
    return exposure
end

mutable struct _ExternalPupilState
    last_end_ns::Int64
    generation::UInt64
    sequence::UInt64
    domain::Union{Nothing,NTuple{16,UInt8}}
    lost_samples::UInt64
    discontinuities::UInt64
    rejected::UInt64
end

"""Prepared graph observations driven by coherent host OPD/amplitude samples.

The graphs use ordinary AOS preparation, stepping and reset. `observations` binds
one intensity output and sampling node per graph; input names bind the two
pupil quantities explicitly. Transport adapters retain ownership of packet
formats and convert their metadata to model timestamps/durations.
"""
struct PreparedExternalPupil{G,E,A,O,I}
    graphs::G
    exposures::E
    pupil_opd::A
    pupil_amplitude::A
    observations::O
    input_names::I
    state::_ExternalPupilState
end

_observation_response(observation) = get(observation, :response, nothing)

_observation_graphs(graph, observation) = haskey(observation, :chopping) ?
    (graph, observation.chopping.unfringed) : (graph,)

function _validate_chopping(graph, observation, input_names)
    haskey(observation, :chopping) || return nothing
    chopping = observation.chopping
    chopping isa NamedTuple || throw(AlgorithmGraphError("chopping must be a named tuple"))
    _require_named_fields(chopping, (:unfringed, :initial_phase),
        (:unfringed, :initial_phase), "chopping")
    chopping.initial_phase in (:fringed, :unfringed) ||
        throw(AlgorithmGraphError("initial_phase must be :fringed or :unfringed"))
    other = chopping.unfringed
    other isa PreparedAlgorithmGraph || throw(AlgorithmGraphError("unfringed must be a prepared graph"))
    compute_device(other) == compute_device(graph) ||
        throw(AlgorithmGraphError("chopping graphs must use the same compute device"))
    output, alternate = graph_output(graph, observation.output), graph_output(other, observation.output)
    size(output) == size(alternate) && eltype(output) == eltype(alternate) ||
        throw(AlgorithmGraphError("chopping output shape and element type must match"))
    for input in values(input_names)
        a, b = graph_input(graph, input), graph_input(other, input)
        size(a) == size(b) && eltype(a) == eltype(b) ||
            throw(AlgorithmGraphError("chopping pupil input contracts must match"))
    end
    return nothing
end

function _validate_chopping_sampling(graph, observation)
    haskey(observation, :chopping) || return nothing
    other = observation.chopping.unfringed
    # Sampling is established by propagation, including all spectral channels.
    try
        step_graph!(graph); step_graph!(other)
        first_sampling = optical_sampling_m(prepared_graph_node(graph, observation.sampling_node))
        other_sampling = optical_sampling_m(prepared_graph_node(other, observation.sampling_node))
        isfinite(first_sampling) && first_sampling > 0 && isfinite(other_sampling) && other_sampling > 0 &&
            isapprox(first_sampling, other_sampling; rtol=8eps(Float32), atol=0) ||
            throw(AlgorithmGraphError("chopping optical sampling must match"))
    finally
        reset_graph!(graph); reset_graph!(other)
    end
    return nothing
end

function _chopped_image_kind(chopping, start, duration)
    even = iseven(fld(start, duration))
    return even == (chopping.initial_phase === :fringed) ? :fringed : :unfringed
end

function _active_observation_graphs(graph, observation, start, finish, duration)
    haskey(observation, :chopping) || return (graph,)
    chopping = observation.chopping
    first = _chopped_image_kind(chopping, start, duration) === :fringed ? graph : chopping.unfringed
    last = _chopped_image_kind(chopping, finish-1, duration) === :fringed ? graph : chopping.unfringed
    return first === last ? (first,) : (first, last)
end

function _validate_observation_response(response, graph, observation, duration)
    response === nothing && return nothing
    response isa NamedTuple || throw(AlgorithmGraphError("exposure response must be a named tuple"))
    _require_named_fields(response, (:graph, :input, :output, :normalized_to_photon_rate,
        :pixel_sampling_m, :exposure_duration_s, :frame_schema),
        (:graph, :input, :output, :normalized_to_photon_rate,
        :pixel_sampling_m, :exposure_duration_s, :frame_schema), "exposure response")
    response.frame_schema isa AbstractString && !isempty(strip(response.frame_schema)) ||
        throw(AlgorithmGraphError("exposure response requires a nonempty frame_schema"))
    compute_device(response.graph) == compute_device(graph) ||
        throw(AlgorithmGraphError("external pupil response must use the optical graph compute device"))
    input = graph_input(response.graph, response.input)
    output = graph_output(response.graph, response.output)
    intensity = graph_output(graph, observation.output)
    input isa AbstractMatrix{<:Union{Float32,Float64}} && size(input) == size(intensity) ||
        throw(AlgorithmGraphError("response photon-rate input must match the optical grid; no resampling is inferred"))
    output isa AbstractMatrix{<:Union{Float32,Float64}} && size(output) == size(intensity) ||
        throw(AlgorithmGraphError("response frame must match the optical grid"))
    Base.mightalias(input, output) && throw(AlgorithmGraphError("response input and output must not alias"))
    scale = response.normalized_to_photon_rate
    scale isa Real && !(scale isa Bool) && isfinite(scale) && scale > 0 ||
        throw(AlgorithmGraphError("normalized_to_photon_rate must be finite and positive"))
    sampling = response.pixel_sampling_m
    sampling isa Tuple && length(sampling) == 2 &&
        all(x -> x isa Real && !(x isa Bool) && isfinite(x) && x > 0, sampling) ||
        throw(AlgorithmGraphError("response pixel_sampling_m requires positive row/column pitches"))
    seconds = response.exposure_duration_s
    seconds isa Real && !(seconds isa Bool) && isfinite(seconds) && seconds > 0 &&
        isapprox(Float64(seconds) * 1e9, model_nanoseconds(duration); rtol=8eps(Float32), atol=0.5) ||
        throw(AlgorithmGraphError("response detector duration must match the integrated exposure"))
    # PROPER reports output sampling only after propagation. Probe once while
    # cold, then reset; this does not consume an acquisition or detector draw.
    step_graph!(graph)
    optical_sampling = optical_sampling_m(prepared_graph_node(graph, observation.sampling_node))
    all(pitch -> isapprox(pitch, optical_sampling; rtol=8eps(Float32), atol=0), sampling) ||
        throw(AlgorithmGraphError("response pixel sampling differs from optical sampling; declare a resampling model explicitly"))
    reset_graph!(graph)
    return nothing
end

"""
    prepare_external_pupil(graphs; exposure_durations, observations,
        input_names=(opd=:pupil_opd,amplitude=:pupil_amplitude))

Bind a nonempty named tuple of prepared AOS graphs to independent intensity
exposures. `exposure_durations` contains `ModelDuration` values. Each observation
is `(output=:intensity,sampling_node=:optics)`, with caller-selected identifiers.
All optical and response graphs must use the same compute device. Pupil ports
must be equally shaped Float32 or Float64 matrices, with one common precision
across every optical branch and chopping variant. Validated external host samples are
staged into each graph; propagation, integration and detector response remain
on its compute device. Emitted images are borrowed arrays on that device.
No detector sampling,
photometry, resampling or pupil registration is inferred.

An observation may additionally declare `response=(graph=..., input=:photon_rate,
output=:frame, normalized_to_photon_rate=..., pixel_sampling_m=(dy,dx),
exposure_duration_s=..., frame_schema=...)`. The explicit scale converts each relative intensity
sample to cell-integrated photons/s; it includes the caller's radiometric and
spatial calibration and is never normalized by the current image sum. A response
uses the same optical grid and runs once per complete exposure. Preparation
probes and resets the optical graph to check sampling. Detector readout and noise
therefore run at the exposure cadence, after optical integration. Responses reset
after an input discontinuity or a new generation. Required `frame_schema` labels
the output quantity without inferring electrons or ADC units.

An observation may declare `chopping=(unfringed=prepared_graph,
initial_phase=:fringed)`. The ordinary graph supplies the fringed intensity;
the separate static graph blocks the reference field in that observation only.
Both variants must have matching pupil/output arrays, compute device and optical
sampling, with independent mutable graph state and boundary storage. The caller
owns their physical optical meaning and calibration. No mask is inferred or
changed by this API. `initial_phase` is `:fringed` or `:unfringed` for the exposure
starting at model time zero. The kind alternates at every model exposure boundary,
including incomplete or discarded intervals, and is published in completed
metadata as `chopped_image_kind`. It is not electrical trigger polarity or pair
completion. A sample crossing a boundary uses the corresponding graph on each
side before integration; other observations retain their original sampling.
One shared detector response runs once per completed chopped exposure.
"""
function prepare_external_pupil(graphs::NamedTuple;exposure_durations::NamedTuple,
        observations::NamedTuple,input_names=(opd=:pupil_opd,amplitude=:pupil_amplitude))
    isempty(graphs) && throw(AlgorithmGraphError("external pupil requires an optical graph"))
    keys(graphs)==keys(exposure_durations)==keys(observations) || throw(AlgorithmGraphError("optical observation names differ"))
    for (name, graph) in pairs(graphs)
        _validate_chopping(graph, getproperty(observations, name), input_names)
    end
    all_graphs = Tuple(variant for (name, graph) in pairs(graphs)
        for variant in _observation_graphs(graph, getproperty(observations, name)))
    # Chopping variants are independent graph owners, not alternate aliases of
    # one graph's mutable execution state or borrowed boundary buffers.
    if any(observation -> haskey(observation, :chopping), values(observations))
        buffers = Tuple((graph_input(variant, input_names.opd),
            graph_input(variant, input_names.amplitude),
            graph_output(variant, getproperty(observations, name).output))
            for (name, graph) in pairs(graphs)
            for variant in _observation_graphs(graph, getproperty(observations, name)))
        for i in eachindex(all_graphs), j in 1:i-1
            a, b = all_graphs[i], all_graphs[j]
            (a === b || a.state === b.state ||
                any(x -> any(y -> x.owner === y.owner, values(b.nodes)), values(a.nodes))) &&
                throw(AlgorithmGraphError("chopping observations must own independent graph state"))
            any(x -> any(y -> Base.mightalias(x, y), buffers[j]), buffers[i]) &&
                throw(AlgorithmGraphError("chopping observation buffers must not alias"))
        end
    end
    opd=graph_input(first(values(graphs)),input_names.opd)
    amplitude=graph_input(first(values(graphs)),input_names.amplitude)
    target=compute_device(first(values(graphs)))
    T = eltype(opd)
    T in (Float32,Float64) && opd isa AbstractMatrix{T} && amplitude isa AbstractMatrix{T} && size(opd)==size(amplitude) ||
        throw(AlgorithmGraphError("external pupil requires equally shaped Float32 or Float64 OPD and amplitude matrices with matching precision"))
    for opd_graph in all_graphs,amplitude_graph in all_graphs
        Base.mightalias(graph_input(opd_graph,input_names.opd),graph_input(amplitude_graph,input_names.amplitude)) &&
            throw(AlgorithmGraphError("external pupil OPD and amplitude buffers must not alias"))
    end
    for (name,base_graph) in pairs(graphs), graph in _observation_graphs(base_graph,getproperty(observations,name))
        intensity=graph_output(graph,getproperty(observations,name).output)
        for input_graph in all_graphs,input in values(input_names)
            Base.mightalias(intensity,graph_input(input_graph,input)) &&
                throw(AlgorithmGraphError("intensity observations must not alias pupil inputs"))
        end
        for (other_name,other_base) in pairs(graphs), other_graph in _observation_graphs(other_base,getproperty(observations,other_name))
            name==other_name && continue
            Base.mightalias(intensity,graph_output(other_graph,getproperty(observations,other_name).output)) &&
                throw(AlgorithmGraphError("independent intensity observations must not alias"))
        end
    end
    for (name,base_graph) in pairs(graphs), graph in _observation_graphs(base_graph,getproperty(observations,name))
        compute_device(graph)==target || throw(AlgorithmGraphError("external pupil graphs must use the same compute device"))
        for input in values(input_names)
            input_buffer=graph_input(graph,input)
            input_buffer isa AbstractMatrix{T} && size(input_buffer)==size(opd) || throw(AlgorithmGraphError("external pupil graph inputs differ in shape or precision"))
        end
        intensity=graph_output(graph,getproperty(observations,name).output)
        ndims(intensity)==2 || throw(AlgorithmGraphError("intensity observation must be a matrix"))
        prepared_graph_node(graph,getproperty(observations,name).sampling_node)
        _validate_observation_response(_observation_response(getproperty(observations,name)),
            graph, getproperty(observations,name), getproperty(exposure_durations,name))
    end
    for (name, graph) in pairs(graphs)
        _validate_chopping_sampling(graph, getproperty(observations, name))
    end
    for (name, observation) in pairs(observations)
        response = _observation_response(observation)
        response === nothing && continue
        buffers = (graph_input(response.graph, response.input), graph_output(response.graph, response.output))
        for other in all_graphs, buffer in buffers
            any(input -> Base.mightalias(buffer, graph_input(other, input)), values(input_names)) &&
                throw(AlgorithmGraphError("response storage must not alias pupil inputs"))
        end
        for (other_name, other_observation) in pairs(observations),
                other_graph in _observation_graphs(getproperty(graphs, other_name), other_observation)
            other_response = _observation_response(other_observation)
            for buffer in buffers
                Base.mightalias(buffer, graph_output(other_graph, other_observation.output)) &&
                    throw(AlgorithmGraphError("response storage must not alias optical outputs"))
                (other_name == name || other_response === nothing) && continue
                any(other -> Base.mightalias(buffer, other),
                    (graph_input(other_response.graph, other_response.input), graph_output(other_response.graph, other_response.output))) &&
                    throw(AlgorithmGraphError("independent responses must not share storage"))
            end
        end
    end
    exposures=NamedTuple{keys(graphs)}(map(keys(graphs)) do name
        IntensityExposure(size(graph_output(getproperty(graphs,name),getproperty(observations,name).output)),
            getproperty(exposure_durations,name);target)
    end)
    # PipeWire packets are host-owned. Keep validation and sanitation outside
    # accelerator scalar indexing, with one retained staging pair per owner.
    if target!=HostComputeDevice()
        opd=zeros(T,size(opd));amplitude=zeros(T,size(amplitude))
    end
    return PreparedExternalPupil(graphs,exposures,opd,amplitude,observations,input_names,
        _ExternalPupilState(-1,0,0,nothing,0,0,0))
end

function _reset_observation_graphs!(prepared)
    for (name, graph) in pairs(prepared.graphs),
            variant in _observation_graphs(graph, getproperty(prepared.observations, name))
        reset_graph!(variant)
    end
    return nothing
end

function _integrate_observation!(emit, prepared, name, graph, observation,
        timestamp, duration, generation, kind=nothing)
    sampling = optical_sampling_m(prepared_graph_node(graph, observation.sampling_node))
    integrate_intensity!(getproperty(prepared.exposures,name),
            graph_output(graph,observation.output),timestamp,duration) do image,pts,exposure,output_sequence
        metadata = (;pts=model_nanoseconds(pts),duration=UInt64(model_nanoseconds(exposure)),
            generation,sequence=output_sequence,sampling_m=sampling)
        kind === nothing || (metadata = merge(metadata,(;chopped_image_kind=kind)))
        response = _observation_response(observation)
        if response === nothing
            emit(name,image,metadata)
        else
            input = graph_input(response.graph, response.input)
            _with_prepared_device_execution_context(response.graph.context) do
                input .= image .* response.normalized_to_photon_rate
            end
            step_graph!(response.graph)
            emit(name,graph_output(response.graph, response.output),
                merge(metadata,(;quantity=:detector_frame,frame_schema=response.frame_schema)))
        end
    end
    return nothing
end

"""
    step_external_pupil!(emit, prepared, opd, amplitude, timestamp, duration;
        domain, generation, sequence, lost_samples=0)

Validate/copy a coherent external host pupil, propagate every ordinary AOS graph,
then integrate intensity. Validate that host values convert to finite values at
the prepared pupil precision before changing exposure or sequence state; blocked
OPD is sanitized to zero. The sample is piecewise constant over its declared
model interval and cannot exceed the shortest exposure. Sequence/time/loss gaps
or a newer generation invalidate partial exposures. Duplicate/overlapping and
older-generation samples are rejected with `false`. Different domains require
repreparation. `emit(name,image,metadata)` runs synchronously; copy its borrowed
image before retaining it. No transport or controller runs in this API.
Chopped observations classify the actual model-time exposure, not receive or
emission parity. A gap discards partial coverage without moving the phase origin.
"""
function step_external_pupil!(emit,prepared::PreparedExternalPupil,opd::AbstractMatrix,
        amplitude::AbstractMatrix,timestamp::ModelTimestamp,duration::ModelDuration;
        domain::NTuple{16,UInt8},generation::UInt64,sequence::UInt64,lost_samples::UInt64=UInt64(0))
    size(opd)==size(amplitude)==size(prepared.pupil_opd) || throw(DimensionMismatch("pupil sample shape"))
    Base.require_one_based_indexing(opd,amplitude)
    compute_device(opd)==compute_device(amplitude)==HostComputeDevice() ||
        throw(AlgorithmGraphError("external pupil samples must be host matrices"))
    (Base.mightalias(opd,prepared.pupil_amplitude) || Base.mightalias(amplitude,prepared.pupil_opd)) &&
        throw(AlgorithmGraphError("external pupil samples must not alias the other bound quantity"))
    span=model_nanoseconds(duration);start=model_nanoseconds(timestamp)
    0<span<=minimum(e.duration_ns for e in values(prepared.exposures)) || throw(AlgorithmGraphError("sample coverage must fit shortest exposure"))
    generation>0 && sequence>0 && any(!iszero,domain) || throw(AlgorithmGraphError("invalid pupil sample identity"))
    start >= 0 || !any(observation -> haskey(observation,:chopping),values(prepared.observations)) ||
        throw(AlgorithmGraphError("chopping requires nonnegative model timestamps"))
    finish=Base.checked_add(start,span);state=prepared.state
    if state.domain!==nothing && (domain!=state.domain || generation<state.generation ||
            (generation==state.generation && (sequence<=state.sequence || start<state.last_end_ns)))
        state.rejected+=1;return false
    end
    T = eltype(prepared.pupil_opd)
    all(eachindex(opd,amplitude)) do index
        a=amplitude[index]
        a isa Real && opd[index] isa Real && isfinite(a) && isfinite(T(a)) && a>=0 &&
            (a==0 || (isfinite(opd[index]) && isfinite(T(opd[index]))))
    end || throw(AlgorithmGraphError("pupil amplitude or transmitting OPD is invalid at the selected optical precision"))
    changed=generation!=state.generation
    if changed || start!=state.last_end_ns || sequence!=state.sequence+1 || lost_samples!=state.lost_samples
        state.discontinuities+=1
        for exposure in values(prepared.exposures)
            reset_intensity_exposure!(exposure,timestamp)
            changed && (exposure.sequence=0)
        end
        changed && _reset_observation_graphs!(prepared)
        for observation in values(prepared.observations)
            response = _observation_response(observation)
            response === nothing || reset_graph!(response.graph)
        end
    end
    for index in eachindex(prepared.pupil_opd)
        a=amplitude[index]
        prepared.pupil_amplitude[index]=a
        prepared.pupil_opd[index]=a==0 ? 0 : opd[index]
    end
    for (name,base_graph) in pairs(prepared.graphs), graph in _active_observation_graphs(
            base_graph,getproperty(prepared.observations,name),start,finish,
            getproperty(prepared.exposures,name).duration_ns)
        _with_prepared_device_execution_context(graph.context) do
            copyto!(graph_input(graph,prepared.input_names.opd),prepared.pupil_opd)
            copyto!(graph_input(graph,prepared.input_names.amplitude),prepared.pupil_amplitude)
        end
        step_graph!(graph)
        intensity=graph_output(graph,getproperty(prepared.observations,name).output)
        _with_prepared_device_execution_context(graph.context) do
            _validate_intensity(intensity)
            response = _observation_response(getproperty(prepared.observations,name))
            if response !== nothing
                target_type = eltype(graph_input(response.graph, response.input))
                scale = response.normalized_to_photon_rate
                _scaled_intensity_valid(target_type,intensity,scale) ||
                    throw(AlgorithmGraphError("scaled photon rate exceeds response input precision"))
            end
        end
    end
    for (name,graph) in pairs(prepared.graphs)
        observation=getproperty(prepared.observations,name)
        if !haskey(observation,:chopping)
            _integrate_observation!(emit,prepared,name,graph,observation,timestamp,duration,generation)
            continue
        end
        exposure=getproperty(prepared.exposures,name)
        cursor=start
        while cursor<finish
            boundary=Base.checked_add(exposure.start_ns,exposure.duration_ns)
            stop=min(finish,boundary)
            kind=_chopped_image_kind(observation.chopping,cursor,exposure.duration_ns)
            selected=kind===:fringed ? graph : observation.chopping.unfringed
            _integrate_observation!(emit,prepared,name,selected,observation,
                ModelTimestamp(cursor),ModelDuration(stop-cursor),generation,kind)
            cursor=stop
        end
    end
    state.domain=domain;state.generation=generation;state.sequence=sequence
    state.last_end_ns=finish;state.lost_samples=lost_samples
    return true
end

"""
    reset_external_pupil!(prepared)

Reset optical and response histories, partial exposures, emitted sequences and
input-domain admission. The sole owner must stop and serialize input processing
before calling this function. The next sample establishes a fresh domain and
generation; no previous partial exposure or detector history participates.
Both chopped graph variants reset; their phase origin remains model time zero.
Transport buffers and source lifecycle remain the caller's responsibility.
"""
function reset_external_pupil!(prepared::PreparedExternalPupil)
    _reset_observation_graphs!(prepared)
    for observation in values(prepared.observations)
        response = _observation_response(observation)
        response === nothing || reset_graph!(response.graph)
    end
    for exposure in values(prepared.exposures)
        reset_intensity_exposure!(exposure, ModelTimestamp(0))
        exposure.sequence = 0
        exposure.incomplete = 0
        Backends._with_compute_device(compute_device(exposure.image)) do
            fill!(exposure.image, 0)
        end
    end
    fill!(prepared.pupil_opd, 0)
    fill!(prepared.pupil_amplitude, 0)
    state = prepared.state
    state.last_end_ns = -1
    state.generation = state.sequence = 0
    state.domain = nothing
    state.lost_samples = state.discontinuities = state.rejected = 0
    return prepared
end
