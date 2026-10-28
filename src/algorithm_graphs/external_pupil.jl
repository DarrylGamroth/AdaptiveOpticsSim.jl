"""
    IntensityExposure

Single-writer CPU storage for duration-weighted relative intensity over exactly
covered model intervals. `image` is borrowed output storage; copy an emitted
image before retaining it across subsequent integrations.
"""
mutable struct IntensityExposure{A<:Matrix{Float64}}
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
    IntensityExposure(shape, duration::ModelDuration)

Create a CPU reference exposure with positive matrix dimensions and a positive
integer-nanosecond duration. Exposure boundaries align to model time zero.
"""
function IntensityExposure(shape::Tuple{Int,Int}, duration::ModelDuration)
    ns=model_nanoseconds(duration)
    ns>0 && all(>(0),shape) || throw(AlgorithmGraphError("exposure duration and shape must be positive"))
    return IntensityExposure(ns,0,0,0,zeros(Float64,shape),zeros(Float64,shape),UInt64(0),UInt64(0))
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
    fill!(exposure.integrated,0)
    return exposure
end

function _validate_intensity(intensity)
    all(intensity) do value
        value isa Real && isfinite(value) && isfinite(Float64(value)) && value >= 0
    end || throw(AlgorithmGraphError("intensity must be real, finite, nonnegative and representable as Float64"))
    return nothing
end

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
    (Base.mightalias(intensity,exposure.integrated) || Base.mightalias(intensity,exposure.image)) &&
        throw(AlgorithmGraphError("intensity input must not alias exposure storage"))
    _validate_intensity(intensity)
    start=model_nanoseconds(timestamp);span=model_nanoseconds(duration)
    span>0 || throw(AlgorithmGraphError("input coverage must be positive"))
    start==exposure.next_ns || throw(AlgorithmGraphError("intensity coverage must be contiguous; reset after a gap"))
    finish=Base.checked_add(start,span);cursor=start
    while cursor<finish
        boundary=Base.checked_add(exposure.start_ns,exposure.duration_ns)
        overlap=min(finish,boundary)-cursor
        overlap>0 || throw(AlgorithmGraphError("nonmonotonic exposure input"))
        exposure.integrated .+= intensity .* overlap
        exposure.covered_ns+=overlap;cursor+=overlap
        exposure.next_ns=cursor
        if cursor==boundary
            if exposure.covered_ns==exposure.duration_ns
                exposure.image .= exposure.integrated ./ exposure.duration_ns
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

"""Prepared CPU graph observations driven by coherent external OPD/amplitude samples.

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

"""
    prepare_external_pupil(graphs; exposure_durations, observations,
        input_names=(opd=:pupil_opd,amplitude=:pupil_amplitude))

Bind a nonempty named tuple of prepared CPU AOS graphs to independent intensity
exposures. `exposure_durations` contains `ModelDuration` values. Each observation
is `(output=:intensity,sampling_node=:optics)`, with caller-selected identifiers.
Pupil ports must be equally shaped Float64 matrices. No detector sampling,
photometry, resampling or pupil registration is inferred.
"""
function prepare_external_pupil(graphs::NamedTuple;exposure_durations::NamedTuple,
        observations::NamedTuple,input_names=(opd=:pupil_opd,amplitude=:pupil_amplitude))
    isempty(graphs) && throw(AlgorithmGraphError("external pupil requires an optical graph"))
    keys(graphs)==keys(exposure_durations)==keys(observations) || throw(AlgorithmGraphError("optical observation names differ"))
    opd=graph_input(first(values(graphs)),input_names.opd)
    amplitude=graph_input(first(values(graphs)),input_names.amplitude)
    opd isa Matrix{Float64} && amplitude isa Matrix{Float64} && size(opd)==size(amplitude) ||
        throw(AlgorithmGraphError("CPU external pupil requires equally shaped Float64 OPD and amplitude matrices"))
    for opd_graph in values(graphs),amplitude_graph in values(graphs)
        Base.mightalias(graph_input(opd_graph,input_names.opd),graph_input(amplitude_graph,input_names.amplitude)) &&
            throw(AlgorithmGraphError("external pupil OPD and amplitude buffers must not alias"))
    end
    for (name,graph) in pairs(graphs)
        intensity=graph_output(graph,getproperty(observations,name).output)
        for input_graph in values(graphs),input in values(input_names)
            Base.mightalias(intensity,graph_input(input_graph,input)) &&
                throw(AlgorithmGraphError("intensity observations must not alias pupil inputs"))
        end
        for (other_name,other_graph) in pairs(graphs)
            name==other_name && continue
            Base.mightalias(intensity,graph_output(other_graph,getproperty(observations,other_name).output)) &&
                throw(AlgorithmGraphError("independent intensity observations must not alias"))
        end
    end
    for (name,graph) in pairs(graphs)
        compute_device(graph)==HostComputeDevice() || throw(AlgorithmGraphError("external pupil reference requires CPU graphs"))
        for input in values(input_names)
            input_buffer=graph_input(graph,input)
            input_buffer isa Matrix{Float64} && size(input_buffer)==size(opd) || throw(AlgorithmGraphError("external pupil graph inputs differ"))
        end
        intensity=graph_output(graph,getproperty(observations,name).output)
        ndims(intensity)==2 || throw(AlgorithmGraphError("intensity observation must be a matrix"))
        prepared_graph_node(graph,getproperty(observations,name).sampling_node)
    end
    exposures=NamedTuple{keys(graphs)}(map(keys(graphs)) do name
        IntensityExposure(size(graph_output(getproperty(graphs,name),getproperty(observations,name).output)),
            getproperty(exposure_durations,name))
    end)
    return PreparedExternalPupil(graphs,exposures,opd,amplitude,observations,input_names,
        _ExternalPupilState(-1,0,0,nothing,0,0,0))
end

"""
    step_external_pupil!(emit, prepared, opd, amplitude, timestamp, duration;
        domain, generation, sequence, lost_samples=0)

Validate/copy a coherent external pupil, propagate every ordinary AOS graph,
then integrate intensity. The sample is piecewise constant over its declared
model interval and cannot exceed the shortest exposure. Sequence/time/loss gaps
or a newer generation invalidate partial exposures. Duplicate/overlapping and
older-generation samples are rejected with `false`. Different domains require
repreparation. `emit(name,image,metadata)` runs synchronously; copy its borrowed
image before retaining it. No transport or controller runs in this API.
"""
function step_external_pupil!(emit,prepared::PreparedExternalPupil,opd::AbstractMatrix,
        amplitude::AbstractMatrix,timestamp::ModelTimestamp,duration::ModelDuration;
        domain::NTuple{16,UInt8},generation::UInt64,sequence::UInt64,lost_samples::UInt64=UInt64(0))
    size(opd)==size(amplitude)==size(prepared.pupil_opd) || throw(DimensionMismatch("pupil sample shape"))
    Base.require_one_based_indexing(opd,amplitude)
    (Base.mightalias(opd,prepared.pupil_amplitude) || Base.mightalias(amplitude,prepared.pupil_opd)) &&
        throw(AlgorithmGraphError("external pupil samples must not alias the other bound quantity"))
    span=model_nanoseconds(duration);start=model_nanoseconds(timestamp)
    0<span<=minimum(e.duration_ns for e in values(prepared.exposures)) || throw(AlgorithmGraphError("sample coverage must fit shortest exposure"))
    generation>0 && sequence>0 && any(!iszero,domain) || throw(AlgorithmGraphError("invalid pupil sample identity"))
    finish=Base.checked_add(start,span);state=prepared.state
    if state.domain!==nothing && (domain!=state.domain || generation<state.generation ||
            (generation==state.generation && (sequence<=state.sequence || start<state.last_end_ns)))
        state.rejected+=1;return false
    end
    all(eachindex(opd,amplitude)) do index
        a=amplitude[index]
        a isa Real && opd[index] isa Real && isfinite(a) && isfinite(Float64(a)) && a>=0 &&
            (a==0 || (isfinite(opd[index]) && isfinite(Float64(opd[index]))))
    end || throw(AlgorithmGraphError("invalid pupil amplitude or transmitting OPD"))
    changed=generation!=state.generation
    if changed || start!=state.last_end_ns || sequence!=state.sequence+1 || lost_samples!=state.lost_samples
        state.discontinuities+=1
        for exposure in values(prepared.exposures)
            reset_intensity_exposure!(exposure,timestamp)
            changed && (exposure.sequence=0)
        end
        changed && foreach(reset_graph!,values(prepared.graphs))
    end
    for index in eachindex(prepared.pupil_opd)
        a=amplitude[index]
        prepared.pupil_amplitude[index]=a
        prepared.pupil_opd[index]=a==0 ? 0 : opd[index]
    end
    for (name,graph) in pairs(prepared.graphs)
        copyto!(graph_input(graph,prepared.input_names.opd),prepared.pupil_opd)
        copyto!(graph_input(graph,prepared.input_names.amplitude),prepared.pupil_amplitude)
        step_graph!(graph)
        intensity=graph_output(graph,getproperty(prepared.observations,name).output)
        _validate_intensity(intensity)
    end
    for (name,graph) in pairs(prepared.graphs)
        observation=getproperty(prepared.observations,name)
        sampling=optical_sampling_m(prepared_graph_node(graph,observation.sampling_node))
        integrate_intensity!(getproperty(prepared.exposures,name),graph_output(graph,observation.output),timestamp,duration) do image,pts,exposure,output_sequence
            emit(name,image,(;pts=model_nanoseconds(pts),duration=UInt64(model_nanoseconds(exposure)),
                generation,sequence=output_sequence,sampling_m=sampling))
        end
    end
    state.domain=domain;state.generation=generation;state.sequence=sequence
    state.last_end_ns=finish;state.lost_samples=lost_samples
    return true
end
