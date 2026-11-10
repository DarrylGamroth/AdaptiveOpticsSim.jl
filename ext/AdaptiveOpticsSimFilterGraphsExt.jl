module AdaptiveOpticsSimFilterGraphsExt

import AdaptiveOpticsSim
import FilterGraphs as FG
import FilterGraphAlgorithms
const FGA = FilterGraphAlgorithms
using FilterGraphAlgorithms: @algorithm
import KernelAbstractions: get_backend
using KernelAbstractions: @atomicreplace, @kernel

const AOG = AdaptiveOpticsSim.AlgorithmGraphs
const B = AdaptiveOpticsSim.Backends
const O = AdaptiveOpticsSim.Optics
proper_extension() = Base.get_extension(AdaptiveOpticsSim, :AdaptiveOpticsSimProperExt)

@enum OpticalFailure::UInt8 InvalidOpticalInput=1

@algorithm OpticalObservationF32 begin
    label = "aos-optical-observation-f32"
    failure = OpticalFailure

    @configuration OpticalObservationConfiguration begin
        observation::AOG.PreparedAlgorithmGraph => graph_rebuild
        single_scc_product::Bool => graph_rebuild
        merge_scc_fields::Bool => graph_rebuild
    end

    ports(c) = (
        input(:pupil_opd, Float32, size(AOG.graph_input(c.observation, :pupil_opd)),
            "org.adaptiveopticssim.pupil-opd/1"),
        input(:pupil_amplitude, Float32, size(AOG.graph_input(c.observation, :pupil_amplitude)),
            "org.adaptiveopticssim.pupil-amplitude/1"),
        output(:intensity, Float32, size(AOG.graph_output(c.observation, :intensity)),
            "org.adaptiveopticssim.optical-intensity/1"; metadata_from=:pupil_opd),
    )

    prepare(c) = (; c.observation, c.single_scc_product, c.merge_scc_fields)
    process!(plan) = propagate!(plan.observation, intensity, pupil_opd, pupil_amplitude)
end

function propagate!(observation, intensity, opd, amplitude)
    copyto!(AOG.graph_input(observation, :pupil_opd), opd)
    copyto!(AOG.graph_input(observation, :pupil_amplitude), amplitude)
    foreach(node -> AOG.step_graph_node!(node.owner), values(observation.nodes))
    copyto!(intensity, AOG.graph_output(observation, :intensity))
    return nothing
end

FGA.capture_support(::Type{OpticalObservationF32}, device) = optical_capture_support(device)

function optical_capture_support(device)
    for (extension, backend, device_type) in (
        (:AdaptiveOpticsSimCUDAExt, :CUDA, :CuDevice),
        (:AdaptiveOpticsSimAMDGPUExt, :AMDGPU, :HIPDevice),
    )
        ext = Base.get_extension(AdaptiveOpticsSim, extension)
        ext === nothing && continue
        device isa getproperty(getproperty(ext, backend), device_type) && return FGA.SupportsCapture()
    end
    return FGA.NoCapture()
end

struct CapturedObservation{G,I,O,S,K,V,F}
    observation::G
    inputs::I
    intensity::O
    status::S
    check!::K
    validation::V
    propagate!::F
    invalid::Int32
end

@kernel function check_optical_status!(status, validation, invalid)
    if any(!iszero, validation)
        @inbounds @atomicreplace status[1] Int32(0) => invalid
    end
end

_validation(owner, proper) =
    owner isa proper.CompiledPathsOwner ? owner.workspace.validation.device : nothing
_validation(owner::AOG._PixelIntegrationOwner, proper) = owner.validation.device

function FGA.prepare_capture(::FGA.CompleteTransition, device,
        prepared::FGA.PreparedAlgorithm{OpticalObservationF32},
        outputs::Tuple, inputs::Tuple, status, index::Integer)
    observation = prepared.plan.observation
    AOG._capture_safe_node_count(values(observation.nodes))
    check! = check_optical_status!(get_backend(status), 1)
    invalid = (Int32(index) << 8) | Int32(InvalidOpticalInput)
    proper = proper_extension()
    validation = map(child -> _validation(child.owner, proper), values(observation.nodes))
    product = AOG.graph_output(observation, :intensity)
    select = prepared.plan.single_scc_product && length(observation.nodes) == 1 &&
        keys(observation.outputs) == (:intensity,)
    propagate! = map(values(observation.nodes)) do child
        capture_observation(child.owner, product, proper, select, prepared.plan.merge_scc_fields)
    end
    return CapturedObservation(observation, inputs, only(outputs), status, check!, validation,
        propagate!, invalid)
end

function capture_observation(owner, product, proper, select, merge_fields)
    if select && owner isa proper.CompiledPathsOwner &&
            owner.propagation.workspace.pixel_mtf === nothing &&
            B.compute_device(product) isa B.AcceleratorComputeDevice{B.CUDABackend}
        if product === owner.outputs.fringed || product === owner.outputs.unfringed
            role = product === owner.outputs.fringed ? Val(:fringed) : Val(:unfringed)
            if merge_fields
                merge = proper._prepare_scc_field_merge(owner.propagation)
                merge === nothing || return Base.Fix2(proper._enqueue_captured_scc_product!, (role, merge))
            end
            return Base.Fix2(proper._enqueue_captured_scc_product!, role)
        end
    end
    return AOG.enqueue_captured_graph_node!
end

FGA.preflight(::CapturedObservation) = nothing
FGA.preserves_inputs(::CapturedObservation) = true
FGA.complete!(node::CapturedObservation, status::Int32) =
    status == node.invalid ? InvalidOpticalInput : nothing

function FGA.reset!(node::CapturedObservation)
    AOG._reset_nodes!(values(node.observation.nodes))
    fill!(node.intensity, 0)
    return nothing
end

function FGA.enqueue!(node::CapturedObservation)
    observation = node.observation
    B.copyto_backend_async!(AOG.graph_input(observation, :pupil_opd), node.inputs[1])
    B.copyto_backend_async!(AOG.graph_input(observation, :pupil_amplitude), node.inputs[2])
    for (child, validation, propagate!) in zip(values(observation.nodes), node.validation, node.propagate!)
        propagate!(child.owner)
        validation === nothing || node.check!(node.status, validation, node.invalid; ndrange=1)
    end
    B.copyto_backend_async!(node.intensity, AOG.graph_output(observation, :intensity))
    return nothing
end

# The pilot borrows exact prepared AOS buffers so DM and detector bindings stay
# in place. The owned-description API copies its configuration, so use the
# spec compiler with explicit bindings in this extension.
function prepare_observations(observations, target, single_scc_product, merge_scc_fields)
    specs = FG.NodeSpec[]
    inputs = FG.Endpoint[]
    outputs = FG.Endpoint[]
    for (name, observation) in pairs(observations)
        config = Dict{String,Union{AOG.PreparedAlgorithmGraph,Bool}}(
            "observation" => observation, "single_scc_product" => single_scc_product,
            "merge_scc_fields" => merge_scc_fields)
        push!(specs, FG.NodeSpec(name, FGA.algorithm_label(OpticalObservationF32),
            config, Dict{String,Nothing}()))
        push!(inputs, FG.Endpoint(name, :pupil_opd, Symbol(name, :_opd)),
            FG.Endpoint(name, :pupil_amplitude, Symbol(name, :_amplitude)))
        push!(outputs, FG.Endpoint(name, :intensity, name))
    end
    spec = FG.GraphSpec(specs, FG.LinkSpec[], inputs, outputs,
        FG.CpuExecutionSpec(0, nothing, :native, false))
    graph = FG._prepare_graph(spec; algorithms=(OpticalObservationF32,))
    return (; graph=prepare_execution(graph, target),
        destinations=map(graph -> AOG.graph_output(graph, :intensity), observations))
end

prepare_execution(graph, ::B.HostComputeDevice) = graph
function prepare_execution(graph, target::B.AcceleratorComputeDevice{B.CUDABackend})
    extension = Base.get_extension(AdaptiveOpticsSim, :AdaptiveOpticsSimCUDAExt)
    return FG.prepare_captured_graph(graph, extension._require_cuda_device(target))
end
function prepare_execution(graph, target::B.AcceleratorComputeDevice{B.AMDGPUBackend})
    extension = Base.get_extension(AdaptiveOpticsSim, :AdaptiveOpticsSimAMDGPUExt)
    return FG.prepare_captured_graph(graph, extension._require_amdgpu_device(target))
end
prepare_execution(graph, target) = throw(ArgumentError("native FilterGraphs capture is unavailable for $target"))

include("filtergraphs/optical_train.jl")
include("filtergraphs/pixel_integration.jl")
include("filtergraphs/llowfs.jl")
include("filtergraphs/pyramid_rate.jl")
include("filtergraphs/pwfs.jl")

mutable struct FilterGraphsPropagation{F,U,C,E,L}
    fringed::F
    unfringed::U
    crossing::C
    context::E
    llowfs::L
    failed::Bool
end

function O.prepare_filtergraphs_external_pupil(prepared::AOG.PreparedExternalPupil;
        single_scc_product::Bool=false, merge_scc_fields::Bool=false, llowfs_definition=nothing)
    merge_scc_fields && !single_scc_product && throw(ArgumentError("SCC field merging requires single_scc_product"))
    keys(prepared.graphs) == (:llowfs, :scc) ||
        throw(ArgumentError("FilterGraphs pilot requires SCC and LLOWFS observations"))
    prepared.state.sequence == 0 && prepared.propagation === nothing ||
        throw(ArgumentError("FilterGraphs pilot requires an unused independent optical owner"))
    scc, llowfs = prepared.graphs.scc, prepared.graphs.llowfs
    chopping = get(prepared.observations.scc, :chopping, nothing)
    chopping !== nothing && !haskey(prepared.observations.llowfs, :chopping) &&
        prepared.input_names == (;opd=:pupil_opd, amplitude=:pupil_amplitude) &&
        all(obs -> obs.output === :intensity, values(prepared.observations)) ||
        throw(ArgumentError("unexpected fitted SPIDERS optical bindings"))
    unfringed = chopping.unfringed
    proper = proper_extension()
    OpticalOwner = Union{proper.CompiledPathsOwner,proper.OpticalTrainOwner,
        proper.RuntimeOPDTrainOwner,AOG._PixelIntegrationOwner}
    for branch in (scc, unfringed, llowfs)
        isempty(branch.delayed_links) && eltype(AOG.graph_output(branch, :intensity)) === Float32 &&
            all(node -> node.owner isa OpticalOwner, values(branch.nodes)) ||
            throw(ArgumentError("FilterGraphs pilot requires stateless Float32 optical observations"))
    end
    target = B.compute_device(scc)
    direct = prepare_direct_llowfs(llowfs_definition, llowfs, target)
    wrapped = direct === nothing ? (; llowfs) : (;)
    propagation = FilterGraphsPropagation(
        prepare_observations((; wrapped..., scc), target, single_scc_product, merge_scc_fields),
        prepare_observations((; wrapped..., scc=unfringed), target, single_scc_product, merge_scc_fields),
        prepare_observations((; wrapped..., scc, scc_unfringed=unfringed), target, single_scc_product, merge_scc_fields),
        scc.context, direct, false,
    )
    return AOG.PreparedExternalPupil(prepared.graphs, prepared.exposures, prepared.pupil_opd,
        prepared.pupil_amplitude, prepared.observations, prepared.input_names, propagation,
        prepared.exposure_execution, prepared.state)
end

function AOG._external_pupil_execution_graph(propagation::FilterGraphsPropagation, graph, name, unfringed)
    name === :llowfs && propagation.llowfs !== nothing && return propagation.llowfs.graph
    return unfringed ? propagation.unfringed.graph : propagation.fringed.graph
end

function AOG._reset_external_pupil_propagation!(propagation::FilterGraphsPropagation)
    FG.reset!(propagation.fringed.graph)
    FG.reset!(propagation.unfringed.graph)
    FG.reset!(propagation.crossing.graph)
    propagation.llowfs === nothing || FG.reset!(propagation.llowfs.graph)
    propagation.failed = false
    return nothing
end

function run_observations!(graph::FG.PreparedGraph, inputs, destinations, context)
    return FG.process!(destinations, graph, inputs)
end

function run_observations!(graph::FG.PreparedCapturedGraph, inputs, destinations, context)
    return FG.process_device!(graph, inputs) do outputs
        B._with_prepared_device_execution_context(context) do
            try
                foreach(B.copyto_backend_async!, values(destinations), values(outputs))
            finally
                B._synchronize_prepared_device_execution_context!(context)
            end
        end
    end
end

function run_observations!(selection, prepared, context)
    input_names = keys(selection.graph.input_formats)
    opd, amplitude = prepared.pupil_opd, prepared.pupil_amplitude
    inputs = NamedTuple{input_names}(ntuple(i -> isodd(i) ? opd : amplitude, length(input_names)))
    result = run_observations!(selection.graph, inputs, selection.destinations, context)
    result isa OpticalFailure && throw(AOG.AlgorithmGraphError("invalid captured optical observation input"))
    return true
end

function AOG._propagate_external_pupil!(propagation::FilterGraphsPropagation, prepared, start, finish)
    propagation.failed && throw(AOG.AlgorithmGraphError(
        "external pupil propagation is stopped after a failed graph step; reset it before reuse"))
    scc = prepared.graphs.scc
    observation = prepared.observations.scc
    active = AOG._active_observation_graphs(scc, observation, start, finish,
        prepared.exposures.scc.duration_ns)
    try
        run_direct_llowfs!(propagation.llowfs, prepared)
        if length(active) == 2
            return run_observations!(propagation.crossing, prepared, propagation.context)
        elseif first(active) === scc
            return run_observations!(propagation.fringed, prepared, propagation.context)
        else
            return run_observations!(propagation.unfringed, prepared, propagation.context)
        end
    catch
        propagation.failed = true
        rethrow()
    end
end

end
