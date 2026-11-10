module AdaptiveOpticsSimFilterGraphAlgorithmsExt

using AdaptiveOpticsSim
using FilterGraphAlgorithms: @algorithm
import FilterGraphAlgorithms
const FGA = FilterGraphAlgorithms
import KernelAbstractions
using KernelAbstractions: @kernel, @index, @atomicreplace

const B = AdaptiveOpticsSim.Backends
const O = AdaptiveOpticsSim.Optics
const D = AdaptiveOpticsSim.Detectors
const AOG = AdaptiveOpticsSim.AlgorithmGraphs
const IMAGE_SCHEMA = "org.adaptiveopticssim.optical-intensity/1"
const PRODUCTS = (:fringed, :unfringed, :reference)

@enum SCCOpticalFailure::UInt8 SCCNonFiniteInput=1 SCCNonFiniteOutput=2

struct OpticalModel{M,R}
    metadata::M
    read_array::R
end

struct PixelResponse{T}
    transfer::T
    padding::Int
end

struct OpticalWorkspace{P,F}
    propagation::P
    propagate!::F
end

struct PixelWorkspace{P}
    response::P
end

const SCCWorkspace = Union{OpticalWorkspace,PixelWorkspace}
FGA.commit!(::SCCWorkspace) = nothing
FGA.discard!(::SCCWorkspace) = nothing
FGA.reset!(::SCCWorkspace) = nothing

function prepare_optics(plan::OpticalModel, target)
    paths = Base.get_extension(AdaptiveOpticsSim, :AdaptiveOpticsSimProperExt)
    owner = paths._prepare_compiled_paths(
        plan.metadata, plan.read_array, nothing; T=Float32, target, summation=:coherent_modes)
    return OpticalWorkspace(owner, paths._enqueue_compiled_paths!)
end

function prepare_optics(configuration)
    metadata, read_array = AOG._read_compiled_optical_model(
        configuration.model_path, configuration.model_sha256)
    Int(metadata["image_size"]) == configuration.image_size &&
        Int(metadata["live"]["dm_size"]) == configuration.dm_size &&
        Int(metadata["live"]["entrance_size"]) == configuration.pupil_size ||
        throw(DimensionMismatch("SCC model and configured grids differ"))
    !metadata["fast"]["apply_pixel_mtf"] ||
        throw(ArgumentError("SCC model must leave pixel MTF to the detector node"))
    plan = OpticalModel(metadata, read_array)
    return (; plan, workspace=prepare_optics(plan, B.HostComputeDevice()))
end

function prepare_response(configuration)
    metadata, read_array = AOG._read_compiled_optical_model(
        configuration.model_path, configuration.model_sha256)
    Int(metadata["image_size"]) == configuration.image_size ||
        throw(DimensionMismatch("SCC model and detector image grids differ"))
    fast = metadata["fast"]
    centered = read_array(fast["pixel_mtf"])
    all(iszero, imag.(centered)) || throw(ArgumentError("pixel MTF transfer must be real"))
    plan = PixelResponse(circshift(real.(centered), map(n -> -fld(n, 2), size(centered))), Int(fast["pixel_pad"]))
    image = zeros(Float32, configuration.image_size, configuration.image_size)
    workspace = PixelWorkspace(D.prepare_pixel_mtf(image, plan.transfer; padding=plan.padding))
    return (; plan, workspace)
end

function run_optics!(outputs, workspace::OpticalWorkspace)
    products = workspace.propagate!(workspace.propagation)
    for (destination, name) in zip(outputs, PRODUCTS)
        B.copyto_backend_async!(destination, getproperty(products, name))
    end
    return nothing
end

function process_optics!(outputs, workspace, inputs)
    all(array -> all(isfinite, array), inputs) || return SCCNonFiniteInput
    foreach(copyto!, values(workspace.propagation.inputs), inputs)
    run_optics!(outputs, workspace)
    return all(array -> all(isfinite, array), outputs) ? nothing : SCCNonFiniteOutput
end

function process_response!(detector, workspace, intensity)
    all(isfinite, intensity) || return SCCNonFiniteInput
    response = workspace.response
    copyto!(response.input, intensity)
    filtered = D.apply_pixel_mtf!(response)
    all(isfinite, filtered) || return SCCNonFiniteOutput
    copyto!(detector, filtered)
    return nothing
end

@algorithm FittedSCCPropagationF32 begin
    label = "aos-fitted-scc-propagation-f32"
    failure = SCCOpticalFailure

    @configuration FittedSCCPropagationConfiguration begin
        model_path::String => graph_rebuild
        model_sha256::String => graph_rebuild
        dm_size::Int => graph_rebuild
        pupil_size::Int => graph_rebuild
        image_size::Int => graph_rebuild
    end

    ports(c) = (
        input(:dm_opd, Float32, (c.dm_size, c.dm_size), "org.adaptiveopticssim.dm-opd/1"),
        input(:pupil_opd, Float32, (c.pupil_size, c.pupil_size), "org.adaptiveopticssim.pupil-opd/1"),
        input(:attenuation, Float32, (c.pupil_size, c.pupil_size), "org.adaptiveopticssim.entrance-attenuation/1"),
        output(:fringed, Float32, (c.image_size, c.image_size), IMAGE_SCHEMA; metadata_from=:pupil_opd),
        output(:unfringed, Float32, (c.image_size, c.image_size), IMAGE_SCHEMA; metadata_from=:pupil_opd),
        output(:reference, Float32, (c.image_size, c.image_size), IMAGE_SCHEMA; metadata_from=:pupil_opd),
    )

    prepare(c) = prepare_optics(c)
    process!(plan, workspace) = process_optics!(
        (fringed, unfringed, reference), workspace, (dm_opd, pupil_opd, attenuation))
end

@algorithm FittedSCCPixelMTFF32 begin
    label = "aos-fitted-scc-pixel-mtf-f32"
    failure = SCCOpticalFailure

    @configuration FittedSCCPixelMTFConfiguration begin
        model_path::String => graph_rebuild
        model_sha256::String => graph_rebuild
        image_size::Int => graph_rebuild
    end

    ports(c) = (
        input(:intensity, Float32, (c.image_size, c.image_size), IMAGE_SCHEMA),
        output(:detector, Float32, (c.image_size, c.image_size), IMAGE_SCHEMA; metadata_from=:intensity),
    )

    prepare(c) = prepare_response(c)
    process!(plan, workspace) = process_response!(detector, workspace, intensity)
end

const SCCDeclarations = Union{FittedSCCPropagationF32,FittedSCCPixelMTFF32}

function FGA.capture_support(::Type{<:SCCDeclarations}, device)
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

struct CapturedSCCNode{W,I,O,S,K}
    workspace::W
    inputs::I
    outputs::O
    status::S
    validate!::K
    invalid_input::Int32
    invalid_output::Int32
end

@kernel function validate_finite!(status, values, failure)
    i = @index(Global, Linear)
    if i <= length(values) && @inbounds(!isfinite(values[i]))
        @inbounds @atomicreplace status[1] Int32(0) => failure
    end
end

@kernel function stage_finite!(destination, source, status, failure)
    i = @index(Global, Linear)
    if i <= length(destination)
        @inbounds value = source[i]
        if !isfinite(value)
            @inbounds @atomicreplace status[1] Int32(0) => failure
            value = zero(value)
        end
        @inbounds destination[i] = value
    end
end

function prepare_captured_scc(workspace, inputs, outputs, status, index)
    0 < index <= (typemax(Int32) >> 8) || throw(ArgumentError("SCC node index cannot encode status"))
    code = Int32(index) << 8
    validate! = validate_finite!(KernelAbstractions.get_backend(status), 256)
    return CapturedSCCNode(workspace, inputs, outputs, status, validate!,
        code | Int32(SCCNonFiniteInput), code | Int32(SCCNonFiniteOutput))
end

function FGA.prepare_capture(::FGA.CompleteTransition, device,
        prepared::FGA.PreparedAlgorithm{FittedSCCPropagationF32},
        outputs::Tuple, inputs::Tuple, status, index::Integer)
    target = B.compute_device(first(inputs))
    workspace = prepare_optics(prepared.plan, target)
    return prepare_captured_scc(workspace, inputs, outputs, status, index)
end

function FGA.prepare_capture(::FGA.CompleteTransition, device,
        prepared::FGA.PreparedAlgorithm{FittedSCCPixelMTFF32},
        outputs::Tuple, inputs::Tuple, status, index::Integer)
    plan = prepared.plan
    response = D.prepare_pixel_mtf(only(inputs), plan.transfer;
        padding=plan.padding, output=only(outputs))
    return prepare_captured_scc(PixelWorkspace(response), inputs, outputs, status, index)
end

FGA.preflight(::CapturedSCCNode) = nothing
FGA.preserves_inputs(::CapturedSCCNode) = true

function FGA.reset!(node::CapturedSCCNode)
    foreach(array -> fill!(array, 0), node.outputs)
    return nothing
end

function FGA.complete!(node::CapturedSCCNode, status::Int32)
    status == node.invalid_input && return SCCNonFiniteInput
    status == node.invalid_output && return SCCNonFiniteOutput
    return nothing
end

function FGA.enqueue!(node::CapturedSCCNode{<:OpticalWorkspace})
    stage! = stage_finite!(KernelAbstractions.get_backend(node.status), 256)
    for (destination, source) in zip(values(node.workspace.propagation.inputs), node.inputs)
        stage!(destination, source, node.status, node.invalid_input; ndrange=length(source))
    end
    run_optics!(node.outputs, node.workspace)
    for output in node.outputs
        node.validate!(node.status, output, node.invalid_output; ndrange=length(output))
    end
    return nothing
end

function FGA.enqueue!(node::CapturedSCCNode{<:PixelWorkspace})
    input, output = only(node.inputs), only(node.outputs)
    node.validate!(node.status, input, node.invalid_input; ndrange=length(input))
    D.apply_pixel_mtf!(node.workspace.response)
    node.validate!(node.status, output, node.invalid_output; ndrange=length(output))
    return nothing
end

function O.fitted_scc_filter_graph(model_path::AbstractString, model_sha256::AbstractString)
    path = abspath(model_path)
    metadata, _ = AOG._read_compiled_optical_model(path, model_sha256)
    detector = Dict("model_path" => path, "model_sha256" => String(model_sha256),
        "image_size" => Int(metadata["image_size"]))
    optics = merge(detector, Dict("dm_size" => Int(metadata["live"]["dm_size"]),
        "pupil_size" => Int(metadata["live"]["entrance_size"])))
    nodes = [
        Dict("type" => "ndarray", "name" => "optics",
            "label" => FGA.algorithm_label(FittedSCCPropagationF32), "config" => optics),
        Dict("type" => "ndarray", "name" => "pixel",
            "label" => FGA.algorithm_label(FittedSCCPixelMTFF32), "config" => detector),
    ]
    inputs = [Dict("name" => name, "endpoint" => "optics:$name")
        for name in ("dm_opd", "pupil_opd", "attenuation")]
    outputs = [Dict("name" => String(name), "endpoint" => "optics:$name") for name in PRODUCTS]
    push!(outputs, Dict("name" => "detector", "endpoint" => "pixel:detector"))
    links = [Dict("output" => "optics:fringed", "input" => "pixel:intensity")]
    return (; document=Dict("nodes" => nodes, "links" => links, "inputs" => inputs, "outputs" => outputs),
        algorithms=(FittedSCCPropagationF32, FittedSCCPixelMTFF32))
end

end
