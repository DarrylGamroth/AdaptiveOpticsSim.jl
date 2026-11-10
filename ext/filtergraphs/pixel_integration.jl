function pixel_integration_config(definition::AOG.AlgorithmNodeDefinition)
    AOG._node_type(definition) === AOG.PixelIntegrationNode{Float32} &&
        isempty(definition.props) ||
        throw(ArgumentError("expected a Float32 AOS pixel integration definition"))
    return definition.config::AOG.PixelIntegrationNodeConfig
end

function prepare_pixel_integration(definition)
    config = pixel_integration_config(definition)
    return AOG._PixelIntegrationPlan(config,
        AOG._pixel_integration_axis(Float32, config, 1),
        AOG._pixel_integration_axis(Float32, config, 2))
end

function process_pixels!(output, plan, input)
    Base.mightalias(input, output) &&
        throw(ArgumentError("pixel integration requires distinct input and output storage"))
    owner = AOG._PixelIntegrationOwner(plan, input, output, nothing)
    try
        AOG.step_graph_node!(owner)
    catch error
        error isa AOG.AlgorithmGraphError || rethrow()
        return InvalidOpticalInput
    end
    return nothing
end

@algorithm PixelIntegrationF32 begin
    label = "aos-pixel-integration-f32"
    failure = OpticalFailure

    @configuration PixelIntegrationConfiguration begin
        definition::AOG.AlgorithmNodeDefinition => graph_rebuild
    end

    ports(c) = (
        input(:input, Float32, pixel_integration_config(c.definition).input_shape,
            pixel_integration_config(c.definition).input_schema),
        output(:output, Float32, pixel_integration_config(c.definition).output_shape,
            pixel_integration_config(c.definition).output_schema; metadata_from=:input),
    )

    prepare(c) = prepare_pixel_integration(c.definition)
    process!(plan) = process_pixels!(output, plan, input)
end

FGA.capture_support(::Type{PixelIntegrationF32}, device) = optical_capture_support(device)

struct CapturedPixelIntegration{O,S,K}
    owner::O
    status::S
    check!::K
    invalid::Int32
end

function FGA.prepare_capture(::FGA.CompleteTransition, device,
        prepared::FGA.PreparedAlgorithm{PixelIntegrationF32},
        outputs::Tuple, inputs::Tuple, status, index::Integer)
    input, output = only(inputs), only(outputs)
    target = B.compute_device(input)
    owner = AOG.prepare_graph_node(AOG.PixelIntegrationNode{Float32}, prepared.plan.config,
        (;), (; input), (; output), (;), target)
    check! = check_optical_status!(get_backend(status), 1)
    invalid = (Int32(index) << 8) | Int32(InvalidOpticalInput)
    return CapturedPixelIntegration(owner, status, check!, invalid)
end

FGA.preflight(::CapturedPixelIntegration) = nothing
FGA.preserves_inputs(::CapturedPixelIntegration) = true
FGA.complete!(node::CapturedPixelIntegration, status::Int32) =
    status == node.invalid ? InvalidOpticalInput : nothing
FGA.reset!(node::CapturedPixelIntegration) = AOG.reset_graph_node!(node.owner)

function FGA.enqueue!(node::CapturedPixelIntegration)
    AOG.enqueue_captured_graph_node!(node.owner)
    node.check!(node.status, node.owner.validation.device, node.invalid; ndrange=1)
    return nothing
end
