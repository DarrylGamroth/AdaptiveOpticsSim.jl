function detector_frame_layout_config(definition::AOG.AlgorithmNodeDefinition)
    config = definition.config
    AOG._node_type(definition) === AOG.DetectorFrameLayoutNode &&
        config isa AOG.DetectorFrameLayoutConfig{UInt16} && isempty(definition.props) ||
        throw(ArgumentError("detector frame layout requires a Float32 input and UInt16 output definition with no properties"))
    return config
end

function bind_detector_frame_layout(config, detector_frame, frame)
    return AOG.prepare_graph_node(AOG.DetectorFrameLayoutNode, config,
        (;), (; detector_frame), (; frame), (;), B.compute_device(frame))
end

function process_detector_frame_layout!(frame, config, detector_frame)
    owner = bind_detector_frame_layout(config, detector_frame, frame)
    return AOG.step_graph_node!(owner)
end

@algorithm DetectorFrameLayoutF32 begin
    label = "aos-detector-frame-layout-f32"

    @configuration DetectorFrameLayoutConfiguration begin
        definition::AOG.AlgorithmNodeDefinition => graph_rebuild
    end

    ports(c) = (
        input(:detector_frame, Float32, detector_frame_layout_config(c.definition).input_size,
            detector_frame_layout_config(c.definition).input_schema),
        output(:frame, UInt16, detector_frame_layout_config(c.definition).output_size,
            detector_frame_layout_config(c.definition).output_schema; metadata_from=:detector_frame),
    )

    prepare(c) = detector_frame_layout_config(c.definition)
    process!(plan) = process_detector_frame_layout!(frame, plan, detector_frame)
end

struct CapturedDetectorFrameLayout{O}
    owner::O
end

FGA.capture_support(::Type{DetectorFrameLayoutF32}, device) = optical_capture_support(device)

function prepare_captured_detector_frame_layout(definition::AOG.AlgorithmNodeDefinition,
        detector_frame, frame)
    return prepare_captured_detector_frame_layout(detector_frame_layout_config(definition),
        detector_frame, frame)
end

function prepare_captured_detector_frame_layout(config::AOG.DetectorFrameLayoutConfig{UInt16},
        detector_frame, frame)
    owner = bind_detector_frame_layout(config, detector_frame, frame)
    AOG.graph_node_capture_capability(owner) isa AOG.GraphNodeCaptureSafe ||
        throw(ArgumentError("detector frame layout does not support native capture on the selected device"))
    return CapturedDetectorFrameLayout(owner)
end

function FGA.prepare_capture(::FGA.CompleteTransition, device,
        prepared::FGA.PreparedAlgorithm{DetectorFrameLayoutF32},
        outputs::Tuple, inputs::Tuple, status, index::Integer)
    return prepare_captured_detector_frame_layout(prepared.plan, only(inputs), only(outputs))
end

FGA.preflight(::CapturedDetectorFrameLayout) = nothing
FGA.preserves_inputs(::CapturedDetectorFrameLayout) = true
FGA.enqueue!(node::CapturedDetectorFrameLayout) = AOG.enqueue_captured_graph_node!(node.owner)
FGA.complete!(::CapturedDetectorFrameLayout, ::Int32) = nothing
FGA.reset!(node::CapturedDetectorFrameLayout) = AOG.reset_graph_node!(node.owner)
