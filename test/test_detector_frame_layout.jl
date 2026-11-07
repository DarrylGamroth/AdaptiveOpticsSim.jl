# Standalone CPU selector for the built-in detector frame layout factory.
using Test
using AdaptiveOpticsSim
using AdaptiveOpticsSim.AlgorithmGraphs

const AdapterAOG = AdaptiveOpticsSim.AlgorithmGraphs
const AdapterTarget = AdaptiveOpticsSim.Backends.HostComputeDevice()

layout_config() = Dict{String,Any}(
    "input_size" => [2, 3],
    "output_size" => [5, 7],
    "origin" => [1, 2],
    "scale" => 2,
    "bias" => -3,
    "input_schema" => "test.output-electrons.f32/1",
    "output_schema" => "test.input-electrons.f32/1",
)

function adapter_document(type, config; props = Dict())
    return Dict{String,Any}(
        "schema_version" => 1,
        "name" => "array_operations",
        "nodes" => [
            Dict(
                "name" => "operation",
                "type" => type,
                "config" => config,
                "props" => props,
            ),
        ],
        "inputs" => [
            Dict(
                "name" => "detector_frame",
                "binding" => "input",
                "destination" => "operation.detector_frame",
            ),
        ],
        "outputs" => [Dict("name" => "frame", "source" => "operation.frame")],
    )
end

function adapter_definition(type, config, input; props = Dict())
    return load_algorithm_graph(
        adapter_document(type, config; props);
        node_types = builtin_graph_node_types(),
        bindings = (; input),
    )
end
adapter_graph(type, config, input; kwargs...) = prepare_algorithm_graph(
    adapter_definition(type, config, input; kwargs...); target = AdapterTarget
)

@testset "UInt16 detector frame output" begin
    input = reshape(Float32[-2, 0.25, 0.75, 1.25, 8191.25, 10000], 2, 3)
    config = merge(
        layout_config(),
        Dict(
            "scale" => 2,
            "bias" => 0.5,
            "range" => [0, 16383],
            "rounding" => "nearest_even",
            "output_element_type" => "uint16",
        ),
    )
    graph = adapter_graph("detector_frame_layout_f32", config, input)
    frame = AdapterAOG.graph_output(graph, :frame)
    @test eltype(frame) === UInt16
    expected = zeros(UInt16, 5, 7)
    expected[2:3, 3:5] .= round.(UInt16, clamp.(2.0f0 .* input .+ 0.5f0, 0.0f0, 16383.0f0))
    step_graph!(graph)
    @test frame == expected
    fill!(frame, typemax(UInt16))
    step_graph!(graph)
    @test frame == expected
    reset_graph!(graph)
    @test all(iszero, frame)
    for changes in (
        Dict("rounding" => "none"),
        Dict("range" => [-1, 16383]),
        Dict("range" => [0, 65536]),
        Dict("output_element_type" => "int16"),
    )
        @test_throws AdapterAOG.AlgorithmGraphError adapter_definition(
            "detector_frame_layout_f32", merge(config, changes), input
        )
    end
end

adapter_owner(
    definition::AlgorithmGraphDefinition, inputs, outputs; target = AdapterTarget
) = adapter_owner(only(definition.nodes), inputs, outputs; target)
function adapter_owner(
    node::AlgorithmNodeDefinition{Name,Node}, inputs, outputs; target = AdapterTarget
) where {Name,Node}
    ports = AdapterAOG.graph_node_ports(Node, node.config)
    @test !isempty(ports)
    return AdapterAOG.prepare_graph_node(
        Node, node.config, node.props, inputs, outputs, NamedTuple(), target
    )
end

# Minimal offset-indexed host wrapper: rejection is the documented contract.
struct AdapterOffsetArray{T,N,A<:AbstractArray{T,N}} <: AbstractArray{T,N}
    parent::A
end
Base.size(array::AdapterOffsetArray) = size(array.parent)
Base.axes(array::AdapterOffsetArray) = map(n -> 0:(n-1), size(array))
Base.getindex(array::AdapterOffsetArray, indices::Vararg{Int}) =
    getindex(array.parent, (indices .+ 1)...)

function adapter_allocations(owner)
    AdapterAOG.step_graph_node!(owner)
    AdapterAOG.reset_graph_node!(owner)
    step = @allocated AdapterAOG.step_graph_node!(owner)
    reset = @allocated AdapterAOG.reset_graph_node!(owner)
    return (; step, reset)
end

@testset "Built-in detector frame layout factory" begin
    registry = builtin_graph_node_types()
    @test hasproperty(registry, :detector_frame_layout_f32)

    @testset "Rectangular placement, scale, bias, and optional clipping" begin
        input = reshape(Float32.(1:6), 2, 3)
        config = layout_config()
        graph = adapter_graph("detector_frame_layout_f32", config, input)
        frame = AdapterAOG.graph_output(graph, :frame)
        expected = zeros(Float32, 5, 7)
        expected[2:3, 3:5] .= input .* 2.0f0 .- 3.0f0
        step_graph!(graph)
        @test frame == expected
        @test frame[2, 3] == -1.0f0 # Omitted range preserves negative Float32 values.
        fill!(frame, 99.0f0)
        step_graph!(graph)
        @test frame == expected # Padding is cleared on every step.
        reset_graph!(graph)
        @test all(iszero, frame)
        @test input == reshape(Float32.(1:6), 2, 3)

        config["range"] = [0, 5]
        clipped = adapter_graph("detector_frame_layout_f32", config, input)
        step_graph!(clipped)
        expected[2:3, 3:5] .= clamp.(input .* 2.0f0 .- 3.0f0, 0.0f0, 5.0f0)
        @test AdapterAOG.graph_output(clipped, :frame) == expected

        config = layout_config()
        merge!(
            config,
            Dict(
                "input_size" => [240, 240],
                "output_size" => [242, 240],
                "origin" => [0, 0],
                "scale" => 1 / 5000,
                "bias" => 0,
            ),
        )
        electrons = fill(20_000.0f0, 240, 240)
        electrons[1, 1] = -5000.0f0
        electrons[2, 1] = 100_000_000.0f0
        gpi = adapter_graph("detector_frame_layout_f32", config, electrons)
        step_graph!(gpi)
        frame = AdapterAOG.graph_output(gpi, :frame)
        @test size(frame) == (242, 240)
        @test eltype(frame) === Float32
        @test frame[1:240, :] == electrons .* Float32(1 / 5000)
        @test frame[2, 1] > 16383.0f0
        @test all(iszero, frame[241:242, :])
    end

    @testset "Optional rounding completes Float32 camera samples" begin
        input = reshape(Float32[0, 0.5, 1, -1, 8191.25, 10000], 2, 3)
        original = copy(input)
        config = merge(
            layout_config(),
            Dict(
                "scale" => 2,
                "bias" => 0.5,
                "range" => [0, 16383],
                "rounding" => "nearest_even",
            ),
        )
        graph = adapter_graph("detector_frame_layout_f32", config, input)
        step_graph!(graph)
        frame = AdapterAOG.graph_output(graph, :frame)
        expected = zeros(Float32, 5, 7)
        expected[2:3, 3:5] .= Float32[0 2 16383; 2 0 16383]
        @test frame == expected
        @test eltype(frame) === Float32
        @test input == original
        fill!(frame, 99.0f0)
        step_graph!(graph)
        @test frame == expected
        reset_graph!(graph)
        @test all(iszero, frame)

        config["rounding"] = "none"
        unrounded = adapter_graph("detector_frame_layout_f32", config, input)
        step_graph!(unrounded)
        expected[2:3, 3:5] .= Float32[0.5 2.5 16383; 1.5 0 16383]
        @test AdapterAOG.graph_output(unrounded, :frame) == expected
        delete!(config, "rounding")
        default = adapter_graph("detector_frame_layout_f32", config, input)
        step_graph!(default)
        @test AdapterAOG.graph_output(default, :frame) == expected

        config["range"] = [0.25, 5.75]
        fractional = adapter_graph("detector_frame_layout_f32", config, input)
        step_graph!(fractional)
        expected[2:3, 3:5] .= Float32[0.5 2.5 5.75; 1.5 0.25 5.75]
        @test AdapterAOG.graph_output(fractional, :frame) == expected
    end

    @testset "Host views, one-based axes, concrete ownership, and allocation" begin
        definition = adapter_definition(
            "detector_frame_layout_f32", layout_config(), zeros(Float32, 2, 3)
        )
        input = view(reshape(Float32.(1:24), 4, 6), 1:2:4, 1:2:6)
        frame = view(fill(99.0f0, 10, 14), 1:2:10, 2:2:14)
        owner = adapter_owner(definition, (; detector_frame = input), (; frame))
        AdapterAOG.step_graph_node!(owner)
        expected = zeros(Float32, 5, 7)
        expected[2:3, 3:5] .= input .* 2.0f0 .- 3.0f0
        @test frame == expected
        @test all(isconcretetype, fieldtypes(typeof(owner)))
        @test adapter_allocations(owner) == (; step = 0, reset = 0)
        @test typeof(AdapterAOG.graph_node_capture_capability(owner)) ===
            AdapterAOG.GraphNodeCaptureUnsupported
        clipped_config = layout_config()
        clipped_config["range"] = [0, 2]
        clipped_definition = adapter_definition(
            "detector_frame_layout_f32", clipped_config, zeros(Float32, 2, 3)
        )
        clipped_owner = adapter_owner(
            clipped_definition, (; detector_frame = input), (; frame)
        )
        @test adapter_allocations(clipped_owner) == (; step = 0, reset = 0)
        rounded_config = merge(
            layout_config(),
            Dict("bias" => 0.5, "range" => [0, 16383], "rounding" => "nearest_even"),
        )
        rounded_definition = adapter_definition(
            "detector_frame_layout_f32", rounded_config, zeros(Float32, 2, 3)
        )
        rounded_owner = adapter_owner(
            rounded_definition, (; detector_frame = input), (; frame)
        )
        AdapterAOG.step_graph_node!(rounded_owner)
        expected[2:3, 3:5] .= input .* 2.0f0
        @test frame == expected
        @test all(isconcretetype, fieldtypes(typeof(rounded_owner)))
        @test adapter_allocations(rounded_owner) == (; step = 0, reset = 0)
        @test_throws ArgumentError adapter_owner(
            definition,
            (; detector_frame = AdapterOffsetArray(zeros(Float32, 2, 3))),
            (; frame),
        )
        @test_throws AlgorithmGraphError adapter_owner(
            definition,
            (; detector_frame = zeros(Float32, 2, 3)),
            (; frame = zeros(Float64, 5, 7)),
        )
        @test_throws AlgorithmGraphError adapter_owner(
            definition,
            (; detector_frame = zeros(Float32, 2, 3)),
            (; frame = zeros(Float32, 5, 7));
            target = AdaptiveOpticsSim.Backends.AcceleratorComputeDevice(
                AdaptiveOpticsSim.Backends.CUDABackend(), 0
            ),
        )
    end

    @testset "Malformed configurations fail through the public loader" begin
        layout_changes = (
            config -> (config["extra"] = 1),
            config -> delete!(config, "scale"),
            config -> (config["input_size"] = [2]),
            config -> (config["output_size"] = [typemax(Int), 2]),
            config -> (config["output_size"] = [0, 7]),
            config -> (config["origin"] = [1]),
            config -> (config["origin"] = [-1, 0]),
            config -> (config["origin"] = [typemax(Int), 0]),
            config -> (config["origin"] = [4, 0]),
            config -> (config["origin"] = [0, 5]),
            config -> (config["origin"] = [0.5, 0]),
            config -> (config["scale"] = Inf),
            config -> (config["scale"] = NaN),
            config -> (config["scale"] = 1e300),
            config -> (config["scale"] = true),
            config -> (config["bias"] = -Inf),
            config -> (config["bias"] = -1e300),
            config -> (config["range"] = [2, 1]),
            config -> (config["range"] = [1 + eps(), 1]),
            config -> (config["range"] = [0, Inf]),
            config -> (config["range"] = [0]),
            config -> (config["range"] = [false, true]),
            config -> (config["range"] = nothing),
            config -> (config["rounding"] = "nearest"),
            config -> (config["rounding"] = true),
            config -> (config["rounding"] = nothing),
            config ->
                merge!(config, Dict("rounding" => "nearest_even", "range" => [0.5, 10])),
            config -> merge!(
                config, Dict("rounding" => "nearest_even", "range" => [1 + eps(), 10])
            ),
            config ->
                merge!(config, Dict("rounding" => "nearest_even", "range" => [0, 10.5])),
        )
        for change in layout_changes
            config = layout_config()
            change(config)
            @test_throws AlgorithmGraphError adapter_definition(
                "detector_frame_layout_f32", config, zeros(Float32, 2, 3)
            )
        end
        @test_throws AlgorithmGraphError adapter_definition(
            "detector_frame_layout_f32",
            layout_config(),
            zeros(Float32, 2, 3);
            props = Dict("extra" => 1),
        )
        @test_throws AlgorithmGraphError adapter_graph(
            "detector_frame_layout_f32", layout_config(), zeros(Float64, 2, 3)
        )
    end
end
