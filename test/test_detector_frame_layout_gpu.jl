# Optional selector: julia --project=<gpu-environment> test/test_detector_frame_layout_gpu.jl cuda|amdgpu
using Test
using LinearAlgebra
using AdaptiveOpticsSim
using AdaptiveOpticsSim.AlgorithmGraphs

length(ARGS) == 1 || error("select cuda or amdgpu")
if only(ARGS) == "cuda"
    using CUDA
    CUDA.allowscalar(false)
    const AdapterGPUArray = CUDA.CuArray
    adapter_gpu_synchronize() = CUDA.synchronize()
elseif only(ARGS) == "amdgpu"
    using AMDGPU
    AMDGPU.allowscalar(false)
    const AdapterGPUArray = AMDGPU.ROCArray
    adapter_gpu_synchronize() = AMDGPU.synchronize()
else
    error("select cuda or amdgpu")
end
BLAS.set_num_threads(1)

function gpu_layout_definition(
    input; range = nothing, rounding = "none", output_element_type = "float32"
)
    config = Dict{String,Any}(
        "input_size" => [2, 3],
        "output_size" => [5, 7],
        "origin" => [1, 2],
        "scale" => 2,
        "bias" => 0.5,
        "rounding" => rounding,
        "output_element_type" => output_element_type,
        "input_schema" => "test.detector-electrons.f32/1",
        "output_schema" => "test.detector-adc.f32/1",
    )
    isnothing(range) || (config["range"] = range)
    document = Dict{String,Any}(
        "schema_version" => 1,
        "name" => "detector_layout",
        "nodes" => [
            Dict(
                "name" => "layout",
                "type" => "detector_frame_layout_f32",
                "config" => config,
            ),
        ],
        "inputs" => [
            Dict(
                "name" => "detector_frame",
                "binding" => "input",
                "destination" => "layout.detector_frame",
            ),
        ],
        "outputs" => [Dict("name" => "frame", "source" => "layout.frame")],
    )
    return load_algorithm_graph(
        document;
        node_types = builtin_graph_node_types(),
        bindings = (; input),
    )
end

@testset "Accelerator detector layout stream and capture" begin
    first = reshape(Float32[0, 0.5, 1, -1, 8191.25, 10000], 2, 3)
    second = reverse(first; dims = 2)
    for options in (
        (;),
        (; range = [0, 16383]),
        (; range = [0, 16383], rounding = "nearest_even"),
        (; range = [0, 16383], rounding = "nearest_even", output_element_type = "uint16"),
    )
        input = AdapterGPUArray(first)
        target = AdaptiveOpticsSim.Backends.compute_device(input)
        definition = gpu_layout_definition(input; options...)
        stream = prepare_algorithm_graph(
            definition; target, execution = StreamGraphExecution()
        )
        captured = prepare_algorithm_graph(
            definition; target, execution = CapturedGraphExecution()
        )
        owner = AlgorithmGraphs.prepared_graph_node(captured, Val(:layout))
        @test typeof(AlgorithmGraphs.graph_node_capture_capability(owner)) ===
            AlgorithmGraphs.GraphNodeCaptureSafe
        @test captured_graph_node_count(captured) == 1
        @test all(isconcretetype, fieldtypes(typeof(owner)))
        addresses = (pointer(owner.detector_frame), pointer(owner.frame))
        first_result = nothing
        for source in (first, second)
            copyto!(input, source)
            fill!(graph_output(captured, Val(:frame)), 99.0f0)
            adapter_gpu_synchronize() # Publish caller input before retained execution contexts read it.
            step_graph!(stream)
            step_graph!(captured)
            T = get(options, :output_element_type, "float32") == "uint16" ? UInt16 : Float32
            expected = zeros(T, 5, 7)
            @test eltype(graph_output(captured, Val(:frame))) === T
            values = source .* 2.0f0 .+ 0.5f0
            haskey(options, :range) && (values = clamp.(values, 0.0f0, 16383.0f0))
            get(options, :rounding, "none") == "nearest_even" &&
                (values = round.(values, RoundNearest))
            expected[2:3, 3:5] .= values
            actual = Array(graph_output(captured, Val(:frame)))
            @test actual == expected
            @test actual == Array(graph_output(stream, Val(:frame)))
            @test all(iszero, actual[[1, 4, 5], :])
            isnothing(first_result) && (first_result = actual)
        end
        reset_graph!(captured)
        @test all(iszero, Array(graph_output(captured, Val(:frame))))
        @test (pointer(owner.detector_frame), pointer(owner.frame)) == addresses
        copyto!(input, first)
        adapter_gpu_synchronize()
        step_graph!(captured)
        @test Array(graph_output(captured, Val(:frame))) == first_result
        @test (pointer(owner.detector_frame), pointer(owner.frame)) == addresses
    end
end
