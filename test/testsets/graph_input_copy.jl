# CPU default; set AOS_TEST_GRAPH_INPUT_COPY_BACKEND=cuda in a CUDA environment.
using Test, AdaptiveOpticsSim
using AdaptiveOpticsSim.AlgorithmGraphs
const InputCopyGraphs = AdaptiveOpticsSim.AlgorithmGraphs
const InputCopyBackends = AdaptiveOpticsSim.Backends
const INPUT_COPY_BACKEND = get(ENV, "AOS_TEST_GRAPH_INPUT_COPY_BACKEND", "cpu")
INPUT_COPY_BACKEND in ("cpu", "cuda") || error("unsupported input-copy test backend")
input_copy_array(::Val{:cpu}, values) = Array(values)
input_copy_execution(::Val{:cpu}) = StreamGraphExecution()
if INPUT_COPY_BACKEND == "cuda"
    using CUDA
    CUDA.allowscalar(false)
    input_copy_array(::Val{:cuda}, values) = CuArray(values)
    input_copy_execution(::Val{:cuda}) = CapturedGraphExecution()
end
function input_copy_graph(::Type{T}, backend) where {T}
    command = input_copy_array(backend, T[1, 2])
    basis = ones(T, 2, 2, 2)
    basis[:, :, 2] .= 2
    definition = algorithm_graph((modal_opd_expansion_node(:modal;
        pupil_rows = 2, pupil_columns = 2, mode_count = 2, T,
        coefficients_schema = "test.hil.command/1", opd_schema = "test.hil.opd/1",
        basis_schema = "test.hil.basis/1", pupil_support_schema = "test.hil.support/1"),);
        inputs = (graph_input(:command, :modal => :coefficients, command),),
        outputs = (graph_output(:frame, :modal => :opd),),
        parameters = (sparse_parameter(:modal => :basis, input_copy_array(backend, basis)),
            sparse_parameter(:modal => :pupil_support, input_copy_array(backend, trues(2, 2)))))
    return prepare_algorithm_graph(definition; target = InputCopyBackends.compute_device(command),
        execution = input_copy_execution(backend))
end
struct FailedInputCopy <: AbstractVector{Float32} end
Base.size(::FailedInputCopy) = (2,)
Base.getindex(::FailedInputCopy, index::Int) = index == 1 ? 9f0 : error("copy failed")

@testset "prepared graph input copy: $(INPUT_COPY_BACKEND)" begin
    backend = Val(Symbol(INPUT_COPY_BACKEND))
    for T in (Float32, Float64)
        graph = input_copy_graph(T, backend)
        initial_input = graph_input(graph, :command)
        address = pointer(initial_input)
        @test InputCopyGraphs.copyto_graph_input!(graph, :command, T[3, -1]) === graph
        @test graph_input(graph, :command) === initial_input
        @test pointer(initial_input) == address
        @test graph_step_sequence(graph) == 0
        step_graph!(graph)
        @test Array(graph_output(graph, :frame)) == ones(T, 2, 2)
        @test_throws AlgorithmGraphError InputCopyGraphs.copyto_graph_input!(graph, :missing, T[3, -1])
        @test_throws AlgorithmGraphError InputCopyGraphs.copyto_graph_input!(graph, :command, zeros(T, 1, 2))
        @test_throws AlgorithmGraphError InputCopyGraphs.copyto_graph_input!(graph, :command, Int[3, -1])
        @test !graph_failed(graph)
        ticket = step_graph_async!(graph)
        @test_throws AlgorithmGraphError InputCopyGraphs.copyto_graph_input!(graph, :command, T[9, 9])
        wait_graph_step!(ticket)
        @test Array(graph_input(graph, :command)) == T[3, -1]
        for command in (T[-2, 4], T[8, 1])
            source = input_copy_array(backend, command)
            InputCopyGraphs.copyto_graph_input!(graph, :command, source)
            step_graph!(graph)
            @test Array(graph_output(graph, :frame)) == fill(command[1] + 2command[2], 2, 2)
            @test Array(source) == command
        end
        if INPUT_COPY_BACKEND == "cuda"
            @test captured_graph_node_count(graph) > 0
            caller_stream, caller_context = CUDA.stream(), CUDA.context()
            other_stream = CUDA.CuStream()
            CUDA.stream!(other_stream) do
                InputCopyGraphs.copyto_graph_input!(graph, :command, T[5, -2])
                @test CUDA.stream() === other_stream
                step_graph!(graph)
                @test CUDA.stream() === other_stream
            end
            @test CUDA.stream() === caller_stream
            @test CUDA.context() == caller_context
            @test Array(graph_output(graph, :frame)) == ones(T, 2, 2)
        end
    end
    if INPUT_COPY_BACKEND == "cpu"
        graph = input_copy_graph(Float32, backend)
        @test_throws ErrorException InputCopyGraphs.copyto_graph_input!(graph, :command, FailedInputCopy())
        @test graph_failed(graph)
        @test_throws AlgorithmGraphError InputCopyGraphs.copyto_graph_input!(graph, :command, Float32[1, 2])
        reset_graph!(graph)
        @test !graph_failed(graph)
        InputCopyGraphs.copyto_graph_input!(graph, :command, Float32[1, 2])
        step_graph!(graph)
        @test graph_output(graph, :frame) == fill(5f0, 2, 2)
    end
end
