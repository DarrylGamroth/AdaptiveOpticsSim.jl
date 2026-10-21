using Test
using CUDA
using AdaptiveOpticsSim
using AdaptiveOpticsSim.AlgorithmGraphs
const HILBackends = AdaptiveOpticsSim.Backends
const HILGraphs = AdaptiveOpticsSim.AlgorithmGraphs

function _hil_copy_graph(::Type{T}; masked=false) where {T}
    command = CuArray(T[1, 2])
    basis = ones(T, 2, 2, 2)
    basis[:, :, 2] .= 2
    support = trues(2, 2)
    masked && (support[1, 2] = false)
    definition = algorithm_graph((modal_opd_expansion_node(:modal;
        pupil_rows=2, pupil_columns=2, mode_count=2, T,
        coefficients_schema="test.hil.command/1", opd_schema="test.hil.opd/1",
        basis_schema="test.hil.basis/1", pupil_support_schema="test.hil.support/1"),);
        inputs=(graph_input(:command, :modal => :coefficients, command),),
        outputs=(graph_output(:frame, :modal => :opd),),
        parameters=(sparse_parameter(:modal => :basis, CuArray(basis)),
            sparse_parameter(:modal => :pupil_support, CuArray(support))))
    return prepare_algorithm_graph(definition; target=HILBackends.compute_device(command),
        execution=CapturedGraphExecution())
end
function _hil_copy_cycle!(boundary)
    sequence = step_hil_frame!(boundary)
    adopt_hil_command!(boundary, sequence)
    return nothing
end
_hil_copy_bytes(boundary) = @allocated _hil_copy_cycle!(boundary)
_exposure_copy_bytes(boundary) = @allocated step_hil_exposure!(boundary)

function _hil_pending_kernel!(output)
    value = zero(eltype(output))
    for _ in 1:2_000_000
        value = muladd(value, 0.99999f0, 0.00001f0)
    end
    output[CUDA.threadIdx().x] = value
    return nothing
end

@testset "CUDA HIL completed host copies" begin
    CUDA.allowscalar(false)
    caller_stream, caller_context = CUDA.stream(), CUDA.context()
    for T in (Float32, Float64)
        graph = _hil_copy_graph(T)
        boundary = prepare_graph_hil_boundary(graph; command_input=:command, frame_output=:frame)
        command, frame = hil_command_buffer(boundary), hil_frame_buffer(boundary)
        addresses = (pointer(command), pointer(frame), pointer(boundary.command_input), pointer(boundary.frame_output))
        @test command == T[1, 2]
        sequence = step_hil_frame!(boundary)
        @test frame == fill(T(5), 2, 2)
        command .= T[3, -1]
        adopt_hil_command!(boundary, sequence)
        sequence = step_hil_frame!(boundary)
        @test frame == ones(T, 2, 2)
        command .= T(NaN)
        @test_throws AlgorithmGraphError adopt_hil_command!(boundary, sequence)
        @test Array(boundary.command_input) == T[3, -1]
        command .= T[3, -1]
        adopt_hil_command!(boundary, sequence)
        for _ in 1:32
            _hil_copy_cycle!(boundary)
        end
        @test _hil_copy_bytes(boundary) == 0
        @test CUDA.stream() == caller_stream
        @test CUDA.context() == caller_context
        reset_hil_boundary!(boundary)
        @test command == T[1, 2]
        @test all(iszero, frame)
        @test addresses == (pointer(command), pointer(frame), pointer(boundary.command_input), pointer(boundary.frame_output))
        step_hil_frame!(boundary)
        @test frame == fill(T(5), 2, 2)

        calibration = prepare_graph_calibration_boundary(_hil_copy_graph(T);
            command_input=:command, frame_output=:frame)
        hil_command_buffer(calibration) .= T[-2, 4]
        adopt_hil_probe!(calibration, UInt64(7))
        for _ in 1:32
            step_hil_exposure!(calibration)
        end
        @test hil_frame_buffer(calibration) == fill(T(6), 2, 2)
        @test _exposure_copy_bytes(calibration) == 0
        @test hil_boundary_status(calibration).active_probe_sequence == UInt64(7)

        # Pending producer work must remain ordered when pointer conversion
        # transfers managed CUDA storage to the prepared graph's stream.
        producer = CUDA.CuStream()
        pending = CUDA.zeros(T, 2)
        CUDA.@cuda threads=2 _hil_pending_kernel!(pending)
        CUDA.synchronize()
        CUDA.stream!(producer) do
            CUDA.@cuda threads=2 _hil_pending_kernel!(pending)
            fill!(pending, T(9))
        end
        destination = zeros(T, 2)
        HILGraphs._copy_hil_buffer!(graph, destination, pending)
        @test destination == fill(T(9), 2)
        @test CUDA.stream() == caller_stream
        @test CUDA.context() == caller_context
        @test_throws AdaptiveOpticsSim.DimensionMismatchError HILGraphs._copy_hil_buffer!(
            graph, zeros(T, 1, 4), CUDA.zeros(T, 2, 2))
        @test CUDA.stream() == caller_stream
        @test CUDA.context() == caller_context
        empty_host = zeros(T, 0, 2)
        @test HILGraphs._copy_hil_buffer!(graph, empty_host, CUDA.zeros(T, 0, 2)) === empty_host
        # Storage kinds outside DeviceMemory retain the ordinary copy path.
        host_storage = CUDA.CuArray{T,2,CUDA.HostMemory}(undef, 2, 2)
        copyto!(host_storage, fill(T(11), 2, 2))
        @test HILGraphs._copy_hil_buffer!(graph, zeros(T, 2, 2), host_storage) == fill(T(11), 2, 2)
    end
    @test CUDA.stream() == caller_stream
    @test CUDA.context() == caller_context
end

@testset "CUDA captured modal OPD numerical edges" begin
    for T in (Float32, Float64)
        graph = _hil_copy_graph(T; masked=true)
        command = graph_input(graph, Val(:command))
        for coefficients in (T[2, -1], T[-3, 0.5], T[NaN, 1], T[Inf, -Inf])
            copyto!(command, coefficients)
            CUDA.synchronize()
            step_graph!(graph)
            expected = coefficients[1] + 2 * coefficients[2]
            actual = Array(graph_output(graph, Val(:frame)))
            @test actual[1, 2] === zero(T)
            @test all(value -> isequal(value, expected), actual[[1, 2, 4]])
        end
    end
end
