# Run explicitly in a project containing CUDA and this checkout:
# JULIA_NUM_THREADS=1 julia --startup-file=no --project=<environment> test/cuda_ingaas_graph_capture.jl
using Test
using AdaptiveOpticsSim
using AdaptiveOpticsSim.Optics
using AdaptiveOpticsSim.Backends
using AdaptiveOpticsSim.Detectors
using AdaptiveOpticsSim.AlgorithmGraphs
using CUDA
using LinearAlgebra

CUDA.allowscalar(false)
BLAS.set_num_threads(1)
Backends.set_fft_provider_threads!(1)

function cuda_ingaas_definition(rate; kwargs...)
    config = merge((rows=size(rate, 1), columns=size(rate, 2),
            pixel_sampling_m=(12e-6, 14e-6), spectral_channel=:capture_test_band,
            exposure_duration_s=0.125, quantum_efficiency=0.5, gain=2,
            rng_seed=0x1984, photon_rate_schema="test.ingaas.cuda-photon-rate.f32/1",
            frame_schema="test.ingaas.cuda-frame.f32/1", photon_noise=false),
        (; kwargs...))
    return algorithm_graph((ingaas_detector_acquisition_node(:detector; config...),);
        name=:cuda_ingaas_capture,
        inputs=(graph_input(:rate, :detector => :photon_rate, rate),),
        outputs=(graph_output(:frame, :detector => :frame),))
end

function compare_cuda_ingaas_capture(rate, target; stochastic=false, kwargs...)
    definition = cuda_ingaas_definition(rate; kwargs...)
    stream = prepare_algorithm_graph(definition; target, execution=StreamGraphExecution())
    captured = prepare_algorithm_graph(definition; target, execution=CapturedGraphExecution())
    @test captured_graph_node_count(captured) == 1
    owner = AlgorithmGraphs.prepared_graph_node(captured, Val(:detector))
    @test owner.prepared.detector.params.sensor isa InGaAsSensor
    @test AlgorithmGraphs.graph_node_capture_capability(owner) isa
        AlgorithmGraphs.GraphNodeCaptureSafe
    initial_rate = Array(rate)
    first_frame = nothing
    first_draw = UInt64(0)
    for frame_index in 1:2
        if frame_index == 2
            copyto!(rate, 1.75f0 .* initial_rate)
            CUDA.synchronize() # Caller publishes input before retained graph streams read it.
        end
        step_graph!(stream)
        step_graph!(captured)
        expected = Array(graph_output(stream, Val(:frame)))
        actual = Array(graph_output(captured, Val(:frame)))
        @test actual ≈ expected rtol=2e-5 atol=1e-6
        @test all(isfinite, actual)
        @test size(actual) == size(initial_rate)
        draw = only(Array(owner.rng.state.draw_sequence))
        if frame_index == 1
            first_frame = actual
            first_draw = draw
        else
            @test actual != first_frame
            if stochastic
                @test draw > first_draw > 0
            else
                @test draw == first_draw == 0
            end
        end
    end
    copyto!(rate, initial_rate)
    CUDA.synchronize()
    reset_graph!(stream)
    reset_graph!(captured)
    @test all(iszero, Array(graph_output(captured, Val(:frame))))
    @test only(Array(owner.rng.state.draw_sequence)) == 0
    step_graph!(stream)
    step_graph!(captured)
    @test Array(graph_output(captured, Val(:frame))) == first_frame
    @test Array(graph_output(captured, Val(:frame))) ≈
        Array(graph_output(stream, Val(:frame))) rtol=2e-5 atol=1e-6
    return nothing
end

@testset "CUDA InGaAs acquisition graph capture" begin
    rate = CuArray(reshape(Float32.(10:57), 8, 6))
    target = compute_device(rate)
    compare_cuda_ingaas_capture(rate, target) # Null glow and deterministic acquisition.
    compare_cuda_ingaas_capture(rate, target;
        stochastic=true, photon_noise=true, readout_noise=true, readout_noise_e=0.7)
    compare_cuda_ingaas_capture(rate, target;
        stochastic=true, photon_noise=true, readout_noise=true, readout_noise_e=0.7,
        dark_current_e_per_pixel_s=2, glow_e_per_pixel_s=3,
        bits=10, full_well_e=100)
    for unsupported in ((persistence_model=ExponentialPersistence(0.1, 0.5),),
            (nonlinearity_model=SaturatingFrameNonlinearity(0.01),),
            (bits=32, full_well_e=100))
        definition = cuda_ingaas_definition(rate; unsupported...)
        @test_throws AlgorithmGraphError prepare_algorithm_graph(definition;
            target, execution=CapturedGraphExecution())
    end
end
