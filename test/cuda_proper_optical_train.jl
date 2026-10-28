# Run explicitly in an environment containing this checkout, Proper and CUDA:
# JULIA_NUM_THREADS=1 julia --startup-file=no --project=<environment> test/cuda_proper_optical_train.jl
using Test
using AdaptiveOpticsSim
using AdaptiveOpticsSim.AlgorithmGraphs
using AdaptiveOpticsSim.Backends
using Proper
using CUDA
using LinearAlgebra

CUDA.allowscalar(false)
BLAS.set_num_threads(1)
Backends.set_fft_provider_threads!(1)

function cuda_proper_definition(config, opd, amplitude; bindings=(;))
    T = eltype(opd)
    node = proper_optical_train_node_factory(; T)(:optics, config, (;), bindings)
    return algorithm_graph((node,); name=:cuda_proper_train,
        inputs=(graph_input(:opd, :optics => :pupil_opd, opd),
            graph_input(:amplitude, :optics => :pupil_amplitude, amplitude)),
        outputs=(graph_output(:intensity, :optics => :output),))
end

function compare_cuda_proper_capture(::Type{T}, config; bindings=(;)) where {T}
    p = config.pupil_resolution
    initial_opd = T[3e-8 * sin(0.7i) * cos(0.4j) for i in 1:p, j in 1:p]
    initial_amplitude = T[0.5 + 0.2cos(0.3i + 0.6j) for i in 1:p, j in 1:p]
    opd, amplitude = CuArray(initial_opd), CuArray(initial_amplitude)
    target = compute_device(opd)
    definition = cuda_proper_definition(config, opd, amplitude; bindings)
    stream = prepare_algorithm_graph(definition; target, execution=StreamGraphExecution())
    captured = prepare_algorithm_graph(definition; target, execution=CapturedGraphExecution())
    host_opd, host_amplitude = copy(initial_opd), copy(initial_amplitude)
    host = prepare_algorithm_graph(cuda_proper_definition(config, host_opd, host_amplitude; bindings);
        target=Backends.HostComputeDevice())
    @test captured_graph_node_count(captured) == 1
    owner = AlgorithmGraphs.prepared_graph_node(captured, Val(:optics))
    @test AlgorithmGraphs.graph_node_capture_capability(owner) isa AlgorithmGraphs.GraphNodeCaptureSafe
    field = owner.workspace.wavefront.field
    fft_workspace = owner.workspace.context.workspace.fft
    fft_scratch = fft_workspace.scratch
    fft_plans = (fft_workspace.forward_plan, fft_workspace.backward_plan)
    masks = owner.workspace.masks
    rtol, atol = T === Float32 ? (5e-4, 5e-6) : (5e-11, 5e-12)
    first_frame = zeros(T, config.resolution, config.resolution)
    first_sampling = optical_sampling_m(owner)
    for frame_index in 1:4
        host_opd .= initial_opd .+ T(frame_index - 1) .* T[4e-8 * (i - 0.4j) / p for i in 1:p, j in 1:p]
        host_amplitude .= initial_amplitude .* T(1 - 0.13 * (frame_index - 1))
        copyto!(opd, host_opd)
        copyto!(amplitude, host_amplitude)
        CUDA.synchronize() # Publish caller input before the retained graph streams consume it.
        step_graph!(host)
        step_graph!(stream)
        step_graph!(captured)
        actual = Array(graph_output(captured, Val(:intensity)))
        expected = Array(graph_output(stream, Val(:intensity)))
        @test actual ≈ expected rtol=rtol atol=atol
        @test actual ≈ graph_output(host, Val(:intensity)) rtol=rtol atol=atol
        @test all(isfinite, actual)
        @test sum(actual) > 0
        @test optical_sampling_m(owner) ≈ optical_sampling_m(AlgorithmGraphs.prepared_graph_node(stream, Val(:optics)))
        @test owner.workspace.wavefront.field === field
        @test fft_workspace.scratch === fft_scratch
        @test fft_workspace.forward_plan === fft_plans[1]
        @test fft_workspace.backward_plan === fft_plans[2]
        @test owner.workspace.masks === masks
        if frame_index == 1
            first_frame .= actual
            first_sampling = optical_sampling_m(owner)
        else
            @test !isapprox(actual, first_frame; rtol=1e-3, atol=atol)
        end
    end
    copyto!(opd, initial_opd)
    copyto!(amplitude, initial_amplitude)
    host_opd .= initial_opd
    host_amplitude .= initial_amplitude
    CUDA.synchronize()
    reset_graph!(stream)
    reset_graph!(captured)
    @test all(iszero, Array(graph_output(captured, Val(:intensity))))
    @test owner.workspace.wavefront.field === field
    step_graph!(stream)
    step_graph!(captured)
    @test Array(graph_output(captured, Val(:intensity))) ≈ first_frame rtol=rtol atol=atol
    @test Array(graph_output(captured, Val(:intensity))) ≈
        Array(graph_output(stream, Val(:intensity))) rtol=rtol atol=atol
    @test optical_sampling_m(owner) == first_sampling
    return nothing
end

@testset "CUDA scalar PROPER optical train graph capture" begin
    for T in (Float32, Float64)
        base = (resolution=32, pupil_resolution=16, diameter_m=0.008, wavelength_um=1.55)
        compare_cuda_proper_capture(T, merge(base, (steps=(),)))
        circles = ((radius_m=0.00012, x_m=0.00003, y_m=-0.00002),)
        mask_step = (op="sampled_mask", name="asymmetric sampled focal mask",
            amplitude="transmission", phase="retardance",
            sampling_m=(0.00001, 0.000014), center_pixel=(2.4, 3.7),
            outside_amplitude=0.2, outside_phase_rad=0.3)
        steps = ((op="circular_aperture", radius_m=0.0035),
            (op="lens", focal_length_m=0.12), (op="propagate", distance_m=0.12),
            mask_step, (op="propagate", distance_m=-0.03),
            (op="lens", focal_length_m=-0.2),
            (op="circular_mask", mode="transmit", circles=circles),
            (op="propagate", distance_m=0.17),
            merge(mask_step, (name="sampled relay mask", sampling_m=(0.0002, 0.0003))),
            (op="multiply", amplitude=0.8),
            (op="circular_obscuration", radius_m=0.00001),
            (op="propagate", distance_m=-0.04), (op="lens", focal_length_m=0.4),
            (op="propagate", distance_m=0.07))
        config = merge(base, (; steps, transmission="mask_amplitude", retardance="mask_phase"))
        mask_amplitude = T[0.4 + 0.3sin(0.7i + 0.2j)^2 for i in 1:5, j in 1:7]
        mask_phase = T[0.8sin(0.3i - 0.5j) for i in 1:5, j in 1:7]
        compare_cuda_proper_capture(T, config; bindings=(; mask_amplitude, mask_phase))
    end
end
