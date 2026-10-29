# Run explicitly with this checkout, Proper and AMDGPU in the selected environment.
# JULIA_NUM_THREADS=1 julia --startup-file=no --project=<environment> test/amdgpu_proper_capture.jl
using Test, AdaptiveOpticsSim, Proper, AMDGPU, LinearAlgebra
using AdaptiveOpticsSim.AlgorithmGraphs
using AdaptiveOpticsSim.Backends

AMDGPU.allowscalar(false)
BLAS.set_num_threads(1)
Backends.set_fft_provider_threads!(1)

function amdgpu_proper_definition(config, opd, amplitude; bindings=(;))
    node = proper_optical_train_node_factory(T=eltype(opd))(:optics, config, (;), bindings)
    return algorithm_graph((node,); name=:amdgpu_proper_capture,
        inputs=(graph_input(:opd, :optics => :pupil_opd, opd),
            graph_input(:amplitude, :optics => :pupil_amplitude, amplitude)),
        outputs=(graph_output(:intensity, :optics => :output),))
end

function compare_amdgpu_proper_capture(::Type{T}, config; bindings=(;)) where {T}
    p = config.pupil_resolution
    initial_opd = T[3e-8sin(0.3i + 0.7j) for i in 1:p, j in 1:p]
    initial_amplitude = T[0.2 + 0.6sin(0.1i + 0.3j)^2 for i in 1:p, j in 1:p]
    host_opd, host_amplitude = copy(initial_opd), copy(initial_amplitude)
    opd, amplitude = ROCArray(initial_opd), ROCArray(initial_amplitude)
    definition = amdgpu_proper_definition(config, opd, amplitude; bindings)
    target = compute_device(opd)
    streamed = prepare_algorithm_graph(definition; target, execution=StreamGraphExecution())
    captured = prepare_algorithm_graph(definition; target, execution=CapturedGraphExecution())
    host = prepare_algorithm_graph(amdgpu_proper_definition(config, host_opd, host_amplitude; bindings))
    owner = AlgorithmGraphs.prepared_graph_node(captured, :optics)
    @test captured_graph_node_count(captured) == 1
    @test AlgorithmGraphs.graph_node_capture_capability(owner) isa AlgorithmGraphs.GraphNodeCaptureSafe
    channels = hasproperty(owner, :channels) ? owner.channels : (owner,)
    retained = map(channels) do channel
        fft = channel.workspace.context.workspace.fft
        (; field=channel.workspace.wavefront.field, scratch=fft.scratch,
            forward=fft.forward_plan, backward=fft.backward_plan, output=channel.output)
    end
    pitch = optical_sampling_m(owner)
    @test isfinite(pitch) && pitch > 0
    if hasproperty(config, :output_sampling_m)
        @test pitch == T(config.output_sampling_m)
    end
    rtol, atol = T === Float32 ? (8e-4, 8e-6) : (8e-11, 8e-12)
    first_image = zeros(eltype(graph_output(captured, :intensity)), config.resolution, config.resolution)
    for frame in 1:4
        host_opd .= initial_opd .+ T(frame-1) .* T[5e-8(i - 0.4j)/p for i in 1:p, j in 1:p]
        host_amplitude .= initial_amplitude .* T(1 - 0.17(frame-1))
        copyto!(opd, host_opd)
        copyto!(amplitude, host_amplitude)
        AMDGPU.synchronize() # Publish inputs before the retained streams consume them.
        step_graph!(host)
        step_graph!(streamed)
        step_graph!(captured)
        actual = Array(graph_output(captured, :intensity))
        @test actual ≈ Array(graph_output(streamed, :intensity)) rtol=rtol atol=atol
        @test actual ≈ graph_output(host, :intensity) rtol=rtol atol=atol
        @test all(isfinite, actual) && sum(abs2, actual) > 0
        @test optical_sampling_m(owner) == pitch
        for (channel, previous) in zip(channels, retained)
            fft = channel.workspace.context.workspace.fft
            @test channel.workspace.wavefront.field === previous.field
            @test fft.scratch === previous.scratch
            @test fft.forward_plan === previous.forward
            @test fft.backward_plan === previous.backward
            @test channel.output === previous.output
        end
        if frame == 1
            first_image .= actual
        else
            @test !isapprox(actual, first_image; rtol=1e-3, atol=atol)
        end
    end
    copyto!(opd, initial_opd)
    copyto!(amplitude, initial_amplitude)
    AMDGPU.synchronize()
    reset_graph!(streamed)
    reset_graph!(captured)
    @test all(iszero, Array(graph_output(captured, :intensity)))
    step_graph!(streamed)
    step_graph!(captured)
    @test Array(graph_output(captured, :intensity)) ≈ first_image rtol=rtol atol=atol
    @test Array(graph_output(captured, :intensity)) ≈ Array(graph_output(streamed, :intensity)) rtol=rtol atol=atol
    @test optical_sampling_m(owner) == pitch
    return nothing
end

@testset "AMDGPU scalar and finite-band Proper HIP capture" begin
    for T in (Float32, Float64)
        base = (resolution=32, pupil_resolution=16, diameter_m=0.008, wavelength_um=1.5)
        compare_amdgpu_proper_capture(T, merge(base, (steps=(),)))
        compare_amdgpu_proper_capture(T, merge(base, (pupil_resolution=32, steps=(),)))
        compare_amdgpu_proper_capture(T, merge(base, (output_quantity="electric_field",
            steps=((op="propagate", distance_m=0.01),))))
        mask = (op="sampled_mask", name="asymmetric focal mask", amplitude="transmission",
            sampling_m=(0.00001, 0.000014), center_pixel=(2.4, 3.7),
            outside_amplitude=0.2, outside_phase_rad=0.3)
        steps = ((op="circular_aperture", radius_m=0.0035),
            (op="lens", focal_length_m=0.12), (op="propagate", distance_m=0.12),
            mask, (op="multiply", amplitude=0.8), (op="propagate", distance_m=-0.03))
        transmission = T[0.4 + 0.3sin(0.7i + 0.2j)^2 for i in 1:5, j in 1:7]
        compare_amdgpu_proper_capture(T, merge(base, (; steps, transmission="transmission")); bindings=(; transmission))
        blue = ((op="multiply", amplitude=0.6), (op="lens", focal_length_m=0.12),
            (op="propagate", distance_m=0.12))
        red = ((op="lens", focal_length_m=0.18), (op="propagate", distance_m=0.18))
        channels = ((wavelength_um=1.0, weight=2.0, transmission=0.35, steps=blue),
            (wavelength_um=2.1, weight=5.0, transmission=0.8, steps=red),
            (wavelength_um=1.6, weight=1.0, transmission=0.0, steps=blue))
        compare_amdgpu_proper_capture(T, merge(base, (steps=blue, output_sampling_m=1.2e-5,)))
        compare_amdgpu_proper_capture(T, merge(base, (steps=blue, channels, output_sampling_m=1.2e-5)))
    end
end
@test !any(p.name == "CUDA" for p in keys(Base.loaded_modules))
