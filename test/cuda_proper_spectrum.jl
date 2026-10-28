# Run explicitly with this checkout, Proper and CUDA in the selected environment.
# JULIA_NUM_THREADS=1 julia --startup-file=no --project=<environment> test/cuda_proper_spectrum.jl
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

function spectrum_definition(config, opd, amplitude)
    node = proper_optical_train_node_factory(T=eltype(opd))(:optics, config, (;))
    return algorithm_graph((node,); name=:cuda_proper_spectrum,
        inputs=(graph_input(:opd, :optics => :pupil_opd, opd),
            graph_input(:amplitude, :optics => :pupil_amplitude, amplitude)),
        outputs=(graph_output(:intensity, :optics => :output),))
end

function compare_spectrum_capture(::Type{T}, config) where {T}
    p = config.pupil_resolution
    initial_opd = T[3e-8sin(0.3i + 0.7j) for i in 1:p, j in 1:p]
    initial_amplitude = T[0.2 + 0.6sin(0.1i + 0.3j)^2 for i in 1:p, j in 1:p]
    host_opd, host_amplitude = copy(initial_opd), copy(initial_amplitude)
    opd, amplitude = CuArray(initial_opd), CuArray(initial_amplitude)
    definition = spectrum_definition(config, opd, amplitude)
    target = compute_device(opd)
    streamed = prepare_algorithm_graph(definition; target, execution=StreamGraphExecution())
    captured = prepare_algorithm_graph(definition; target, execution=CapturedGraphExecution())
    host = prepare_algorithm_graph(spectrum_definition(config, host_opd, host_amplitude))
    owner = AlgorithmGraphs.prepared_graph_node(captured, Val(:optics))
    @test captured_graph_node_count(captured) == 1
    @test AlgorithmGraphs.graph_node_capture_capability(owner) isa AlgorithmGraphs.GraphNodeCaptureSafe
    @test all(channel -> AlgorithmGraphs.graph_node_capture_capability(channel) isa
        AlgorithmGraphs.GraphNodeCaptureSafe, owner.channels)
    scratch, overlaps = owner.scratch, owner.overlaps
    channel_outputs = map(channel -> channel.output, owner.channels)
    pitch = optical_sampling_m(owner)
    if hasproperty(config, :output_sampling_m)
        @test pitch === T(config.output_sampling_m)
    end
    rtol, atol = T === Float32 ? (8e-4, 8e-6) : (8e-11, 8e-12)
    first_image = zeros(T, config.resolution, config.resolution)
    for frame in 1:4
        host_opd .= initial_opd .+ T(frame-1) .* T[5e-8(i - 0.4j)/p for i in 1:p, j in 1:p]
        host_amplitude .= initial_amplitude .* T(1 - 0.17(frame-1))
        copyto!(opd, host_opd)
        copyto!(amplitude, host_amplitude)
        CUDA.synchronize()
        step_graph!(host)
        step_graph!(streamed)
        step_graph!(captured)
        actual = Array(graph_output(captured, Val(:intensity)))
        @test actual ≈ Array(graph_output(streamed, Val(:intensity))) rtol=rtol atol=atol
        @test actual ≈ graph_output(host, Val(:intensity)) rtol=rtol atol=atol
        @test all(isfinite, actual)
        @test sum(actual) > 0
        @test optical_sampling_m(owner) == pitch
        @test owner.scratch === scratch
        @test owner.overlaps === overlaps
        @test all(i -> owner.channels[i].output === channel_outputs[i], eachindex(owner.channels))
        if frame == 1
            first_image .= actual
        else
            @test !isapprox(actual, first_image; rtol=1e-3, atol=atol)
        end
    end
    copyto!(opd, initial_opd)
    copyto!(amplitude, initial_amplitude)
    CUDA.synchronize()
    reset_graph!(streamed)
    reset_graph!(captured)
    @test all(iszero, Array(graph_output(captured, Val(:intensity))))
    @test all(iszero, Array(owner.scratch))
    step_graph!(streamed)
    step_graph!(captured)
    @test Array(graph_output(captured, Val(:intensity))) ≈ first_image rtol=rtol atol=atol
    @test Array(graph_output(captured, Val(:intensity))) ≈
        Array(graph_output(streamed, Val(:intensity))) rtol=rtol atol=atol
    @test optical_sampling_m(owner) == pitch
    return nothing
end

@testset "CUDA finite-band PROPER capture" begin
    for T in (Float32, Float64)
        base = (resolution=32, pupil_resolution=16, diameter_m=0.008, wavelength_um=1.5,
            steps=((op="lens", focal_length_m=0.14), (op="propagate", distance_m=0.14)))
        # A single default-transmission channel retains the scalar grid and path.
        single = merge(base, (channels=((wavelength_um=1.5, weight=3.0, steps=base.steps),),))
        compare_spectrum_capture(T, single)
        # Fixed physical pitch clips the broader red support. Different optical
        # trains, photon fractions and filter attenuation all survive replay.
        blue = ((op="multiply", amplitude=0.6), (op="lens", focal_length_m=0.12),
            (op="propagate", distance_m=0.12))
        red = ((op="lens", focal_length_m=0.18), (op="propagate", distance_m=0.18),
            (op="multiply", amplitude=0.85))
        channels = ((wavelength_um=1.0, weight=2.0, transmission=0.35, steps=blue),
            (wavelength_um=2.1, weight=5.0, transmission=0.8, steps=red),
            (wavelength_um=1.6, weight=1.0, transmission=0.0, steps=base.steps))
        compare_spectrum_capture(T, merge(base, (; channels, output_sampling_m=1.2e-5)))
        # Fixed-grid single-channel selection also uses the gather path.
        compare_spectrum_capture(T, merge(base, (output_sampling_m=1.2e-5,)))
    end
end

@testset "CUDA cell-integrated overlap conservation and clipping" begin
    extension = Base.get_extension(AdaptiveOpticsSim, :AdaptiveOpticsSimProperExt)
    for T in (Float32, Float64), (source_pitch, output_pitch) in ((1, 1), (1, 2), (2, 1), (1.25, 1))
        image = reshape(T.(1:64), 8, 8)
        overlaps = extension._cell_overlaps(8, T(source_pitch), T(output_pitch))
        expected = zeros(T, 8, 8)
        extension._add_channel_intensity!(expected, image, similar(expected), overlaps, T(0.3))
        intensity = CuArray(image)
        target = compute_device(intensity)
        gather = extension._prepare_cell_overlaps(overlaps, 8, target)
        output, scratch = CUDA.zeros(T, 8, 8), CUDA.zeros(T, 8, 8)
        extension._add_channel_intensity!(output, intensity, scratch, gather, T(0.3))
        actual = Array(output)
        @test actual ≈ expected rtol=50eps(T)
        if output_pitch >= source_pitch
            @test sum(actual) ≈ T(0.3)sum(image) rtol=50eps(T)
        else
            @test sum(actual) < T(0.3)sum(image)
        end
        if gather !== nothing
            @test length(gather.offsets) == 9
            @test length(gather.sources) == length(overlaps)
            @test length(gather.sources) <= 15
        end
    end
    edge = zeros(Float32, 8, 8); edge[1, 1] = 1
    intensity = CuArray(edge)
    overlaps = extension._cell_overlaps(8, 2f0, 1f0)
    gather = extension._prepare_cell_overlaps(overlaps, 8, compute_device(intensity))
    output, scratch = CUDA.zeros(Float32, 8, 8), CUDA.zeros(Float32, 8, 8)
    extension._add_channel_intensity!(output, intensity, scratch, gather, 1f0)
    @test iszero(sum(Array(output)))
end
