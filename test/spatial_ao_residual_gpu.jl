# Focused optional-backend selector. No GPU package is loaded without an explicit
# target: julia --project=<validation-env> test/spatial_ao_residual_gpu.jl cuda|amdgpu
using Test
using AdaptiveOpticsSim
using AdaptiveOpticsSim.AlgorithmGraphs
using AdaptiveOpticsSim.Atmospheres
using AdaptiveOpticsSim.Backends
using LinearAlgebra

length(ARGS) == 1 && only(ARGS) in ("cuda", "amdgpu") ||
    error("expected one target: cuda or amdgpu")
if only(ARGS) == "cuda"
    import CUDA
    CUDA.functional() || error("CUDA is not functional")
    CUDA.device!(0)
    target = Backends.AcceleratorComputeDevice(CUDABackend(), 0)
else
    import AMDGPU
    AMDGPU.functional() || error("AMDGPU is not functional")
    target = Backends.AcceleratorComputeDevice(AMDGPUBackend(), 1)
end
BLAS.set_num_threads(1)
Backends.set_fft_provider_threads!(1)

# Reuse the same CPU reference equations and declarations, then exercise the
# actual device FFT at reset and the existing graph replay paths independently.
include("testsets/spatial_ao_residual.jl")

@testset "spatial AO residual device Fourier response and reset" begin
    rejection = SpatialAORejection(control_cutoff_cycles_per_m=1,
        residual_amplitude=0.25)
    graph = spatial_ao_test_graph(rejection; target)
    owner = SpatialAOGraphs.prepared_graph_node(graph, Val(:residual))
    filter = owner.rejection
    @test compute_device(graph) == target
    @test compute_device(filter.workspace.buffer) == target
    @test compute_device(filter.transfer) == target
    n = size(first(filter.screens), 1)
    for (kx, ky, amplitude) in ((0, 0, 0.0), (1, 2, 0.25), (0, 3, 0.25), (4, 1, 1.0))
        wave = Float32[cos(2pi * (kx * (i - 1) + ky * (j - 1)) / n)
            for i in 1:n, j in 1:n]
        Backends._with_prepared_device_execution_context(graph.context) do
            for screen in filter.screens
                copyto!(screen, wave)
            end
            SpatialAOAtmospheres._apply_spatial_ao_rejection!(filter)
            Backends._synchronize_prepared_device_execution_context!(graph.context)
        end
        for screen in filter.screens
            @test Array(screen) ≈ Float32(amplitude) .* wave atol=2e-5 rtol=2e-5
        end
    end
    reset_graph!(graph)
    first_screens = map(Array, filter.screens)
    step_graph!(graph)
    first_output = Array(graph_output(graph, Val(:opd)))
    @test all(isfinite, first_output)
    @test !all(iszero, first_output)
    step_graph!(graph)
    @test Array(graph_output(graph, Val(:opd))) != first_output
    reset_graph!(graph)
    @test map(Array, filter.screens) == first_screens
    @test all(iszero, Array(graph_output(graph, Val(:opd))))
    step_graph!(graph)
    @test Array(graph_output(graph, Val(:opd))) == first_output
end

@testset "spatial AO residual captured replay" begin
    rejection = SpatialAORejection(control_cutoff_cycles_per_m=1,
        residual_amplitude=0.25)
    ordinary = spatial_ao_test_graph(rejection; target)
    captured = spatial_ao_test_graph(rejection; target, execution=CapturedGraphExecution())
    @test captured_graph_node_count(captured) == 1
    expected = Matrix{Float32}[]
    for _ in 1:5
        step_graph!(ordinary)
        step_graph!(captured)
        frame = Array(graph_output(captured, Val(:opd)))
        @test frame ≈ Array(graph_output(ordinary, Val(:opd))) atol=2e-12 rtol=2e-5
        push!(expected, frame)
    end
    reset_graph!(captured)
    for frame in expected
        step_graph!(captured)
        @test Array(graph_output(captured, Val(:opd))) == frame
    end
end
println("SPATIAL_AO_GPU_PASS target=", only(ARGS))
