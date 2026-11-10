# Run with AdaptiveOpticsSim, Proper, FilterGraphAlgorithms and FilterGraphs.
# The first argument selects cpu (default), cuda, or amdgpu.
using AdaptiveOpticsSim, Proper, FilterGraphs, FilterGraphAlgorithms, LinearAlgebra, Test
using AdaptiveOpticsSim.AlgorithmGraphs
const A = AdaptiveOpticsSim.AlgorithmGraphs
const B = AdaptiveOpticsSim.Backends
const O = AdaptiveOpticsSim.Optics

BLAS.set_num_threads(1)
B.set_fft_provider_threads!(1)

selection = get(ARGS, 1, "cpu")
if selection == "cuda"
    using CUDA
    @test CUDA.functional()
    CUDA.allowscalar(false)
    target = B.AcceleratorComputeDevice(B.CUDABackend(), 0)
elseif selection == "amdgpu"
    using AMDGPU
    @test AMDGPU.functional()
    AMDGPU.allowscalar(false)
    target = B.AcceleratorComputeDevice(B.AMDGPUBackend(), AMDGPU.device_id(AMDGPU.device()))
else
    @assert selection == "cpu"
    target = B.HostComputeDevice()
end
policy = selection == "cpu" ? StreamGraphExecution() : CapturedGraphExecution()

@testset "Direct and AOS-bound FilterGraphAlgorithms PWFS" begin
    node = pyramid_rate_node(:pyramid; T=Float32, resolution=8, telescope_diameter_m=1.22,
        pupil_samples=4, modulation=2, modulation_points=4, n_pix_separation=2, n_pix_edge=1,
        source_wavelength_m=0.55e-6, source_photon_irradiance_m2_s=1.0,
        opd_schema="test.pupil-opd/1", photon_rate_schema="test.photon-rate/1")
    opd = fill!(B.allocate_device_array(target, Float32, 8, 8), 0)
    definition = algorithm_graph((node,); inputs=(graph_input(:opd, :pyramid=>:opd, opd),),
        outputs=(graph_output(:photon_rate, :pyramid=>:photon_rate),))
    baseline = prepare_algorithm_graph(definition; target, execution=policy)
    substituted = O.with_filtergraphs_pwfs(definition)
    @test A.graph_input_contract(substituted, :opd) == A.graph_input_contract(definition, :opd)
    @test A.graph_output_contract(substituted, :photon_rate) == A.graph_output_contract(definition, :photon_rate)
    @test_throws ArgumentError O.with_filtergraphs_pwfs(definition; node=:missing)
    bridge = prepare_algorithm_graph(substituted; target, execution=policy)
    recipe = O.pwfs_filter_graph(node)
    checked = FilterGraphs.prepare_graph(recipe.document; algorithms=recipe.algorithms)
    direct = if selection == "cuda"
        result = FilterGraphs.prepare_captured_graph(checked, CUDA.device())
        @test result.executable.executable isa CUDA.CuGraphExec
        result
    elseif selection == "amdgpu"
        FilterGraphs.prepare_captured_graph(checked, AMDGPU.device())
    else
        checked
    end
    if selection != "cpu"
        @test direct isa FilterGraphs.PreparedCapturedGraph
        @test captured_graph_node_count(bridge) == 1
    end
    photon_rate = zeros(Float32, size(graph_output(baseline, :photon_rate)))
    cases = (
        zeros(Float32, 8, 8),
        Float32[(i + 2j) * 1e-8 for i in 1:8, j in 1:8],
        Float32[2e-8 * sin(i + j) for i in 1:8, j in 1:8],
    )
    rates = map(cases) do host_opd
        A.copyto_graph_input!(baseline, :opd, host_opd)
        step_graph!(baseline)
        A.copyto_graph_input!(bridge, :opd, host_opd)
        step_graph!(bridge)
        @test FilterGraphs.process!((; photon_rate), direct, (; opd=host_opd)) === (; photon_rate=true)
        expected = Array(graph_output(baseline, :photon_rate))
        @test all(isfinite, photon_rate) && any(>(0), photon_rate)
        @test photon_rate == expected
        @test Array(graph_output(bridge, :photon_rate)) == expected
        copy(photon_rate)
    end
    @test rates[1] != rates[2] != rates[3]
    reset_graph!(bridge)
    FilterGraphs.reset!(direct)
    @test FilterGraphs.process!((; photon_rate), direct, (; opd=first(cases))) === (; photon_rate=true)
    @test photon_rate == first(rates)
end
