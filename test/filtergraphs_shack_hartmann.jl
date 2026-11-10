# Run with AdaptiveOpticsSim, Proper, FilterGraphAlgorithms and FilterGraphs.
# The first argument selects cpu (default), cuda, or amdgpu.
using AdaptiveOpticsSim, Proper, FilterGraphs, FilterGraphAlgorithms, LinearAlgebra, Test
using AdaptiveOpticsSim.AlgorithmGraphs

const A = AdaptiveOpticsSim.AlgorithmGraphs
const B = AdaptiveOpticsSim.Backends
const Ext = Base.get_extension(AdaptiveOpticsSim, :AdaptiveOpticsSimFilterGraphsExt)

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
policy = selection == "cpu" ? A.StreamGraphExecution() : A.CapturedGraphExecution()

@testset "FilterGraphs Shack–Hartmann optical rate" begin
    node = A.shack_hartmann_rate_node(:shwfs; T=Float32,
        resolution=8, telescope_diameter_m=1.22,
        n_lenslets=2, n_pix_subap=4, pixel_scale_arcsec=0.1,
        source_wavelength_m=0.55e-6, source_photon_irradiance_m2_s=1e6,
        opd_schema="test.pupil-opd/1", photon_rate_schema="test.photon-rate/1")
    opd_shape = (8, 8)
    photon_rate_shape = (8, 8)
    buffer(shape) = fill!(B.allocate_device_array(target, Float32, shape...), 0)
    opd = buffer(opd_shape)
    photon_rate = zeros(Float32, photon_rate_shape)
    definition = A.algorithm_graph((node,);
        inputs=(A.graph_input(:opd, :shwfs=>:opd, opd),),
        outputs=(A.graph_output(:photon_rate, :shwfs=>:photon_rate),))

    baseline = A.prepare_algorithm_graph(definition; target, execution=policy)
    bound = A.prepare_algorithm_graph(Ext.replace_camera_nodes(definition, (:shwfs,));
        target, execution=policy)
    recipe = AdaptiveOpticsSim.Optics.shack_hartmann_rate_filter_graph(node)
    checked = FilterGraphs.prepare_graph(recipe.document; algorithms=recipe.algorithms)
    direct = if selection == "cuda"
        captured = FilterGraphs.prepare_captured_graph(checked, CUDA.device())
        @test captured.executable.executable isa CUDA.CuGraphExec
        captured
    elseif selection == "amdgpu"
        FilterGraphs.prepare_captured_graph(checked, AMDGPU.device())
    else
        checked
    end
    if selection != "cpu"
        @test direct isa FilterGraphs.PreparedCapturedGraph
        @test A.captured_graph_node_count(baseline) == A.captured_graph_node_count(bound) == 1
    end
    @test recipe.algorithms == (Ext.ShackHartmannRateF32,)
    @test A.graph_input_contract(definition, :opd).shape == opd_shape
    @test A.graph_output_contract(definition, :photon_rate).shape == photon_rate_shape

    cases = (
        zeros(Float32, opd_shape),
        Float32[(i - 1) * 1e-8 for i in 1:8, j in 1:8],
        Float32[2e-8 * sin(i + 2j) for i in 1:8, j in 1:8],
    )
    frames = map(enumerate(cases)) do (index, input)
        metadata = (; opd=FilterGraphAlgorithms.SampleMetadata(index;
            timestamp_nanoseconds=(index - 1) * 500_000,
            exposure_duration_nanoseconds=500_000, terminal=true))
        A.copyto_graph_input!(baseline, :opd, input)
        A.step_graph!(baseline)
        A.copyto_graph_input!(bound, :opd, input)
        A.step_graph!(bound)
        @test FilterGraphs.process!((; photon_rate), direct, (; opd=input), metadata) ===
            (; photon_rate=true)
        expected = Array(A.graph_output(baseline, :photon_rate))
        @test photon_rate == expected == Array(A.graph_output(bound, :photon_rate))
        @test all(isfinite, photon_rate) && any(>(0), photon_rate)
        @test FilterGraphs.output_metadata(direct).photon_rate == metadata.opd
        copy(photon_rate)
    end
    @test frames[1] != frames[2] != frames[3]

    FilterGraphs.reset!(direct)
    @test FilterGraphs.process!((; photon_rate), direct, (; opd=first(cases))) ===
        (; photon_rate=true)
    @test photon_rate == first(frames)

    input_binding = buffer(opd_shape)
    output_binding = buffer(photon_rate_shape)
    bound = Ext.prepare_host_shack_hartmann(node, input_binding, output_binding)
    @test bound.workspace.opd === input_binding
    @test bound.workspace.photon_rate === output_binding
end
