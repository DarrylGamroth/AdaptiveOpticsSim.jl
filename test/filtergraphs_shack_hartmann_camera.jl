# Run with the AOS, Proper, FilterGraphs and FilterGraphAlgorithms environment.
using AdaptiveOpticsSim, Proper, FilterGraphs, FilterGraphAlgorithms, LinearAlgebra, Test
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
policy = selection == "cpu" ? A.StreamGraphExecution() : A.CapturedGraphExecution()

@testset "Full Shack–Hartmann CMOS camera" begin
    schema = "test.shwfs/"
    pupil = A.pupil_opd_composition_node(:pupil; resolution=8,
        uncompensated_opd_schema=schema*"uncompensated/1",
        surface_opd_schema=schema*"dm/1", pupil_opd_schema=schema*"opd/1")
    optics = A.shack_hartmann_rate_node(:shwfs; T=Float32,
        resolution=8, telescope_diameter_m=1.22, n_lenslets=2, n_pix_subap=4,
        pixel_scale_arcsec=0.1, source_wavelength_m=0.55e-6,
        source_photon_irradiance_m2_s=2e8, opd_schema=schema*"opd/1",
        photon_rate_schema=schema*"photon-rate/1")
    detector = A.cmos_detector_acquisition_node(:camera; rows=8, columns=8,
        pixel_scale_arcsec=0.1, wavelength_m=0.55e-6, exposure_duration_s=0.001,
        quantum_efficiency=0.7, gain=2, dark_current_e_per_pixel_s=0.2,
        bits=12, full_well_e=90_000, photon_noise=true, readout_noise=true,
        readout_noise_e=2, column_readout_noise_e=0.3, row_readout_noise_e=0.5,
        rng_seed=121, photon_rate_schema=schema*"photon-rate/1",
        frame_schema=schema*"electrons/1")
    layout = A.builtin_graph_node_types().detector_frame_layout_f32(:layout,
        (; input_size=(8,8), output_size=(10,12), origin=(1,2),
            scale=1.5, bias=3, output_element_type="uint16", range=(0,65535),
            rounding="nearest_even", input_schema=schema*"electrons/1",
            output_schema=schema*"raw/1"), (;))
    host_opd = Float32[(i+2j)*1e-8 for i in 1:8, j in 1:8]
    host_surface = Float32[1e-8*sin(i-j) for i in 1:8, j in 1:8]
    native(values) = copyto!(B.allocate_device_array(target, Float32, size(values)...), values)
    definition = A.algorithm_graph((pupil, optics, detector, layout);
        inputs=(A.graph_input(:opd, :pupil=>:uncompensated_opd, native(host_opd)),
            A.graph_input(:surface, :pupil=>:surface_opd, native(host_surface))),
        links=(A.link(:pupil=>:pupil_opd, :shwfs=>:opd),
            A.link(:shwfs=>:photon_rate, :camera=>:photon_rate),
            A.link(:camera=>:frame, :layout=>:detector_frame)),
        outputs=(A.graph_output(:raw_camera, :layout=>:frame),
            A.graph_output(:photon_rate, :shwfs=>:photon_rate),
            A.graph_output(:detector_frame, :camera=>:frame),
            A.graph_output(:pupil_opd, :pupil=>:pupil_opd)))
    rewritten = O.with_filtergraphs_shack_hartmann_camera(definition; frame_output=:raw_camera)
    @test length(rewritten.nodes) == length(definition.nodes) == 4
    @test rewritten.links == definition.links
    for name in (:opd, :surface)
        @test A.graph_input_contract(rewritten, name) == A.graph_input_contract(definition, name)
    end
    for name in (:raw_camera, :photon_rate, :detector_frame, :pupil_opd)
        @test A.graph_output_contract(rewritten, name) == A.graph_output_contract(definition, name)
    end
    @test_throws ArgumentError O.shack_hartmann_camera_filter_graph(definition)
    @test_throws ArgumentError O.with_filtergraphs_shack_hartmann_camera(definition; frame_output=:missing)
    baseline = A.prepare_algorithm_graph(definition; target, execution=policy)
    bound = A.prepare_algorithm_graph(rewritten; target, execution=policy)
    recipe = O.shack_hartmann_camera_filter_graph(definition; frame_output=:raw_camera)
    checked = FilterGraphs.prepare_graph(recipe.document; algorithms=recipe.algorithms)
    direct = if selection == "cuda"
        owner = FilterGraphs.prepare_captured_graph(checked, CUDA.device())
        @test owner.executable.executable isa CUDA.CuGraphExec
        owner
    elseif selection == "amdgpu"
        FilterGraphs.prepare_captured_graph(checked, AMDGPU.device())
    else
        checked
    end
    @test length(direct.nodes) == 3
    if selection != "cpu"
        @test A.captured_graph_node_count(baseline) == A.captured_graph_node_count(bound) == 4
        @test direct isa FilterGraphs.PreparedCapturedGraph
    end
    schedule = A.PeriodicSchedule(; period_ns=2_000_000, phase_ns=400_000)
    drivers = (A.FixedStepModelTimeDriver(schedule), A.FixedStepModelTimeDriver(schedule))
    cases = (host_opd, host_opd, -host_opd)
    frame = zeros(UInt16, 10, 12)
    frames = Matrix{UInt16}[]
    for (index, opd) in enumerate(cases)
        timestamps = map((baseline, bound), drivers) do graph, driver
            A.copyto_graph_input!(graph, :opd, opd)
            A.step_graph_at!(graph, driver)
        end
        time = 400_000+(index-1)*2_000_000
        @test timestamps[1] == timestamps[2]
        @test A.model_nanoseconds(first(timestamps)) == time
        metadata = FilterGraphAlgorithms.SampleMetadata(index;
            timestamp_nanoseconds=time, exposure_duration_nanoseconds=1_000_000, terminal=true)
        @test FilterGraphs.process!((; frame), direct, (; opd=opd+host_surface), (; opd=metadata)) ==
            (; frame=true)
        @test frame == Array(A.graph_output(baseline, :raw_camera)) ==
            Array(A.graph_output(bound, :raw_camera))
        @test FilterGraphs.output_metadata(direct).frame == metadata
        @test all(iszero, frame[1,:]) && all(iszero, frame[10,:])
        @test all(iszero, frame[:,1:2]) && all(iszero, frame[:,11:12])
        for name in (:photon_rate, :detector_frame, :pupil_opd)
            @test Array(A.graph_output(bound, name)) == Array(A.graph_output(baseline, name))
        end
        push!(frames, copy(frame))
    end
    @test frames[1] != frames[2]
    if selection == "cpu"
        inputs = ((; opd=copy(host_opd+host_surface)), (; opd=copy(-host_opd+host_surface)))
        outputs = ((; frame=fill(typemax(UInt16), 10, 12)),
            (; frame=fill(typemax(UInt16), 10, 12)))
        loans = FilterGraphs.prepare_cpu_loans(recipe.document; algorithms=recipe.algorithms,
            input_slots=inputs, output_slots=outputs)
        first_slot, second_slot = loans.nodes_by_input
        @test first_slot[1].prepared.workspace.opd === inputs[1].opd
        @test second_slot[1].prepared.workspace.opd === inputs[2].opd
        @test first_slot[1].prepared.workspace.owner.prepared.workspace ===
            second_slot[1].prepared.workspace.owner.prepared.workspace
        @test first_slot[2].prepared === second_slot[2].prepared
        @test first_slot[1].prepared.workspace.photon_rate ===
            first_slot[2].prepared.workspace.owner.prepared.input.values
        for (index, input_slot, output_slot) in ((1,1,2), (2,1,1), (3,2,2))
            @test FilterGraphs.process_cpu_loan!(loans, input_slot, output_slot) == (; frame=true)
            @test outputs[output_slot].frame == frames[index]
        end
        saved = copy(outputs[2].frame)
        @test_throws BoundsError FilterGraphs.process_cpu_loan!(loans, 3, 1)
        @test outputs[2].frame == saved
        FilterGraphs.reset!(loans)
        @test FilterGraphs.process_cpu_loan!(loans, 1, 1) == (; frame=true)
        @test outputs[1].frame == first(frames)
    end
    prior = copy(frame)
    prior_metadata = FilterGraphs.output_metadata(direct)
    @test_throws DimensionMismatch FilterGraphs.process!((; frame), direct, (; opd=zeros(Float32, 7, 8)))
    @test frame == prior
    @test FilterGraphs.output_metadata(direct) == prior_metadata

    for (graph, driver) in zip((baseline, bound), drivers)
        A.reset_graph!(graph)
        A.reset_model_time!(driver)
    end
    FilterGraphs.reset!(direct)
    for (index, opd) in enumerate(cases)
        for (graph, driver) in zip((baseline, bound), drivers)
            A.copyto_graph_input!(graph, :opd, opd)
            @test A.model_nanoseconds(A.step_graph_at!(graph, driver)) ==
                400_000+(index-1)*2_000_000
        end
        @test FilterGraphs.process!((; frame), direct, (; opd=opd+host_surface)) == (; frame=true)
        @test frame == frames[index] == Array(A.graph_output(baseline, :raw_camera)) ==
            Array(A.graph_output(bound, :raw_camera))
    end
end
