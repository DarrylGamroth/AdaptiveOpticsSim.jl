using YAML

function hgcdte_avalanche_test_graph(input; kwargs...)
    config = merge((rows=size(input, 1), columns=size(input, 2),
            pixel_sampling_m=(24e-6, 25e-6), spectral_channel=:test_hgcdte_band,
            exposure_duration_s=0.5, quantum_efficiency=0.75, gain=2,
            avalanche_gain=4, rng_seed=0xabc1,
            photon_rate_schema="test.hgcdte.photon-rate/1",
            frame_schema="test.hgcdte.frame/1", photon_noise=false, T=eltype(input)),
        (; kwargs...))
    node = hgcdte_avalanche_detector_acquisition_node(:detector; config...)
    return algorithm_graph((node,); name=:hgcdte_avalanche_test,
        inputs=(graph_input(:rate, :detector => :photon_rate, input),),
        outputs=(graph_output(:frame, :detector => :frame),))
end

@testset "HgCdTe linear-avalanche complete-frame graph node" begin
    graphs = AdaptiveOpticsSim.AlgorithmGraphs
    for T in (Float32, Float64)
        input = reshape(T.(1:15), 3, 5)
        graph = prepare_algorithm_graph(hgcdte_avalanche_test_graph(input))
        owner = prepared_graph_node(graph, Val(:detector))
        output = graph_output(graph, Val(:frame))
        @test owner.prepared.input.values === input
        @test owner.prepared.detector.params.sensor isa HgCdTeAvalancheArraySensor
        @test owner.prepared.input.metadata.sampling == (T(24e-6), T(25e-6))
        @test owner.prepared.input.metadata.spectral ==
            IntegratedSpectralChannel(:test_hgcdte_band)
        @test @inferred(step_graph!(graph)) === graph
        @test output == input .* T(3) # rate × exposure × QE × avalanche × gain
        input .*= T(2)
        step_graph!(graph)
        @test output == input .* T(3)
        saturated = prepare_algorithm_graph(hgcdte_avalanche_test_graph(
            fill(T(20), 3, 5); full_well_e=8))
        step_graph!(saturated)
        @test graph_output(saturated, Val(:frame)) == fill(T(16), 3, 5)

        # Oracle is the independently prepared numerical sensor, including its
        # stochastic multiplication, saturation, persistence and nonlinearity.
        for model in (ConditionalGammaAvalancheMultiplication(),
                ClippedGaussianAvalancheMultiplicationApproximation())
            photons = fill(T(200), 3, 5)
            graph = prepare_algorithm_graph(hgcdte_avalanche_test_graph(photons;
                photon_noise=true, readout_noise=true, readout_noise_e=0.7,
                dark_current_e_per_pixel_s=2, glow_e_per_pixel_s=3,
                excess_noise_factor=1.2, multiplication_model=model,
                persistence_model=ExponentialPersistence(0.1, 0.5),
                nonlinearity_model=SaturatingFrameNonlinearity(0.001)))
            owner = prepared_graph_node(graph, Val(:detector))
            detector = Detector(exposure_duration=T(0.5), qe=T(0.75), gain=T(2),
                dark_current=T(2), noise=NoisePhotonReadout(T(0.7)),
                sensor=HgCdTeAvalancheArraySensor(avalanche_gain=4,
                    excess_noise_factor=1.2, multiplication_model=model,
                    glow_rate=3, persistence_model=ExponentialPersistence(0.1, 0.5), T=T),
                nonlinearity_model=SaturatingFrameNonlinearity(T(0.001)),
                response_model=NullFrameResponse(), T=T)
            reference = prepare_detector_acquisition(detector, owner.prepared.input)
            rng = AdaptiveOpticsSim.Backends.runtime_rng(0xabc1)
            first = nothing
            for factor in (1, 2, 0)
                fill!(photons, T(200factor))
                step_graph!(graph)
                @test graph_output(graph, Val(:frame)) == capture!(reference, rng)
                factor == 1 && (first = copy(graph_output(graph, Val(:frame))))
            end
            @test any(>(0), owner.latent_buffer)
            reset_graph!(graph)
            @test all(iszero, owner.latent_buffer)
            @test all(iszero, graph_output(graph, Val(:frame)))
            @test graph_step_sequence(graph) == 0
            fill!(photons, T(200))
            step_graph!(graph)
            @test graph_output(graph, Val(:frame)) == first
            @test graphs.graph_node_capture_capability(owner) isa
                graphs.GraphNodeCaptureUnsupported
        end
    end
    input = ones(Float32, 3, 5)
    for bad in ((rows=0,), (columns=-1,), (quantum_efficiency=1.1,),
            (exposure_duration_s=0,), (gain=0,), (avalanche_gain=0.9,),
            (avalanche_gain=Inf,), (excess_noise_factor=0.9,),
            (excess_noise_factor=NaN,), (dark_current_e_per_pixel_s=-1,),
            (glow_e_per_pixel_s=-1,), (pixel_sampling_m=(0, 1e-6),),
            (spectral_channel=Symbol(""),), (bits=8,), (bits=0, full_well_e=10),
            (readout_noise_e=2,), (readout_noise=true, readout_noise_e=0),
            (rng_seed=-1,), (T=Float16,), (photon_rate_schema="",),
            (frame_schema="",))
        @test_throws AlgorithmGraphError hgcdte_avalanche_test_graph(input; bad...)
    end
    @test_throws MethodError hgcdte_avalanche_test_graph(input; sampling_mode=SingleRead())
    @test Base.isexported(graphs, :hgcdte_avalanche_detector_acquisition_node)
    @test Base.ispublic(graphs, :HgCdTeAvalancheDetectorAcquisitionNode)
    @test !Base.ispublic(AdaptiveOpticsSim, :HgCdTeAvalancheDetectorAcquisitionNode)
end

@testset "Declarative HgCdTe linear-avalanche acquisition" begin
    input = fill(200f0, 3, 5)
    path = joinpath(dirname(@__DIR__), "graph_files", "hgcdte_avalanche.yaml")
    document = YAML.load_file(path)
    config = document["nodes"][1]["config"]
    bindings = (; photons=input)
    from_file = prepare_algorithm_graph(load_algorithm_graph(path; bindings))
    from_dict = prepare_algorithm_graph(load_algorithm_graph(document; bindings))
    step_graph!(from_file)
    step_graph!(from_dict)
    @test graph_output(from_file, Val(:frame)) == graph_output(from_dict, Val(:frame))
    for model in ("conditional_gamma", "clipped_gaussian_approximation")
        config["multiplication_model"] = model
        actual = prepare_algorithm_graph(load_algorithm_graph(document; bindings))
        numerical_model = model == "conditional_gamma" ?
            ConditionalGammaAvalancheMultiplication() :
            ClippedGaussianAvalancheMultiplicationApproximation()
        expected = prepare_algorithm_graph(hgcdte_avalanche_test_graph(input;
            excess_noise_factor=1.2, multiplication_model=numerical_model,
            photon_noise=true, readout_noise=true, readout_noise_e=0.7,
            dark_current_e_per_pixel_s=2, glow_e_per_pixel_s=3, bits=12, full_well_e=1000,
            persistence_model=ExponentialPersistence(0.1, 0.5),
            nonlinearity_model=SaturatingFrameNonlinearity(0.001)))
        for _ in 1:2
            step_graph!(actual)
            step_graph!(expected)
            @test graph_output(actual, Val(:frame)) == graph_output(expected, Val(:frame))
        end
    end
    for (field, invalid) in (("multiplication_model", "unknown"),
            ("avalanche_gain", 0.5), ("excess_noise_factor", Inf),
            ("pixel_sampling_m", [24e-6]), ("shutter", "rolling"),
            ("readout_noise", 1), ("bits", true), ("unknown_setting", 1),
            ("read_duration_s", 0.01), ("sampling_mode", "fowler"))
        bad = deepcopy(document)
        bad["nodes"][1]["config"][field] = invalid
        @test_throws AlgorithmGraphError load_algorithm_graph(bad; bindings)
    end
end
