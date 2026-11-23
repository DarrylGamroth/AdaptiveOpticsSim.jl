using AdaptiveOpticsSim.AlgorithmGraphs: warmup_external_pupil!

# Use the same optical-train/pixel-integration topology in both chopper phases.
function warmup_test_camera(mask)
    source = proper_camera_test_source()
    camera = AlgorithmGraphs.compose_proper_camera_graph(source; bindings=(; mask),
        output_shape=(8, 8), output_sampling_m=(1e-5, 1e-5),
        output_schema="test.camera-cell-flux/1", name="camera")
    return prepare_algorithm_graph(load_algorithm_graph(camera; bindings=(
        mask=mask, opd=zeros(Float32, 16, 16), amplitude=ones(Float32, 16, 16))))
end

function warmup_test_owner()
    fringed = warmup_test_camera(Float32[0.6 + 0.2cos(i + 2j) for i in 1:16, j in 1:16])
    unfringed = warmup_test_camera(fill(0.5f0, 16, 16))
    llowfs = warmup_test_camera(ones(Float32, 16, 16))
    plant = prepare_external_pupil((; llowfs, scc=fringed);
        exposure_durations=(llowfs=ModelDuration(2), scc=ModelDuration(10)),
        observations=(llowfs=(output=:intensity, sampling_node=:camera_pixels),
            scc=(output=:intensity, sampling_node=:camera_pixels,
                chopping=(; unfringed, initial_phase=:fringed))))
    return plant, (llowfs, fringed, unfringed)
end

function test_warmup_reset(plant, graphs)
    state = plant.state
    @test state.domain === nothing && state.last_end_ns == -1
    @test (state.generation, state.sequence, state.lost_samples,
        state.discontinuities, state.rejected) == (0, 0, 0, 0, 0)
    @test all(iszero, plant.pupil_opd) && all(iszero, plant.pupil_amplitude)
    @test all(exposure -> exposure_state(exposure) == (0, 0, 0,
        zeros(size(exposure.image)), zeros(size(exposure.image)), 0, 0), values(plant.exposures))
    @test all(graph -> graph_step_sequence(graph) == 0 && !graph_failed(graph) &&
        !graph_step_pending(graph) && all(iszero, graph_output(graph, :intensity)), graphs)
end

@testset "Fixed external pupil warmup" begin
    reference, reference_graphs = warmup_test_owner()
    warmed, warmed_graphs = warmup_test_owner()
    packet = zeros(Float32, 2, 16, 16)
    opd = view(packet, 1, :, :)
    amplitude = view(packet, 2, :, :)
    fill!(amplitude, 1)
    domain = ntuple(UInt8, 16)
    live_domain = ntuple(i -> UInt8(i + 16), 16)
    expected = NamedTuple[]
    observed = NamedTuple[]
    steps = Ref((0, 0, 0))
    checks = Ref(0)
    record!(records, name, image, metadata) = push!(records,
        (; name, bytes=copy(reinterpret(UInt8, vec(image))), metadata))
    record_reference(name, image, metadata) = record!(expected, name, image, metadata)
    function record_warmup(name, image, metadata)
        steps[] = map(graph_step_sequence, warmed_graphs)
        record!(observed, name, image, metadata)
    end
    function feed!(emit, plant, sequence; domain=domain, generation=UInt64(1))
        step_external_pupil!(emit, plant, opd, amplitude,
            ModelTimestamp(2(sequence - 1)), ModelDuration(2);
            domain, generation, sequence=UInt64(sequence))
    end
    for sequence in 1:10
        @test feed!(record_reference, reference, sequence)
    end
    @test map(graph_step_sequence, reference_graphs) == (10, 5, 5)
    @test warmup_external_pupil!(record_warmup, warmed, opd, amplitude;
        domain, check=() -> (checks[] += 1)) === nothing
    @test checks[] == 10
    @test steps[] == (1, 1, 1)
    @test length(observed) == 12
    @test isequal(observed, expected)
    chopped = filter(record -> record.name === :scc, observed)
    @test [record.metadata.chopped_image_kind for record in chopped] == [:fringed, :unfringed]
    @test chopped[1].bytes != chopped[2].bytes
    test_warmup_reset(warmed, warmed_graphs)

    for cancel_before in (1, 4)
        empty!(observed)
        checks[] = 0
        check = () -> (checks[] += 1; checks[] == cancel_before && throw(InterruptException()))
        @test_throws InterruptException warmup_external_pupil!(record_warmup,
            warmed, opd, amplitude; domain, check)
        @test length(observed) == cancel_before - 1
        test_warmup_reset(warmed, warmed_graphs)
    end
    @test_throws ErrorException warmup_external_pupil!((_, _, _) -> error("emitter failure"),
        warmed, opd, amplitude; domain)
    test_warmup_reset(warmed, warmed_graphs)

    # Ordinary acquisition must recompute both phases for a changing pupil and
    # start with a fresh identity, just like an independently reset plant.
    reset_external_pupil!(reference)
    empty!(expected)
    empty!(observed)
    for sequence in 1:10
        for column in axes(opd, 2), row in axes(opd, 1)
            opd[row, column] = 2.5f-8 * (mod(row + 2column + 3sequence, 13) - 6)
            amplitude[row, column] = 0.84f0 + 0.01f0 * mod(row + column + sequence, 9)
        end
        @test feed!(record_reference, reference, sequence; domain=live_domain, generation=UInt64(7))
        @test feed!(record_warmup, warmed, sequence; domain=live_domain, generation=UInt64(7))
    end
    @test isequal(observed, expected)
    @test steps[] == (10, 5, 5)
    @test all(record -> record.metadata.generation == 7, observed)
    @test_throws ArgumentError warmup_external_pupil!(record_warmup, warmed, opd, amplitude; domain)
    @test warmed.state.generation == 7 && warmed.state.sequence == 10

    # Mutating the source in check invalidates retained optical intensities even
    # when the pupil changes back to its initial value on the following sample.
    reset_external_pupil!(reference)
    reset_external_pupil!(warmed)
    empty!(expected)
    empty!(observed)
    fill!(opd, 0)
    amplitude_for(sequence) = sequence == 3 ? 0.5f0 : 1.0f0
    for sequence in 1:10
        fill!(amplitude, amplitude_for(sequence))
        @test feed!(record_reference, reference, sequence)
    end
    checks[] = 0
    check = () -> fill!(amplitude, amplitude_for(checks[] += 1))
    warmup_external_pupil!(record_warmup, warmed, opd, amplitude; domain, check)
    @test isequal(observed, expected)
    @test steps[] == (9, 4, 5)
    test_warmup_reset(warmed, warmed_graphs)

    fill!(warmed.pupil_amplitude, 1)
    warmup_external_pupil!(record_warmup, warmed, warmed.pupil_opd, warmed.pupil_amplitude; domain)
    @test steps[] == (10, 5, 5)
    test_warmup_reset(warmed, warmed_graphs)
end

@testset "External pupil warmup with stateful optical and response owners" begin
    intensity = Ref(0.0)
    optical = external_intensity_test_graph(intensity)
    response = external_response_graph()
    specification = (graph=response, input=:photon_rate, output=:frame,
        normalized_to_photon_rate=1e9, pixel_sampling_m=(1.0, 1.0),
        exposure_duration_s=4e-9, frame_schema="test.frame/1")
    owner = prepare_external_pupil((sensor=optical,);
        exposure_durations=(sensor=ModelDuration(4),),
        observations=(sensor=(output=:intensity, sampling_node=:optics, response=specification),))
    output = Tuple{Float32, UInt64, UInt64}[]
    function emit(name, image, metadata)
        push!(output, (first(image), graph_step_sequence(optical), graph_step_sequence(response)))
    end
    warmup_external_pupil!(emit, owner, zeros(4, 4), ones(4, 4);
        domain=ntuple(UInt8, 16), check=() -> (intensity[] += 2))
    @test output == [(4f0, 1, 1), (8f0, 2, 2)]
    @test map(graph_step_sequence, (optical, response)) == (0, 0)
    @test all(iszero, graph_output(response, :frame))

    static_optics = external_test_graph()
    pitch = optical_sampling_m(AlgorithmGraphs.prepared_graph_node(static_optics, :optics))
    sampling = (pitch, pitch)
    detector = external_response_graph(; sampling)
    detector_owner = prepare_external_pupil((sensor=static_optics,);
        exposure_durations=(sensor=ModelDuration(4),),
        observations=(sensor=(output=:intensity, sampling_node=:optics,
            response=merge(specification, (graph=detector, pixel_sampling_m=sampling))),))
    steps = Tuple{UInt64, UInt64}[]
    warmup_external_pupil!((_, _, _) ->
            push!(steps, map(graph_step_sequence, (static_optics, detector))),
        detector_owner, zeros(4, 4), ones(4, 4); domain=ntuple(UInt8, 16))
    @test steps == [(1, 1), (2, 2)]
end
