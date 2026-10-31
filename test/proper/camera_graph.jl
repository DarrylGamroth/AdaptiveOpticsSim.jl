using AdaptiveOpticsSim, Proper, Test
using AdaptiveOpticsSim.AlgorithmGraphs

function proper_camera_test_source()
    n = 16
    return Dict("schema_version" => 1, "name" => "optical_source",
        "nodes" => [Dict("name" => "optics", "type" => "proper_optical_train_f32",
            "config" => Dict("resolution" => n, "pupil_resolution" => n,
                "diameter_m" => 0.016, "wavelength_um" => 1.55,
                "output_schema" => "test.raw-cell-flux/1",
                "steps" => [Dict("op" => "sampled_mask", "name" => "entrance",
                    "amplitude" => "entrance", "sampling_m" => [0.001, 0.001],
                    "center_pixel" => [9.0, 9.0], "outside_amplitude" => 0.0,
                    "outside_phase_rad" => 0.0),
                    Dict("op" => "lens", "focal_length_m" => 0.1),
                    Dict("op" => "propagate", "distance_m" => 0.07)]),
            "bindings" => Dict("entrance" => "mask"))],
        "inputs" => [Dict("name" => "pupil_opd", "destination" => "optics.pupil_opd", "binding" => "opd"),
            Dict("name" => "pupil_amplitude", "destination" => "optics.pupil_amplitude", "binding" => "amplitude")],
        "outputs" => [Dict("name" => "intensity", "source" => "optics.output")])
end

# Independent rectangular-cell area oracle in source-pixel coordinates.
function proper_camera_test_integral(raw, shape, ratio, source_center, destination_center)
    result = zeros(Float64, shape)
    for j in axes(result, 2), i in axes(result, 1)
        ylo = (i - destination_center[1] - 0.5)*ratio[1] + source_center[1]
        yhi = ylo + ratio[1]
        xlo = (j - destination_center[2] - 0.5)*ratio[2] + source_center[2]
        xhi = xlo + ratio[2]
        for col in axes(raw, 2), row in axes(raw, 1)
            dy = max(0.0, min(yhi, row + 0.5) - max(ylo, row - 0.5))
            dx = max(0.0, min(xhi, col + 0.5) - max(xlo, col - 0.5))
            result[i, j] += raw[row, col]*dy*dx
        end
    end
    return result
end

@testset "Cold Proper camera graph composition" begin
    AG = AdaptiveOpticsSim.AlgorithmGraphs
    source = proper_camera_test_source()
    snapshot = deepcopy(source)
    mask = Float32[0.6 + 0.2cos(i + 2j) for i in 1:16, j in 1:16]
    bindings = (; mask, opd=zeros(Float32, 16, 16), amplitude=ones(Float32, 16, 16))
    raw = prepare_algorithm_graph(load_algorithm_graph(source; bindings))
    pitch = Float64(optical_sampling_m(AG.prepared_graph_node(raw, :optics)))
    step_graph!(raw)
    for (shape, ratio, center_pixels) in (((7, 9), (1.0, 1.0), (0.5, -0.5)),
            ((8, 10), (2.0, 0.5), (0.0, 0.0)))
        sampling = pitch .* ratio
        center = sampling .* center_pixels
        graph = AG.compose_proper_camera_graph(source; bindings=(;mask),
            output_shape=collect(shape), output_sampling_m=sampling,
            array_center_m=center, output_schema="test.camera-cell-flux/1", name="camera")
        @test source == snapshot
        @test graph["inputs"] == source["inputs"]
        @test graph["nodes"][1] == source["nodes"][1]
        @test graph["nodes"][1] !== source["nodes"][1]
        @test graph["outputs"] == [Dict("name" => "intensity", "source" => "camera_pixels.output")]
        config = graph["nodes"][2]["config"]
        expected_center = ((shape .+ 1) ./ 2) .- center_pixels
        @test config["input_sampling_m"] == [pitch, pitch]
        @test config["input_center_pixel"] == [9.0, 9.0]
        @test config["output_center_pixel"] == collect(expected_center)
        @test config["input_schema"] == "test.raw-cell-flux/1"
        camera = prepare_algorithm_graph(load_algorithm_graph(graph; bindings))
        step_graph!(camera)
        expected = proper_camera_test_integral(graph_output(raw, :intensity),
            shape, ratio, (9.0, 9.0), expected_center)
        @test graph_output(camera, :intensity) ≈ expected rtol=2e-6
    end
    compose(graph=source; kwargs...) = AG.compose_proper_camera_graph(graph;
        bindings=(;mask), output_shape=(7, 9), output_sampling_m=(pitch, pitch),
        output_schema="test.camera-cell-flux/1", name="camera", kwargs...)
    @test compose()["nodes"][2]["config"]["output_center_pixel"] == [4.0, 5.0]
    @test compose(;pixel_node_name="pixels")["outputs"][1]["source"] == "pixels.output"
    @test_throws AlgorithmGraphError compose(;pixel_node_name="optics")
    @test_throws AlgorithmGraphError AG.compose_proper_camera_graph(source;
        output_shape=(7, 9), output_sampling_m=(pitch, pitch),
        output_schema="test.camera-cell-flux/1", name="camera")
    for mutate in (g -> push!(g["nodes"], deepcopy(g["nodes"][1])),
            g -> (g["nodes"][1]["type"] = "compiled_optical_paths_f32"),
            g -> (g["nodes"][1]["type"] = "proper_optical_train_f64"),
            g -> (g["nodes"][1]["config"]["channels"] = []),
            g -> (g["nodes"][1]["config"]["output_sampling_m"] = pitch),
            g -> (g["nodes"][1]["config"]["output_quantity"] = "electric_field"),
            g -> (g["outputs"] = []),
            g -> push!(g["outputs"], deepcopy(g["outputs"][1])),
            g -> (g["outputs"][1]["source"] = "missing.output"),
            g -> (g["outputs"][1]["source"] = "optics.pupil_opd"),
            g -> (g["outputs"][1]["name"] = "field"),
            g -> (g["outputs"][1]["unknown"] = "invalid"),
            g -> (g["inputs"][1]["destination"] = "missing.pupil_opd"),
            g -> (g["inputs"][1]["destination"] = "optics.output"),
            g -> (g["inputs"][1]["unknown"] = "invalid"),
            g -> pop!(g["inputs"]),
            g -> push!(g["inputs"], deepcopy(g["inputs"][1])),
            g -> (g["links"] = [Dict("source" => "optics.output", "destination" => "optics.pupil_opd")]),
            g -> (g["parameters"] = [Dict("destination" => "optics.pupil_opd", "binding" => "opd")]),
            g -> (g["schema_version"] = 2),
            g -> (g["unknown"] = true))
        graph = deepcopy(source)
        mutate(graph)
        @test_throws AlgorithmGraphError compose(graph)
    end
    for geometry in ((;output_shape=(7,)), (;output_shape=(7, 0)),
            (;output_shape=(true, 9)), (;output_sampling_m=(pitch,)),
            (;output_sampling_m=(pitch, 0)), (;array_center_m=(NaN, 0)),
            (;array_center_m=(0, 0, 0)))
        kwargs = merge((;bindings=(;mask), output_shape=(7, 9),
            output_sampling_m=(pitch, pitch), output_schema="test.camera-cell-flux/1", name="camera"), geometry)
        @test_throws AlgorithmGraphError AG.compose_proper_camera_graph(source; kwargs...)
    end
end
