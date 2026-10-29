# Run after test_zemax.jl; its synthetic prescription has no optional dependencies.
@testset "authored Zemax optical precision" begin
    mktempdir() do root
        write_zemax_fixture(root)
        Z.prepare(root; yamlfmt="true")
        original = read(joinpath(root, "generated/test.yaml"))
        path = joinpath(root, "simulation/definition.yaml")
        definition = YAML.load_file(path)
        definition["defaults"]["precision"] = "Float64"
        definition["bindings"] = Dict(name => Dict("kind" => "zeros", "type" => "Float64",
            "shape" => ["pupil_resolution", "pupil_resolution"])
            for name in ("pupil_opd", "pupil_amplitude"))
        YAML.write_file(path, definition)
        Z.prepare(root; yamlfmt="true")
        @test read(joinpath(root, "generated/test.yaml")) == original
        definition["defaults"]["precision"] = "Float32"
        for binding in values(definition["bindings"])
            binding["type"] = "Float32"
        end
        YAML.write_file(path, definition)
        authored = authored_fixture_bytes(root)
        Z.prepare(root; yamlfmt="true")
        graph = YAML.load_file(joinpath(root, "generated/test.yaml"))
        @test only(graph["nodes"])["type"] == "proper_optical_train_f32"
        expected = YAML.load(String(original))
        only(expected["nodes"])["type"] = "proper_optical_train_f32"
        @test graph == expected # Precision changes neither equations nor physical settings.
        execution = YAML.load_file(joinpath(root, "generated/simulation.yaml"))
        @test all(b["type"] == "Float32" for b in values(execution["bindings"]))
        @test Z.execution_graphs(root)["test"] == graph
        @test bounded_zemax_storage(Z.execution_graphs(root))
        @test Z.validate_provenance(root)["schema_version"] == 1
        @test Z.execution_graphs(root; settings=Dict("precision" => "Float32"))["test"] == graph
        @test_throws ArgumentError Z.execution_graphs(root; settings=Dict("precision" => "Float64"))
        @test_throws ArgumentError Z.execution_graphs(root; settings=Dict("optical_precision" => "Float32"))
        @test authored_fixture_bytes(root) == authored
        generated = Dict(p => read(joinpath(root, p)) for p in
            ("generated/test.yaml", "generated/simulation.yaml", "generated/optics/provenance.yaml"))
        for precision in ("F32", "float32", 32, true, nothing)
            definition["defaults"]["precision"] = precision
            YAML.write_file(path, definition)
            @test_throws ArgumentError Z.prepare(root; yamlfmt="true")
            @test all(read(joinpath(root, p)) == bytes for (p, bytes) in generated)
        end
        definition["defaults"]["precision"] = "Float32"
        definition["bindings"]["pupil_opd"]["type"] = "Float64"
        YAML.write_file(path, definition)
        @test_throws ArgumentError Z.prepare(root; yamlfmt="true")
        @test all(read(joinpath(root, p)) == bytes for (p, bytes) in generated)
        definition["bindings"]["pupil_opd"]["type"] = "Float32"
        definition["preparation"]["variants"] = Dict("other" => Dict("branch" => "test",
            "settings" => Dict("precision" => "Float64")))
        YAML.write_file(path, definition)
        @test_throws ArgumentError Z.prepare(root; yamlfmt="true")
        @test all(read(joinpath(root, p)) == bytes for (p, bytes) in generated)
    end
end
