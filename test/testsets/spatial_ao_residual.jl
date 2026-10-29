using YAML

# This selector exercises the declared ideal spatial model, not a facility fit.
const SpatialAOAtmospheres = AdaptiveOpticsSim.Atmospheres
const SpatialAOGraphs = AdaptiveOpticsSim.AlgorithmGraphs

function spatial_ao_test_keywords(; T=Float32)
    return (; resolution=8, telescope_diameter_m=1.0, r0=0.15,
        reference_wavelength_m=500e-9, fractional_cn2=(0.7, 0.3),
        wind_speed=(5.0, 9.0), wind_direction_deg=(0.0, 90.0),
        altitude=(0.0, 5000.0), layer_ids=(:ground, :high),
        atmosphere_step=0.001, rng_seed=17, T)
end

function spatial_ao_test_graph(rejection; T=Float32, target=Backends.HostComputeDevice(),
        execution=StreamGraphExecution())
    node = spatial_ao_residual_node(:residual;
        spatial_ao_test_keywords(; T)..., rejection,
        residual_opd_schema="test.spatial-ao-residual-opd-m/1")
    definition = algorithm_graph((node,); name=:spatial_ao_residual,
        outputs=(graph_output(:opd, :residual => :residual_opd),))
    return prepare_algorithm_graph(definition; target, execution)
end

@testset "spatial AO rejection admission" begin
    rejection = SpatialAORejection(control_cutoff_cycles_per_m=1, residual_amplitude=0.25)
    @test rejection.control_cutoff_cycles_per_m == 1
    @test rejection.residual_amplitude == 0.25
    for cutoff in (0, -1, Inf, NaN, true), amplitude in (0.25,)
        @test_throws InvalidConfiguration SpatialAORejection(
            control_cutoff_cycles_per_m=cutoff, residual_amplitude=amplitude)
    end
    for amplitude in (-0.1, 1.1, Inf, NaN, true)
        @test_throws InvalidConfiguration SpatialAORejection(
            control_cutoff_cycles_per_m=1, residual_amplitude=amplitude)
    end
    for cutoff in (nextfloat(4.0), 1e100, 1e-100)
        model = SpatialAORejection(control_cutoff_cycles_per_m=cutoff, residual_amplitude=0.25)
        @test_throws InvalidConfiguration spatial_ao_test_graph(model)
    end
    @test spatial_ao_test_graph(SpatialAORejection(
        control_cutoff_cycles_per_m=4, residual_amplitude=0)) isa PreparedAlgorithmGraph
    @test_throws AlgorithmGraphError spatial_ao_residual_node(:invalid_schema;
        spatial_ao_test_keywords()..., rejection,
        residual_opd_schema="residual/1", atmosphere_opd_schema="raw/1")
    @test_throws AlgorithmGraphError spatial_ao_residual_node(:invalid;
        spatial_ao_test_keywords()..., rejection,
        residual_opd_schema="")
end

@testset "full-screen spatial AO Fourier response" begin
    for T in (Float32, Float64)
        graph = spatial_ao_test_graph(SpatialAORejection(
            control_cutoff_cycles_per_m=1, residual_amplitude=0.25); T)
        owner = SpatialAOGraphs.prepared_graph_node(graph, Val(:residual))
        filter = owner.rejection
        n = size(first(filter.screens), 1)
        @test n == 3 * spatial_ao_test_keywords().resolution
        @test filter.screens[1] === owner.atmosphere.layers[1].generator.state.phase_rad
        @test !Base.mightalias(filter.workspace.buffer, filter.screens[1])
        @test filter.transfer[1, 1] == 0
        tolerance = 100eps(T)
        for (kx, ky, amplitude) in ((0, 0, 0.0), (1, 0, 0.25),
                (0, 3, 0.25), (1, 2, 0.25), (4, 1, 1.0), (0, 12, 1.0))
            # The full-screen physical width is 3 m: k=3 is exactly the cutoff.
            wave = T[cos(2pi * (kx * (i - 1) + ky * (j - 1)) / n)
                for i in 1:n, j in 1:n]
            for screen in filter.screens
                copyto!(screen, wave)
            end
            SpatialAOAtmospheres._apply_spatial_ao_rejection!(filter)
            for screen in filter.screens
                @test screen ≈ T(amplitude) .* wave atol=tolerance rtol=tolerance
            end
        end
        # No per-frame normalization: multiplying the input multiplies output.
        wave = T[sin(2pi * (i - 1) / n) + cos(2pi * (4j - 4) / n)
            for i in 1:n, j in 1:n]
        expected = T[0.25sin(2pi * (i - 1) / n) + cos(2pi * (4j - 4) / n)
            for i in 1:n, j in 1:n]
        for screen in filter.screens
            screen .= T(2) .* wave
        end
        SpatialAOAtmospheres._apply_spatial_ao_rejection!(filter)
        @test first(filter.screens) ≈ T(2) .* expected atol=tolerance rtol=tolerance
    end
end

@testset "spatial AO timeline, RNG and reset" begin
    raw_node = multilayer_atmosphere_opd_node(:raw; spatial_ao_test_keywords()...,
        atmosphere_opd_schema="test.raw-atmosphere-opd-m/1")
    raw = prepare_algorithm_graph(algorithm_graph((raw_node,); name=:raw_atmosphere,
        outputs=(graph_output(:opd, :raw => :atmosphere_opd),)))
    raw_owner = SpatialAOGraphs.prepared_graph_node(raw, Val(:raw))
    residual = spatial_ao_test_graph(SpatialAORejection(
        control_cutoff_cycles_per_m=1, residual_amplitude=0.25))
    owner = SpatialAOGraphs.prepared_graph_node(residual, Val(:residual))
    @test raw_owner.rejection === nothing
    @test typeof(owner).name === typeof(raw_owner).name
    @test owner.rng == raw_owner.rng
    screen_snapshot = map(copy, owner.rejection.screens)
    @test any(screen_snapshot[i] != raw_owner.atmosphere.layers[i].generator.state.phase_rad
        for i in eachindex(screen_snapshot))
    for _ in 1:3
        step_graph!(raw)
        @test @inferred(step_graph!(residual)) === residual
        @test owner.rng == raw_owner.rng
        @test epoch_time(current_epoch(owner.atmosphere)) ==
            epoch_time(current_epoch(raw_owner.atmosphere))
        @test epoch_sequence(current_epoch(owner.atmosphere)) ==
            epoch_sequence(current_epoch(raw_owner.atmosphere))
        @test owner.rejection.screens == screen_snapshot
    end
    @test epoch_time(current_epoch(owner.atmosphere)) ≈ 3 * Float32(0.001)
    output = graph_output(residual, Val(:opd))
    @test all(isfinite, output)
    @test !all(iszero, output)
    @test all(iszero, output[.!owner.renderer.pupil])
    saved_third = copy(output)
    fill!(owner.rejection.workspace.buffer, ComplexF32(NaN, NaN))
    step_graph!(residual) # warm exact entrypoint
    @test all(isnan, owner.rejection.workspace.buffer) # no per-step FFT or filter read
    @test all(isfinite, output)
    @test @allocated(step_graph!(residual)) == 0
    reset_graph!(residual)
    @test all(iszero, output)
    @test owner.rejection.screens == screen_snapshot # not cumulatively filtered
    @test epoch_time(current_epoch(owner.atmosphere)) == 0
    reset_graph!(raw)
    @test owner.rng == raw_owner.rng
    for _ in 1:3
        step_graph!(residual)
    end
    @test output == saved_third
    # With a=1 only full-screen DC is removed; raw RNG and OPD scale are retained.
    unity = spatial_ao_test_graph(SpatialAORejection(
        control_cutoff_cycles_per_m=1, residual_amplitude=1))
    reset_graph!(raw)
    step_graph!(raw)
    step_graph!(unity)
    @test graph_output(unity, Val(:opd)) ≈ graph_output(raw, Val(:opd)) rtol=2e-6 atol=1e-12
end

@testset "spatial AO residual YAML" begin
    config = Dict{String,Any}(String(key) => value for (key, value) in
        pairs(Base.structdiff(spatial_ao_test_keywords(), NamedTuple{(:T,)})))
    for (key, value) in config
        value isa Tuple && (config[key] = collect(value))
    end
    config["layer_ids"] = ["ground", "high"]
    config["residual_opd_schema"] = "test.spatial-ao-residual-opd-m/1"
    config["control_cutoff_cycles_per_m"] = 1.0
    config["residual_amplitude"] = 0.25
    document = Dict("schema_version" => 1, "name" => "spatial_ao_residual",
        "nodes" => [Dict("name" => "residual", "type" => "spatial_ao_residual_f32",
            "config" => config)],
        "outputs" => [Dict("name" => "opd", "source" => "residual.residual_opd")])
    parsed = mktemp() do path, io
        YAML.write(io, document)
        close(io)
        prepare_algorithm_graph(load_algorithm_graph(path))
    end
    native = spatial_ao_test_graph(SpatialAORejection(
        control_cutoff_cycles_per_m=1, residual_amplitude=0.25))
    step_graph!(parsed)
    step_graph!(native)
    @test graph_output(parsed, Val(:opd)) == graph_output(native, Val(:opd))
    for key in ("control_cutoff_cycles_per_m", "residual_amplitude", "residual_opd_schema")
        value = pop!(config, key)
        @test_throws AlgorithmGraphError load_algorithm_graph(document)
        config[key] = value
    end
    for (key, value) in (("unknown", 1), ("atmosphere_opd_schema", "raw/1"))
        config[key] = value
        @test_throws AlgorithmGraphError load_algorithm_graph(document)
        delete!(config, key)
    end
    document["nodes"][1]["props"] = Dict("gain" => 0.1)
    @test_throws AlgorithmGraphError load_algorithm_graph(document)
end
