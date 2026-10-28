using Test
using SHA
using YAML
using AdaptiveOpticsSim

const ZCOMP = AdaptiveOpticsSim.Optics.Zemax

function component_zmx(path, text)
    open(path, "w") do io
        write(io, UInt8[0xff, 0xfe])
        for unit in transcode(UInt16, text)
            write(io, UInt8(unit & 0xff), UInt8(unit >> 8))
        end
    end
    return path
end

function component_yaml(path, document)
    open(io -> YAML.write(io, document), path, "w")
    return path
end

component_hash(path) = bytes2hex(sha256(read(path)))
component_reference(path) = Dict("path" => basename(path), "sha256" => component_hash(path))

function component_fixture(root, text; first_surface=1, last_surface=2,
        initial_direction=1, powered_mirrors=[1, 2], catalogs=Dict{String,String}[])
    source = component_zmx(joinpath(root, "source.ZMX"), text)
    instrument = Dict("optical_prescription" => merge(component_reference(source), Dict("format" => "zemax")))
    component_yaml(joinpath(root, "instrument.yaml"), instrument)
    recipe = Dict("schema_version" => 2, "catalogs" => catalogs,
        "conversion" => Dict("first_surface" => first_surface,
            "initial_direction" => initial_direction, "powered_mirrors" => powered_mirrors),
        "branches" => Dict("relay" => Dict("configuration_id" => 1, "last_surface" => last_surface)),
        "assumptions" => ["Synthetic centered relay for verification."])
    component_yaml(joinpath(root, "recipe.yaml"), recipe)
    return (; source, instrument, recipe)
end

const COMPONENT_SETTINGS = Dict("resolution" => 32, "pupil_resolution" => 24,
    "wavelength_um" => 0.6, "temperature_c" => 20.0)

const COMPONENT_RELAY = """
UNIT MM
SURF 0
  TYPE STANDARD
SURF 1
  TYPE STANDARD
  CURV -0.001
  GLAS MIRROR
  DISZ -750
  DIAM 100 0
SURF 2
  TYPE STANDARD
  CURV 0.002
  GLAS MIRROR
  DISZ 999
  CLAP 0 10 0
MNUM 1 1
"""

@testset "Composable centered ZEMAX mirror subsystem" begin
    mktempdir() do root
        fixture = component_fixture(root, COMPONENT_RELAY)
        before = Dict(name => read(joinpath(root, name)) for name in readdir(root))
        prepared = ZCOMP.prepare_component(root, "recipe.yaml", "relay";
            settings=COMPONENT_SETTINGS, diameter_m=0.04)
        @test prepared.config isa ZCOMP.TrainConfiguration
        @test Set(keys(prepared.config)) == Set(("resolution", "pupil_resolution", "wavelength_um", "diameter_m", "steps"))
        steps = prepared.config["steps"]
        @test [s["focal_length_m"] for s in steps if s["op"] == "lens"] == [0.5, 0.25]
        @test [s["distance_m"] for s in steps if s["op"] == "propagate"] == [0.75]
        @test [s["radius_m"] for s in steps if s["op"] == "circular_aperture"] == [0.01]
        # Two concave mirrors separated by f₁ + f₂ form an afocal relay:
        # A = -f₂/f₁, C = 0, D = -f₁/f₂ in unfolded coordinates.
        @test prepared.diagnostics["2"] ≈ [-0.5 0.75; 0.0 -2.0]
        @test Set(keys(prepared.sources)) == Set(("instrument.yaml", "source.ZMX", "recipe.yaml"))
        @test all(component_hash(joinpath(root, name)) == digest for (name, digest) in prepared.sources)
        @test first(prepared.assumptions) == only(fixture.recipe["assumptions"])
        @test all(read(joinpath(root, name)) == bytes for (name, bytes) in before)
        @test Set(readdir(root)) == Set(keys(before))
        @test !ispath(joinpath(root, "generated"))

        fixture.recipe["branches"]["first_mirror"] = Dict("configuration_id" => 1, "last_surface" => 1)
        component_yaml(joinpath(root, "recipe.yaml"), fixture.recipe)
        subset = ZCOMP.prepare_component(root, "recipe.yaml", "first_mirror";
            settings=COMPONENT_SETTINGS, diameter_m=0.04)
        @test length(subset.config["steps"]) == 1
        @test only(subset.config["steps"])["focal_length_m"] == 0.5

        # Reversing the incoming axis and all signed curvatures/thicknesses
        # preserves the unfolded train, including both reflection sign flips.
        reverse = replace(COMPONENT_RELAY, "CURV -0.001" => "CURV 0.001",
            "CURV 0.002" => "CURV -0.002", "DISZ -750" => "DISZ 750")
        component_fixture(root, reverse; initial_direction=-1)
        reversed = ZCOMP.prepare_component(root, "recipe.yaml", "relay";
            settings=COMPONENT_SETTINGS, diameter_m=0.04)
        @test reversed.config == prepared.config
        @test reversed.diagnostics == prepared.diagnostics

        # Convex incidence gives negative optical power rather than focusing.
        component_fixture(root, replace(COMPONENT_RELAY, "CURV -0.001" => "CURV 0.001"))
        convex = ZCOMP.prepare_component(root, "recipe.yaml", "relay";
            settings=COMPONENT_SETTINGS, diameter_m=0.04)
        @test first(convex.config["steps"])["focal_length_m"] == -0.5
    end
end

@testset "Component catalog ordering, refraction and provenance" begin
    mktempdir() do root
        first_catalog = joinpath(root, "first.agf")
        second_catalog = joinpath(root, "second.agf")
        write(first_catalog, "NM G 2 0 1.5 0\nCD 1.25 0 0 0 0 0\nLD 0.3 2.5\n")
        write(second_catalog, "NM G 2 0 2.0 0\nCD 3 0 0 0 0 0\nLD 0.3 2.5\n")
        text = """
        UNIT MM
        SURF 0
          TYPE STANDARD
        SURF 1
          TYPE STANDARD
          CURV 0.001
          GLAS G
          DISZ 10
          DIAM 5 1
        SURF 2
          TYPE STANDARD
          DISZ 999
        MNUM 1 1
        """
        fixture = component_fixture(root, text; powered_mirrors=Int[],
            catalogs=[component_reference(first_catalog), component_reference(second_catalog)])
        prepared = ZCOMP.prepare_component(root, "recipe.yaml", "relay";
            settings=COMPONENT_SETTINGS, diameter_m=0.04)
        steps = prepared.config["steps"]
        @test first(steps)["focal_length_m"] == 2.0 # φ = (1.5 - 1) × 1 m⁻¹
        @test only(s["distance_m"] for s in steps if s["op"] == "propagate") ≈ 0.01 / 1.5
        @test only(s["radius_m"] for s in steps if s["op"] == "circular_aperture") == 0.005
        @test Set(keys(prepared.sources)) == Set(("instrument.yaml", "source.ZMX", "recipe.yaml", "first.agf", "second.agf"))
        write(first_catalog, read(first_catalog, String) * "CC Changed catalog\n")
        @test_throws ArgumentError ZCOMP.prepare_component(root, "recipe.yaml", "relay";
            settings=COMPONENT_SETTINGS, diameter_m=0.04)
        # A subsystem cut after a glass entrance cannot implicitly reset n to 1.
        component_fixture(root, text; first_surface=2, powered_mirrors=Int[],
            catalogs=[component_reference(first_catalog)])
        @test_throws ArgumentError ZCOMP.prepare_component(root, "recipe.yaml", "relay";
            settings=COMPONENT_SETTINGS, diameter_m=0.04)
    end
end

@testset "Component recipe and selected-surface rejection" begin
    mktempdir() do root
        prepare() = ZCOMP.prepare_component(root, "recipe.yaml", "relay";
            settings=COMPONENT_SETTINGS, diameter_m=0.04)
        for text in (replace(COMPONENT_RELAY, "TYPE STANDARD\n  CURV -0.001" => "TYPE EVENASPH\n  CURV -0.001"),
                replace(COMPONENT_RELAY, "GLAS MIRROR" => "GLAS UNKNOWN"),
                replace(COMPONENT_RELAY, "CLAP 0 10 0" => "CLAP 1 10 0"))
            component_fixture(root, text)
            @test_throws ArgumentError prepare()
        end
        fixture = component_fixture(root, COMPONENT_RELAY; powered_mirrors=[1])
        @test_throws ArgumentError prepare()
        fixture = component_fixture(root, COMPONENT_RELAY; powered_mirrors=[1, 1])
        @test_throws ArgumentError prepare()
        fixture = component_fixture(root, COMPONENT_RELAY)
        fixture.recipe["conversion"]["unknown"] = 1
        component_yaml(joinpath(root, "recipe.yaml"), fixture.recipe)
        @test_throws ArgumentError prepare()
        component_fixture(root, COMPONENT_RELAY)
        open(io -> write(io, "schema_version: 2\n"), joinpath(root, "recipe.yaml"), "a")
        @test_throws ArgumentError prepare()
        component_fixture(root, COMPONENT_RELAY)
        @test_throws ArgumentError ZCOMP.prepare_component(root, "recipe.yaml", "relay";
            settings=merge(COMPONENT_SETTINGS, Dict("unknown" => 1)), diameter_m=0.04)
        @test_throws ArgumentError ZCOMP.prepare_component(root, "recipe.yaml", "relay";
            settings=COMPONENT_SETTINGS, diameter_m=0.0)
        @test_throws ArgumentError ZCOMP.prepare_component(root, "recipe.yaml", "relay";
            settings=COMPONENT_SETTINGS, diameter_m=true)
        for override in (Dict("resolution" => true), Dict("resolution" => 31), Dict("wavelength_um" => true),
                Dict("temperature_c" => Inf), Dict("pupil_resolution" => 64))
            @test_throws ArgumentError ZCOMP.prepare_component(root, "recipe.yaml", "relay";
                settings=merge(COMPONENT_SETTINGS, override), diameter_m=0.04)
        end
        @test_throws ArgumentError ZCOMP.prepare_component(root, "../recipe.yaml", "relay";
            settings=COMPONENT_SETTINGS, diameter_m=0.04)
        mkpath(joinpath(root, "generated"))
        cp(joinpath(root, "recipe.yaml"), joinpath(root, "generated", "recipe.yaml"))
        @test_throws ArgumentError ZCOMP.prepare_component(root, "generated/recipe.yaml", "relay";
            settings=COMPONENT_SETTINGS, diameter_m=0.04)
        write(joinpath(root, "source.ZMX"), "changed")
        @test_throws ArgumentError prepare()
    end
end

function authored_component_fixture(root; catalogs=Dict{String,String}[])
    fixture = component_fixture(root, COMPONENT_RELAY; catalogs)
    mkpath(joinpath(root, "simulation"))
    imported = Dict("template" => "proper_optical_train_f32",
        "preparation" => Dict("kind" => "zemax", "recipe" => "recipe.yaml", "branch" => "relay"),
        "config" => Dict("resolution" => Dict("\$ref" => "/settings/resolution")))
    graph = Dict("schema_version" => 1, "nodes" => [Dict("id" => "relay", "component" => "relay")],
        "outputs" => Dict("image" => "relay.image"))
    component_yaml(joinpath(root, "simulation", "production.yaml"), graph)
    definition = Dict("schema_version" => 1, "instrument" => "../instrument.yaml",
        "graphs" => Dict("production" => "simulation/production.yaml", "optical-reference" => graph),
        "components" => Dict("relay" => imported, "relay_copy" => deepcopy(imported)),
        "defaults" => COMPONENT_SETTINGS,
        "profiles" => Dict("test" => Dict("template" => "test")),
        "responses" => Dict("mirror" => Dict("kind" => "second_order")),
        "bindings" => Dict("pupil" => Dict("kind" => "pupil_amplitude")),
        "hil" => Dict("command" => ["mirror.command"], "frame" => "relay.image"))
    component_yaml(joinpath(root, "simulation", "definition.yaml"), definition)
    return (; fixture, definition, graph)
end

@testset "Authored component graphs generated preparation" begin
    mktempdir() do root
        (; definition, graph) = authored_component_fixture(root)
        authored_paths = ["instrument.yaml", "source.ZMX", "recipe.yaml", "simulation/definition.yaml", "simulation/production.yaml"]
        authored = Dict(path => read(joinpath(root, path)) for path in authored_paths)
        in_memory = ZCOMP.prepare_component(root, "recipe.yaml", "relay";
            settings=COMPONENT_SETTINGS, diameter_m=0.04)
        @test !ispath(joinpath(root, "generated"))
        @test_throws ArgumentError ZCOMP.prepare_component(root, "recipe.yaml", "relay";
            settings=COMPONENT_SETTINGS, diameter_m=0.04, frozen=true)
        ZCOMP.prepare(root; yamlfmt="true")
        receipt = ZCOMP.validate_provenance(root)
        @test all(read(joinpath(root, path)) == bytes for (path, bytes) in authored)
        generated = ["generated/production.yaml", "generated/optical-reference.yaml", "generated/simulation.yaml",
            "generated/optics/resolved.yaml", "generated/optics/glasses.yaml", "generated/optics/provenance.yaml"]
        before = Dict(path => read(joinpath(root, path)) for path in generated)
        @test YAML.load_file(joinpath(root, "generated", "production.yaml")) == graph
        @test YAML.load_file(joinpath(root, "generated", "optical-reference.yaml")) == graph
        expected = deepcopy(definition)
        expected["graphs"] = Dict("production" => "generated/production.yaml", "optical-reference" => "generated/optical-reference.yaml")
        expected["optical_model"] = Dict("format" => "zemax", "source" => "../simulation/definition.yaml")
        @test YAML.load_file(joinpath(root, "generated", "simulation.yaml")) == expected
        withenv("YAMLFMT" => "/does/not/exist") do
            ZCOMP.prepare(root; yamlfmt=nothing)
        end
        @test all(read(joinpath(root, path)) == bytes for (path, bytes) in before)
        frozen = ZCOMP.prepare_component(root, "recipe.yaml", "relay";
            settings=COMPONENT_SETTINGS, diameter_m=0.04, frozen=true)
        @test frozen == in_memory
        changed_settings = merge(COMPONENT_SETTINGS, Dict("resolution" => 64, "pupil_resolution" => 48, "wavelength_um" => 0.8, "temperature_c" => 5.0))
        changed = ZCOMP.prepare_component(root, "recipe.yaml", "relay";
            settings=changed_settings, diameter_m=0.08, frozen=true)
        @test changed == ZCOMP.prepare_component(root, "recipe.yaml", "relay";
            settings=changed_settings, diameter_m=0.08)
        @test changed.config != frozen.config
        @test changed.config["resolution"] == 64
        @test changed.config["diameter_m"] == 0.08
        @test all(read(joinpath(root, path)) == bytes for (path, bytes) in before)

        # Updating receipt hashes alone cannot authorize a changed generated graph
        # or a changed execution definition: both must match authored content.
        for (relative, key) in (("generated/production.yaml", "graphs"), ("generated/simulation.yaml", "execution"))
            document = YAML.load_file(joinpath(root, relative))
            document["unexpected"] = "tampered"
            component_yaml(joinpath(root, relative), document)
            @test_throws ArgumentError ZCOMP.validate_provenance(root)
            rewritten = deepcopy(receipt)
            entry = key == "execution" ? rewritten[key] : only(filter(e -> e["path"] == relative, rewritten[key]))
            entry["sha256"] = component_hash(joinpath(root, relative))
            @test_throws ArgumentError ZCOMP.validate_provenance(root, rewritten)
            write(joinpath(root, relative), before[relative])
        end
        path = joinpath(root, "simulation", "production.yaml")
        write(path, read(path, String) * "# source drift\n")
        @test_throws ArgumentError ZCOMP.validate_provenance(root)
        @test_throws ArgumentError ZCOMP.prepare_component(root, "recipe.yaml", "relay";
            settings=COMPONENT_SETTINGS, diameter_m=0.04, frozen=true)
        write(path, authored["simulation/production.yaml"])
        @test ZCOMP.validate_provenance(root) == receipt
    end
end

@testset "Generated component catalog snapshots and output guards" begin
    mktempdir() do root
        catalog = joinpath(root, "glass.agf")
        write(catalog, "NM G 2 0 1.5 0\nCD 1.25 0 0 0 0 0\nLD 0.3 2.5\nTD 0 0 0 0 0 0 20\nNM UNUSED 2 0 1.5 0\nTD\n")
        (; definition) = authored_component_fixture(root; catalogs=[component_reference(catalog)])
        # Two recipes can order duplicate glass names independently while
        # sharing the same prescription and selected configuration.
        second_catalog = joinpath(root, "second.agf")
        write(second_catalog, "NM G 2 0 2.0 0\nCD 3 0 0 0 0 0\nLD 0.3 2.5\n")
        component_fixture(root, """
            UNIT MM
            SURF 0
              TYPE STANDARD
            SURF 1
              TYPE STANDARD
              CURV 0.001
              GLAS G
              DISZ 10
              DIAM 5 1
            SURF 2
              TYPE STANDARD
              DISZ 999
            MNUM 1 1
            """; powered_mirrors=Int[], catalogs=[component_reference(catalog), component_reference(second_catalog)])
        other = YAML.load_file(joinpath(root, "recipe.yaml"))
        reverse!(other["catalogs"])
        component_yaml(joinpath(root, "other.yaml"), other)
        definition["components"]["relay_copy"]["preparation"]["recipe"] = "other.yaml"
        component_yaml(joinpath(root, "simulation", "definition.yaml"), definition)
        ZCOMP.prepare(root; yamlfmt="true")
        snapshot = YAML.load_file(joinpath(root, "generated", "optics", "glasses.yaml"))
        @test Set(keys(snapshot["recipes"])) == Set(("recipe.yaml", "other.yaml"))
        extension = Base.get_extension(AdaptiveOpticsSim, :AdaptiveOpticsSimYAMLExt)
        empty_glass = extension.frozen_glass(snapshot["recipes"]["recipe.yaml"]["UNUSED"])
        @test empty_glass["TD"] isa Vector{Float64}
        @test isempty(empty_glass["TD"])
        # Empty strings and coefficient sequences must retain their declared
        # element types even when YAML supplied Vector{Any}.
        @test extension.frozen_glass(Dict("raw_records" => Any[]))["raw_records"] isa Vector{String}
        for (recipe, expected_focal) in (("recipe.yaml", 2.0), ("other.yaml", 1.0))
            frozen = ZCOMP.prepare_component(root, recipe, "relay"; settings=COMPONENT_SETTINGS,
                diameter_m=0.04, frozen=true)
            @test frozen == ZCOMP.prepare_component(root, recipe, "relay";
                settings=COMPONENT_SETTINGS, diameter_m=0.04)
            @test first(frozen.config["steps"])["focal_length_m"] == expected_focal
        end
        write(catalog, read(catalog, String) * "CC catalog changed\n")
        @test_throws ArgumentError ZCOMP.validate_provenance(root)
    end
    for alias in (:symlink, :hardlink)
        mktempdir() do root
            authored_component_fixture(root)
            mkpath(joinpath(root, "generated"))
            source = joinpath(root, "simulation", "production.yaml")
            output = joinpath(root, "generated", "production.yaml")
            alias == :symlink ? symlink(source, output) : hardlink(source, output)
            before = read(source)
            @test_throws ArgumentError ZCOMP.prepare(root; yamlfmt="true")
            @test read(source) == before
            @test !isfile(joinpath(root, "generated", "simulation.yaml"))
        end
    end
    for change in (d -> (d["graphs"]["production"] = "generated/production.yaml"),
            d -> (d["components"]["relay"]["preparation"]["recipe"] = "generated/recipe.yaml"),
            d -> (d["optical_model"] = Dict("format" => "zemax", "source" => "../generated/simulation.yaml")))
        mktempdir() do root
            (; definition) = authored_component_fixture(root)
            change(definition)
            component_yaml(joinpath(root, "simulation", "definition.yaml"), definition)
            @test_throws ArgumentError ZCOMP.prepare(root; yamlfmt="true")
            @test !ispath(joinpath(root, "generated"))
        end
    end
end

@testset "Runtime OPD placement uses active surface identities" begin
    mktempdir() do root
        component_fixture(root, COMPONENT_RELAY)
        baseline = ZCOMP.prepare_component(root,"recipe.yaml","relay";
            settings=COMPONENT_SETTINGS,diameter_m=0.04)
        screen = Dict("surface"=>1,"input"=>"woofer_opd","size"=>[12,8],"axes"=>[2,1],
            "diameter_m"=>0.04,"schema"=>"test.opd-metres/1","outside_opd_m"=>0)
        last = merge(screen,Dict("surface"=>2,"input"=>"tweeter_opd","diameter_m"=>0.02))
        result = ZCOMP.prepare_component(root,"recipe.yaml","relay";
            settings=COMPONENT_SETTINGS,diameter_m=0.04,phase_screens=[screen,last])
        steps = result.config["steps"]
        @test filter(s->s["op"]!="add_opd",steps)==baseline.config["steps"]
        @test result.diagnostics==baseline.diagnostics
        @test result.sources==baseline.sources
        @test steps[1]["name"]=="surface_1_power"
        @test steps[2]["op"]=="add_opd" && steps[2]["input"]=="woofer_opd"
        @test steps[3]["name"]=="surface_1_propagation"
        @test steps[end]["op"]=="add_opd" && steps[end]["input"]=="tweeter_opd"
        @test steps[2]["size"]==[12,8] && steps[2]["axes"]==[2,1]
        @test steps[2]["sampling_m"]==[0.04/8,0.04/12]
        @test steps[2]["center_pixel"]==[4.5,6.5]
        @test steps[2]["outside_opd_m"]==0.0
        for change in (Dict("surface"=>0),Dict("surface"=>3),Dict("surface"=>1.0),
                Dict("input"=>"pupil_opd"),Dict("input"=>"a.b"),Dict("size"=>[0,8]),
                Dict("size"=>[8]),Dict("size"=>4),Dict("axes"=>[1,1]),Dict("axes"=>4),
                Dict("diameter_m"=>0),Dict("diameter_m"=>true),Dict("diameter_m"=>Inf),
                Dict("schema"=>""),Dict("outside_opd_m"=>NaN),Dict("outside_opd_m"=>true))
            @test_throws ArgumentError ZCOMP.prepare_component(root,"recipe.yaml","relay";
                settings=COMPONENT_SETTINGS,diameter_m=0.04,phase_screens=[merge(screen,change)])
        end
        grid = Dict("size"=>[12,8],"axes"=>[2,1],"diameter_m"=>0.04)
        mapped = ZCOMP.prepare_component(root,"recipe.yaml","relay";
            settings=COMPONENT_SETTINGS,diameter_m=0.04,pupil_grid=grid)
        @test mapped.config["pupil_grid"] == Dict("size"=>[12,8],"axes"=>[2,1],
            "sampling_m"=>[0.04/8,0.04/12],"center_pixel"=>[4.5,6.5])
        @test mapped.config["steps"] == baseline.config["steps"]
        for change in (Dict("size"=>[0,8]),Dict("axes"=>[2,2]),
                Dict("diameter_m"=>NaN),Dict("extra"=>1))
            @test_throws ArgumentError ZCOMP.prepare_component(root,"recipe.yaml","relay";
                settings=COMPONENT_SETTINGS,diameter_m=0.04,pupil_grid=merge(grid,change))
        end
        @test_throws ArgumentError ZCOMP.prepare_component(root,"recipe.yaml","relay";
            settings=COMPONENT_SETTINGS,diameter_m=0.04,phase_screens=[screen,screen])
    end
end
