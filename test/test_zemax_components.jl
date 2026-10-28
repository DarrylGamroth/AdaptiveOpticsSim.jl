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
