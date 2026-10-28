using Test
using SHA
using YAML
using AdaptiveOpticsSim
const Z = AdaptiveOpticsSim.Optics.Zemax

function utf16_fixture(path, text)
    open(path, "w") do io
        write(io, UInt8[0xff, 0xfe])
        for unit in transcode(UInt16, text)
            write(io, UInt8(unit & 0xff), UInt8(unit >> 8))
        end
    end
    return path
end

const SYNTHETIC_ZMX = """
UNIT MM
SURF 0
  COMM object
  TYPE STANDARD
  DISZ INFINITY
SURF 1
  TYPE COORDBRK
  CURV 0.02 0
  DISZ -2
  DIAM 4 0
  PARM 1 2
  PPAR 1 0 -1 0 0
  SQAP 3 4 0
SURF 2
  TYPE STANDARD
  DISZ 3
SURF 3
  TYPE STANDARD
  GLAS MIRROR 0
SURF 4
  TYPE STANDARD
MNUM 5 3
MOFF 0 3 "SCC" 0
MOFF 0 3 "" 0
MOFF 0 5 "LLOWFS" 0
IGNM 1 3 1 0 2
IGNM 2 3 0 0 3
IGNM 4 3 1 0 4
THIC 1 5 -12 0 0
CRVT 1 5 0.1 0 0
SDIA 1 5 4.5 0 0
PRAM 1 5 5 0 1
XFIE 1 5 0.25 0 0
"""

@testset "bounded sequential ZMX inspection" begin
    mktempdir() do dir
        path = utf16_fixture(joinpath(dir, "synthetic.ZMX"), SYNTHETIC_ZMX)
        result = Z.inspect(path; configurations=[3, 5])
        @test result["source"]["sha256"] == bytes2hex(sha256(read(path)))
        @test result["source"]["surface_count"] == 5
        @test result["saved_configuration_id"] == 3
        @test Z.inspect(path)["selected_configuration_ids"] == [3]
        a, b = result["configurations"]
        @test a["name"] == "SCC" # Empty trailing MOFF must not erase the name.
        @test a["ignored_surface_ids"] == [1, 2, 4]
        @test [s["id"] for s in a["active_surfaces"]] == [0, 3]
        @test a["active_surfaces"][1]["thickness_m"] == Inf
        s = only(filter(s -> s["id"] == 1, b["active_surfaces"]))
        @test s["thickness_m"] == -0.012
        @test s["curvature_inv_m"] == 100
        @test s["semi_diameter_m"] ≈ 0.0045
        @test s["parameters_native"]["1"] == 5
        @test "PPAR 1 0 -1 0 0" in s["raw_records"]
        @test s["square_aperture_native_record"] == "SQAP 3 4 0"
        @test b["field_overrides_native"][1]["value_native"] == 0.25
        @test_throws ErrorException Z.inspect(path; configurations=[6])
        @test_throws ErrorException Z.inspect(path; configurations=[3, 3])
        utf16_fixture(path, SYNTHETIC_ZMX * "WAVE 1 5 1 0 0\n")
        @test_throws ErrorException Z.inspect(path; configurations=[3, 5])
        utf16_fixture(path, replace(SYNTHETIC_ZMX, "SDIA 1 5 4.5 0" => "SDIA 1 5 4.5 6"))
        @test_throws ErrorException Z.inspect(path; configurations=[3, 5])
        utf16_fixture(path, SYNTHETIC_ZMX * "THIC 1 9 1 0 0\n")
        @test_throws ErrorException Z.inspect(path; configurations=[3, 5])
        utf16_fixture(path, replace(SYNTHETIC_ZMX, "UNIT MM" => "UNIT IN"))
        @test_throws ErrorException Z.inspect(path; configurations=[3, 5])
        write(path, "not UTF-16LE")
        @test_throws ErrorException Z.inspect(path; configurations=[3, 5])
    end
end

# Inspect the actual stored types, including generated child documents. This
# also catches accidental widening caused by heterogeneous constructor inputs.
function bounded_zemax_storage(value)
    return true
end
function bounded_zemax_storage(value::AbstractDict)
    return (valtype(value) === Union{} || all(isconcretetype, Base.uniontypes(valtype(value)))) &&
        all(bounded_zemax_storage, values(value))
end
function bounded_zemax_storage(value::AbstractArray)
    return (eltype(value) === Union{} || all(isconcretetype, Base.uniontypes(eltype(value)))) &&
        all(bounded_zemax_storage, value)
end

@testset "bounded optical conversion documents" begin
    mktempdir() do dir
        path = utf16_fixture(joinpath(dir, "synthetic.ZMX"), SYNTHETIC_ZMX)
        inspected = Z.inspect(path; configurations=[3, 5])
        @test bounded_zemax_storage(inspected)
        branch = last(inspected["configurations"])
        model = Dict("conversion" => Dict(
            "first_surface" => 1, "input_lens_surface" => 1,
            "input_focal_parameter" => "1", "oap_surface" => 14,
            "apodizer_surface" => 0, "focal_mask_surface" => 10))
        selection = Dict("configuration_id" => 5, "detector_surface" => 4, "lyot_surface" => 2)
        settings = Dict("resolution" => 16, "pupil_resolution" => 8,
            "wavelength_um" => 1.55, "temperature_c" => 20)
        masks = Dict("lyot" => Dict("circles" => [
            Dict("radius_m" => 0.001, "x_m" => 0, "y_m" => 0)]))
        convert_train(b) = Z.train(b, "test", selection, model, settings,
            masks, Dict(), 0.012; lyot_mode="transmit")
        graph, diagnostics = convert_train(branch)
        @test_throws UndefKeywordError Z.train(branch, "test", selection, model, settings, masks, Dict(), 0.012)
        @test_throws ArgumentError Z.train(branch, "test", selection, model, settings, masks, Dict(), 0.012; lyot_mode="absorb")
        reflected, _ = Z.train(branch, "test", selection, model, settings, masks, Dict(), 0.012; lyot_mode="reflect")
        @test only(filter(s -> s["op"] == "circular_mask", first(reflected["nodes"])["config"]["steps"]))["mode"] == "reflect"
        @test bounded_zemax_storage(graph)
        @test bounded_zemax_storage(diagnostics)
        @test first(graph["nodes"])["config"]["steps"][1]["focal_length_m"] == 0.005
        # Explicitly selected input-lens idealization is independent of its
        # native surface type. Ordinary unmodeled surface types must reject.
        only(s for s in branch["active_surfaces"] if s["id"] == 1)["type"] = "PARAXIAL"
        @test first(convert_train(branch)) == graph
        flat_tilt = deepcopy(branch)
        tilted_surface = only(s for s in flat_tilt["active_surfaces"] if s["id"] == 2)
        tilted_surface["type"] = "TILTSURF"
        @test first(convert_train(flat_tilt)) == graph
        tilted_surface["curvature_inv_m"] = 10.0
        @test_throws ErrorException convert_train(flat_tilt)
        unsupported = deepcopy(branch)
        only(s for s in unsupported["active_surfaces"] if s["id"] == 2)["type"] = "BICONICX"
        @test_throws ErrorException convert_train(unsupported)
        curved_mirror = deepcopy(branch)
        only(s for s in curved_mirror["active_surfaces"] if s["id"] == 3)["curvature_inv_m"] = 10.0
        @test_throws ErrorException convert_train(curved_mirror)
    end
end

@testset "bounded Zemax catalog and YAML preparation" begin
    extension = Base.get_extension(AdaptiveOpticsSim, :AdaptiveOpticsSimYAMLExt)
    mktempdir() do dir
        archive = joinpath(dir, "source.zar")
        write(archive, "synthetic archive identity")
        write(joinpath(dir, "first.AGF"), """
        NM SAMPLE 2
        CD 1.0 0.01 0.2 0.02 0.3 0.03
        LD 0.4 2.0
        TD 0 0 0 0 0 0 20
        """)
        write(joinpath(dir, "second.AGF"), """
        NM SAMPLE 2
        CD 2.0 0.01 0.2 0.02 0.3 0.03
        NM OTHER 2
        CD 1.2 0.01 0.2 0.02 0.3 0.03
        """)
        catalog = extension.read_catalogs(dir,
            Dict("catalog_order" => ["first.AGF", "second.AGF"]), archive)
        @test bounded_zemax_storage(catalog)
        @test catalog["glasses"]["SAMPLE"]["catalog"] == "first.AGF"
        @test catalog["glasses"]["SAMPLE"]["CD"][1] == 1.0
        @test catalog["glasses"]["OTHER"]["CD"][1] == 1.2
        @test catalog["archive"]["sha256"] == bytes2hex(sha256(read(archive)))
        ordered = extension.ordered(catalog)
        @test bounded_zemax_storage(ordered)
        serialized = YAML.write(ordered)
        @test YAML.load(serialized) == catalog
        @test startswith(serialized, "archive:")
        @test bounded_zemax_storage(extension.ordered(Dict{String,String}()))
        # Dynamically parsed YAML is an external boundary; ordering its values
        # must still produce concrete stored types before handing it to YAML.write.
        parsed = YAML.load(serialized)
        @test bounded_zemax_storage(extension.ordered(parsed))
    end
end

function zemax_model_fixture()
    return Dict("schema_version" => 1,
        "catalog_order" => ["test.AGF"], "assumptions" => ["Synthetic centered train."],
        "conversion" => Dict("first_surface" => 1, "input_lens_surface" => 1,
            "input_focal_parameter" => "1", "focal_mask_surface" => 2,
            "oap_surface" => 3, "apodizer_surface" => 1),
        "branches" => Dict("test" => Dict{String,Any}("configuration_id" => 5,
            "detector_surface" => 4, "lyot_surface" => 3)))
end

function write_zemax_fixture(root)
    optics = joinpath(root, "simulation")
    mkpath(optics)
    mkpath(joinpath(root, "generated", "optics"))
    source = utf16_fixture(joinpath(root, "source.ZMX"), SYNTHETIC_ZMX)
    write(joinpath(root, "source.zar"), "synthetic archive identity")
    YAML.write_file(joinpath(optics, "zemax.yaml"), zemax_model_fixture())
    YAML.write_file(joinpath(root, "generated", "optics", "resolved.yaml"), Z.inspect(source; configurations=[5]))
    YAML.write_file(joinpath(root, "generated", "optics", "glasses.yaml"), Dict("schema_version" => 1,
        "archive" => Dict("sha256" => bytes2hex(sha256(read(joinpath(root, "source.zar"))))),
        "catalogs" => [Dict("basename" => "test.AGF", "sha256" => repeat("0", 64))],
        "glasses" => Dict{String,String}()))
    YAML.write_file(joinpath(root, "instrument.yaml"), Dict("id" => "synthetic",
        "optical_input" => Dict("diameter_m" => 0.012),
        "wavefront_sensors" => Dict("camera" => Dict("kind" => "test")),
        "optical_prescription" => Dict("format" => "zemax", "path" => "source.ZMX",
            "sha256" => bytes2hex(sha256(read(source))),
            "catalog_archive" => Dict("path" => "source.zar",
                "sha256" => bytes2hex(sha256(read(joinpath(root, "source.zar")))))),
        "optical_paths" => Dict("test" => Dict("sensor" => "camera", "lyot_stop_port" => "transmitted"))))
    YAML.write_file(joinpath(root, "simulation", "definition.yaml"), Dict(
        "instrument" => "../instrument.yaml",
        "preparation" => Dict("kind" => "zemax", "recipe" => "zemax.yaml"),
        "defaults" => Dict("resolution" => 16, "pupil_resolution" => 8,
            "wavelength_um" => 1.55, "temperature_c" => 20,
            "masks" => Dict("focal_plane" => Dict("radius_m" => 0.001),
                "lyot" => Dict("circles" => [Dict("radius_m" => 0.001, "x_m" => 0.0, "y_m" => 0.0)])))))
    return optics
end

@testset "strict authored ZEMAX conversion recipe schema" begin
    mktempdir() do root
        path = joinpath(root, "zemax.yaml")
        model = zemax_model_fixture()
        YAML.write_file(path, model)
        @test Z.load_model(path) == model
        @test bounded_zemax_storage(Z.load_model(path))
        for change in (
            d -> (d["unexpected"] = 1),
            d -> (d["schema_version"] = true),
            d -> (d["schema_version"] = 2),
            d -> (d["conversion"]["lyot_surfaces"] = Dict("test" => 2)),
            d -> (d["conversion"]["first_surface"] = -1),
            d -> (d["conversion"]["input_lens_surface"] = 5),
            d -> (d["conversion"]["input_focal_parameter"] = 1),
            d -> (d["branches"]["test"]["lyot_mode"] = "transmit"),
            d -> (d["source_archive"] = "../../source.zar"),
            d -> (d["branches"]["test"]["configuration_id"] = false),
            d -> (d["branches"]["test"]["lyot_surface"] = 5),
            d -> (d["branches"]["test"]["lyot_surface"] = 2),
            d -> (d["conversion"]["apodizer_surface"] = 2),
            d -> (d["assumptions"] = fill("assumption", 257)),
            d -> (d["resolved_source"] = repeat("a", 4097)),
            d -> (d["branches"] = Dict("branch$i" => deepcopy(d["branches"]["test"]) for i in 1:17)),
            d -> (d["branches"]["test"]["extra"] = 1),
            d -> delete!(d["branches"]["test"], "lyot_surface"),
            d -> (d["glass_source"] = "/tmp/outside.yaml"),
            d -> (d["catalog_order"] = ["../test.AGF"]),
            d -> (d["catalog_order"] = ["test.AGF", "test.AGF"]),
            d -> (d["assumptions"] = [1]),
        )
            invalid = deepcopy(model)
            change(invalid)
            YAML.write_file(path, invalid)
            @test_throws ArgumentError Z.load_model(path)
        end
        text = YAML.write(model)
        write(path, text * "schema_version: 1\n")
        @test_throws ArgumentError Z.load_model(path)
        write(path, replace(text, "lyot_surface: 3" => "lyot_surface: 3\n    lyot_surface: 3"))
        @test_throws ArgumentError Z.load_model(path)
        # Merge keys must not disable explicit duplicate detection.
        write(path, "<<: {schema_version: 1}\n" * text * "schema_version: 1\n")
        @test_throws ArgumentError Z.load_model(path)
        write(path, text * "\n---\n" * text)
        @test_throws ArgumentError Z.load_model(path)
        write(path, "1: value\n")
        @test_throws ArgumentError Z.load_model(path)
    end
end

@testset "invalid ideal mask geometry cannot replace generated graphs" begin
    for change in (
        masks -> (masks["focal_plane"]["radius_m"] = -0.001),
        masks -> (masks["lyot"]["circles"][1]["radius_m"] = -0.001),
        masks -> (masks["lyot"]["circles"][1]["x_m"] = Inf),
    )
        mktempdir() do root
            write_zemax_fixture(root)
            path=joinpath(root,"simulation","definition.yaml")
            definition=YAML.load_file(path)
            change(definition["defaults"]["masks"])
            YAML.write_file(path,definition)
            output=joinpath(root,"generated","test.yaml")
            write(output,"preserve existing graph")
            @test_throws ArgumentError Z.prepare(root;yamlfmt="true")
            @test read(output,String)=="preserve existing graph"
        end
    end
end

function authored_fixture_bytes(root)
    paths = ("instrument.yaml", "simulation/definition.yaml", "simulation/zemax.yaml",
        "source.ZMX", "source.zar")
    return Dict(path => read(joinpath(root, path)) for path in paths)
end

@testset "Zemax spectra replan without changing source or snapshots" begin
    mktempdir() do root
        write_zemax_fixture(root)
        source = joinpath(root,"simulation","definition.yaml")
        definition = YAML.load_file(source)
        delete!(definition["defaults"],"wavelength_um")
        definition["defaults"]["spectra"] = Dict("test"=>Dict("wavelength_um"=>1.55,
            "fractional_bandwidth"=>0.02,"samples"=>3,"filter"=>"attenuator"))
        YAML.write_file(source,definition)
        instrument = YAML.load_file(joinpath(root,"instrument.yaml"))
        instrument["optical_filters"] = Dict("attenuator"=>Dict(
            "wavelength_m"=>[1.0e-6,2.0e-6],"transmission"=>[0.25,0.25]))
        instrument["optical_paths"]["test"]["filters"] = ["attenuator"]
        YAML.write_file(joinpath(root,"instrument.yaml"),instrument)
        Z.prepare(root;yamlfmt="true")
        paths = ("instrument.yaml","simulation/definition.yaml","simulation/zemax.yaml",
            "generated/simulation.yaml","generated/test.yaml","generated/optics/provenance.yaml")
        original = Dict(path=>read(joinpath(root,path)) for path in paths)
        graphs = Z.execution_graphs(root)
        channels = only(graphs["test"]["nodes"])["config"]["channels"]
        @test bounded_zemax_storage(graphs)
        @test length(channels) == 3
        @test [c["wavelength_um"] for c in channels] ≈ [1.55*(1-0.02/3),1.55,1.55*(1+0.02/3)]
        @test sum(c["weight"] for c in channels) ≈ 1
        @test all(c->c["transmission"]==0.25,channels)
        overrides = Dict("spectra"=>Dict("test"=>Dict("wavelength_um"=>1.7,
            "fractional_bandwidth"=>0.0,"samples"=>1)))
        changed = Z.execution_graphs(root;settings=overrides)
        selected = only(only(changed["test"]["nodes"])["config"]["channels"])
        @test selected["wavelength_um"] == 1.7
        @test selected["transmission"] == 0.25
        @test Dict(path=>read(joinpath(root,path)) for path in paths) == original
        @test_throws ArgumentError Z.execution_graphs(root;settings=Dict("unknown"=>1))
        @test_throws ArgumentError Z.execution_graphs(root;settings=Dict("spectra"=>Dict("test"=>Dict("filter"=>"none"))))
        @test_throws ArgumentError Z.execution_graphs(root;settings=Dict("spectra"=>Dict("test"=>Dict("samples"=>1))))
        outside = Z.execution_graphs(root;settings=Dict("spectra"=>Dict("test"=>Dict(
            "wavelength_um"=>2.5,"fractional_bandwidth"=>0.0,"samples"=>1))))
        @test only(only(outside["test"]["nodes"])["config"]["channels"])["transmission"] == 0
    end
end

function change_physical_fixture(root, change)
    path = joinpath(root, "instrument.yaml")
    instrument = YAML.load_file(path)
    change(instrument)
    YAML.write_file(path, instrument)
end

@testset "Authored observation graphs and calibration identities are prepared together" begin
    mktempdir() do root
        write_zemax_fixture(root)
        definition_path = joinpath(root,"simulation","definition.yaml")
        definition = YAML.load_file(definition_path)
        # This unit boundary hashes source identities. RTC separately validates
        # graph contracts and decodes the FITS binding in integration tests.
        response_path = "simulation/detector.yaml"
        write(joinpath(root,response_path),"schema_version: 1\nname: authored_response\n")
        calibration_path = "transfer.fits"
        write(joinpath(root,calibration_path),"calibration identity bytes")
        definition["observations"] = Dict("test"=>Dict("graph"=>response_path,
            "detector"=>"detector","normalized_to_photon_rate"=>1e6))
        definition["bindings"] = Dict("transfer"=>Dict("kind"=>"fits","path"=>calibration_path,
            "sha256"=>bytes2hex(sha256(read(joinpath(root,calibration_path))))))
        YAML.write_file(definition_path,definition)
        Z.prepare(root;yamlfmt="true")
        execution = YAML.load_file(joinpath(root,"generated","simulation.yaml"))
        @test execution["observations"]["test"]["graph"]=="test_response"
        @test execution["graphs"]["test_response"]==response_path
        receipt = Z.validate_provenance(root)
        @test response_path in [entry["path"] for entry in receipt["sources"]]
        @test calibration_path in [entry["path"] for entry in receipt["sources"]]
        @test length(receipt["graphs"])==1
        write(joinpath(root,response_path),"changed authored graph")
        @test_throws ArgumentError Z.validate_provenance(root)
        # A changed calibration rejects before it can replace any output.
        output = read(joinpath(root,"generated","test.yaml"))
        write(joinpath(root,calibration_path),"stale calibration bytes")
        @test_throws ArgumentError Z.prepare(root;yamlfmt="true")
        @test read(joinpath(root,"generated","test.yaml"))==output
    end
end

@testset "physical prescription preparation and exact provenance" begin
    extension = Base.get_extension(AdaptiveOpticsSim, :AdaptiveOpticsSimYAMLExt)
    mktempdir() do root
        optics = write_zemax_fixture(root)
        authored = authored_fixture_bytes(root)
        @test Z.prepare(root; yamlfmt="true") === nothing
        @test authored_fixture_bytes(root) == authored
        @test readdir(optics) == ["definition.yaml", "zemax.yaml"]
        @test readdir(joinpath(root, "simulation")) == ["definition.yaml", "zemax.yaml"]
        @test isfile(joinpath(root, "generated", "optics", "provenance.yaml"))
        receipt = Z.validate_provenance(root)
        @test Set(e["path"] for e in receipt["sources"]) == Set([
            "instrument.yaml", "simulation/definition.yaml", "simulation/zemax.yaml",
            "generated/optics/resolved.yaml", "generated/optics/glasses.yaml", "source.ZMX", "source.zar"])
        execution = YAML.load_file(joinpath(root, "generated", "simulation.yaml"))
        expected = YAML.load_file(joinpath(root, "simulation", "definition.yaml"))
        delete!(expected, "preparation")
        expected["graphs"] = Dict("test" => "generated/test.yaml")
        expected["instrument"] = "../instrument.yaml"
        expected["optical_model"] = Dict("format"=>"zemax", "source"=>"../simulation/definition.yaml")
        @test execution == expected
        @test receipt["execution"]["path"] == "generated/simulation.yaml"
        @test !occursin("generated/", read(joinpath(root, "simulation", "definition.yaml"), String))
        @test !occursin("generated/", read(joinpath(optics, "zemax.yaml"), String))
        @test bounded_zemax_storage(extension.preparation_inputs(root).execution)
        original_graph, _ = Z.train(only(Z.inspect(joinpath(root, "source.ZMX"); configurations=[5])["configurations"]),
            "test", zemax_model_fixture()["branches"]["test"], zemax_model_fixture(), expected["defaults"],
            expected["defaults"]["masks"], Dict{String,String}(), 0.012; lyot_mode="transmit", instrument_id="synthetic")
        graph = YAML.load_file(joinpath(root, "generated", "test.yaml"))
        @test graph == original_graph
        @test first(graph["nodes"])["config"]["steps"][1]["name"] == "ideal_unity_apodizer"
        steps = first(graph["nodes"])["config"]["steps"]
        @test count(s -> s["op"] == "circular_obscuration", steps) == 1
        @test count(s -> s["op"] == "circular_mask", steps) == 1
        @test only(filter(s -> s["op"] == "circular_mask", steps))["mode"] == "transmit"
        @test startswith(read(joinpath(root, "generated", "test.yaml"), String), "# Generated")
        # Physical routing is the sole authority for the mask mode.
        change_physical_fixture(root, d -> (d["optical_paths"]["test"]["lyot_stop_port"] = "reflected"))
        Z.prepare(root; yamlfmt="true")
        reflected = YAML.load_file(joinpath(root, "generated", "test.yaml"))
        reflected_steps = first(reflected["nodes"])["config"]["steps"]
        @test only(filter(s -> s["op"] == "circular_mask", reflected_steps))["mode"] == "reflect"
        only(filter(s -> s["op"] == "circular_mask", reflected_steps))["mode"] = "transmit"
        @test reflected == graph
        receipt = Z.validate_provenance(root)
        # Additional hash-valid entries cannot expand the source contract.
        write(joinpath(root, "extra.txt"), "extra")
        extra = deepcopy(receipt)
        push!(extra["sources"], Dict("path" => "extra.txt", "sha256" => extension.hashfile(joinpath(root, "extra.txt"))))
        @test_throws ArgumentError Z.validate_provenance(root, extra)
        missing = deepcopy(receipt)
        pop!(missing["sources"])
        @test_throws ArgumentError Z.validate_provenance(root, missing)
        missing_execution = deepcopy(receipt)
        delete!(missing_execution, "execution")
        @test_throws ArgumentError Z.validate_provenance(root, missing_execution)
        wrong_execution = deepcopy(receipt)
        wrong_execution["execution"]["path"] = "generated/test.yaml"
        @test_throws ArgumentError Z.validate_provenance(root, wrong_execution)
        execution_path = joinpath(root, "generated", "simulation.yaml")
        execution_bytes = read(execution_path)
        write(execution_path, read(execution_path, String) * "# drift\n")
        @test_throws ArgumentError Z.validate_provenance(root)
        write(execution_path, execution_bytes)
        # Validation reads frozen artifacts and never invokes a formatter.
        frozen = Dict(path => read(joinpath(root, path)) for path in
            ("generated/simulation.yaml", "generated/test.yaml", "generated/optics/resolved.yaml", "generated/optics/glasses.yaml", "generated/optics/provenance.yaml"))
        withenv("YAMLFMT" => "/does/not/exist") do
            @test Z.validate_provenance(root) == receipt
        end
        @test all(read(joinpath(root, path)) == bytes for (path, bytes) in frozen)
        duplicate = deepcopy(receipt)
        push!(duplicate["graphs"], first(duplicate["graphs"]))
        @test_throws ArgumentError Z.validate_provenance(root, duplicate)
        write(joinpath(root, "generated", "test.yaml"), "stale")
        @test_throws ArgumentError Z.validate_provenance(root)
    end
    # Each rejection occurs before replacing even an already-existing graph.
    physical_changes = (
        d -> (d["optical_prescription"]["format"] = "other"),
        d -> (d["optical_prescription"]["extra"] = "unknown"),
        d -> delete!(d["optical_prescription"], "catalog_archive"),
        d -> (d["optical_prescription"]["sha256"] = repeat("0", 64)),
        d -> (d["optical_prescription"]["sha256"] = repeat("A", 64)),
        d -> (d["optical_prescription"]["path"] = "../outside.ZMX"),
        d -> (d["optical_prescription"]["path"] = "/outside.ZMX"),
        d -> (d["optical_prescription"]["path"] = "./source.ZMX"),
        d -> (d["optical_prescription"]["path"] = "design/../source.ZMX"),
        d -> (d["optical_prescription"]["catalog_archive"]["sha256"] = repeat("0", 64)),
        d -> (d["optical_prescription"]["catalog_archive"]["path"] = "../outside.zar"),
        d -> (d["optical_prescription"]["catalog_archive"]["extra"] = "unknown"),
        d -> (d["optical_paths"]["test"]["sensor"] = "unknown"),
        d -> (d["optical_paths"]["duplicate"] = deepcopy(d["optical_paths"]["test"])),
        d -> (d["optical_paths"]["test"]["sensor"] = ["camera"]),
        d -> (d["optical_paths"]["test"]["lyot_stop_port"] = "transmit"),
        d -> (d["optical_paths"]["test"]["lyot_mode"] = "transmit"),
        d -> (d["optical_paths"]["other"] = pop!(d["optical_paths"], "test")),
    )
    source_changes = (
        d -> (d["graphs"] = Dict("test" => "generated/test.yaml")),
        d -> delete!(d, "preparation"),
        d -> (d["preparation"]["kind"] = "unknown"),
        d -> (d["preparation"]["extra"] = 1),
        d -> delete!(d["preparation"], "recipe"),
        d -> (d["preparation"]["recipe"] = "../../outside.yaml"),
        d -> (d["preparation"]["recipe"] = "../generated/optics/resolved.yaml"),
        d -> (d["preparation"]["recipe"] = "/tmp/zemax.yaml"),
        d -> (d["instrument"] = "zemax.yaml"),
    )
    changes = (
        ((root, optics) -> change_physical_fixture(root, change) for change in physical_changes)...,
        ((root, optics) -> begin
            path = joinpath(root, "simulation", "definition.yaml")
            definition = YAML.load_file(path); change(definition)
            YAML.write_file(path, definition)
        end for change in source_changes)...,
        (root, optics) -> begin
            model = zemax_model_fixture(); model["branches"]["other"] = pop!(model["branches"], "test")
            YAML.write_file(joinpath(optics, "zemax.yaml"), model)
        end,
        ((root, optics) -> begin
            model = zemax_model_fixture(); model["branches"]["test"]["configuration_id"] = id
            YAML.write_file(joinpath(optics, "zemax.yaml"), model)
        end for id in (3, 6))...,
        (root, optics) -> begin
            path = joinpath(root, "generated", "optics", "glasses.yaml"); catalog = YAML.load_file(path)
            catalog["archive"]["sha256"] = repeat("0", 64)
            YAML.write_file(path, catalog)
        end,
        (root, optics) -> begin
            path = joinpath(root, "generated", "optics", "glasses.yaml"); catalog = YAML.load_file(path)
            catalog["schema_version"] = 2
            YAML.write_file(path, catalog)
        end,
        ((root, optics) -> begin
            mv(joinpath(root, source), joinpath(root, "generated", source))
            change_physical_fixture(root, d -> begin
                reference = source == "source.ZMX" ? d["optical_prescription"] : d["optical_prescription"]["catalog_archive"]
                reference["path"] = "generated/$source"
            end)
        end for source in ("source.ZMX", "source.zar"))...,
        (root, optics) -> begin
            path = joinpath(root, "instrument.yaml")
            write(path, read(path, String) * "id: repeated\n")
        end,
    )
    for change in changes
        mktempdir() do root
            optics = write_zemax_fixture(root)
            graph_path = joinpath(root, "generated", "test.yaml")
            write(graph_path, "preserve existing graph")
            change(root, optics)
            snapshots = Dict(filename => read(joinpath(root, "generated", "optics", filename))
                for filename in ("resolved.yaml", "glasses.yaml"))
            @test_throws Union{ArgumentError,ErrorException} Z.prepare(root; yamlfmt="true")
            @test all(read(joinpath(root, "generated", "optics", filename)) == bytes for (filename, bytes) in snapshots)
            @test read(graph_path, String) == "preserve existing graph"
            @test !isfile(joinpath(root, "generated", "optics", "provenance.yaml"))
        end
    end
    mktempdir() do root
        optics = write_zemax_fixture(root)
        original = read(joinpath(root, "generated", "optics", "glasses.yaml"))
        write(joinpath(root, "test.AGF"), "NM UNUSED 2\nCD 1 0.01 0.2 0.02 0.3 0.03\n")
        model = zemax_model_fixture()
        model["branches"]["test"]["configuration_id"] = 3
        YAML.write_file(joinpath(optics, "zemax.yaml"), model)
        @test_throws Union{ArgumentError,ErrorException} Z.prepare(root; catalog_dir=root, yamlfmt="true")
        @test read(joinpath(root, "generated", "optics", "glasses.yaml")) == original
        model["branches"]["test"]["configuration_id"] = 5
        YAML.write_file(joinpath(optics, "zemax.yaml"), model)
        authored = authored_fixture_bytes(root)
        @test Z.prepare(root; catalog_dir=root, yamlfmt="true") === nothing
        @test authored_fixture_bytes(root) == authored
        @test readdir(optics) == ["definition.yaml", "zemax.yaml"]
        @test readdir(joinpath(root, "simulation")) == ["definition.yaml", "zemax.yaml"]
        @test startswith(read(joinpath(root, "generated", "optics", "glasses.yaml"), String), "# Generated")
        @test isempty(YAML.load_file(joinpath(root, "generated", "optics", "glasses.yaml"))["glasses"])
        @test Z.validate_provenance(root)["schema_version"] == 1
    end
    mktempdir() do root
        write_zemax_fixture(root)
        mkpath(joinpath(root, "design"))
        for source in ("source.ZMX", "source.zar")
            mv(joinpath(root, source), joinpath(root, "design", source))
        end
        change_physical_fixture(root, d -> begin
            d["optical_prescription"]["path"] = "design/source.ZMX"
            d["optical_prescription"]["catalog_archive"]["path"] = "design/source.zar"
        end)
        @test Z.prepare(root; yamlfmt="true") === nothing
        sources = Set(e["path"] for e in Z.validate_provenance(root)["sources"])
        @test "design/source.ZMX" in sources && "design/source.zar" in sources
        @test "source.ZMX" ∉ sources && "source.zar" ∉ sources
    end
    for source in ("source.ZMX", "source.zar")
        mktempdir() do outside
            mktempdir() do root
                write_zemax_fixture(root)
                mv(joinpath(root, source), joinpath(outside, source))
                symlink(joinpath(outside, source), joinpath(root, source))
                @test_throws ArgumentError Z.prepare(root; yamlfmt="true")
                @test !isfile(joinpath(root, "generated", "test.yaml"))
            end
        end
    end
    # A lexical generated path cannot target authored or supplied data.
    for filename in ("resolved.yaml", "glasses.yaml")
        mktempdir() do root
            optics = write_zemax_fixture(root)
            asset = joinpath(root, "generated", "optics", filename)
            target = joinpath(optics, filename)
            mv(asset, target)
            symlink(target, asset)
            authored = authored_fixture_bytes(root)
            original = read(target)
            @test_throws ArgumentError Z.prepare(root; yamlfmt="true")
            @test read(target) == original
            @test authored_fixture_bytes(root) == authored
            @test !isfile(joinpath(root, "generated", "test.yaml"))
        end
    end
    for relative in ("generated/test.yaml", "generated/simulation.yaml", "generated/optics/provenance.yaml")
        mktempdir() do root
            write_zemax_fixture(root)
            target = joinpath(root, "authored.yaml")
            write(target, "preserve authored output target")
            symlink(target, joinpath(root, relative))
            authored = authored_fixture_bytes(root)
            @test_throws ArgumentError Z.prepare(root; yamlfmt="true")
            @test read(target, String) == "preserve authored output target"
            @test authored_fixture_bytes(root) == authored
            if relative != "generated/test.yaml"
                @test !isfile(joinpath(root, "generated", "test.yaml"))
            end
        end
    end
    # Nor can an authored or supplied path resolve into generated.
    for relative in ("instrument.yaml", "simulation/definition.yaml", "simulation/zemax.yaml", "source.ZMX", "source.zar")
        mktempdir() do root
            write_zemax_fixture(root)
            path = joinpath(root, relative)
            target = joinpath(root, "generated", basename(relative))
            mv(path, target)
            symlink(target, path)
            original = read(target)
            @test_throws ArgumentError Z.prepare(root; yamlfmt="true")
            @test read(target) == original
            @test !isfile(joinpath(root, "generated", "test.yaml"))
        end
    end
    for destination in (".", "simulation", "supplied")
        mktempdir() do root
            write_zemax_fixture(root)
            mv(joinpath(root, "generated"), joinpath(root, "supplied"))
            target = normpath(joinpath(root, destination))
            symlink(target, joinpath(root, "generated"))
            graph = joinpath(target, "test.yaml")
            write(graph, "preserve aliased data")
            authored = authored_fixture_bytes(root)
            @test_throws ArgumentError Z.prepare(root; yamlfmt="true")
            @test read(graph, String) == "preserve aliased data"
            @test authored_fixture_bytes(root) == authored
            @test !isfile(joinpath(target, "optics", "provenance.yaml"))
        end
    end
    mktempdir() do outside
        mktempdir() do root
            write_zemax_fixture(root)
            mv(joinpath(root, "generated"), joinpath(outside, "generated"))
            symlink(joinpath(outside, "generated"), joinpath(root, "generated"))
            @test_throws ArgumentError Z.prepare(root; yamlfmt="true")
            @test !isfile(joinpath(outside, "generated", "test.yaml"))
        end
    end
    mktempdir() do outside
        mktempdir() do root
            optics = write_zemax_fixture(root)
            asset = joinpath(root, "generated", "optics", "resolved.yaml")
            mv(asset, joinpath(outside, "resolved.yaml"))
            symlink(joinpath(outside, "resolved.yaml"), asset)
            @test_throws ArgumentError Z.prepare(root; yamlfmt="true")
            @test !isfile(joinpath(root, "generated", "test.yaml"))
        end
    end
end

@testset "fresh authored prescription preparation" begin
    mktempdir() do root
        write_zemax_fixture(root)
        rm(joinpath(root, "generated"); recursive=true)
        write(joinpath(root, "test.AGF"), "NM UNUSED 2\nCD 1 0.01 0.2 0.02 0.3 0.03\n")
        @test !ispath(joinpath(root, "generated"))
        @test_throws ArgumentError Z.validate_provenance(root)
        @test_throws ArgumentError Z.prepare(root; yamlfmt="true")
        @test !ispath(joinpath(root, "generated"))
        @test Z.prepare(root; catalog_dir=root, yamlfmt="true") === nothing
        for relative in ("simulation.yaml", "test.yaml", "optics/resolved.yaml", "optics/glasses.yaml", "optics/provenance.yaml")
            @test isfile(joinpath(root, "generated", relative))
        end
        @test Z.validate_provenance(root)["schema_version"] == 1
        # The declared recipe path is configurable and relative to the source definition.
        mv(joinpath(root, "simulation", "zemax.yaml"), joinpath(root, "recipe.yaml"))
        path = joinpath(root, "simulation", "definition.yaml")
        definition = YAML.load_file(path); definition["preparation"]["recipe"] = "../recipe.yaml"
        YAML.write_file(path, definition)
        @test Z.prepare(root; yamlfmt="true") === nothing
        @test "recipe.yaml" in [entry["path"] for entry in Z.validate_provenance(root)["sources"]]
    end
end

@testset "fresh invalid recipe and output aliases have no side effects" begin
    for change in (
        model -> (model["unexpected"] = 1),
        model -> (model["branches"]["simulation"] = pop!(model["branches"], "test")),
        model -> (model["branches"]["../test"] = pop!(model["branches"], "test")),
    )
        mktempdir() do root
            optics = write_zemax_fixture(root)
            rm(joinpath(root, "generated"); recursive=true)
            write(joinpath(root, "test.AGF"), "NM UNUSED 2\n")
            model = zemax_model_fixture(); change(model)
            YAML.write_file(joinpath(optics, "zemax.yaml"), model)
            if haskey(model["branches"], "simulation")
                change_physical_fixture(root, d -> (d["optical_paths"]["simulation"] = pop!(d["optical_paths"], "test")))
            end
            @test_throws ArgumentError Z.prepare(root; catalog_dir=root, yamlfmt="true")
            @test !ispath(joinpath(root, "generated"))
        end
    end
    # Supplied AGF inputs cannot alias a generated output, even inside generated.
    mktempdir() do root
        optics = write_zemax_fixture(root)
        model = zemax_model_fixture(); model["catalog_order"] = ["simulation.yaml"]
        YAML.write_file(joinpath(optics, "zemax.yaml"), model)
        path = joinpath(root, "generated", "simulation.yaml")
        write(path, "NM UNUSED 2\n")
        @test_throws ArgumentError Z.prepare(root; catalog_dir=dirname(path), yamlfmt="true")
        @test read(path, String) == "NM UNUSED 2\n"
        @test !isfile(joinpath(root, "generated", "test.yaml"))
    end
end

@testset "frozen catalog precedence matches the authored recipe" begin
    mktempdir() do root
        optics = write_zemax_fixture(root)
        write(joinpath(root, "test.AGF"), "NM UNUSED 2\n")
        write(joinpath(root, "other.AGF"), "NM UNUSED 2\n")
        @test Z.prepare(root; catalog_dir=root, yamlfmt="true") === nothing
        model = zemax_model_fixture(); model["catalog_order"] = ["other.AGF", "test.AGF"]
        YAML.write_file(joinpath(optics, "zemax.yaml"), model)
        paths = ("generated/simulation.yaml", "generated/test.yaml", "generated/optics/resolved.yaml",
            "generated/optics/glasses.yaml", "generated/optics/provenance.yaml")
        frozen = Dict(path => read(joinpath(root, path)) for path in paths)
        @test_throws ArgumentError Z.prepare(root; yamlfmt="true")
        @test all(read(joinpath(root, path)) == bytes for (path, bytes) in frozen)
        @test Z.prepare(root; catalog_dir=root, yamlfmt="true") === nothing
        @test Z.validate_provenance(root)["schema_version"] == 1
        # Reversing the same files also changes precedence and requires regeneration.
        reverse!(model["catalog_order"])
        YAML.write_file(joinpath(optics, "zemax.yaml"), model)
        @test_throws ArgumentError Z.prepare(root; yamlfmt="true")
        @test Z.prepare(root; catalog_dir=root, yamlfmt="true") === nothing
        @test Z.validate_provenance(root)["schema_version"] == 1
    end
end

@testset "generated output hard links cannot replace authored inputs" begin
    mktempdir() do root
        write_zemax_fixture(root)
        source = joinpath(root, "simulation", "zemax.yaml")
        output = joinpath(root, "generated", "simulation.yaml")
        Base.Filesystem.hardlink(source, output)
        original = read(source)
        @test_throws ArgumentError Z.prepare(root; yamlfmt="true")
        @test read(source) == original
        @test !isfile(joinpath(root, "generated", "test.yaml"))
    end
    mktempdir() do root
        write_zemax_fixture(root)
        output = joinpath(root, "generated", "test.yaml")
        receipt = joinpath(root, "generated", "optics", "provenance.yaml")
        write(output, "preserve aliased outputs")
        Base.Filesystem.hardlink(output, receipt)
        @test_throws ArgumentError Z.prepare(root; yamlfmt="true")
        @test read(output, String) == "preserve aliased outputs"
    end
end
