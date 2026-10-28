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
            "apodizer_surface" => 0, "focal_mask_surface" => 10,
            "lyot_surfaces" => Dict("test" => 2),
            "lyot_modes" => Dict("test" => "transmit")))
        selection = Dict("detector_surface" => 4)
        settings = Dict("resolution" => 16, "pupil_resolution" => 8,
            "wavelength_um" => 1.55, "temperature_c" => 20)
        masks = Dict("lyot" => Dict("circles" => [
            Dict("radius_m" => 0.001, "x_m" => 0, "y_m" => 0)]))
        convert_train(b) = Z.train(b, "test", selection, model, settings,
            masks, Dict(), 0.012)
        graph, diagnostics = convert_train(branch)
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
