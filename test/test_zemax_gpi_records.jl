using Test
using SHA
using AdaptiveOpticsSim

const ZGPI = AdaptiveOpticsSim.Optics.Zemax

function gpi_zmx_fixture(path, text)
    open(path, "w") do io
        write(io, UInt8[0xff, 0xfe])
        for unit in transcode(UInt16, text)
            write(io, UInt8(unit & 0xff), UInt8(unit >> 8))
        end
    end
    return path
end

const GPI_MCE_RECORDS = """
UNIT MM
NAME bounded GPI operand fixture
NOTE 0 System stop and saved surface diameters are different quantities.
ENVD 10 0.75 0
SURF 0
  TYPE STANDARD
  DISZ 100
SURF 1
  TYPE STANDARD
  DISZ 999
  DIAM 3 0
SURF 2
  TYPE STANDARD
  DISZ 0
SURF 3
  TYPE STANDARD
  GLAS F_SILICA
SURF 4
  TYPE STANDARD
MNUM 2 1
MOFF 0 1 "NO ADC" 0
MOFF 0 2 "WITH ADC" 0
THIC 0 1 100 0 0
THIC 0 2 200 0 0
THIC 1 1 999 2 0 0 1 2 -1 0 0
THIC 1 2 999 3 0 0 1 3 1 10 0
SATP 0 1 3 0 0
SATP 0 2 3 0 0
STPS 0 1 1 0 0
STPS 0 2 1 0 0
APER 0 1 2 0 0
APER 0 2 2.5 0 0
FLTP 0 1 1 0 0
FLTP 0 2 1 0 0
GLSS 3 1 "" 0 0
GLSS 3 2 "N-SK2" 0 0
IGNR 4 1 0 0 0
IGNR 4 2 1 0 0
XFIE 2 1 0.62 0 0
XFIE 2 2 0.62 0 0
YFIE 2 1 -0.62 0 0
YFIE 2 2 -0.62 0 0
"""

@testset "GPI system operands, materials and bounded pickups" begin
    mktempdir() do dir
        path = gpi_zmx_fixture(joinpath(dir, "gpi.ZMX"), GPI_MCE_RECORDS)
        d = ZGPI.inspect(path; configurations=[1, 2])
        a, b = d["configurations"]
        surface(branch, id) = only(s for s in branch["active_surfaces"] if s["id"] == id)
        @test surface(a, 1)["thickness_m"] == -0.1
        # The reference is configuration 1, row 3 (MOFF counts as row 1),
        # not configuration 2's THIC 0 or the stale saved value 999.
        @test surface(b, 1)["thickness_m"] == -0.09
        @test surface(a, 3)["glass_after"] == ""
        @test surface(b, 3)["glass_after"] == "N-SK2"
        @test b["ignored_surface_ids"] == [4]
        @test a["ignored_surface_ids"] == Int[]
        @test b["system_overrides_native"] == Dict("SATP" => 3, "STPS" => 1, "APER" => 2.5, "FLTP" => 1)
        @test surface(b, 1)["semi_diameter_m"] == 0.003
        @test b["field_overrides_native"] == [
            ZGPI.FieldDocument("operand" => "XFIE", "field_id" => 2, "value_native" => 0.62),
            ZGPI.FieldDocument("operand" => "YFIE", "field_id" => 2, "value_native" => -0.62)]
        @test first(a["raw_mce_records"]) == "MOFF 0 1 \"NO ADC\" 0"
        @test length(d["notes"]) == 1
        @test "ENVD 10 0.75 0" in d["raw_system_records"]
        @test "THIC 1 2 999 3 0 0 1 3 1 10 0" in b["raw_mce_records"]
        @test d["source"]["sha256"] == bytes2hex(sha256(read(path)))

        reject(text; configs=[2]) = begin
            gpi_zmx_fixture(path, text)
            @test_throws ErrorException ZGPI.inspect(path; configurations=configs)
        end
        reject(replace(GPI_MCE_RECORDS, "THIC 1 2 999 3" => "THIC 1 2 999 1"))
        reject(replace(GPI_MCE_RECORDS, "1 3 1 10 0" => "2 3 1 10 0")) # cyclic pickup
        reject(replace(GPI_MCE_RECORDS, "1 3 1 10 0" => "1 99 1 10 0"))
        reject(replace(GPI_MCE_RECORDS, "1 3 1 10 0" => "1 1 1 10 0")) # MOFF source
        reject(replace(GPI_MCE_RECORDS, "MOFF 0 2 \"WITH ADC\" 0\n" => ""))
        reject(replace(GPI_MCE_RECORDS, "THIC 1 2 999 3" => "CRVT 1 2 999 2"))
        reject(replace(GPI_MCE_RECORDS, "  DISZ 999" => "  CURV 0.1 0\n  DISZ 999"))
        reject(replace(GPI_MCE_RECORDS, "  DISZ 999" => "  CURV 0 1\n  DISZ 999"))
        reject(replace(GPI_MCE_RECORDS, "ENVD 10 0.75 0\n" => ""))
        reject(replace(GPI_MCE_RECORDS, "ENVD 10 0.75 0" => "ENVD 10 0.75 1"))
        reject(GPI_MCE_RECORDS * "TEMP 0 1 20 0 0\n") # invalidates the shared environment proof
        reject(replace(GPI_MCE_RECORDS, "  DISZ 999" => "  PZUP 3 0 0\n  DISZ 999"))
        reject(replace(GPI_MCE_RECORDS, "GLSS 3 2 \"N-SK2\" 0" => "GLSS 3 2 \"N-SK2\" 2"))
        reject(replace(GPI_MCE_RECORDS, "GLSS 3 2 \"N-SK2\"" => "GLSS 3 2 N-SK2"))
        reject(replace(GPI_MCE_RECORDS, "SATP 0 2 3" => "SATP 0 2 0"))
        reject(replace(GPI_MCE_RECORDS, "FLTP 0 2 1" => "FLTP 0 2 0"))
        reject(replace(GPI_MCE_RECORDS, "STPS 0 2 1" => "STPS 0 2 4"))
        reject(replace(GPI_MCE_RECORDS, "APER 0 2 2.5" => "APER 0 2 -1"))
        reject(replace(GPI_MCE_RECORDS, r"SATP[^\n]*\n" => ""))
        reject(replace(GPI_MCE_RECORDS, "IGNR 4 2 1" => "IGNR 4 2 2"))
        reject(replace(GPI_MCE_RECORDS, "IGNR 4" => "IGNR 1"))
        reject(replace(GPI_MCE_RECORDS, "IGNR 4" => "IGNR 2"))
        reject(GPI_MCE_RECORDS * "IGNM 1 1 0 0 2\nIGNM 1 2 1 0 2\n")
        reject(GPI_MCE_RECORDS * "WAVE 1 2 1 0 0\n")
    end
end

# Supply the unchanged vendor source for integration validation without copying
# its contents into this package or making routine tests depend on local paths.
if haskey(ENV, "AOS_GPI_ZMX")
    @testset "GPI 2008 source configurations 1 and 3" begin
        path = ENV["AOS_GPI_ZMX"]
        bytes = read(path)
        d = ZGPI.inspect(path; configurations=[1, 3])
        @test d["source"]["surface_count"] == 109 # SURF 0 through SURF 108
        @test d["configuration_count"] == 5
        @test d["saved_configuration_id"] == 1
        @test d["source"]["sha256"] == bytes2hex(sha256(bytes))
        a, b = d["configurations"]
        @test a["name"] == "NO ADC"
        @test b["name"] == "WITH ADC"
        for branch in (a, b)
            @test only(s["thickness_m"] for s in branch["active_surfaces"] if s["id"] == 5) ≈ 16.542574776
            @test branch["system_overrides_native"] == Dict("SATP" => 3, "STPS" => 4, "APER" => 511.75, "FLTP" => 1)
            @test branch["ignored_surface_ids"] == Int[]
        end
        for id in (21, 24)
            @test only(s["glass_after"] for s in a["active_surfaces"] if s["id"] == id) == ""
        end
        @test only(s["glass_after"] for s in b["active_surfaces"] if s["id"] == 21) == "S-NPH2"
        @test only(s["glass_after"] for s in b["active_surfaces"] if s["id"] == 24) == "S-BAL42"
        @test read(path) == bytes
    end
end
