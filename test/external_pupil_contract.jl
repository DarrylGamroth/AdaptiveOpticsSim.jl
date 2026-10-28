using Test
using AdaptiveOpticsSim.AlgorithmGraphs

@testset "external pupil contract is shared and preserves existing identities" begin
    arguments = (; pupil_resolution=128, diameter_m=0.012,
        plane="upstream_of_apodizer", curvature="collimated")
    contract = AlgorithmGraphs.external_pupil_contract(; arguments...)
    @test contract.schema == "org.adaptiveopticssim.pupil-field/1;sha256=77300590e4a9770a5ff1666e4d2be17a41d19cce7a8f5e031a07877fef3f750c"
    @test contract.shape == (2,128,128)
    @test contract.properties isa Dict{String,String}
    @test contract.properties["optics.sampling-m"] == string(0.012/128)
    @test contract.properties["optics.spectral-model"] == "grey_pupil_amplitude"
    @test contract.properties["optics.axis-order"] == "quantity,y,x"
    @test contract.properties["optics.quantities"] == "opd_m,amplitude"
    @test contract.properties["optics.registration"] == "provisional"
    @test AlgorithmGraphs.external_pupil_contract(; arguments..., registration="declared").schema == contract.schema
    for change in ((; pupil_resolution=240), (; diameter_m=0.013), (; plane="woofer_pupil"))
        @test AlgorithmGraphs.external_pupil_contract(; merge(arguments, change)...).schema != contract.schema
    end
    legacy = AlgorithmGraphs.external_pupil_contract(; arguments..., wavelength_um=1.55)
    @test legacy.schema == "org.adaptiveopticssim.pupil-field/1;sha256=193d0b331073706ed4937ee4ec3245a37e504bf2952b74b7c12700324811cabb"
    @test legacy.properties["optics.wavelength-m"] == string(1.55e-6)
    @test !haskey(legacy.properties, "optics.spectral-model")
    for change in ((; pupil_resolution=0), (; pupil_resolution=true), (; pupil_resolution=2.5),
            (; diameter_m=0), (; diameter_m=NaN), (; diameter_m=Inf), (; diameter_m=true),
            (; plane=""), (; plane="pupil\nother"), (; curvature="converging"),
            (; registration=""), (; spectral_model="chromatic"),
            (; wavelength_um=0), (; wavelength_um=Inf), (; wavelength_um=true))
        @test_throws ArgumentError AlgorithmGraphs.external_pupil_contract(; merge(arguments, change)...)
    end
    @test all(id -> id.name ∉ ("CUDA", "AMDGPU", "Proper"), keys(Base.loaded_modules))
end
