using Test
using LinearAlgebra
using TOML

@testset "AppleAccelerate remains application-owned" begin
    @test Sys.isapple()
    @test Sys.ARCH === :aarch64
    @test !any(mod -> nameof(mod) === :AppleAccelerate,
        values(Base.loaded_modules))
end

using AdaptiveOpticsSim
import FFTW

@testset "Backend-neutral package load" begin
    @test !any(mod -> nameof(mod) === :AppleAccelerate,
        values(Base.loaded_modules))

    project = TOML.parsefile(joinpath(pkgdir(AdaptiveOpticsSim),
        "Project.toml"))
    @test !haskey(project["deps"], "AppleAccelerate")
    @test haskey(project["weakdeps"], "AppleAccelerate")
    @test project["extensions"]["AdaptiveOpticsSimAppleAccelerateExt"] ==
        "AppleAccelerate"
    @test Base.get_extension(AdaptiveOpticsSim,
        :AdaptiveOpticsSimAppleAccelerateExt) === nothing

    buffer = zeros(ComplexF64, 8, 8)
    plan = AdaptiveOpticsSim.plan_fft_backend!(buffer)
    @test parentmodule(typeof(plan)) === FFTW
end
