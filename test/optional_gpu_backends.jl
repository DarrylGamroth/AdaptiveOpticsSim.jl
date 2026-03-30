@testset "Optional AMDGPU smoke" begin
    amdgpu_path = Base.find_package("AMDGPU")
    if amdgpu_path === nothing
        @info "Skipping AMDGPU smoke: AMDGPU.jl is not available in this environment"
        @test true
    else
        @eval import AMDGPU
        if !AMDGPU.functional()
            @info "Skipping AMDGPU smoke: ROCm runtime/device is not functional on this host"
            @test true
        else
            AdaptiveOpticsSim.disable_scalar_backend!(AMDGPUBackendTag)
            backend = AdaptiveOpticsSim.gpu_backend_array_type(AMDGPUBackendTag)
            @test backend !== nothing

            if get(ENV, "ADAPTIVEOPTICS_TEST_FULL_AMDGPU", "0") == "1"
                include(joinpath(dirname(@__DIR__), "scripts", "gpu_smoke_contract.jl"))
                run_gpu_smoke_matrix(AMDGPUBackendTag)
                @test true
            else
                T = Float32
                rng = MersenneTwister(7)
                BackendArray = backend

                tel = Telescope(resolution=16, diameter=8.0f0, sampling_time=1.0f-3,
                    central_obstruction=0.0f0, T=T, backend=BackendArray)
                src = Source(band=:I, magnitude=0.0, T=T)

                atm = MultiLayerAtmosphere(tel;
                    r0=T(0.2),
                    L0=T(25.0),
                    fractional_cn2=T[0.7, 0.3],
                    wind_speed=T[8.0, 4.0],
                    wind_direction=T[0.0, 90.0],
                    altitude=T[0.0, 5000.0],
                    T=T,
                    backend=BackendArray,
                )
                advance!(atm, tel; rng=rng)
                propagate!(atm, tel, src)
                @test atm.state.opd isa BackendArray
                @test tel.state.opd isa BackendArray

                inf_atm = InfiniteMultiLayerAtmosphere(tel;
                    r0=T(0.2),
                    L0=T(25.0),
                    fractional_cn2=T[0.7, 0.3],
                    wind_speed=T[8.0, 4.0],
                    wind_direction=T[0.0, 90.0],
                    altitude=T[0.0, 5000.0],
                    screen_resolution=33,
                    stencil_size=35,
                    T=T,
                    backend=BackendArray,
                )
                advance!(inf_atm, tel; rng=rng)
                propagate!(inf_atm, tel, src)
                @test inf_atm.state.opd isa BackendArray
                @test inf_atm.layers[1].screen.state.screen isa BackendArray

                prop = AtmosphericFieldPropagation(atm, tel, src;
                    model=GeometricAtmosphericPropagation(T=T),
                    zero_padding=2,
                    T=T)
                field = propagate_atmosphere_field!(prop, atm, tel, src)
                @test field.state.field isa BackendArray
                intensity = atmospheric_intensity!(prop, atm, tel, src)
                @test intensity isa BackendArray

                bundle = SpectralBundle(T[0.9 * wavelength(src), 1.1 * wavelength(src)], T[0.4, 0.6]; T=T)
                poly = with_spectrum(src, bundle)
                sh = ShackHartmann(tel; n_subap=4, mode=Diffractive(), T=T, backend=BackendArray)
                slopes = measure!(sh, tel, poly)
                @test slopes isa BackendArray

                curv = CurvatureWFS(tel; n_subap=4, T=T, backend=BackendArray)
                curv_slopes = measure!(curv, tel, src, atm)
                @test curv_slopes isa BackendArray
            end
        end
    end
end
