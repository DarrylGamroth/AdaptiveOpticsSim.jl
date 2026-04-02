backend_package_name(::Type{CUDABackendTag}) = "CUDA"
backend_package_name(::Type{AMDGPUBackendTag}) = "AMDGPU"

backend_label(::Type{CUDABackendTag}) = "CUDA"
backend_label(::Type{AMDGPUBackendTag}) = "AMDGPU"

backend_full_smoke_env(::Type{CUDABackendTag}) = "ADAPTIVEOPTICS_TEST_FULL_CUDA"
backend_full_smoke_env(::Type{AMDGPUBackendTag}) = "ADAPTIVEOPTICS_TEST_FULL_AMDGPU"

function import_backend_package!(::Type{CUDABackendTag})
    @eval import CUDA
    return nothing
end

function import_backend_package!(::Type{AMDGPUBackendTag})
    @eval import AMDGPU
    return nothing
end

backend_functional(::Type{CUDABackendTag}) = CUDA.functional()
backend_functional(::Type{AMDGPUBackendTag}) = AMDGPU.functional()

run_optional_backend_plan_checks(::Type{<:GPUBackendTag}, tel, backend) = nothing

function run_optional_backend_plan_checks(::Type{AMDGPUBackendTag}, tel, backend)
    T = Float32
    pyr = PyramidWFS(tel; n_subap=4, modulation=T(1.0), mode=Diffractive(), T=T, backend=backend)
    bio = BioEdgeWFS(tel; n_subap=4, modulation=T(1.0), mode=Diffractive(), T=T, backend=backend)
    @test AdaptiveOpticsSim.grouped_accumulation_plan(AdaptiveOpticsSim.execution_style(pyr.state.intensity), pyr) isa AdaptiveOpticsSim.GroupedStaged2DPlan
    @test AdaptiveOpticsSim.grouped_accumulation_plan(AdaptiveOpticsSim.execution_style(bio.state.intensity), bio) isa AdaptiveOpticsSim.GroupedStaged2DPlan
    return nothing
end

function run_optional_backend_plan_checks(::Type{CUDABackendTag}, tel, backend)
    T = Float32
    pyr = PyramidWFS(tel; n_subap=4, modulation=T(1.0), mode=Diffractive(), T=T, backend=backend)
    bio = BioEdgeWFS(tel; n_subap=4, modulation=T(1.0), mode=Diffractive(), T=T, backend=backend)
    @test AdaptiveOpticsSim.grouped_accumulation_plan(AdaptiveOpticsSim.execution_style(pyr.state.intensity), pyr) isa AdaptiveOpticsSim.GroupedStackReducePlan
    @test AdaptiveOpticsSim.grouped_accumulation_plan(AdaptiveOpticsSim.execution_style(bio.state.intensity), bio) isa AdaptiveOpticsSim.GroupedStackReducePlan
    return nothing
end

function run_optional_backend_smoke(::Type{B}) where {B<:GPUBackendTag}
    pkg = backend_package_name(B)
    pkg_path = Base.find_package(pkg)
    if pkg_path === nothing
        @info "Skipping $(backend_label(B)) smoke: $(pkg).jl is not available in this environment"
        @test true
        return nothing
    end

    import_backend_package!(B)
    if !backend_functional(B)
        @info "Skipping $(backend_label(B)) smoke: backend runtime/device is not functional on this host"
        @test true
        return nothing
    end

    AdaptiveOpticsSim.disable_scalar_backend!(B)
    backend = AdaptiveOpticsSim.gpu_backend_array_type(B)
    @test backend !== nothing

    if get(ENV, backend_full_smoke_env(B), "0") == "1"
        include(joinpath(dirname(@__DIR__), "scripts", "gpu_smoke_contract.jl"))
        run_gpu_smoke_matrix(B)
        @test true
        return nothing
    end

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

    run_optional_backend_plan_checks(B, tel, BackendArray)

    curv = CurvatureWFS(tel; n_subap=4, T=T, backend=BackendArray)
    curv_slopes = measure!(curv, tel, src, atm)
    @test curv_slopes isa BackendArray
    return nothing
end
