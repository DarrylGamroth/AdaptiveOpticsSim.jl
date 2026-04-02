using AdaptiveOpticsSim
using LinearAlgebra

function run_gpu_builder_smoke(::Type{B}) where {B<:GPUBackendTag}
    disable_scalar_backend!(B)
    BackendArray = gpu_backend_array_type(B)
    BackendArray === nothing && error("GPU backend $(B) is not available")

    T = Float32
    build_backend = AdaptiveOpticsSim.GPUArrayBuildBackend(B)

    A = backend_rand(B, T, 8, 4)
    imat = InteractionMatrix(A, T(0.1))
    cpu_A = Array(A)
    vault = CalibrationVault(A; build_backend=build_backend)
    @assert vault.M isa BackendArray

    recon = ModalReconstructor(imat; build_backend=build_backend)
    @assert recon.reconstructor isa BackendArray
    recon_cpu = ModalReconstructor(InteractionMatrix(cpu_A, T(0.1)); build_backend=AdaptiveOpticsSim.CPUBuildBackend())
    slopes_modal = reshape(T.(1:8), 8)
    @assert isapprox(
        Array(reconstruct(recon, AdaptiveOpticsSim.materialize_build(build_backend, slopes_modal))),
        reconstruct(recon_cpu, slopes_modal);
        rtol=1f-5,
        atol=1f-6,
    )

    atm = TomographyAtmosphereParams(
        zenith_angle_deg=T(0.0),
        altitude_km=T[0.0],
        L0=T(25.0),
        r0_zenith=T(0.2),
        fractional_r0=T[1.0],
        wavelength=T(500e-9),
        wind_direction_deg=T[0.0],
        wind_speed=T[10.0],
    )
    lgs = LGSAsterismParams(
        radius_arcsec=T(7.6),
        wavelength=T(589e-9),
        base_height_m=T(90_000.0),
        n_lgs=1,
    )
    wfs = LGSWFSParams(
        diameter=T(8.0),
        n_lenslet=1,
        n_px=8,
        field_stop_size_arcsec=T(2.0),
        valid_lenslet_map=trues(1, 1),
        lenslet_rotation_rad=zeros(T, 1),
        lenslet_offset=zeros(T, 2, 1),
    )
    tomo = TomographyParams(
        n_fit_src=1,
        fov_optimization_arcsec=T(0.0),
        fit_src_height_m=T(Inf),
    )
    dm = TomographyDMParams(
        heights_m=T[0.0],
        pitch_m=T[0.5],
        cross_coupling=T(0.2),
        n_actuators=[1],
        valid_actuators=trues(1, 1),
    )
    grid_mask = trues(1, 1)
    imat_t = AdaptiveOpticsSim.materialize_build(build_backend, reshape(T[1.0, 0.5], 2, 1))
    imat_t_cpu = reshape(T[1.0, 0.5], 2, 1)
    noise = RelativeSignalNoise(T(0.1))

    tr = build_reconstructor(
        InteractionMatrixTomography(),
        imat_t,
        grid_mask,
        atm,
        lgs,
        wfs,
        tomo,
        dm;
        noise_model=noise,
        build_backend=build_backend,
    )
    @assert tr.reconstructor isa BackendArray
    @assert tr.operators.cxx isa BackendArray
    @assert tr.operators.cox isa BackendArray
    @assert tr.operators.cnz isa BackendArray
    tr_cpu = build_reconstructor(
        InteractionMatrixTomography(),
        imat_t_cpu,
        grid_mask,
        atm,
        lgs,
        wfs,
        tomo,
        dm;
        noise_model=noise,
        build_backend=AdaptiveOpticsSim.CPUBuildBackend(),
    )
    slopes_tomo = T[0.25, -0.5]
    slopes_tomo_tr_gpu = AdaptiveOpticsSim.materialize_build(build_backend, convert.(eltype(tr.reconstructor), slopes_tomo))
    @assert isapprox(
        Array(reconstruct_wavefront(tr, slopes_tomo_tr_gpu)),
        reconstruct_wavefront(tr_cpu, slopes_tomo);
        rtol=1f-5,
        atol=1f-6,
    )
    @assert isapprox(
        Array(dm_commands(tr, slopes_tomo_tr_gpu)),
        dm_commands(tr_cpu, slopes_tomo);
        rtol=1f-5,
        atol=1f-6,
    )

    mr = build_reconstructor(
        ModelBasedTomography(),
        atm,
        lgs,
        wfs,
        tomo,
        dm;
        noise_model=noise,
        build_backend=build_backend,
    )
    @assert mr.reconstructor isa BackendArray
    @assert mr.operators.cxx isa BackendArray
    @assert mr.operators.cox isa BackendArray
    @assert mr.operators.cnz isa BackendArray
    mr_cpu = build_reconstructor(
        ModelBasedTomography(),
        atm,
        lgs,
        wfs,
        tomo,
        dm;
        noise_model=noise,
        build_backend=AdaptiveOpticsSim.CPUBuildBackend(),
    )
    slopes_tomo_mr = convert.(eltype(mr_cpu.reconstructor), slopes_tomo)
    slopes_tomo_mr_gpu = AdaptiveOpticsSim.materialize_build(build_backend, slopes_tomo_mr)
    @assert isapprox(
        Array(reconstruct_wavefront(mr, slopes_tomo_mr_gpu)),
        reconstruct_wavefront(mr_cpu, slopes_tomo_mr);
        rtol=1f-5,
        atol=1f-6,
    )

    println("gpu_builder_smoke complete")
    return nothing
end
