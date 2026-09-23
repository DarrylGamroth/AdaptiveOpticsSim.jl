using AdaptiveOpticsSim
using AdaptiveOpticsSim.Optics
using AdaptiveOpticsSim.Calibration
using AdaptiveOpticsSim.Tomography
using LinearAlgebra
using Test: @inferred

function run_gpu_covariance_geometry_contract(::Type{B}) where
    {B<:AdaptiveOpticsSim.Backends.GPUBackendTag}
    T = Float32
    build_backend = Calibration.GPUArrayBuildBackend(B)
    BackendArray = AdaptiveOpticsSim.Backends.gpu_backend_array_type(B)
    atmosphere = TomographyAtmosphereParams(
        zenith_angle_deg=T(7),
        layer_altitudes_m=T[0, 10_000],
        L0=T(25),
        r0_zenith=T(0.2),
        fractional_cn2=T[0.6, 0.4],
        reference_wavelength_m=T(500e-9),
        wind_direction_deg=T[0, 90],
        wind_speed=T[5, 10],
    )
    asterism = LGSAsterismParams(
        radius_arcsec=T(7.6),
        wavelength_m=T(589e-9),
        base_height_m=T(90_000),
        n_lgs=2,
    )
    wfs = LGSWFSParams(
        pupil_diameter_m=T(8),
        n_lenslets=2,
        n_px=8,
        field_stop_size_arcsec=T(2),
        valid_lenslet_map=trues(2, 2),
        lenslet_grid_rotations_rad=T[0.06, -0.04],
        lenslet_grid_offsets_fraction=T[0.03 -0.02; -0.01 0.04],
    )
    tomography = TomographyParams(
        n_fit_src=2,
        fov_optimization_arcsec=T(15),
        fit_src_height_m=T(120_000),
    )
    grid_mask = Bool[true false true; false true false; true false true]
    cxx_gpu = @inferred AdaptiveOpticsSim.Tomography.auto_correlation(
        build_backend, atmosphere, asterism, wfs, grid_mask)
    cox_gpu = @inferred AdaptiveOpticsSim.Tomography.cross_correlation(
        build_backend, atmosphere, asterism, wfs, tomography;
        grid_mask=grid_mask)
    cxx_cpu = AdaptiveOpticsSim.Tomography.auto_correlation(
        atmosphere, asterism, wfs, grid_mask)
    cox_cpu = AdaptiveOpticsSim.Tomography.cross_correlation(
        atmosphere, asterism, wfs, tomography; grid_mask=grid_mask)
    @assert cxx_gpu isa BackendArray
    @assert cox_gpu isa BackendArray
    @assert size(cxx_gpu) == (10, 10)
    @assert size(cox_gpu) == (4, 5, 10)
    @assert norm(Array(cxx_gpu) - cxx_cpu) / norm(cxx_cpu) <= 2f-4
    @assert norm(Array(cox_gpu) - cox_cpu) / norm(cox_cpu) <= 2f-4

    empty_mask = falses(size(grid_mask))
    empty_cxx = @inferred AdaptiveOpticsSim.Tomography.auto_correlation(
        build_backend, atmosphere, asterism, wfs, empty_mask)
    empty_cox = @inferred AdaptiveOpticsSim.Tomography.cross_correlation(
        build_backend, atmosphere, asterism, wfs, tomography;
        grid_mask=empty_mask)
    @assert empty_cxx isa BackendArray && size(empty_cxx) == (0, 0)
    @assert empty_cox isa BackendArray && size(empty_cox) == (4, 0, 0)

    projection_host = T[1 0; 0 1]
    phase_host = T[2 0.25; 0.25 3]
    fit_host = reshape(T[1, 2], 1, 2)
    noise_host = T[0.1 0; 0 0.2]
    projection = Calibration.materialize_build(build_backend, projection_host)
    phase = Calibration.materialize_build(build_backend, phase_host)
    fit = Calibration.materialize_build(build_backend, fit_host)
    noise = Calibration.materialize_build(build_backend, noise_host)
    recstat = @inferred AdaptiveOpticsSim.Tomography._tomographic_covariance_reconstructor(
        AdaptiveOpticsSim.Backends.execution_style(phase), build_backend,
        projection, phase, fit, noise)
    expected = fit_host * transpose(projection_host) /
        (projection_host * phase_host * transpose(projection_host) + noise_host)
    @assert recstat isa BackendArray
    @assert Array(recstat) ≈ expected rtol=3f-5 atol=3f-5
    return nothing
end

function run_gpu_builder_smoke(::Type{B}) where {B<:AdaptiveOpticsSim.Backends.GPUBackendTag}
    AdaptiveOpticsSim.Backends.disable_scalar_backend!(B)
    BackendArray = AdaptiveOpticsSim.Backends.gpu_backend_array_type(B)
    BackendArray === nothing && error("GPU backend $(B) is not available")

    T = Float32
    build_backend = Calibration.GPUArrayBuildBackend(B)
    atm = TomographyAtmosphereParams(
        zenith_angle_deg=T(0.0),
        layer_altitudes_m=T[0.0],
        L0=T(25.0),
        r0_zenith=T(0.2),
        fractional_cn2=T[1.0],
        reference_wavelength_m=T(500e-9),
        wind_direction_deg=T[0.0],
        wind_speed=T[10.0],
    )
    lgs = LGSAsterismParams(
        radius_arcsec=T(7.6),
        wavelength_m=T(589e-9),
        base_height_m=T(90_000.0),
        n_lgs=1,
    )
    wfs = LGSWFSParams(
        pupil_diameter_m=T(8.0),
        n_lenslets=1,
        n_px=8,
        field_stop_size_arcsec=T(2.0),
        valid_lenslet_map=trues(1, 1),
        lenslet_grid_rotations_rad=zeros(T, 1),
        lenslet_grid_offsets_fraction=zeros(T, 2, 1),
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
    imat_t = Calibration.materialize_build(build_backend,
        reshape(T[1.0, 0.5], 2, 1))
    imat_t_cpu = reshape(T[1.0, 0.5], 2, 1)
    noise = AdaptiveOpticsSim.Tomography.RelativeSignalNoise(T(0.1))

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
        build_backend=Calibration.CPUBuildBackend(),
    )
    relative_interaction_matrix_error =
        norm(Array(tr.reconstructor) - tr_cpu.reconstructor) /
        norm(tr_cpu.reconstructor)
    @assert relative_interaction_matrix_error <= 2f-2 "GPU/CPU interaction-matrix reconstructor error: $relative_interaction_matrix_error"

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
        build_backend=Calibration.CPUBuildBackend(),
    )
    cxx_direct = @inferred AdaptiveOpticsSim.Tomography.auto_correlation(
        build_backend, atm, lgs, wfs, grid_mask)
    cox_direct = @inferred AdaptiveOpticsSim.Tomography.cross_correlation(
        build_backend, atm, lgs, wfs, tomo; grid_mask=grid_mask)
    @assert cxx_direct isa BackendArray
    @assert cox_direct isa BackendArray
    cxx_direct_cpu = AdaptiveOpticsSim.Tomography.auto_correlation(
        atm, lgs, wfs, grid_mask)
    cox_direct_cpu = AdaptiveOpticsSim.Tomography.cross_correlation(
        atm, lgs, wfs, tomo; grid_mask=grid_mask)
    @assert Array(cxx_direct) ≈ cxx_direct_cpu rtol=2f-2
    @assert Array(cox_direct) ≈ cox_direct_cpu rtol=2f-2
    cxx_gpu = Array(mr.operators.cxx)
    cox_gpu = Array(mr.operators.cox)
    relative_grid_covariance_error = norm(cox_gpu - cxx_gpu) / norm(cxx_gpu)
    relative_cpu_covariance_error = norm(cxx_gpu - mr_cpu.operators.cxx) /
        norm(mr_cpu.operators.cxx)
    relative_reconstructor_error = norm(Array(mr.reconstructor) - mr_cpu.reconstructor) /
        norm(mr_cpu.reconstructor)
    @assert relative_grid_covariance_error <= 1f-3 "GPU Cox/Cxx grid error: $relative_grid_covariance_error"
    @assert relative_cpu_covariance_error <= 2f-2 "GPU/CPU Cxx error: $relative_cpu_covariance_error"
    @assert relative_reconstructor_error <= 2f-2 "GPU/CPU reconstructor error: $relative_reconstructor_error"
    run_gpu_covariance_geometry_contract(B)
    println("gpu_builder_smoke complete")
    return nothing
end
