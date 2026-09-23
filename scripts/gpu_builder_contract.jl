using AdaptiveOpticsSim
using AdaptiveOpticsSim.Optics
using AdaptiveOpticsSim.Calibration
using AdaptiveOpticsSim.Tomography
using LinearAlgebra

function run_gpu_builder_smoke(::Type{B}) where {B<:AdaptiveOpticsSim.Backends.GPUBackendTag}
    AdaptiveOpticsSim.Backends.disable_scalar_backend!(B)
    BackendArray = AdaptiveOpticsSim.Backends.gpu_backend_array_type(B)
    BackendArray === nothing && error("GPU backend $(B) is not available")

    T = Float32
    build_backend = Calibration.GPUArrayBuildBackend(B)

    A = AdaptiveOpticsSim.Backends.backend_rand(B, T, 8, 4)
    control_matrix = ControlMatrix(A; build_backend=build_backend)
    @assert control_matrix.M isa BackendArray

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
    println("gpu_builder_smoke complete")
    return nothing
end
