function _sync_backend!(x)
    AdaptiveOpticsSim.synchronize_backend!(AdaptiveOpticsSim.execution_style(x))
    return x
end

function _time_phase(f)
    t0 = time_ns()
    value = f()
    dt = time_ns() - t0
    return value, dt
end

function run_gpu_model_tomography_phase_profile(::Type{B}) where {B<:GPUBackendTag}
    AdaptiveOpticsSim.disable_scalar_backend!(B)
    BackendArray = AdaptiveOpticsSim.gpu_backend_array_type(B)
    BackendArray === nothing && error("GPU backend $(B) is not available")

    policy = AdaptiveOpticsSim.default_gpu_precision_policy(B)
    TB = AdaptiveOpticsSim.gpu_build_type(policy)
    build_backend = AdaptiveOpticsSim.GPUArrayBuildBackend(B)

    n_lenslet = 3
    n_lgs = 2
    n_fit_src = 2
    n_dm = 2
    grid_side = 2 * n_lenslet

    atmosphere = TomographyAtmosphereParams(
        zenith_angle_deg=TB(0.0),
        altitude_km=TB[0.0, 10.0],
        L0=TB(25.0),
        r0_zenith=TB(0.2),
        fractional_r0=TB[0.6, 0.4],
        wavelength=TB(500e-9),
        wind_direction_deg=TB[0.0, 45.0],
        wind_speed=TB[10.0, 20.0],
    )
    asterism = LGSAsterismParams(
        radius_arcsec=TB(7.6),
        wavelength=TB(589e-9),
        base_height_m=TB(90_000.0),
        n_lgs=n_lgs,
    )
    wfs = LGSWFSParams(
        diameter=TB(8.0),
        n_lenslet=n_lenslet,
        n_px=8,
        field_stop_size_arcsec=TB(2.0),
        valid_lenslet_map=trues(n_lenslet, n_lenslet),
        lenslet_rotation_rad=zeros(TB, n_lenslet^2),
        lenslet_offset=zeros(TB, 2, n_lenslet^2),
    )
    tomography = TomographyParams(
        n_fit_src=n_fit_src,
        fov_optimization_arcsec=TB(15.0),
        fit_src_height_m=TB(Inf),
    )
    dm = TomographyDMParams(
        heights_m=collect(TB, range(TB(0.0), length=n_dm, step=TB(6000.0))),
        pitch_m=fill(TB(8.0 / n_lenslet), n_dm),
        cross_coupling=TB(0.2),
        n_actuators=fill(grid_side, n_dm),
        valid_actuators=trues(grid_side, grid_side),
    )
    noise_model = RelativeSignalNoise(TB(0.1))

    gamma_single, t_gamma_single = _time_phase() do
        AdaptiveOpticsSim.sparse_gradient_matrix(AdaptiveOpticsSim.valid_lenslet_support(wfs); over_sampling=2)
    end
    gamma_base, grid_mask = gamma_single
    gamma, t_blockdiag = _time_phase() do
        blockdiag(ntuple(_ -> gamma_base, asterism.n_lgs)...)
    end

    cxx, t_cxx = _time_phase() do
        value = AdaptiveOpticsSim.auto_correlation(build_backend, atmosphere, asterism, wfs, grid_mask)
        _sync_backend!(value)
    end
    cross, t_cross = _time_phase() do
        value = AdaptiveOpticsSim.cross_correlation(build_backend, atmosphere, asterism, wfs, tomography)
        _sync_backend!(value)
    end

    weights = AdaptiveOpticsSim._equal_fit_source_weights(tomography)
    cox_full, t_fit_average = _time_phase() do
        value = AdaptiveOpticsSim._fit_source_average(cross, weights)
        _sync_backend!(value)
    end

    row_positions = findall(vec(grid_mask))
    col_positions = findall(repeat(vec(grid_mask), asterism.n_lgs))
    cox, t_extract = _time_phase() do
        value = AdaptiveOpticsSim._extract_submatrix(cox_full, row_positions, col_positions, build_backend)
        _sync_backend!(value)
    end

    gamma_native, t_gamma_native = _time_phase() do
        value = AdaptiveOpticsSim.materialize_build(build_backend, gamma, gamma)
        _sync_backend!(value)
    end
    cxx_native, t_cxx_native = _time_phase() do
        value = AdaptiveOpticsSim.materialize_build(build_backend, gamma_native, cxx)
        _sync_backend!(value)
    end
    cox_native, t_cox_native = _time_phase() do
        value = AdaptiveOpticsSim.materialize_build(build_backend, gamma_native, cox)
        _sync_backend!(value)
    end
    native_mask, t_mask_native = _time_phase() do
        value = AdaptiveOpticsSim.materialize_build(build_backend, gamma_native, grid_mask)
        _sync_backend!(value)
    end

    css_signal, t_css_signal = _time_phase() do
        value = gamma_native * cxx_native * transpose(gamma_native)
        _sync_backend!(value)
    end
    cnz, t_cnz = _time_phase() do
        value = AdaptiveOpticsSim.tomography_noise_covariance(build_backend, noise_model, diag(css_signal))
        _sync_backend!(value)
    end
    css, t_css = _time_phase() do
        value = css_signal .+ cnz
        _sync_backend!(value)
    end

    rhs, t_rhs = _time_phase() do
        value = cox_native * transpose(gamma_native)
        _sync_backend!(value)
    end
    recstat, t_recstat = _time_phase() do
        value = AdaptiveOpticsSim.stable_hermitian_right_division(build_backend, rhs, css)
        _sync_backend!(value)
    end

    d = AdaptiveOpticsSim.support_diameter(wfs) / size(AdaptiveOpticsSim.valid_lenslet_support(wfs), 1)
    wavefront_to_meter = asterism.wavelength / d / 2
    recon, t_recon = _time_phase() do
        value = d * wavefront_to_meter .* recstat
        _sync_backend!(value)
    end

    total_ns = t_gamma_single + t_blockdiag + t_cxx + t_cross + t_fit_average + t_extract +
               t_gamma_native + t_cxx_native + t_cox_native + t_mask_native + t_css_signal +
               t_cnz + t_css + t_rhs + t_recstat + t_recon

    println("GPU model tomography phase profile")
    println("  backend: ", string(something(AdaptiveOpticsSim.gpu_backend_name(B), B)))
    println("  case: medium")
    println("  gamma_single_ns: ", t_gamma_single)
    println("  gamma_blockdiag_ns: ", t_blockdiag)
    println("  auto_correlation_ns: ", t_cxx)
    println("  cross_correlation_ns: ", t_cross)
    println("  fit_source_average_ns: ", t_fit_average)
    println("  extract_submatrix_ns: ", t_extract)
    println("  materialize_gamma_ns: ", t_gamma_native)
    println("  materialize_cxx_ns: ", t_cxx_native)
    println("  materialize_cox_ns: ", t_cox_native)
    println("  materialize_mask_ns: ", t_mask_native)
    println("  css_signal_ns: ", t_css_signal)
    println("  cnz_ns: ", t_cnz)
    println("  css_sum_ns: ", t_css)
    println("  rhs_ns: ", t_rhs)
    println("  recstat_ns: ", t_recstat)
    println("  recon_scale_ns: ", t_recon)
    println("  total_timed_ns: ", total_ns)
    println("  reconstructor_type: ", typeof(recon))
    println("  native_mask_type: ", typeof(native_mask))
    return nothing
end
