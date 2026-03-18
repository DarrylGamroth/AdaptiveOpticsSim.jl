using AdaptiveOpticsSim

try
    using CUDA
catch err
    error("gpu_profile_auto_correlation_cuda.jl requires CUDA.jl: $(sprint(showerror, err))")
end

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

function _profile_case(::Type{B}) where {B<:GPUBackendTag}
    AdaptiveOpticsSim.disable_scalar_backend!(B)
    BackendArray = AdaptiveOpticsSim.gpu_backend_array_type(B)
    BackendArray === nothing && error("GPU backend $(B) is not available")

    policy = AdaptiveOpticsSim.default_gpu_precision_policy(B)
    T = AdaptiveOpticsSim.gpu_build_type(policy)
    backend = GPUArrayBuildBackend(B)

    n_lenslet = 3
    n_lgs = 2

    atmosphere = TomographyAtmosphereParams(
        zenith_angle_deg=T(0.0),
        altitude_km=T[0.0, 10.0],
        L0=T(25.0),
        r0_zenith=T(0.2),
        fractional_r0=T[0.6, 0.4],
        wavelength=T(500e-9),
        wind_direction_deg=T[0.0, 45.0],
        wind_speed=T[10.0, 20.0],
    )
    asterism = LGSAsterismParams(
        radius_arcsec=T(7.6),
        wavelength=T(589e-9),
        base_height_m=T(90_000.0),
        n_lgs=n_lgs,
    )
    wfs = LGSWFSParams(
        diameter=T(8.0),
        n_lenslet=n_lenslet,
        n_px=8,
        field_stop_size_arcsec=T(2.0),
        valid_lenslet_map=trues(n_lenslet, n_lenslet),
        lenslet_rotation_rad=zeros(T, n_lenslet^2),
        lenslet_offset=zeros(T, 2, n_lenslet^2),
    )

    _, grid_mask = AdaptiveOpticsSim.sparse_gradient_matrix(AdaptiveOpticsSim.valid_lenslet_support(wfs); over_sampling=2)
    sampling = size(grid_mask, 1)
    mask_vec = vec(grid_mask)
    valid_positions = findall(mask_vec)
    n_valid = count(mask_vec)
    altitude = AdaptiveOpticsSim.layer_altitude_m(atmosphere)
    r0 = AdaptiveOpticsSim._fried_parameter(atmosphere)
    support_d = AdaptiveOpticsSim.support_diameter(wfs)
    lgs_dir = AdaptiveOpticsSim.lgs_directions(asterism)
    directions = AdaptiveOpticsSim.direction_vectors(view(lgs_dir, :, 1), view(lgs_dir, :, 2))
    source_height = AdaptiveOpticsSim.lgs_height_m(asterism, atmosphere)

    result = AdaptiveOpticsSim._backend_array(B, T, n_lgs * n_valid, n_lgs * n_valid)
    fill!(result, zero(T))
    valid_positions_native = AdaptiveOpticsSim._backend_array(B, Int, length(valid_positions))
    copyto!(valid_positions_native, valid_positions)
    style = AdaptiveOpticsSim.execution_style(result)

    guide_xy, t_guides = _time_phase() do
        gx, gy = AdaptiveOpticsSim._guide_star_grids(
            backend,
            sampling,
            support_d,
            wfs.lenslet_rotation_rad,
            view(wfs.lenslet_offset, 1, :),
            view(wfs.lenslet_offset, 2, :),
        )
        _sync_backend!(gx)
        _sync_backend!(gy)
        (gx, gy)
    end
    guide_x, guide_y = guide_xy

    shifted, t_shift = _time_phase() do
        value = AdaptiveOpticsSim._scaled_shifted_coord_stack(
            backend,
            guide_x,
            guide_y,
            directions,
            altitude,
            source_height,
        )
        _sync_backend!(value)
    end

    cst, var_term, inv_L0 = AdaptiveOpticsSim._covariance_constants(r0, atmosphere.L0)
    block = AdaptiveOpticsSim._backend_array(B, T, n_valid, n_valid)
    cov = AdaptiveOpticsSim._backend_array(B, T, sampling * sampling, sampling * sampling)

    t_cov = 0
    t_accumulate = 0
    t_scatter = 0

    for jgs in 1:n_lgs
        for igs in 1:jgs
            fill!(block, zero(T))
            for layer in eachindex(altitude)
                iz = @view shifted[:, :, igs, layer]
                jz = @view shifted[:, :, jgs, layer]

                _, dt_cov = _time_phase() do
                    AdaptiveOpticsSim._covariance_matrix!(
                        backend,
                        cov,
                        vec(iz),
                        vec(jz),
                        cst,
                        var_term,
                        inv_L0,
                        atmosphere.fractional_r0[layer],
                    )
                    _sync_backend!(cov)
                end
                t_cov += dt_cov

                _, dt_accumulate = _time_phase() do
                    AdaptiveOpticsSim.launch_kernel_async!(style, AdaptiveOpticsSim.accumulate_selected_block_kernel!, block, cov,
                        valid_positions_native, n_valid; ndrange=size(block))
                    _sync_backend!(block)
                end
                t_accumulate += dt_accumulate
            end

            _, dt_scatter = _time_phase() do
                rows = (igs - 1) * n_valid + 1:igs * n_valid
                cols = (jgs - 1) * n_valid + 1:jgs * n_valid
                result[rows, cols] .= block
                if igs != jgs
                    result[cols, rows] .= transpose(block)
                end
                _sync_backend!(result)
            end
            t_scatter += dt_scatter
        end
    end

    println("GPU auto_correlation phase profile")
    println("  backend: ", string(something(AdaptiveOpticsSim.gpu_backend_name(B), B)))
    println("  case: medium")
    println("  guide_grids_ns: ", t_guides)
    println("  scaled_shifted_coords_ns: ", t_shift)
    println("  covariance_matrix_ns: ", t_cov)
    println("  selected_block_accumulate_ns: ", t_accumulate)
    println("  block_scatter_ns: ", t_scatter)
    println("  total_timed_ns: ", t_guides + t_shift + t_cov + t_accumulate + t_scatter)
    println("  result_type: ", typeof(result))
    return nothing
end

_profile_case(CUDABackendTag)
