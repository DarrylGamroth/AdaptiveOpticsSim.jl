using AdaptiveOpticsSim
using AdaptiveOpticsSim.Optics
using AdaptiveOpticsSim.Tomography

try
    using CUDA
catch err
    error("gpu_profile_auto_correlation_cuda.jl requires CUDA.jl: $(sprint(showerror, err))")
end

function _sync_backend!(x)
    AdaptiveOpticsSim.Backends.synchronize_backend!(AdaptiveOpticsSim.Backends.execution_style(x))
    return x
end

function _time_phase(f)
    t0 = time_ns()
    value = f()
    dt = time_ns() - t0
    return value, dt
end

function _profile_case(::Type{B}) where {B<:AdaptiveOpticsSim.Backends.GPUBackendTag}
    AdaptiveOpticsSim.Backends.disable_scalar_backend!(B)
    BackendArray = AdaptiveOpticsSim.Backends.gpu_backend_array_type(B)
    BackendArray === nothing && error("GPU backend $(B) is not available")

    policy = AdaptiveOpticsSim.Backends.default_gpu_precision_policy(B)
    T = AdaptiveOpticsSim.Backends.gpu_build_type(policy)
    backend = AdaptiveOpticsSim.Calibration.GPUArrayBuildBackend(B)

    n_lenslets = 3
    n_lgs = 2

    atmosphere = TomographyAtmosphereParams(
        zenith_angle_deg=T(0.0),
        layer_altitudes_m=T[0.0, 10_000.0],
        L0=T(25.0),
        r0_zenith=T(0.2),
        fractional_cn2=T[0.6, 0.4],
        reference_wavelength_m=T(500e-9),
        wind_direction_deg=T[0.0, 45.0],
        wind_speed=T[10.0, 20.0],
    )
    asterism = LGSAsterismParams(
        radius_arcsec=T(7.6),
        wavelength_m=T(589e-9),
        base_height_m=T(90_000.0),
        n_lgs=n_lgs,
    )
    wfs = LGSWFSParams(
        pupil_diameter_m=T(8.0),
        n_lenslets=n_lenslets,
        n_px=8,
        field_stop_size_arcsec=T(2.0),
        valid_lenslet_map=trues(n_lenslets, n_lenslets),
        lenslet_grid_rotations_rad=zeros(T, n_lgs),
        lenslet_grid_offsets_fraction=zeros(T, 2, n_lgs),
    )

    _, grid_mask = AdaptiveOpticsSim.Tomography.sparse_gradient_matrix(
        AdaptiveOpticsSim.Tomography.valid_lenslet_support(wfs);
        over_sampling=2)
    sampling = size(grid_mask, 1)
    mask_vec = vec(grid_mask)
    valid_positions = findall(mask_vec)
    n_valid = count(mask_vec)
    slant_ranges_m = AdaptiveOpticsSim.Tomography.layer_slant_ranges_m(atmosphere)
    r0 = AdaptiveOpticsSim.Tomography._fried_parameter(atmosphere)
    support_diameter_m =
        AdaptiveOpticsSim.Tomography.lenslet_grid_support_diameter_m(wfs)
    lgs_dir = AdaptiveOpticsSim.Tomography.lgs_directions(asterism)
    directions = AdaptiveOpticsSim.Tomography.direction_vectors(
        view(lgs_dir, :, 1), view(lgs_dir, :, 2))
    source_height_m =
        AdaptiveOpticsSim.Tomography.lgs_height_m(asterism, atmosphere)
    rotations_rad, offset_fractions_x, offset_fractions_y =
        AdaptiveOpticsSim.Tomography._active_guide_grid_params(
        wfs.lenslet_grid_rotations_rad,
        view(wfs.lenslet_grid_offsets_fraction, 1, :),
        view(wfs.lenslet_grid_offsets_fraction, 2, :),
        n_lgs,
    )

    result = AdaptiveOpticsSim.Calibration._backend_array(
        B, T, n_lgs * n_valid, n_lgs * n_valid)
    fill!(result, zero(T))
    valid_positions_native = AdaptiveOpticsSim.Calibration._backend_array(
        B, Int, length(valid_positions))
    copyto!(valid_positions_native, valid_positions)
    style = AdaptiveOpticsSim.Backends.execution_style(result)

    guide_xy, t_guides = _time_phase() do
        gx, gy = AdaptiveOpticsSim.Tomography._guide_star_grids(
            backend,
            sampling,
            support_diameter_m,
            rotations_rad,
            offset_fractions_x,
            offset_fractions_y,
        )
        _sync_backend!(gx)
        _sync_backend!(gy)
        (gx, gy)
    end
    guide_x, guide_y = guide_xy

    shifted, t_shift = _time_phase() do
        value = AdaptiveOpticsSim.Tomography._scaled_shifted_coord_stack(
            backend,
            guide_x,
            guide_y,
            directions,
            slant_ranges_m,
            source_height_m,
        )
        _sync_backend!(value)
    end
    shifted_flat = reshape(shifted, :, n_lgs, length(slant_ranges_m))

    cst, var_term, inv_L0 =
        AdaptiveOpticsSim.Tomography._covariance_constants(
            r0, atmosphere.L0)
    block = AdaptiveOpticsSim.Calibration._backend_array(
        B, T, n_valid, n_valid)
    fractional_cn2_native = AdaptiveOpticsSim.Calibration.materialize_build(
        backend, atmosphere.fractional_cn2)

    t_selected = 0
    t_scatter = 0

    for jgs in 1:n_lgs
        for igs in 1:jgs
            _, dt_selected = _time_phase() do
                AdaptiveOpticsSim.Backends.launch_kernel_async!(
                    style,
                    AdaptiveOpticsSim.Tomography.selected_covariance_block_kernel!,
                    block,
                    shifted_flat,
                    valid_positions_native,
                    fractional_cn2_native,
                    igs,
                    jgs,
                    cst,
                    var_term,
                    inv_L0,
                    n_valid,
                    length(slant_ranges_m);
                    ndrange=size(block),
                )
                _sync_backend!(block)
            end
            t_selected += dt_selected

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
    println("  backend: ", string(something(AdaptiveOpticsSim.Backends.gpu_backend_name(B), B)))
    println("  case: medium")
    println("  guide_grids_ns: ", t_guides)
    println("  scaled_shifted_coords_ns: ", t_shift)
    println("  selected_block_kernel_ns: ", t_selected)
    println("  block_scatter_ns: ", t_scatter)
    println("  total_timed_ns: ", t_guides + t_shift + t_selected + t_scatter)
    println("  result_type: ", typeof(result))
    return nothing
end

_profile_case(AdaptiveOpticsSim.Backends.CUDABackendTag)
