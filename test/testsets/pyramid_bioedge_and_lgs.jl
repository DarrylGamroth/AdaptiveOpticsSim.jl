@testset "LGS convolution normalization" begin
    n = 8
    expected = reshape(collect(range(0.25, 2.0; length=n * n)), n, n)
    identity_kernel_fft = ones(ComplexF64, n, n)

    fft_buffer = zeros(ComplexF64, n, n)
    fft_plan = AdaptiveOpticsSim.plan_fft_backend!(fft_buffer)
    ifft_plan = AdaptiveOpticsSim.plan_ifft_backend!(fft_buffer)
    same_buffer = copy(expected)
    AdaptiveOpticsSim.apply_lgs_convolution!(
        same_buffer, identity_kernel_fft, fft_buffer, fft_plan, ifft_plan)
    @test same_buffer ≈ expected rtol=1e-12 atol=1e-12
    @test sum(same_buffer) ≈ sum(expected) rtol=1e-12

    split_fft_buffer = zeros(ComplexF64, n, n)
    split_ifft_buffer = similar(split_fft_buffer)
    split_fft_plan = AdaptiveOpticsSim.plan_fft_backend!(split_fft_buffer)
    split_ifft_plan = AdaptiveOpticsSim.plan_ifft_backend!(split_ifft_buffer)
    split_buffer = copy(expected)
    AdaptiveOpticsSim.apply_lgs_convolution!(split_buffer, identity_kernel_fft,
        split_fft_buffer, split_fft_plan, split_ifft_buffer, split_ifft_plan)
    @test split_buffer ≈ expected rtol=1e-12 atol=1e-12
    @test sum(split_buffer) ≈ sum(expected) rtol=1e-12

    expected_stack = cat(expected, reverse(expected; dims=2); dims=3)
    identity_kernel_stack_fft = ones(ComplexF64, n, n, 2)

    scalar_fft_stack = zeros(ComplexF64, n, n, 2)
    scalar_fft_plan = AdaptiveOpticsSim.plan_fft_backend!(scalar_fft_stack, (1, 2))
    scalar_ifft_plan = AdaptiveOpticsSim.plan_ifft_backend!(scalar_fft_stack, (1, 2))
    scalar_stack = copy(expected_stack)
    AdaptiveOpticsSim.apply_lgs_convolution_stack!(scalar_stack,
        identity_kernel_stack_fft, scalar_fft_stack, scalar_fft_plan, scalar_ifft_plan)
    @test scalar_stack ≈ expected_stack rtol=1e-12 atol=1e-12
    @test vec(sum(scalar_stack; dims=(1, 2))) ≈ vec(sum(expected_stack; dims=(1, 2))) rtol=1e-12

    accelerator_fft_stack = zeros(ComplexF64, n, n, 2)
    accelerator_fft_plan = AdaptiveOpticsSim.plan_fft_backend!(accelerator_fft_stack, (1, 2))
    accelerator_ifft_plan = AdaptiveOpticsSim.plan_ifft_backend!(accelerator_fft_stack, (1, 2))
    accelerator_stack = copy(expected_stack)
    AdaptiveOpticsSim._apply_lgs_convolution_stack!(KA_CPU_STYLE, accelerator_stack,
        identity_kernel_stack_fft, accelerator_fft_stack, accelerator_fft_plan, accelerator_ifft_plan)
    @test accelerator_stack ≈ scalar_stack rtol=1e-12 atol=1e-12
    @test vec(sum(accelerator_stack; dims=(1, 2))) ≈
          vec(sum(expected_stack; dims=(1, 2))) rtol=1e-12
end

@testset "Sodium-profile LGS kernel cache invalidation" begin
    tel = Telescope(resolution=16, diameter=8.0,
        central_obstruction=0.0)
    signature_source = LGSSource(
        na_profile=[80000.0 90000.0 100000.0; 0.2 0.6 0.2],
        laser_coordinates=(1.0, -0.5),
        fwhm_spot_up=0.8,
        photon_irradiance=1.0,
    )
    signature_args = (tel, signature_source, 16, 4, 0.1, Float64)
    base_signature = AdaptiveOpticsSim.lgs_kernel_signature(
        signature_args...; model=:subaperture_average)
    @test AdaptiveOpticsSim.lgs_kernel_signature(
        signature_args...; model=:per_subaperture) != base_signature
    @test AdaptiveOpticsSim.lgs_kernel_signature(
        signature_args...; model=:per_subaperture,
        threshold=0.2) != AdaptiveOpticsSim.lgs_kernel_signature(
        signature_args...; model=:per_subaperture, threshold=0.1)
    @test AdaptiveOpticsSim.lgs_kernel_signature(
        tel, signature_source, 32, 4, 0.1, Float64;
        model=:subaperture_average) != base_signature
    @test AdaptiveOpticsSim.lgs_kernel_signature(
        tel, signature_source, 16, 8, 0.1, Float64;
        model=:subaperture_average) != base_signature
    @test AdaptiveOpticsSim.lgs_kernel_signature(
        tel, signature_source, 16, 4, 0.2, Float64;
        model=:subaperture_average) != base_signature
    @test AdaptiveOpticsSim.lgs_kernel_signature(
        tel, signature_source, 16, 4, 0.1, Float32;
        model=:subaperture_average) != base_signature
    @test AdaptiveOpticsSim.lgs_kernel_signature(
        Telescope(resolution=16, diameter=10.0), signature_source,
        16, 4, 0.1, Float64;
        model=:subaperture_average) != base_signature
    changed_wavelength = LGSSource(
        wavelength=600e-9,
        na_profile=signature_source.params.na_profile,
        laser_coordinates=signature_source.params.laser_coordinates,
        fwhm_spot_up=signature_source.params.fwhm_spot_up,
        photon_irradiance=1.0,
    )
    @test AdaptiveOpticsSim.lgs_kernel_signature(
        tel, changed_wavelength, 16, 4, 0.1, Float64;
        model=:subaperture_average) != base_signature

    for family in (:shack_hartmann, :pyramid, :bioedge)
        src = LGSSource(
            na_profile=[80000.0 90000.0 100000.0; 0.2 0.6 0.2],
            laser_coordinates=(1.0, -0.5),
            fwhm_spot_up=0.8,
            photon_irradiance=1.0,
        )
        wfs = if family === :shack_hartmann
            sensor = ShackHartmannWFS(tel; n_lenslets=4,
                mode=Diffractive(), n_pix_subap=4)
            prepare_sampling!(sensor, tel, src)
            sensor
        elseif family === :pyramid
            PyramidWFS(tel; pupil_samples=4, mode=Diffractive(),
                modulation=0.0)
        else
            BioEdgeWFS(tel; pupil_samples=4, mode=Diffractive(),
                modulation=0.0)
        end
        ensure_kernel! = family === :shack_hartmann ?
            AdaptiveOpticsSim.ensure_lgs_kernels! :
            AdaptiveOpticsSim.ensure_lgs_kernel!
        kernel_state = family === :shack_hartmann ?
            wfs.front_end.propagation : wfs.front_end.propagation

        ensure_kernel!(wfs, tel, src)
        original_tag = kernel_state.lgs_kernel_tag
        original_kernel = copy(kernel_state.lgs_kernel_fft)
        ensure_kernel!(wfs, tel, src)
        @test kernel_state.lgs_kernel_tag == original_tag

        src.params.na_profile[2, :] .= [0.8, 0.1, 0.1]
        ensure_kernel!(wfs, tel, src)
        @test kernel_state.lgs_kernel_tag != original_tag
        @test !isapprox(kernel_state.lgs_kernel_fft, original_kernel;
            rtol=1e-12, atol=1e-14)
    end
end

function pyramid_signal_allocation_bytes(wfs, tel, frame, src, scale)
    pyramid_signal!(wfs, tel, frame, src, scale)
    return @allocated pyramid_signal!(wfs, tel, frame, src, scale)
end

function detector_calibration_signature_allocation_bytes(det, sig)
    detector_calibration_signature(det, sig)
    return @allocated detector_calibration_signature(det, sig)
end

@testset "Pyramid, BioEdge, and LGS" begin
    tel = Telescope(resolution=32, diameter=8.0, central_obstruction=0.0)
    for i in 1:tel.params.resolution, j in 1:tel.params.resolution
        tel.state.opd[i, j] = i
    end

    pyr = PyramidWFS(tel; pupil_samples=4, modulation=1.0)
    pyr_slopes = measure!(pyr, tel)
    @test length(pyr_slopes) == 2 * 4 * 4

    bio = BioEdgeWFS(tel; pupil_samples=4)
    bio_slopes = measure!(bio, tel)
    @test length(bio_slopes) == 2 * 4 * 4

    ngs = Source(band=:I, magnitude=0.0)

    pyr_direct = PyramidWFS(tel; pupil_samples=2, mode=Diffractive(), modulation=0.0)
    pyr_direct.acquisition.state.nominal_detector_resolution = 4
    AdaptiveOpticsSim.resize_pyramid_signal_buffers!(pyr_direct, 4)
    pyr_direct.estimator.state.valid_i4q .= Bool[1 0; 1 1]
    AdaptiveOpticsSim.update_pyramid_valid_signal!(pyr_direct)
    pyr_direct.estimator.state.valid_signal_indices = Int[]
    pyr_direct.estimator.state.valid_signal_indices_host = Int[]
    @test AdaptiveOpticsSim.update_pyramid_valid_signal_indices!(pyr_direct) == 3
    AdaptiveOpticsSim.resize_pyramid_slope_buffers!(pyr_direct)
    fill!(pyr_direct.estimator.state.reference_signal_2d, 0.0)
    pyr_frame = [4.0 4.0 1.0 1.0;
                 4.0 4.0 1.0 1.0;
                 3.0 3.0 2.0 2.0;
                 3.0 3.0 2.0 2.0]
    pyr_direct_slopes = AdaptiveOpticsSim.pyramid_signal!(pyr_direct, tel, pyr_frame)
    @test length(pyr_direct_slopes) == 6
    @test pyr_direct_slopes[1:3] ≈ fill(0.4, 3)
    @test pyr_direct_slopes[4:6] ≈ zeros(3)
    copyto!(pyr_direct.acquisition.state.camera_frame, pyr_frame)
    @test AdaptiveOpticsSim.pyramid_slopes!(pyr_direct, tel) ≈ pyr_direct_slopes
    @test AdaptiveOpticsSim.pyramid_slopes!(pyr_direct, tel, pyr_frame) ≈ pyr_direct_slopes
    pyr_direct_accel = PyramidWFS(tel; pupil_samples=2, mode=Diffractive(), modulation=0.0)
    pyr_direct_accel.acquisition.state.nominal_detector_resolution = 4
    AdaptiveOpticsSim.resize_pyramid_signal_buffers!(pyr_direct_accel, 4)
    pyr_direct_accel.estimator.state.valid_i4q .= pyr_direct.estimator.state.valid_i4q
    AdaptiveOpticsSim.update_pyramid_valid_signal!(pyr_direct_accel)
    @test AdaptiveOpticsSim.update_pyramid_valid_signal_indices!(pyr_direct_accel) == 3
    AdaptiveOpticsSim.resize_pyramid_slope_buffers!(pyr_direct_accel)
    fill!(pyr_direct_accel.estimator.state.reference_signal_2d, 0.0)
    @test AdaptiveOpticsSim.pyramid_signal!(KA_CPU_STYLE, pyr_direct_accel, tel, pyr_frame, nothing) ≈ pyr_direct_slopes

    pyr_unit_mask = PyramidWFS(tel; pupil_samples=2, mode=Diffractive(),
        modulation=0.0, mask_scale=1.0)
    pyr_scaled_mask = PyramidWFS(tel; pupil_samples=2, mode=Diffractive(),
        modulation=0.0, mask_scale=1.5)
    @test pyr_scaled_mask.front_end.propagation.pyramid_mask != pyr_unit_mask.front_end.propagation.pyramid_mask
    pyr_old_unit_mask = PyramidWFS(tel; pupil_samples=2, mode=Diffractive(),
        modulation=0.0, old_mask=true, mask_scale=1.0)
    pyr_old_scaled_mask = PyramidWFS(tel; pupil_samples=2,
        mode=Diffractive(), modulation=0.0, old_mask=true, mask_scale=1.5)
    @test pyr_old_scaled_mask.front_end.propagation.pyramid_mask !=
          pyr_old_unit_mask.front_end.propagation.pyramid_mask
    for invalid_scale in (0.0, -1.0, NaN, Inf)
        @test_throws InvalidConfiguration PyramidWFS(tel; pupil_samples=2,
            mode=Diffractive(), mask_scale=invalid_scale)
    end

    # Separated pupil images occupy the configured centered regions, not the
    # outer corners of a padded detector frame.
    @test_throws InvalidConfiguration PyramidWFS(tel; pupil_samples=2,
        mode=Diffractive(), n_pix_separation=-1)
    @test_throws InvalidConfiguration PyramidWFS(tel; pupil_samples=2,
        mode=Diffractive(), n_pix_separation=2, n_pix_edge=-1)
    @test_throws InvalidConfiguration PyramidWFS(tel; pupil_samples=2,
        mode=Diffractive(), n_pix_separation=2, n_pix_edge=1, binning=2)
    pyr_separated = PyramidWFS(tel; pupil_samples=2, mode=Diffractive(),
        modulation=0.0, n_pix_separation=2, n_pix_edge=1)
    pyr_separated.acquisition.state.nominal_detector_resolution = 8
    @test_throws InvalidConfiguration AdaptiveOpticsSim.resize_pyramid_signal_buffers!(
        pyr_separated, 0)
    @test_throws InvalidConfiguration AdaptiveOpticsSim.resize_pyramid_signal_buffers!(
        pyr_separated, 4)
    AdaptiveOpticsSim.resize_pyramid_signal_buffers!(pyr_separated, 8)
    @test_throws DimensionMismatchError AdaptiveOpticsSim.pyramid_signal!(
        pyr_separated, tel, zeros(0, 0))
    fill!(pyr_separated.estimator.state.valid_i4q, true)
    AdaptiveOpticsSim.update_pyramid_valid_signal!(pyr_separated)
    @test AdaptiveOpticsSim.update_pyramid_valid_signal_indices!(pyr_separated) == 4
    AdaptiveOpticsSim.resize_pyramid_slope_buffers!(pyr_separated)
    fill!(pyr_separated.estimator.state.reference_signal_2d, 0.0)
    separated_frame = zeros(8, 8)
    separated_frame[[1, 8], :] .= 100.0
    separated_frame[:, [1, 8]] .= 100.0
    separated_frame[2:3, 2:3] .= 4.0
    separated_frame[2:3, 6:7] .= 1.0
    separated_frame[6:7, 6:7] .= 2.0
    separated_frame[6:7, 2:3] .= 3.0
    separated_slopes = AdaptiveOpticsSim.pyramid_signal!(
        pyr_separated, tel, separated_frame)
    @test separated_slopes[1:4] ≈ fill(0.4, 4)
    @test separated_slopes[5:8] ≈ zeros(4)

    pyr_separated_accel = PyramidWFS(tel; pupil_samples=2,
        mode=Diffractive(), modulation=0.0, n_pix_separation=2,
        n_pix_edge=1)
    pyr_separated_accel.acquisition.state.nominal_detector_resolution = 8
    AdaptiveOpticsSim.resize_pyramid_signal_buffers!(pyr_separated_accel, 8)
    fill!(pyr_separated_accel.estimator.state.valid_i4q, true)
    AdaptiveOpticsSim.update_pyramid_valid_signal!(pyr_separated_accel)
    @test AdaptiveOpticsSim.update_pyramid_valid_signal_indices!(
        pyr_separated_accel) == 4
    AdaptiveOpticsSim.resize_pyramid_slope_buffers!(pyr_separated_accel)
    fill!(pyr_separated_accel.estimator.state.reference_signal_2d, 0.0)
    @test AdaptiveOpticsSim.pyramid_signal!(KA_CPU_STYLE,
        pyr_separated_accel, tel, separated_frame, nothing) ≈ separated_slopes

    zero_slopes = fill(1.0, 8)
    AdaptiveOpticsSim._pyramid_slopes!(AdaptiveOpticsSim.ScalarCPUStyle(), zero_slopes, zeros(4, 4), trues(2, 2),
        2, 2, 4, 4, 0, 0, 0, 0, 0, 0, 0, 0, (0, 0, 0, 0), (0, 0, 0, 0))
    @test all(iszero, zero_slopes)
    invalid_slopes = fill(1.0, 8)
    AdaptiveOpticsSim._pyramid_slopes!(AdaptiveOpticsSim.ScalarCPUStyle(), invalid_slopes, ones(4, 4), falses(2, 2),
        2, 2, 4, 4, 0, 0, 0, 0, 0, 0, 0, 0, (0, 0, 0, 0), (0, 0, 0, 0))
    @test all(iszero, invalid_slopes)
    AdaptiveOpticsSim.apply_shift_wfs!(pyr_direct; sx=1.2, sy=-1.8)
    @test pyr_direct.estimator.state.shift_x == (1, 1, 1, 1)
    @test pyr_direct.estimator.state.shift_y == (-2, -2, -2, -2)
    shift_x, shift_y = AdaptiveOpticsSim.pyramid_shift_components([0, 1, 2, 3], [3, 2, 1, 0])
    @test shift_x == (0, 1, 2, 3)
    @test shift_y == (3, 2, 1, 0)
    @test_throws InvalidConfiguration AdaptiveOpticsSim.pyramid_shift_components([1, 2, 3], [1, 2, 3, 4])
    AdaptiveOpticsSim.set_optical_gain!(pyr_direct, 2.0)
    @test all(==(2.0), pyr_direct.estimator.state.optical_gain)
    AdaptiveOpticsSim.set_optical_gain!(pyr_direct, collect(1.0:length(pyr_direct.estimator.state.optical_gain)))
    @test pyr_direct.estimator.state.optical_gain[end] == length(pyr_direct.estimator.state.optical_gain)
    @test_throws InvalidConfiguration AdaptiveOpticsSim.set_optical_gain!(pyr_direct, [1.0])

    pyr_invalid = PyramidWFS(tel; pupil_samples=2, mode=Diffractive(), modulation=0.0)
    pyr_invalid.acquisition.state.nominal_detector_resolution = 4
    AdaptiveOpticsSim.resize_pyramid_signal_buffers!(pyr_invalid, 4)
    fill!(pyr_invalid.estimator.state.valid_i4q, false)
    AdaptiveOpticsSim.update_pyramid_valid_signal!(pyr_invalid)
    @test AdaptiveOpticsSim.update_pyramid_valid_signal_indices!(pyr_invalid) == 0
    @test_throws InvalidConfiguration AdaptiveOpticsSim.resize_pyramid_slope_buffers!(pyr_invalid)

    pyr_incidence = PyramidWFS(tel; pupil_samples=2, mode=Diffractive(), modulation=0.0,
        normalization=IncidenceFluxNormalization())
    expected_pyr_norm = AdaptiveOpticsSim.photon_irradiance(ngs) *
                        (tel.params.diameter /
                         pyr_incidence.estimator.params.pupil_samples)^2
    @test AdaptiveOpticsSim.pyramid_normalization(
        pyr_incidence.estimator.params.normalization,
        pyr_incidence, tel, ngs, 3, 10.0) ≈ expected_pyr_norm
    @test AdaptiveOpticsSim.pyramid_normalization(
        pyr_incidence.estimator.params.normalization,
        pyr_incidence, tel, nothing, 3, 10.0) == 1.0

    bio_direct = BioEdgeWFS(tel; pupil_samples=2, mode=Diffractive())
    bio_direct.acquisition.state.nominal_detector_resolution = 4
    AdaptiveOpticsSim.resize_bioedge_signal_buffers!(bio_direct, 4)
    bio_direct.estimator.state.valid_i4q .= Bool[1 0; 1 1]
    AdaptiveOpticsSim.update_bioedge_valid_signal!(bio_direct)
    bio_direct.estimator.state.valid_signal_indices = Int[]
    bio_direct.estimator.state.valid_signal_indices_host = Int[]
    @test AdaptiveOpticsSim.update_bioedge_valid_signal_indices!(bio_direct) == 3
    AdaptiveOpticsSim.resize_bioedge_slope_buffers!(bio_direct)
    fill!(bio_direct.estimator.state.reference_signal_2d, 0.0)
    fill!(bio_direct.estimator.state.optical_gain, 2.0)
    bio_frame = copy(pyr_frame)
    bio_direct_slopes = AdaptiveOpticsSim.bioedge_signal!(bio_direct, tel, bio_frame)
    @test length(bio_direct_slopes) == 6
    @test bio_direct_slopes[1:3] ≈ fill(0.8, 3)
    @test bio_direct_slopes[4:6] ≈ zeros(3)
    @test AdaptiveOpticsSim.bioedge_slopes_intensity!(bio_direct, tel, bio_frame) ≈ bio_direct_slopes
    bio_phase = reshape(collect(range(0.0, 1.0; length=32 * 32)), 32, 32)
    bio_edge_mask = falses(32, 32)
    bio_edge_mask[1, :] .= true
    bio_edge_mask[end, :] .= true
    bio_edge_mask[:, 1] .= true
    bio_edge_mask[:, end] .= true
    @test length(AdaptiveOpticsSim.bioedge_slopes!(bio, bio_phase, bio_edge_mask)) == 2 * 4 * 4
    bio_direct_accel = BioEdgeWFS(tel; pupil_samples=2, mode=Diffractive())
    bio_direct_accel.acquisition.state.nominal_detector_resolution = 4
    AdaptiveOpticsSim.resize_bioedge_signal_buffers!(bio_direct_accel, 4)
    bio_direct_accel.estimator.state.valid_i4q .= bio_direct.estimator.state.valid_i4q
    AdaptiveOpticsSim.update_bioedge_valid_signal!(bio_direct_accel)
    @test AdaptiveOpticsSim.update_bioedge_valid_signal_indices!(bio_direct_accel) == 3
    AdaptiveOpticsSim.resize_bioedge_slope_buffers!(bio_direct_accel)
    fill!(bio_direct_accel.estimator.state.reference_signal_2d, 0.0)
    fill!(bio_direct_accel.estimator.state.optical_gain, 2.0)
    @test AdaptiveOpticsSim.bioedge_signal!(KA_CPU_STYLE, bio_direct_accel, tel, bio_frame, nothing) ≈ bio_direct_slopes
    AdaptiveOpticsSim.set_optical_gain!(bio_direct, 3.0)
    @test all(==(3.0), bio_direct.estimator.state.optical_gain)
    AdaptiveOpticsSim.set_optical_gain!(bio_direct, collect(1.0:length(bio_direct.estimator.state.optical_gain)))
    @test bio_direct.estimator.state.optical_gain[end] == length(bio_direct.estimator.state.optical_gain)
    @test_throws InvalidConfiguration AdaptiveOpticsSim.set_optical_gain!(bio_direct, [1.0])

    bio_invalid = BioEdgeWFS(tel; pupil_samples=2, mode=Diffractive())
    bio_invalid.acquisition.state.nominal_detector_resolution = 4
    AdaptiveOpticsSim.resize_bioedge_signal_buffers!(bio_invalid, 4)
    fill!(bio_invalid.estimator.state.valid_i4q, false)
    AdaptiveOpticsSim.update_bioedge_valid_signal!(bio_invalid)
    @test AdaptiveOpticsSim.update_bioedge_valid_signal_indices!(bio_invalid) == 0
    @test_throws InvalidConfiguration AdaptiveOpticsSim.resize_bioedge_slope_buffers!(bio_invalid)

    bio_incidence = BioEdgeWFS(tel; pupil_samples=2, mode=Diffractive(),
        normalization=IncidenceFluxNormalization())
    expected_bio_norm = AdaptiveOpticsSim.photon_irradiance(ngs) *
                        (tel.params.diameter /
                         bio_incidence.estimator.params.pupil_samples)^2
    @test AdaptiveOpticsSim.bioedge_normalization(
        bio_incidence.estimator.params.normalization,
        bio_incidence, tel, ngs, 3, 10.0) ≈ expected_bio_norm
    @test AdaptiveOpticsSim.bioedge_normalization(
        bio_incidence.estimator.params.normalization,
        bio_incidence, tel, nothing, 3, 10.0) == 1.0
    bio_flux_select = BioEdgeWFS(tel; pupil_samples=2, mode=Diffractive(), light_ratio=0.25)
    @test AdaptiveOpticsSim.select_bioedge_valid_i4q!(
        AdaptiveOpticsSim.ScalarCPUStyle(), bio_flux_select, tel, ngs) === bio_flux_select
    @test bio_flux_select.estimator.state.valid_signal_count > 0
    bio_flux_select_accel = BioEdgeWFS(tel; pupil_samples=2, mode=Diffractive(), light_ratio=0.25)
    @test AdaptiveOpticsSim.select_bioedge_valid_i4q!(
        KA_CPU_STYLE, bio_flux_select_accel, tel, ngs) === bio_flux_select_accel
    @test bio_flux_select_accel.estimator.state.valid_signal_count > 0
    sh = ShackHartmannWFS(tel; n_lenslets=4)
    lgs = LGSSource(elongation_factor=2.0)
    slopes_ngs = measure!(sh, tel, ngs)
    slopes_lgs = measure!(sh, tel, lgs)
    n = microlens_array(sh).params.n_lenslets * microlens_array(sh).params.n_lenslets
    @test slopes_lgs[n+1:end] ≈ slopes_ngs[n+1:end] .* 2.0

    bio_lgs = BioEdgeWFS(tel; pupil_samples=4, mode=Diffractive())
    @test AdaptiveOpticsSim.ensure_lgs_kernel!(bio_lgs, tel, lgs) === bio_lgs
    na_profile = [80000.0 90000.0 100000.0; 0.2 0.6 0.2]
    lgs_profile = LGSSource(elongation_factor=1.2, na_profile=na_profile, fwhm_spot_up=1.0)
    bio_lgs_profile = BioEdgeWFS(tel; pupil_samples=4, mode=Diffractive())
    @test AdaptiveOpticsSim.ensure_lgs_kernel!(bio_lgs_profile, tel, lgs_profile) === bio_lgs_profile
    cached_tag = bio_lgs_profile.front_end.propagation.lgs_kernel_tag
    @test AdaptiveOpticsSim.ensure_lgs_kernel!(bio_lgs_profile, tel, lgs_profile) === bio_lgs_profile
    @test bio_lgs_profile.front_end.propagation.lgs_kernel_tag == cached_tag
end

@testset "Diffractive WFS" begin
    tel = Telescope(resolution=32, diameter=8.0, central_obstruction=0.0)
    for i in 1:tel.params.resolution, j in 1:tel.params.resolution
        tel.state.opd[i, j] = i
    end
    ngs = Source(band=:I, magnitude=0.0)
    lgs = LGSSource(elongation_factor=1.5, photon_irradiance=1.0)

    sh = ShackHartmannWFS(tel; n_lenslets=4, mode=Diffractive())
    @test_throws InvalidConfiguration measure!(sh, tel)
    sh_slopes = measure!(sh, tel, ngs)
    @test length(sh_slopes) == 2 * 4 * 4
    @test all(isfinite, sh_slopes)
    sh_lgs = measure!(sh, tel, lgs)
    @test all(isfinite, sh_lgs)

    na_profile = [80000.0 90000.0 100000.0; 0.2 0.6 0.2]
    lgs_profile = LGSSource(elongation_factor=1.2, na_profile=na_profile,
        fwhm_spot_up=1.0, photon_irradiance=1.0)
    sh_profile = ShackHartmannWFS(tel; n_lenslets=4, mode=Diffractive())
    sh_profile_slopes = measure!(sh_profile, tel, lgs_profile)
    @test all(isfinite, sh_profile_slopes)

    sh_sampled = ShackHartmannWFS(tel; n_lenslets=4, mode=Diffractive(), pixel_scale_arcsec=0.06, n_pix_subap=8)
    sh_sampled_slopes = measure!(sh_sampled, tel, ngs)
    @test length(sh_sampled_slopes) == 2 * 4 * 4

    pyr_sampled = PyramidWFS(tel; pupil_samples=4, mode=Diffractive(),
        n_pix_separation=4, n_pix_edge=2, binning=2)
    pyr_sampled_slopes = measure!(pyr_sampled, tel, ngs)
    @test length(pyr_sampled_slopes) == 2 * count(pyr_sampled.estimator.state.valid_i4q)
    pyr_intensity = reshape(
        Float64.(1:length(pyr_sampled.front_end.propagation.intensity)),
        size(pyr_sampled.front_end.propagation.intensity),
    )
    pyr_frame = copy(AdaptiveOpticsSim.sample_pyramid_intensity!(pyr_sampled, tel, pyr_intensity))
    pyr_camera = zeros(Float64, 16, 16)
    pyr_manual = zeros(Float64, 8, 8)
    AdaptiveOpticsSim.bin2d!(pyr_camera, pyr_intensity, 8)
    AdaptiveOpticsSim.bin2d!(pyr_manual, pyr_camera, 2)
    @test pyr_frame == pyr_manual

    bio_sampled = BioEdgeWFS(tel; pupil_samples=4, mode=Diffractive(), binning=2)
    bio_sampled_slopes = measure!(bio_sampled, tel, ngs)
    @test length(bio_sampled_slopes) == 2 * count(bio_sampled.estimator.state.valid_i4q)
    bio_intensity = reshape(Float64.(1:size(tel.state.opd, 1)^2), size(tel.state.opd))
    bio_frame = copy(AdaptiveOpticsSim.sample_bioedge_intensity!(bio_sampled, tel, bio_intensity))
    bio_camera = zeros(Float64, 4, 4)
    bio_manual = similar(bio_frame)
    AdaptiveOpticsSim.bin2d!(bio_camera, bio_intensity, 8)
    AdaptiveOpticsSim.bin2d!(bio_manual, bio_camera, div(size(bio_camera, 1), size(bio_frame, 1)))
    @test bio_frame == bio_manual

    pyr_profile = PyramidWFS(tel; pupil_samples=4, mode=Diffractive())
    pyr_profile_slopes = measure!(pyr_profile, tel, lgs_profile)
    @test all(isfinite, pyr_profile_slopes)

    bio_profile = BioEdgeWFS(tel; pupil_samples=4, mode=Diffractive())
    bio_profile_slopes = measure!(bio_profile, tel, lgs_profile)
    @test all(isfinite, bio_profile_slopes)

    pyr = PyramidWFS(tel; pupil_samples=4, mode=Diffractive())
    @test_throws InvalidConfiguration measure!(pyr, tel)
    pyr_slopes = measure!(pyr, tel, ngs)
    @test length(pyr_slopes) == 2 * 4 * 4

    bio = BioEdgeWFS(tel; pupil_samples=4, mode=Diffractive())
    @test_throws InvalidConfiguration measure!(bio, tel)
    bio_slopes = measure!(bio, tel, ngs)
    @test length(bio_slopes) == 2 * 4 * 4

    det = Detector(noise=NoiseNone(), binning=1)
    sh_det = ShackHartmannWFS(tel; n_lenslets=4, mode=Diffractive())
    sh_det_slopes = measure!(sh_det, tel, ngs, det)
    @test length(sh_det_slopes) == 2 * 4 * 4
    sh_det_image = wfs_detector_image(sh_det, det; gap=1)
    @test ndims(sh_det_image) == 2
    sh_adu_det = Detector(noise=NoiseNone(), binning=1, full_well=30_000.0, bits=12, output_type=UInt16)
    sh_adu = ShackHartmannWFS(tel; n_lenslets=4, mode=Diffractive())
    measure!(sh_adu, tel, ngs, sh_adu_det; rng=MersenneTwister(15))
    sh_adu_image = wfs_detector_image(sh_adu, sh_adu_det; gap=1)
    @test sh_adu_image isa Matrix{UInt16}
    @test maximum(sh_adu_image) <= 0x0fff
    pyr_det = PyramidWFS(tel; pupil_samples=4, mode=Diffractive())
    pyr_det_slopes = measure!(pyr_det, tel, ngs, det)
    @test length(pyr_det_slopes) == 2 * 4 * 4
    @test wfs_detector_image(pyr_det, det) === output_frame(det)
    bio_det = BioEdgeWFS(tel; pupil_samples=4, mode=Diffractive())
    bio_det_slopes = measure!(bio_det, tel, ngs, det)
    @test length(bio_det_slopes) == 2 * 4 * 4
    @test wfs_detector_image(bio_det, det) === output_frame(det)

    ast = Asterism([ngs, Source(band=:I, magnitude=0.0, coordinates=(0.0, 0.0))])
    sh_ast = ShackHartmannWFS(tel; n_lenslets=4, mode=Diffractive())
    sh_ast_slopes = copy(measure!(sh_ast, tel, ast))
    @test length(sh_ast_slopes) == 2 * 4 * 4
    sh_ast_serial = ShackHartmannWFS(tel; n_lenslets=4, mode=Diffractive())
    AdaptiveOpticsSim.prepare_sampling!(sh_ast_serial, tel, ast.sources[1])
    AdaptiveOpticsSim.ensure_sh_calibration!(sh_ast_serial, tel, ast.sources[1])
    fill!(sh_ast_serial.acquisition.detector_noise_cube, zero(eltype(sh_ast_serial.acquisition.detector_noise_cube)))
    for src in ast.sources
        AdaptiveOpticsSim.sampled_spots_peak!(sh_ast_serial, tel, src)
        sh_ast_serial.acquisition.detector_noise_cube .+= sh_ast_serial.acquisition.spot_cube
    end
    copyto!(sh_ast_serial.acquisition.spot_cube, sh_ast_serial.acquisition.detector_noise_cube)
    sh_ast_serial_peak = maximum(sh_ast_serial.acquisition.spot_cube)
    AdaptiveOpticsSim.sh_signal_from_spots!(sh_ast_serial, sh_ast_serial_peak, slope_extraction_model(sh_ast_serial))
    AdaptiveOpticsSim.subtract_reference_and_scale!(sh_ast_serial)
    sh_ast_serial_slopes = copy(slopes(sh_ast_serial))
    @test norm(sh_ast_slopes - sh_ast_serial_slopes) / norm(sh_ast_slopes) < 0.07
    mixed_ngs = Source(wavelength=wavelength(lgs_profile), magnitude=0.0, coordinates=(0.0, 0.0))
    mixed_ast = Asterism([mixed_ngs, lgs_profile])
    sh_mixed_det = ShackHartmannWFS(tel; n_lenslets=4, mode=Diffractive())
    @test_throws InvalidConfiguration measure!(sh_mixed_det, tel,
        mixed_ast, det; rng=MersenneTwister(14))
    @test !sh_mixed_det.calibration.calibrated
    pyr_ast = PyramidWFS(tel; pupil_samples=4, mode=Diffractive())
    pyr_ast_slopes = copy(measure!(pyr_ast, tel, ast))
    @test length(pyr_ast_slopes) == 2 * 4 * 4
    pyr_ast_serial = PyramidWFS(tel; pupil_samples=4, mode=Diffractive())
    AdaptiveOpticsSim.ensure_pyramid_calibration!(pyr_ast_serial, tel, ast.sources[1])
    pyr_ast_stack = @view AdaptiveOpticsSim.ensure_pyramid_asterism_stack!(pyr_ast_serial, length(ast.sources))[:, :, 1:length(ast.sources)]
    fill!(pyr_ast_serial.front_end.propagation.intensity, zero(eltype(pyr_ast_serial.front_end.propagation.intensity)))
    for (src_idx, src) in pairs(ast.sources)
        AdaptiveOpticsSim.pyramid_intensity!(@view(pyr_ast_stack[:, :, src_idx]), pyr_ast_serial, tel, src)
        pyr_ast_serial.front_end.propagation.intensity .+= @view(pyr_ast_stack[:, :, src_idx])
    end
    pyr_ast_intensity = AdaptiveOpticsSim.sample_pyramid_intensity!(pyr_ast_serial, tel, pyr_ast_serial.front_end.propagation.intensity)
    AdaptiveOpticsSim.pyramid_signal!(pyr_ast_serial, tel, pyr_ast_intensity)
    slopes(pyr_ast_serial) .*= pyr_ast_serial.estimator.state.optical_gain
    @test pyr_ast_slopes ≈ slopes(pyr_ast_serial)
    bio_ast = BioEdgeWFS(tel; pupil_samples=4, mode=Diffractive())
    bio_ast_slopes = copy(measure!(bio_ast, tel, ast))
    @test length(bio_ast_slopes) == 2 * 4 * 4
    bio_ast_serial = BioEdgeWFS(tel; pupil_samples=4, mode=Diffractive())
    AdaptiveOpticsSim.ensure_bioedge_calibration!(bio_ast_serial, tel, ast.sources[1])
    fill!(bio_ast_serial.acquisition.state.binned_intensity, zero(eltype(bio_ast_serial.acquisition.state.binned_intensity)))
    for src in ast.sources
        AdaptiveOpticsSim.bioedge_intensity!(bio_ast_serial.front_end.propagation.intensity, bio_ast_serial, tel, src)
        bio_ast_serial.acquisition.state.binned_intensity .+= bio_ast_serial.front_end.propagation.intensity
    end
    bio_ast_intensity = AdaptiveOpticsSim.sample_bioedge_intensity!(bio_ast_serial, tel, bio_ast_serial.acquisition.state.binned_intensity)
    AdaptiveOpticsSim.bioedge_signal!(bio_ast_serial, tel, bio_ast_intensity)
    @test bio_ast_slopes ≈ slopes(bio_ast_serial)

    pyr_ast_det = PyramidWFS(tel; pupil_samples=4, mode=Diffractive())
    pyr_ast_det_slopes = copy(measure!(pyr_ast_det, tel, ast, det))
    pyr_ast_det_serial = PyramidWFS(tel; pupil_samples=4, mode=Diffractive())
    AdaptiveOpticsSim.ensure_pyramid_calibration!(pyr_ast_det_serial, tel, ast.sources[1])
    pyr_ast_det_stack = @view AdaptiveOpticsSim.ensure_pyramid_asterism_stack!(pyr_ast_det_serial, length(ast.sources))[:, :, 1:length(ast.sources)]
    fill!(pyr_ast_det_serial.front_end.propagation.intensity, zero(eltype(pyr_ast_det_serial.front_end.propagation.intensity)))
    for (src_idx, src) in pairs(ast.sources)
        AdaptiveOpticsSim.pyramid_intensity!(@view(pyr_ast_det_stack[:, :, src_idx]), pyr_ast_det_serial, tel, src)
        pyr_ast_det_serial.front_end.propagation.intensity .+= @view(pyr_ast_det_stack[:, :, src_idx])
    end
    pyr_ast_det_intensity = AdaptiveOpticsSim.sample_pyramid_intensity!(pyr_ast_det_serial, tel, pyr_ast_det_serial.front_end.propagation.intensity)
    pyr_ast_det_frame = capture!(det, pyr_ast_det_intensity,
        first(ast.sources); rng=MersenneTwister(12))
    AdaptiveOpticsSim.resize_pyramid_signal_buffers!(pyr_ast_det_serial, size(pyr_ast_det_frame, 1))
    AdaptiveOpticsSim.pyramid_signal!(pyr_ast_det_serial, tel, pyr_ast_det_frame)
    slopes(pyr_ast_det_serial) .*= pyr_ast_det_serial.estimator.state.optical_gain
    @test pyr_ast_det_slopes ≈ slopes(pyr_ast_det_serial)

    bio_ast_det = BioEdgeWFS(tel; pupil_samples=4, mode=Diffractive())
    bio_ast_det_slopes = copy(measure!(bio_ast_det, tel, ast, det; rng=MersenneTwister(13)))
    bio_ast_det_serial = BioEdgeWFS(tel; pupil_samples=4, mode=Diffractive())
    AdaptiveOpticsSim.ensure_bioedge_calibration!(bio_ast_det_serial, tel, ast.sources[1])
    fill!(bio_ast_det_serial.acquisition.state.binned_intensity, zero(eltype(bio_ast_det_serial.acquisition.state.binned_intensity)))
    for src in ast.sources
        AdaptiveOpticsSim.bioedge_intensity!(bio_ast_det_serial.front_end.propagation.intensity, bio_ast_det_serial, tel, src)
        bio_ast_det_serial.acquisition.state.binned_intensity .+= bio_ast_det_serial.front_end.propagation.intensity
    end
    bio_ast_det_intensity = AdaptiveOpticsSim.sample_bioedge_intensity!(bio_ast_det_serial, tel, bio_ast_det_serial.acquisition.state.binned_intensity)
    bio_ast_det_frame = capture!(det, bio_ast_det_intensity,
        first(ast.sources); rng=MersenneTwister(13))
    AdaptiveOpticsSim.resize_bioedge_signal_buffers!(bio_ast_det_serial, size(bio_ast_det_frame, 1))
    AdaptiveOpticsSim.bioedge_signal!(bio_ast_det_serial, tel, bio_ast_det_frame)
    @test bio_ast_det_slopes ≈ slopes(bio_ast_det_serial)

    common_qe_wavelength = 550e-9
    common_qe_ast = Asterism([
        Source(band=:custom, wavelength=common_qe_wavelength,
            photon_irradiance=1.0, coordinates=(0.0, 0.0)),
        Source(band=:custom, wavelength=common_qe_wavelength,
            photon_irradiance=2.0, coordinates=(0.5, 90.0)),
    ])
    wavelength_dependent_qe = SampledQuantumEfficiency(
        [500e-9, common_qe_wavelength, 600e-9], [0.1, 0.35, 0.9])

    pyr_sampled_qe = PyramidWFS(tel; pupil_samples=4, mode=Diffractive())
    pyr_sampled_qe_det = Detector(noise=NoiseNone(),
        qe=wavelength_dependent_qe, integration_time=1.0, binning=1)
    measure!(pyr_sampled_qe, tel, common_qe_ast, pyr_sampled_qe_det;
        rng=MersenneTwister(21))
    pyr_sampled_qe_frame = copy(output_frame(pyr_sampled_qe_det))
    pyr_scalar_qe = PyramidWFS(tel; pupil_samples=4, mode=Diffractive())
    pyr_scalar_qe_det = Detector(noise=NoiseNone(), qe=0.35,
        integration_time=1.0, binning=1)
    measure!(pyr_scalar_qe, tel, common_qe_ast, pyr_scalar_qe_det;
        rng=MersenneTwister(21))
    @test sum(pyr_sampled_qe_frame) > 0
    @test pyr_sampled_qe_frame ≈ output_frame(pyr_scalar_qe_det)

    bio_sampled_qe = BioEdgeWFS(tel; pupil_samples=4, mode=Diffractive())
    bio_sampled_qe_det = Detector(noise=NoiseNone(),
        qe=wavelength_dependent_qe, integration_time=1.0, binning=1)
    measure!(bio_sampled_qe, tel, common_qe_ast, bio_sampled_qe_det;
        rng=MersenneTwister(22))
    bio_sampled_qe_frame = copy(output_frame(bio_sampled_qe_det))
    bio_scalar_qe = BioEdgeWFS(tel; pupil_samples=4, mode=Diffractive())
    bio_scalar_qe_det = Detector(noise=NoiseNone(), qe=0.35,
        integration_time=1.0, binning=1)
    measure!(bio_scalar_qe, tel, common_qe_ast, bio_scalar_qe_det;
        rng=MersenneTwister(22))
    @test sum(bio_sampled_qe_frame) > 0
    @test bio_sampled_qe_frame ≈ output_frame(bio_scalar_qe_det)

    mixed_qe_ast = Asterism([
        Source(band=:custom, wavelength=common_qe_wavelength,
            photon_irradiance=1.0),
        Source(band=:custom, wavelength=600e-9,
            photon_irradiance=1.0),
    ])
    @test_throws InvalidConfiguration measure!(
        PyramidWFS(tel; pupil_samples=4, mode=Diffractive()),
        tel, mixed_qe_ast, pyr_sampled_qe_det)
    @test_throws InvalidConfiguration measure!(
        BioEdgeWFS(tel; pupil_samples=4, mode=Diffractive()),
        tel, mixed_qe_ast, bio_sampled_qe_det)
end

@testset "Pyramid and BioEdge incidence-normalization contracts" begin
    tel = Telescope(resolution=32, diameter=8.0, central_obstruction=0.0)
    @inbounds for j in axes(tel.state.opd, 2), i in axes(tel.state.opd, 1)
        tel.state.opd[i, j] = 2e-8 * sinpi(2 * i / 32) * cospi(2 * j / 32)
    end
    src = Source(band=:custom, wavelength=0.75e-6,
        photon_irradiance=10.0)

    for family in (:pyramid, :bioedge)
        make_wfs = () -> family === :pyramid ?
            PyramidWFS(tel; pupil_samples=4, mode=Diffractive(),
                modulation=1.0,
                normalization=IncidenceFluxNormalization()) :
            BioEdgeWFS(tel; pupil_samples=4, mode=Diffractive(),
                modulation=1.0,
                normalization=IncidenceFluxNormalization())

        asymmetric_response = [0.0 0.0 0.0;
                               0.0 0.2 0.8;
                               0.0 0.0 0.0]
        reference = copy(measure!(make_wfs(), tel, src,
            Detector(noise=NoiseNone(), integration_time=1.0, qe=1.0,
                response_model=SampledFrameResponse(asymmetric_response))))
        scaled = copy(measure!(make_wfs(), tel, src,
            Detector(noise=NoiseNone(), integration_time=0.5, qe=0.25,
                response_model=SampledFrameResponse(asymmetric_response))))
        @test norm(reference) > 1e-6
        @test scaled ≈ reference atol=1e-12 rtol=1e-12

        gained = copy(measure!(make_wfs(), tel, src,
            Detector(noise=NoiseNone(), integration_time=1.0, qe=1.0,
                gain=4.0,
                response_model=SampledFrameResponse(asymmetric_response))))
        @test gained ≈ reference atol=1e-12 rtol=1e-12

        hgcdte_reference = copy(measure!(make_wfs(), tel, src,
            Detector(noise=NoiseNone(), integration_time=1.0, qe=1.0,
                gain=1.0,
                sensor=HgCdTeAvalancheArraySensor(avalanche_gain=1.0,
                    sampling_mode=CorrelatedDoubleSampling()),
                response_model=SampledFrameResponse(asymmetric_response))))
        hgcdte_gained = copy(measure!(make_wfs(), tel, src,
            Detector(noise=NoiseNone(), integration_time=1.0, qe=1.0,
                gain=3.0,
                sensor=HgCdTeAvalancheArraySensor(avalanche_gain=2.0,
                    sampling_mode=CorrelatedDoubleSampling()),
                response_model=SampledFrameResponse(asymmetric_response))))
        @test norm(hgcdte_reference) > 1e-6
        @test hgcdte_gained ≈ hgcdte_reference atol=1e-12 rtol=1e-12

        ast = Asterism([src, src])
        single_source = copy(measure!(make_wfs(), tel, src))
        two_sources = copy(measure!(make_wfs(), tel, ast))
        @test norm(single_source) > 1e-6
        @test two_sources ≈ single_source atol=1e-12 rtol=1e-12

        expected_ast_norm = 2 * photon_irradiance(src) *
            (tel.params.diameter / 4)^2
        ast_wfs = make_wfs()
        normalization = ast_wfs.estimator.params.normalization
        actual_ast_norm = family === :pyramid ?
            pyramid_normalization(normalization, ast_wfs, tel, ast, 16, 1.0) :
            bioedge_normalization(normalization, ast_wfs, tel, ast, 16, 1.0)
        @test actual_ast_norm ≈ expected_ast_norm
    end

    zero_src = Source(band=:custom, wavelength=wavelength(src),
        photon_irradiance=0.0)
    for normalization in (MeanValidFluxNormalization(),
            IncidenceFluxNormalization())
        zero_pyramid = PyramidWFS(tel; pupil_samples=4,
            mode=Diffractive(), modulation=0.0,
            normalization=normalization)
        zero_bioedge = BioEdgeWFS(tel; pupil_samples=4,
            mode=Diffractive(), normalization=normalization)
        @test all(iszero, measure!(zero_pyramid, tel, zero_src))
        @test all(isfinite, slopes(zero_pyramid))
        @test all(iszero, measure!(zero_bioedge, tel, zero_src))
        @test all(isfinite, slopes(zero_bioedge))
    end
    undetectable = PyramidWFS(tel; pupil_samples=4, mode=Diffractive(),
        modulation=0.0, normalization=IncidenceFluxNormalization())
    @test all(iszero, measure!(undetectable, tel, src,
        Detector(noise=NoiseNone(), integration_time=0.5, qe=0.0)))

    flat_tel = Telescope(resolution=32, diameter=8.0,
        central_obstruction=0.0)
    lgs = LGSSource(wavelength=589e-9, photon_irradiance=10.0,
        elongation_factor=1.8)
    spectral_lgs = with_spectrum(lgs,
        SpectralBundle([0.57e-6, 0.61e-6], [0.4, 0.6]))
    spectral_flat = PyramidWFS(flat_tel; pupil_samples=4,
        mode=Diffractive(), modulation=1.0,
        normalization=IncidenceFluxNormalization())
    @test measure!(spectral_flat, flat_tel, spectral_lgs) ≈
        zeros(length(slopes(spectral_flat))) atol=1e-12

    sampled_qe = SampledQuantumEfficiency(
        [0.50e-6, 0.57e-6, 0.61e-6, 0.70e-6],
        [0.1, 0.2, 0.9, 0.7])
    spectral_detector = Detector(noise=NoiseNone(), integration_time=0.4,
        qe=sampled_qe)
    spectral_detected = PyramidWFS(flat_tel; pupil_samples=4,
        mode=Diffractive(), modulation=1.0,
        normalization=IncidenceFluxNormalization())
    @test measure!(spectral_detected, flat_tel, spectral_lgs,
        spectral_detector) ≈ zeros(length(slopes(spectral_detected))) atol=1e-12
    spectral_scale = wfs_detector_incidence_scale(spectral_detector,
        spectral_lgs, eltype(output_frame(spectral_detector)))
    @test pyramid_signal_allocation_bytes(spectral_detected, flat_tel,
        output_frame(spectral_detector), spectral_lgs, spectral_scale) == 0

    aberrated_qe_tel = Telescope(resolution=32, diameter=8.0,
        central_obstruction=0.0)
    @inbounds for j in axes(aberrated_qe_tel.state.opd, 2),
            i in axes(aberrated_qe_tel.state.opd, 1)
        aberrated_qe_tel.state.opd[i, j] =
            2e-8 * sinpi(2 * i / 32) * cospi(2 * j / 32)
    end
    qe_reference_wfs = PyramidWFS(aberrated_qe_tel; pupil_samples=4,
        mode=Diffractive(), modulation=1.0,
        normalization=IncidenceFluxNormalization())
    qe_reference = copy(measure!(qe_reference_wfs, aberrated_qe_tel,
        spectral_lgs, Detector(noise=NoiseNone(), integration_time=1.0,
            qe=sampled_qe)))
    scaled_qe = SampledQuantumEfficiency(
        [0.50e-6, 0.57e-6, 0.61e-6, 0.70e-6],
        0.5 .* [0.1, 0.2, 0.9, 0.7])
    qe_scaled_wfs = PyramidWFS(aberrated_qe_tel; pupil_samples=4,
        mode=Diffractive(), modulation=1.0,
        normalization=IncidenceFluxNormalization())
    qe_scaled = copy(measure!(qe_scaled_wfs, aberrated_qe_tel,
        spectral_lgs, Detector(noise=NoiseNone(), integration_time=0.5,
            qe=scaled_qe)))
    @test norm(qe_reference) > 1e-6
    @test qe_scaled ≈ qe_reference atol=1e-12 rtol=1e-12
end


@testset "Pyramid and BioEdge detector-response references" begin
    tel = Telescope(resolution=32, diameter=8.0, central_obstruction=0.0)
    src = Source(band=:custom, wavelength=0.75e-6,
        photon_irradiance=10.0)
    asymmetric_response = [0.0 0.0 0.0;
                           0.0 0.2 0.8;
                           0.0 0.0 0.0]
    responses = (
        GaussianPixelResponse(response_width_px=0.75),
        SampledFrameResponse(asymmetric_response),
    )

    for response in responses, family in (:pyramid, :bioedge)
        wfs = family === :pyramid ?
            PyramidWFS(tel; pupil_samples=4, mode=Diffractive(),
                modulation=1.0) :
            BioEdgeWFS(tel; pupil_samples=4, mode=Diffractive(),
                modulation=1.0)
        det = Detector(noise=NoiseNone(), integration_time=0.4, qe=0.3,
            response_model=response)
        flat = copy(measure!(wfs, tel, src, det))
        @test all(iszero, flat)
        @test wfs.estimator.state.calibration_signature ==
            detector_calibration_signature(det,
                telescope_aperture_calibration_signature(tel,
                    calibration_signature(src)))
        @test detector_calibration_signature_allocation_bytes(det,
            telescope_aperture_calibration_signature(tel,
                calibration_signature(src))) == 0
        @test measure!(wfs, tel, src, det) == flat
    end

    make_wfs = family -> family === :pyramid ?
        PyramidWFS(tel; pupil_samples=4, mode=Diffractive(),
            modulation=1.0, light_ratio=0.0) :
        BioEdgeWFS(tel; pupil_samples=4, mode=Diffractive(),
            modulation=1.0, light_ratio=0.0)

    for family in (:pyramid, :bioedge)
        probe_wfs = make_wfs(family)
        probe_detector = Detector(noise=NoiseNone(), sensor=CMOSSensor())
        measure!(probe_wfs, tel, src, probe_detector;
            rng=MersenneTwister(701))
        detector_size = size(output_frame(probe_detector))

        gain_map = ones(detector_size)
        gain_map[1:2:end] .= 0.65
        bad_mask = falses(detector_size)
        bad_mask[2, 2] = true
        prnu = PixelResponseNonuniformity(gain_map)
        bad_pixels = BadPixelMask(bad_mask; throughput=0.0)
        detector = Detector(noise=NoiseNone(), sensor=CMOSSensor(),
            defect_model=CompositeDetectorDefectModel(prnu, bad_pixels))
        wfs = make_wfs(family)

        flat = copy(measure!(wfs, tel, src, detector;
            rng=MersenneTwister(702)))
        @test all(iszero, flat)
        @test all(iszero, measure!(wfs, tel, src, detector;
            rng=MersenneTwister(703)))
        signature = wfs.estimator.state.calibration_signature
        @test detector_calibration_signature_allocation_bytes(detector,
            telescope_aperture_calibration_signature(tel,
                calibration_signature(src))) == 0

        # Public models and Detector each own their parameter storage. Changes
        # to the model used to build a run do not mutate or invalidate that
        # run's frozen detector configuration.
        prnu.gain_map[1, 1] = 0.2
        bad_pixels.mask[end, end] = true
        @test detector_calibration_signature(detector,
            telescope_aperture_calibration_signature(tel,
                calibration_signature(src))) == signature

        replacement_gain = copy(gain_map)
        replacement_gain[1, 1] = 0.4
        replacement_mask = copy(bad_mask)
        replacement_mask[end, end] = true
        replacement = Detector(noise=NoiseNone(), sensor=CMOSSensor(),
            defect_model=CompositeDetectorDefectModel(
                PixelResponseNonuniformity(replacement_gain),
                BadPixelMask(replacement_mask; throughput=0.0)))
        replacement_signature = detector_calibration_signature(replacement,
            telescope_aperture_calibration_signature(tel,
                calibration_signature(src)))
        @test replacement_signature != signature
        @test all(iszero, measure!(wfs, tel, src, replacement;
            rng=MersenneTwister(704)))
        @test wfs.estimator.state.calibration_signature == replacement_signature
    end

    hgcdte_probe_wfs = make_wfs(:pyramid)
    hgcdte_probe_detector = Detector(noise=NoiseNone(),
        sensor=HgCdTeAvalancheArraySensor(
            sampling_mode=CorrelatedDoubleSampling()))
    measure!(hgcdte_probe_wfs, tel, src, hgcdte_probe_detector;
        rng=MersenneTwister(705))
    hgcdte_size = size(output_frame(hgcdte_probe_detector))
    hgcdte_gain_map = ones(hgcdte_size)
    hgcdte_gain_map[1:2:end] .= 0.7
    hgcdte_bad_mask = falses(hgcdte_size)
    hgcdte_bad_mask[2, 2] = true
    hgcdte_correction = ReferencePixelCommonModeCorrection(1, 1)
    hgcdte_detector = Detector(noise=NoiseNone(),
        sensor=HgCdTeAvalancheArraySensor(
            sampling_mode=CorrelatedDoubleSampling()),
        correction_model=hgcdte_correction,
        defect_model=CompositeDetectorDefectModel(
            PixelResponseNonuniformity(hgcdte_gain_map),
            BadPixelMask(hgcdte_bad_mask; throughput=0.0)))
    hgcdte_wfs = make_wfs(:pyramid)
    @test all(iszero, measure!(hgcdte_wfs, tel, src, hgcdte_detector;
        rng=MersenneTwister(706)))
    hgcdte_seed = telescope_aperture_calibration_signature(tel,
        calibration_signature(src))
    @test detector_calibration_signature_allocation_bytes(hgcdte_detector,
        hgcdte_seed) == 0

    correction_one = Detector(noise=NoiseNone(),
        sensor=HgCdTeAvalancheArraySensor(
            sampling_mode=CorrelatedDoubleSampling()),
        correction_model=ReferencePixelCommonModeCorrection(1, 1))
    correction_two = Detector(noise=NoiseNone(),
        sensor=HgCdTeAvalancheArraySensor(
            sampling_mode=CorrelatedDoubleSampling()),
        correction_model=ReferencePixelCommonModeCorrection(2, 1))
    @test detector_calibration_signature(correction_one, hgcdte_seed) !=
        detector_calibration_signature(correction_two, hgcdte_seed)
    composite_correction = Detector(noise=NoiseNone(),
        sensor=HgCdTeAvalancheArraySensor(
            sampling_mode=CorrelatedDoubleSampling()),
        correction_model=CompositeFrameReadoutCorrection((
            ReferenceRowCommonModeCorrection(1),
            ReferenceColumnCommonModeCorrection(1),
        )))
    @test detector_calibration_signature_allocation_bytes(
        composite_correction, hgcdte_seed) == 0

    ramp_correction = ReferencePixelCommonModeCorrection(1, 1)
    ramp_detector = Detector(noise=NoiseNone(), integration_time=1.0,
        gain=3.0,
        sensor=HgCdTeAvalancheArraySensor(avalanche_gain=2.0,
            read_time=0.1, sampling_mode=UpTheRampSampling(4)),
        correction_model=ramp_correction)
    @test detector_calibration_signature_allocation_bytes(
        ramp_detector, hgcdte_seed) == 0
    ramp_flat = measure!(make_wfs(:pyramid), tel, src,
        ramp_detector; rng=MersenneTwister(707))
    @test maximum(abs, ramp_flat) <= 16eps(eltype(ramp_flat))

    gain_signature_detector = Detector(noise=NoiseNone(), gain=2.0,
        sensor=HgCdTeAvalancheArraySensor(
            sampling_mode=CorrelatedDoubleSampling()),
        correction_model=hgcdte_correction)
    avalanche_signature_detector = Detector(noise=NoiseNone(),
        sensor=HgCdTeAvalancheArraySensor(avalanche_gain=2.0,
            sampling_mode=CorrelatedDoubleSampling()),
        correction_model=hgcdte_correction)
    ramp_read_time_detector = Detector(noise=NoiseNone(),
        integration_time=1.0, gain=3.0,
        sensor=HgCdTeAvalancheArraySensor(avalanche_gain=2.0,
            read_time=0.2, sampling_mode=UpTheRampSampling(4)),
        correction_model=ramp_correction)
    @test detector_calibration_signature(gain_signature_detector,
        hgcdte_seed) != detector_calibration_signature(correction_one,
        hgcdte_seed)
    @test detector_calibration_signature(avalanche_signature_detector,
        hgcdte_seed) != detector_calibration_signature(correction_one,
        hgcdte_seed)
    @test detector_calibration_signature(ramp_read_time_detector,
        hgcdte_seed) != detector_calibration_signature(ramp_detector,
        hgcdte_seed)

    windowed = Detector(noise=NoiseNone(), integration_time=1.0, qe=1.0,
        readout_window=FrameWindow(1:2, 1:2))
    @test_throws InvalidConfiguration measure!(
        PyramidWFS(tel; pupil_samples=4, mode=Diffractive()),
        tel, src, windowed)

    unsupported_calibration_detectors = (
        Detector(noise=NoiseNone(),
            charge_coupling_model=InterpixelCapacitance(
                [0.0 0.1 0.0; 0.1 0.6 0.1; 0.0 0.1 0.0])),
        Detector(noise=NoiseNone(), sensor=InGaAsSensor(),
            response_model=NullFrameResponse(),
            nonlinearity_model=SaturatingFrameNonlinearity(0.1)),
        Detector(noise=NoiseNone(), sensor=CMOSSensor(),
            defect_model=DarkSignalNonuniformity(ones(8, 8))),
        Detector(noise=NoiseNone(),
            sensor=HgCdTeAvalancheArraySensor(
                read_time=0.4, sampling_mode=UpTheRampSampling(4))),
    )
    for detector in unsupported_calibration_detectors
        @test_throws InvalidConfiguration measure!(
            PyramidWFS(tel; pupil_samples=4, mode=Diffractive()),
            tel, src, detector)
        @test_throws InvalidConfiguration measure!(
            BioEdgeWFS(tel; pupil_samples=4, mode=Diffractive()),
            tel, src, detector)
    end

    fill!(tel.state.opd, 1e-8)
    opd_before_failure = copy(tel.state.opd)
    for family in (:pyramid, :bioedge)
        wfs = family === :pyramid ?
            PyramidWFS(tel; pupil_samples=4, mode=Diffractive()) :
            BioEdgeWFS(tel; pupil_samples=4, mode=Diffractive())
        invalid_sampling = Detector(noise=NoiseNone(), psf_sampling=3)
        @test_throws DimensionMismatchError measure!(
            wfs, tel, src, invalid_sampling)
        @test tel.state.opd == opd_before_failure
    end
end

@testset "BioEdge source-composition support boundary" begin
    tel = Telescope(resolution=16, diameter=8.0,
        central_obstruction=0.0)
    src = Source(band=:custom, wavelength=0.75e-6,
        photon_irradiance=1.0)
    spectral = with_spectrum(src,
        SpectralBundle([0.70e-6, 0.80e-6], [0.5, 0.5]))
    extended = with_extended_source(src,
        PointCloudSourceModel([(0.0, 0.0)], [1.0]))

    for expanded in (spectral, extended)
        @test_throws UnsupportedAlgorithm measure!(
            BioEdgeWFS(tel; pupil_samples=2, mode=Diffractive()),
            tel, expanded)
        @test_throws UnsupportedAlgorithm measure!(
            BioEdgeWFS(tel; pupil_samples=2, mode=Diffractive()),
            tel, expanded,
            Detector(noise=NoiseNone(), integration_time=1.0, qe=1.0))
    end
end

@testset "Pyramid and BioEdge pupil-reflectivity throughput" begin
    transmission = 0.25
    src = Source(band=:custom, wavelength=0.75e-6,
        photon_irradiance=1.0)
    full_tel = Telescope(resolution=16, diameter=8.0,
        central_obstruction=0.0)
    attenuated_tel = Telescope(resolution=16, diameter=8.0,
        central_obstruction=0.0, pupil_reflectivity=transmission)

    for style in (ScalarCPUStyle(), KA_CPU_STYLE)
        full_wfs = PyramidWFS(full_tel; pupil_samples=2,
            mode=Diffractive(), modulation=0.0)
        attenuated_wfs = PyramidWFS(attenuated_tel; pupil_samples=2,
            mode=Diffractive(), modulation=0.0)
        pyramid_intensity_core!(style, full_wfs.front_end.propagation.intensity,
            full_wfs, full_tel, src, pyramid_operating_modulation(full_wfs))
        pyramid_intensity_core!(style, attenuated_wfs.front_end.propagation.intensity,
            attenuated_wfs, attenuated_tel, src,
            pyramid_operating_modulation(attenuated_wfs))
        full_rate = sum(full_wfs.front_end.propagation.intensity)
        @test full_rate > 0
        @test sum(attenuated_wfs.front_end.propagation.intensity) ≈ transmission * full_rate rtol=1e-12
    end

    full_pyramid = PyramidWFS(full_tel; pupil_samples=2,
        mode=Diffractive(), modulation=0.0)
    attenuated_pyramid = PyramidWFS(attenuated_tel; pupil_samples=2,
        mode=Diffractive(), modulation=0.0)
    full_modulation_frame = similar(full_pyramid.front_end.propagation.intensity)
    attenuated_modulation_frame = similar(attenuated_pyramid.front_end.propagation.intensity)
    pyramid_modulation_frame!(full_modulation_frame, full_pyramid, full_tel, src)
    pyramid_modulation_frame!(attenuated_modulation_frame,
        attenuated_pyramid, attenuated_tel, src)
    @test sum(attenuated_modulation_frame) ≈
        transmission * sum(full_modulation_frame) rtol=1e-12

    full_bioedge = BioEdgeWFS(full_tel; pupil_samples=2,
        mode=Diffractive(), modulation=0.0)
    attenuated_bioedge = BioEdgeWFS(attenuated_tel; pupil_samples=2,
        mode=Diffractive(), modulation=0.0)
    bioedge_intensity_core!(full_bioedge.front_end.propagation.intensity,
        full_bioedge, full_tel, src)
    bioedge_intensity_core!(attenuated_bioedge.front_end.propagation.intensity,
        attenuated_bioedge, attenuated_tel, src)
    full_bioedge_rate = sum(full_bioedge.front_end.propagation.intensity)
    @test full_bioedge_rate > 0
    @test sum(attenuated_bioedge.front_end.propagation.intensity) ≈
        transmission * full_bioedge_rate rtol=1e-12
end

@testset "Shack-Hartmann subapertures" begin
    tel = Telescope(resolution=24, diameter=8.0, central_obstruction=0.1)
    src = Source(band=:I, magnitude=0.0)
    sh = ShackHartmannWFS(tel; n_lenslets=6, mode=Diffractive(), pixel_scale_arcsec=0.06, n_pix_subap=8, threshold_cog=0.02)

    layout = subaperture_layout(sh)
    calibration = subaperture_calibration(sh)
    @test layout isa SubapertureLayout
    @test calibration isa SubapertureCalibration
    @test layout.n_subap == 6
    @test layout.subap_pixels == 4
    @test layout.pitch_m ≈ tel.params.diameter / 6
    @test !calibration.calibrated
    @test slope_extraction_model(sh) isa CenterOfGravityExtraction
    @test slope_extraction_model(sh).threshold ≈ 0.02
    @test n_valid_subapertures(layout) == count(layout.valid_mask_host)
    @test valid_subaperture_indices(layout) == findall(layout.valid_mask_host)

    prepare_runtime_wfs!(sh, tel, src)
    @test calibration.calibrated
    @test calibration.centroid_response == sh.calibration.centroid_response
    @test calibration.wavelength == sh.calibration.wavelength
    @test calibration.signature == sh.calibration.signature
    @test calibration.reference_signal_2d ===
        sh.calibration.reference_signal_2d
    @test calibration.reference_signal_host ===
        sh.calibration.reference_signal_host
    @test length(valid_subaperture_indices(layout)) == n_valid_subapertures(layout)

    slopes = measure!(sh, tel, src)
    @test all(isfinite, slopes)
    meta = AdaptiveOpticsSim.wfs_output_metadata(sh)
    @test meta.n_valid_subap == n_valid_subapertures(layout)
    @test meta.subap_pixels == layout.subap_pixels
    @test meta.calibrated

    dm = DeformableMirror(tel; n_act=5)
    imat = interaction_matrix(dm, sh, tel, src; amplitude=1e-8)
    @test size(imat.matrix, 1) == length(AdaptiveOpticsSim.slopes(sh))
    @test size(imat.matrix, 2) == length(dm.state.coefs)
end
