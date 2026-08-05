@testset "Detector response, MTF, and batched capture" begin
    function impulse_transfer_magnitude(frame, spatial_frequency_x,
        spatial_frequency_y)
        center_i = fld(size(frame, 1), 2) + 1
        center_j = fld(size(frame, 2), 2) + 1
        response = zero(Complex{eltype(frame)})
        @inbounds for j in axes(frame, 2), i in axes(frame, 1)
            phase = -2pi * (spatial_frequency_y * (i - center_i) +
                            spatial_frequency_x * (j - center_j))
            response += frame[i, j] * cis(phase)
        end
        return abs(response) / sum(frame)
    end

    det_mtf = Detector(exposure_duration=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        response_model=GaussianPixelResponse(response_width_px=0.75))
    impulse = zeros(9, 9)
    impulse[5, 5] = 1.0
    frame_mtf = capture!(det_mtf, impulse; rng=MersenneTwister(3))
    @test sum(frame_mtf) ≈ 1.0 atol=1e-6
    @test frame_mtf[5, 5] < 1.0
    @test frame_mtf[5, 4] > 0
    @test supports_detector_mtf(det_mtf)
    @test detector_mtf(det_mtf, 0.0, 0.0) ≈ 1.0
    @test detector_mtf(det_mtf, 0.5, 0.0) < 1.0
    @test detector_mtf(det_mtf, 0.5, 0.0) ≈
        impulse_transfer_magnitude(frame_mtf, 0.5, 0.0)
    mtf_meta = detector_export_metadata(det_mtf)
    @test mtf_meta.frame_response == :gaussian
    @test mtf_meta.response_width_px == 0.75
    @test mtf_meta.response_application_domain == :image
    @test mtf_meta.response_is_separable
    @test !mtf_meta.response_is_shift_invariant
    @test mtf_meta.response_support_rows == mtf_meta.response_support_cols
    @test mtf_meta.pitch_x_px === nothing
    @test mtf_meta.aperture_shape === nothing
    @test_throws InvalidConfiguration GaussianPixelResponse(response_width_px=0.0)

    sampled_kernel = [0.0 0.1 0.0; 0.1 0.6 0.1; 0.0 0.1 0.0]
    sampled_det = Detector(exposure_duration=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        response_model=SampledFrameResponse(sampled_kernel))
    sampled_frame = capture!(sampled_det, impulse; rng=MersenneTwister(6))
    @test sum(sampled_frame) ≈ 1.0 atol=1e-6
    @test sampled_frame[5, 5] ≈ 0.6 atol=1e-6
    @test sampled_frame[5, 4] ≈ 0.1 atol=1e-6
    sampled_meta = detector_export_metadata(sampled_det)
    @test sampled_meta.frame_response == :sampled
    @test sampled_meta.response_application_domain == :image
    @test !sampled_meta.response_is_separable
    @test sampled_meta.response_support_rows == 3
    @test sampled_meta.response_support_cols == 3
    @test sampled_meta.aperture_shape == :sampled
    @test supports_detector_mtf(sampled_det)
    @test detector_mtf(sampled_det, 0.0, 0.0) ≈ 1.0
    @test_throws InvalidConfiguration SampledFrameResponse(zeros(3, 3))
    @test_throws InvalidConfiguration SampledFrameResponse(ones(2, 3))
    negative_kernel = copy(sampled_kernel)
    negative_kernel[1, 1] = -0.1
    @test_throws InvalidConfiguration SampledFrameResponse(negative_kernel)
    nan_kernel = copy(sampled_kernel)
    nan_kernel[1, 1] = NaN
    @test_throws InvalidConfiguration SampledFrameResponse(nan_kernel)
    infinite_kernel = copy(sampled_kernel)
    infinite_kernel[1, 1] = Inf
    @test_throws InvalidConfiguration SampledFrameResponse(infinite_kernel)
    @test_throws InvalidConfiguration SampledFrameResponse(
        fill(0.2, 3, 3); normalize=false)

    asymmetric_kernel = [0.0 0.0 0.0; 0.1 0.2 0.7; 0.0 0.0 0.0]
    asymmetric_det = Detector(exposure_duration=1.0, noise=NoiseNone(),
        qe=1.0, binning=1,
        response_model=SampledFrameResponse(asymmetric_kernel))
    center_impulse = zeros(9, 9)
    center_impulse[5, 5] = 1.0
    center_frame = copy(capture!(asymmetric_det, center_impulse;
        rng=MersenneTwister(61)))
    left_impulse = zeros(9, 9)
    left_impulse[5, 1] = 1.0
    left_frame = copy(capture!(asymmetric_det, left_impulse;
        rng=MersenneTwister(62)))
    right_impulse = zeros(9, 9)
    right_impulse[5, end] = 1.0
    right_frame = copy(capture!(asymmetric_det, right_impulse;
        rng=MersenneTwister(63)))
    @test sum(center_frame) ≈ 1.0
    @test center_frame[5, 4:6] ≈ asymmetric_kernel[2, :]
    @test sum(left_frame) ≈ 0.9
    @test sum(right_frame) ≈ 0.3
    @test sum(left_frame) <= sum(left_impulse)
    @test sum(right_frame) <= sum(right_impulse)
    @test minimum(left_frame) >= 0
    @test minimum(right_frame) >= 0
    @test detector_mtf(asymmetric_det, 0.37, 0.0) ≈
        impulse_transfer_magnitude(center_frame, 0.37, 0.0)
    @test !detector_export_metadata(asymmetric_det).response_is_shift_invariant

    cube_mtf = Array{Float64}(undef, 2, size(impulse, 1), size(impulse, 2))
    cube_mtf[1, :, :] .= impulse
    cube_mtf[2, :, :] .= impulse
    scratch_mtf = similar(cube_mtf)
    stack_mtf = AdaptiveOpticsSim.Detectors.capture_stack!(det_mtf, cube_mtf, scratch_mtf; rng=MersenneTwister(10))
    @test size(stack_mtf) == size(cube_mtf)
    @test all(isfinite, stack_mtf)
    cube_sampled = Array{Float64}(undef, 2, size(impulse, 1), size(impulse, 2))
    cube_sampled[1, :, :] .= impulse
    cube_sampled[2, :, :] .= impulse
    scratch_sampled = similar(cube_sampled)
    stack_sampled = AdaptiveOpticsSim.Detectors.capture_stack!(sampled_det, cube_sampled, scratch_sampled; rng=MersenneTwister(10))
    @test size(stack_sampled) == size(cube_sampled)
    @test all(isfinite, stack_sampled)
    @test stack_sampled[1, :, :] ≈ sampled_frame atol=1e-6
    @test stack_sampled[2, :, :] ≈ sampled_frame atol=1e-6
    det_stack_adc = Detector(exposure_duration=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        bits=8, full_well=10.0, output_type=UInt16)
    cube_stack_adc = fill(10.0, 2, 4, 4)
    stack_adc = AdaptiveOpticsSim.Detectors.capture_stack!(det_stack_adc, cube_stack_adc, similar(cube_stack_adc);
        rng=MersenneTwister(10))
    @test stack_adc === cube_stack_adc
    @test all(stack_adc .== 255.0)

    rect_det = Detector(exposure_duration=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        response_model=RectangularPixelAperture(pitch_x_px=2.0, pitch_y_px=2.0,
            fill_factor_x=0.6, fill_factor_y=0.8))
    rect_frame = capture!(rect_det, impulse; rng=MersenneTwister(4))
    @test sum(rect_frame) ≈ 1.0 atol=1e-6
    @test rect_frame[5, 5] < 1.0
    rect_meta = detector_export_metadata(rect_det)
    @test rect_meta.frame_response == :rectangular_aperture
    @test rect_meta.pitch_x_px == 2.0
    @test rect_meta.pitch_y_px == 2.0
    @test rect_meta.fill_factor_x == 0.6
    @test rect_meta.fill_factor_y == 0.8
    @test rect_meta.aperture_shape == :rectangular
    @test rect_meta.response_application_domain == :image
    @test supports_detector_mtf(rect_det)
    @test detector_mtf(rect_det, 0.5, 0.0) ≈
        impulse_transfer_magnitude(rect_frame, 0.5, 0.0)
    unit_rectangular = RectangularPixelAperture()
    half_fill_rectangular = RectangularPixelAperture(fill_factor_x=0.5,
        fill_factor_y=0.5)
    @test detector_mtf(unit_rectangular, 0.5, 0.0) ≈ 1.0
    @test unit_rectangular.kernel_x == half_fill_rectangular.kernel_x
    @test unit_rectangular.kernel_y == half_fill_rectangular.kernel_y
    @test detector_mtf(unit_rectangular, 0.5, 0.0) ==
        detector_mtf(half_fill_rectangular, 0.5, 0.0)
    @test !AdaptiveOpticsSim.Detectors.supports_subpixel_geometry(unit_rectangular)
    @test !AdaptiveOpticsSim.Detectors.supports_subpixel_geometry(
        half_fill_rectangular)

    @test_throws InvalidConfiguration RectangularPixelAperture(fill_factor_x=0.0)
    @test_throws InvalidConfiguration RectangularPixelAperture(fill_factor_y=1.5)
    @test_throws InvalidConfiguration RectangularPixelAperture(pitch_x_px=Inf)
    @test_throws InvalidConfiguration RectangularPixelAperture(fill_factor_y=NaN)
    @test_throws InvalidConfiguration GaussianPixelResponse(response_width_px=Inf)
    @test_throws InvalidConfiguration GaussianPixelResponse(truncate_at=Inf)

    ipc_kernel = [0.0 0.01 0.0; 0.01 0.96 0.01; 0.0 0.01 0.0]
    ipc_det = Detector(exposure_duration=1.0, noise=NoisePhoton(), qe=1.0,
        response_model=NullFrameResponse(),
        charge_coupling_model=InterpixelCapacitance(ipc_kernel))
    ipc_input = zeros(9, 9)
    ipc_input[5, 5] = 100.0
    ipc_frame = capture!(ipc_det, ipc_input; rng=MersenneTwister(45))
    @test ipc_frame[5, 4] > 0
    @test !isinteger(ipc_frame[5, 4])
    @test sum(ipc_frame) ≈ round(sum(ipc_frame)) atol=1e-10
    ipc_meta = detector_export_metadata(ipc_det)
    @test ipc_meta.charge_coupling == :interpixel_capacitance
    @test ipc_meta.charge_coupling_support_rows == 3
    @test ipc_meta.charge_coupling_support_cols == 3

    det_window_stack = Detector(exposure_duration=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        readout_window=FrameWindow(2:8, 2:8))
    cube_window = Array{Float64}(undef, 2, size(impulse, 1), size(impulse, 2))
    cube_window[1, :, :] .= impulse
    cube_window[2, :, :] .= impulse
    scratch_window = similar(cube_window)
    @test_throws InvalidConfiguration AdaptiveOpticsSim.Detectors.capture_stack!(det_window_stack, cube_window, scratch_window; rng=MersenneTwister(10))
    input_window_stack = copy(cube_window)
    output_window_stack = Array{Float64}(undef, 2, 7, 7)
    generalized_window = AdaptiveOpticsSim.Detectors.capture_stack!(det_window_stack, output_window_stack, input_window_stack; rng=MersenneTwister(10))
    @test size(generalized_window) == (2, 7, 7)
    @test generalized_window[1, :, :] ≈ capture!(det_window_stack, impulse; rng=MersenneTwister(10))

    corrected_stack_models = (
        ReferencePixelCommonModeCorrection(1, 1),
        ReferenceRowCommonModeCorrection(1),
        ReferenceColumnCommonModeCorrection(1),
        ReferenceOutputCommonModeCorrection(2; edge_rows=1, edge_cols=1),
        CompositeFrameReadoutCorrection((
            ReferenceRowCommonModeCorrection(1),
            ReferenceColumnCommonModeCorrection(1),
        )),
    )
    for correction_model in corrected_stack_models
        corrected_stack_det = Detector(exposure_duration=1.0, noise=NoiseNone(), qe=1.0, binning=1,
            sensor=HgCdTeSensor(sampling_mode=SingleRead()),
            correction_model=correction_model)
        corrected_stack_in = Array{Float64}(undef, 2, 5, 5)
        corrected_stack_in[1, :, :] .= reshape(collect(1.0:25.0), 5, 5)
        corrected_stack_in[2, :, :] .= reshape(collect(26.0:50.0), 5, 5)
        corrected_stack_ref = copy(corrected_stack_in)
        corrected_stack = AdaptiveOpticsSim.Detectors.capture_stack!(corrected_stack_det, corrected_stack_in,
            similar(corrected_stack_in); rng=MersenneTwister(10))
        @test size(corrected_stack) == size(corrected_stack_in)
        corrected_frame_1 = capture!(Detector(exposure_duration=1.0, noise=NoiseNone(), qe=1.0, binning=1,
                sensor=HgCdTeSensor(sampling_mode=SingleRead()),
                correction_model=correction_model),
            @view(corrected_stack_ref[1, :, :]); rng=MersenneTwister(10))
        corrected_frame_2 = capture!(Detector(exposure_duration=1.0, noise=NoiseNone(), qe=1.0, binning=1,
                sensor=HgCdTeSensor(sampling_mode=SingleRead()),
                correction_model=correction_model),
            @view(corrected_stack_ref[2, :, :]); rng=MersenneTwister(10))
        @test corrected_stack[1, :, :] ≈ corrected_frame_1 atol=1e-12 rtol=1e-12
        @test corrected_stack[2, :, :] ≈ corrected_frame_2 atol=1e-12 rtol=1e-12
    end

    det_cmos_batched = Detector(exposure_duration=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        sensor=CMOSSensor(column_readout_sigma=1.0))
    @test_throws InvalidConfiguration AdaptiveOpticsSim.Detectors.capture_stack!(det_cmos_batched, cube_mtf, scratch_mtf; rng=MersenneTwister(10))

    det_generalized = Detector(exposure_duration=1.0, noise=NoiseNone(), qe=1.0, psf_sampling=2, binning=2,
        bits=8, full_well=10.0, output_type=UInt8)
    input_generalized = zeros(Float64, 2, 8, 8)
    input_generalized[1, 4, 4] = 10.0
    input_generalized[2, 5, 5] = 10.0
    output_generalized = Array{UInt8}(undef, 2, 2, 2)
    generalized_stack = AdaptiveOpticsSim.Detectors.capture_stack!(det_generalized, output_generalized, input_generalized; rng=MersenneTwister(10))
    @test size(generalized_stack) == (2, 2, 2)
    @test readout_ready(det_generalized)
    @test iszero(det_generalized.state.integrated_time)
    @test generalized_stack[1, :, :] == capture!(det_generalized, @view(input_generalized[1, :, :]); rng=MersenneTwister(10))
    @test generalized_stack[2, :, :] == capture!(det_generalized, @view(input_generalized[2, :, :]); rng=MersenneTwister(10))

end
