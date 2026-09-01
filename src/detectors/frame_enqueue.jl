@kernel function _prepare_enqueued_detector_frame_kernel!(
    frame,
    photon_rate,
    scale,
    count::Int,
)
    index = @index(Global, Linear)
    if index <= count
        @inbounds frame[index] = photon_rate[index] * scale
    end
end

@kernel function _fill_enqueued_detector_scratch_kernel!(
    scratch,
    value,
    count::Int,
)
    index = @index(Global, Linear)
    if index <= count
        @inbounds scratch[index] = value
    end
end

@kernel function _add_enqueued_detector_scratch_kernel!(
    frame,
    scratch,
    scale,
    count::Int,
)
    index = @index(Global, Linear)
    if index <= count
        @inbounds frame[index] += scale * scratch[index]
    end
end

@kernel function _scale_enqueued_detector_frame_kernel!(
    frame,
    scale,
    count::Int,
)
    index = @index(Global, Linear)
    if index <= count
        @inbounds frame[index] *= scale
    end
end

@kernel function _clamp_enqueued_detector_frame_kernel!(
    frame,
    lower,
    upper,
    count::Int,
)
    index = @index(Global, Linear)
    if index <= count
        @inbounds frame[index] = clamp(frame[index], lower, upper)
    end
end

@kernel function _multiply_enqueued_emccd_frame_kernel!(
    frame,
    scratch,
    gain,
    factor,
    count::Int,
)
    index = @index(Global, Linear)
    if index <= count
        value = max(@inbounds(frame[index]), zero(eltype(frame)))
        noise = @inbounds scratch[index]
        @inbounds frame[index] = gain *
            max(value + factor * sqrt(value) * noise, zero(value))
    end
end

@kernel function _quantize_enqueued_detector_frame_kernel!(
    frame,
    scale,
    upper,
    count::Int,
)
    index = @index(Global, Linear)
    if index <= count
        value = @inbounds(frame[index]) * scale
        @inbounds frame[index] = clamp(value, zero(value), upper)
    end
end

@kernel function _write_enqueued_detector_integer_output_kernel!(
    output,
    frame,
    count::Int,
)
    index = @index(Global, Linear)
    if index <= count
        # Quantization has already clamped the floating frame to the exact ADC
        # range represented by `output`. Avoid checked float-to-integer
        # conversion in device code: the range proof belongs to graph
        # admission and the immediately preceding quantization kernel.
        value = round(@inbounds frame[index])
        @inbounds output[index] = unsafe_trunc(eltype(output), value)
    end
end

@inline _supports_enqueued_detector_sensor(::AbstractFrameSensor) = false

@inline _supports_enqueued_detector_sensor(
    ::CCDSensor{T,SingleRead},
) where {T<:AbstractFloat} = true

@inline _supports_enqueued_detector_sensor(
    ::CMOSSensor{
        T,
        NullCMOSReadNoise,
        NullCMOSOutputModel,
        GlobalShutter,
    },
) where {T<:AbstractFloat} = true

@inline _supports_enqueued_detector_sensor(
    ::EMCCDSensor{
        T,
        M,
        LinearEMMode,
        EMOutput,
        SequentialAcquisition,
    },
) where {
    T<:AbstractFloat,
    M<:ClippedGaussianMultiplicationApproximation,
} = true

@inline _supports_enqueued_detector_models(args...) = false

@inline function _supports_enqueued_detector_models(
    ::NullFrameResponse,
    ::NullChargeCoupling,
    ::NullDetectorDefectModel,
    ::GlobalShutter,
    ::NullFrameReadoutCorrection,
    ::NullFrameNonlinearity,
    ::NullDetectorThermalModel,
    ::NoBackground,
    ::NoBackground,
    ::NoFrameReadoutProducts,
)
    return true
end

@inline _supports_enqueued_detector_output(
    ::Nothing,
    bits::Union{Nothing,Int},
    ::Type{<:AbstractFloat},
) = true
@inline _supports_enqueued_detector_output(
    ::AbstractMatrix{T},
    bits::Int,
    ::Type{F},
) where {T<:Integer,F<:AbstractFloat} =
    bits <= 8 * sizeof(T) && bits <= precision(F)
@inline _supports_enqueued_detector_output(output, bits, frame_type) = false

@inline function _supports_enqueued_detector_acquisition(
    prepared::PreparedDetectorAcquisition,
    rng::AbstractRNG,
    output,
)
    return false
end

@inline function _supports_enqueued_detector_acquisition(
    prepared::PreparedDetectorAcquisition,
    rng::_PreparedCounterRNG,
    output,
)
    detector = prepared.detector
    return _supports_enqueued_detector_acquisition(
        execution_style(detector.products.frame),
        prepared,
        output,
    )
end

@inline function _supports_enqueued_detector_acquisition(
    ::ExecutionStyle,
    prepared::PreparedDetectorAcquisition,
    output,
)
    return false
end

@inline function _supports_enqueued_detector_acquisition(
    ::AcceleratorStyle,
    prepared::PreparedDetectorAcquisition,
    output,
)
    detector = prepared.detector
    params = detector.params
    frame = detector.products.frame
    input = prepared.input.values
    detector_output = detector.products.output_buffer
    return params.psf_sampling == 1 &&
        params.binning == 1 &&
        params.readout_window === nothing &&
        size(input) == size(frame) == size(output) &&
        eltype(input) === eltype(frame) === eltype(output) &&
        compute_device(input) == compute_device(frame) ==
            compute_device(output) &&
        _supports_enqueued_detector_sensor(params.sensor) &&
        _supports_enqueued_detector_models(
            params.response_model,
            params.charge_coupling_model,
            params.defect_model,
            params.timing_model,
            params.correction_model,
            params.nonlinearity_model,
            params.thermal_model,
            detector.background_flux,
            detector.background_map,
            detector.products.readout,
        ) &&
        _supports_enqueued_detector_output(
            detector_output,
            params.bits,
            eltype(frame),
        ) &&
        (detector_output === nothing ||
            compute_device(detector_output) == compute_device(frame))
end

@inline function _enqueue_detector_poisson_rate!(
    style::AcceleratorStyle,
    detector::Detector,
    rng::_PreparedCounterRNG,
    rate,
)
    frame = detector.products.frame
    value = eltype(frame)(rate)
    value <= zero(value) && return frame
    scratch = detector.workspace.noise_buffer
    launch_kernel_async!(
        style,
        _fill_enqueued_detector_scratch_kernel!,
        scratch,
        value,
        length(scratch);
        ndrange=length(scratch),
    )
    poisson_noise_async!(style, rng, scratch)
    launch_kernel_async!(
        style,
        _add_enqueued_detector_scratch_kernel!,
        frame,
        scratch,
        one(value),
        length(frame);
        ndrange=length(frame),
    )
    return frame
end

@inline _enqueue_detector_photon_noise!(
    style::AcceleratorStyle,
    detector::Detector{NoiseNone},
    rng::_PreparedCounterRNG,
) = detector.products.frame

@inline _enqueue_detector_photon_noise!(
    style::AcceleratorStyle,
    detector::Detector{<:NoiseReadout},
    rng::_PreparedCounterRNG,
) = detector.products.frame

@inline function _enqueue_detector_photon_noise!(
    style::AcceleratorStyle,
    detector::Detector{NoisePhoton},
    rng::_PreparedCounterRNG,
)
    return poisson_noise_async!(style, rng, detector.products.frame)
end

@inline function _enqueue_detector_photon_noise!(
    style::AcceleratorStyle,
    detector::Detector{<:NoisePhotonReadout},
    rng::_PreparedCounterRNG,
)
    return poisson_noise_async!(style, rng, detector.products.frame)
end

@inline function _enqueue_detector_sensor_statistics!(
    style::AcceleratorStyle,
    sensor::CCDSensor,
    detector::Detector,
    rng::_PreparedCounterRNG,
)
    return _enqueue_detector_poisson_rate!(
        style,
        detector,
        rng,
        effective_cic_per_frame(detector),
    )
end

@inline _enqueue_detector_sensor_statistics!(
    style::AcceleratorStyle,
    sensor::CMOSSensor,
    detector::Detector,
    rng::_PreparedCounterRNG,
) = detector.products.frame

@inline function _enqueue_detector_sensor_statistics!(
    style::AcceleratorStyle,
    sensor::EMCCDSensor,
    detector::Detector,
    rng::_PreparedCounterRNG,
)
    return _enqueue_detector_poisson_rate!(
        style,
        detector,
        rng,
        effective_cic_per_frame(detector),
    )
end

@inline function _enqueue_detector_clamp!(
    style::AcceleratorStyle,
    detector::Detector,
)
    frame = detector.products.frame
    upper = sensor_saturation_limit(detector)
    upper === nothing && return frame
    launch_kernel_async!(
        style,
        _clamp_enqueued_detector_frame_kernel!,
        frame,
        zero(eltype(frame)),
        eltype(frame)(upper),
        length(frame);
        ndrange=length(frame),
    )
    return frame
end

@inline function _enqueue_detector_scale!(
    style::AcceleratorStyle,
    detector::Detector,
    scale,
)
    frame = detector.products.frame
    value = eltype(frame)(scale)
    isone(value) && return frame
    launch_kernel_async!(
        style,
        _scale_enqueued_detector_frame_kernel!,
        frame,
        value,
        length(frame);
        ndrange=length(frame),
    )
    return frame
end

@inline _enqueue_detector_pre_readout_gain!(
    style::AcceleratorStyle,
    sensor::Union{CCDSensor,CMOSSensor},
    detector::Detector,
    rng::_PreparedCounterRNG,
) = detector.products.frame

@inline function _enqueue_detector_pre_readout_gain!(
    style::AcceleratorStyle,
    sensor::EMCCDSensor{
        T,
        M,
        LinearEMMode,
        EMOutput,
        SequentialAcquisition,
    },
    detector::Detector,
    rng::_PreparedCounterRNG,
) where {
    T<:AbstractFloat,
    M<:ClippedGaussianMultiplicationApproximation,
}
    frame = detector.products.frame
    gain = eltype(frame)(detector.params.gain)
    factor = eltype(frame)(em_conditional_noise_factor(
        sensor.multiplication_model,
        sensor,
    ))
    factor <= zero(factor) &&
        return _enqueue_detector_scale!(style, detector, gain)
    scratch = detector.workspace.noise_buffer
    randn_backend_async!(style, rng, scratch)
    launch_kernel_async!(
        style,
        _multiply_enqueued_emccd_frame_kernel!,
        frame,
        scratch,
        gain,
        factor,
        length(frame);
        ndrange=length(frame),
    )
    return frame
end

@inline function _enqueue_detector_gaussian_noise!(
    style::AcceleratorStyle,
    detector::Detector,
    rng::_PreparedCounterRNG,
    sigma,
)
    frame = detector.products.frame
    value = eltype(frame)(sigma)
    value <= zero(value) && return frame
    scratch = detector.workspace.noise_buffer
    randn_backend_async!(style, rng, scratch)
    launch_kernel_async!(
        style,
        _add_enqueued_detector_scratch_kernel!,
        frame,
        scratch,
        value,
        length(frame);
        ndrange=length(frame),
    )
    return frame
end

@inline _enqueue_detector_readout_noise!(
    style::AcceleratorStyle,
    detector::Detector{NoiseNone},
    rng::_PreparedCounterRNG,
) = detector.products.frame

@inline _enqueue_detector_readout_noise!(
    style::AcceleratorStyle,
    detector::Detector{NoisePhoton},
    rng::_PreparedCounterRNG,
) = detector.products.frame

@inline function _enqueue_detector_readout_noise!(
    style::AcceleratorStyle,
    detector::Detector{<:NoiseReadout},
    rng::_PreparedCounterRNG,
)
    sigma = effective_readout_sigma(
        detector.params.sensor,
        detector.noise.sigma,
    )
    return _enqueue_detector_gaussian_noise!(
        style,
        detector,
        rng,
        sigma,
    )
end

@inline function _enqueue_detector_readout_noise!(
    style::AcceleratorStyle,
    detector::Detector{<:NoisePhotonReadout},
    rng::_PreparedCounterRNG,
)
    sigma = effective_readout_sigma(
        detector.params.sensor,
        detector.noise.sigma,
    )
    return _enqueue_detector_gaussian_noise!(
        style,
        detector,
        rng,
        sigma,
    )
end

@inline _enqueue_detector_sensor_readout_noise!(
    style::AcceleratorStyle,
    sensor::Union{CCDSensor,EMCCDSensor},
    detector::Detector,
    rng::_PreparedCounterRNG,
) = detector.products.frame

@inline function _enqueue_detector_sensor_readout_noise!(
    style::AcceleratorStyle,
    sensor::CMOSSensor,
    detector::Detector,
    rng::_PreparedCounterRNG,
)
    frame = detector.products.frame
    scratch = detector.workspace.noise_buffer
    n, m = size(frame)
    if sensor.column_readout_sigma > zero(sensor.column_readout_sigma)
        randn_backend_async!(style, rng, scratch)
        launch_kernel_async!(
            style,
            add_column_noise_kernel!,
            frame,
            scratch,
            sensor.column_readout_sigma,
            n,
            m;
            ndrange=(n, m),
        )
    end
    if sensor.row_readout_sigma > zero(sensor.row_readout_sigma)
        randn_backend_async!(style, rng, scratch)
        launch_kernel_async!(
            style,
            add_row_noise_kernel!,
            frame,
            scratch,
            sensor.row_readout_sigma,
            n,
            m;
            ndrange=(n, m),
        )
    end
    return frame
end

@inline function _enqueue_detector_post_readout_gain!(
    style::AcceleratorStyle,
    sensor::Union{CCDSensor,CMOSSensor},
    detector::Detector,
)
    return _enqueue_detector_scale!(style, detector, detector.params.gain)
end

@inline _enqueue_detector_post_readout_gain!(
    style::AcceleratorStyle,
    sensor::EMCCDSensor,
    detector::Detector,
) = detector.products.frame

@inline function _enqueue_detector_quantization!(
    style::AcceleratorStyle,
    detector::Detector,
)
    bits = detector.params.bits
    bits === nothing && return detector.products.frame
    frame = detector.products.frame
    levels = exp2(eltype(frame)(bits))
    upper = levels - one(levels)
    scale = upper / something(detector.params.full_well)
    launch_kernel_async!(
        style,
        _quantize_enqueued_detector_frame_kernel!,
        frame,
        scale,
        upper,
        length(frame);
        ndrange=length(frame),
    )
    return frame
end

@inline function _enqueue_detector_output!(
    style::AcceleratorStyle,
    ::Nothing,
    frame,
    output,
)
    copyto_backend_async!(output, frame)
    return output
end

@inline function _enqueue_detector_output!(
    style::AcceleratorStyle,
    detector_output::AbstractMatrix{T},
    frame,
    output,
) where {T<:Integer}
    launch_kernel_async!(
        style,
        _write_enqueued_detector_integer_output_kernel!,
        detector_output,
        frame,
        length(frame);
        ndrange=length(frame),
    )
    copyto_backend_async!(output, detector_output)
    return output
end

"""
Queue one admitted, fixed-address frame-detector acquisition without a device
completion boundary or host-state publication. This internal seam is narrower
than `capture!`: callers must first establish
`_supports_enqueued_detector_acquisition` and must publish completion only
after the owning device context reports success.
"""
function _enqueue_detector_acquisition!(
    prepared::PreparedDetectorAcquisition,
    rng::_PreparedCounterRNG,
    output,
)
    detector = prepared.detector
    params = detector.params
    frame = detector.products.frame
    style = execution_style(frame)
    scale = prepared.plan.quantum_efficiency * prepared.plan.rate_scale *
        params.exposure_duration
    launch_kernel_async!(
        style,
        _prepare_enqueued_detector_frame_kernel!,
        frame,
        prepared.input.values,
        scale,
        length(frame);
        ndrange=length(frame),
    )
    _enqueue_detector_photon_noise!(style, detector, rng)
    dark_signal = effective_dark_current(detector) *
        effective_dark_current_duration(
            params.sensor,
            params.exposure_duration,
        )
    _enqueue_detector_poisson_rate!(style, detector, rng, dark_signal)
    _enqueue_detector_sensor_statistics!(style, params.sensor, detector, rng)
    _enqueue_detector_clamp!(style, detector)
    _enqueue_detector_pre_readout_gain!(style, params.sensor, detector, rng)
    _enqueue_detector_readout_noise!(style, detector, rng)
    _enqueue_detector_sensor_readout_noise!(
        style,
        params.sensor,
        detector,
        rng,
    )
    _enqueue_detector_post_readout_gain!(style, params.sensor, detector)
    _enqueue_detector_quantization!(style, detector)
    return _enqueue_detector_output!(
        style,
        detector.products.output_buffer,
        frame,
        output,
    )
end

@inline function _preflight_enqueued_detector_acquisition!(
    prepared::PreparedDetectorAcquisition,
)
    _require_prepared_whole_acquisition(prepared)
    return nothing
end

@inline function _complete_enqueued_detector_acquisition!(
    prepared::PreparedDetectorAcquisition,
)
    state = prepared.detector.state
    state.integrated_time = zero(state.integrated_time)
    state.readout_ready = true
    return nothing
end
