struct CCDSensor{T<:AbstractFloat,M<:FrameSamplingMode} <: FrameSensorType
    clock_induced_charge_per_frame::T
    read_time::T
    sampling_mode::M
end

function CCDSensor(;
    clock_induced_charge_per_frame::Real=0.0,
    read_time::Real=0.0,
    sampling_mode::FrameSamplingMode=SingleRead(),
    T::Type{<:AbstractFloat}=Float64,
)
    clock_induced_charge_per_frame >= 0 || throw(InvalidConfiguration(
        "CCDSensor clock_induced_charge_per_frame must be >= 0"))
    read_time >= 0 ||
        throw(InvalidConfiguration("CCDSensor read_time must be >= 0"))
    mode = validate_ccd_sampling_mode(sampling_mode)
    return CCDSensor{T,typeof(mode)}(
        T(clock_induced_charge_per_frame), T(read_time), mode)
end

validate_ccd_sampling_mode(mode::SingleRead) = validate_frame_sampling_mode(mode)
validate_ccd_sampling_mode(mode::SkipperSampling) = validate_frame_sampling_mode(mode)
validate_ccd_sampling_mode(mode::FrameSamplingMode) = throw(InvalidConfiguration(
    "CCDSensor sampling_mode must be SingleRead or SkipperSampling"))

detector_sensor_symbol(::CCDSensor) = :ccd
supports_clock_induced_charge(::CCDSensor) = true
supports_nondestructive_reads(sensor::CCDSensor) =
    supports_skipper_sampling(sensor.sampling_mode)
supports_multi_read_readout_products(sensor::CCDSensor) =
    supports_skipper_sampling(sensor.sampling_mode)
supports_skipper_sampling(::FrameSamplingMode) = false
supports_skipper_sampling(::SkipperSampling) = true

multi_read_sampling_mode(sensor::CCDSensor) = sensor.sampling_mode
configured_cic_rate(sensor::CCDSensor, ::Type{T}) where {T<:AbstractFloat} =
    T(sensor.clock_induced_charge_per_frame)
effective_readout_sigma(sensor::CCDSensor, sigma) =
    effective_readout_sigma(sensor.sampling_mode, sigma)

sampling_read_time(sensor::CCDSensor, ::Type{T}) where {T<:AbstractFloat} =
    ccd_sampling_read_time(sensor.sampling_mode, sensor.read_time, T)
ccd_sampling_read_time(::SingleRead, read_time, ::Type{T}) where {T<:AbstractFloat} =
    nothing
ccd_sampling_read_time(::SkipperSampling, read_time, ::Type{T}) where {T<:AbstractFloat} =
    T(read_time)

function sampling_wallclock_time(sensor::CCDSensor, integration_time,
    ::Type{T}) where {T<:AbstractFloat}
    return ccd_sampling_wallclock_time(sensor.sampling_mode, integration_time,
        sensor.read_time, T)
end

ccd_sampling_wallclock_time(::SingleRead, integration_time, read_time,
    ::Type{T}) where {T<:AbstractFloat} = T(integration_time)
ccd_sampling_wallclock_time(mode::SkipperSampling, integration_time, read_time,
    ::Type{T}) where {T<:AbstractFloat} =
    T(integration_time) + T(mode.n_samples) * T(read_time)

function apply_sensor_statistics!(sensor::CCDSensor, det::Detector,
    rng::AbstractRNG)
    rate = effective_cic_rate(det)
    rate <= zero(rate) && return det.state.frame
    fill!(det.state.noise_buffer, rate)
    poisson_noise_frame!(det, rng, det.state.noise_buffer)
    det.state.frame .+= det.state.noise_buffer
    return det.state.frame
end

function apply_post_readout_gain!(::CCDSensor, det::Detector)
    det.state.frame .*= det.params.gain
    return det.state.frame
end

function _batched_sensor_statistics!(sensor::CCDSensor, det::Detector,
    cube::AbstractArray, scratch::AbstractArray, rng::AbstractRNG)
    rate = effective_cic_rate(det)
    rate <= zero(rate) && return cube
    fill!(scratch, rate)
    poisson_noise_frame!(det, rng, scratch)
    cube .+= scratch
    return cube
end

_batched_post_readout_gain!(::CCDSensor, det::Detector,
    cube::AbstractArray) = (cube .*= det.params.gain; cube)

_require_batched_sensor_compat(sensor::CCDSensor) =
    require_batched_ccd_sampling(sensor.sampling_mode)
require_batched_ccd_sampling(::SingleRead) = nothing
require_batched_ccd_sampling(::SkipperSampling) = throw(InvalidConfiguration(
    "batched detector capture does not retain Skipper read samples; use repeated capture! calls"))

function detector_readout_products_type(
    ::CCDSensor{<:AbstractFloat,<:SkipperSampling}, frame::A,
    ::Type{T}) where {T<:AbstractFloat,A<:AbstractMatrix{T}}
    return Union{NoFrameReadoutProducts,SkipperReadoutProducts{A}}
end

function _skipper_products_ready(products::SkipperReadoutProducts,
    det::Detector, n_samples::Int)
    return size(products.mean_frame) == readout_product_shape(det) &&
        products.sample_count == n_samples
end

_skipper_products_ready(::FrameReadoutProducts, det::Detector,
    n_samples::Int) = false

function ensure_skipper_products!(det::Detector, n_samples::Int)
    current = readout_products(det)
    _skipper_products_ready(current, det, n_samples) && return current
    mean_frame = similar(det.state.frame, readout_product_shape(det)...)
    products = SkipperReadoutProducts(mean_frame, n_samples)
    det.state.readout_products = products
    return products
end

function _finalize_capture!(sensor::CCDSensor, det::Detector,
    rng::AbstractRNG, exposure_time::Real)
    return finalize_ccd_capture!(sensor.sampling_mode, sensor, det, rng,
        exposure_time)
end

finalize_ccd_capture!(::SingleRead, sensor::CCDSensor, det::Detector,
    rng::AbstractRNG, exposure_time::Real) =
    finalize_readout_pipeline!(det, rng, exposure_time)

function finalize_ccd_capture!(mode::SkipperSampling, sensor::CCDSensor,
    det::Detector, rng::AbstractRNG, exposure_time::Real)
    finalize_charge_generation!(det, rng, exposure_time)
    finalize_charge_transport!(det, rng)

    copyto!(det.state.response_buffer, det.state.frame)
    fill!(det.state.accum_buffer, zero(eltype(det.state.accum_buffer)))
    sigma = _raw_sampling_sigma(det)
    for _ in 1:mode.n_samples
        copyto!(det.state.frame, det.state.response_buffer)
        add_gaussian_noise!(det.state.frame, det, rng, sigma)
        det.state.accum_buffer .+= det.state.frame
    end
    det.state.frame .= det.state.accum_buffer ./ mode.n_samples

    apply_post_readout_gain!(sensor, det)
    apply_readout_correction!(det.params.correction_model, det.state.frame, det)
    apply_quantization!(det)
    subtract_background_map!(det.background_map, det)
    update_sensor_persistence!(sensor, det, exposure_time)

    products = ensure_skipper_products!(det, mode.n_samples)
    _copy_windowed_frame!(products.mean_frame, det.state.frame, det)
    return det.state.frame
end
