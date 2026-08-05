"""
    CCDSensor(; clock_induced_charge_per_frame=0,
        sample_duration=0, sampling_mode=SingleRead(), T=Float64)

CCD area-sensor model. `sample_duration` is the duration of one configured
Skipper nondestructive sample and must remain zero for conventional
`SingleRead` operation. Whole-acquisition readout and readiness timing belongs
to the Plant acquisition definition.
"""
struct CCDSensor{T<:AbstractFloat,M<:FrameSamplingMode} <: AbstractFrameSensor
    clock_induced_charge_per_frame::T
    sample_duration::T
    sampling_mode::M
end

function CCDSensor(;
    clock_induced_charge_per_frame::Real=0.0,
    sample_duration::Real=0.0,
    sampling_mode::FrameSamplingMode=SingleRead(),
    T::Type{<:AbstractFloat}=Float64,
)
    mode = validate_ccd_sampling_mode(sampling_mode)
    cic = T(clock_induced_charge_per_frame)
    duration = T(sample_duration)
    isfinite(cic) && cic >= zero(T) || throw(InvalidConfiguration(
        "CCDSensor clock_induced_charge_per_frame must be finite and >= 0"))
    isfinite(duration) && duration >= zero(T) ||
        throw(InvalidConfiguration(
            "CCDSensor sample_duration must be finite and >= 0"))
    validate_ccd_sample_duration(mode, duration)
    return CCDSensor{T,typeof(mode)}(cic, duration, mode)
end

validate_ccd_sampling_mode(mode::SingleRead) = validate_frame_sampling_mode(mode)
validate_ccd_sampling_mode(mode::SkipperSampling) = validate_frame_sampling_mode(mode)
validate_ccd_sampling_mode(mode::FrameSamplingMode) = throw(InvalidConfiguration(
    "CCDSensor sampling_mode must be SingleRead or SkipperSampling"))

function validate_ccd_sample_duration(::SingleRead, sample_duration)
    iszero(sample_duration) || throw(InvalidConfiguration(
        "CCDSensor sample_duration applies only to SkipperSampling; " *
        "Plant acquisition definitions own single-read readout and readiness timing"))
    return nothing
end

validate_ccd_sample_duration(::SkipperSampling, sample_duration) = nothing

detector_sensor_symbol(::CCDSensor) = :ccd
supports_clock_induced_charge(::CCDSensor) = true
supports_nondestructive_reads(sensor::CCDSensor) =
    supports_skipper_sampling(sensor.sampling_mode)
supports_multi_read_readout_products(sensor::CCDSensor) =
    supports_skipper_sampling(sensor.sampling_mode)
supports_skipper_sampling(::FrameSamplingMode) = false
supports_skipper_sampling(::SkipperSampling) = true

multi_read_sampling_mode(sensor::CCDSensor) = sensor.sampling_mode
configured_cic_per_frame(sensor::CCDSensor,
    ::Type{T}) where {T<:AbstractFloat} =
    T(sensor.clock_induced_charge_per_frame)
effective_readout_sigma(sensor::CCDSensor, sigma) =
    effective_readout_sigma(sensor.sampling_mode, sigma)

sampling_read_time(sensor::CCDSensor, ::Type{T}) where {T<:AbstractFloat} =
    ccd_sampling_read_time(sensor.sampling_mode, sensor.sample_duration, T)
ccd_sampling_read_time(::SingleRead, sample_duration,
    ::Type{T}) where {T<:AbstractFloat} =
    nothing
ccd_sampling_read_time(::SkipperSampling, sample_duration,
    ::Type{T}) where {T<:AbstractFloat} =
    T(sample_duration)

function sampling_wallclock_time(sensor::CCDSensor, integration_time,
    ::Type{T}) where {T<:AbstractFloat}
    return ccd_sampling_wallclock_time(sensor.sampling_mode, integration_time,
        sensor.sample_duration, T)
end

ccd_sampling_wallclock_time(::SingleRead, integration_time, sample_duration,
    ::Type{T}) where {T<:AbstractFloat} = T(integration_time)
ccd_sampling_wallclock_time(mode::SkipperSampling, integration_time,
    sample_duration, ::Type{T}) where {T<:AbstractFloat} =
    T(integration_time) + T(mode.n_samples) * T(sample_duration)

function apply_sensor_statistics!(sensor::CCDSensor, det::Detector,
    rng::AbstractRNG, exposure_time::Real)
    mean_per_frame = effective_cic_per_frame(det)
    mean_per_frame <= zero(mean_per_frame) && return det.products.frame
    fill!(det.workspace.noise_buffer, mean_per_frame)
    poisson_noise_frame!(det, rng, det.workspace.noise_buffer)
    det.products.frame .+= det.workspace.noise_buffer
    return det.products.frame
end

function apply_post_readout_gain!(::CCDSensor, det::Detector)
    det.products.frame .*= det.params.gain
    return det.products.frame
end

function _batched_sensor_statistics!(sensor::CCDSensor, det::Detector,
    cube::AbstractArray, scratch::AbstractArray, rng::AbstractRNG,
    exposure_time::Real)
    mean_per_frame = effective_cic_per_frame(det)
    mean_per_frame <= zero(mean_per_frame) && return cube
    fill!(scratch, mean_per_frame)
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

function detector_readout_workspace_type(
    ::CCDSensor{<:AbstractFloat,<:SkipperSampling}, frame::A,
    ::Type{T}) where {T<:AbstractFloat,A<:AbstractMatrix{T}}
    return Union{NoFrameReadoutWorkspace,SkipperReadoutWorkspace{A}}
end

function _skipper_products_ready(products::SkipperReadoutProducts,
    det::Detector, n_samples::Int)
    return size(products.mean_frame) == readout_product_shape(det) &&
        products.sample_count == n_samples
end

_skipper_products_ready(::FrameReadoutProducts, det::Detector,
    n_samples::Int) = false

function _skipper_workspace_ready(workspace::SkipperReadoutWorkspace,
    det::Detector)
    frame_shape = size(det.products.frame)
    return size(workspace.baseline_frame) == frame_shape &&
        size(workspace.sample_sum) == frame_shape
end

_skipper_workspace_ready(::FrameReadoutWorkspace, ::Detector) = false

function ensure_skipper_products!(det::Detector, n_samples::Int)
    current = readout_products(det)
    _skipper_products_ready(current, det, n_samples) &&
        _skipper_workspace_ready(det.workspace.readout, det) && return current
    mean_frame = similar(det.products.frame, readout_product_shape(det)...)
    baseline_frame = similar(det.products.frame, size(det.products.frame)...)
    sample_sum = similar(det.products.frame, size(det.products.frame)...)
    fill!(mean_frame, zero(eltype(mean_frame)))
    fill!(baseline_frame, zero(eltype(baseline_frame)))
    fill!(sample_sum, zero(eltype(sample_sum)))
    products = SkipperReadoutProducts(mean_frame, n_samples)
    workspace = SkipperReadoutWorkspace(baseline_frame, sample_sum)
    det.products.readout = products
    det.workspace.readout = workspace
    return products
end

@inline function prepare_frame_readout_state!(sensor::CCDSensor,
    det::Detector)
    return prepare_ccd_readout_state!(sensor.sampling_mode, det)
end

@inline prepare_ccd_readout_state!(::SingleRead, det::Detector) = det

@inline function prepare_ccd_readout_state!(mode::SkipperSampling,
    det::Detector)
    ensure_skipper_products!(det, mode.n_samples)
    return det
end

function _finalize_capture!(sensor::CCDSensor, det::Detector,
    rng::AbstractRNG, exposure_time::Real)
    return finalize_ccd_capture!(sensor.sampling_mode, sensor, det, rng,
        exposure_time, exposure_time)
end

function _finalize_incremental_capture!(sensor::CCDSensor, det::Detector,
    rng::AbstractRNG, exposure_time::Real)
    return finalize_ccd_capture!(sensor.sampling_mode, sensor, det, rng,
        exposure_time, zero(exposure_time))
end

finalize_ccd_capture!(::SingleRead, sensor::CCDSensor, det::Detector,
    rng::AbstractRNG, exposure_time::Real, charge_exposure_time::Real) =
    finalize_readout_pipeline!(det, rng, exposure_time,
        charge_exposure_time)

function finalize_ccd_capture!(mode::SkipperSampling, sensor::CCDSensor,
    det::Detector, rng::AbstractRNG, exposure_time::Real,
    charge_exposure_time::Real)
    finalize_charge_generation!(det, rng, charge_exposure_time)
    finalize_charge_transport!(det, rng)

    products = ensure_skipper_products!(det, mode.n_samples)
    workspace = det.workspace.readout
    copyto!(workspace.baseline_frame, det.products.frame)
    fill!(workspace.sample_sum, zero(eltype(workspace.sample_sum)))
    sigma = _raw_sampling_sigma(det)
    for _ in 1:mode.n_samples
        copyto!(det.products.frame, workspace.baseline_frame)
        add_gaussian_noise!(det.products.frame, det, rng, sigma)
        workspace.sample_sum .+= det.products.frame
    end
    det.products.frame .= workspace.sample_sum ./ mode.n_samples

    apply_post_readout_gain!(sensor, det)
    apply_readout_correction!(det.params.correction_model, det.products.frame, det)
    apply_quantization!(det)
    subtract_background_map!(det.background_map, det)
    update_sensor_persistence!(sensor, det, exposure_time)

    _copy_windowed_frame!(products.mean_frame, det.products.frame, det)
    return det.products.frame
end
