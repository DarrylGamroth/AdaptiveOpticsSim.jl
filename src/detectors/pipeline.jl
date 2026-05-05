abstract type AbstractDetectorExecutionPlan end

struct DetectorDirectPlan <: AbstractDetectorExecutionPlan end
struct DetectorHostMirrorPlan <: AbstractDetectorExecutionPlan end

@inline detector_execution_plan(style::ExecutionStyle, det::Detector) =
    detector_execution_plan(typeof(style), typeof(det))

@inline detector_execution_plan(::Type{<:ExecutionStyle}, ::Type{<:Detector}) = DetectorDirectPlan()

@inline photon_noise_enabled(::Detector{NoiseNone}) = false
@inline photon_noise_enabled(::Detector{NoisePhoton}) = true
@inline photon_noise_enabled(::Detector{<:NoiseReadout}) = false
@inline photon_noise_enabled(::Detector{<:NoisePhotonReadout}) = true

@inline function prepare_signal_frame!(det::Detector, psf::AbstractMatrix, exposure_time::Real, qe=det.params.qe)
    fill_frame!(det, psf, exposure_time, qe)
    apply_signal_defects!(det.params.defect_model, det, exposure_time)
    apply_sensor_persistence!(det.params.sensor, det, exposure_time)
    return det.state.frame
end

function add_poisson_rate!(dest::AbstractMatrix{T}, det::Detector, rng::AbstractRNG, rate) where {T<:AbstractFloat}
    rate_t = T(rate)
    rate_t <= zero(T) && return dest
    fill!(det.state.noise_buffer, rate_t)
    poisson_noise_frame!(det, rng, det.state.noise_buffer)
    dest .+= det.state.noise_buffer
    return dest
end

function add_gaussian_noise!(dest::AbstractMatrix{T}, det::Detector, rng::AbstractRNG, sigma) where {T<:AbstractFloat}
    sigma_t = T(sigma)
    sigma_t <= zero(T) && return dest
    randn_frame_noise!(det, rng, det.state.noise_buffer)
    dest .+= sigma_t .* det.state.noise_buffer
    return dest
end

function capture_signal_pipeline!(det::Detector, psf::AbstractMatrix, rng::AbstractRNG, exposure_time::Real, qe=det.params.qe)
    prepare_signal_frame!(det, psf, exposure_time, qe)
    photon_noise_enabled(det) && poisson_noise_frame!(det, rng, det.state.frame)
    apply_background_flux!(det.background_flux, det, rng, exposure_time)
    return det.state.frame
end

function finalize_readout_pipeline!(det::Detector, rng::AbstractRNG, exposure_time::Real)
    apply_dark_current!(det, rng, exposure_time)
    apply_dark_defects!(det.params.defect_model, det, exposure_time)
    apply_sensor_statistics!(det.params.sensor, det, rng)
    apply_frame_nonlinearity!(det.params.nonlinearity_model, det)
    apply_saturation!(det)
    apply_pre_readout_gain!(det.params.sensor, det, rng)
    apply_readout_noise!(det, rng)
    apply_post_readout_gain!(det.params.sensor, det)
    finalize_readout_products!(det.params.sensor, det, rng, exposure_time)
    apply_readout_correction!(det.params.correction_model, det.state.frame, det)
    apply_quantization!(det)
    subtract_background_map!(det.background_map, det)
    update_sensor_persistence!(det.params.sensor, det, exposure_time)
    return det.state.frame
end
