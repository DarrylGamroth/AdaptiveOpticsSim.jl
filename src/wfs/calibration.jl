@inline calibration_wavelength(src::AbstractSource, ::Type{T}) where {T<:AbstractFloat} = T(wavelength(src))

@inline lgs_optical_signature(::AbstractSource, sig::UInt) = sig

function lgs_optical_signature(src::LGSSource, sig::UInt)
    params = src.params
    sig = hash(params.altitude, sig)
    sig = hash(params.elongation_factor, sig)
    sig = hash(params.laser_coordinates, sig)
    sig = hash(params.fwhm_spot_up, sig)
    profile = params.na_profile
    if profile !== nothing
        sig = hash(size(profile), sig)
        @inbounds for value in profile
            sig = hash(value, sig)
        end
    end
    return sig
end

@inline lgs_optical_signature(src::SpectralSource, sig::UInt) =
    lgs_optical_signature(src.source, sig)

@inline lgs_optical_signature(src::ExtendedSource, sig::UInt) =
    lgs_optical_signature(src.source, sig)

@inline calibration_signature(src::AbstractSource) =
    lgs_optical_signature(src, source_measurement_signature(src))

@inline pupil_aperture_calibration_signature(pupil::PupilFunction,
    sig::UInt) = hash(aperture_revision(pupil), sig)

"""
    common_wfs_calibration_source(asterism, sensor_name)

Return the source that defines the shared reference calibration for a flat
asterism measurement. Every source must have the same optical calibration
signature; otherwise a reference made from the first source would be
incorrect for one or more accumulated components.
"""
function common_wfs_calibration_source(ast::Asterism,
    sensor_name::AbstractString)
    isempty(ast.sources) && throw(InvalidConfiguration(
        "asterism must contain at least one source"))
    common_source = first(ast.sources)
    common_signature = calibration_signature(common_source)
    @inbounds for i in 2:length(ast.sources)
        calibration_signature(ast.sources[i]) == common_signature || throw(
            InvalidConfiguration(
                "$(sensor_name) asterism sources must share a common optical calibration signature"))
    end
    return common_source
end

@inline function calibration_storage_signature(values::AbstractArray,
    sig::UInt)
    sig = hash(size(values), sig)
    return hash(objectid(values), sig)
end

function detector_response_calibration_signature(
    model::AbstractFrameResponse, ::UInt)
    throw(InvalidConfiguration(
        "detector-coupled WFS calibration does not support response model " *
        "$(typeof(model)); define an instance-complete " *
        "detector_response_calibration_signature overload for this model"))
end

@inline detector_response_calibration_signature(::NullFrameResponse,
    sig::UInt) = hash(NullFrameResponse, sig)

function detector_response_calibration_signature(model::GaussianPixelResponse,
    sig::UInt)
    sig = hash(typeof(model), sig)
    sig = hash(model.response_width_px, sig)
    return calibration_storage_signature(model.kernel, sig)
end

function detector_response_calibration_signature(model::SampledFrameResponse,
    sig::UInt)
    sig = hash(typeof(model), sig)
    return calibration_storage_signature(model.kernel, sig)
end

function detector_response_calibration_signature(model::RectangularPixelAperture,
    sig::UInt)
    sig = hash(typeof(model), sig)
    sig = hash(model.pitch_x_px, sig)
    sig = hash(model.pitch_y_px, sig)
    sig = hash(model.fill_factor_x, sig)
    sig = hash(model.fill_factor_y, sig)
    sig = calibration_storage_signature(model.kernel_x, sig)
    return calibration_storage_signature(model.kernel_y, sig)
end

function detector_qe_calibration_signature(model::ScalarQuantumEfficiency,
    sig::UInt)
    sig = hash(typeof(model), sig)
    return hash(model.value, sig)
end

function detector_qe_calibration_signature(model::SampledQuantumEfficiency,
    sig::UInt)
    sig = hash(typeof(model), sig)
    sig = calibration_storage_signature(model.wavelengths, sig)
    sig = calibration_storage_signature(model.values, sig)
    return hash(model.out_of_band, sig)
end

@inline detector_readout_correction_calibration_signature(
    ::NullFrameReadoutCorrection, sig::UInt) = sig

@inline detector_sampling_calibration_signature(mode::FrameSamplingMode,
    sig::UInt) = hash(typeof(mode), sig)

function detector_sampling_calibration_signature(
    mode::AveragedNonDestructiveReads, sig::UInt)
    sig = hash(typeof(mode), sig)
    return hash(mode.n_reads, sig)
end

function detector_sampling_calibration_signature(mode::FowlerSampling,
    sig::UInt)
    sig = hash(typeof(mode), sig)
    return hash(mode.n_pairs, sig)
end

function detector_sampling_calibration_signature(mode::UpTheRampSampling,
    sig::UInt)
    sig = hash(typeof(mode), sig)
    return hash(mode.n_reads, sig)
end

function detector_sampling_calibration_signature(mode::SkipperSampling,
    sig::UInt)
    sig = hash(typeof(mode), sig)
    return hash(mode.n_samples, sig)
end

@inline detector_sensor_calibration_signature(sensor::FrameSensorType,
    sig::UInt) = hash(typeof(sensor), sig)

function detector_sensor_calibration_signature(sensor::CCDSensor,
    sig::UInt)
    sig = hash(typeof(sensor), sig)
    sig = hash(sensor.read_time, sig)
    return detector_sampling_calibration_signature(sensor.sampling_mode, sig)
end

function detector_sensor_calibration_signature(
    sensor::HgCdTeAvalancheArraySensor, sig::UInt)
    sig = hash(typeof(sensor), sig)
    sig = hash(sensor.avalanche_gain, sig)
    sig = hash(sensor.read_time, sig)
    return detector_sampling_calibration_signature(sensor.sampling_mode, sig)
end

function detector_readout_correction_calibration_signature(
    model::ReferencePixelCommonModeCorrection, sig::UInt)
    sig = hash(typeof(model), sig)
    sig = hash(model.edge_rows, sig)
    return hash(model.edge_cols, sig)
end

function detector_readout_correction_calibration_signature(
    model::ReferenceRowCommonModeCorrection, sig::UInt)
    sig = hash(typeof(model), sig)
    return hash(model.edge_cols, sig)
end

function detector_readout_correction_calibration_signature(
    model::ReferenceColumnCommonModeCorrection, sig::UInt)
    sig = hash(typeof(model), sig)
    return hash(model.edge_rows, sig)
end

function detector_readout_correction_calibration_signature(
    model::ReferenceOutputCommonModeCorrection, sig::UInt)
    sig = hash(typeof(model), sig)
    sig = hash(model.output_cols, sig)
    sig = hash(model.edge_rows, sig)
    return hash(model.edge_cols, sig)
end

@inline detector_readout_correction_calibration_signature(
    ::Tuple{}, sig::UInt) = sig

@inline function detector_readout_correction_calibration_signature(
    stages::Tuple, sig::UInt)
    next_sig = detector_readout_correction_calibration_signature(
        first(stages), sig)
    return detector_readout_correction_calibration_signature(
        Base.tail(stages), next_sig)
end

@inline function detector_readout_correction_calibration_signature(
    model::CompositeFrameReadoutCorrection, sig::UInt)
    sig = hash(typeof(model), sig)
    return detector_readout_correction_calibration_signature(model.stages,
        sig)
end

@inline detector_signal_defect_calibration_signature(
    ::NullDetectorDefectModel, sig::UInt) = sig

@inline detector_signal_defect_calibration_signature(
    ::DarkSignalNonuniformity, sig::UInt) = sig

function detector_signal_defect_calibration_signature(
    model::PixelResponseNonuniformity, sig::UInt)
    sig = hash(typeof(model), sig)
    return calibration_storage_signature(model.gain_map, sig)
end

function detector_signal_defect_calibration_signature(model::BadPixelMask,
    sig::UInt)
    sig = hash(typeof(model), sig)
    sig = hash(model.throughput, sig)
    return calibration_storage_signature(model.mask, sig)
end

@inline detector_signal_defect_calibration_signature(
    ::Tuple{}, sig::UInt) = sig

@inline function detector_signal_defect_calibration_signature(stages::Tuple,
    sig::UInt)
    next_sig = detector_signal_defect_calibration_signature(first(stages),
        sig)
    return detector_signal_defect_calibration_signature(Base.tail(stages),
        next_sig)
end

@inline detector_signal_defect_calibration_signature(
    model::CompositeDetectorDefectModel, sig::UInt) =
    detector_signal_defect_calibration_signature(model.stages, sig)

@inline require_wfs_calibration_defects(::NullDetectorDefectModel) = nothing
@inline require_wfs_calibration_defects(::PixelResponseNonuniformity) = nothing
@inline require_wfs_calibration_defects(::BadPixelMask) = nothing

function require_wfs_calibration_defects(::DarkSignalNonuniformity)
    throw(InvalidConfiguration(
        "detector-coupled WFS calibration does not support dark-signal nonuniformity"))
end

@inline require_wfs_calibration_defects(::Tuple{}) = nothing

@inline function require_wfs_calibration_defects(stages::Tuple)
    require_wfs_calibration_defects(first(stages))
    return require_wfs_calibration_defects(Base.tail(stages))
end

@inline require_wfs_calibration_defects(
    model::CompositeDetectorDefectModel) =
    require_wfs_calibration_defects(model.stages)

@inline require_wfs_calibration_charge_coupling(::NullChargeCoupling) = nothing

function require_wfs_calibration_charge_coupling(
    ::AbstractChargeCouplingModel)
    throw(InvalidConfiguration(
        "detector-coupled WFS calibration does not support post-collection charge coupling"))
end

@inline require_wfs_calibration_background_map(::NoBackground) = nothing

function require_wfs_calibration_background_map(::BackgroundModel)
    throw(InvalidConfiguration(
        "detector-coupled WFS calibration does not support deterministic background-map subtraction"))
end

@inline require_wfs_calibration_correction(
    ::NullFrameReadoutCorrection) = nothing
@inline require_wfs_calibration_correction(
    ::ReferencePixelCommonModeCorrection) = nothing
@inline require_wfs_calibration_correction(
    ::ReferenceRowCommonModeCorrection) = nothing
@inline require_wfs_calibration_correction(
    ::ReferenceColumnCommonModeCorrection) = nothing
@inline require_wfs_calibration_correction(
    ::ReferenceOutputCommonModeCorrection) = nothing
@inline require_wfs_calibration_correction(::Tuple{}) = nothing

@inline function require_wfs_calibration_correction(stages::Tuple)
    require_wfs_calibration_correction(first(stages))
    return require_wfs_calibration_correction(Base.tail(stages))
end

@inline require_wfs_calibration_correction(
    model::CompositeFrameReadoutCorrection) =
    require_wfs_calibration_correction(model.stages)

function require_wfs_calibration_correction(
    ::FrameReadoutCorrectionModel)
    throw(InvalidConfiguration(
        "detector-coupled WFS calibration does not support this readout correction"))
end

@inline require_wfs_calibration_sensor(::CCDSensor) = nothing
@inline require_wfs_calibration_sensor(::InGaAsSensor) = nothing

function require_wfs_calibration_sensor(sensor::CMOSSensor)
    is_null_cmos_output_model(sensor.output_model) || throw(
        InvalidConfiguration(
            "detector-coupled WFS calibration does not support a CMOS output pattern"))
    return nothing
end

@inline require_wfs_hgcdte_sampling(::SingleRead,
    ::HgCdTeAvalancheArraySensor, ::Detector) = nothing
@inline require_wfs_hgcdte_sampling(::AveragedNonDestructiveReads,
    ::HgCdTeAvalancheArraySensor, ::Detector) = nothing
@inline require_wfs_hgcdte_sampling(::CorrelatedDoubleSampling,
    ::HgCdTeAvalancheArraySensor, ::Detector) = nothing
@inline require_wfs_hgcdte_sampling(::FowlerSampling,
    ::HgCdTeAvalancheArraySensor, ::Detector) = nothing

@inline function require_wfs_hgcdte_sampling(mode::UpTheRampSampling,
    sensor::HgCdTeAvalancheArraySensor, det::Detector)
    validate_up_the_ramp_schedule(sensor, det, mode,
        det.params.integration_time)
    return nothing
end

function require_wfs_hgcdte_sampling(::FrameSamplingMode,
    ::HgCdTeAvalancheArraySensor, ::Detector)
    throw(InvalidConfiguration(
        "detector-coupled WFS calibration does not support this HgCdTe sampling mode"))
end

@inline require_wfs_calibration_sensor(
    ::HgCdTeAvalancheArraySensor) = nothing

@inline require_wfs_calibration_schedule(::FrameSensorType,
    ::Detector) = nothing

@inline function require_wfs_calibration_schedule(
    sensor::HgCdTeAvalancheArraySensor, det::Detector)
    require_wfs_hgcdte_sampling(sensor.sampling_mode, sensor, det)
    return nothing
end

@inline require_wfs_emccd_output(::ConventionalOutput) = nothing

function require_wfs_emccd_output(::AbstractEMCCDOutputPath)
    throw(InvalidConfiguration(
        "detector-coupled WFS calibration does not support EM-register output"))
end

@inline require_wfs_emccd_mode(::LinearEMMode) = nothing

function require_wfs_emccd_mode(::AbstractEMCCDOperatingMode)
    throw(InvalidConfiguration(
        "detector-coupled WFS calibration only supports linear EMCCD output"))
end

function require_wfs_calibration_sensor(sensor::EMCCDSensor)
    require_wfs_emccd_output(sensor.output_path)
    require_wfs_emccd_mode(sensor.operating_mode)
    return nothing
end

function require_wfs_calibration_sensor(::FrameSensorType)
    throw(InvalidConfiguration(
        "detector-coupled WFS calibration does not support this frame sensor"))
end

"""
    detector_calibration_signature(detector, signature)

Extend a source calibration signature with the deterministic optical-signal
configuration used before WFS signal extraction. This includes presampling
response, sampling, QE, exposure, binning, PRNU, bad-pixel throughput, and
the configured deterministic gain, sensor sampling/read timing, and supported
homogeneous readout correction. Noise, generated charge,
persistence, and other stateful or stochastic electronics are intentionally
excluded. Parameter arrays are run-owned immutable configuration; rebuild the
model or detector instead of mutating them in place.
"""
@inline function require_wfs_detector_calibration_output(det::Detector)
    params = det.params
    params.readout_window === nothing || throw(InvalidConfiguration(
        "detector-coupled WFS calibration does not support a readout_window"))
    params.output_type === nothing || throw(InvalidConfiguration(
        "detector-coupled WFS calibration requires floating-point detector output"))
    params.bits === nothing || throw(InvalidConfiguration(
        "detector-coupled WFS calibration does not support quantization"))
    params.full_well === nothing || throw(InvalidConfiguration(
        "detector-coupled WFS calibration does not support saturation"))
    is_null_frame_nonlinearity(params.nonlinearity_model) || throw(
        InvalidConfiguration(
            "detector-coupled WFS calibration does not support frame nonlinearity"))
    is_null_persistence(persistence_model(params.sensor)) || throw(
        InvalidConfiguration(
            "detector-coupled WFS calibration does not support persistence"))
    require_wfs_calibration_charge_coupling(params.charge_coupling_model)
    require_wfs_calibration_defects(params.defect_model)
    require_wfs_calibration_background_map(det.background_map)
    require_wfs_calibration_correction(params.correction_model)
    require_wfs_calibration_sensor(params.sensor)
    require_wfs_calibration_schedule(params.sensor, det)
    return nothing
end

function detector_calibration_signature(det::Detector, sig::UInt)
    require_wfs_detector_calibration_output(det)
    params = det.params
    sig = hash(params.integration_time, sig)
    sig = hash(params.psf_sampling, sig)
    sig = hash(params.binning, sig)
    sig = hash(params.gain, sig)
    sig = detector_sensor_calibration_signature(params.sensor, sig)
    sig = detector_qe_calibration_signature(params.quantum_efficiency_model,
        sig)
    sig = detector_response_calibration_signature(params.response_model, sig)
    sig = detector_signal_defect_calibration_signature(params.defect_model,
        sig)
    return detector_readout_correction_calibration_signature(
        params.correction_model, sig)
end

@inline function apply_detector_calibration_readout!(det::Detector,
    frame::AbstractMatrix)
    T = eltype(frame)
    gain = deterministic_frame_readout_gain(det.params.sensor,
        det.params.gain, T)
    gain == one(T) || (frame .*= gain)
    return apply_readout_correction!(det.params.correction_model, frame, det)
end

"""
    detector_calibration_frame!(detector, photon_rate, quantum_efficiency)

Apply the same deterministic optical-signal path as frame acquisition:
presampling response, detector sampling, exposure and QE, detector binning,
PRNU, bad-pixel throughput, deterministic sensor and detector gain, and
supported homogeneous readout correction.
No noise, generated charge, persistence, saturation, quantization, or thermal
state transition is applied. Supported HgCdTe single, averaged,
correlated-double, Fowler, and valid up-the-ramp reads have this same noiseless
transform. Signal/reference sampling collapses to the fully integrated frame;
up-the-ramp calibration validates the read schedule without producing or
mutating acquisition readout products.
"""
function detector_calibration_frame!(det::Detector,
    photon_rate::AbstractMatrix, quantum_efficiency::Real)
    require_wfs_detector_calibration_output(det)
    require_whole_capture_idle(det)
    frame = prepare_signal_frame!(det, photon_rate,
        det.params.integration_time, quantum_efficiency, false,
        det.params.integration_time)
    return apply_detector_calibration_readout!(det, frame)
end

@inline function detector_calibration_frame!(det::Detector,
    photon_rate::AbstractMatrix, src::AbstractSource)
    T = eltype(det.state.frame)
    return detector_calibration_frame!(det, photon_rate,
        effective_qe(det, src, T))
end

@inline calibration_matches(calibrated::Bool, stored_λ, λ) = calibrated && stored_λ == λ

@inline calibration_matches(calibrated::Bool, stored_λ, λ, stored_sig::UInt, sig::UInt) =
    calibrated && stored_λ == λ && stored_sig == sig

function save_zero_opd!(pupil::PupilFunction)
    saved = copy(pupil.opd)
    reset_opd!(pupil)
    return saved
end

@inline function restore_opd!(pupil::PupilFunction,
    saved_opd::AbstractMatrix)
    copyto!(pupil.opd, saved_opd)
    return pupil
end

@inline function store_reference_signal!(reference::AbstractMatrix, signal::AbstractMatrix, slopes::AbstractVector)
    copyto!(reference, signal)
    fill!(slopes, zero(eltype(slopes)))
    return reference
end
