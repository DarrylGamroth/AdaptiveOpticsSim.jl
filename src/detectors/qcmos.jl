abstract type AbstractQCMOSCameraModel end
abstract type AbstractQCMOSScanMode end
abstract type AbstractQCMOSPhotonNumberQuantization end

struct ORCAQuest <: AbstractQCMOSCameraModel end
struct ORCAQuest2 <: AbstractQCMOSCameraModel end
struct ORCAQuestIQ <: AbstractQCMOSCameraModel end

struct QCMOSStandardScan <: AbstractQCMOSScanMode end
struct QCMOSUltraQuietScan <: AbstractQCMOSScanMode end
struct QCMOSRawScan <: AbstractQCMOSScanMode end

struct NoQCMOSPhotonNumberQuantization <: AbstractQCMOSPhotonNumberQuantization end

struct QCMOSPhotonNumberQuantization{
    T<:AbstractFloat,
    G,
    O,
} <: AbstractQCMOSPhotonNumberQuantization
    gain::G
    offset::O
    max_electrons::Union{Nothing,T}
end

struct QCMOSPhotonNumberResolvingScan{Q<:AbstractQCMOSPhotonNumberQuantization} <: AbstractQCMOSScanMode
    photon_quantization::Q
end

function QCMOSPhotonNumberQuantization(; gain::Union{Real,AbstractMatrix}=1.0,
    offset::Union{Real,AbstractMatrix}=0.0, max_electrons::Union{Nothing,Real}=nothing,
    T::Type{<:AbstractFloat}=Float64, backend::AbstractArrayBackend=CPUBackend())
    backend = _resolve_array_backend(backend)
    converted_gain = _qcmos_photon_calibration_value(gain, T, backend)
    converted_offset = _qcmos_photon_calibration_value(offset, T, backend)
    max_value = max_electrons === nothing ? nothing : T(max_electrons)
    return validate_qcmos_photon_quantization(
        QCMOSPhotonNumberQuantization{T,typeof(converted_gain),typeof(converted_offset)}(
            converted_gain, converted_offset, max_value))
end

function QCMOSPhotonNumberResolvingScan(;
    photon_quantization::Union{Nothing,AbstractQCMOSPhotonNumberQuantization}=nothing,
    gain::Union{Nothing,Real,AbstractMatrix}=nothing,
    offset::Union{Real,AbstractMatrix}=0.0,
    max_electrons::Union{Nothing,Real}=nothing,
    T::Type{<:AbstractFloat}=Float64,
    backend::AbstractArrayBackend=CPUBackend(),
)
    quantization = _resolve_qcmos_photon_quantization(photon_quantization, gain, offset, max_electrons, T, backend)
    converted = convert_qcmos_photon_quantization(quantization, T, backend)
    return QCMOSPhotonNumberResolvingScan(validate_qcmos_photon_quantization(converted))
end

function _resolve_qcmos_photon_quantization(::Nothing, ::Nothing, offset, max_electrons, ::Type{T}, backend) where {T<:AbstractFloat}
    _is_default_qcmos_photon_offset(offset) ||
        throw(InvalidConfiguration("QCMOSPhotonNumberResolvingScan offset requires photon-number gain calibration"))
    max_electrons === nothing ||
        throw(InvalidConfiguration("QCMOSPhotonNumberResolvingScan max_electrons requires photon-number gain calibration"))
    return NoQCMOSPhotonNumberQuantization()
end

function _resolve_qcmos_photon_quantization(::Nothing, gain::Union{Real,AbstractMatrix}, offset, max_electrons,
    ::Type{T}, backend) where {T<:AbstractFloat}
    quantization = QCMOSPhotonNumberQuantization(; gain, offset, max_electrons, T, backend)
    return quantization
end

function _resolve_qcmos_photon_quantization(photon_quantization::AbstractQCMOSPhotonNumberQuantization, gain, offset,
    max_electrons, ::Type{T}, backend) where {T<:AbstractFloat}
    gain === nothing ||
        throw(InvalidConfiguration("QCMOSPhotonNumberResolvingScan cannot combine photon_quantization with gain"))
    _is_default_qcmos_photon_offset(offset) ||
        throw(InvalidConfiguration("QCMOSPhotonNumberResolvingScan cannot combine photon_quantization with offset"))
    max_electrons === nothing ||
        throw(InvalidConfiguration("QCMOSPhotonNumberResolvingScan cannot combine photon_quantization with max_electrons"))
    return photon_quantization
end

_is_default_qcmos_photon_offset(offset::Real) = iszero(offset)
_is_default_qcmos_photon_offset(offset::AbstractMatrix) = all(iszero, offset)

struct QCMOSSensor{
    T<:AbstractFloat,
    C<:AbstractQCMOSCameraModel,
    S<:AbstractQCMOSScanMode,
    O<:AbstractCMOSOutputModel,
    M<:AbstractFrameTimingModel,
} <: FrameSensorType
    camera_model::C
    scan_mode::S
    readout_noise_sigma::T
    readout_noise_median::T
    dsnu_sigma::T
    prnu_sigma::T
    linearity_error::T
    conversion_factor_e_per_count::T
    pixel_size_um::T
    peak_qe::T
    full_well::T
    dark_current::T
    frame_size::Tuple{Int,Int}
    frame_rate_hz::T
    min_exposure_time::T
    output_model::O
    timing_model::M
end

qcmos_camera_symbol(::ORCAQuest) = :orca_quest
qcmos_camera_symbol(::ORCAQuest2) = :orca_quest_2
qcmos_camera_symbol(::ORCAQuestIQ) = :orca_quest_iq

qcmos_scan_mode_symbol(::QCMOSStandardScan) = :standard
qcmos_scan_mode_symbol(::QCMOSUltraQuietScan) = :ultra_quiet
qcmos_scan_mode_symbol(::QCMOSPhotonNumberResolvingScan) = :photon_number_resolving
qcmos_scan_mode_symbol(::QCMOSRawScan) = :raw
qcmos_photon_quantization(sensor::QCMOSSensor) = qcmos_photon_quantization(sensor.scan_mode)
qcmos_photon_quantization(::AbstractQCMOSScanMode) = NoQCMOSPhotonNumberQuantization()
qcmos_photon_quantization(mode::QCMOSPhotonNumberResolvingScan) = mode.photon_quantization

detector_sensor_symbol(::QCMOSSensor) = :qcmos
supports_detector_defect_maps(::QCMOSSensor) = true
supports_shutter_timing(::QCMOSSensor) = true
supports_photon_number_resolving(::QCMOSSensor) = true
supports_calibrated_photon_number_output(sensor::QCMOSSensor) =
    supports_calibrated_photon_number_output(qcmos_photon_quantization(sensor))
supports_calibrated_photon_number_output(::AbstractQCMOSPhotonNumberQuantization) = false
supports_calibrated_photon_number_output(::QCMOSPhotonNumberQuantization) = true
supports_raw_digital_output(sensor::QCMOSSensor) = supports_raw_digital_output(sensor.scan_mode)
supports_raw_digital_output(::AbstractQCMOSScanMode) = false
supports_raw_digital_output(::QCMOSRawScan) = true

qcmos_camera_model(sensor::QCMOSSensor) = sensor.camera_model
qcmos_scan_mode(sensor::QCMOSSensor) = sensor.scan_mode
qcmos_readout_noise_sigma(sensor::QCMOSSensor) = sensor.readout_noise_sigma
qcmos_readout_noise_median(sensor::QCMOSSensor) = sensor.readout_noise_median
qcmos_dsnu_sigma(sensor::QCMOSSensor) = sensor.dsnu_sigma
qcmos_prnu_sigma(sensor::QCMOSSensor) = sensor.prnu_sigma
qcmos_linearity_error(sensor::QCMOSSensor) = sensor.linearity_error
qcmos_conversion_factor_e_per_count(sensor::QCMOSSensor) = sensor.conversion_factor_e_per_count
qcmos_pixel_size_um(sensor::QCMOSSensor) = sensor.pixel_size_um
qcmos_peak_qe(sensor::QCMOSSensor) = sensor.peak_qe
qcmos_full_well(sensor::QCMOSSensor) = sensor.full_well
qcmos_dark_current(sensor::QCMOSSensor) = sensor.dark_current
qcmos_frame_size(sensor::QCMOSSensor) = sensor.frame_size
qcmos_frame_rate_hz(sensor::QCMOSSensor) = sensor.frame_rate_hz
qcmos_min_exposure_time(sensor::QCMOSSensor) = sensor.min_exposure_time

_qcmos_photon_calibration_value(value::Real, ::Type{T}, backend) where {T<:AbstractFloat} = T(value)
_qcmos_photon_calibration_value(value::AbstractMatrix, ::Type{T}, backend) where {T<:AbstractFloat} =
    _to_backend_matrix(T.(Array(value)), backend)

convert_qcmos_scan_mode(mode::AbstractQCMOSScanMode, ::Type{T}, backend) where {T<:AbstractFloat} = mode

function convert_qcmos_scan_mode(mode::QCMOSPhotonNumberResolvingScan, ::Type{T}, backend) where {T<:AbstractFloat}
    quantization = convert_qcmos_photon_quantization(mode.photon_quantization, T, backend)
    return QCMOSPhotonNumberResolvingScan(validate_qcmos_photon_quantization(quantization))
end

convert_qcmos_photon_quantization(::NoQCMOSPhotonNumberQuantization, ::Type{T}, backend) where {T<:AbstractFloat} =
    NoQCMOSPhotonNumberQuantization()

function convert_qcmos_photon_quantization(model::QCMOSPhotonNumberQuantization, ::Type{T}, backend) where {T<:AbstractFloat}
    gain = _qcmos_photon_calibration_value(model.gain, T, backend)
    offset = _qcmos_photon_calibration_value(model.offset, T, backend)
    max_value = model.max_electrons === nothing ? nothing : T(model.max_electrons)
    return QCMOSPhotonNumberQuantization{T,typeof(gain),typeof(offset)}(gain, offset, max_value)
end

validate_qcmos_photon_quantization(::NoQCMOSPhotonNumberQuantization) = NoQCMOSPhotonNumberQuantization()

function validate_qcmos_photon_quantization(model::QCMOSPhotonNumberQuantization)
    _validate_qcmos_photon_gain(model.gain)
    _validate_qcmos_photon_offset(model.offset)
    model.max_electrons === nothing || model.max_electrons >= zero(model.max_electrons) ||
        throw(InvalidConfiguration("QCMOSPhotonNumberQuantization max_electrons must be >= 0"))
    return model
end

function _validate_qcmos_photon_gain(gain::Real)
    gain > zero(gain) || throw(InvalidConfiguration("QCMOSPhotonNumberQuantization gain must be > 0"))
    isfinite(gain) || throw(InvalidConfiguration("QCMOSPhotonNumberQuantization gain must be finite"))
    return nothing
end

function _validate_qcmos_photon_gain(gain::AbstractMatrix)
    isempty(gain) && throw(InvalidConfiguration("QCMOSPhotonNumberQuantization gain map must not be empty"))
    minimum(Array(gain)) > zero(eltype(gain)) ||
        throw(InvalidConfiguration("QCMOSPhotonNumberQuantization gain map values must be > 0"))
    all(isfinite, Array(gain)) ||
        throw(InvalidConfiguration("QCMOSPhotonNumberQuantization gain map values must be finite"))
    return nothing
end

function _validate_qcmos_photon_offset(offset::Real)
    isfinite(offset) || throw(InvalidConfiguration("QCMOSPhotonNumberQuantization offset must be finite"))
    return nothing
end

function _validate_qcmos_photon_offset(offset::AbstractMatrix)
    isempty(offset) && throw(InvalidConfiguration("QCMOSPhotonNumberQuantization offset map must not be empty"))
    all(isfinite, Array(offset)) ||
        throw(InvalidConfiguration("QCMOSPhotonNumberQuantization offset map values must be finite"))
    return nothing
end

qcmos_frame_size(::AbstractQCMOSCameraModel) = (2304, 4096)
qcmos_pixel_size_um(::AbstractQCMOSCameraModel, ::Type{T}) where {T<:AbstractFloat} = T(4.6)
qcmos_peak_qe(::AbstractQCMOSCameraModel, ::Type{T}) where {T<:AbstractFloat} = T(0.85)
qcmos_full_well(::AbstractQCMOSCameraModel, ::Type{T}) where {T<:AbstractFloat} = T(7000)
qcmos_dsnu_sigma(::AbstractQCMOSCameraModel, ::Type{T}) where {T<:AbstractFloat} = T(0.06)
qcmos_prnu_sigma(::AbstractQCMOSCameraModel, ::Type{T}) where {T<:AbstractFloat} = T(0.001)
qcmos_linearity_error(::AbstractQCMOSCameraModel, ::Type{T}) where {T<:AbstractFloat} = T(0.005)
qcmos_conversion_factor_e_per_count(::AbstractQCMOSCameraModel, ::Type{T}) where {T<:AbstractFloat} = T(0.107)

qcmos_readout_noise_sigma(::ORCAQuest, ::QCMOSStandardScan, ::Type{T}) where {T<:AbstractFloat} = T(0.43)
qcmos_readout_noise_sigma(::ORCAQuest, ::QCMOSUltraQuietScan, ::Type{T}) where {T<:AbstractFloat} = T(0.27)
qcmos_readout_noise_sigma(::ORCAQuest, ::QCMOSPhotonNumberResolvingScan, ::Type{T}) where {T<:AbstractFloat} = T(0.27)
qcmos_readout_noise_sigma(::ORCAQuest, ::QCMOSRawScan, ::Type{T}) where {T<:AbstractFloat} = T(0.27)

qcmos_readout_noise_sigma(::ORCAQuest2, ::QCMOSStandardScan, ::Type{T}) where {T<:AbstractFloat} = T(0.43)
qcmos_readout_noise_sigma(::ORCAQuest2, ::QCMOSUltraQuietScan, ::Type{T}) where {T<:AbstractFloat} = T(0.30)
qcmos_readout_noise_sigma(::ORCAQuest2, ::QCMOSPhotonNumberResolvingScan, ::Type{T}) where {T<:AbstractFloat} = T(0.30)
qcmos_readout_noise_sigma(::ORCAQuest2, ::QCMOSRawScan, ::Type{T}) where {T<:AbstractFloat} = T(0.30)

qcmos_readout_noise_sigma(::ORCAQuestIQ, mode::AbstractQCMOSScanMode, ::Type{T}) where {T<:AbstractFloat} =
    qcmos_readout_noise_sigma(ORCAQuest2(), mode, T)

qcmos_readout_noise_median(camera::AbstractQCMOSCameraModel, mode::AbstractQCMOSScanMode, ::Type{T}) where {T<:AbstractFloat} =
    qcmos_readout_noise_sigma(camera, mode, T)
qcmos_readout_noise_median(::ORCAQuest2, ::QCMOSStandardScan, ::Type{T}) where {T<:AbstractFloat} = T(0.39)
qcmos_readout_noise_median(::ORCAQuest2, ::QCMOSUltraQuietScan, ::Type{T}) where {T<:AbstractFloat} = T(0.25)
qcmos_readout_noise_median(::ORCAQuest2, ::QCMOSPhotonNumberResolvingScan, ::Type{T}) where {T<:AbstractFloat} = T(0.25)
qcmos_readout_noise_median(::ORCAQuest2, ::QCMOSRawScan, ::Type{T}) where {T<:AbstractFloat} = T(0.25)
qcmos_readout_noise_median(::ORCAQuestIQ, mode::AbstractQCMOSScanMode, ::Type{T}) where {T<:AbstractFloat} =
    qcmos_readout_noise_median(ORCAQuest2(), mode, T)

qcmos_dark_current(::ORCAQuest, ::Type{T}) where {T<:AbstractFloat} = T(0.016)
qcmos_dark_current(::ORCAQuest2, ::Type{T}) where {T<:AbstractFloat} = T(0.016)
qcmos_dark_current(::ORCAQuestIQ, ::Type{T}) where {T<:AbstractFloat} = T(0.032)

qcmos_frame_rate_hz(::ORCAQuest, ::QCMOSStandardScan, ::Type{T}) where {T<:AbstractFloat} = T(120)
qcmos_frame_rate_hz(::ORCAQuest, ::QCMOSUltraQuietScan, ::Type{T}) where {T<:AbstractFloat} = T(5)
qcmos_frame_rate_hz(::ORCAQuest, ::QCMOSPhotonNumberResolvingScan, ::Type{T}) where {T<:AbstractFloat} = T(5)
qcmos_frame_rate_hz(::ORCAQuest, ::QCMOSRawScan, ::Type{T}) where {T<:AbstractFloat} = T(5)

qcmos_frame_rate_hz(::ORCAQuest2, ::QCMOSStandardScan, ::Type{T}) where {T<:AbstractFloat} = T(120)
qcmos_frame_rate_hz(::ORCAQuest2, ::QCMOSUltraQuietScan, ::Type{T}) where {T<:AbstractFloat} = T(25.4)
qcmos_frame_rate_hz(::ORCAQuest2, ::QCMOSPhotonNumberResolvingScan, ::Type{T}) where {T<:AbstractFloat} = T(25.4)
qcmos_frame_rate_hz(::ORCAQuest2, ::QCMOSRawScan, ::Type{T}) where {T<:AbstractFloat} = T(25.4)
qcmos_frame_rate_hz(::ORCAQuestIQ, mode::AbstractQCMOSScanMode, ::Type{T}) where {T<:AbstractFloat} =
    qcmos_frame_rate_hz(ORCAQuest2(), mode, T)

qcmos_min_exposure_time(::ORCAQuest, ::QCMOSStandardScan, ::Type{T}) where {T<:AbstractFloat} = T(7.2e-6)
qcmos_min_exposure_time(::ORCAQuest, ::QCMOSUltraQuietScan, ::Type{T}) where {T<:AbstractFloat} = T(199.9e-3)
qcmos_min_exposure_time(::ORCAQuest, ::QCMOSPhotonNumberResolvingScan, ::Type{T}) where {T<:AbstractFloat} = T(199.9e-3)
qcmos_min_exposure_time(::ORCAQuest, ::QCMOSRawScan, ::Type{T}) where {T<:AbstractFloat} = T(199.9e-3)

qcmos_min_exposure_time(::ORCAQuest2, ::QCMOSStandardScan, ::Type{T}) where {T<:AbstractFloat} = T(7.2e-6)
qcmos_min_exposure_time(::ORCAQuest2, ::QCMOSUltraQuietScan, ::Type{T}) where {T<:AbstractFloat} = T(33.9e-6)
qcmos_min_exposure_time(::ORCAQuest2, ::QCMOSPhotonNumberResolvingScan, ::Type{T}) where {T<:AbstractFloat} = T(33.9e-6)
qcmos_min_exposure_time(::ORCAQuest2, ::QCMOSRawScan, ::Type{T}) where {T<:AbstractFloat} = T(33.9e-6)
qcmos_min_exposure_time(::ORCAQuestIQ, mode::AbstractQCMOSScanMode, ::Type{T}) where {T<:AbstractFloat} =
    qcmos_min_exposure_time(ORCAQuest2(), mode, T)

qcmos_line_time(camera::AbstractQCMOSCameraModel, scan_mode::AbstractQCMOSScanMode, ::Type{T}) where {T<:AbstractFloat} =
    qcmos_min_exposure_time(camera, scan_mode, T)
qcmos_line_time(::ORCAQuest, ::QCMOSUltraQuietScan, ::Type{T}) where {T<:AbstractFloat} = T(172.8e-6)
qcmos_line_time(::ORCAQuest, ::QCMOSPhotonNumberResolvingScan, ::Type{T}) where {T<:AbstractFloat} = T(172.8e-6)
qcmos_line_time(::ORCAQuest, ::QCMOSRawScan, ::Type{T}) where {T<:AbstractFloat} = T(172.8e-6)
qcmos_default_timing_model(camera::AbstractQCMOSCameraModel, scan_mode::AbstractQCMOSScanMode, ::Type{T}) where {T<:AbstractFloat} =
    RollingShutter{T}(qcmos_line_time(camera, scan_mode, T), 2)

function QCMOSSensor(;
    camera_model::AbstractQCMOSCameraModel=ORCAQuest2(),
    scan_mode::AbstractQCMOSScanMode=QCMOSUltraQuietScan(),
    readout_noise_sigma::Union{Nothing,Real}=nothing,
    readout_noise_median::Union{Nothing,Real}=nothing,
    dsnu_sigma::Union{Nothing,Real}=nothing,
    prnu_sigma::Union{Nothing,Real}=nothing,
    linearity_error::Union{Nothing,Real}=nothing,
    conversion_factor_e_per_count::Union{Nothing,Real}=nothing,
    pixel_size_um::Union{Nothing,Real}=nothing,
    peak_qe::Union{Nothing,Real}=nothing,
    full_well::Union{Nothing,Real}=nothing,
    dark_current::Union{Nothing,Real}=nothing,
    frame_size::Union{Nothing,Tuple{Int,Int}}=nothing,
    frame_rate_hz::Union{Nothing,Real}=nothing,
    min_exposure_time::Union{Nothing,Real}=nothing,
    output_model::AbstractCMOSOutputModel=NullCMOSOutputModel(),
    timing_model::Union{Nothing,AbstractFrameTimingModel}=nothing,
    T::Type{<:AbstractFloat}=Float64,
    backend::AbstractArrayBackend=CPUBackend(),
)
    backend = _resolve_array_backend(backend)
    resolved_readout_noise = T(something(readout_noise_sigma, qcmos_readout_noise_sigma(camera_model, scan_mode, T)))
    resolved_readout_median = T(something(readout_noise_median, qcmos_readout_noise_median(camera_model, scan_mode, T)))
    resolved_dsnu = T(something(dsnu_sigma, qcmos_dsnu_sigma(camera_model, T)))
    resolved_prnu = T(something(prnu_sigma, qcmos_prnu_sigma(camera_model, T)))
    resolved_linearity = T(something(linearity_error, qcmos_linearity_error(camera_model, T)))
    resolved_conversion = T(something(conversion_factor_e_per_count, qcmos_conversion_factor_e_per_count(camera_model, T)))
    resolved_pixel_size = T(something(pixel_size_um, qcmos_pixel_size_um(camera_model, T)))
    resolved_peak_qe = T(something(peak_qe, qcmos_peak_qe(camera_model, T)))
    resolved_full_well = T(something(full_well, qcmos_full_well(camera_model, T)))
    resolved_dark_current = T(something(dark_current, qcmos_dark_current(camera_model, T)))
    resolved_frame_size = something(frame_size, qcmos_frame_size(camera_model))
    resolved_frame_rate = T(something(frame_rate_hz, qcmos_frame_rate_hz(camera_model, scan_mode, T)))
    resolved_min_exposure = T(something(min_exposure_time, qcmos_min_exposure_time(camera_model, scan_mode, T)))
    resolved_readout_noise >= zero(T) || throw(InvalidConfiguration("QCMOSSensor readout_noise_sigma must be >= 0"))
    resolved_readout_median >= zero(T) || throw(InvalidConfiguration("QCMOSSensor readout_noise_median must be >= 0"))
    resolved_dsnu >= zero(T) || throw(InvalidConfiguration("QCMOSSensor dsnu_sigma must be >= 0"))
    resolved_prnu >= zero(T) || throw(InvalidConfiguration("QCMOSSensor prnu_sigma must be >= 0"))
    resolved_linearity >= zero(T) || throw(InvalidConfiguration("QCMOSSensor linearity_error must be >= 0"))
    resolved_conversion > zero(T) || throw(InvalidConfiguration("QCMOSSensor conversion_factor_e_per_count must be > 0"))
    resolved_pixel_size > zero(T) || throw(InvalidConfiguration("QCMOSSensor pixel_size_um must be > 0"))
    zero(T) <= resolved_peak_qe <= one(T) || throw(InvalidConfiguration("QCMOSSensor peak_qe must lie in [0, 1]"))
    resolved_full_well > zero(T) || throw(InvalidConfiguration("QCMOSSensor full_well must be > 0"))
    resolved_dark_current >= zero(T) || throw(InvalidConfiguration("QCMOSSensor dark_current must be >= 0"))
    all(>(0), resolved_frame_size) || throw(InvalidConfiguration("QCMOSSensor frame_size dimensions must be > 0"))
    resolved_frame_rate > zero(T) || throw(InvalidConfiguration("QCMOSSensor frame_rate_hz must be > 0"))
    resolved_min_exposure >= zero(T) || throw(InvalidConfiguration("QCMOSSensor min_exposure_time must be >= 0"))
    converted_scan_mode = convert_qcmos_scan_mode(scan_mode, T, backend)
    converted_output = convert_cmos_output_model(output_model, T, backend)
    validated_output = validate_cmos_output_model(converted_output)
    resolved_timing = something(timing_model, qcmos_default_timing_model(camera_model, scan_mode, T))
    validated_timing = validate_frame_timing_model(convert_frame_timing_model(resolved_timing, T))
    return QCMOSSensor{
        T,
        typeof(camera_model),
        typeof(converted_scan_mode),
        typeof(validated_output),
        typeof(validated_timing),
    }(
        camera_model,
        converted_scan_mode,
        resolved_readout_noise,
        resolved_readout_median,
        resolved_dsnu,
        resolved_prnu,
        resolved_linearity,
        resolved_conversion,
        resolved_pixel_size,
        resolved_peak_qe,
        resolved_full_well,
        resolved_dark_current,
        resolved_frame_size,
        resolved_frame_rate,
        resolved_min_exposure,
        validated_output,
        validated_timing,
    )
end

ORCAQuestSensor(; kwargs...) = QCMOSSensor(; camera_model=ORCAQuest(), kwargs...)
ORCAQuest2Sensor(; kwargs...) = QCMOSSensor(; camera_model=ORCAQuest2(), kwargs...)
ORCAQuestIQSensor(; kwargs...) = QCMOSSensor(; camera_model=ORCAQuestIQ(), kwargs...)

default_frame_timing_model(sensor::QCMOSSensor; T::Type{<:AbstractFloat}=Float64) = sensor.timing_model
default_response_model(::QCMOSSensor; T::Type{<:AbstractFloat}=Float64, backend::AbstractArrayBackend=CPUBackend()) =
    RectangularPixelAperture(T=T, backend=backend)

function sampling_wallclock_time(sensor::QCMOSSensor, integration_time, frame_size::Tuple{Int,Int},
    window::Union{Nothing,FrameWindow}, ::Type{T}) where {T<:AbstractFloat}
    readout_time = inv(T(sensor.frame_rate_hz))
    if window !== nothing
        readout_time *= T(length(window.rows)) / T(sensor.frame_size[1])
    end
    return max(T(integration_time), readout_time)
end

function apply_post_readout_gain!(::QCMOSSensor, det::Detector)
    det.state.frame .*= det.params.gain
    apply_output_model!(execution_style(det.state.frame), det.params.sensor.output_model, det.state.frame)
    apply_qcmos_photon_number_quantization!(det.params.sensor.scan_mode, det.state.frame)
    return det.state.frame
end

apply_qcmos_photon_number_quantization!(::AbstractQCMOSScanMode, frame) = frame
apply_qcmos_photon_number_quantization!(mode::QCMOSPhotonNumberResolvingScan, frame) =
    apply_qcmos_photon_number_quantization!(mode.photon_quantization, frame)
apply_qcmos_photon_number_quantization!(::NoQCMOSPhotonNumberQuantization, frame) = frame

function apply_qcmos_photon_number_quantization!(model::QCMOSPhotonNumberQuantization, frame)
    _validate_qcmos_photon_calibration_shape(model.gain, frame, "gain")
    _validate_qcmos_photon_calibration_shape(model.offset, frame, "offset")
    @. frame = round(max(frame * model.gain + model.offset, zero(frame)))
    _apply_qcmos_photon_number_limit!(model.max_electrons, frame)
    return frame
end

_validate_qcmos_photon_calibration_shape(value::Real, frame, name) = nothing

function _validate_qcmos_photon_calibration_shape(value::AbstractMatrix, frame::AbstractMatrix, name)
    size(value) == size(frame) ||
        throw(DimensionMismatchError("QCMOSPhotonNumberQuantization $(name) map must match detector frame size"))
    return nothing
end

function _validate_qcmos_photon_calibration_shape(value::AbstractMatrix, cube::AbstractArray{T,3}, name) where {T}
    size(value) == (size(cube, 2), size(cube, 3)) ||
        throw(DimensionMismatchError("QCMOSPhotonNumberQuantization $(name) map must match detector frame size"))
    return nothing
end

_apply_qcmos_photon_number_limit!(::Nothing, frame) = frame

function _apply_qcmos_photon_number_limit!(limit, frame)
    clamp!(frame, zero(eltype(frame)), eltype(frame)(limit))
    return frame
end

function _require_batched_sensor_compat(sensor::QCMOSSensor)
    is_null_cmos_output_model(sensor.output_model) ||
        throw(InvalidConfiguration("batched detector capture does not yet support QCMOSSensor output-group patterns"))
    is_global_shutter(sensor.timing_model) ||
        throw(InvalidConfiguration("batched detector capture does not yet support QCMOSSensor rolling-shutter timing"))
    return nothing
end

function _batched_post_readout_gain!(::QCMOSSensor, det::Detector, cube::AbstractArray)
    isone(det.params.gain) || (cube .*= det.params.gain)
    apply_qcmos_photon_number_quantization!(det.params.sensor.scan_mode, cube)
    return cube
end

function qcmos_snr(signal_photons::Real; quantum_efficiency::Real=0.85,
    dark_electrons::Real=0, readout_noise::Real=0.3, T::Type{<:AbstractFloat}=Float64)
    signal = T(signal_photons)
    qe = T(quantum_efficiency)
    dark = T(dark_electrons)
    read_noise = T(readout_noise)
    signal >= zero(T) || throw(InvalidConfiguration("qcmos_snr signal_photons must be >= 0"))
    zero(T) <= qe <= one(T) || throw(InvalidConfiguration("qcmos_snr quantum_efficiency must lie in [0, 1]"))
    dark >= zero(T) || throw(InvalidConfiguration("qcmos_snr dark_electrons must be >= 0"))
    read_noise >= zero(T) || throw(InvalidConfiguration("qcmos_snr readout_noise must be >= 0"))
    detected_signal = qe * signal
    variance = detected_signal + dark + read_noise^2
    return _qcmos_snr_ratio(detected_signal, variance)
end

function _qcmos_snr_ratio(signal, variance)
    variance > zero(variance) && return signal / sqrt(variance)
    return signal > zero(signal) ? oftype(signal / one(signal), Inf) : zero(signal)
end

function qcmos_snr(sensor::QCMOSSensor, signal_photons::Real; exposure_time::Real=1,
    dark_electrons::Union{Nothing,Real}=nothing, readout_noise::Union{Nothing,Real}=nothing,
    quantum_efficiency::Union{Nothing,Real}=nothing)
    T = typeof(sensor.readout_noise_sigma)
    dark = something(dark_electrons, sensor.dark_current * T(exposure_time))
    sigma = something(readout_noise, sensor.readout_noise_sigma)
    qe = something(quantum_efficiency, sensor.peak_qe)
    return qcmos_snr(signal_photons; quantum_efficiency=qe, dark_electrons=dark, readout_noise=sigma, T=T)
end

function qcmos_relative_snr(signal_photons::Real; kwargs...)
    signal_photons >= 0 || throw(InvalidConfiguration("qcmos_relative_snr signal_photons must be >= 0"))
    signal_photons == 0 && return zero(float(signal_photons))
    return qcmos_snr(signal_photons; kwargs...) / sqrt(float(signal_photons))
end

function qcmos_relative_snr(sensor::QCMOSSensor, signal_photons::Real; kwargs...)
    signal_photons >= 0 || throw(InvalidConfiguration("qcmos_relative_snr signal_photons must be >= 0"))
    signal_photons == 0 && return zero(typeof(sensor.readout_noise_sigma))
    return qcmos_snr(sensor, signal_photons; kwargs...) / sqrt(typeof(sensor.readout_noise_sigma)(signal_photons))
end

function _default_qcmos_noise(sensor::QCMOSSensor, noise::Nothing)
    return NoisePhotonReadout(sensor.readout_noise_sigma)
end

_default_qcmos_noise(::QCMOSSensor, noise::NoiseModel) = noise
_qcmos_float_type(::QCMOSSensor{T}) where {T<:AbstractFloat} = T

function QCMOSDetector(;
    sensor::Union{Nothing,QCMOSSensor}=nothing,
    integration_time::Union{Nothing,Real}=nothing,
    noise::Union{Nothing,NoiseModel}=nothing,
    qe::Union{Nothing,Real}=nothing,
    dark_current::Union{Nothing,Real}=nothing,
    full_well::Union{Nothing,Real}=nothing,
    bits::Union{Nothing,Int,Symbol}=:default,
    output_type::Union{Nothing,DataType}=UInt16,
    T::Union{Nothing,Type{<:AbstractFloat}}=nothing,
    backend::AbstractArrayBackend=CPUBackend(),
    kwargs...,
)
    detector_T = something(T, sensor === nothing ? Float64 : _qcmos_float_type(sensor))
    resolved_sensor = sensor === nothing ? QCMOSSensor(T=detector_T, backend=backend) : sensor
    resolved_integration_time = something(integration_time, max(resolved_sensor.min_exposure_time, inv(resolved_sensor.frame_rate_hz)))
    resolved_bits = resolve_qcmos_detector_bits(resolved_sensor.scan_mode, bits)
    return Detector(;
        integration_time=resolved_integration_time,
        noise=_default_qcmos_noise(resolved_sensor, noise),
        qe=something(qe, resolved_sensor.peak_qe),
        dark_current=something(dark_current, resolved_sensor.dark_current),
        full_well=something(full_well, resolved_sensor.full_well),
        bits=resolved_bits,
        output_type=output_type,
        sensor=resolved_sensor,
        T=detector_T,
        backend=backend,
        kwargs...,
    )
end

resolve_qcmos_detector_bits(mode::AbstractQCMOSScanMode, bits::Nothing) = nothing
resolve_qcmos_detector_bits(mode::AbstractQCMOSScanMode, bits::Integer) = Int(bits)

function resolve_qcmos_detector_bits(mode::AbstractQCMOSScanMode, bits::Symbol)
    bits === :default || throw(InvalidConfiguration("QCMOSDetector bits symbol must be :default"))
    return default_qcmos_detector_bits(mode)
end

default_qcmos_detector_bits(::AbstractQCMOSScanMode) = 16
default_qcmos_detector_bits(mode::QCMOSPhotonNumberResolvingScan) =
    default_qcmos_detector_bits(mode.photon_quantization)
default_qcmos_detector_bits(::NoQCMOSPhotonNumberQuantization) = 16
default_qcmos_detector_bits(::QCMOSPhotonNumberQuantization) = nothing

function _qcmos_detector_for_camera(camera_model::AbstractQCMOSCameraModel;
    scan_mode::AbstractQCMOSScanMode=QCMOSUltraQuietScan(),
    readout_noise_sigma::Union{Nothing,Real}=nothing,
    readout_noise_median::Union{Nothing,Real}=nothing,
    dsnu_sigma::Union{Nothing,Real}=nothing,
    prnu_sigma::Union{Nothing,Real}=nothing,
    linearity_error::Union{Nothing,Real}=nothing,
    conversion_factor_e_per_count::Union{Nothing,Real}=nothing,
    pixel_size_um::Union{Nothing,Real}=nothing,
    peak_qe::Union{Nothing,Real}=nothing,
    full_well::Union{Nothing,Real}=nothing,
    dark_current::Union{Nothing,Real}=nothing,
    frame_size::Union{Nothing,Tuple{Int,Int}}=nothing,
    frame_rate_hz::Union{Nothing,Real}=nothing,
    min_exposure_time::Union{Nothing,Real}=nothing,
    output_model::AbstractCMOSOutputModel=NullCMOSOutputModel(),
    timing_model::Union{Nothing,AbstractFrameTimingModel}=nothing,
    T::Type{<:AbstractFloat}=Float64,
    backend::AbstractArrayBackend=CPUBackend(),
    kwargs...,
)
    sensor = QCMOSSensor(;
        camera_model,
        scan_mode,
        readout_noise_sigma,
        readout_noise_median,
        dsnu_sigma,
        prnu_sigma,
        linearity_error,
        conversion_factor_e_per_count,
        pixel_size_um,
        peak_qe,
        full_well,
        dark_current,
        frame_size,
        frame_rate_hz,
        min_exposure_time,
        output_model,
        timing_model,
        T,
        backend,
    )
    return QCMOSDetector(; sensor, T, backend, kwargs...)
end

ORCAQuestDetector(; kwargs...) = _qcmos_detector_for_camera(ORCAQuest(); kwargs...)
ORCAQuest2Detector(; kwargs...) = _qcmos_detector_for_camera(ORCAQuest2(); kwargs...)
ORCAQuestIQDetector(; kwargs...) = _qcmos_detector_for_camera(ORCAQuestIQ(); kwargs...)
