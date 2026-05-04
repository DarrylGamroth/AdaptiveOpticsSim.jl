abstract type AbstractQCMOSCameraModel end
abstract type AbstractQCMOSScanMode end

struct ORCAQuest <: AbstractQCMOSCameraModel end
struct ORCAQuest2 <: AbstractQCMOSCameraModel end
struct ORCAQuestIQ <: AbstractQCMOSCameraModel end

struct QCMOSStandardScan <: AbstractQCMOSScanMode end
struct QCMOSUltraQuietScan <: AbstractQCMOSScanMode end
struct QCMOSPhotonNumberResolvingScan <: AbstractQCMOSScanMode end
struct QCMOSRawScan <: AbstractQCMOSScanMode end

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

detector_sensor_symbol(::QCMOSSensor) = :qcmos
supports_detector_defect_maps(::QCMOSSensor) = true
supports_shutter_timing(::QCMOSSensor) = true
supports_photon_number_resolving(::QCMOSSensor) = true
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
    converted_output = convert_cmos_output_model(output_model, T, backend)
    validated_output = validate_cmos_output_model(converted_output)
    resolved_timing = something(timing_model, qcmos_default_timing_model(camera_model, scan_mode, T))
    validated_timing = validate_frame_timing_model(convert_frame_timing_model(resolved_timing, T))
    return QCMOSSensor{
        T,
        typeof(camera_model),
        typeof(scan_mode),
        typeof(validated_output),
        typeof(validated_timing),
    }(
        camera_model,
        scan_mode,
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
    return det.state.frame
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
    return cube
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
    bits::Union{Nothing,Int}=16,
    output_type::Union{Nothing,DataType}=UInt16,
    T::Union{Nothing,Type{<:AbstractFloat}}=nothing,
    backend::AbstractArrayBackend=CPUBackend(),
    kwargs...,
)
    detector_T = something(T, sensor === nothing ? Float64 : _qcmos_float_type(sensor))
    resolved_sensor = sensor === nothing ? QCMOSSensor(T=detector_T, backend=backend) : sensor
    resolved_integration_time = something(integration_time, max(resolved_sensor.min_exposure_time, inv(resolved_sensor.frame_rate_hz)))
    return Detector(;
        integration_time=resolved_integration_time,
        noise=_default_qcmos_noise(resolved_sensor, noise),
        qe=something(qe, resolved_sensor.peak_qe),
        dark_current=something(dark_current, resolved_sensor.dark_current),
        full_well=something(full_well, resolved_sensor.full_well),
        bits=bits,
        output_type=output_type,
        sensor=resolved_sensor,
        T=detector_T,
        backend=backend,
        kwargs...,
    )
end

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
