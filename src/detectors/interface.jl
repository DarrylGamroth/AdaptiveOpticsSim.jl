abstract type NoiseModel end
abstract type SensorType end
abstract type FrameSensorType <: SensorType end
abstract type CountingSensorType <: SensorType end

abstract type AbstractFrameDetector <: AbstractDetector end
abstract type AbstractCountingDetector <: AbstractDetector end
abstract type CountingDeadTimeModel end
abstract type AbstractCountingGateModel end
abstract type AbstractCountingCorrelationModel end
abstract type AbstractDetectorResponse end
abstract type AbstractFrameResponse <: AbstractDetectorResponse end
abstract type AbstractChargeCouplingModel end
const FrameResponseModel = AbstractFrameResponse
abstract type BackgroundModel end
abstract type AbstractDetectorDefectModel end
abstract type FrameSamplingMode end
abstract type AbstractFrameTimingModel end
abstract type FrameReadoutCorrectionModel end
abstract type FrameReadoutProducts end
abstract type AbstractFrameNonlinearityModel end
abstract type AbstractPersistenceModel end
abstract type AbstractDetectorThermalModel end
abstract type AbstractDetectorThermalState end
abstract type AbstractTemperatureLaw end
abstract type AbstractQuantumEfficiencyModel end
abstract type AbstractTemporalFrameSource end
abstract type AbstractRollingShutterExposureMode end
abstract type AvalancheFrameSensorType <: FrameSensorType end
abstract type HgCdTeAvalancheArraySensorType <: AvalancheFrameSensorType end
abstract type SPADArraySensorType <: CountingSensorType end
abstract type MKIDArraySensorType <: CountingSensorType end

supports_detector_thermal_model(::AbstractFrameDetector) = false
supports_detector_thermal_model(det::AbstractCountingDetector) = !is_null_thermal_model(thermal_model(det))
supports_counting_noise(::AbstractCountingDetector) = true
supports_dead_time(det::AbstractCountingDetector) = supports_dead_time(counting_dead_time_model(det))
supports_channel_gain_map(det::AbstractCountingDetector) = !isnothing(counting_channel_gain_map(det))
supports_counting_gating(det::AbstractCountingDetector) = !is_null_counting_gate(counting_gate_model(det))
supports_afterpulsing(det::AbstractCountingDetector) = _supports_afterpulsing(counting_correlation_model(det))
supports_channel_crosstalk(det::AbstractCountingDetector) = _supports_channel_crosstalk(counting_correlation_model(det))
supports_paralyzable_dead_time(det::AbstractCountingDetector) = is_paralyzable_dead_time(counting_dead_time_model(det))

detector_sensor_symbol(sensor::SensorType) =
    throw(InvalidConfiguration("missing detector_sensor_symbol overload for $(typeof(sensor))"))

supports_clock_induced_charge(::FrameSensorType) = false
supports_column_readout_noise(::FrameSensorType) = false
supports_avalanche_gain(::FrameSensorType) = false
supports_avalanche_gain(::AvalancheFrameSensorType) = true
supports_sensor_glow(::FrameSensorType) = false
supports_detector_defect_maps(::FrameSensorType) = false
supports_detector_persistence(::FrameSensorType) = false
supports_detector_nonlinearity(::FrameSensorType) = false
supports_shutter_timing(::FrameSensorType) = false
supports_photon_number_resolving(::SensorType) = false
supports_energy_resolving(::SensorType) = false
supports_raw_digital_output(::SensorType) = false
supports_nondestructive_reads(::FrameSensorType) = false
supports_up_the_ramp(::FrameSensorType) = false
supports_reference_read_subtraction(::FrameSensorType) = false
supports_readout_correction(::FrameSensorType) = false
supports_read_cube(::FrameSensorType) = false
supports_detector_response(::SensorType, ::AbstractDetectorResponse) = false
supports_detector_response(::FrameSensorType, ::AbstractFrameResponse) = true
supports_multi_read_readout_products(::FrameSensorType) = false

struct SingleRead <: FrameSamplingMode end

struct AveragedNonDestructiveReads <: FrameSamplingMode
    n_reads::Int
end

struct CorrelatedDoubleSampling <: FrameSamplingMode end

struct FowlerSampling <: FrameSamplingMode
    n_pairs::Int
end

"""
    UpTheRampSampling(n_reads)

Fit a line to `n_reads` evenly spaced nondestructive reads spanning the
integration. The detector output remains in integrated-frame units and is the
fitted slope multiplied by the integration time. At least two reads are
required so that both an intercept and slope can be estimated.
"""
struct UpTheRampSampling <: FrameSamplingMode
    n_reads::Int
end

"""
    SkipperSampling(n_samples)

Repeated nondestructive sampling of a CCD charge packet. The detector keeps a
streaming mean and does not retain a rows-by-columns-by-samples cube.
"""
struct SkipperSampling <: FrameSamplingMode
    n_samples::Int
end

struct FrameWindow
    rows::UnitRange{Int}
    cols::UnitRange{Int}
    function FrameWindow(rows::UnitRange{Int}, cols::UnitRange{Int})
        window = new(rows, cols)
        validate_readout_window(window)
        return window
    end
end

frame_sampling_symbol(sensor::FrameSensorType) = frame_sampling_symbol(multi_read_sampling_mode(sensor))
frame_sampling_reads(sensor::FrameSensorType) = frame_sampling_reads(multi_read_sampling_mode(sensor))
frame_sampling_reference_reads(sensor::FrameSensorType) = frame_sampling_reference_reads(multi_read_sampling_mode(sensor))
frame_sampling_signal_reads(sensor::FrameSensorType) = frame_sampling_signal_reads(multi_read_sampling_mode(sensor))
sampling_read_time(::FrameSensorType, ::Type{T}) where {T<:AbstractFloat} = nothing
sampling_read_time(sensor::FrameSensorType, frame_size::Tuple{Int,Int}, window::Union{Nothing,FrameWindow}, ::Type{T}) where {T<:AbstractFloat} =
    sampling_read_time(sensor, T)
sampling_wallclock_time(::FrameSensorType, integration_time, ::Type{T}) where {T<:AbstractFloat} = T(integration_time)
sampling_wallclock_time(sensor::FrameSensorType, integration_time, frame_size::Tuple{Int,Int},
    window::Union{Nothing,FrameWindow}, ::Type{T}) where {T<:AbstractFloat} =
    sampling_wallclock_time(sensor, integration_time, T)
acquisition_mode_symbol(::FrameSensorType) = :standard
frame_transfer_time(::FrameSensorType, ::Type{T}) where {T<:AbstractFloat} = nothing
steady_state_frame_period(sensor::FrameSensorType, integration_time,
    frame_size::Tuple{Int,Int}, window::Union{Nothing,FrameWindow},
    ::Type{T}) where {T<:AbstractFloat} =
    sampling_wallclock_time(sensor, integration_time, frame_size, window, T)
frame_sampling_symbol(::SingleRead) = :single_read
frame_sampling_symbol(::AveragedNonDestructiveReads) = :averaged_non_destructive_reads
frame_sampling_symbol(::CorrelatedDoubleSampling) = :correlated_double_sampling
frame_sampling_symbol(::FowlerSampling) = :fowler_sampling
frame_sampling_symbol(::UpTheRampSampling) = :up_the_ramp
frame_sampling_symbol(::SkipperSampling) = :skipper

frame_sampling_reads(::SingleRead) = 1
frame_sampling_reads(mode::AveragedNonDestructiveReads) = mode.n_reads
frame_sampling_reads(::CorrelatedDoubleSampling) = 2
frame_sampling_reads(mode::FowlerSampling) = 2 * mode.n_pairs
frame_sampling_reads(mode::UpTheRampSampling) = mode.n_reads
frame_sampling_reads(mode::SkipperSampling) = mode.n_samples

frame_sampling_reference_reads(::SingleRead) = 0
frame_sampling_reference_reads(::AveragedNonDestructiveReads) = 0
frame_sampling_reference_reads(::CorrelatedDoubleSampling) = 1
frame_sampling_reference_reads(mode::FowlerSampling) = mode.n_pairs
frame_sampling_reference_reads(::UpTheRampSampling) = 0
frame_sampling_reference_reads(::SkipperSampling) = 0

frame_sampling_signal_reads(::SingleRead) = 1
frame_sampling_signal_reads(mode::AveragedNonDestructiveReads) = mode.n_reads
frame_sampling_signal_reads(::CorrelatedDoubleSampling) = 1
frame_sampling_signal_reads(mode::FowlerSampling) = mode.n_pairs
frame_sampling_signal_reads(mode::UpTheRampSampling) = mode.n_reads
frame_sampling_signal_reads(mode::SkipperSampling) = mode.n_samples

effective_readout_sigma(::FrameSamplingMode, sigma) = sigma
effective_readout_sigma(mode::AveragedNonDestructiveReads, sigma) = sigma / sqrt(mode.n_reads)
effective_readout_sigma(::CorrelatedDoubleSampling, sigma) = sigma * sqrt(2)
effective_readout_sigma(mode::FowlerSampling, sigma) = sigma * sqrt(2 / mode.n_pairs)
effective_readout_sigma(mode::UpTheRampSampling, sigma) =
    sigma * sqrt(12 * (mode.n_reads - 1) / (mode.n_reads * (mode.n_reads + 1)))
effective_readout_sigma(mode::SkipperSampling, sigma) = sigma / sqrt(mode.n_samples)

validate_frame_sampling_mode(::SingleRead) = SingleRead()

function validate_frame_sampling_mode(mode::AveragedNonDestructiveReads)
    mode.n_reads >= 1 || throw(InvalidConfiguration("AveragedNonDestructiveReads n_reads must be >= 1"))
    return mode
end

validate_frame_sampling_mode(::CorrelatedDoubleSampling) = CorrelatedDoubleSampling()

function validate_frame_sampling_mode(mode::FowlerSampling)
    mode.n_pairs >= 1 || throw(InvalidConfiguration("FowlerSampling n_pairs must be >= 1"))
    return mode
end

function validate_frame_sampling_mode(mode::UpTheRampSampling)
    mode.n_reads >= 2 ||
        throw(InvalidConfiguration("UpTheRampSampling n_reads must be >= 2"))
    return mode
end


function validate_frame_sampling_mode(mode::SkipperSampling)
    mode.n_samples >= 1 ||
        throw(InvalidConfiguration("SkipperSampling n_samples must be >= 1"))
    return mode
end

struct NullFrameReadoutCorrection <: FrameReadoutCorrectionModel end

struct ReferencePixelCommonModeCorrection <: FrameReadoutCorrectionModel
    edge_rows::Int
    edge_cols::Int
    function ReferencePixelCommonModeCorrection(edge_rows::Integer=4, edge_cols::Integer=4)
        edge_rows >= 0 || throw(InvalidConfiguration("ReferencePixelCommonModeCorrection edge_rows must be >= 0"))
        edge_cols >= 0 || throw(InvalidConfiguration("ReferencePixelCommonModeCorrection edge_cols must be >= 0"))
        (edge_rows > 0 || edge_cols > 0) ||
            throw(InvalidConfiguration("ReferencePixelCommonModeCorrection requires at least one reference-pixel edge"))
        return new(Int(edge_rows), Int(edge_cols))
    end
end

struct ReferenceRowCommonModeCorrection <: FrameReadoutCorrectionModel
    edge_cols::Int
    function ReferenceRowCommonModeCorrection(edge_cols::Integer=4)
        edge_cols > 0 || throw(InvalidConfiguration("ReferenceRowCommonModeCorrection edge_cols must be > 0"))
        return new(Int(edge_cols))
    end
end

struct ReferenceColumnCommonModeCorrection <: FrameReadoutCorrectionModel
    edge_rows::Int
    function ReferenceColumnCommonModeCorrection(edge_rows::Integer=4)
        edge_rows > 0 || throw(InvalidConfiguration("ReferenceColumnCommonModeCorrection edge_rows must be > 0"))
        return new(Int(edge_rows))
    end
end

struct ReferenceOutputCommonModeCorrection <: FrameReadoutCorrectionModel
    output_cols::Int
    edge_rows::Int
    edge_cols::Int
    function ReferenceOutputCommonModeCorrection(output_cols::Integer; edge_rows::Integer=4, edge_cols::Integer=4)
        output_cols > 0 || throw(InvalidConfiguration("ReferenceOutputCommonModeCorrection output_cols must be > 0"))
        edge_rows >= 0 || throw(InvalidConfiguration("ReferenceOutputCommonModeCorrection edge_rows must be >= 0"))
        edge_cols >= 0 || throw(InvalidConfiguration("ReferenceOutputCommonModeCorrection edge_cols must be >= 0"))
        (edge_rows > 0 || edge_cols > 0) ||
            throw(InvalidConfiguration("ReferenceOutputCommonModeCorrection requires at least one reference-pixel edge"))
        return new(Int(output_cols), Int(edge_rows), Int(edge_cols))
    end
end

struct CompositeFrameReadoutCorrection{M<:Tuple} <: FrameReadoutCorrectionModel
    stages::M
    function CompositeFrameReadoutCorrection(stages::Tuple{Vararg{FrameReadoutCorrectionModel}})
        isempty(stages) && throw(InvalidConfiguration("CompositeFrameReadoutCorrection requires at least one stage"))
        return new{typeof(stages)}(stages)
    end
end

CompositeFrameReadoutCorrection(stages::Tuple) =
    throw(InvalidConfiguration("CompositeFrameReadoutCorrection stages must be FrameReadoutCorrectionModel values"))
CompositeFrameReadoutCorrection(stages::FrameReadoutCorrectionModel...) = CompositeFrameReadoutCorrection(tuple(stages...))

readout_correction_symbol(::NullFrameReadoutCorrection) = :none
readout_correction_symbol(::ReferencePixelCommonModeCorrection) = :reference_pixel_common_mode
readout_correction_symbol(::ReferenceRowCommonModeCorrection) = :reference_row_common_mode
readout_correction_symbol(::ReferenceColumnCommonModeCorrection) = :reference_column_common_mode
readout_correction_symbol(::ReferenceOutputCommonModeCorrection) = :reference_output_common_mode
readout_correction_symbol(::CompositeFrameReadoutCorrection) = :composite

is_null_readout_correction(::FrameReadoutCorrectionModel) = false
is_null_readout_correction(::NullFrameReadoutCorrection) = true

correction_edge_rows(::FrameReadoutCorrectionModel) = nothing
correction_edge_cols(::FrameReadoutCorrectionModel) = nothing
correction_edge_rows(model::ReferencePixelCommonModeCorrection) = model.edge_rows
correction_edge_cols(model::ReferencePixelCommonModeCorrection) = model.edge_cols
correction_edge_rows(model::ReferenceColumnCommonModeCorrection) = model.edge_rows
correction_edge_cols(model::ReferenceRowCommonModeCorrection) = model.edge_cols
correction_edge_rows(model::ReferenceOutputCommonModeCorrection) = model.edge_rows
correction_edge_cols(model::ReferenceOutputCommonModeCorrection) = model.edge_cols

correction_group_rows(::FrameReadoutCorrectionModel) = nothing
correction_group_cols(::FrameReadoutCorrectionModel) = nothing
correction_group_cols(model::ReferenceOutputCommonModeCorrection) = model.output_cols

correction_stage_count(::FrameReadoutCorrectionModel) = 1
correction_stage_count(model::CompositeFrameReadoutCorrection) = length(model.stages)

supports_batched_readout_correction(::FrameReadoutCorrectionModel) = false
supports_batched_readout_correction(::NullFrameReadoutCorrection) = true
supports_batched_readout_correction(::ReferencePixelCommonModeCorrection) = true
supports_batched_readout_correction(::ReferenceRowCommonModeCorrection) = true
supports_batched_readout_correction(::ReferenceColumnCommonModeCorrection) = true
supports_batched_readout_correction(::ReferenceOutputCommonModeCorrection) = true
supports_batched_readout_correction(model::CompositeFrameReadoutCorrection) =
    all(supports_batched_readout_correction, model.stages)

struct NoFrameReadoutProducts <: FrameReadoutProducts end

struct SkipperReadoutProducts{A<:AbstractMatrix} <: FrameReadoutProducts
    mean_frame::A
    sample_count::Int
end

# Sampled detector parameters are copied into run-owned storage at every
# public construction boundary. The token is reserved for package code that
# has already allocated and populated a fresh owned buffer.
struct OwnedDetectorParameterToken end
const OWNED_DETECTOR_PARAMETER = OwnedDetectorParameterToken()

struct NullDetectorDefectModel <: AbstractDetectorDefectModel end

struct PixelResponseNonuniformity{T<:AbstractFloat,A<:AbstractMatrix{T}} <: AbstractDetectorDefectModel
    gain_map::A
    function PixelResponseNonuniformity{T,A}(
        ::OwnedDetectorParameterToken, gain_map::A) where {
        T<:AbstractFloat,A<:AbstractMatrix{T}}
        return new{T,A}(gain_map)
    end

    function PixelResponseNonuniformity{T,A}(gain_map::A) where {
        T<:AbstractFloat,A<:AbstractMatrix{T}}
        owned = copy(gain_map)
        return validate_detector_defect_model(
            PixelResponseNonuniformity{T,typeof(owned)}(
                OWNED_DETECTOR_PARAMETER, owned))
    end
end

struct DarkSignalNonuniformity{T<:AbstractFloat,A<:AbstractMatrix{T}} <: AbstractDetectorDefectModel
    dark_map::A
    function DarkSignalNonuniformity{T,A}(
        ::OwnedDetectorParameterToken, dark_map::A) where {
        T<:AbstractFloat,A<:AbstractMatrix{T}}
        return new{T,A}(dark_map)
    end

    function DarkSignalNonuniformity{T,A}(dark_map::A) where {
        T<:AbstractFloat,A<:AbstractMatrix{T}}
        owned = copy(dark_map)
        return validate_detector_defect_model(
            DarkSignalNonuniformity{T,typeof(owned)}(
                OWNED_DETECTOR_PARAMETER, owned))
    end
end

struct BadPixelMask{T<:AbstractFloat,A<:AbstractMatrix{Bool}} <: AbstractDetectorDefectModel
    mask::A
    throughput::T
    function BadPixelMask{T,A}(::OwnedDetectorParameterToken, mask::A,
        throughput::T) where {T<:AbstractFloat,A<:AbstractMatrix{Bool}}
        return new{T,A}(mask, throughput)
    end

    function BadPixelMask{T,A}(mask::A, throughput::T) where {
        T<:AbstractFloat,A<:AbstractMatrix{Bool}}
        owned = copy(mask)
        return validate_detector_defect_model(
            BadPixelMask{T,typeof(owned)}(OWNED_DETECTOR_PARAMETER,
                owned, throughput))
    end
end

struct CompositeDetectorDefectModel{M<:Tuple} <: AbstractDetectorDefectModel
    stages::M
    function CompositeDetectorDefectModel(::OwnedDetectorParameterToken,
        stages::Tuple{Vararg{AbstractDetectorDefectModel}})
        isempty(stages) && throw(InvalidConfiguration("CompositeDetectorDefectModel requires at least one stage"))
        return new{typeof(stages)}(stages)
    end
end

@inline owned_detector_defect_model(model::AbstractDetectorDefectModel) = model
@inline owned_detector_defect_model(::NullDetectorDefectModel) =
    NullDetectorDefectModel()
@inline owned_detector_defect_model(model::PixelResponseNonuniformity) =
    PixelResponseNonuniformity(model.gain_map)
@inline owned_detector_defect_model(model::DarkSignalNonuniformity) =
    DarkSignalNonuniformity(model.dark_map)
@inline owned_detector_defect_model(model::BadPixelMask) =
    BadPixelMask(model.mask; throughput=model.throughput)

function owned_detector_defect_model(model::CompositeDetectorDefectModel)
    return CompositeDetectorDefectModel(model.stages)
end

function CompositeDetectorDefectModel(
    stages::Tuple{Vararg{AbstractDetectorDefectModel}})
    owned_stages = map(owned_detector_defect_model, stages)
    return CompositeDetectorDefectModel(OWNED_DETECTOR_PARAMETER,
        owned_stages)
end

CompositeDetectorDefectModel(stages::Tuple) =
    throw(InvalidConfiguration("CompositeDetectorDefectModel stages must be AbstractDetectorDefectModel values"))
CompositeDetectorDefectModel(stages::AbstractDetectorDefectModel...) = CompositeDetectorDefectModel(tuple(stages...))

detector_defect_symbol(::NullDetectorDefectModel) = :none
detector_defect_symbol(::PixelResponseNonuniformity) = :prnu
detector_defect_symbol(::DarkSignalNonuniformity) = :dsnu
detector_defect_symbol(::BadPixelMask) = :bad_pixel_mask
detector_defect_symbol(::CompositeDetectorDefectModel) = :composite
default_detector_defect_model(::FrameSensorType; T::Type{<:AbstractFloat}=Float64, backend::AbstractArrayBackend=CPUBackend()) = NullDetectorDefectModel()

struct GlobalShutter <: AbstractFrameTimingModel end

struct RollingExposure <: AbstractRollingShutterExposureMode end
struct GlobalResetExposure <: AbstractRollingShutterExposureMode end

struct RollingShutter{T<:AbstractFloat,M<:AbstractRollingShutterExposureMode} <: AbstractFrameTimingModel
    line_time::T
    row_group_size::Int
    exposure_mode::M
end

RollingShutter(line_time::Real; row_group_size::Integer=1,
    exposure_mode::AbstractRollingShutterExposureMode=RollingExposure()) =
    RollingShutter{Float64,typeof(exposure_mode)}(float(line_time), Int(row_group_size), exposure_mode)
RollingShutter{T}(line_time::Real, row_group_size::Integer=1;
    exposure_mode::AbstractRollingShutterExposureMode=RollingExposure()) where {T<:AbstractFloat} =
    RollingShutter{T,typeof(exposure_mode)}(T(line_time), Int(row_group_size), exposure_mode)

timing_model_symbol(::GlobalShutter) = :global_shutter
timing_model_symbol(::RollingShutter) = :rolling_shutter
is_global_shutter(::AbstractFrameTimingModel) = false
is_global_shutter(::GlobalShutter) = true
default_frame_timing_model(::FrameSensorType; T::Type{<:AbstractFloat}=Float64) = GlobalShutter()

struct FunctionFrameSource{F} <: AbstractTemporalFrameSource
    f::F
end

struct InPlaceFrameSource{F} <: AbstractTemporalFrameSource
    f::F
    frame_size::Tuple{Int,Int}
end

struct FunctionExposureFrameSource{F} <: AbstractTemporalFrameSource
    f::F
end

struct InPlaceExposureFrameSource{F} <: AbstractTemporalFrameSource
    f::F
    frame_size::Tuple{Int,Int}
end

function sample_frame!(dest::AbstractMatrix, source::FunctionFrameSource, time)
    frame = source.f(time)
    size(frame) == size(dest) || throw(DimensionMismatchError("temporal frame source returned an unexpected frame size"))
    copyto!(dest, frame)
    return dest
end

function sample_frame!(dest::AbstractMatrix, source::InPlaceFrameSource, time)
    size(dest) == source.frame_size || throw(DimensionMismatchError("temporal frame source destination has an unexpected frame size"))
    source.f(dest, time)
    return dest
end

sample_exposure_frame!(dest::AbstractMatrix, source::AbstractTemporalFrameSource, start_time, exposure_time) =
    sample_frame!(dest, source, start_time)

function sample_exposure_frame!(dest::AbstractMatrix, source::FunctionExposureFrameSource, start_time, exposure_time)
    frame = source.f(start_time, exposure_time)
    size(frame) == size(dest) || throw(DimensionMismatchError("temporal exposure frame source returned an unexpected frame size"))
    copyto!(dest, frame)
    return dest
end

function sample_exposure_frame!(dest::AbstractMatrix, source::InPlaceExposureFrameSource, start_time, exposure_time)
    size(dest) == source.frame_size ||
        throw(DimensionMismatchError("temporal exposure frame source destination has an unexpected frame size"))
    source.f(dest, start_time, exposure_time)
    return dest
end

struct NullFrameNonlinearity <: AbstractFrameNonlinearityModel end

struct SaturatingFrameNonlinearity{T<:AbstractFloat} <: AbstractFrameNonlinearityModel
    coefficient::T
end

SaturatingFrameNonlinearity(coefficient::Real) = SaturatingFrameNonlinearity{Float64}(float(coefficient))

nonlinearity_symbol(::NullFrameNonlinearity) = :none
nonlinearity_symbol(::SaturatingFrameNonlinearity) = :saturating
is_null_frame_nonlinearity(::AbstractFrameNonlinearityModel) = false
is_null_frame_nonlinearity(::NullFrameNonlinearity) = true
default_frame_nonlinearity_model(::FrameSensorType; T::Type{<:AbstractFloat}=Float64) = NullFrameNonlinearity()

struct NullPersistence <: AbstractPersistenceModel end

struct NullDetectorThermalModel <: AbstractDetectorThermalModel end

struct NoThermalState <: AbstractDetectorThermalState end

mutable struct DetectorThermalState{T<:AbstractFloat} <: AbstractDetectorThermalState
    temperature_K::T
end

struct NullTemperatureLaw <: AbstractTemperatureLaw end

struct ArrheniusRateLaw{T<:AbstractFloat} <: AbstractTemperatureLaw
    reference_temperature_K::T
    activation_temperature_K::T
end

ArrheniusRateLaw(reference_temperature_K::Real, activation_temperature_K::Real) =
    ArrheniusRateLaw{Float64}(float(reference_temperature_K), float(activation_temperature_K))

struct LinearTemperatureLaw{T<:AbstractFloat} <: AbstractTemperatureLaw
    reference_temperature_K::T
    slope_per_K::T
end

LinearTemperatureLaw(reference_temperature_K::Real, slope_per_K::Real) =
    LinearTemperatureLaw{Float64}(float(reference_temperature_K), float(slope_per_K))

struct ExponentialTemperatureLaw{T<:AbstractFloat} <: AbstractTemperatureLaw
    reference_temperature_K::T
    exponent_per_K::T
end

ExponentialTemperatureLaw(reference_temperature_K::Real, exponent_per_K::Real) =
    ExponentialTemperatureLaw{Float64}(float(reference_temperature_K), float(exponent_per_K))

struct FixedTemperature{
    T<:AbstractFloat,
    DC<:AbstractTemperatureLaw,
    G<:AbstractTemperatureLaw,
    DCR<:AbstractTemperatureLaw,
    CIC<:AbstractTemperatureLaw,
} <: AbstractDetectorThermalModel
    temperature_K::T
    dark_current_law::DC
    glow_rate_law::G
    dark_count_law::DCR
    cic_rate_law::CIC
end

function FixedTemperature(; temperature_K::Real,
    dark_current_law::AbstractTemperatureLaw=NullTemperatureLaw(),
    glow_rate_law::AbstractTemperatureLaw=NullTemperatureLaw(),
    dark_count_law::AbstractTemperatureLaw=NullTemperatureLaw(),
    cic_rate_law::AbstractTemperatureLaw=NullTemperatureLaw(),
    T::Type{<:AbstractFloat}=Float64)
    temperature_K > 0 || throw(InvalidConfiguration("FixedTemperature temperature_K must be > 0"))
    return FixedTemperature{
        T,
        typeof(dark_current_law),
        typeof(glow_rate_law),
        typeof(dark_count_law),
        typeof(cic_rate_law),
    }(
        T(temperature_K),
        dark_current_law,
        glow_rate_law,
        dark_count_law,
        cic_rate_law,
    )
end

struct FirstOrderThermalModel{
    T<:AbstractFloat,
    DC<:AbstractTemperatureLaw,
    G<:AbstractTemperatureLaw,
    DCR<:AbstractTemperatureLaw,
    CIC<:AbstractTemperatureLaw,
} <: AbstractDetectorThermalModel
    ambient_temperature_K::T
    setpoint_temperature_K::T
    initial_temperature_K::T
    time_constant_s::T
    min_temperature_K::T
    max_temperature_K::T
    dark_current_law::DC
    glow_rate_law::G
    dark_count_law::DCR
    cic_rate_law::CIC
end

function FirstOrderThermalModel(; ambient_temperature_K::Real,
    setpoint_temperature_K::Real,
    time_constant_s::Real,
    initial_temperature_K::Real=setpoint_temperature_K,
    min_temperature_K::Real=min(ambient_temperature_K, setpoint_temperature_K, initial_temperature_K),
    max_temperature_K::Real=max(ambient_temperature_K, setpoint_temperature_K, initial_temperature_K),
    dark_current_law::AbstractTemperatureLaw=NullTemperatureLaw(),
    glow_rate_law::AbstractTemperatureLaw=NullTemperatureLaw(),
    dark_count_law::AbstractTemperatureLaw=NullTemperatureLaw(),
    cic_rate_law::AbstractTemperatureLaw=NullTemperatureLaw(),
    T::Type{<:AbstractFloat}=Float64)
    ambient_temperature_K > 0 || throw(InvalidConfiguration("FirstOrderThermalModel ambient_temperature_K must be > 0"))
    setpoint_temperature_K > 0 || throw(InvalidConfiguration("FirstOrderThermalModel setpoint_temperature_K must be > 0"))
    initial_temperature_K > 0 || throw(InvalidConfiguration("FirstOrderThermalModel initial_temperature_K must be > 0"))
    time_constant_s > 0 || throw(InvalidConfiguration("FirstOrderThermalModel time_constant_s must be > 0"))
    min_temperature_K > 0 || throw(InvalidConfiguration("FirstOrderThermalModel min_temperature_K must be > 0"))
    max_temperature_K > 0 || throw(InvalidConfiguration("FirstOrderThermalModel max_temperature_K must be > 0"))
    min_temperature_K <= max_temperature_K ||
        throw(InvalidConfiguration("FirstOrderThermalModel min_temperature_K must be <= max_temperature_K"))
    min_temperature_K <= initial_temperature_K <= max_temperature_K ||
        throw(InvalidConfiguration("FirstOrderThermalModel initial_temperature_K must lie within [min_temperature_K, max_temperature_K]"))
    min_temperature_K <= setpoint_temperature_K <= max_temperature_K ||
        throw(InvalidConfiguration("FirstOrderThermalModel setpoint_temperature_K must lie within [min_temperature_K, max_temperature_K]"))
    min_temperature_K <= ambient_temperature_K <= max_temperature_K ||
        throw(InvalidConfiguration("FirstOrderThermalModel ambient_temperature_K must lie within [min_temperature_K, max_temperature_K]"))
    return FirstOrderThermalModel{
        T,
        typeof(dark_current_law),
        typeof(glow_rate_law),
        typeof(dark_count_law),
        typeof(cic_rate_law),
    }(
        T(ambient_temperature_K),
        T(setpoint_temperature_K),
        T(initial_temperature_K),
        T(time_constant_s),
        T(min_temperature_K),
        T(max_temperature_K),
        dark_current_law,
        glow_rate_law,
        dark_count_law,
        cic_rate_law,
    )
end

struct ExponentialPersistence{T<:AbstractFloat} <: AbstractPersistenceModel
    coupling::T
    decay::T
end

ExponentialPersistence(coupling::Real, decay::Real) = ExponentialPersistence{Float64}(float(coupling), float(decay))

persistence_symbol(::NullPersistence) = :none
persistence_symbol(::ExponentialPersistence) = :exponential
is_null_persistence(::AbstractPersistenceModel) = false
is_null_persistence(::NullPersistence) = true

thermal_model_symbol(::NullDetectorThermalModel) = :none
thermal_model_symbol(::FixedTemperature) = :fixed_temperature
thermal_model_symbol(::FirstOrderThermalModel) = :first_order
is_null_thermal_model(::AbstractDetectorThermalModel) = false
is_null_thermal_model(::NullDetectorThermalModel) = true

temperature_law_symbol(::NullTemperatureLaw) = :none
temperature_law_symbol(::ArrheniusRateLaw) = :arrhenius
temperature_law_symbol(::LinearTemperatureLaw) = :linear
temperature_law_symbol(::ExponentialTemperatureLaw) = :exponential
is_null_temperature_law(::AbstractTemperatureLaw) = false
is_null_temperature_law(::NullTemperatureLaw) = true

ambient_temperature_K(::AbstractDetectorThermalModel, ::Type{T}) where {T<:AbstractFloat} = nothing
cooling_setpoint_K(::AbstractDetectorThermalModel, ::Type{T}) where {T<:AbstractFloat} = nothing
thermal_time_constant_s(::AbstractDetectorThermalModel, ::Type{T}) where {T<:AbstractFloat} = nothing
cooling_setpoint_K(model::FixedTemperature, ::Type{T}) where {T<:AbstractFloat} = T(model.temperature_K)
ambient_temperature_K(model::FirstOrderThermalModel, ::Type{T}) where {T<:AbstractFloat} = T(model.ambient_temperature_K)
cooling_setpoint_K(model::FirstOrderThermalModel, ::Type{T}) where {T<:AbstractFloat} = T(model.setpoint_temperature_K)
thermal_time_constant_s(model::FirstOrderThermalModel, ::Type{T}) where {T<:AbstractFloat} = T(model.time_constant_s)
detector_temperature_K(::NullDetectorThermalModel, ::AbstractDetectorThermalState, ::Type{T}) where {T<:AbstractFloat} = nothing
detector_temperature_K(model::FixedTemperature, ::AbstractDetectorThermalState, ::Type{T}) where {T<:AbstractFloat} = T(model.temperature_K)
detector_temperature_K(::NullDetectorThermalModel, ::DetectorThermalState, ::Type{T}) where {T<:AbstractFloat} = nothing
detector_temperature_K(model::FixedTemperature, ::DetectorThermalState, ::Type{T}) where {T<:AbstractFloat} = T(model.temperature_K)
detector_temperature_K(model::AbstractDetectorThermalModel, state::DetectorThermalState, ::Type{T}) where {T<:AbstractFloat} = T(state.temperature_K)

advance_thermal!(::NullDetectorThermalModel, ::AbstractDetectorThermalState, dt) = nothing
advance_thermal!(::FixedTemperature, ::AbstractDetectorThermalState, dt) = nothing
function advance_thermal!(model::FirstOrderThermalModel, state::DetectorThermalState, dt)
    dt >= 0 || throw(InvalidConfiguration("thermal evolution step dt must be >= 0"))
    iszero(dt) && return nothing
    equilibrium = clamp(model.setpoint_temperature_K, model.min_temperature_K, model.max_temperature_K)
    state.temperature_K = clamp(equilibrium + (state.temperature_K - equilibrium) * exp(-dt / model.time_constant_s),
        model.min_temperature_K, model.max_temperature_K)
    return nothing
end

evaluate_temperature_law(::NullTemperatureLaw, base_value, detector_temperature) = base_value

function evaluate_temperature_law(law::ArrheniusRateLaw{T}, base_value, detector_temperature) where {T<:AbstractFloat}
    isnothing(detector_temperature) && return base_value
    detector_temperature > zero(detector_temperature) ||
        throw(InvalidConfiguration("detector temperature must be > 0 K"))
    scale = exp(law.activation_temperature_K * (inv(law.reference_temperature_K) - inv(detector_temperature)))
    return base_value * scale
end

function evaluate_temperature_law(law::LinearTemperatureLaw, base_value, detector_temperature)
    isnothing(detector_temperature) && return base_value
    scaled = base_value * (one(base_value) + law.slope_per_K * (detector_temperature - law.reference_temperature_K))
    return max(zero(base_value), scaled)
end

function evaluate_temperature_law(law::ExponentialTemperatureLaw, base_value, detector_temperature)
    isnothing(detector_temperature) && return base_value
    return base_value * exp(law.exponent_per_K * (detector_temperature - law.reference_temperature_K))
end

default_thermal_model(::SensorType; T::Type{<:AbstractFloat}=Float64) = NullDetectorThermalModel()
supports_dynamic_thermal_state(::AbstractDetectorThermalModel) = false
supports_dynamic_thermal_state(::FirstOrderThermalModel) = true

convert_temperature_law(::NullTemperatureLaw, ::Type{T}) where {T<:AbstractFloat} = NullTemperatureLaw()
convert_temperature_law(model::ArrheniusRateLaw, ::Type{T}) where {T<:AbstractFloat} =
    ArrheniusRateLaw{T}(T(model.reference_temperature_K), T(model.activation_temperature_K))
convert_temperature_law(model::LinearTemperatureLaw, ::Type{T}) where {T<:AbstractFloat} =
    LinearTemperatureLaw{T}(T(model.reference_temperature_K), T(model.slope_per_K))
convert_temperature_law(model::ExponentialTemperatureLaw, ::Type{T}) where {T<:AbstractFloat} =
    ExponentialTemperatureLaw{T}(T(model.reference_temperature_K), T(model.exponent_per_K))

validate_temperature_law(::NullTemperatureLaw) = NullTemperatureLaw()

function validate_temperature_law(model::ArrheniusRateLaw)
    model.reference_temperature_K > zero(model.reference_temperature_K) ||
        throw(InvalidConfiguration("ArrheniusRateLaw reference_temperature_K must be > 0"))
    return model
end

function validate_temperature_law(model::LinearTemperatureLaw)
    model.reference_temperature_K > zero(model.reference_temperature_K) ||
        throw(InvalidConfiguration("LinearTemperatureLaw reference_temperature_K must be > 0"))
    return model
end

function validate_temperature_law(model::ExponentialTemperatureLaw)
    model.reference_temperature_K > zero(model.reference_temperature_K) ||
        throw(InvalidConfiguration("ExponentialTemperatureLaw reference_temperature_K must be > 0"))
    return model
end

convert_thermal_model(::NullDetectorThermalModel, ::Type{T}) where {T<:AbstractFloat} = NullDetectorThermalModel()

function convert_thermal_model(model::FixedTemperature, ::Type{T}) where {T<:AbstractFloat}
    return FixedTemperature(
        temperature_K=T(model.temperature_K),
        dark_current_law=validate_temperature_law(convert_temperature_law(model.dark_current_law, T)),
        glow_rate_law=validate_temperature_law(convert_temperature_law(model.glow_rate_law, T)),
        dark_count_law=validate_temperature_law(convert_temperature_law(model.dark_count_law, T)),
        cic_rate_law=validate_temperature_law(convert_temperature_law(model.cic_rate_law, T)),
        T=T,
    )
end

function convert_thermal_model(model::FirstOrderThermalModel, ::Type{T}) where {T<:AbstractFloat}
    return FirstOrderThermalModel(
        ambient_temperature_K=T(model.ambient_temperature_K),
        setpoint_temperature_K=T(model.setpoint_temperature_K),
        initial_temperature_K=T(model.initial_temperature_K),
        time_constant_s=T(model.time_constant_s),
        min_temperature_K=T(model.min_temperature_K),
        max_temperature_K=T(model.max_temperature_K),
        dark_current_law=validate_temperature_law(convert_temperature_law(model.dark_current_law, T)),
        glow_rate_law=validate_temperature_law(convert_temperature_law(model.glow_rate_law, T)),
        dark_count_law=validate_temperature_law(convert_temperature_law(model.dark_count_law, T)),
        cic_rate_law=validate_temperature_law(convert_temperature_law(model.cic_rate_law, T)),
        T=T,
    )
end

validate_thermal_model(::NullDetectorThermalModel) = NullDetectorThermalModel()

function validate_thermal_model(model::FixedTemperature)
    model.temperature_K > zero(model.temperature_K) ||
        throw(InvalidConfiguration("FixedTemperature temperature_K must be > 0"))
    validate_temperature_law(model.dark_current_law)
    validate_temperature_law(model.glow_rate_law)
    validate_temperature_law(model.dark_count_law)
    validate_temperature_law(model.cic_rate_law)
    return model
end

function validate_thermal_model(model::FirstOrderThermalModel)
    model.ambient_temperature_K > zero(model.ambient_temperature_K) ||
        throw(InvalidConfiguration("FirstOrderThermalModel ambient_temperature_K must be > 0"))
    model.setpoint_temperature_K > zero(model.setpoint_temperature_K) ||
        throw(InvalidConfiguration("FirstOrderThermalModel setpoint_temperature_K must be > 0"))
    model.initial_temperature_K > zero(model.initial_temperature_K) ||
        throw(InvalidConfiguration("FirstOrderThermalModel initial_temperature_K must be > 0"))
    model.time_constant_s > zero(model.time_constant_s) ||
        throw(InvalidConfiguration("FirstOrderThermalModel time_constant_s must be > 0"))
    model.min_temperature_K > zero(model.min_temperature_K) ||
        throw(InvalidConfiguration("FirstOrderThermalModel min_temperature_K must be > 0"))
    model.max_temperature_K > zero(model.max_temperature_K) ||
        throw(InvalidConfiguration("FirstOrderThermalModel max_temperature_K must be > 0"))
    model.min_temperature_K <= model.initial_temperature_K <= model.max_temperature_K ||
        throw(InvalidConfiguration("FirstOrderThermalModel initial_temperature_K must lie within [min_temperature_K, max_temperature_K]"))
    model.min_temperature_K <= model.setpoint_temperature_K <= model.max_temperature_K ||
        throw(InvalidConfiguration("FirstOrderThermalModel setpoint_temperature_K must lie within [min_temperature_K, max_temperature_K]"))
    model.min_temperature_K <= model.ambient_temperature_K <= model.max_temperature_K ||
        throw(InvalidConfiguration("FirstOrderThermalModel ambient_temperature_K must lie within [min_temperature_K, max_temperature_K]"))
    validate_temperature_law(model.dark_current_law)
    validate_temperature_law(model.glow_rate_law)
    validate_temperature_law(model.dark_count_law)
    validate_temperature_law(model.cic_rate_law)
    return model
end

resolve_thermal_model(sensor::SensorType, ::Nothing; T::Type{<:AbstractFloat}) =
    default_thermal_model(sensor; T=T)
resolve_thermal_model(sensor::SensorType, model::AbstractDetectorThermalModel; T::Type{<:AbstractFloat}) = model

thermal_state_from_model(::NullDetectorThermalModel, ::Type{T}) where {T<:AbstractFloat} = NoThermalState()
thermal_state_from_model(::FixedTemperature, ::Type{T}) where {T<:AbstractFloat} = NoThermalState()
thermal_state_from_model(model::FirstOrderThermalModel, ::Type{T}) where {T<:AbstractFloat} =
    DetectorThermalState{T}(T(model.initial_temperature_K))

struct NullCountingGate <: AbstractCountingGateModel end

struct DutyCycleGate{T<:AbstractFloat} <: AbstractCountingGateModel
    duty_cycle::T
end

DutyCycleGate(duty_cycle::Real) = DutyCycleGate{Float64}(float(duty_cycle))

struct NullCountingCorrelation <: AbstractCountingCorrelationModel end

struct AfterpulsingModel{T<:AbstractFloat} <: AbstractCountingCorrelationModel
    probability::T
end

AfterpulsingModel(probability::Real) = AfterpulsingModel{Float64}(float(probability))

struct ChannelCrosstalkModel{T<:AbstractFloat} <: AbstractCountingCorrelationModel
    coupling::T
end

ChannelCrosstalkModel(coupling::Real) = ChannelCrosstalkModel{Float64}(float(coupling))

struct CompositeCountingCorrelation{M<:Tuple} <: AbstractCountingCorrelationModel
    stages::M
    function CompositeCountingCorrelation(stages::Tuple{Vararg{AbstractCountingCorrelationModel}})
        isempty(stages) && throw(InvalidConfiguration("CompositeCountingCorrelation requires at least one stage"))
        return new{typeof(stages)}(stages)
    end
end

CompositeCountingCorrelation(stages::Tuple) =
    throw(InvalidConfiguration("CompositeCountingCorrelation stages must be AbstractCountingCorrelationModel values"))
CompositeCountingCorrelation(stages::AbstractCountingCorrelationModel...) = CompositeCountingCorrelation(tuple(stages...))

counting_gate_symbol(::NullCountingGate) = :none
counting_gate_symbol(::DutyCycleGate) = :duty_cycle

counting_correlation_symbol(::NullCountingCorrelation) = :none
counting_correlation_symbol(::AfterpulsingModel) = :afterpulsing
counting_correlation_symbol(::ChannelCrosstalkModel) = :channel_crosstalk
counting_correlation_symbol(::CompositeCountingCorrelation) = :composite

function PixelResponseNonuniformity(gain_map::AbstractMatrix; T::Type{<:AbstractFloat}=Float64, backend::AbstractArrayBackend=CPUBackend())
    backend = _resolve_array_backend(backend)
    backend_map = _to_backend_matrix(T.(gain_map), backend)
    return validate_detector_defect_model(
        PixelResponseNonuniformity{T,typeof(backend_map)}(
            OWNED_DETECTOR_PARAMETER, backend_map))
end

function DarkSignalNonuniformity(dark_map::AbstractMatrix; T::Type{<:AbstractFloat}=Float64, backend::AbstractArrayBackend=CPUBackend())
    backend = _resolve_array_backend(backend)
    backend_map = _to_backend_matrix(T.(dark_map), backend)
    return validate_detector_defect_model(
        DarkSignalNonuniformity{T,typeof(backend_map)}(
            OWNED_DETECTOR_PARAMETER, backend_map))
end

function BadPixelMask(mask::AbstractMatrix{Bool}; throughput::Real=0.0, T::Type{<:AbstractFloat}=Float64, backend::AbstractArrayBackend=CPUBackend())
    backend = _resolve_array_backend(backend)
    backend_mask = _to_backend_bool_matrix(mask, backend)
    return validate_detector_defect_model(
        BadPixelMask{T,typeof(backend_mask)}(OWNED_DETECTOR_PARAMETER,
            backend_mask, T(throughput)))
end

struct SampledFrameReadoutProducts{A<:AbstractMatrix,C} <: FrameReadoutProducts
    reference_frame::Union{Nothing,A}
    signal_frame::A
    read_cube::Union{Nothing,C}
end

struct MultiReadFrameReadoutProducts{A<:AbstractMatrix,C,V} <: FrameReadoutProducts
    reference_frame::Union{Nothing,A}
    signal_frame::A
    combined_frame::A
    reference_cube::Union{Nothing,C}
    signal_cube::Union{Nothing,C}
    read_cube::Union{Nothing,C}
    read_times::Union{Nothing,V}

    function MultiReadFrameReadoutProducts{A,C,V}(reference_frame::Union{Nothing,A}, signal_frame::A, combined_frame::A,
        reference_cube::Union{Nothing,C}, signal_cube::Union{Nothing,C}, read_cube::Union{Nothing,C},
        read_times::Union{Nothing,V}) where {A<:AbstractMatrix,C,V}
        return new{A,C,V}(reference_frame, signal_frame, combined_frame,
            reference_cube, signal_cube, read_cube, read_times)
    end
end

"""
    UpTheRampReadoutProducts

Preallocated products from an up-the-ramp fit. `integrated_frame` is the slope
multiplied by the integration time and therefore matches ordinary detector
output units. The `workspace_*` arrays are detector-owned full-frame storage;
the public products may be windowed views copied into reusable arrays.
"""
struct UpTheRampReadoutProducts{A<:AbstractMatrix,C<:AbstractArray,V<:AbstractVector} <:
    FrameReadoutProducts
    slope_frame::A
    intercept_frame::A
    integrated_frame::A
    read_cube::C
    read_times::V
    workspace_slope::A
    workspace_intercept::A
    workspace_integrated::A
    workspace_cube::C
end

@inline function _multi_read_cube_param(reference_cube, signal_cube, read_cube)
    if !isnothing(reference_cube)
        return typeof(reference_cube)
    elseif !isnothing(signal_cube)
        return typeof(signal_cube)
    else
        return typeof(read_cube)
    end
end

function MultiReadFrameReadoutProducts(reference_frame::Union{Nothing,A}, signal_frame::A, combined_frame::A,
    reference_cube, signal_cube, read_cube, read_times) where {A<:AbstractMatrix}
    C = _multi_read_cube_param(reference_cube, signal_cube, read_cube)
    V = typeof(read_times)
    return MultiReadFrameReadoutProducts{A,C,V}(reference_frame, signal_frame, combined_frame,
        reference_cube, signal_cube, read_cube, read_times)
end

const HgCdTeReadoutProducts = MultiReadFrameReadoutProducts

function SampledFrameReadoutProducts(reference_frame::Union{Nothing,A}, signal_frame::A,
    read_cube::Nothing) where {A<:AbstractMatrix}
    return SampledFrameReadoutProducts{A,Nothing}(reference_frame, signal_frame, read_cube)
end

function SampledFrameReadoutProducts(reference_frame::Union{Nothing,A}, signal_frame::A,
    read_cube::C) where {A<:AbstractMatrix,C<:AbstractArray}
    return SampledFrameReadoutProducts{A,C}(reference_frame, signal_frame, read_cube)
end

struct NoiseNone <: NoiseModel end
struct NoisePhoton <: NoiseModel end

struct NoiseReadout{T<:AbstractFloat} <: NoiseModel
    sigma::T
end

struct NoisePhotonReadout{T<:AbstractFloat} <: NoiseModel
    sigma::T
end

detector_noise_symbol(::NoiseNone) = :none
detector_noise_symbol(::NoisePhoton) = :photon
detector_noise_symbol(::NoiseReadout) = :readout
detector_noise_symbol(::NoisePhotonReadout) = :photon_readout

struct NullFrameResponse <: AbstractFrameResponse end

struct GaussianPixelResponse{T<:AbstractFloat,V<:AbstractVector{T}} <: AbstractFrameResponse
    response_width_px::T
    kernel::V
    function GaussianPixelResponse{T,V}(::OwnedDetectorParameterToken,
        response_width_px::T,
        kernel::V) where {T<:AbstractFloat,V<:AbstractVector{T}}
        isfinite(response_width_px) && response_width_px > zero(T) ||
            throw(InvalidConfiguration(
                "GaussianPixelResponse response_width_px must be finite and > 0"))
        isempty(kernel) && throw(InvalidConfiguration(
            "GaussianPixelResponse kernel must not be empty"))
        isodd(length(kernel)) || throw(InvalidConfiguration(
            "GaussianPixelResponse kernel length must be odd"))
        _validate_physical_response_kernel(kernel, "GaussianPixelResponse")
        return new{T,V}(response_width_px, kernel)
    end

    function GaussianPixelResponse{T,V}(response_width_px::T,
        kernel::V) where {T<:AbstractFloat,V<:AbstractVector{T}}
        owned = copy(kernel)
        return GaussianPixelResponse{T,typeof(owned)}(
            OWNED_DETECTOR_PARAMETER, response_width_px, owned)
    end
end

struct SampledFrameResponse{T<:AbstractFloat,A<:AbstractMatrix{T}} <: AbstractFrameResponse
    kernel::A
    function SampledFrameResponse{T,A}(::OwnedDetectorParameterToken,
        kernel::A) where {T<:AbstractFloat,A<:AbstractMatrix{T}}
        all(size(kernel) .> 0) || throw(InvalidConfiguration("SampledFrameResponse kernel must not be empty"))
        isodd(size(kernel, 1)) || throw(InvalidConfiguration("SampledFrameResponse kernel row count must be odd"))
        isodd(size(kernel, 2)) || throw(InvalidConfiguration("SampledFrameResponse kernel column count must be odd"))
        _validate_physical_response_kernel(kernel, "SampledFrameResponse")
        return new{T,A}(kernel)
    end

    function SampledFrameResponse{T,A}(kernel::A) where {
        T<:AbstractFloat,A<:AbstractMatrix{T}}
        owned = copy(kernel)
        return SampledFrameResponse{T,typeof(owned)}(
            OWNED_DETECTOR_PARAMETER, owned)
    end
end

struct NullChargeCoupling <: AbstractChargeCouplingModel end

"""
    InterpixelCapacitance(kernel; normalize=true, T=Float64, backend=CPUBackend())

Post-collection capacitive coupling between neighboring detector nodes. Unlike a
`FrameResponseModel`, this stage is applied after photon and generated-charge
statistics, so it correlates the collected charge without smoothing the
Poisson expectation before sampling.
"""
struct InterpixelCapacitance{R<:SampledFrameResponse} <: AbstractChargeCouplingModel
    response::R
end

struct RectangularPixelAperture{T<:AbstractFloat,VX<:AbstractVector{T},VY<:AbstractVector{T}} <: AbstractFrameResponse
    pitch_x_px::T
    pitch_y_px::T
    fill_factor_x::T
    fill_factor_y::T
    kernel_x::VX
    kernel_y::VY
    function RectangularPixelAperture{T,VX,VY}(
        ::OwnedDetectorParameterToken, pitch_x_px::T,
        pitch_y_px::T, fill_factor_x::T, fill_factor_y::T, kernel_x::VX,
        kernel_y::VY) where {T<:AbstractFloat,VX<:AbstractVector{T},VY<:AbstractVector{T}}
        isfinite(pitch_x_px) && pitch_x_px > zero(T) ||
            throw(InvalidConfiguration(
                "RectangularPixelAperture pitch_x_px must be finite and > 0"))
        isfinite(pitch_y_px) && pitch_y_px > zero(T) ||
            throw(InvalidConfiguration(
                "RectangularPixelAperture pitch_y_px must be finite and > 0"))
        isfinite(fill_factor_x) && zero(T) < fill_factor_x <= one(T) ||
            throw(InvalidConfiguration(
                "RectangularPixelAperture fill_factor_x must be finite and lie in (0, 1]"))
        isfinite(fill_factor_y) && zero(T) < fill_factor_y <= one(T) ||
            throw(InvalidConfiguration(
                "RectangularPixelAperture fill_factor_y must be finite and lie in (0, 1]"))
        isempty(kernel_x) && throw(InvalidConfiguration(
            "RectangularPixelAperture kernel_x must not be empty"))
        isempty(kernel_y) && throw(InvalidConfiguration(
            "RectangularPixelAperture kernel_y must not be empty"))
        isodd(length(kernel_x)) || throw(InvalidConfiguration(
            "RectangularPixelAperture kernel_x length must be odd"))
        isodd(length(kernel_y)) || throw(InvalidConfiguration(
            "RectangularPixelAperture kernel_y length must be odd"))
        _validate_physical_response_kernel(kernel_x,
            "RectangularPixelAperture kernel_x")
        _validate_physical_response_kernel(kernel_y,
            "RectangularPixelAperture kernel_y")
        return new{T,VX,VY}(pitch_x_px, pitch_y_px, fill_factor_x,
            fill_factor_y, kernel_x, kernel_y)
    end


    function RectangularPixelAperture{T,VX,VY}(pitch_x_px::T,
        pitch_y_px::T, fill_factor_x::T, fill_factor_y::T, kernel_x::VX,
        kernel_y::VY) where {T<:AbstractFloat,VX<:AbstractVector{T},VY<:AbstractVector{T}}
        owned_x = copy(kernel_x)
        owned_y = copy(kernel_y)
        return RectangularPixelAperture{T,typeof(owned_x),typeof(owned_y)}(
            OWNED_DETECTOR_PARAMETER, pitch_x_px, pitch_y_px,
            fill_factor_x, fill_factor_y, owned_x, owned_y)
    end
end

const SeparableGaussianPixelResponse = GaussianPixelResponse

@inline _frame_response_kernel_sum(kernel) = _frame_response_kernel_sum(execution_style(kernel), kernel)
@inline _frame_response_kernel_sum(::ScalarCPUStyle, kernel) = sum(kernel)
@inline _frame_response_kernel_sum(::ExecutionStyle, kernel) = sum(Array(kernel))

function _validate_physical_response_kernel(kernel::AbstractArray{T},
    label::AbstractString; require_nonamplifying::Bool=true) where {T<:AbstractFloat}
    host_kernel = Array(kernel)
    isempty(host_kernel) && throw(InvalidConfiguration(
        "$(label) kernel must not be empty"))
    kernel_sum = zero(T)
    @inbounds for weight in host_kernel
        isfinite(weight) || throw(InvalidConfiguration(
            "$(label) kernel weights must be finite"))
        weight >= zero(T) || throw(InvalidConfiguration(
            "$(label) kernel weights must be nonnegative"))
        kernel_sum += weight
    end
    isfinite(kernel_sum) && kernel_sum > zero(T) ||
        throw(InvalidConfiguration("$(label) kernel must have finite positive sum"))
    if require_nonamplifying
        tolerance = T(16) * eps(T) * T(max(length(host_kernel), 1))
        kernel_sum <= one(T) + tolerance || throw(InvalidConfiguration(
            "$(label) kernel sum must not exceed one"))
    end
    return kernel_sum
end

response_family(::NullFrameResponse) = :none
response_family(::GaussianPixelResponse) = :gaussian
response_family(::SampledFrameResponse) = :sampled
response_family(::RectangularPixelAperture) = :rectangular_aperture

frame_response_symbol(model::AbstractFrameResponse) = response_family(model)

response_application_domain(::AbstractFrameResponse) = :image

# The sampled kernels are shift invariant only in the interior. Their finite
# frame realization uses zero extension, so the full applied operator is not
# shift invariant at a boundary.
is_shift_invariant(::AbstractFrameResponse) = false
is_shift_invariant(::NullFrameResponse) = true
supports_frequency_domain_application(::AbstractFrameResponse) = false
supports_separable_application(::AbstractFrameResponse) = false
supports_separable_application(::NullFrameResponse) = true
supports_separable_application(::GaussianPixelResponse) = true
supports_separable_application(::RectangularPixelAperture) = true
supports_batched_response_application(::ScalarCPUStyle, ::AbstractFrameResponse) = true
supports_batched_response_application(::AcceleratorStyle, ::AbstractFrameResponse) = false
supports_batched_response_application(::ScalarCPUStyle, ::NullFrameResponse) = true
supports_batched_response_application(::AcceleratorStyle, ::NullFrameResponse) = true
supports_batched_response_application(::ScalarCPUStyle, ::GaussianPixelResponse) = true
supports_batched_response_application(::AcceleratorStyle, ::GaussianPixelResponse) = true
supports_batched_response_application(::ScalarCPUStyle, ::SampledFrameResponse) = true
supports_batched_response_application(::AcceleratorStyle, ::SampledFrameResponse) = true
supports_batched_response_application(::ScalarCPUStyle, ::RectangularPixelAperture) = true
supports_batched_response_application(::AcceleratorStyle, ::RectangularPixelAperture) = true
supports_subpixel_geometry(::AbstractFrameResponse) = false

response_support(::NullFrameResponse) = nothing, nothing
response_support(model::GaussianPixelResponse) = length(model.kernel), length(model.kernel)
response_support(model::SampledFrameResponse) = size(model.kernel)
response_support(model::RectangularPixelAperture) = length(model.kernel_y), length(model.kernel_x)

response_width_px(::NullFrameResponse, ::Type{T}) where {T<:AbstractFloat} = nothing
response_width_px(model::GaussianPixelResponse, ::Type{T}) where {T<:AbstractFloat} = T(model.response_width_px)
response_width_px(::SampledFrameResponse, ::Type{T}) where {T<:AbstractFloat} = nothing
response_width_px(::RectangularPixelAperture, ::Type{T}) where {T<:AbstractFloat} = nothing

response_pitch_x_px(::AbstractFrameResponse, ::Type{T}) where {T<:AbstractFloat} = nothing
response_pitch_y_px(::AbstractFrameResponse, ::Type{T}) where {T<:AbstractFloat} = nothing
response_fill_factor_x(::AbstractFrameResponse, ::Type{T}) where {T<:AbstractFloat} = nothing
response_fill_factor_y(::AbstractFrameResponse, ::Type{T}) where {T<:AbstractFloat} = nothing
response_aperture_shape(::AbstractFrameResponse) = nothing
response_aperture_shape(::SampledFrameResponse) = :sampled

response_pitch_x_px(model::RectangularPixelAperture, ::Type{T}) where {T<:AbstractFloat} = T(model.pitch_x_px)
response_pitch_y_px(model::RectangularPixelAperture, ::Type{T}) where {T<:AbstractFloat} = T(model.pitch_y_px)
response_fill_factor_x(model::RectangularPixelAperture, ::Type{T}) where {T<:AbstractFloat} = T(model.fill_factor_x)
response_fill_factor_y(model::RectangularPixelAperture, ::Type{T}) where {T<:AbstractFloat} = T(model.fill_factor_y)
response_aperture_shape(::RectangularPixelAperture) = :rectangular

supports_detector_mtf(::AbstractFrameResponse) = false
supports_detector_mtf(::GaussianPixelResponse) = true
supports_detector_mtf(::SampledFrameResponse) = true
supports_detector_mtf(::RectangularPixelAperture) = true

"""
    detector_mtf(response, spatial_frequency_x, spatial_frequency_y)

Evaluate the normalized, interior modulation transfer function of the sampled
response kernel applied by detector acquisition, at spatial frequencies in
cycles per detector pixel. This is the infinite-grid transfer function of the
realized discrete kernel. A finite detector frame uses non-amplifying zero
extension, so edge pixels do not share one global transfer function with the
interior. This diagnostic does not imply that the response is applied through
an FFT in the capture hot path. Continuous subpixel aperture MTF requires an
explicitly prepared oversampled optical-grid mapping.
"""
detector_mtf(::NullFrameResponse, spatial_frequency_x::Real,
    spatial_frequency_y::Real) = one(float(spatial_frequency_x + spatial_frequency_y))

function _sampled_axis_transfer(kernel::AbstractVector,
    spatial_frequency::Real)
    host_kernel = Array(kernel)
    T = promote_type(eltype(host_kernel),
        typeof(float(spatial_frequency)))
    frequency = T(spatial_frequency)
    center = fld(length(host_kernel), 2) + 1
    response = zero(Complex{T})
    normalization = zero(T)
    @inbounds for index in eachindex(host_kernel)
        weight = T(host_kernel[index])
        phase = -T(2pi) * frequency * T(index - center)
        response += weight * cis(phase)
        normalization += weight
    end
    return response / normalization
end

function detector_mtf(model::GaussianPixelResponse, spatial_frequency_x::Real,
    spatial_frequency_y::Real)
    response_x = _sampled_axis_transfer(model.kernel,
        spatial_frequency_x)
    response_y = _sampled_axis_transfer(model.kernel,
        spatial_frequency_y)
    return abs(response_x * response_y)
end

function detector_mtf(model::SampledFrameResponse, spatial_frequency_x::Real,
    spatial_frequency_y::Real)
    kernel = Array(model.kernel)
    T = promote_type(eltype(kernel), typeof(float(spatial_frequency_x)),
        typeof(float(spatial_frequency_y)))
    fx = T(spatial_frequency_x)
    fy = T(spatial_frequency_y)
    center_i = fld(size(kernel, 1), 2) + 1
    center_j = fld(size(kernel, 2), 2) + 1
    response = zero(Complex{T})
    normalization = zero(T)
    @inbounds for j in axes(kernel, 2), i in axes(kernel, 1)
        weight = T(kernel[i, j])
        phase = -T(2pi) * (fy * T(i - center_i) + fx * T(j - center_j))
        response += weight * cis(phase)
        normalization += weight
    end
    return abs(response) / normalization
end

function detector_mtf(model::RectangularPixelAperture,
    spatial_frequency_x::Real, spatial_frequency_y::Real)
    response_x = _sampled_axis_transfer(model.kernel_x,
        spatial_frequency_x)
    response_y = _sampled_axis_transfer(model.kernel_y,
        spatial_frequency_y)
    return abs(response_x * response_y)
end

default_response_model(::FrameSensorType; T::Type{<:AbstractFloat}=Float64, backend::AbstractArrayBackend=CPUBackend()) =
    NullFrameResponse()
default_charge_coupling_model(::FrameSensorType; T::Type{<:AbstractFloat}=Float64,
    backend::AbstractArrayBackend=CPUBackend()) = NullChargeCoupling()

charge_coupling_symbol(::NullChargeCoupling) = :none
charge_coupling_symbol(::InterpixelCapacitance) = :interpixel_capacitance
charge_coupling_support(::NullChargeCoupling) = (nothing, nothing)
charge_coupling_support(model::InterpixelCapacitance) = size(model.response.kernel)

struct NoBackground <: BackgroundModel end

struct ScalarBackground{T<:AbstractFloat} <: BackgroundModel
    level::T
end

struct BackgroundFrame{T<:AbstractFloat,A<:AbstractMatrix{T}} <: BackgroundModel
    map::A
    function BackgroundFrame{T,A}(::OwnedDetectorParameterToken,
        map::A) where {T<:AbstractFloat,A<:AbstractMatrix{T}}
        return new{T,A}(map)
    end

    function BackgroundFrame{T,A}(map::A) where {
        T<:AbstractFloat,A<:AbstractMatrix{T}}
        owned = copy(map)
        return BackgroundFrame{T,typeof(owned)}(
            OWNED_DETECTOR_PARAMETER, owned)
    end
end

BackgroundFrame(map::A) where {T<:AbstractFloat,A<:AbstractMatrix{T}} =
    BackgroundFrame{T,A}(map)

struct ScalarQuantumEfficiency{T<:AbstractFloat} <: AbstractQuantumEfficiencyModel
    value::T
    function ScalarQuantumEfficiency{T}(value::T) where {T<:AbstractFloat}
        zero(T) <= value <= one(T) ||
            throw(InvalidConfiguration("ScalarQuantumEfficiency value must lie in [0, 1]"))
        return new{T}(value)
    end
end

struct SampledQuantumEfficiency{T<:AbstractFloat,V<:AbstractVector{T}} <: AbstractQuantumEfficiencyModel
    wavelengths::V
    values::V
    out_of_band::T
    function SampledQuantumEfficiency{T,V}(
        ::OwnedDetectorParameterToken, wavelengths::V, values::V,
        out_of_band::T) where {T<:AbstractFloat,V<:AbstractVector{T}}
        length(wavelengths) >= 2 ||
            throw(InvalidConfiguration("SampledQuantumEfficiency requires at least two samples"))
        length(wavelengths) == length(values) ||
            throw(DimensionMismatchError("SampledQuantumEfficiency wavelengths and values must have the same length"))
        zero(T) <= out_of_band <= one(T) ||
            throw(InvalidConfiguration("SampledQuantumEfficiency out_of_band must lie in [0, 1]"))
        @inbounds for i in eachindex(wavelengths)
            wavelengths[i] > zero(T) ||
                throw(InvalidConfiguration("SampledQuantumEfficiency wavelengths must be > 0"))
            zero(T) <= values[i] <= one(T) ||
                throw(InvalidConfiguration("SampledQuantumEfficiency values must lie in [0, 1]"))
            if i > firstindex(wavelengths)
                wavelengths[i] > wavelengths[i - 1] ||
                    throw(InvalidConfiguration("SampledQuantumEfficiency wavelengths must be strictly increasing"))
            end
        end
        return new{T,V}(wavelengths, values, out_of_band)
    end


    function SampledQuantumEfficiency{T,V}(wavelengths::V, values::V,
        out_of_band::T) where {T<:AbstractFloat,V<:AbstractVector{T}}
        owned_wavelengths = copy(wavelengths)
        owned_values = copy(values)
        return SampledQuantumEfficiency{T,typeof(owned_wavelengths)}(
            OWNED_DETECTOR_PARAMETER, owned_wavelengths, owned_values,
            out_of_band)
    end
end

struct DetectorExportMetadata{T<:AbstractFloat}
    integration_time::T
    qe::T
    psf_sampling::Int
    binning::Int
    gain::T
    dark_current::T
    bits::Union{Nothing,Int}
    full_well::Union{Nothing,T}
    sensor::Symbol
    noise::Symbol
    readout_sigma::Union{Nothing,T}
    output_type::Union{Nothing,DataType}
    frame_size::Tuple{Int,Int}
    output_size::Tuple{Int,Int}
    frame_response::Symbol
    response_width_px::Union{Nothing,T}
    response_application_domain::Symbol
    response_is_separable::Bool
    response_is_shift_invariant::Bool
    response_support_rows::Union{Nothing,Int}
    response_support_cols::Union{Nothing,Int}
    pitch_x_px::Union{Nothing,T}
    pitch_y_px::Union{Nothing,T}
    fill_factor_x::Union{Nothing,T}
    fill_factor_y::Union{Nothing,T}
    aperture_shape::Union{Nothing,Symbol}
    charge_coupling::Symbol
    charge_coupling_support_rows::Union{Nothing,Int}
    charge_coupling_support_cols::Union{Nothing,Int}
    detector_defects::Symbol
    has_prnu::Bool
    has_dsnu::Bool
    has_bad_pixels::Bool
    window_rows::Union{Nothing,Tuple{Int,Int}}
    window_cols::Union{Nothing,Tuple{Int,Int}}
    timing_model::Symbol
    timing_line_time::Union{Nothing,T}
    acquisition_mode::Symbol
    frame_transfer_time::Union{Nothing,T}
    steady_state_frame_period::Union{Nothing,T}
    thermal_model::Symbol
    detector_temperature_K::Union{Nothing,T}
    ambient_temperature_K::Union{Nothing,T}
    cooling_setpoint_K::Union{Nothing,T}
    thermal_time_constant_s::Union{Nothing,T}
    dark_current_law::Symbol
    glow_rate_law::Symbol
    cic_rate_law::Symbol
    sampling_mode::Symbol
    sampling_reads::Union{Nothing,Int}
    sampling_reference_reads::Union{Nothing,Int}
    sampling_signal_reads::Union{Nothing,Int}
    sampling_read_time::Union{Nothing,T}
    sampling_wallclock_time::Union{Nothing,T}
    readout_correction::Symbol
    correction_edge_rows::Union{Nothing,Int}
    correction_edge_cols::Union{Nothing,Int}
    correction_group_rows::Union{Nothing,Int}
    correction_group_cols::Union{Nothing,Int}
    correction_stage_count::Int
    nonlinearity_model::Symbol
    persistence_model::Symbol
    provides_reference_frame::Bool
    provides_signal_frame::Bool
    provides_combined_frame::Bool
    provides_reference_cube::Bool
    provides_signal_cube::Bool
    provides_read_cube::Bool
    reference_cube_reads::Union{Nothing,Int}
    signal_cube_reads::Union{Nothing,Int}
    read_cube_reads::Union{Nothing,Int}
end

struct CountingReadoutMetadata
    layout::Symbol
    output_size::Tuple{Int,Int}
    n_channels::Int
end

struct CountingDetectorExportMetadata{T<:AbstractFloat}
    integration_time::T
    qe::T
    fill_factor::Union{Nothing,T}
    gain::T
    dark_count_rate::T
    dead_time_model::Symbol
    dead_time::Union{Nothing,T}
    gate_model::Symbol
    duty_cycle::Union{Nothing,T}
    correlation_model::Symbol
    afterpulse_probability::Union{Nothing,T}
    crosstalk::Union{Nothing,T}
    thermal_model::Symbol
    detector_temperature_K::Union{Nothing,T}
    ambient_temperature_K::Union{Nothing,T}
    cooling_setpoint_K::Union{Nothing,T}
    thermal_time_constant_s::Union{Nothing,T}
    dark_count_law::Symbol
    sensor::Symbol
    noise::Symbol
    output_type::Union{Nothing,DataType}
    readout::CountingReadoutMetadata
    energy_resolution::Union{Nothing,T}
    timing_jitter_s::Union{Nothing,T}
    wavelength_min_m::Union{Nothing,T}
    wavelength_max_m::Union{Nothing,T}
end

NoiseReadout(sigma::Real) = NoiseReadout{Float64}(float(sigma))
NoisePhotonReadout(sigma::Real) = NoisePhotonReadout{Float64}(float(sigma))

function _to_backend_vector(host_kernel::AbstractVector{T}, backend) where {T<:AbstractFloat}
    backend = _resolve_array_backend(backend)
    kernel = backend{T}(undef, length(host_kernel))
    copyto!(kernel, host_kernel)
    return kernel
end

function _to_backend_matrix(host_data::AbstractMatrix{T}, backend) where {T<:AbstractFloat}
    backend = _resolve_array_backend(backend)
    data = backend{T}(undef, size(host_data)...)
    copyto!(data, host_data)
    return data
end

function _to_backend_bool_matrix(host_data::AbstractMatrix{Bool}, backend)
    backend = _resolve_array_backend(backend)
    host_matrix = Matrix{Bool}(host_data)
    data = backend{Bool}(undef, size(host_data)...)
    copyto!(data, host_matrix)
    return data
end

function _gaussian_kernel(response_width_px::Real, truncate_at::Real, ::Type{T}) where {T<:AbstractFloat}
    width = T(response_width_px)
    truncation = T(truncate_at)
    isfinite(width) && width > zero(T) || throw(InvalidConfiguration(
        "GaussianPixelResponse response_width_px must be finite and > 0"))
    isfinite(truncation) && truncation > zero(T) || throw(InvalidConfiguration(
        "GaussianPixelResponse truncate_at must be finite and > 0"))
    support_radius = truncation * width
    isfinite(support_radius) || throw(InvalidConfiguration(
        "GaussianPixelResponse response support must be finite"))
    radius = max(1, ceil(Int, support_radius))
    host_kernel = Vector{T}(undef, 2 * radius + 1)
    inv_sigma2 = inv(width^2)
    for (idx, offset) in enumerate(-radius:radius)
        host_kernel[idx] = exp(-T(0.5) * T(offset * offset) * inv_sigma2)
    end
    host_kernel ./= sum(host_kernel)
    return host_kernel
end

function _pixel_aperture_kernel(pitch_px::Real, fill_factor::Real, ::Type{T}) where {T<:AbstractFloat}
    pitch = T(pitch_px)
    fill = T(fill_factor)
    isfinite(pitch) && pitch > zero(T) ||
        throw(InvalidConfiguration("pixel pitch must be finite and > 0"))
    isfinite(fill) && zero(T) < fill <= one(T) ||
        throw(InvalidConfiguration("pixel fill factor must be finite and lie in (0, 1]"))
    width = pitch * fill
    radius = max(1, ceil(Int, width))
    host_kernel = Vector{T}(undef, 2 * radius + 1)
    half_width = width / T(2)
    for (idx, offset) in enumerate(-radius:radius)
        pixel_left = T(offset) - T(0.5)
        pixel_right = T(offset) + T(0.5)
        left = max(pixel_left, -half_width)
        right = min(pixel_right, half_width)
        host_kernel[idx] = max(right - left, zero(T))
    end
    sum(host_kernel) > zero(T) || throw(InvalidConfiguration("pixel aperture kernel must have positive support"))
    host_kernel ./= sum(host_kernel)
    return host_kernel
end

function GaussianPixelResponse(; response_width_px::Real=0.5, truncate_at::Real=3.0,
    T::Type{<:AbstractFloat}=Float64, backend::AbstractArrayBackend=CPUBackend())
    backend = _resolve_array_backend(backend)
    kernel = _to_backend_vector(_gaussian_kernel(response_width_px, truncate_at, T), backend)
    return GaussianPixelResponse{T,typeof(kernel)}(
        OWNED_DETECTOR_PARAMETER, T(response_width_px), kernel)
end

function RectangularPixelAperture(; pitch_x_px::Real=1.0, pitch_y_px::Real=1.0,
    fill_factor_x::Real=1.0, fill_factor_y::Real=1.0,
    T::Type{<:AbstractFloat}=Float64, backend::AbstractArrayBackend=CPUBackend())
    backend = _resolve_array_backend(backend)
    kernel_x = _to_backend_vector(_pixel_aperture_kernel(pitch_x_px, fill_factor_x, T), backend)
    kernel_y = _to_backend_vector(_pixel_aperture_kernel(pitch_y_px, fill_factor_y, T), backend)
    return RectangularPixelAperture{T,typeof(kernel_x),typeof(kernel_y)}(
        OWNED_DETECTOR_PARAMETER, T(pitch_x_px), T(pitch_y_px),
        T(fill_factor_x), T(fill_factor_y), kernel_x, kernel_y)
end

function SampledFrameResponse(kernel::AbstractMatrix; normalize::Bool=true,
    T::Type{<:AbstractFloat}=Float64, backend::AbstractArrayBackend=CPUBackend())
    backend = _resolve_array_backend(backend)
    isempty(kernel) && throw(InvalidConfiguration("SampledFrameResponse kernel must not be empty"))
    host_kernel = T.(kernel)
    kernel_sum = _validate_physical_response_kernel(host_kernel,
        "SampledFrameResponse"; require_nonamplifying=false)
    if normalize
        host_kernel ./= kernel_sum
    end
    kernel_backend = backend{T}(undef, size(host_kernel)...)
    copyto!(kernel_backend, host_kernel)
    model = SampledFrameResponse{T,typeof(kernel_backend)}(
        OWNED_DETECTOR_PARAMETER, kernel_backend)
    return validate_frame_response_model(model)
end

function InterpixelCapacitance(kernel::AbstractMatrix; normalize::Bool=true,
    T::Type{<:AbstractFloat}=Float64, backend::AbstractArrayBackend=CPUBackend())
    response = SampledFrameResponse(kernel; normalize=normalize, T=T, backend=backend)
    return InterpixelCapacitance{typeof(response)}(response)
end

@kernel function separable_response_rows_kernel!(out, img, kernel, radius::Int, n::Int, m::Int, klen::Int)
    i, j = @index(Global, NTuple)
    if i <= n && j <= m
        acc = zero(eltype(out))
        @inbounds for kk in 1:klen
            jj = j + kk - radius - 1
            if 1 <= jj <= m
                acc += kernel[kk] * img[i, jj]
            end
        end
        @inbounds out[i, j] = acc
    end
end

@kernel function separable_response_cols_kernel!(out, img, kernel, radius::Int, n::Int, m::Int, klen::Int)
    i, j = @index(Global, NTuple)
    if i <= n && j <= m
        acc = zero(eltype(out))
        @inbounds for kk in 1:klen
            ii = i + kk - radius - 1
            if 1 <= ii <= n
                acc += kernel[kk] * img[ii, j]
            end
        end
        @inbounds out[i, j] = acc
    end
end

@kernel function sampled_response_kernel!(out, img, kernel, radius_i::Int, radius_j::Int,
    n::Int, m::Int, kn::Int, km::Int)
    i, j = @index(Global, NTuple)
    if i <= n && j <= m
        acc = zero(eltype(out))
        @inbounds for ki in 1:kn
            ii = i + ki - radius_i - 1
            if 1 <= ii <= n
                for kj in 1:km
                    jj = j + kj - radius_j - 1
                    if 1 <= jj <= m
                        acc += kernel[ki, kj] * img[ii, jj]
                    end
                end
            end
        end
        @inbounds out[i, j] = acc
    end
end

@kernel function sampled_response_stack_kernel!(out, img, kernel, radius_i::Int, radius_j::Int,
    n_batch::Int, n::Int, m::Int, kn::Int, km::Int)
    b, i, j = @index(Global, NTuple)
    if b <= n_batch && i <= n && j <= m
        acc = zero(eltype(out))
        @inbounds for ki in 1:kn
            ii = i + ki - radius_i - 1
            if 1 <= ii <= n
                for kj in 1:km
                    jj = j + kj - radius_j - 1
                    if 1 <= jj <= m
                        acc += kernel[ki, kj] * img[b, ii, jj]
                    end
                end
            end
        end
        @inbounds out[b, i, j] = acc
    end
end

@kernel function separable_response_stack_rows_kernel!(out, img, kernel, radius::Int,
    n_batch::Int, n::Int, m::Int, klen::Int)
    b, i, j = @index(Global, NTuple)
    if b <= n_batch && i <= n && j <= m
        acc = zero(eltype(out))
        @inbounds for kk in 1:klen
            jj = j + kk - radius - 1
            if 1 <= jj <= m
                acc += kernel[kk] * img[b, i, jj]
            end
        end
        @inbounds out[b, i, j] = acc
    end
end

@kernel function separable_response_stack_cols_kernel!(out, img, kernel, radius::Int,
    n_batch::Int, n::Int, m::Int, klen::Int)
    b, i, j = @index(Global, NTuple)
    if b <= n_batch && i <= n && j <= m
        acc = zero(eltype(out))
        @inbounds for kk in 1:klen
            ii = i + kk - radius - 1
            if 1 <= ii <= n
                acc += kernel[kk] * img[b, ii, j]
            end
        end
        @inbounds out[b, i, j] = acc
    end
end

@kernel function add_column_noise_kernel!(frame, noise, sigma, n::Int, m::Int)
    i, j = @index(Global, NTuple)
    if i <= n && j <= m
        @inbounds frame[i, j] += sigma * noise[1, j]
    end
end

struct DetectorParams{T<:AbstractFloat,S<:SensorType,QE<:AbstractQuantumEfficiencyModel,R<:AbstractFrameResponse,
    CC<:AbstractChargeCouplingModel,
    D<:AbstractDetectorDefectModel,FT<:AbstractFrameTimingModel,
    C<:FrameReadoutCorrectionModel,NL<:AbstractFrameNonlinearityModel,
    TM<:AbstractDetectorThermalModel}
    integration_time::T
    qe::T
    psf_sampling::Int
    binning::Int
    gain::T
    dark_current::T
    bits::Union{Nothing,Int}
    full_well::Union{Nothing,T}
    sensor::S
    quantum_efficiency_model::QE
    response_model::R
    charge_coupling_model::CC
    defect_model::D
    timing_model::FT
    correction_model::C
    nonlinearity_model::NL
    thermal_model::TM
    readout_window::Union{Nothing,FrameWindow}
    output_type::Union{Nothing,DataType}
end

mutable struct DetectorState{T<:AbstractFloat,A<:AbstractMatrix{T},O,OH,P,
    TS<:AbstractDetectorThermalState}
    frame::A
    presampling_buffer::A
    presampling_scratch::A
    response_buffer::A
    bin_buffer::A
    temporal_buffer::A
    noise_buffer::A
    noise_buffer_host::Matrix{T}
    batched_buffer_host::Array{T,3}
    accum_buffer::A
    latent_buffer::A
    output_buffer::O
    output_buffer_host::OH
    readout_products::P
    thermal_state::TS
    integrated_time::T
    readout_ready::Bool
end

struct Detector{N<:NoiseModel,P<:DetectorParams,S<:DetectorState,BF<:BackgroundModel,BM<:BackgroundModel,B<:AbstractArrayBackend} <: AbstractFrameDetector
    noise::N
    params::P
    state::S
    background_flux::BF
    background_map::BM
end

@inline backend(::Detector{<:Any,<:Any,<:Any,<:Any,<:Any,B}) where {B} = B()
