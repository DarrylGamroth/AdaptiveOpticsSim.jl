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
abstract type AbstractFrameMTF <: AbstractFrameResponse end
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
abstract type AvalancheFrameSensorType <: FrameSensorType end
abstract type HgCdTeAvalancheArraySensorType <: AvalancheFrameSensorType end

struct FrameWindow
    rows::UnitRange{Int}
    cols::UnitRange{Int}
    function FrameWindow(rows::UnitRange{Int}, cols::UnitRange{Int})
        window = new(rows, cols)
        validate_readout_window(window)
        return window
    end
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

struct NoFrameReadoutProducts <: FrameReadoutProducts end

struct NullDetectorDefectModel <: AbstractDetectorDefectModel end

struct PixelResponseNonuniformity{T<:AbstractFloat,A<:AbstractMatrix{T}} <: AbstractDetectorDefectModel
    gain_map::A
end

struct DarkSignalNonuniformity{T<:AbstractFloat,A<:AbstractMatrix{T}} <: AbstractDetectorDefectModel
    dark_map::A
end

struct BadPixelMask{T<:AbstractFloat,A<:AbstractMatrix{Bool}} <: AbstractDetectorDefectModel
    mask::A
    throughput::T
end

struct CompositeDetectorDefectModel{M<:Tuple} <: AbstractDetectorDefectModel
    stages::M
    function CompositeDetectorDefectModel(stages::Tuple{Vararg{AbstractDetectorDefectModel}})
        isempty(stages) && throw(InvalidConfiguration("CompositeDetectorDefectModel requires at least one stage"))
        return new{typeof(stages)}(stages)
    end
end

CompositeDetectorDefectModel(stages::Tuple) =
    throw(InvalidConfiguration("CompositeDetectorDefectModel stages must be AbstractDetectorDefectModel values"))
CompositeDetectorDefectModel(stages::AbstractDetectorDefectModel...) = CompositeDetectorDefectModel(tuple(stages...))

struct GlobalShutter <: AbstractFrameTimingModel end

struct RollingShutter{T<:AbstractFloat} <: AbstractFrameTimingModel
    line_time::T
end

RollingShutter(line_time::Real) = RollingShutter{Float64}(float(line_time))

struct NullFrameNonlinearity <: AbstractFrameNonlinearityModel end

struct SaturatingFrameNonlinearity{T<:AbstractFloat} <: AbstractFrameNonlinearityModel
    coefficient::T
end

SaturatingFrameNonlinearity(coefficient::Real) = SaturatingFrameNonlinearity{Float64}(float(coefficient))

struct NullPersistence <: AbstractPersistenceModel end

struct NullDetectorThermalModel <: AbstractDetectorThermalModel end

struct NoThermalState <: AbstractDetectorThermalState end

struct DetectorThermalState{T<:AbstractFloat} <: AbstractDetectorThermalState
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

struct ExponentialPersistence{T<:AbstractFloat} <: AbstractPersistenceModel
    coupling::T
    decay::T
end

ExponentialPersistence(coupling::Real, decay::Real) = ExponentialPersistence{Float64}(float(coupling), float(decay))

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

function PixelResponseNonuniformity(gain_map::AbstractMatrix; T::Type{<:AbstractFloat}=Float64, backend=Array)
    backend_map = _to_backend_matrix(T.(gain_map), backend)
    return validate_detector_defect_model(PixelResponseNonuniformity{T,typeof(backend_map)}(backend_map))
end

function DarkSignalNonuniformity(dark_map::AbstractMatrix; T::Type{<:AbstractFloat}=Float64, backend=Array)
    backend_map = _to_backend_matrix(T.(dark_map), backend)
    return validate_detector_defect_model(DarkSignalNonuniformity{T,typeof(backend_map)}(backend_map))
end

function BadPixelMask(mask::AbstractMatrix{Bool}; throughput::Real=0.0, T::Type{<:AbstractFloat}=Float64, backend=Array)
    backend_mask = _to_backend_bool_matrix(mask, backend)
    return validate_detector_defect_model(BadPixelMask{T,typeof(backend_mask)}(backend_mask, T(throughput)))
end

struct SampledFrameReadoutProducts{A<:AbstractMatrix,C} <: FrameReadoutProducts
    reference_frame::Union{Nothing,A}
    signal_frame::A
    read_cube::Union{Nothing,C}
end

struct HgCdTeReadoutProducts{A<:AbstractMatrix,C,V} <: FrameReadoutProducts
    reference_frame::Union{Nothing,A}
    signal_frame::A
    combined_frame::A
    reference_cube::Union{Nothing,C}
    signal_cube::Union{Nothing,C}
    read_cube::Union{Nothing,C}
    read_times::Union{Nothing,V}
end

function HgCdTeReadoutProducts(reference_frame::Union{Nothing,A}, signal_frame::A, combined_frame::A,
    reference_cube::Nothing, signal_cube::Nothing, read_cube::Nothing, read_times::Nothing) where {A<:AbstractMatrix}
    return HgCdTeReadoutProducts{A,Nothing,Nothing}(reference_frame, signal_frame, combined_frame,
        reference_cube, signal_cube, read_cube, read_times)
end

function HgCdTeReadoutProducts(reference_frame::Union{Nothing,A}, signal_frame::A, combined_frame::A,
    reference_cube::Union{Nothing,C}, signal_cube::Union{Nothing,C}, read_cube::Union{Nothing,C}, read_times::Nothing) where
    {A<:AbstractMatrix,C<:AbstractArray}
    return HgCdTeReadoutProducts{A,C,Nothing}(reference_frame, signal_frame, combined_frame,
        reference_cube, signal_cube, read_cube, read_times)
end

function HgCdTeReadoutProducts(reference_frame::Union{Nothing,A}, signal_frame::A, combined_frame::A,
    reference_cube::Union{Nothing,C}, signal_cube::Union{Nothing,C}, read_cube::Union{Nothing,C}, read_times::Union{Nothing,V}) where
    {A<:AbstractMatrix,C<:AbstractArray,V<:AbstractVector}
    return HgCdTeReadoutProducts{A,C,V}(reference_frame, signal_frame, combined_frame,
        reference_cube, signal_cube, read_cube, read_times)
end

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

struct NullFrameResponse <: AbstractFrameResponse end

struct GaussianPixelResponse{T<:AbstractFloat,V<:AbstractVector{T}} <: AbstractFrameResponse
    response_width_px::T
    kernel::V
end

struct SampledFrameResponse{T<:AbstractFloat,A<:AbstractMatrix{T}} <: AbstractFrameResponse
    kernel::A
    function SampledFrameResponse{T,A}(kernel::A) where {T<:AbstractFloat,A<:AbstractMatrix{T}}
        all(size(kernel) .> 0) || throw(InvalidConfiguration("SampledFrameResponse kernel must not be empty"))
        isodd(size(kernel, 1)) || throw(InvalidConfiguration("SampledFrameResponse kernel row count must be odd"))
        isodd(size(kernel, 2)) || throw(InvalidConfiguration("SampledFrameResponse kernel column count must be odd"))
        sum(kernel) > zero(T) || throw(InvalidConfiguration("SampledFrameResponse kernel must have positive sum"))
        return new{T,A}(kernel)
    end
end

struct RectangularPixelAperture{T<:AbstractFloat,VX<:AbstractVector{T},VY<:AbstractVector{T}} <: AbstractFrameMTF
    pitch_x_px::T
    pitch_y_px::T
    fill_factor_x::T
    fill_factor_y::T
    kernel_x::VX
    kernel_y::VY
end

struct SeparablePixelMTF{T<:AbstractFloat,VX<:AbstractVector{T},VY<:AbstractVector{T}} <: AbstractFrameMTF
    pitch_x_px::T
    pitch_y_px::T
    fill_factor_x::T
    fill_factor_y::T
    kernel_x::VX
    kernel_y::VY
end

const SeparableGaussianPixelResponse = GaussianPixelResponse

struct NoBackground <: BackgroundModel end

struct ScalarBackground{T<:AbstractFloat} <: BackgroundModel
    level::T
end

struct BackgroundFrame{T<:AbstractFloat,A<:AbstractMatrix{T}} <: BackgroundModel
    map::A
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
    output_precision::Union{Nothing,DataType}
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
    detector_defects::Symbol
    has_prnu::Bool
    has_dsnu::Bool
    has_bad_pixels::Bool
    window_rows::Union{Nothing,Tuple{Int,Int}}
    window_cols::Union{Nothing,Tuple{Int,Int}}
    timing_model::Symbol
    timing_line_time::Union{Nothing,T}
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
    output_precision::Union{Nothing,DataType}
    readout::CountingReadoutMetadata
end

NoiseReadout(sigma::Real) = NoiseReadout{Float64}(float(sigma))
NoisePhotonReadout(sigma::Real) = NoisePhotonReadout{Float64}(float(sigma))

function _to_backend_vector(host_kernel::AbstractVector{T}, backend) where {T<:AbstractFloat}
    kernel = backend{T}(undef, length(host_kernel))
    copyto!(kernel, host_kernel)
    return kernel
end

function _to_backend_matrix(host_data::AbstractMatrix{T}, backend) where {T<:AbstractFloat}
    data = backend{T}(undef, size(host_data)...)
    copyto!(data, host_data)
    return data
end

function _to_backend_bool_matrix(host_data::AbstractMatrix{Bool}, backend)
    data = backend{Bool}(undef, size(host_data)...)
    copyto!(data, host_data)
    return data
end

function _gaussian_kernel(response_width_px::Real, truncate_at::Real, ::Type{T}) where {T<:AbstractFloat}
    response_width_px > 0 || throw(InvalidConfiguration("GaussianPixelResponse response_width_px must be > 0"))
    truncate_at > 0 || throw(InvalidConfiguration("GaussianPixelResponse truncate_at must be > 0"))
    radius = max(1, ceil(Int, truncate_at * response_width_px))
    host_kernel = Vector{T}(undef, 2 * radius + 1)
    inv_sigma2 = inv(T(response_width_px)^2)
    for (idx, offset) in enumerate(-radius:radius)
        host_kernel[idx] = exp(-T(0.5) * T(offset * offset) * inv_sigma2)
    end
    host_kernel ./= sum(host_kernel)
    return host_kernel
end

function _pixel_aperture_kernel(pitch_px::Real, fill_factor::Real, ::Type{T}) where {T<:AbstractFloat}
    pitch_px > 0 || throw(InvalidConfiguration("pixel pitch must be > 0"))
    (zero(fill_factor) < fill_factor <= one(fill_factor)) ||
        throw(InvalidConfiguration("pixel fill factor must lie in (0, 1]"))
    width = T(pitch_px * fill_factor)
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
    T::Type{<:AbstractFloat}=Float64, backend=Array)
    kernel = _to_backend_vector(_gaussian_kernel(response_width_px, truncate_at, T), backend)
    return GaussianPixelResponse{T,typeof(kernel)}(T(response_width_px), kernel)
end

function RectangularPixelAperture(; pitch_x_px::Real=1.0, pitch_y_px::Real=1.0,
    fill_factor_x::Real=1.0, fill_factor_y::Real=1.0,
    T::Type{<:AbstractFloat}=Float64, backend=Array)
    kernel_x = _to_backend_vector(_pixel_aperture_kernel(pitch_x_px, fill_factor_x, T), backend)
    kernel_y = _to_backend_vector(_pixel_aperture_kernel(pitch_y_px, fill_factor_y, T), backend)
    return RectangularPixelAperture{T,typeof(kernel_x),typeof(kernel_y)}(
        T(pitch_x_px), T(pitch_y_px), T(fill_factor_x), T(fill_factor_y), kernel_x, kernel_y)
end

function SeparablePixelMTF(; pitch_x_px::Real=1.0, pitch_y_px::Real=1.0,
    fill_factor_x::Real=1.0, fill_factor_y::Real=1.0,
    T::Type{<:AbstractFloat}=Float64, backend=Array)
    kernel_x = _to_backend_vector(_pixel_aperture_kernel(pitch_x_px, fill_factor_x, T), backend)
    kernel_y = _to_backend_vector(_pixel_aperture_kernel(pitch_y_px, fill_factor_y, T), backend)
    return SeparablePixelMTF{T,typeof(kernel_x),typeof(kernel_y)}(
        T(pitch_x_px), T(pitch_y_px), T(fill_factor_x), T(fill_factor_y), kernel_x, kernel_y)
end

function SampledFrameResponse(kernel::AbstractMatrix; normalize::Bool=true,
    T::Type{<:AbstractFloat}=Float64, backend=Array)
    isempty(kernel) && throw(InvalidConfiguration("SampledFrameResponse kernel must not be empty"))
    host_kernel = T.(kernel)
    if normalize
        kernel_sum = sum(host_kernel)
        kernel_sum > zero(T) || throw(InvalidConfiguration("SampledFrameResponse normalized kernel must have positive sum"))
        host_kernel ./= kernel_sum
    end
    kernel_backend = backend{T}(undef, size(host_kernel)...)
    copyto!(kernel_backend, host_kernel)
    model = SampledFrameResponse{T,typeof(kernel_backend)}(kernel_backend)
    return validate_frame_response_model(model)
end

@kernel function separable_response_rows_kernel!(out, img, kernel, radius::Int, n::Int, m::Int, klen::Int)
    i, j = @index(Global, NTuple)
    if i <= n && j <= m
        acc = zero(eltype(out))
        @inbounds for kk in 1:klen
            jj = clamp(j + kk - radius - 1, 1, m)
            acc += kernel[kk] * img[i, jj]
        end
        @inbounds out[i, j] = acc
    end
end

@kernel function separable_response_cols_kernel!(out, img, kernel, radius::Int, n::Int, m::Int, klen::Int)
    i, j = @index(Global, NTuple)
    if i <= n && j <= m
        acc = zero(eltype(out))
        @inbounds for kk in 1:klen
            ii = clamp(i + kk - radius - 1, 1, n)
            acc += kernel[kk] * img[ii, j]
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
            ii = clamp(i + ki - radius_i - 1, 1, n)
            for kj in 1:km
                jj = clamp(j + kj - radius_j - 1, 1, m)
                acc += kernel[ki, kj] * img[ii, jj]
            end
        end
        @inbounds out[i, j] = acc
    end
end

@kernel function separable_response_stack_rows_kernel!(out, img, kernel, radius::Int,
    n_batch::Int, n::Int, m::Int, klen::Int)
    b, i, j = @index(Global, NTuple)
    if b <= n_batch && i <= n && j <= m
        acc = zero(eltype(out))
        @inbounds for kk in 1:klen
            jj = clamp(j + kk - radius - 1, 1, m)
            acc += kernel[kk] * img[b, i, jj]
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
            ii = clamp(i + kk - radius - 1, 1, n)
            acc += kernel[kk] * img[b, ii, j]
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

struct DetectorParams{T<:AbstractFloat,S<:SensorType,R<:AbstractFrameResponse,TM<:AbstractDetectorThermalModel}
    integration_time::T
    qe::T
    psf_sampling::Int
    binning::Int
    gain::T
    dark_current::T
    bits::Union{Nothing,Int}
    full_well::Union{Nothing,T}
    sensor::S
    response_model::R
    defect_model::AbstractDetectorDefectModel
    timing_model::AbstractFrameTimingModel
    correction_model::FrameReadoutCorrectionModel
    nonlinearity_model::AbstractFrameNonlinearityModel
    thermal_model::TM
    readout_window::Union{Nothing,FrameWindow}
    output_precision::Union{Nothing,DataType}
end

mutable struct DetectorState{T<:AbstractFloat,A<:AbstractMatrix{T},O,P<:FrameReadoutProducts,TS<:AbstractDetectorThermalState}
    frame::A
    response_buffer::A
    bin_buffer::A
    noise_buffer::A
    accum_buffer::A
    latent_buffer::A
    output_buffer::O
    readout_products::P
    thermal_state::TS
    integrated_time::T
    readout_ready::Bool
end

struct Detector{N<:NoiseModel,P<:DetectorParams,S<:DetectorState,BF<:BackgroundModel,BM<:BackgroundModel} <: AbstractFrameDetector
    noise::N
    params::P
    state::S
    background_flux::BF
    background_map::BM
end
