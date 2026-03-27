abstract type NoiseModel end
abstract type SensorType end
abstract type FrameSensorType <: SensorType end
abstract type CountingSensorType <: SensorType end
abstract type AbstractFrameDetector <: AbstractDetector end
abstract type AbstractCountingDetector <: AbstractDetector end
abstract type CountingDeadTimeModel end
abstract type AbstractDetectorResponse end
abstract type AbstractFrameResponse <: AbstractDetectorResponse end
abstract type AbstractFrameMTF <: AbstractFrameResponse end
const FrameResponseModel = AbstractFrameResponse
abstract type BackgroundModel end
abstract type FrameSamplingMode end
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
    window_rows::Union{Nothing,Tuple{Int,Int}}
    window_cols::Union{Nothing,Tuple{Int,Int}}
    sampling_mode::Symbol
    sampling_reads::Union{Nothing,Int}
    sampling_reference_reads::Union{Nothing,Int}
    sampling_signal_reads::Union{Nothing,Int}
    sampling_read_time::Union{Nothing,T}
    sampling_wallclock_time::Union{Nothing,T}
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

struct DetectorParams{T<:AbstractFloat,S<:SensorType,R<:AbstractFrameResponse}
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
    readout_window::Union{Nothing,FrameWindow}
    output_precision::Union{Nothing,DataType}
end

mutable struct DetectorState{T<:AbstractFloat,A<:AbstractMatrix{T},O}
    frame::A
    response_buffer::A
    bin_buffer::A
    noise_buffer::A
    accum_buffer::A
    output_buffer::O
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
