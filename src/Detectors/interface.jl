abstract type NoiseModel end
abstract type SensorType end
abstract type FrameSensorType <: SensorType end
abstract type CountingSensorType <: SensorType end
abstract type AbstractFrameDetector <: AbstractDetector end
abstract type AbstractCountingDetector <: AbstractDetector end
abstract type CountingDeadTimeModel end
abstract type FrameResponseModel end
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

struct NullFrameResponse <: FrameResponseModel end

struct SeparableGaussianPixelResponse{T<:AbstractFloat,V<:AbstractVector{T}} <: FrameResponseModel
    response_width_px::T
    kernel::V
end

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

function SeparableGaussianPixelResponse(; response_width_px::Real=0.5, truncate_at::Real=3.0,
    T::Type{<:AbstractFloat}=Float64, backend=Array)
    response_width_px > 0 || throw(InvalidConfiguration("SeparableGaussianPixelResponse response_width_px must be > 0"))
    truncate_at > 0 || throw(InvalidConfiguration("SeparableGaussianPixelResponse truncate_at must be > 0"))
    radius = max(1, ceil(Int, truncate_at * response_width_px))
    host_kernel = Vector{T}(undef, 2 * radius + 1)
    inv_sigma2 = inv(T(response_width_px)^2)
    for (idx, offset) in enumerate(-radius:radius)
        host_kernel[idx] = exp(-T(0.5) * T(offset * offset) * inv_sigma2)
    end
    host_kernel ./= sum(host_kernel)
    kernel = backend{T}(undef, length(host_kernel))
    copyto!(kernel, host_kernel)
    return SeparableGaussianPixelResponse{T,typeof(kernel)}(T(response_width_px), kernel)
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

@kernel function add_column_noise_kernel!(frame, noise, sigma, n::Int, m::Int)
    i, j = @index(Global, NTuple)
    if i <= n && j <= m
        @inbounds frame[i, j] += sigma * noise[1, j]
    end
end

struct DetectorParams{T<:AbstractFloat,S<:SensorType,R<:FrameResponseModel}
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
