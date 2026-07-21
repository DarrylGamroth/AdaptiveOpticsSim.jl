module AdaptiveOpticsSimAppleAccelerateExt

using AbstractFFTs
using AdaptiveOpticsSim
using AppleAccelerate
using LinearAlgebra

@static if Sys.isapple()

abstract type AppleFFTDirection end
struct AppleFFTForward <: AppleFFTDirection end
struct AppleFFTInverse <: AppleFFTDirection end

"""
Type-stable Apple FFT plan with a vDSP fast path and an FFTW shape fallback.

AppleAccelerate 0.7's general AbstractFFTs adapter allocates split-complex arrays
while executing an in-place plan. AdaptiveOpticsSim owns these buffers in the
prepared plan so repeated optical propagation remains allocation-free. The
concrete fallback field also lets resized prepared state cross the vDSP shape
boundary without changing its Julia field type.
"""
struct AdaptiveOpticsAppleFFTPlan{
    T<:Union{Float32,Float64},N,D<:AppleFFTDirection,P,
} <: AbstractFFTs.Plan{Complex{T}}
    setup::Union{Nothing,AppleAccelerate.FFTSetup{T}}
    real_buffer::Array{T,N}
    imag_buffer::Array{T,N}
    fftw_plan::Union{Nothing,P}
    shape::NTuple{N,Int}
end

Base.size(plan::AdaptiveOpticsAppleFFTPlan) = plan.shape

@inline _vdsp_shape_supported(buffer::Vector) =
    !isempty(buffer) && ispow2(length(buffer))
@inline _vdsp_shape_supported(buffer::Matrix) =
    !isempty(buffer) && ispow2(size(buffer, 1)) && ispow2(size(buffer, 2))

@inline function _full_transform_region(dims, ::Val{N}) where {N}
    dims isa Integer && return N == 1 && dims == 1
    return Tuple(dims) == ntuple(identity, N)
end

@inline function _fallback_plan(
    ::Type{AppleFFTForward}, buffer::Array{Complex{T},N}, dims,
) where {T,N}
    return AdaptiveOpticsSim._plan_fftw_fft!(buffer, dims)
end

@inline function _fallback_plan(
    ::Type{AppleFFTInverse}, buffer::Array{Complex{T},N}, dims,
) where {T,N}
    scale = one(T) / AdaptiveOpticsSim._fft_region_length(buffer, dims)
    return scale * AdaptiveOpticsSim._plan_fftw_bfft!(buffer, dims)
end

function _prepare_apple_fft_plan(
    ::Type{D}, buffer::Array{Complex{T},N}, dims, use_vdsp::Bool,
) where {T<:Union{Float32,Float64},N,D<:AppleFFTDirection}
    # Preparing the concrete fallback provides a type witness that keeps this
    # wrapper stable when optical state is resized across the vDSP capability
    # boundary. Supported transforms do not retain or execute that FFTW plan.
    fallback_plan = _fallback_plan(D, buffer, dims)
    fftw_plan = use_vdsp ? nothing : fallback_plan
    setup = use_vdsp ? AppleAccelerate.plan_fft(buffer) : nothing
    buffer_shape = use_vdsp ? size(buffer) : ntuple(_ -> 0, N)
    real_buffer = Array{T,N}(undef, buffer_shape)
    imag_buffer = similar(real_buffer)
    return AdaptiveOpticsAppleFFTPlan{T,N,D,typeof(fallback_plan)}(
        setup, real_buffer, imag_buffer, fftw_plan, size(buffer))
end

function _plan_forward(buffer::Array{Complex{T},N}, dims) where {
    T<:Union{Float32,Float64},N,
}
    if _vdsp_shape_supported(buffer) &&
            _full_transform_region(dims, Val(N))
        return _prepare_apple_fft_plan(AppleFFTForward, buffer, dims, true)
    end
    return _prepare_apple_fft_plan(AppleFFTForward, buffer, dims, false)
end

function _plan_inverse(buffer::Array{Complex{T},N}, dims) where {
    T<:Union{Float32,Float64},N,
}
    if _vdsp_shape_supported(buffer) &&
            _full_transform_region(dims, Val(N))
        return _prepare_apple_fft_plan(AppleFFTInverse, buffer, dims, true)
    end
    return _prepare_apple_fft_plan(AppleFFTInverse, buffer, dims, false)
end

function AdaptiveOpticsSim.plan_fft_backend!(
    buffer::Vector{Complex{T}},
) where {T<:Union{Float32,Float64}}
    return _plan_forward(buffer, (1,))
end

function AdaptiveOpticsSim.plan_fft_backend!(
    buffer::Matrix{Complex{T}},
) where {T<:Union{Float32,Float64}}
    return _plan_forward(buffer, (1, 2))
end

function AdaptiveOpticsSim.plan_fft_backend!(
    buffer::Vector{Complex{T}}, dims,
) where {T<:Union{Float32,Float64}}
    return _plan_forward(buffer, dims)
end

function AdaptiveOpticsSim.plan_fft_backend!(
    buffer::Matrix{Complex{T}}, dims,
) where {T<:Union{Float32,Float64}}
    return _plan_forward(buffer, dims)
end

function AdaptiveOpticsSim.plan_ifft_backend!(
    buffer::Vector{Complex{T}},
) where {T<:Union{Float32,Float64}}
    return _plan_inverse(buffer, (1,))
end

function AdaptiveOpticsSim.plan_ifft_backend!(
    buffer::Matrix{Complex{T}},
) where {T<:Union{Float32,Float64}}
    return _plan_inverse(buffer, (1, 2))
end

function AdaptiveOpticsSim.plan_ifft_backend!(
    buffer::Vector{Complex{T}}, dims,
) where {T<:Union{Float32,Float64}}
    return _plan_inverse(buffer, dims)
end

function AdaptiveOpticsSim.plan_ifft_backend!(
    buffer::Matrix{Complex{T}}, dims,
) where {T<:Union{Float32,Float64}}
    return _plan_inverse(buffer, dims)
end

@inline _vdsp_direction(::AdaptiveOpticsAppleFFTPlan{T,N,AppleFFTForward}) where {T,N} =
    AppleAccelerate.FFT_FORWARD
@inline _vdsp_direction(::AdaptiveOpticsAppleFFTPlan{T,N,AppleFFTInverse}) where {T,N} =
    AppleAccelerate.FFT_INVERSE
@inline _output_scale(::AdaptiveOpticsAppleFFTPlan{T,N,AppleFFTForward}) where {T,N} =
    one(T)
@inline _output_scale(plan::AdaptiveOpticsAppleFFTPlan{T,N,AppleFFTInverse}) where {T,N} =
    one(T) / length(plan.real_buffer)

@inline function _load_split_buffers!(
    plan::AdaptiveOpticsAppleFFTPlan{T,N}, buffer::Array{Complex{T},N},
) where {T,N}
    real_buffer = plan.real_buffer
    imag_buffer = plan.imag_buffer
    @inbounds @simd for index in eachindex(buffer)
        value = buffer[index]
        real_buffer[index] = real(value)
        imag_buffer[index] = imag(value)
    end
    return nothing
end

@inline function _store_split_buffers!(
    buffer::Array{Complex{T},N}, plan::AdaptiveOpticsAppleFFTPlan{T,N},
) where {T,N}
    real_buffer = plan.real_buffer
    imag_buffer = plan.imag_buffer
    scale = _output_scale(plan)
    @inbounds @simd for index in eachindex(buffer)
        buffer[index] = Complex{T}(
            real_buffer[index] * scale, imag_buffer[index] * scale)
    end
    return buffer
end

@inline function _execute_vdsp!(
    plan::AdaptiveOpticsAppleFFTPlan{Float32,1},
)
    real_buffer = plan.real_buffer
    imag_buffer = plan.imag_buffer
    setup = plan.setup::AppleAccelerate.FFTSetup{Float32}
    GC.@preserve real_buffer imag_buffer setup begin
        split = AppleAccelerate.DSPSplitComplex(
            pointer(real_buffer), pointer(imag_buffer))
        AppleAccelerate.LibAccelerate.vDSP_fft_zip(
            setup.plan, Ref(split), AppleAccelerate.SIGNAL_STRIDE,
            trailing_zeros(length(real_buffer)), _vdsp_direction(plan))
    end
    return nothing
end

@inline function _execute_vdsp!(
    plan::AdaptiveOpticsAppleFFTPlan{Float64,1},
)
    real_buffer = plan.real_buffer
    imag_buffer = plan.imag_buffer
    setup = plan.setup::AppleAccelerate.FFTSetup{Float64}
    GC.@preserve real_buffer imag_buffer setup begin
        split = AppleAccelerate.DSPDoubleSplitComplex(
            pointer(real_buffer), pointer(imag_buffer))
        AppleAccelerate.LibAccelerate.vDSP_fft_zipD(
            setup.plan, Ref(split), AppleAccelerate.SIGNAL_STRIDE,
            trailing_zeros(length(real_buffer)), _vdsp_direction(plan))
    end
    return nothing
end

@inline function _execute_vdsp!(
    plan::AdaptiveOpticsAppleFFTPlan{Float32,2},
)
    real_buffer = plan.real_buffer
    imag_buffer = plan.imag_buffer
    setup = plan.setup::AppleAccelerate.FFTSetup{Float32}
    rows, columns = size(real_buffer)
    GC.@preserve real_buffer imag_buffer setup begin
        split = AppleAccelerate.DSPSplitComplex(
            pointer(real_buffer), pointer(imag_buffer))
        AppleAccelerate.LibAccelerate.vDSP_fft2d_zip(
            setup.plan, Ref(split), AppleAccelerate.SIGNAL_STRIDE, 0,
            trailing_zeros(rows), trailing_zeros(columns),
            _vdsp_direction(plan))
    end
    return nothing
end

@inline function _execute_vdsp!(
    plan::AdaptiveOpticsAppleFFTPlan{Float64,2},
)
    real_buffer = plan.real_buffer
    imag_buffer = plan.imag_buffer
    setup = plan.setup::AppleAccelerate.FFTSetup{Float64}
    rows, columns = size(real_buffer)
    GC.@preserve real_buffer imag_buffer setup begin
        split = AppleAccelerate.DSPDoubleSplitComplex(
            pointer(real_buffer), pointer(imag_buffer))
        AppleAccelerate.LibAccelerate.vDSP_fft2d_zipD(
            setup.plan, Ref(split), AppleAccelerate.SIGNAL_STRIDE, 0,
            trailing_zeros(rows), trailing_zeros(columns),
            _vdsp_direction(plan))
    end
    return nothing
end

@inline function _execute_plan!(
    buffer::Array{Complex{T},N}, plan::AdaptiveOpticsAppleFFTPlan{T,N},
) where {T,N}
    size(buffer) == size(plan) || throw(DimensionMismatch(
        "Apple FFT plan size $(size(plan)) does not match buffer size $(size(buffer))"))
    if isnothing(plan.setup)
        fallback_plan = plan.fftw_plan
        isnothing(fallback_plan) && throw(ArgumentError(
            "Apple FFT plan has neither a vDSP setup nor an FFTW fallback"))
        return AdaptiveOpticsSim.execute_fft_plan!(buffer, fallback_plan)
    end
    _load_split_buffers!(plan, buffer)
    _execute_vdsp!(plan)
    return _store_split_buffers!(buffer, plan)
end

function AdaptiveOpticsSim.execute_fft_plan!(
    buffer::Array{Complex{T},N}, plan::AdaptiveOpticsAppleFFTPlan{T,N},
) where {T,N}
    return _execute_plan!(buffer, plan)
end

function LinearAlgebra.mul!(
    destination::Array{Complex{T},N},
    plan::AdaptiveOpticsAppleFFTPlan{T,N},
    source::Array{Complex{T},N},
) where {T,N}
    size(source) == size(plan) || throw(DimensionMismatch(
        "Apple FFT plan size $(size(plan)) does not match source size $(size(source))"))
    destination === source || copyto!(destination, source)
    return _execute_plan!(destination, plan)
end

end # Sys.isapple()

end
