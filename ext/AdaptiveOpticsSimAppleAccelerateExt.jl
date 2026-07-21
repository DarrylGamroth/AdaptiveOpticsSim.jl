module AdaptiveOpticsSimAppleAccelerateExt

using AbstractFFTs
using AdaptiveOpticsSim
using AppleAccelerate
using LinearAlgebra

@static if Sys.isapple()

abstract type VDSPTransformDirection end
struct VDSPForward <: VDSPTransformDirection end
struct VDSPInverse <: VDSPTransformDirection end

"""
Reusable Apple vDSP FFT setup with package-owned split-complex work buffers.

AppleAccelerate 0.7's general AbstractFFTs adapter allocates split-complex arrays
while executing an in-place plan. AdaptiveOpticsSim owns these buffers in the
prepared plan so repeated optical propagation remains allocation-free.
"""
struct AdaptiveOpticsVDSPPlan{
    T<:Union{Float32,Float64},N,D<:VDSPTransformDirection,
} <: AbstractFFTs.Plan{Complex{T}}
    setup::AppleAccelerate.FFTSetup{T}
    real_buffer::Array{T,N}
    imag_buffer::Array{T,N}
end

Base.size(plan::AdaptiveOpticsVDSPPlan) = size(plan.real_buffer)

@inline _vdsp_shape_supported(buffer::Vector) =
    !isempty(buffer) && ispow2(length(buffer))
@inline _vdsp_shape_supported(buffer::Matrix) =
    !isempty(buffer) && ispow2(size(buffer, 1)) && ispow2(size(buffer, 2))

@inline function _full_transform_region(dims, ::Val{N}) where {N}
    dims isa Integer && return N == 1 && dims == 1
    return Tuple(dims) == ntuple(identity, N)
end

function _prepare_vdsp_plan(
    ::Type{D}, buffer::Array{Complex{T},N},
) where {T<:Union{Float32,Float64},N,D<:VDSPTransformDirection}
    setup = AppleAccelerate.plan_fft(buffer)
    real_buffer = Array{T,N}(undef, size(buffer))
    imag_buffer = similar(real_buffer)
    return AdaptiveOpticsVDSPPlan{T,N,D}(
        setup, real_buffer, imag_buffer)
end

function _plan_forward(buffer::Array{Complex{T},N}, dims) where {
    T<:Union{Float32,Float64},N,
}
    if _vdsp_shape_supported(buffer) &&
            _full_transform_region(dims, Val(N))
        return _prepare_vdsp_plan(VDSPForward, buffer)
    end
    return AdaptiveOpticsSim._plan_fftw_fft!(buffer, dims)
end

function _plan_inverse(buffer::Array{Complex{T},N}, dims) where {
    T<:Union{Float32,Float64},N,
}
    if _vdsp_shape_supported(buffer) &&
            _full_transform_region(dims, Val(N))
        return _prepare_vdsp_plan(VDSPInverse, buffer)
    end
    scale = one(T) / AdaptiveOpticsSim._fft_region_length(buffer, dims)
    return scale * AdaptiveOpticsSim._plan_fftw_bfft!(buffer, dims)
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

@inline _vdsp_direction(::AdaptiveOpticsVDSPPlan{T,N,VDSPForward}) where {T,N} =
    AppleAccelerate.FFT_FORWARD
@inline _vdsp_direction(::AdaptiveOpticsVDSPPlan{T,N,VDSPInverse}) where {T,N} =
    AppleAccelerate.FFT_INVERSE
@inline _output_scale(::AdaptiveOpticsVDSPPlan{T,N,VDSPForward}) where {T,N} =
    one(T)
@inline _output_scale(plan::AdaptiveOpticsVDSPPlan{T,N,VDSPInverse}) where {T,N} =
    one(T) / length(plan.real_buffer)

@inline function _load_split_buffers!(
    plan::AdaptiveOpticsVDSPPlan{T,N}, buffer::Array{Complex{T},N},
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
    buffer::Array{Complex{T},N}, plan::AdaptiveOpticsVDSPPlan{T,N},
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
    plan::AdaptiveOpticsVDSPPlan{Float32,1},
)
    real_buffer = plan.real_buffer
    imag_buffer = plan.imag_buffer
    setup = plan.setup
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
    plan::AdaptiveOpticsVDSPPlan{Float64,1},
)
    real_buffer = plan.real_buffer
    imag_buffer = plan.imag_buffer
    setup = plan.setup
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
    plan::AdaptiveOpticsVDSPPlan{Float32,2},
)
    real_buffer = plan.real_buffer
    imag_buffer = plan.imag_buffer
    setup = plan.setup
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
    plan::AdaptiveOpticsVDSPPlan{Float64,2},
)
    real_buffer = plan.real_buffer
    imag_buffer = plan.imag_buffer
    setup = plan.setup
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
    buffer::Array{Complex{T},N}, plan::AdaptiveOpticsVDSPPlan{T,N},
) where {T,N}
    size(buffer) == size(plan) || throw(DimensionMismatch(
        "vDSP FFT plan size $(size(plan)) does not match buffer size $(size(buffer))"))
    _load_split_buffers!(plan, buffer)
    _execute_vdsp!(plan)
    return _store_split_buffers!(buffer, plan)
end

function AdaptiveOpticsSim.execute_fft_plan!(
    buffer::Array{Complex{T},N}, plan::AdaptiveOpticsVDSPPlan{T,N},
) where {T,N}
    return _execute_plan!(buffer, plan)
end

function LinearAlgebra.mul!(
    destination::Array{Complex{T},N},
    plan::AdaptiveOpticsVDSPPlan{T,N},
    source::Array{Complex{T},N},
) where {T,N}
    size(source) == size(plan) || throw(DimensionMismatch(
        "vDSP FFT plan size $(size(plan)) does not match source size $(size(source))"))
    destination === source || copyto!(destination, source)
    return _execute_plan!(destination, plan)
end

end # Sys.isapple()

end
