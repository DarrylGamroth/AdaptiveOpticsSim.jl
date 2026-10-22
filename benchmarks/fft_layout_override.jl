# Benchmark-only FFT decomposition. Include before preparing any optical node,
# then include profile_physical_pyramid_cuda.jl and call physical_profile_main.
# AOS_FFT_LAYOUT=native|axis|transpose never changes production configuration.
include("fft_transpose_cuda.jl")
module FFTLayoutExperiment

using AdaptiveOpticsSim
import AdaptiveOpticsSim.Backends
import CUDA
import ..FFTTransposeExperiment
using LinearAlgebra: mul!, permutedims!
import AdaptiveOpticsSim.WavefrontSensors as WFS
function WFS._prepare_pyramid_modulation_fft(
    ::Backends.AcceleratorStyle{CUDA.CUDABackend},
    stack::CUDA.CuArray{ComplexF32,3}, points::Int)
    return WFS._prepare_pyramid_native_fft(stack)
end

const LAYOUT = get(ENV, "AOS_FFT_LAYOUT", "native")
LAYOUT in ("native", "axis", "transpose") ||
    error("AOS_FFT_LAYOUT must be native, axis, or transpose")

struct AxisFFTPlan{P1,P2}
    axis_1::P1
    axis_2::P2
end

struct TransposeFFTPlan{P1,P2,A}
    axis_1::P1
    transposed_axis_1::P2
    temporary::A
end

function prepare_layout_plan(buffer, dims, native_plan)
    (LAYOUT == "native" || dims != (1, 2)) && return native_plan(buffer, dims)
    axis_1 = native_plan(buffer, (1,))
    if LAYOUT == "axis"
        return AxisFFTPlan(axis_1, native_plan(buffer, (2,)))
    end
    temporary = similar(buffer, (size(buffer, 2), size(buffer, 1), size(buffer, 3)))
    return TransposeFFTPlan(axis_1, native_plan(temporary, (1,)), temporary)
end

function Backends.plan_repeated_fft_backend!(
    buffer::CUDA.CuArray{T,3}, dims::Tuple{Int,Int},
) where {T<:Union{ComplexF32,ComplexF64}}
    return prepare_layout_plan(buffer, dims, Backends.plan_fft_backend!)
end

function Backends.plan_repeated_bfft_backend!(
    buffer::CUDA.CuArray{T,3}, dims::Tuple{Int,Int},
) where {T<:Union{ComplexF32,ComplexF64}}
    # The inverse remains unnormalized, matching the production bfft path.
    return prepare_layout_plan(buffer, dims, Backends.plan_bfft_backend!)
end

@inline function Backends.enqueue_fft_plan!(buffer, plan::AxisFFTPlan)
    mul!(buffer, plan.axis_1, buffer)
    mul!(buffer, plan.axis_2, buffer)
    return buffer
end

@inline function Backends.enqueue_fft_plan!(buffer, plan::TransposeFFTPlan)
    mul!(buffer, plan.axis_1, buffer)
    FFTTransposeExperiment.transpose_planes!(plan.temporary, buffer)
    mul!(plan.temporary, plan.transposed_axis_1, plan.temporary)
    FFTTransposeExperiment.transpose_planes!(buffer, plan.temporary)
    return buffer
end

@inline Backends.execute_fft_plan!(buffer, plan::AxisFFTPlan) =
    Backends.enqueue_fft_plan!(buffer, plan)
@inline Backends.execute_fft_plan!(buffer, plan::TransposeFFTPlan) =
    Backends.enqueue_fft_plan!(buffer, plan)

# Cold reporting only: cuFFT's internal work area is not included. Forward and
# backward transpose plans own separate arrays, so report and add both plans.
scratch_bytes(plan) = 0
scratch_bytes(plan::TransposeFFTPlan) = sizeof(eltype(plan.temporary)) * length(plan.temporary)

layout_name(plan) = "native"
layout_name(plan::AxisFFTPlan) = "axis"
layout_name(plan::TransposeFFTPlan) = "transpose"

end
