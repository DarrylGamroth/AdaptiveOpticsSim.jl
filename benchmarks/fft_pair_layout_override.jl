# Benchmark-only full FFT/mask/inverse pipeline. Include fft_layout_override.jl
# first, then this file, then the physical profile harness. Set
# AOS_FFT_LAYOUT=native AOS_FFT_PAIR=1. The temporary holds transposed frequency
# coordinates through the mask and first inverse transform, avoiding two of the
# four transposes in the independent forward/backward transpose experiment.
module FFTPairLayoutExperiment

using AdaptiveOpticsSim
import AdaptiveOpticsSim.Backends
import AdaptiveOpticsSim.WavefrontSensors as WFS
import CUDA
import ..FFTLayoutExperiment
import ..FFTTransposeExperiment
using KernelAbstractions: @kernel, @index
using LinearAlgebra: mul!, permutedims!

const ENABLED = get(ENV, "AOS_FFT_PAIR", "0") == "1"
ENABLED && FFTLayoutExperiment.LAYOUT != "native" &&
    error("AOS_FFT_PAIR=1 requires AOS_FFT_LAYOUT=native")

struct PairFFTPlan{P1,P2,P3,A}
    axis_1::P1
    transposed_axis_1::P2
    transposed_inverse_axis_1::P3
    temporary::A
end

@kernel function transposed_mask_kernel!(stack, mask, pad::Int, batch_size::Int)
    i, j = @index(Global, NTuple)
    if i <= pad && j <= pad
        mask_value = @inbounds mask[j, i]
        @inbounds for batch_index in 1:batch_size
            stack[i, j, batch_index] *= mask_value
        end
    end
end

if ENABLED
    function Backends.plan_repeated_fft_backend!(
        buffer::CUDA.CuArray{T,3}, dims::Tuple{Int,Int},
    ) where {T<:Union{ComplexF32,ComplexF64}}
        dims == (1, 2) || return Backends.plan_fft_backend!(buffer, dims)
        size(buffer, 1) == size(buffer, 2) || error("FFT pair benchmark requires square planes")
        temporary = similar(buffer)
        return PairFFTPlan(Backends.plan_fft_backend!(buffer, (1,)),
            Backends.plan_fft_backend!(temporary, (1,)),
            Backends.plan_bfft_backend!(temporary, (1,)), temporary)
    end

    function Backends.plan_repeated_bfft_backend!(
        buffer::CUDA.CuArray{T,3}, dims::Tuple{Int,Int},
    ) where {T<:Union{ComplexF32,ComplexF64}}
        return Backends.plan_bfft_backend!(buffer, dims == (1, 2) ? (1,) : dims)
    end
end

# This owner is a pipeline experiment, not a replacement FFT with an implicit
# output layout. Reject accidental use outside the explicit pair execution.
Backends.enqueue_fft_plan!(buffer, ::PairFFTPlan) =
    error("PairFFTPlan requires the complete pyramid FFT/mask/inverse pipeline")
Backends.execute_fft_plan!(buffer, ::PairFFTPlan) =
    error("PairFFTPlan requires the complete pyramid FFT/mask/inverse pipeline")

@inline function execute_pair!(batch, mask, style, launch)
    forward = batch.fft_plan
    stack = batch.field_stack
    temporary = forward.temporary
    mul!(stack, forward.axis_1, stack)
    FFTTransposeExperiment.transpose_planes!(temporary, stack)
    mul!(temporary, forward.transposed_axis_1, temporary)
    launch(style, transposed_mask_kernel!, temporary, mask,
        size(stack, 1), batch.batch_size; ndrange=size(mask))
    mul!(temporary, forward.transposed_inverse_axis_1, temporary)
    FFTTransposeExperiment.transpose_planes!(stack, temporary)
    mul!(stack, batch.bfft_plan, stack)
    return stack
end

@inline function WFS._propagate_pyramid_modulation_batch!(out,
    front_end::WFS.PyramidOpticalFrontEnd,
    batch::WFS.PyramidModulationBatchWorkspace{C,P,V,Pf,Pb},
    style::Backends.AcceleratorStyle,
) where {C,P,V,Pf<:PairFFTPlan,Pb}
    pad = size(out, 1)
    intensity_scale = WFS.pyramid_bfft_intensity_scale(eltype(out), pad)
    execute_pair!(batch, WFS.pyramid_propagation_workspace(front_end).pyramid_mask,
        style, WFS.launch_kernel!)
    WFS.launch_kernel!(style, WFS.pyramid_modulation_batch_intensity_kernel!,
        out, batch.field_stack, intensity_scale, pad, batch.batch_size;
        ndrange=size(out))
    return nothing
end

@inline function WFS._enqueue_pyramid_rate_batch!(
    style::Backends.AcceleratorStyle,
    out,
    front_end::WFS.PyramidOpticalFrontEnd,
    input::WFS.PupilFunction,
    source::WFS.Source,
    batch::WFS.PyramidModulationBatchWorkspace{C,P,V,Pf,Pb},
) where {C,P,V,Pf<:PairFFTPlan,Pb}
    propagation = WFS.pyramid_propagation_workspace(front_end)
    resolution = input.metadata.dimensions[1]
    pad = size(propagation.field, 1)
    offset = div(pad - resolution, 2)
    T = eltype(propagation.intensity)
    amplitude_scale = sqrt(T(WFS.photon_irradiance(source)) *
        T(input.metadata.sampling[1] * input.metadata.sampling[2]))
    opd_to_cycles = T(2) / T(WFS.wavelength(source))
    point_count = WFS.modulation_point_count(front_end.modulation)
    intensity_scale = WFS.pyramid_bfft_intensity_scale(T, pad)

    fill!(out, zero(T))
    WFS.launch_kernel_async!(style, WFS.pyramid_pupil_phase_kernel!, batch.pupil_phase,
        input.opd, opd_to_cycles, resolution; ndrange=size(batch.pupil_phase))
    phases = AdaptiveOpticsSim.Optics._modulation_phase_storage(front_end.modulation)
    @inbounds for first_point in 1:batch.batch_size:point_count
        WFS.launch_kernel_async!(style, WFS.pyramid_pupil_modulation_batch_kernel!,
            batch.field_stack, input.amplitude, batch.pupil_phase,
            phases, batch.operating_weights,
            propagation.phasor, amplitude_scale, offset,
            resolution, first_point, pad, batch.batch_size;
            ndrange=size(batch.field_stack))
        execute_pair!(batch, propagation.pyramid_mask, style, WFS.launch_kernel_async!)
        WFS.launch_kernel_async!(style,
            WFS.pyramid_modulation_batch_intensity_kernel!, out,
            batch.field_stack, intensity_scale, pad, batch.batch_size;
            ndrange=size(out))
    end
    return out
end

FFTLayoutExperiment.layout_name(::PairFFTPlan) = "pair"
FFTLayoutExperiment.scratch_bytes(plan::PairFFTPlan) =
    sizeof(eltype(plan.temporary)) * length(plan.temporary)

end
