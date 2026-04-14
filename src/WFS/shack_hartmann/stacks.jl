function compute_intensity!(wfs::ShackHartmann, tel::Telescope, src::AbstractSource,
    xs::Int, ys::Int, xe::Int, ye::Int, ox::Int, oy::Int, sub::Int)
    opd_to_cycles = eltype(wfs.state.intensity)(2) / wavelength(src)
    amp_scale = sqrt(eltype(wfs.state.intensity)(photon_flux(src) * tel.params.sampling_time *
        (tel.params.diameter / tel.params.resolution)^2))
    fill!(wfs.state.field, zero(eltype(wfs.state.field)))
    @views @. wfs.state.field[ox+1:ox+sub, oy+1:oy+sub] =
        amp_scale * tel.state.pupil[xs:xe, ys:ye] * cispi(opd_to_cycles * tel.state.opd[xs:xe, ys:ye])
    @. wfs.state.field *= wfs.state.phasor
    copyto!(wfs.state.fft_buffer, wfs.state.field)
    execute_fft_plan!(wfs.state.fft_buffer, wfs.state.fft_plan)
    @. wfs.state.intensity = abs2(wfs.state.fft_buffer)
    return wfs.state.intensity
end

@inline function compute_intensity_safe!(style::ExecutionStyle, wfs::ShackHartmann, tel::Telescope, src::AbstractSource,
    xs::Int, ys::Int, xe::Int, ye::Int, ox::Int, oy::Int, sub::Int)
    return compute_intensity!(wfs, tel, src, xs, ys, xe, ye, ox, oy, sub)
end

@kernel function sh_field_stack_kernel!(fft_stack, valid_mask, pupil, opd, phasor,
    amp_scale, opd_to_cycles, n_sub::Int, sub::Int, ox::Int, oy::Int, n::Int, pad::Int)
    x, y, idx = @index(Global, NTuple)
    if x <= pad && y <= pad && idx <= n_sub * n_sub
        i = (idx - 1) ÷ n_sub + 1
        j = idx - (i - 1) * n_sub
        T = eltype(phasor)
        val = zero(T)
        if @inbounds valid_mask[i, j]
            if ox + 1 <= x <= ox + sub && oy + 1 <= y <= oy + sub
                px = (i - 1) * sub + (x - ox)
                py = (j - 1) * sub + (y - oy)
                if px <= n && py <= n
                    amp = amp_scale * pupil[px, py]
                    val = amp * cispi(opd_to_cycles * opd[px, py]) * phasor[x, y]
                end
            end
        end
        @inbounds fft_stack[x, y, idx] = val
    end
end

@kernel function complex_abs2_stack_kernel!(intensity_stack, fft_stack, pad::Int, n_spots::Int)
    x, y, idx = @index(Global, NTuple)
    if x <= pad && y <= pad && idx <= n_spots
        @inbounds intensity_stack[x, y, idx] = abs2(fft_stack[x, y, idx])
    end
end

@kernel function sh_field_asterism_stack_kernel!(fft_stack, valid_mask, pupil, opd, phasor,
    amp_scales, opd_to_cycles, n_sub::Int, sub::Int, ox::Int, oy::Int, n::Int, pad::Int, n_spots::Int)
    x, y, idx_total = @index(Global, NTuple)
    if x <= pad && y <= pad && idx_total <= size(fft_stack, 3)
        src_idx = (idx_total - 1) ÷ n_spots + 1
        idx = idx_total - (src_idx - 1) * n_spots
        i = (idx - 1) ÷ n_sub + 1
        j = idx - (i - 1) * n_sub
        T = eltype(phasor)
        val = zero(T)
        if @inbounds valid_mask[i, j]
            if ox + 1 <= x <= ox + sub && oy + 1 <= y <= oy + sub
                px = (i - 1) * sub + (x - ox)
                py = (j - 1) * sub + (y - oy)
                if px <= n && py <= n
                    amp = (@inbounds amp_scales[src_idx]) * pupil[px, py]
                    val = amp * cispi((@inbounds opd_to_cycles[src_idx]) * opd[px, py]) * phasor[x, y]
                end
            end
        end
        @inbounds fft_stack[x, y, idx_total] = val
    end
end

@kernel function sh_sample_spot_stack_kernel!(spot_cube, intensity_stack, valid_mask,
    binning::Int, n_sub::Int, n_binned::Int, n_out::Int)
    idx, u, v = @index(Global, NTuple)
    n_spots = n_sub * n_sub
    if idx <= n_spots && u <= n_out && v <= n_out
        i = (idx - 1) ÷ n_sub + 1
        j = idx - (i - 1) * n_sub
        T = eltype(spot_cube)
        val = zero(T)
        if @inbounds valid_mask[i, j]
            ox = div(n_out - n_binned, 2)
            oy = div(n_out - n_binned, 2)
            bu = u - ox
            bv = v - oy
            if 1 <= bu <= n_binned && 1 <= bv <= n_binned
                if binning == 1
                    @inbounds val = intensity_stack[bu, bv, idx]
                else
                    acc = zero(T)
                    @inbounds for ii in 1:binning, jj in 1:binning
                        acc += intensity_stack[(bu - 1) * binning + ii, (bv - 1) * binning + jj, idx]
                    end
                    val = acc
                end
            end
        end
        @inbounds spot_cube[idx, u, v] = val
    end
end

function compute_intensity_stack!(style::AcceleratorStyle, wfs::ShackHartmann, tel::Telescope, src::AbstractSource)
    pad = size(wfs.state.fft_stack, 1)
    n_sub = wfs.params.n_subap
    sub = div(tel.params.resolution, n_sub)
    ox = div(pad - sub, 2)
    oy = div(pad - sub, 2)
    T = eltype(wfs.state.intensity)
    opd_to_cycles = T(2) / wavelength(src)
    amp_scale = sqrt(T(photon_flux(src) * tel.params.sampling_time *
        (tel.params.diameter / tel.params.resolution)^2))
    launch_kernel_async!(style, sh_field_stack_kernel!, wfs.state.fft_stack, wfs.state.valid_mask,
        tel.state.pupil, tel.state.opd, wfs.state.phasor, amp_scale, opd_to_cycles, n_sub, sub, ox, oy,
        tel.params.resolution, pad; ndrange=size(wfs.state.fft_stack))
    synchronize_backend!(style)
    execute_fft_plan!(wfs.state.fft_stack, wfs.state.fft_stack_plan)
    launch_kernel!(style, complex_abs2_stack_kernel!, wfs.state.intensity_stack, wfs.state.fft_stack,
        pad, n_sub * n_sub; ndrange=size(wfs.state.intensity_stack))
    return wfs.state.intensity_stack
end

@inline sh_stacked_asterism_compatible(ast::Asterism) =
    !isempty(ast.sources) && all(src -> !(src isa LGSSource), ast.sources)

function compute_intensity_asterism_stack!(style::AcceleratorStyle, wfs::ShackHartmann, tel::Telescope, ast::Asterism)
    n_src = length(ast.sources)
    ensure_sh_asterism_buffers!(wfs, n_src)
    pad = size(wfs.state.fft_stack, 1)
    n_sub = wfs.params.n_subap
    n_spots = n_sub * n_sub
    sub = div(tel.params.resolution, n_sub)
    ox = div(pad - sub, 2)
    oy = div(pad - sub, 2)
    T = eltype(wfs.state.intensity)
    amp_scales = wfs.state.amp_scales
    host_amp_scales = wfs.state.amp_scales_host
    opd_to_cycles = wfs.state.opd_to_cycles
    host_opd_to_cycles = wfs.state.opd_to_cycles_host
    @inbounds for i in eachindex(ast.sources)
        src = ast.sources[i]
        host_amp_scales[i] = sqrt(T(photon_flux(src) * tel.params.sampling_time *
            (tel.params.diameter / tel.params.resolution)^2))
        host_opd_to_cycles[i] = T(2) / wavelength(src)
    end
    copyto!(amp_scales, host_amp_scales)
    copyto!(opd_to_cycles, host_opd_to_cycles)
    total = n_spots * n_src
    fft_view = @view wfs.state.fft_asterism_stack[:, :, 1:total]
    intensity_view = @view wfs.state.intensity_tmp_stack[:, :, 1:total]
    launch_kernel_async!(style, sh_field_asterism_stack_kernel!, fft_view, wfs.state.valid_mask,
        tel.state.pupil, tel.state.opd, wfs.state.phasor, amp_scales, opd_to_cycles,
        n_sub, sub, ox, oy, tel.params.resolution, pad, n_spots; ndrange=size(fft_view))
    synchronize_backend!(style)
    execute_fft_plan!(fft_view, wfs.state.fft_asterism_plan)
    phase = begin_kernel_phase(style)
    queue_kernel!(phase, complex_abs2_stack_kernel!, intensity_view, fft_view, pad, total; ndrange=size(intensity_view))
    queue_kernel!(phase, reduce_grouped_blocks_kernel!, wfs.state.intensity_stack, intensity_view,
        n_spots, n_src, size(wfs.state.intensity_stack, 1), size(wfs.state.intensity_stack, 2); ndrange=size(wfs.state.intensity_stack))
    finish_kernel_phase!(phase)
    return wfs.state.intensity_stack
end

function compute_intensity_spectral_stack!(style::AcceleratorStyle, wfs::ShackHartmann, tel::Telescope, src::SpectralSource)
    bundle = spectral_bundle(src)
    n_src = length(bundle)
    ensure_sh_asterism_buffers!(wfs, n_src)
    pad = size(wfs.state.fft_stack, 1)
    n_sub = wfs.params.n_subap
    n_spots = n_sub * n_sub
    sub = div(tel.params.resolution, n_sub)
    ox = div(pad - sub, 2)
    oy = div(pad - sub, 2)
    T = eltype(wfs.state.intensity)
    amp_scales = wfs.state.amp_scales
    host_amp_scales = wfs.state.amp_scales_host
    opd_to_cycles = wfs.state.opd_to_cycles
    host_opd_to_cycles = wfs.state.opd_to_cycles_host
    base_flux = T(photon_flux(src) * tel.params.sampling_time * (tel.params.diameter / tel.params.resolution)^2)
    @inbounds for i in eachindex(bundle.samples)
        sample = bundle.samples[i]
        host_amp_scales[i] = sqrt(base_flux * sample.weight)
        host_opd_to_cycles[i] = T(2) / sample.wavelength
    end
    copyto!(amp_scales, host_amp_scales)
    copyto!(opd_to_cycles, host_opd_to_cycles)
    total = n_spots * n_src
    fft_view = @view wfs.state.fft_asterism_stack[:, :, 1:total]
    intensity_view = @view wfs.state.intensity_tmp_stack[:, :, 1:total]
    launch_kernel_async!(style, sh_field_asterism_stack_kernel!, fft_view, wfs.state.valid_mask,
        tel.state.pupil, tel.state.opd, wfs.state.phasor, amp_scales, opd_to_cycles,
        n_sub, sub, ox, oy, tel.params.resolution, pad, n_spots; ndrange=size(fft_view))
    synchronize_backend!(style)
    execute_fft_plan!(fft_view, wfs.state.fft_asterism_plan)
    phase = begin_kernel_phase(style)
    queue_kernel!(phase, complex_abs2_stack_kernel!, intensity_view, fft_view, pad, total; ndrange=size(intensity_view))
    queue_kernel!(phase, reduce_grouped_blocks_kernel!, wfs.state.intensity_stack, intensity_view,
        n_spots, n_src, size(wfs.state.intensity_stack, 1), size(wfs.state.intensity_stack, 2); ndrange=size(wfs.state.intensity_stack))
    finish_kernel_phase!(phase)
    return wfs.state.intensity_stack
end

function sample_spot_stack!(::ScalarCPUStyle, wfs::ShackHartmann)
    n_spots = size(wfs.state.sampled_spot_cube, 1)
    @inbounds for idx in 1:n_spots
        sample_spot!(wfs, @view(wfs.state.intensity_stack[:, :, idx]))
        copyto!(@view(wfs.state.sampled_spot_cube[idx, :, :]), wfs.state.spot)
    end
    return wfs.state.sampled_spot_cube
end

function sample_spot_stack!(style::AcceleratorStyle, wfs::ShackHartmann)
    pad = size(wfs.state.intensity_stack, 1)
    binning = wfs.state.binning_pixel_scale
    @assert pad % binning == 0
    n_binned = div(pad, binning)
    n_out = size(wfs.state.sampled_spot_cube, 2)
    launch_kernel!(style, sh_sample_spot_stack_kernel!, wfs.state.sampled_spot_cube, wfs.state.intensity_stack,
        wfs.state.valid_mask, binning, wfs.params.n_subap, n_binned, n_out; ndrange=size(wfs.state.sampled_spot_cube))
    return wfs.state.sampled_spot_cube
end

function sampled_spots_peak_asterism_stacked!(::ScalarCPUStyle, wfs::ShackHartmann, tel::Telescope, ast::Asterism)
    fill!(wfs.state.detector_noise_cube, zero(eltype(wfs.state.detector_noise_cube)))
    for src in ast.sources
        sampled_spots_peak!(ScalarCPUStyle(), wfs, tel, src)
        wfs.state.detector_noise_cube .+= wfs.state.spot_cube
    end
    copyto!(wfs.state.spot_cube, wfs.state.detector_noise_cube)
    return sh_safe_peak_value(wfs.state.spot_cube)
end

function sampled_spots_peak_asterism_stacked!(::ScalarCPUStyle, wfs::ShackHartmann, tel::Telescope, ast::Asterism,
    det::AbstractDetector, rng::AbstractRNG)
    fill!(wfs.state.detector_noise_cube, zero(eltype(wfs.state.detector_noise_cube)))
    for src in ast.sources
        sampled_spots_peak!(ScalarCPUStyle(), wfs, tel, src, det, rng)
        wfs.state.detector_noise_cube .+= wfs.state.spot_cube
    end
    copyto!(wfs.state.spot_cube, wfs.state.detector_noise_cube)
    return sh_safe_peak_value(wfs.state.spot_cube)
end

function sampled_spots_peak_asterism_stacked!(style::AcceleratorStyle, wfs::ShackHartmann, tel::Telescope, ast::Asterism)
    compute_intensity_asterism_stack!(style, wfs, tel, ast)
    sample_spot_stack!(style, wfs)
    sync_signal_spots_from_sampled!(wfs)
    return sh_safe_peak_value(wfs.state.spot_cube)
end

function sampled_spots_peak_asterism_stacked!(style::AcceleratorStyle, wfs::ShackHartmann, tel::Telescope, ast::Asterism,
    det::AbstractDetector, rng::AbstractRNG)
    compute_intensity_asterism_stack!(style, wfs, tel, ast)
    sample_spot_stack!(style, wfs)
    n_sub = wfs.params.n_subap
    capture_sampled_spot_stack!(wfs, det, rng)
    launch_kernel!(style, zero_invalid_spots_kernel!, wfs.state.spot_cube, wfs.state.valid_mask,
        n_sub, size(wfs.state.spot_cube, 2), size(wfs.state.spot_cube, 3); ndrange=size(wfs.state.spot_cube))
    return sh_safe_peak_value(wfs.state.spot_cube)
end

@kernel function sh_spot_centroid_stats_kernel!(stats, spot_cube, valid_mask, threshold, n_sub::Int, n1::Int, n2::Int)
    idx = @index(Global, Linear)
    n_spots = n_sub * n_sub
    if idx <= n_spots
        i = (idx - 1) ÷ n_sub + 1
        j = idx - (i - 1) * n_sub
        base = 3 * (idx - 1)
        T = eltype(stats)
        peak = zero(T)
        total = zero(T)
        sx = zero(T)
        sy = zero(T)
        if @inbounds valid_mask[i, j]
            @inbounds for x in 1:n1, y in 1:n2
                val = spot_cube[idx, x, y]
                peak = max(peak, val)
            end
            cutoff = threshold * peak
            @inbounds for x in 1:n1, y in 1:n2
                val = spot_cube[idx, x, y]
                if val >= cutoff
                    total += val
                    sx += T(x - 1) * val
                    sy += T(y - 1) * val
                end
            end
            if total > 0
                sx /= total
                sy /= total
            end
        end
        @inbounds begin
            stats[base + 1] = total
            stats[base + 2] = sx
            stats[base + 3] = sy
        end
    end
end

@kernel function accumulate_spot_stats_kernel!(accum, stats, n_spots::Int)
    idx = @index(Global, Linear)
    if idx <= 3 * n_spots
        @inbounds accum[idx] += stats[idx]
    end
end

@kernel function sh_finalize_asterism_slopes_kernel!(slopes, stats_accum, reference, valid_mask,
    slopes_units, n_src::Int, n_sub::Int, offset::Int)
    idx = @index(Global, Linear)
    n_spots = n_sub * n_sub
    if idx <= n_spots
        i = (idx - 1) ÷ n_sub + 1
        j = idx - (i - 1) * n_sub
        base = 3 * (idx - 1)
        T = eltype(slopes)
        if @inbounds valid_mask[i, j]
            total_sum = stats_accum[base + 1]
            if total_sum <= 0
                @inbounds begin
                    slopes[idx] = zero(T)
                    slopes[idx + offset] = zero(T)
                end
            else
                sy_sum = stats_accum[base + 3]
                sx_sum = stats_accum[base + 2]
                inv_nsrc = inv(T(n_src))
                @inbounds begin
                    slopes[idx] = ((sy_sum * inv_nsrc) - reference[idx]) / slopes_units
                    slopes[idx + offset] = ((sx_sum * inv_nsrc) - reference[idx + offset]) / slopes_units
                end
            end
        else
            @inbounds begin
                slopes[idx] = zero(T)
                slopes[idx + offset] = zero(T)
            end
        end
    end
end

@inline function _sampled_spots_peak_source_batched!(style::AcceleratorStyle, wfs::ShackHartmann, tel::Telescope, src::AbstractSource)
    return sampled_spots_peak!(style, wfs, tel, src)
end

@inline function _sampled_spots_peak_source_batched!(style::AcceleratorStyle, wfs::ShackHartmann, tel::Telescope, src::LGSSource)
    return sampled_spots_peak!(style, wfs, tel, src)
end

@inline function _sampled_spots_peak_source_batched!(style::AcceleratorStyle, wfs::ShackHartmann, tel::Telescope,
    src::AbstractSource, det::AbstractDetector, rng::AbstractRNG)
    return sampled_spots_peak!(style, wfs, tel, src, det, rng)
end

@inline function _sampled_spots_peak_source_batched!(style::AcceleratorStyle, wfs::ShackHartmann, tel::Telescope,
    src::LGSSource, det::AbstractDetector, rng::AbstractRNG)
    return sampled_spots_peak!(style, wfs, tel, src, det, rng)
end

function measure_sh_asterism_batched!(style::AcceleratorStyle, wfs::ShackHartmann, tel::Telescope, ast::Asterism)
    n_src = length(ast.sources)
    n_sub = wfs.params.n_subap
    n_spots = n_sub * n_sub
    fill!(wfs.state.spot_cube_accum, zero(eltype(wfs.state.spot_cube_accum)))
    fill!(wfs.state.spot_stats_accum, zero(eltype(wfs.state.spot_stats_accum)))
    for src in ast.sources
        _sampled_spots_peak_source_batched!(style, wfs, tel, src)
        wfs.state.spot_cube_accum .+= wfs.state.spot_cube
        phase = begin_kernel_phase(style)
        queue_kernel!(phase, sh_spot_centroid_stats_kernel!, wfs.state.spot_stats, wfs.state.spot_cube,
            wfs.state.valid_mask, centroid_threshold(wfs), n_sub, size(wfs.state.spot_cube, 2),
            size(wfs.state.spot_cube, 3); ndrange=n_spots)
        queue_kernel!(phase, accumulate_spot_stats_kernel!, wfs.state.spot_stats_accum,
            wfs.state.spot_stats, n_spots; ndrange=3 * n_spots)
        finish_kernel_phase!(phase)
    end
    copyto!(wfs.state.spot_cube, wfs.state.spot_cube_accum)
    sync_exported_spots!(wfs)
    offset = n_spots
    reference = vec(wfs.state.reference_signal_2d)
    launch_kernel!(style, sh_finalize_asterism_slopes_kernel!, wfs.state.slopes, wfs.state.spot_stats_accum,
        reference, wfs.state.valid_mask, wfs.state.slopes_units, n_src, n_sub, offset; ndrange=n_spots)
    return wfs.state.slopes
end

function measure_sh_asterism_batched!(style::AcceleratorStyle, wfs::ShackHartmann, tel::Telescope,
    ast::Asterism, det::AbstractDetector, rng::AbstractRNG)
    n_src = length(ast.sources)
    n_sub = wfs.params.n_subap
    n_spots = n_sub * n_sub
    fill!(wfs.state.spot_cube_accum, zero(eltype(wfs.state.spot_cube_accum)))
    fill!(wfs.state.spot_stats_accum, zero(eltype(wfs.state.spot_stats_accum)))
    for src in ast.sources
        _sampled_spots_peak_source_batched!(style, wfs, tel, src, det, rng)
        wfs.state.spot_cube_accum .+= wfs.state.spot_cube
        phase = begin_kernel_phase(style)
        queue_kernel!(phase, sh_spot_centroid_stats_kernel!, wfs.state.spot_stats, wfs.state.spot_cube,
            wfs.state.valid_mask, centroid_threshold(wfs), n_sub, size(wfs.state.spot_cube, 2),
            size(wfs.state.spot_cube, 3); ndrange=n_spots)
        queue_kernel!(phase, accumulate_spot_stats_kernel!, wfs.state.spot_stats_accum,
            wfs.state.spot_stats, n_spots; ndrange=3 * n_spots)
        finish_kernel_phase!(phase)
    end
    copyto!(wfs.state.spot_cube, wfs.state.spot_cube_accum)
    sync_exported_spots!(wfs)
    offset = n_spots
    reference = vec(wfs.state.reference_signal_2d)
    launch_kernel!(style, sh_finalize_asterism_slopes_kernel!, wfs.state.slopes, wfs.state.spot_stats_accum,
        reference, wfs.state.valid_mask, wfs.state.slopes_units, n_src, n_sub, offset; ndrange=n_spots)
    return wfs.state.slopes
end

@kernel function sh_spot_centroid_kernel!(slopes, spot_cube, valid_mask, cutoff, n_sub::Int, offset::Int, n1::Int, n2::Int)
    idx = @index(Global, Linear)
    n_spots = n_sub * n_sub
    if idx <= n_spots
        i = (idx - 1) ÷ n_sub + 1
        j = idx - (i - 1) * n_sub
        T = eltype(slopes)
        total = zero(T)
        sx = zero(T)
        sy = zero(T)
        if @inbounds valid_mask[i, j]
            @inbounds for x in 1:n1, y in 1:n2
                val = spot_cube[idx, x, y]
                if val < cutoff
                    spot_cube[idx, x, y] = zero(T)
                else
                    total += val
                    sx += T(x - 1) * val
                    sy += T(y - 1) * val
                end
            end
            if total <= 0
                slopes[idx] = zero(T)
                slopes[idx + offset] = zero(T)
            else
                slopes[idx] = sy / total
                slopes[idx + offset] = sx / total
            end
        else
            @inbounds begin
                slopes[idx] = zero(T)
                slopes[idx + offset] = zero(T)
            end
        end
    end
end

@kernel function sh_spot_cutoff_stats_kernel!(stats, spot_cube, valid_mask, cutoff, n_sub::Int, n1::Int, n2::Int)
    idx = @index(Global, Linear)
    n_spots = n_sub * n_sub
    if idx <= n_spots
        i = (idx - 1) ÷ n_sub + 1
        j = idx - (i - 1) * n_sub
        base = 3 * (idx - 1)
        T = eltype(stats)
        total = zero(T)
        sx = zero(T)
        sy = zero(T)
        if @inbounds valid_mask[i, j]
            @inbounds for x in 1:n1, y in 1:n2
                val = spot_cube[idx, x, y]
                if val >= cutoff
                    total += val
                    sx += T(x - 1) * val
                    sy += T(y - 1) * val
                end
            end
        end
        @inbounds begin
            stats[base + 1] = total
            stats[base + 2] = sx
            stats[base + 3] = sy
        end
    end
end

@kernel function sh_finalize_spot_slopes_kernel!(slopes, stats, valid_mask, n_sub::Int, offset::Int)
    idx = @index(Global, Linear)
    n_spots = n_sub * n_sub
    if idx <= n_spots
        i = (idx - 1) ÷ n_sub + 1
        j = idx - (i - 1) * n_sub
        base = 3 * (idx - 1)
        T = eltype(slopes)
        if @inbounds valid_mask[i, j]
            total = @inbounds stats[base + 1]
            if total > zero(T)
                sx = @inbounds stats[base + 2] / total
                sy = @inbounds stats[base + 3] / total
                @inbounds begin
                    slopes[idx] = sy
                    slopes[idx + offset] = sx
                end
            else
                @inbounds begin
                    slopes[idx] = zero(T)
                    slopes[idx + offset] = zero(T)
                end
            end
        else
            @inbounds begin
                slopes[idx] = zero(T)
                slopes[idx + offset] = zero(T)
            end
        end
    end
end

@kernel function zero_invalid_spots_kernel!(spot_cube, valid_mask, n_sub::Int, n1::Int, n2::Int)
    idx, x, y = @index(Global, NTuple)
    n_spots = n_sub * n_sub
    if idx <= n_spots && x <= n1 && y <= n2
        i = (idx - 1) ÷ n_sub + 1
        j = idx - (i - 1) * n_sub
        if !(@inbounds valid_mask[i, j])
            @inbounds spot_cube[idx, x, y] = zero(eltype(spot_cube))
        end
    end
end

@kernel function zero_invalid_sh_slopes_kernel!(slopes, valid_mask, n_sub::Int, offset::Int)
    idx = @index(Global, Linear)
    n_spots = n_sub * n_sub
    if idx <= n_spots
        i = (idx - 1) ÷ n_sub + 1
        j = idx - (i - 1) * n_sub
        if !(@inbounds valid_mask[i, j])
            T = eltype(slopes)
            @inbounds begin
                slopes[idx] = zero(T)
                slopes[idx + offset] = zero(T)
            end
        end
    end
end
