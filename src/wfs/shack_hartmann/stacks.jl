@inline sh_fft_intensity_scale(::Type{T}, pad::Integer) where {T} =
    inv(T(pad) * T(pad))

function compute_intensity!(wfs::ShackHartmannWFS, tel::Telescope, src::AbstractSource,
    xs::Int, ys::Int, xe::Int, ye::Int, ox::Int, oy::Int, sub::Int)
    opd_to_cycles = eltype(wfs.front_end.propagation.intensity)(2) / wavelength(src)
    amp_scale = sqrt(eltype(wfs.front_end.propagation.intensity)(photon_irradiance(src) *
        (tel.params.diameter / tel.params.resolution)^2))
    fill!(wfs.front_end.propagation.field, zero(eltype(wfs.front_end.propagation.field)))
    reflectivity = pupil_reflectivity(tel)
    @views @. wfs.front_end.propagation.field[ox+1:ox+sub, oy+1:oy+sub] =
        amp_scale * sqrt(reflectivity[xs:xe, ys:ye]) *
        cispi(opd_to_cycles * tel.state.opd[xs:xe, ys:ye])
    @. wfs.front_end.propagation.field *= wfs.front_end.propagation.phasor
    copyto!(wfs.front_end.propagation.fft_buffer, wfs.front_end.propagation.field)
    execute_fft_plan!(wfs.front_end.propagation.fft_buffer, wfs.front_end.propagation.fft_plan)
    intensity_scale = sh_fft_intensity_scale(eltype(wfs.front_end.propagation.intensity),
        size(wfs.front_end.propagation.fft_buffer, 1))
    @. wfs.front_end.propagation.intensity = abs2(wfs.front_end.propagation.fft_buffer) * intensity_scale
    return wfs.front_end.propagation.intensity
end

@inline function compute_intensity_safe!(style::ExecutionStyle, wfs::ShackHartmannWFS, tel::Telescope, src::AbstractSource,
    xs::Int, ys::Int, xe::Int, ye::Int, ox::Int, oy::Int, sub::Int)
    return compute_intensity!(wfs, tel, src, xs, ys, xe, ye, ox, oy, sub)
end

@kernel function sh_field_stack_kernel!(fft_stack, valid_mask, reflectivity, opd, phasor,
    amp_scale, opd_to_cycles, n_sub::Int, sub::Int, ox::Int, oy::Int, n::Int, pad::Int)
    x, y, i, j = @index(Global, NTuple)
    if x <= pad && y <= pad && i <= n_sub && j <= n_sub
        idx = sh_lenslet_index(i, j, n_sub)
        xi = x - ox
        yi = y - oy
        val = zero(eltype(fft_stack))
        if 1 <= xi <= sub && 1 <= yi <= sub && @inbounds(valid_mask[i, j])
            px = (i - 1) * sub + xi
            py = (j - 1) * sub + yi
            if px <= n && py <= n
                @inbounds val = amp_scale * sqrt(reflectivity[px, py]) *
                    cispi(opd_to_cycles * opd[px, py]) * phasor[x, y]
            end
        end
        @inbounds fft_stack[x, y, idx] = val
    end
end

@kernel function complex_abs2_stack_kernel!(intensity_stack, fft_stack,
    intensity_scale, pad::Int, n_spots::Int)
    x, y, idx = @index(Global, NTuple)
    if x <= pad && y <= pad && idx <= n_spots
        @inbounds intensity_stack[x, y, idx] =
            abs2(fft_stack[x, y, idx]) * intensity_scale
    end
end

@kernel function sh_field_asterism_stack_kernel!(fft_stack, valid_mask, reflectivity, opd, phasor,
    amp_scales, opd_to_cycles, n_sub::Int, sub::Int, ox::Int, oy::Int, n::Int,
    pad::Int, n_spots::Int, n_src::Int)
    x, y, src_idx, i, j = @index(Global, NTuple)
    if x <= pad && y <= pad && src_idx <= n_src && i <= n_sub && j <= n_sub
        idx = sh_lenslet_index(i, j, n_sub)
        idx_total = (src_idx - 1) * n_spots + idx
        xi = x - ox
        yi = y - oy
        val = zero(eltype(fft_stack))
        if 1 <= xi <= sub && 1 <= yi <= sub && @inbounds(valid_mask[i, j])
            px = (i - 1) * sub + xi
            py = (j - 1) * sub + yi
            if px <= n && py <= n
                @inbounds val = amp_scales[src_idx] * sqrt(reflectivity[px, py]) *
                    cispi(opd_to_cycles[src_idx] * opd[px, py]) * phasor[x, y]
            end
        end
        @inbounds fft_stack[x, y, idx_total] = val
    end
end

@kernel function sh_sample_spot_stack_kernel!(spot_cube, intensity_stack, valid_mask,
    binning::Int, n_sub::Int, n_binned::Int, n_out::Int, ox::Int, oy::Int)
    i, j, u, v = @index(Global, NTuple)
    if i <= n_sub && j <= n_sub && u <= n_out && v <= n_out
        idx = sh_lenslet_index(i, j, n_sub)
        T = eltype(spot_cube)
        val = zero(T)
        if @inbounds valid_mask[i, j]
            bu = u - ox
            bv = v - oy
            if 1 <= bu <= n_binned && 1 <= bv <= n_binned
                if binning == 1
                    @inbounds val = intensity_stack[bu, bv, idx]
                else
                    acc = zero(T)
                    @inbounds for jj in 1:binning, ii in 1:binning
                        acc += intensity_stack[(bu - 1) * binning + ii, (bv - 1) * binning + jj, idx]
                    end
                    val = acc
                end
            end
        end
        @inbounds spot_cube[idx, u, v] = val
    end
end

function compute_intensity_stack!(::ScalarCPUStyle, wfs::ShackHartmannWFS, tel::Telescope, src::AbstractSource)
    pad = size(wfs.front_end.propagation.fft_stack, 1)
    n_sub = n_lenslets(wfs)
    sub = div(tel.params.resolution, n_sub)
    ox = div(pad - sub, 2)
    oy = div(pad - sub, 2)
    T = eltype(wfs.front_end.propagation.intensity)
    opd_to_cycles = T(2) / wavelength(src)
    amp_scale = sqrt(T(photon_irradiance(src) *
        (tel.params.diameter / tel.params.resolution)^2))
    reflectivity = pupil_reflectivity(tel)
    fill!(wfs.front_end.propagation.fft_stack, zero(eltype(wfs.front_end.propagation.fft_stack)))
    idx = 1
    @inbounds for j in 1:n_sub, i in 1:n_sub
        if wfs.front_end.layout.valid_mask[i, j]
            for y in 1:sub, x in 1:sub
                px = (i - 1) * sub + x
                py = (j - 1) * sub + y
                wfs.front_end.propagation.fft_stack[ox + x, oy + y, idx] = amp_scale * sqrt(reflectivity[px, py]) *
                    cispi(opd_to_cycles * tel.state.opd[px, py]) * wfs.front_end.propagation.phasor[ox + x, oy + y]
            end
        end
        idx += 1
    end
    execute_fft_plan!(wfs.front_end.propagation.fft_stack, wfs.front_end.propagation.fft_stack_plan)
    intensity_scale = sh_fft_intensity_scale(T, pad)
    @inbounds for idx in axes(wfs.front_end.propagation.fft_stack, 3), y in axes(wfs.front_end.propagation.fft_stack, 2), x in axes(wfs.front_end.propagation.fft_stack, 1)
        wfs.front_end.propagation.intensity_stack[x, y, idx] =
            abs2(wfs.front_end.propagation.fft_stack[x, y, idx]) * intensity_scale
    end
    return wfs.front_end.propagation.intensity_stack
end

function compute_intensity_stack!(style::AcceleratorStyle, wfs::ShackHartmannWFS, tel::Telescope, src::AbstractSource)
    pad = size(wfs.front_end.propagation.fft_stack, 1)
    n_sub = n_lenslets(wfs)
    sub = div(tel.params.resolution, n_sub)
    ox = div(pad - sub, 2)
    oy = div(pad - sub, 2)
    T = eltype(wfs.front_end.propagation.intensity)
    opd_to_cycles = T(2) / wavelength(src)
    amp_scale = sqrt(T(photon_irradiance(src) *
        (tel.params.diameter / tel.params.resolution)^2))
    launch_kernel_async!(style, sh_field_stack_kernel!, wfs.front_end.propagation.fft_stack, wfs.front_end.layout.valid_mask,
        pupil_reflectivity(tel), tel.state.opd, wfs.front_end.propagation.phasor, amp_scale, opd_to_cycles, n_sub, sub, ox, oy,
        tel.params.resolution, pad; ndrange=(pad, pad, n_sub, n_sub))
    synchronize_backend!(style)
    execute_fft_plan!(wfs.front_end.propagation.fft_stack, wfs.front_end.propagation.fft_stack_plan)
    intensity_scale = sh_fft_intensity_scale(T, pad)
    launch_kernel!(style, complex_abs2_stack_kernel!, wfs.front_end.propagation.intensity_stack, wfs.front_end.propagation.fft_stack,
        intensity_scale, pad, n_sub * n_sub;
        ndrange=size(wfs.front_end.propagation.intensity_stack))
    return wfs.front_end.propagation.intensity_stack
end

@inline sh_stacked_asterism_compatible(ast::Asterism) =
    !isempty(ast.sources) && all(!is_lgs_source, ast.sources)

function compute_intensity_asterism_stack!(style::AcceleratorStyle, wfs::ShackHartmannWFS, tel::Telescope, ast::Asterism)
    require_sh_asterism_common_wavelength(ast)
    n_src = length(ast.sources)
    ensure_sh_asterism_buffers!(wfs.front_end, n_src)
    pad = size(wfs.front_end.propagation.fft_stack, 1)
    n_sub = n_lenslets(wfs)
    n_spots = n_sub * n_sub
    sub = div(tel.params.resolution, n_sub)
    ox = div(pad - sub, 2)
    oy = div(pad - sub, 2)
    T = eltype(wfs.front_end.propagation.intensity)
    amp_scales = wfs.front_end.propagation.amp_scales
    host_amp_scales = wfs.front_end.propagation.amp_scales_host
    opd_to_cycles = wfs.front_end.propagation.opd_to_cycles
    host_opd_to_cycles = wfs.front_end.propagation.opd_to_cycles_host
    @inbounds for i in eachindex(ast.sources)
        src = ast.sources[i]
        host_amp_scales[i] = sqrt(T(photon_irradiance(src) *
            (tel.params.diameter / tel.params.resolution)^2))
        host_opd_to_cycles[i] = T(2) / wavelength(src)
    end
    copyto!(amp_scales, host_amp_scales)
    copyto!(opd_to_cycles, host_opd_to_cycles)
    total = n_spots * n_src
    fft_view = @view wfs.front_end.propagation.fft_asterism_stack[:, :, 1:total]
    intensity_view = @view wfs.front_end.propagation.intensity_tmp_stack[:, :, 1:total]
    launch_kernel_async!(style, sh_field_asterism_stack_kernel!, fft_view, wfs.front_end.layout.valid_mask,
        pupil_reflectivity(tel), tel.state.opd, wfs.front_end.propagation.phasor, amp_scales, opd_to_cycles,
        n_sub, sub, ox, oy, tel.params.resolution, pad, n_spots, n_src;
        ndrange=(pad, pad, n_src, n_sub, n_sub))
    synchronize_backend!(style)
    execute_fft_plan!(fft_view, wfs.front_end.propagation.fft_asterism_plan)
    intensity_scale = sh_fft_intensity_scale(T, pad)
    phase = begin_kernel_phase(style)
    queue_kernel!(phase, complex_abs2_stack_kernel!, intensity_view, fft_view,
        intensity_scale, pad, total; ndrange=size(intensity_view))
    queue_kernel!(phase, reduce_grouped_blocks_kernel!, wfs.front_end.propagation.intensity_stack, intensity_view,
        n_spots, n_src, size(wfs.front_end.propagation.intensity_stack, 1), size(wfs.front_end.propagation.intensity_stack, 2); ndrange=size(wfs.front_end.propagation.intensity_stack))
    finish_kernel_phase!(phase)
    return wfs.front_end.propagation.intensity_stack
end

@inline sh_spectral_component_qe(::Nothing, sample,
    ::Type{T}) where {T<:AbstractFloat} = one(T)

@inline sh_spectral_component_qe(model::AbstractQuantumEfficiencyModel,
    sample, ::Type{T}) where {T<:AbstractFloat} =
    T(qe_at(model, sample.wavelength))

function compute_intensity_spectral_stack!(style::AcceleratorStyle,
    wfs::ShackHartmannWFS, tel::Telescope, src::SpectralSource,
    qe_model::Union{Nothing,AbstractQuantumEfficiencyModel})
    grid_wavelength = require_sh_common_spectral_grid(wfs, src)
    bundle = spectral_bundle(src)
    n_src = length(bundle)
    ensure_sh_asterism_buffers!(wfs.front_end, n_src)
    pad = size(wfs.front_end.propagation.fft_stack, 1)
    n_sub = n_lenslets(wfs)
    n_spots = n_sub * n_sub
    sub = div(tel.params.resolution, n_sub)
    ox = div(pad - sub, 2)
    oy = div(pad - sub, 2)
    T = eltype(wfs.front_end.propagation.intensity)
    amp_scales = wfs.front_end.propagation.amp_scales
    host_amp_scales = wfs.front_end.propagation.amp_scales_host
    opd_to_cycles = wfs.front_end.propagation.opd_to_cycles
    host_opd_to_cycles = wfs.front_end.propagation.opd_to_cycles_host
    base_photon_rate = T(photon_irradiance(src) * (tel.params.diameter / tel.params.resolution)^2)
    @inbounds for i in eachindex(bundle.samples)
        sample = bundle.samples[i]
        channel_qe = sh_spectral_component_qe(qe_model, sample, T)
        host_amp_scales[i] = sqrt(base_photon_rate * sample.weight *
            channel_qe)
        host_opd_to_cycles[i] = T(2) / grid_wavelength
    end
    copyto!(amp_scales, host_amp_scales)
    copyto!(opd_to_cycles, host_opd_to_cycles)
    total = n_spots * n_src
    fft_view = @view wfs.front_end.propagation.fft_asterism_stack[:, :, 1:total]
    intensity_view = @view wfs.front_end.propagation.intensity_tmp_stack[:, :, 1:total]
    launch_kernel_async!(style, sh_field_asterism_stack_kernel!, fft_view, wfs.front_end.layout.valid_mask,
        pupil_reflectivity(tel), tel.state.opd, wfs.front_end.propagation.phasor, amp_scales, opd_to_cycles,
        n_sub, sub, ox, oy, tel.params.resolution, pad, n_spots, n_src;
        ndrange=(pad, pad, n_src, n_sub, n_sub))
    synchronize_backend!(style)
    execute_fft_plan!(fft_view, wfs.front_end.propagation.fft_asterism_plan)
    intensity_scale = sh_fft_intensity_scale(T, pad)
    phase = begin_kernel_phase(style)
    queue_kernel!(phase, complex_abs2_stack_kernel!, intensity_view, fft_view,
        intensity_scale, pad, total; ndrange=size(intensity_view))
    queue_kernel!(phase, reduce_grouped_blocks_kernel!, wfs.front_end.propagation.intensity_stack, intensity_view,
        n_spots, n_src, size(wfs.front_end.propagation.intensity_stack, 1), size(wfs.front_end.propagation.intensity_stack, 2); ndrange=size(wfs.front_end.propagation.intensity_stack))
    finish_kernel_phase!(phase)
    return wfs.front_end.propagation.intensity_stack
end

compute_intensity_spectral_stack!(style::AcceleratorStyle,
    wfs::ShackHartmannWFS, tel::Telescope, src::SpectralSource) =
    compute_intensity_spectral_stack!(style, wfs, tel, src, nothing)

@inline function copy_stack_plane_to_matrix!(dest::AbstractMatrix{T}, stack::AbstractArray{T,3}, idx::Int) where {T}
    @inbounds for y in axes(dest, 2), x in axes(dest, 1)
        dest[x, y] = stack[x, y, idx]
    end
    return dest
end

@inline function copy_matrix_to_stack_plane!(stack::AbstractArray{T,3}, idx::Int, src::AbstractMatrix{T}) where {T}
    @inbounds for y in axes(src, 2), x in axes(src, 1)
        stack[idx, x, y] = src[x, y]
    end
    return stack
end

function sample_spot_stack!(::ScalarCPUStyle,
    front_end::ShackHartmannOpticalFrontEnd)
    propagation = front_end.propagation
    n_spots = size(propagation.sampled_spot_cube, 1)
    @inbounds for idx in 1:n_spots
        copy_stack_plane_to_matrix!(propagation.intensity,
            propagation.intensity_stack, idx)
        sample_spot!(front_end, propagation.intensity)
        copy_matrix_to_stack_plane!(propagation.sampled_spot_cube, idx,
            propagation.spot)
    end
    return propagation.sampled_spot_cube
end

function sample_spot_stack!(style::AcceleratorStyle,
    front_end::ShackHartmannOpticalFrontEnd)
    propagation = front_end.propagation
    pad = size(propagation.intensity_stack, 1)
    binning = propagation.binning_pixel_scale
    @assert pad % binning == 0
    n_binned = div(pad, binning)
    n_out = size(propagation.sampled_spot_cube, 2)
    ox = div(n_out - n_binned, 2)
    oy = div(n_out - n_binned, 2)
    n_sub = n_lenslets(front_end)
    launch_kernel!(style, sh_sample_spot_stack_kernel!,
        propagation.sampled_spot_cube, propagation.intensity_stack,
        front_end.layout.valid_mask, binning, n_sub, n_binned, n_out, ox,
        oy;
        ndrange=(n_sub, n_sub, n_out, n_out))
    return propagation.sampled_spot_cube
end

@inline require_sh_asterism_common_wavelength(ast::Asterism) = wavelength(ast)

function accumulate_sh_asterism_spots!(style::ExecutionStyle, wfs::ShackHartmannWFS,
    tel::Telescope, ast::Asterism)
    require_sh_asterism_common_wavelength(ast)
    fill!(wfs.front_end.propagation.spot_cube_accum, zero(eltype(wfs.front_end.propagation.spot_cube_accum)))
    @inbounds for src in ast.sources
        sampled_spots_peak!(style, wfs, tel, src)
        wfs.front_end.propagation.spot_cube_accum .+= wfs.acquisition.spot_cube
    end
    copyto!(wfs.acquisition.spot_cube, wfs.front_end.propagation.spot_cube_accum)
    return sh_safe_peak_value(wfs.acquisition.spot_cube)
end

function capture_sh_asterism_spots!(style::ExecutionStyle, wfs::ShackHartmannWFS,
    ast::Asterism, det::AbstractDetector, rng::AbstractRNG)
    capture_stack!(det, wfs.acquisition.spot_cube, wfs.acquisition.detector_noise_cube,
        first(ast.sources), rng)
    zero_invalid_sh_spot_cube!(style, wfs.acquisition.spot_cube,
        wfs.front_end.layout.valid_mask)
    return sh_safe_peak_value(wfs.acquisition.spot_cube)
end

function capture_sh_source_spots!(style::ExecutionStyle, wfs::ShackHartmannWFS,
    src::AbstractSource, det::AbstractDetector, rng::AbstractRNG)
    capture_stack!(det, wfs.acquisition.spot_cube, wfs.acquisition.detector_noise_cube,
        src, rng)
    zero_invalid_sh_spot_cube!(style, wfs.acquisition.spot_cube,
        wfs.front_end.layout.valid_mask)
    return sh_safe_peak_value(wfs.acquisition.spot_cube)
end

function capture_sh_qe_weighted_spots!(style::ExecutionStyle,
    wfs::ShackHartmannWFS, det::AbstractDetector, rng::AbstractRNG)
    capture_stack_with_quantum_efficiency!(det, wfs.acquisition.spot_cube,
        wfs.acquisition.detector_noise_cube, one(eltype(wfs.acquisition.spot_cube)), rng)
    zero_invalid_sh_spot_cube!(style, wfs.acquisition.spot_cube,
        wfs.front_end.layout.valid_mask)
    return sh_safe_peak_value(wfs.acquisition.spot_cube)
end

function finalize_sh_asterism_slopes!(wfs::ShackHartmannWFS, peak)
    sync_exported_spots!(wfs)
    sh_signal_from_spots!(wfs, peak, slope_extraction_model(wfs))
    subtract_reference_and_scale!(wfs)
    return wfs.estimator.slopes
end

function sampled_spots_peak_asterism_stacked!(::ScalarCPUStyle, wfs::ShackHartmannWFS, tel::Telescope, ast::Asterism)
    require_sh_asterism_common_wavelength(ast)
    fill!(wfs.acquisition.detector_noise_cube, zero(eltype(wfs.acquisition.detector_noise_cube)))
    for src in ast.sources
        sampled_spots_peak!(ScalarCPUStyle(), wfs, tel, src)
        wfs.acquisition.detector_noise_cube .+= wfs.acquisition.spot_cube
    end
    copyto!(wfs.acquisition.spot_cube, wfs.acquisition.detector_noise_cube)
    return sh_safe_peak_value(wfs.acquisition.spot_cube)
end

function sampled_spots_peak_asterism_stacked!(::ScalarCPUStyle, wfs::ShackHartmannWFS, tel::Telescope, ast::Asterism,
    det::AbstractDetector, rng::AbstractRNG)
    sampled_spots_peak_asterism_stacked!(ScalarCPUStyle(), wfs, tel, ast)
    return capture_sh_asterism_spots!(ScalarCPUStyle(), wfs, ast, det, rng)
end

function sampled_spots_peak_asterism_stacked!(style::AcceleratorStyle, wfs::ShackHartmannWFS, tel::Telescope, ast::Asterism)
    compute_intensity_asterism_stack!(style, wfs, tel, ast)
    sample_spot_stack!(style, wfs.front_end)
    sync_signal_spots_from_sampled!(wfs)
    return sh_safe_peak_value(wfs.acquisition.spot_cube)
end

function sampled_spots_peak_asterism_stacked!(style::AcceleratorStyle, wfs::ShackHartmannWFS, tel::Telescope, ast::Asterism,
    det::AbstractDetector, rng::AbstractRNG)
    sampled_spots_peak_asterism_stacked!(style, wfs, tel, ast)
    return capture_sh_asterism_spots!(style, wfs, ast, det, rng)
end

@kernel function sh_spot_centroid_stats_kernel!(stats, spot_cube, valid_mask, threshold, n_sub::Int, n1::Int, n2::Int)
    i, j = @index(Global, NTuple)
    if i <= n_sub && j <= n_sub
        idx = sh_lenslet_index(i, j, n_sub)
        base = 3 * (idx - 1)
        T = eltype(stats)
        peak = zero(T)
        total = zero(T)
        sx = zero(T)
        sy = zero(T)
        if @inbounds valid_mask[i, j]
            @inbounds for y in 1:n2, x in 1:n1
                val = spot_cube[idx, x, y]
                peak = max(peak, val)
            end
            cutoff = threshold * peak
            @inbounds for y in 1:n2, x in 1:n1
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
    centroid_response, n_src::Int, n_sub::Int, offset::Int)
    i, j = @index(Global, NTuple)
    if i <= n_sub && j <= n_sub
        idx = sh_lenslet_index(i, j, n_sub)
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
                sx_sum = stats_accum[base + 2]
                sy_sum = stats_accum[base + 3]
                inv_nsrc = inv(T(n_src))
                @inbounds begin
                    slopes[idx] = ((sx_sum * inv_nsrc) - reference[idx]) / centroid_response
                    slopes[idx + offset] = ((sy_sum * inv_nsrc) - reference[idx + offset]) / centroid_response
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

@inline function _sampled_spots_peak_source_batched!(style::AcceleratorStyle, wfs::ShackHartmannWFS, tel::Telescope, src::AbstractSource)
    return sampled_spots_peak!(style, wfs, tel, src)
end

@inline function _sampled_spots_peak_source_batched!(style::AcceleratorStyle, wfs::ShackHartmannWFS, tel::Telescope, src::LGSSource)
    return sampled_spots_peak!(style, wfs, tel, src)
end

function measure_sh_asterism_batched!(style::AcceleratorStyle, wfs::ShackHartmannWFS, tel::Telescope, ast::Asterism)
    peak = accumulate_sh_asterism_spots!(style, wfs, tel, ast)
    return finalize_sh_asterism_slopes!(wfs, peak)
end

function measure_sh_asterism_batched!(style::AcceleratorStyle, wfs::ShackHartmannWFS, tel::Telescope,
    ast::Asterism, det::AbstractDetector, rng::AbstractRNG)
    accumulate_sh_asterism_spots!(style, wfs, tel, ast)
    peak = capture_sh_asterism_spots!(style, wfs, ast, det, rng)
    return finalize_sh_asterism_slopes!(wfs, peak)
end

@kernel function sh_spot_centroid_kernel!(slopes, spot_cube, valid_mask, cutoff, n_sub::Int, offset::Int, n1::Int, n2::Int)
    i, j = @index(Global, NTuple)
    if i <= n_sub && j <= n_sub
        idx = sh_lenslet_index(i, j, n_sub)
        T = eltype(slopes)
        total = zero(T)
        sx = zero(T)
        sy = zero(T)
        if @inbounds valid_mask[i, j]
            @inbounds for y in 1:n2, x in 1:n1
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
                slopes[idx] = sx / total
                slopes[idx + offset] = sy / total
            end
        else
            @inbounds begin
                slopes[idx] = zero(T)
                slopes[idx + offset] = zero(T)
            end
        end
    end
end

@kernel function sh_spot_centroid_reference_scale_kernel!(slopes, spot_cube, reference, valid_mask,
    cutoff, centroid_response, n_sub::Int, offset::Int, n1::Int, n2::Int)
    i, j = @index(Global, NTuple)
    if i <= n_sub && j <= n_sub
        idx = sh_lenslet_index(i, j, n_sub)
        T = eltype(slopes)
        total = zero(T)
        sx = zero(T)
        sy = zero(T)
        inv_units = inv(centroid_response)
        ref_x = @inbounds reference[idx]
        ref_y = @inbounds reference[idx + offset]
        if @inbounds valid_mask[i, j]
            @inbounds for y in 1:n2, x in 1:n1
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
                slopes[idx] = -ref_x * inv_units
                slopes[idx + offset] = -ref_y * inv_units
            else
                slopes[idx] = (sx / total - ref_x) * inv_units
                slopes[idx + offset] = (sy / total - ref_y) * inv_units
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
    i, j = @index(Global, NTuple)
    if i <= n_sub && j <= n_sub
        idx = sh_lenslet_index(i, j, n_sub)
        base = 3 * (idx - 1)
        T = eltype(stats)
        total = zero(T)
        sx = zero(T)
        sy = zero(T)
        if @inbounds valid_mask[i, j]
            @inbounds for y in 1:n2, x in 1:n1
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
    i, j = @index(Global, NTuple)
    if i <= n_sub && j <= n_sub
        idx = sh_lenslet_index(i, j, n_sub)
        base = 3 * (idx - 1)
        T = eltype(slopes)
        if @inbounds valid_mask[i, j]
            total = @inbounds stats[base + 1]
            if total > zero(T)
                sx = @inbounds stats[base + 2] / total
                sy = @inbounds stats[base + 3] / total
                @inbounds begin
                    slopes[idx] = sx
                    slopes[idx + offset] = sy
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

@kernel function sh_finalize_spot_slopes_reference_scale_kernel!(slopes, stats, reference, valid_mask,
    centroid_response, n_sub::Int, offset::Int)
    i, j = @index(Global, NTuple)
    if i <= n_sub && j <= n_sub
        idx = sh_lenslet_index(i, j, n_sub)
        base = 3 * (idx - 1)
        T = eltype(slopes)
        inv_units = inv(centroid_response)
        ref_x = @inbounds reference[idx]
        ref_y = @inbounds reference[idx + offset]
        if @inbounds valid_mask[i, j]
            total = @inbounds stats[base + 1]
            if total > zero(T)
                sx = @inbounds stats[base + 2] / total
                sy = @inbounds stats[base + 3] / total
                @inbounds begin
                    slopes[idx] = (sx - ref_x) * inv_units
                    slopes[idx + offset] = (sy - ref_y) * inv_units
                end
            else
                @inbounds begin
                    slopes[idx] = -ref_x * inv_units
                    slopes[idx + offset] = -ref_y * inv_units
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
    i, j, x, y = @index(Global, NTuple)
    if i <= n_sub && j <= n_sub && x <= n1 && y <= n2
        idx = sh_lenslet_index(i, j, n_sub)
        if !(@inbounds valid_mask[i, j])
            @inbounds spot_cube[idx, x, y] = zero(eltype(spot_cube))
        end
    end
end

@kernel function zero_invalid_sh_slopes_kernel!(slopes, valid_mask, n_sub::Int, offset::Int)
    i, j = @index(Global, NTuple)
    if i <= n_sub && j <= n_sub
        idx = sh_lenslet_index(i, j, n_sub)
        if !(@inbounds valid_mask[i, j])
            T = eltype(slopes)
            @inbounds begin
                slopes[idx] = zero(T)
                slopes[idx + offset] = zero(T)
            end
        end
    end
end
