@kernel function apply_centering_phase_kernel!(field, phase_shift, n::Int)
    i, j = @index(Global, NTuple)
    if i <= n && j <= n
        @inbounds field[i, j] *= cis(phase_shift * (i + j - 2))
    end
end

function apply_centering_phase!(::ScalarCPUStyle, field::AbstractMatrix{Complex{T}}, phase_shift::T) where {T<:AbstractFloat}
    n, m = size(field)
    @inbounds for j in 1:m, i in 1:n
        field[i, j] *= cis(phase_shift * (i + j - 2))
    end
    return field
end

function apply_centering_phase!(style::AcceleratorStyle, field::AbstractMatrix{Complex{T}}, phase_shift::T) where {T<:AbstractFloat}
    launch_kernel!(style, apply_centering_phase_kernel!, field, phase_shift, size(field, 1); ndrange=size(field))
    return field
end

function ensure_psf_state!(tel::Telescope, n::Int)
    if size(tel.state.psf) != (n, n)
        tel.state.psf = similar(tel.state.psf, n, n)
    end
    return tel
end

function ensure_psf_workspace!(tel::Telescope, n::Int)
    ensure_psf_buffers!(tel.state.psf_workspace, n)
    return tel.state.psf_workspace
end

function compute_psf_centered!(tel::Telescope, src::Source, ws::Workspace, zero_padding::Int=1)
    n = tel.params.resolution
    if zero_padding < 1
        throw(InvalidConfiguration("zero_padding must be >= 1"))
    end
    n_pad = n * zero_padding
    T = eltype(tel.state.opd)

    ensure_psf_buffers!(ws, n_pad)

    fill!(ws.pupil_field, zero(eltype(ws.pupil_field)))
    opd_to_cycles = T(2) / wavelength(src)
    amp_scale = sqrt(T(photon_flux(src) * tel.params.sampling_time * (tel.params.diameter / tel.params.resolution)^2))
    ox = div(n_pad - n, 2)
    oy = div(n_pad - n, 2)
    @views @. ws.pupil_field[ox+1:ox+n, oy+1:oy+n] = amp_scale * sqrt(tel.state.pupil_reflectivity) * cispi(opd_to_cycles * tel.state.opd)
    if iseven(n_pad)
        phase_shift = -T(pi) * (T(n_pad) + one(T)) / T(n_pad)
        apply_centering_phase!(execution_style(ws.pupil_field), ws.pupil_field, phase_shift)
    end

    copyto!(ws.fft_buffer, ws.pupil_field)
    execute_fft_plan!(ws.fft_buffer, ws.fft_plan)

    fft_scale = inv(T(n_pad))
    @. ws.psf_buffer = abs2(ws.fft_buffer) * (fft_scale * fft_scale)

    ensure_psf_state!(tel, n_pad)
    copyto!(tel.state.psf, ws.psf_buffer)
    return tel.state.psf
end

function compute_psf!(tel::Telescope, src::Source, ws::Workspace, zero_padding::Int=1)
    psf = compute_psf_centered!(tel, src, ws, zero_padding)
    if iszero(src.params.coordinates[1])
        return psf
    end
    scale = psf_pixel_scale_arcsec(tel, src, zero_padding)
    dx_arcsec, dy_arcsec = coordinates_xy_arcsec(src)
    shift_psf!(tel.state.psf, ws.psf_buffer, dx_arcsec / scale, dy_arcsec / scale)
    return tel.state.psf
end

function compute_psf!(tel::Telescope, src::Source, zero_padding::Int)
    if zero_padding < 1
        throw(InvalidConfiguration("zero_padding must be >= 1"))
    end
    n_pad = tel.params.resolution * zero_padding
    ws = ensure_psf_workspace!(tel, n_pad)
    return compute_psf!(tel, src, ws, zero_padding)
end

function compute_psf!(tel::Telescope, src::Source; zero_padding::Int=1, ws::Union{Workspace,Nothing}=nothing)
    if ws === nothing
        return compute_psf!(tel, src, zero_padding)
    end
    return compute_psf!(tel, src, ws, zero_padding)
end
