using FFTW

function ensure_psf_state!(tel::Telescope, n::Int)
    if size(tel.state.psf) != (n, n)
        tel.state.psf = similar(tel.state.psf, n, n)
    end
    return tel
end

function compute_psf!(tel::Telescope, src::Source, ws::Workspace, zero_padding::Int=1)
    n = tel.params.resolution
    if zero_padding < 1
        throw(InvalidConfiguration("zero_padding must be >= 1"))
    end
    n_pad = n * zero_padding
    T = eltype(tel.state.opd)

    ensure_psf_buffers!(ws, n_pad)

    fill!(ws.pupil_field, zero(eltype(ws.pupil_field)))
    phase_scale = T(2 * pi) / wavelength(src)
    ox = div(n_pad - n, 2)
    oy = div(n_pad - n, 2)
    @views @. ws.pupil_field[ox+1:ox+n, oy+1:oy+n] = tel.state.pupil * cis(phase_scale * tel.state.opd)

    copyto!(ws.fft_buffer, ws.pupil_field)
    mul!(ws.fft_buffer, ws.fft_plan, ws.fft_buffer)

    @. ws.psf_buffer = abs2(ws.fft_buffer)

    ensure_psf_state!(tel, n_pad)
    FFTW.fftshift!(tel.state.psf, ws.psf_buffer)
    return tel.state.psf
end

function compute_psf!(tel::Telescope, src::Source; zero_padding::Int=1, ws::Union{Workspace,Nothing}=nothing)
    n = tel.params.resolution
    if ws === nothing
        n_pad = n * zero_padding
        T = eltype(tel.state.opd)
        ws = Workspace(tel.state.opd, n_pad; T=T)
    end
    return compute_psf!(tel, src, ws, zero_padding)
end
