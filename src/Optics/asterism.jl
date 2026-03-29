struct Asterism{S<:AbstractSource,V<:AbstractVector{S}} <: AbstractSource
    sources::V
end

Base.length(ast::Asterism) = length(ast.sources)

function wavelength(ast::Asterism)
    if isempty(ast.sources)
        throw(InvalidConfiguration("asterism must contain at least one source"))
    end
    w0 = wavelength(ast.sources[1])
    for src in ast.sources[2:end]
        if wavelength(src) != w0
            throw(InvalidConfiguration("asterism sources must share a common wavelength"))
        end
    end
    return w0
end

function psf_pixel_scale_arcsec(tel::Telescope, src::Source, zero_padding::Int)
    return (180 * 3600 / π) * wavelength(src) / tel.params.diameter / zero_padding
end

function shift_psf!(out::AbstractMatrix, psf::AbstractMatrix, dx::Real, dy::Real)
    circshift2d!(out, psf, (round(Int, dy), round(Int, dx)))
    return out
end

function compute_psf!(tel::Telescope, ast::Asterism; zero_padding::Int=1, ws::Union{Workspace,Nothing}=nothing)
    if isempty(ast.sources)
        throw(InvalidConfiguration("asterism must contain at least one source"))
    end
    n_pad = tel.params.resolution * zero_padding
    T = eltype(tel.state.opd)
    combined = similar(tel.state.psf, T, n_pad, n_pad)
    fill!(combined, zero(T))
    tel.state.psf_stack = similar(tel.state.psf, T, n_pad, n_pad, length(ast.sources))
    shifted = similar(tel.state.psf, T, n_pad, n_pad)

    for (i, src) in enumerate(ast.sources)
        psf = compute_psf_centered!(tel, src, something(ws, ensure_psf_workspace!(tel, n_pad)), zero_padding)
        scale = psf_pixel_scale_arcsec(tel, src, zero_padding)
        dx_arcsec, dy_arcsec = coordinates_xy_arcsec(src)
        shift_psf!(shifted, psf, dx_arcsec / scale, dy_arcsec / scale)
        @views tel.state.psf_stack[:, :, i] .= shifted
        combined .+= shifted
    end
    ensure_psf_state!(tel, n_pad)
    tel.state.psf .= combined
    return tel.state.psf
end
