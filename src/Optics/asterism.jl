struct Asterism{S<:AbstractSource,V<:AbstractVector{S}} <: AbstractSource
    sources::V
end

Base.length(ast::Asterism) = length(ast.sources)

function coordinates_xy_arcsec(src::Source)
    r = src.params.coordinates[1]
    theta = src.params.coordinates[2]
    x = r * cosd(theta)
    y = r * sind(theta)
    return x, y
end

function psf_pixel_scale_arcsec(tel::Telescope, src::Source, zero_padding::Int)
    return 206265 * wavelength(src) / tel.params.diameter / zero_padding
end

function shift_psf(psf::AbstractMatrix, dx::Real, dy::Real)
    return circshift(psf, (round(Int, dy), round(Int, dx)))
end

function compute_psf!(tel::Telescope, ast::Asterism; zero_padding::Int=1, ws::Union{Workspace,Nothing}=nothing)
    if isempty(ast.sources)
        throw(InvalidConfiguration("asterism must contain at least one source"))
    end
    n_pad = tel.params.resolution * zero_padding
    T = eltype(tel.state.opd)
    combined = zeros(T, n_pad, n_pad)
    tel.state.psf_list = Vector{typeof(tel.state.psf)}(undef, length(ast.sources))

    for (i, src) in enumerate(ast.sources)
        psf = compute_psf!(tel, src; zero_padding=zero_padding, ws=ws)
        scale = psf_pixel_scale_arcsec(tel, src, zero_padding)
        dx_arcsec, dy_arcsec = coordinates_xy_arcsec(src)
        shifted = shift_psf(psf, dx_arcsec / scale, dy_arcsec / scale)
        tel.state.psf_list[i] = shifted
        combined .+= shifted
    end
    ensure_psf_state!(tel, n_pad)
    tel.state.psf .= combined
    return tel.state.psf
end
