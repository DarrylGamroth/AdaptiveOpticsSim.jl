"""Immutable physical definition of a pyramid focal-plane phase mask."""
struct PyramidPhaseMask{T<:AbstractFloat}
    old_mask::Bool
    rooftop::T
    theta_rotation::T
    mask_scale::T
    diffraction_padding::Int
    psf_centering::Bool
    n_pix_separation::Union{Int,Nothing}
    n_pix_edge::Union{Int,Nothing}
end

"""Immutable physical definition of the four BioEdge amplitude masks."""
struct BioEdgeAmplitudeMask{T<:AbstractFloat}
    grey_width::T
    grey_length::Union{Bool,T}
    diffraction_padding::Int
    psf_centering::Bool
    n_pix_separation::Union{Int,Nothing}
    n_pix_edge::Union{Int,Nothing}
end

"""Immutable focal-plane phase-shifting spot for a Zernike WFS."""
struct ZernikePhaseSpot{T<:AbstractFloat}
    phase_shift_pi::T
    radius_lambda_over_d::T
    diffraction_padding::Int
end

"""Immutable equal-and-opposite defocus pair used by a Curvature WFS."""
struct CurvatureDefocusPair{T<:AbstractFloat}
    rms_nm::T
    diffraction_padding::Int
end
