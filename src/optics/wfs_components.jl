"""
Immutable physical definition of a pyramid focal-plane phase mask.
`rotation_rad` is the rotation angle in radians applied by the mask-coordinate
transform. `pupil_shift_x_pixels` and `pupil_shift_y_pixels` are the
per-face displacement terms, ordered q1 through q4. With `old_mask=false`,
they use native, unbinned detector pixels; final binning divides the displacement
by `binning`. `old_mask=true` retains the historical quadrant phase normalization
and does not calibrate these terms to the same pixel units.
"""
struct PyramidPhaseMask{T<:AbstractFloat}
    old_mask::Bool
    rooftop::T
    rotation_rad::T
    mask_scale::T
    diffraction_padding::Int
    psf_centering::Bool
    n_pix_separation::Union{Int,Nothing}
    n_pix_edge::Union{Int,Nothing}
    pupil_shift_x_pixels::NTuple{4,T}
    pupil_shift_y_pixels::NTuple{4,T}
end

"""Immutable physical definition of the four Bi-O-edge amplitude masks."""
struct BiOEdgeAmplitudeMask{T<:AbstractFloat}
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
