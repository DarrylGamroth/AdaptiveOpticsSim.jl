"""Base type for optical elements."""
abstract type AbstractOpticalElement end

"""Sources must provide wavelength(src)."""
abstract type AbstractSource <: AbstractOpticalElement end

"""Atmospheres implement advance!(atm, tel; rng) and propagate!(atm, tel)."""
abstract type AbstractAtmosphere <: AbstractOpticalElement end

"""Wavefront sensors implement measure!(wfs, tel[, src]) and update_valid_mask!(wfs, tel)."""
abstract type AbstractWFS <: AbstractOpticalElement end

function apply_shift_wfs!(::AbstractWFS; sx, sy)
    throw(InvalidConfiguration("apply_shift_wfs! is not supported for this WFS"))
end

"""Detectors implement capture!(det, psf; rng)."""
abstract type AbstractDetector <: AbstractOpticalElement end

"""Deformable mirrors implement build_influence_functions!(dm, tel) and apply!(dm, tel, mode)."""
abstract type AbstractDeformableMirror <: AbstractOpticalElement end
