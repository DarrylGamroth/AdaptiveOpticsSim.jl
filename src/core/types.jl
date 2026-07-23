"""Base type for optical elements."""
abstract type AbstractOpticalElement end

"""Telescopes own immutable pupil geometry used to prepare optical paths."""
abstract type AbstractTelescope <: AbstractOpticalElement end

"""Sources must provide wavelength(src)."""
abstract type AbstractSource <: AbstractOpticalElement end

"""Base type for atmospheric models."""
abstract type AbstractAtmosphere <: AbstractOpticalElement end

"""
Atmospheres whose shared evolution is advanced with explicit model time and
published as immutable current-state epoch identity tokens.
"""
abstract type AbstractTimedAtmosphere <: AbstractAtmosphere end

# Extension seam retained for untimed/static test atmospheres and user models.
# Maintained stochastic atmospheres use `advance_by!` / `advance_to!` instead.
function advance! end

"""
Wavefront sensors implement `measure!(wfs, pupil[, src])` and
`update_valid_mask!(wfs, pupil)`.

Optional behavior such as detector-coupled measurement, runtime preparation,
stacked-source support, and grouped execution is expressed through the
capability queries documented in the runtime/control interface layer.
"""
abstract type AbstractWFS <: AbstractOpticalElement end

function apply_shift_wfs!(::AbstractWFS; sx, sy)
    throw(InvalidConfiguration("apply_shift_wfs! is not supported for this WFS"))
end

"""Detectors implement capture!(det, psf; rng)."""
abstract type AbstractDetector <: AbstractOpticalElement end

"""Controllable optics form surfaces that are applied to explicit pupil paths."""
abstract type AbstractControllableOptic <: AbstractOpticalElement end

abstract type DMApplyMode end
struct DMAdditive <: DMApplyMode end
struct DMReplace <: DMApplyMode end

"""Deformable mirrors implement prepared influence functions and surface formation."""
abstract type AbstractDeformableMirror <: AbstractControllableOptic end
