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
