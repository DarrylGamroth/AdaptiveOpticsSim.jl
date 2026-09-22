"""Prepared-contract violation at one semantic wavefront-sensor stage."""
struct WFSPreparationError <: AdaptiveOpticsSimError
    stage::Symbol
    reason::Symbol
    msg::String
end

"""
Wavefront sensors expose the physical stages they own and may extend the
prepared WFS optics, acquisition, and estimation protocols. Sensor families
that also own an in-package estimator implement
`measure!(wfs, pupil[, source])`; a physical-only front end need not do so.

Optional detector coupling, runtime preparation, stacked-source support, and
grouped execution are expressed through capability queries rather than
subtype-specific conditionals.
"""
abstract type AbstractWFS <: AbstractOpticalElement end

function apply_shift_wfs!(::AbstractWFS; sx, sy)
    throw(InvalidConfiguration(
        "apply_shift_wfs! is not supported for this WFS"))
end
