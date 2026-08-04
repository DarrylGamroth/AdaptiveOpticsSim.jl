"""Prepared-contract violation at one semantic wavefront-sensor stage."""
struct WFSPreparationError <: AdaptiveOpticsSimError
    stage::Symbol
    reason::Symbol
    msg::String
end

"""
Wavefront sensors implement `measure!(wfs, pupil[, source])` and may extend the
prepared WFS optics, acquisition, and estimation protocols.

Optional detector coupling, runtime preparation, stacked-source support, and
grouped execution are expressed through capability queries rather than
subtype-specific conditionals.
"""
abstract type AbstractWFS <: AbstractOpticalElement end

function apply_shift_wfs!(::AbstractWFS; sx, sy)
    throw(InvalidConfiguration(
        "apply_shift_wfs! is not supported for this WFS"))
end
