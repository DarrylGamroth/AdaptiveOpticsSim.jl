"""Prepared-contract violation at one semantic wavefront-sensor stage."""
struct WFSPreparationError <: AdaptiveOpticsSimError
    stage::Symbol
    reason::Symbol
    msg::String
end

"""
Wavefront sensors expose the physical stages they own. Built-in families
provide optical formation and detector acquisition; operational estimators
belong to the composing package.

Optional detector coupling, runtime preparation, stacked-source support, and
grouped execution are expressed through capability queries rather than
subtype-specific conditionals.
"""
abstract type AbstractWFS <: AbstractOpticalElement end
