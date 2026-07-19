"""
    OPDMap(opd)

Wrap a fixed OPD array as an optical element.

This is the simplest way to inject a precomputed phase screen or static
aberration map into the telescope pipeline.
"""
struct OPDMap{T<:AbstractFloat,A<:AbstractMatrix{T}} <: AbstractOpticalElement
    opd::A
end


@inline surface_opd(map::OPDMap) = map.opd
