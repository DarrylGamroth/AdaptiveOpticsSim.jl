"""
    NCPA(opd)

An explicit sampled non-common-path aberration represented as optical path
difference in metres. The array remains on its existing backend and is applied
to a selected optical branch through `apply_surface!`.

Basis selection, coefficient generation, and atmosphere- or DM-derived
synthesis are calibration policy and are intentionally not stored by this
physical optical element.
"""
struct NCPA{T<:AbstractFloat,A<:AbstractMatrix{T}} <: AbstractOpticalElement
    opd::A
end

@inline surface_opd(ncpa::NCPA) = ncpa.opd
