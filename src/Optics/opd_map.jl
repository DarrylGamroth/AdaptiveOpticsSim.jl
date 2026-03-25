"""
    OPDMap(opd)

Wrap a fixed OPD array as an optical element.

This is the simplest way to inject a precomputed phase screen or static
aberration map into the telescope pipeline.
"""
struct OPDMap{T<:AbstractFloat,A<:AbstractMatrix{T}} <: AbstractOpticalElement
    opd::A
end

"""
    apply!(map, tel, mode)

Apply a fixed OPD map to the telescope either additively or by replacement,
following the same semantics as DM application modes.
"""
function apply!(map::OPDMap, tel::Telescope, ::DMAdditive)
    if size(map.opd) != size(tel.state.opd)
        throw(DimensionMismatchError("OPD map size does not match telescope resolution"))
    end
    tel.state.opd .+= map.opd
    return tel
end

function apply!(map::OPDMap, tel::Telescope, ::DMReplace)
    if size(map.opd) != size(tel.state.opd)
        throw(DimensionMismatchError("OPD map size does not match telescope resolution"))
    end
    tel.state.opd .= map.opd
    return tel
end
