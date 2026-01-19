struct OPDMap{T<:AbstractFloat,A<:AbstractMatrix{T}} <: AbstractOpticalElement
    opd::A
end

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
