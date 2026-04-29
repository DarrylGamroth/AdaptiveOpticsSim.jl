struct AOSimulation{TEL<:AbstractTelescope,S<:AbstractSource,A<:AbstractAtmosphere,O<:AbstractControllableOptic,W<:AbstractWFS,B<:AbstractArrayBackend}
    tel::TEL
    src::S
    atm::A
    optic::O
    wfs::W
end

@inline backend(::AOSimulation{<:Any,<:Any,<:Any,<:Any,<:Any,B}) where {B} = B()

function AOSimulation(tel::AbstractTelescope, src::AbstractSource, atm::AbstractAtmosphere, optic::AbstractControllableOptic, wfs::AbstractWFS)
    selector = require_same_backend(tel, atm, optic, wfs)
    return AOSimulation{typeof(tel), typeof(src), typeof(atm), typeof(optic), typeof(wfs), typeof(selector)}(tel, src, atm, optic, wfs)
end

function AOSimulation(; telescope::AbstractTelescope, source::AbstractSource,
    atmosphere::AbstractAtmosphere, optic::AbstractControllableOptic,
    sensor::AbstractWFS)
    return AOSimulation(telescope, source, atmosphere, optic, sensor)
end
