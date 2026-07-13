struct AOSimulation{TEL<:AbstractTelescope,S<:AbstractSource,SS<:AbstractSource,A<:AbstractAtmosphere,O<:AbstractControllableOptic,W<:AbstractWFS,B<:AbstractArrayBackend}
    tel::TEL
    src::S
    science_src::SS
    atm::A
    optic::O
    wfs::W
end

@inline backend(::AOSimulation{<:Any,<:Any,<:Any,<:Any,<:Any,<:Any,B}) where {B} = B()
@inline wfs_source(simulation::AOSimulation) = simulation.src
@inline science_source(simulation::AOSimulation) = simulation.science_src

function AOSimulation(tel::AbstractTelescope, src::AbstractSource,
    science_src::AbstractSource, atm::AbstractAtmosphere,
    optic::AbstractControllableOptic, wfs::AbstractWFS)
    selector = require_same_backend(tel, atm, optic, wfs)
    return AOSimulation{
        typeof(tel),
        typeof(src),
        typeof(science_src),
        typeof(atm),
        typeof(optic),
        typeof(wfs),
        typeof(selector),
    }(tel, src, science_src, atm, optic, wfs)
end

function AOSimulation(tel::AbstractTelescope, src::AbstractSource,
    atm::AbstractAtmosphere, optic::AbstractControllableOptic, wfs::AbstractWFS;
    science_source::AbstractSource=src)
    return AOSimulation(tel, src, science_source, atm, optic, wfs)
end

function AOSimulation(; telescope::AbstractTelescope, source::AbstractSource,
    atmosphere::AbstractAtmosphere, optic::AbstractControllableOptic,
    sensor::AbstractWFS, science_source::AbstractSource=source)
    return AOSimulation(telescope, source, science_source, atmosphere, optic, sensor)
end
