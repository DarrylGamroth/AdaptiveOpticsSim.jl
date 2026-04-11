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

function initialize_ao_pyramid(; resolution::Int, diameter::Real, sampling_time::Real,
    r0::Real, L0::Real=25.0, fractional_cn2::AbstractVector=[1.0],
    wind_speed::AbstractVector=[0.0], wind_direction::AbstractVector=[0.0], altitude::AbstractVector=[0.0],
    central_obstruction::Real=0.0, band::Symbol=:I, magnitude::Real=0.0,
    n_act::Int, influence_width::Real=0.2, n_subap::Int, threshold::Real=0.1,
    modulation::Real=2.0, T::Type{<:AbstractFloat}=Float64, backend::AbstractArrayBackend=CPUBackend())

    selector = backend

    tel = Telescope(resolution=resolution, diameter=diameter, sampling_time=sampling_time,
        central_obstruction=central_obstruction, T=T, backend=selector)
    src = Source(band=band, magnitude=magnitude, T=T)
    atm = MultiLayerAtmosphere(tel; r0=r0, L0=L0, fractional_cn2=fractional_cn2,
        wind_speed=wind_speed, wind_direction=wind_direction, altitude=altitude, T=T, backend=selector)
    dm = DeformableMirror(tel; n_act=n_act, influence_width=influence_width, T=T, backend=selector)
    wfs = PyramidWFS(tel; n_subap=n_subap, threshold=threshold, modulation=modulation, T=T, backend=selector)
    return AOSimulation(tel, src, atm, dm, wfs)
end

function initialize_ao_shwfs(; resolution::Int, diameter::Real, sampling_time::Real,
    r0::Real, L0::Real=25.0, fractional_cn2::AbstractVector=[1.0],
    wind_speed::AbstractVector=[0.0], wind_direction::AbstractVector=[0.0], altitude::AbstractVector=[0.0],
    central_obstruction::Real=0.0, band::Symbol=:I, magnitude::Real=0.0,
    n_act::Int, influence_width::Real=0.2, n_subap::Int, threshold::Real=0.1,
    T::Type{<:AbstractFloat}=Float64, backend::AbstractArrayBackend=CPUBackend())

    selector = backend

    tel = Telescope(resolution=resolution, diameter=diameter, sampling_time=sampling_time,
        central_obstruction=central_obstruction, T=T, backend=backend)
    src = Source(band=band, magnitude=magnitude, T=T)
    atm = MultiLayerAtmosphere(tel; r0=r0, L0=L0, fractional_cn2=fractional_cn2,
        wind_speed=wind_speed, wind_direction=wind_direction, altitude=altitude, T=T, backend=selector)
    dm = DeformableMirror(tel; n_act=n_act, influence_width=influence_width, T=T, backend=selector)
    wfs = ShackHartmann(tel; n_subap=n_subap, threshold=threshold, T=T, backend=selector)
    return AOSimulation(tel, src, atm, dm, wfs)
end
