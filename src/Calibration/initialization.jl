struct AOSimulation{A<:AbstractAtmosphere,S<:AbstractSource,O<:AbstractControllableOptic,W<:AbstractWFS}
    tel::Telescope
    atm::A
    src::S
    optic::O
    wfs::W
end

function initialize_ao_pyramid(; resolution::Int, diameter::Real, sampling_time::Real,
    r0::Real, L0::Real=25.0, fractional_cn2::AbstractVector=[1.0],
    wind_speed::AbstractVector=[0.0], wind_direction::AbstractVector=[0.0], altitude::AbstractVector=[0.0],
    central_obstruction::Real=0.0, band::Symbol=:I, magnitude::Real=0.0,
    n_act::Int, influence_width::Real=0.2, n_subap::Int, threshold::Real=0.1,
    modulation::Real=2.0, T::Type{<:AbstractFloat}=Float64, backend=CPUBackend())

    backend = _resolve_array_backend(backend)

    tel = Telescope(resolution=resolution, diameter=diameter, sampling_time=sampling_time,
        central_obstruction=central_obstruction, T=T, backend=backend)
    src = Source(band=band, magnitude=magnitude, T=T)
    atm = MultiLayerAtmosphere(tel; r0=r0, L0=L0, fractional_cn2=fractional_cn2,
        wind_speed=wind_speed, wind_direction=wind_direction, altitude=altitude, T=T, backend=backend)
    dm = DeformableMirror(tel; n_act=n_act, influence_width=influence_width, T=T, backend=backend)
    wfs = PyramidWFS(tel; n_subap=n_subap, threshold=threshold, modulation=modulation, T=T, backend=backend)
    return AOSimulation(tel, atm, src, dm, wfs)
end

function initialize_ao_shwfs(; resolution::Int, diameter::Real, sampling_time::Real,
    r0::Real, L0::Real=25.0, fractional_cn2::AbstractVector=[1.0],
    wind_speed::AbstractVector=[0.0], wind_direction::AbstractVector=[0.0], altitude::AbstractVector=[0.0],
    central_obstruction::Real=0.0, band::Symbol=:I, magnitude::Real=0.0,
    n_act::Int, influence_width::Real=0.2, n_subap::Int, threshold::Real=0.1,
    T::Type{<:AbstractFloat}=Float64, backend=CPUBackend())

    backend = _resolve_array_backend(backend)

    tel = Telescope(resolution=resolution, diameter=diameter, sampling_time=sampling_time,
        central_obstruction=central_obstruction, T=T, backend=backend)
    src = Source(band=band, magnitude=magnitude, T=T)
    atm = MultiLayerAtmosphere(tel; r0=r0, L0=L0, fractional_cn2=fractional_cn2,
        wind_speed=wind_speed, wind_direction=wind_direction, altitude=altitude, T=T, backend=backend)
    dm = DeformableMirror(tel; n_act=n_act, influence_width=influence_width, T=T, backend=backend)
    wfs = ShackHartmann(tel; n_subap=n_subap, threshold=threshold, T=T, backend=backend)
    return AOSimulation(tel, atm, src, dm, wfs)
end
