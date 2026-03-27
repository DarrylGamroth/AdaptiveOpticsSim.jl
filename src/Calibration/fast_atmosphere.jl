function fast_atmosphere(tel::Telescope, atm::KolmogorovAtmosphere; speed_factor::Real)
    if speed_factor <= 0
        throw(InvalidConfiguration("speed_factor must be > 0"))
    end
    params = tel.params
    tel_fast = Telescope(
        resolution=params.resolution,
        diameter=params.diameter,
        sampling_time=params.sampling_time / speed_factor,
        central_obstruction=params.central_obstruction,
        fov_arcsec=params.fov_arcsec,
        T=eltype(tel.state.opd),
    )
    atm_fast = KolmogorovAtmosphere(tel_fast; r0=atm.params.r0, L0=atm.params.L0, T=eltype(tel.state.opd))
    return tel_fast, atm_fast
end

function fast_atmosphere(tel::Telescope, atm::MultiLayerAtmosphere; speed_factor::Real)
    if speed_factor <= 0
        throw(InvalidConfiguration("speed_factor must be > 0"))
    end
    params = tel.params
    tel_fast = Telescope(
        resolution=params.resolution,
        diameter=params.diameter,
        sampling_time=params.sampling_time / speed_factor,
        central_obstruction=params.central_obstruction,
        fov_arcsec=params.fov_arcsec,
        T=eltype(tel.state.opd),
    )
    atm_fast = MultiLayerAtmosphere(tel_fast;
        r0=atm.params.r0,
        L0=atm.params.L0,
        fractional_cn2=atm.params.cn2_fractions,
        wind_speed=atm.params.wind_speed,
        wind_direction=atm.params.wind_direction,
        altitude=atm.params.altitude,
        T=eltype(tel.state.opd),
    )
    return tel_fast, atm_fast
end
