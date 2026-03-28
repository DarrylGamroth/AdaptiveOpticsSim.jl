using AdaptiveOpticsSim
using Random
using TOML

const REVOLT_MODELS = (:pwfs, :pwfs_unmod, :shwfs)
const REVOLT_PARAMETER_DIR = joinpath(@__DIR__, "parameter_files")

# Parameters

function model_file(kind::Symbol)
    kind === :pwfs && return "pwfs.toml"
    kind === :pwfs_unmod && return "pwfs_unmod.toml"
    kind === :shwfs && return "shwfs.toml"
    throw(ArgumentError("unsupported REVOLT model '$kind'; use :pwfs, :pwfs_unmod, or :shwfs"))
end

read_toml(filename::AbstractString) = TOML.parsefile(joinpath(REVOLT_PARAMETER_DIR, filename))

function load_revolt_parameters(kind::Symbol)
    return (
        common=read_toml("common.toml"),
        cameras=read_toml("cameras.toml"),
        model=read_toml(model_file(kind)),
    )
end

# Sources

@inline function scaled_magnitude(base_magnitude::Real, flux_scale::Real)
    flux_scale > 0 || throw(ArgumentError("flux_scale must be > 0"))
    return base_magnitude - 2.5 * log10(flux_scale)
end

function source_from_flux_scale(; band::Symbol, wavelength=nothing, magnitude::Real,
    flux_scale::Real, T::Type{<:AbstractFloat}=Float64)
    effective_magnitude = scaled_magnitude(magnitude, flux_scale)
    if wavelength === nothing
        return Source(band=band, magnitude=effective_magnitude, T=T)
    end
    return Source(band=band, wavelength=wavelength, magnitude=effective_magnitude, T=T)
end

# Cameras and detectors

function sensor_factory(name::AbstractString)
    lowered = lowercase(name)
    lowered == "cmos" && return T -> CMOSSensor(T=T)
    lowered == "ingaas" && return T -> InGaAsSensor(T=T)
    lowered == "emccd" && return T -> EMCCDSensor(T=T)
    lowered == "ccd" && return T -> CCDSensor(T=T)
    throw(ArgumentError("unsupported sensor '$name'"))
end

function camera_configs(parameters, ::Type{T}) where {T<:AbstractFloat}
    configs = Dict{Symbol,NamedTuple}()
    for (name, camera) in parameters.cameras["cameras"]
        configs[Symbol(name)] = (
            resolution=Int(camera["resolution"]),
            bits=Int(camera["bits"]),
            full_well=T(camera["full_well"]),
            sensor=sensor_factory(camera["sensor"]),
            qe=T(camera["qe"]),
            dark_current=T(camera["dark_current"]),
            readout_noise=T(camera["readout_noise"]),
            photon_noise=Bool(camera["photon_noise"]),
            psf_sampling=T(camera["psf_sampling"]),
        )
    end
    return configs
end

function noise_model(camera)
    if camera.photon_noise
        return NoisePhotonReadout(camera.readout_noise)
    elseif iszero(camera.readout_noise)
        return NoiseNone()
    end
    return NoiseReadout(camera.readout_noise)
end

function detector_from_camera(camera; integration_time::Real, T::Type{<:AbstractFloat}=Float64,
    backend=Array, response_mode::Symbol=:default)
    response_model = response_mode === :default ? nothing :
        response_mode === :null ? NullFrameResponse() :
        throw(ArgumentError("unsupported response_mode '$response_mode'; use :default or :null"))
    return Detector(
        noise_model(camera);
        integration_time=integration_time,
        qe=camera.qe,
        gain=1.0,
        dark_current=camera.dark_current,
        bits=camera.bits,
        full_well=camera.full_well,
        sensor=camera.sensor(T),
        response_model=response_model,
        T=T,
        backend=backend,
    )
end

# Optical train

function misregistration(telescope_resolution::Integer, dm_cfg; T::Type{<:AbstractFloat}=Float64)
    pixel_scale = T(2 / telescope_resolution)
    return Misregistration(
        shift_x=T(dm_cfg["shift_x_px"]) * pixel_scale,
        shift_y=T(dm_cfg["shift_y_px"]) * pixel_scale,
        rotation_deg=T(dm_cfg["rotation_deg"]),
        anamorphosis_angle=T(dm_cfg["anamorphosis_angle"]),
        tangential_scaling=T(dm_cfg["tangential_scaling"]),
        radial_scaling=T(dm_cfg["radial_scaling"]),
        T=T,
    )
end

function wavefront_sensor(wfs_cfg::Dict, telescope; T::Type{<:AbstractFloat}=Float64, backend=Array)
    family = Symbol(wfs_cfg["family"])
    if family === :pyramid
        return PyramidWFS(
            telescope;
            n_subap=Int(wfs_cfg["n_subap"]),
            threshold=T(wfs_cfg["threshold"]),
            modulation=T(wfs_cfg["modulation"]),
            modulation_points=Int(wfs_cfg["modulation_points"]),
            light_ratio=T(wfs_cfg["light_ratio"]),
            n_pix_separation=Int(wfs_cfg["n_pix_separation"]),
            n_pix_edge=Int(wfs_cfg["n_pix_edge"]),
            psf_centering=Bool(wfs_cfg["psf_centering"]),
            mode=Diffractive(),
            T=T,
            backend=backend,
        )
    end
    if family === :shack_hartmann
        n_pix_subap = haskey(wfs_cfg, "n_pix_subap") ? Int(wfs_cfg["n_pix_subap"]) : nothing
        return ShackHartmann(
            telescope;
            n_subap=Int(wfs_cfg["n_subap"]),
            threshold=T(wfs_cfg["threshold"]),
            shannon_sampling=Bool(wfs_cfg["shannon_sampling"]),
            n_pix_subap=n_pix_subap,
            mode=Diffractive(),
            T=T,
            backend=backend,
        )
    end
    throw(ArgumentError("unsupported WFS family '$family'"))
end

# Runtime

function revolt_setup(kind::Symbol; T::Type{<:AbstractFloat}=Float64, backend=Array,
    response_mode::Symbol=:default)
    parameters = load_revolt_parameters(kind)
    common = parameters.common
    model = parameters.model
    cameras = camera_configs(parameters, T)

    model_info = model["model"]
    telescope_info = common["telescope"]
    telescope_override = get(model, "telescope", Dict{String,Any}())
    source_info = model["sources"]
    wfs_info = model["wfs"]
    detector_info = model["detectors"]
    atmosphere_info = common["atmosphere"]
    dm_info = common["dm"]

    telescope = Telescope(
        resolution=Int(get(telescope_override, "resolution", telescope_info["resolution"])),
        diameter=T(telescope_info["diameter"]),
        sampling_time=T(telescope_info["sampling_time"]),
        central_obstruction=T(telescope_info["central_obstruction"]),
        T=T,
        backend=backend,
    )
    atmosphere = MultiLayerAtmosphere(
        telescope;
        r0=T(atmosphere_info["r0"]),
        L0=T(atmosphere_info["L0"]),
        fractional_cn2=T.(atmosphere_info["fractional_r0"]),
        wind_speed=T.(atmosphere_info["wind_speed"]),
        wind_direction=T.(atmosphere_info["wind_direction"]),
        altitude=T.(atmosphere_info["altitude"]),
        T=T,
        backend=backend,
    )
    dm = DeformableMirror(
        telescope;
        n_act=Int(dm_info["n_act"]),
        influence_width=T(dm_info["influence_width"]),
        misregistration=misregistration(telescope.params.resolution, dm_info; T=T),
        T=T,
        backend=backend,
    )
    wfs = wavefront_sensor(wfs_info, telescope; T=T, backend=backend)

    calibration_source = source_from_flux_scale(
        band=:R,
        magnitude=source_info["calibration_magnitude"],
        flux_scale=source_info["calibration_flux_scale"],
        T=T,
    )
    sky_source = source_from_flux_scale(
        band=:R,
        magnitude=source_info["star_magnitude_vis"],
        flux_scale=source_info["vis_flux_scale"],
        T=T,
    )
    science_source = source_from_flux_scale(
        band=:IR1310,
        wavelength=1.310e-6,
        magnitude=source_info["star_magnitude_ir"],
        flux_scale=1.0,
        T=T,
    )

    wfs_detector_name = Symbol(detector_info["wfs"])
    gain_detector_name = haskey(detector_info, "gain") ? Symbol(detector_info["gain"]) : nothing
    science_detector_name = Symbol(get(detector_info, "science", "cred2"))

    wfs_detector = detector_from_camera(cameras[wfs_detector_name];
        integration_time=telescope.params.sampling_time,
        T=T,
        backend=backend,
        response_mode=response_mode,
    )
    gain_detector = isnothing(gain_detector_name) ? nothing :
        detector_from_camera(cameras[gain_detector_name];
            integration_time=telescope.params.sampling_time,
            T=T,
            backend=backend,
            response_mode=response_mode,
        )
    science_detector = detector_from_camera(cameras[science_detector_name];
        integration_time=telescope.params.sampling_time,
        T=T,
        backend=backend,
        response_mode=response_mode,
    )

    return (
        name=kind,
        label=model_info["label"],
        wfs_family=Symbol(wfs_info["family"]),
        n_subap=Int(wfs_info["n_subap"]),
        wfs_detector_name=wfs_detector_name,
        gain_detector_name=gain_detector_name,
        science_detector_name=science_detector_name,
        telescope=telescope,
        atmosphere=atmosphere,
        dm=dm,
        wfs=wfs,
        calibration_source=calibration_source,
        sky_source=sky_source,
        science_source=science_source,
        wfs_detector=wfs_detector,
        gain_detector=gain_detector,
        science_detector=science_detector,
        camera_configs=cameras,
        rng=MersenneTwister(0),
    )
end
