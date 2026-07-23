module Gate2SerialPlantBenchmark

using AdaptiveOpticsSim
using AdaptiveOpticsSim.Plant
using SHA

const AOS = AdaptiveOpticsSim
const AOSPlant = AdaptiveOpticsSim.Plant

struct DirectSciencePathModel{R}
    zero_padding::Int
    revision::R
end

struct ShackHartmannPathModel{R}
    n_lenslets::Int
    n_pix_subap::Int
    revision::R
end

struct PyramidPathModel{T<:AbstractFloat,R}
    pupil_samples::Int
    modulation::T
    modulation_points::Int
    revision::R
end

struct FrameAcquisitionModel{T<:AbstractFloat}
    exposure_s::T
    quantum_efficiency::T
end

for model in (
    DirectSciencePathModel,
    ShackHartmannPathModel,
    PyramidPathModel,
    FrameAcquisitionModel,
)
    @eval AOSPlant.plant_model_definition_style(::Type{<:$model}) =
        AOSPlant.ColdPlantModelDefinition()
end

function AOSPlant.prepare_path_executor(
    model::DirectSciencePathModel,
    definition::AOSPlant.OpticalPathDefinition,
    source::AOS.AbstractSource,
    telescope::AOS.Telescope,
    atmosphere::AOS.AbstractTimedAtmosphere,
)
    pupil = AOS.PupilFunction(telescope)
    imaging = AOS.prepare_direct_imaging(pupil, source;
        zero_padding=model.zero_padding)
    return AOSPlant.PreparedPathExecutor(
        definition,
        source,
        telescope,
        atmosphere,
        pupil,
        AOS.direct_imaging_output(imaging),
        imaging;
        materialization=AOSPlant.prepare_pupil_opd_materialization(
            atmosphere, telescope, source, pupil),
        optical_model=(kind=:direct_imaging,
            zero_padding=model.zero_padding),
        propagation_model=:fraunhofer_fft,
        model_revisions=model.revision,
    )
end

function AOSPlant.prepare_path_executor(
    model::ShackHartmannPathModel,
    definition::AOSPlant.OpticalPathDefinition,
    source::AOS.AbstractSource,
    telescope::AOS.Telescope,
    atmosphere::AOS.AbstractTimedAtmosphere,
)
    T = eltype(AOS.pupil_reflectivity(telescope))
    pupil = AOS.PupilFunction(telescope; T=T)
    sensor = AOS.ShackHartmannWFS(telescope;
        n_lenslets=model.n_lenslets,
        n_pix_subap=model.n_pix_subap,
        mode=AOS.Diffractive(),
        T=T,
        backend=AOS.backend(telescope),
    )
    front_end = AOS.ShackHartmannOpticalFrontEnd(sensor.front_end, source)
    result = AOS.shack_hartmann_rate_map(front_end, pupil)
    plan = AOS.prepare_wfs_optical_formation(front_end, pupil, result)
    execution = AOSPlant.WFSOpticalPathExecution(plan)
    return AOSPlant.PreparedPathExecutor(
        definition,
        source,
        telescope,
        atmosphere,
        pupil,
        result,
        execution;
        materialization=AOSPlant.prepare_pupil_opd_materialization(
            atmosphere, telescope, source, pupil),
        optical_model=(kind=:shack_hartmann,
            n_lenslets=model.n_lenslets,
            n_pix_subap=model.n_pix_subap),
        propagation_model=:microlens_fraunhofer,
        model_revisions=(definition=model.revision,
            layout=AOS.subaperture_layout_revision(front_end.layout)),
    )
end

function AOSPlant.prepare_path_executor(
    model::PyramidPathModel,
    definition::AOSPlant.OpticalPathDefinition,
    source::AOS.AbstractSource,
    telescope::AOS.Telescope,
    atmosphere::AOS.AbstractTimedAtmosphere,
)
    T = eltype(AOS.pupil_reflectivity(telescope))
    pupil = AOS.PupilFunction(telescope; T=T)
    sensor = AOS.PyramidWFS(telescope;
        pupil_samples=model.pupil_samples,
        modulation=model.modulation,
        modulation_points=model.modulation_points,
        mode=AOS.Diffractive(),
        T=T,
        backend=AOS.backend(telescope),
    )
    front_end = AOS.PyramidOpticalFrontEnd(sensor, source)
    result = AOS.pyramid_rate_map(front_end, pupil)
    plan = AOS.prepare_wfs_optical_formation(front_end, pupil, result)
    execution = AOSPlant.WFSOpticalPathExecution(plan)
    return AOSPlant.PreparedPathExecutor(
        definition,
        source,
        telescope,
        atmosphere,
        pupil,
        result,
        execution;
        materialization=AOSPlant.prepare_pupil_opd_materialization(
            atmosphere, telescope, source, pupil),
        optical_model=(kind=:pyramid,
            pupil_samples=model.pupil_samples,
            modulation=model.modulation,
            modulation_points=model.modulation_points),
        propagation_model=:modulated_pyramid_fraunhofer,
        model_revisions=model.revision,
    )
end

function AOSPlant.prepare_acquisition_provider(
    model::FrameAcquisitionModel,
    ::AOSPlant.AcquisitionDefinition,
    path::AOSPlant.PreparedPathExecutor,
)
    AOSPlant.require_path_result(path)
    T = eltype(AOSPlant._first_path_result(path.result).values)
    detector = AOS.Detector(
        integration_time=T(model.exposure_s),
        noise=AOS.NoiseNone(),
        qe=T(model.quantum_efficiency),
        response_model=AOS.NullFrameResponse(),
        T=T,
        backend=path.key.backend,
    )
    execution = AOSPlant.FrameAcquisitionExecution(detector, path.result)
    metadata = (
        kind=:detector_frame,
        units=:detected_electrons,
        geometry=path.result.metadata,
        detector=AOS.detector_export_metadata(detector),
        semantics=:complete_acquisition,
    )
    products = AOSPlant.AcquisitionProducts(execution.observation; metadata)
    return AOSPlant.prepare_full_optical_provider(execution, products)
end

mutable struct SerialPlantOperation{S,T<:AbstractFloat}
    selection::S
    step_index::UInt64
    model_step_s::T
end

@inline function (operation::SerialPlantOperation)()
    operation.step_index += UInt64(1)
    next_model_time_s = operation.step_index * operation.model_step_s
    AOSPlant.execute_acquisition_selection_at!(operation.selection,
        next_model_time_s)
    return nothing
end

function serial_plant_definition(raw::AbstractDict;
    reverse_declarations::Bool=false)
    T = Float64
    resolution = Int(raw["resolution"])
    telescope = AOS.Telescope(
        resolution=resolution,
        diameter=T(raw["diameter_m"]),
        central_obstruction=T(raw["central_obstruction"]),
        T=T,
    )
    atmosphere = AOS.MultiLayerAtmosphere(telescope;
        r0=T(raw["r0_m"]),
        L0=T(raw["outer_scale_m"]),
        fractional_cn2=T.(raw["fractional_cn2"]),
        wind_speed=T.(raw["wind_speed_m_per_s"]),
        wind_direction=T.(raw["wind_direction_deg"]),
        altitude=T.(raw["layer_altitude_m"]),
        layer_ids=Tuple(Symbol.(raw["layer_ids"])),
        T=T,
    )

    science_source = AOS.Source(
        band=:custom,
        wavelength=T(raw["science_wavelength_m"]),
        photon_irradiance=T(raw["science_photon_irradiance"]),
        coordinates=(T(raw["science_radius_arcsec"]),
            T(raw["science_azimuth_deg"])),
        T=T,
    )
    ngs_source = AOS.Source(
        band=:custom,
        wavelength=T(raw["ngs_wavelength_m"]),
        photon_irradiance=T(raw["ngs_photon_irradiance"]),
        coordinates=(T(raw["ngs_radius_arcsec"]),
            T(raw["ngs_azimuth_deg"])),
        T=T,
    )
    lgs_source = AOS.LGSSource(
        wavelength=T(raw["lgs_wavelength_m"]),
        photon_irradiance=T(raw["lgs_photon_irradiance"]),
        coordinates=(T(raw["lgs_radius_arcsec"]),
            T(raw["lgs_azimuth_deg"])),
        altitude=T(raw["lgs_altitude_m"]),
        T=T,
    )

    paths = (
        AOSPlant.OpticalPathDefinition(:science, science_source,
            DirectSciencePathModel(Int(raw["science_zero_padding"]),
                UInt(1))),
        AOSPlant.OpticalPathDefinition(:ngs_shack_hartmann, ngs_source,
            ShackHartmannPathModel(Int(raw["sh_n_lenslets"]),
                Int(raw["sh_n_pix_subap"]), UInt(2))),
        AOSPlant.OpticalPathDefinition(:lgs_pyramid, lgs_source,
            PyramidPathModel(Int(raw["pyramid_pupil_samples"]),
                T(raw["pyramid_modulation_lambda_over_d"]),
                Int(raw["pyramid_modulation_points"]), UInt(3))),
    )
    acquisitions = (
        AOSPlant.AcquisitionDefinition(:science_fast, :science,
            FrameAcquisitionModel(T(raw["science_fast_exposure_s"]),
                T(raw["quantum_efficiency"]))),
        AOSPlant.AcquisitionDefinition(:science_slow, :science,
            FrameAcquisitionModel(T(raw["science_slow_exposure_s"]),
                T(raw["quantum_efficiency"]))),
        AOSPlant.AcquisitionDefinition(:ngs_frame, :ngs_shack_hartmann,
            FrameAcquisitionModel(T(raw["wfs_exposure_s"]),
                T(raw["quantum_efficiency"]))),
        AOSPlant.AcquisitionDefinition(:lgs_frame, :lgs_pyramid,
            FrameAcquisitionModel(T(raw["wfs_exposure_s"]),
                T(raw["quantum_efficiency"]))),
    )
    ordered_paths = reverse_declarations ? reverse(paths) : paths
    ordered_acquisitions = reverse_declarations ? reverse(acquisitions) :
        acquisitions
    return AOSPlant.PlantDefinition(; telescope, atmosphere,
        paths=ordered_paths, acquisitions=ordered_acquisitions)
end

function prepare_serial_plant_operation(raw::AbstractDict;
    reverse_declarations::Bool=false)
    definition = serial_plant_definition(raw; reverse_declarations)
    plant = AOSPlant.prepare_plant(definition;
        run_seed=UInt64(raw["run_seed"]),
        rng_derivation_version=Int(raw["rng_derivation_version"]),
    )
    selection = AOSPlant.prepare_acquisition_selection(plant,
        (:lgs_frame, :science_slow, :ngs_frame, :science_fast))
    model_step_s = Float64(raw["model_step_s"])
    return SerialPlantOperation(selection, UInt64(0), model_step_s)
end

@inline observation_values(values::AbstractArray) = values
@inline observation_values(observation::AOS.WFSObservation) =
    observation.storage

function observation_snapshot(operation::SerialPlantOperation)
    acquisitions = AOSPlant.prepared_acquisitions(operation.selection)
    return map(acquisitions) do owner
        values = observation_values(AOSPlant.acquisition_observation(owner))
        Array(values)
    end
end

function validate_replay_and_reordering(raw::AbstractDict)
    canonical = prepare_serial_plant_operation(raw)
    reordered = prepare_serial_plant_operation(raw;
        reverse_declarations=true)
    canonical()
    reordered()
    canonical_snapshot = observation_snapshot(canonical)
    reordered_snapshot = observation_snapshot(reordered)
    canonical_snapshot == reordered_snapshot || error(
        "Gate 2 serial plant changed outputs after declaration reordering")

    canonical_paths = AOSPlant.prepared_paths(canonical.selection)
    path_ids = collect(map(
        path -> String(AOSPlant.path_id(path.definition).name), canonical_paths))
    length(unique(path_ids)) == 3 || error(
        "Gate 2 serial plant must materialize three unique optical paths")
    science_fast = AOSPlant.prepared_acquisition(canonical.selection.plant,
        :science_fast)
    science_slow = AOSPlant.prepared_acquisition(canonical.selection.plant,
        :science_slow)
    science_fast.path_result === science_slow.path_result || error(
        "Gate 2 science acquisitions must reuse one optical result")
    exposure_ratio = Float64(raw["science_slow_exposure_s"]) /
        Float64(raw["science_fast_exposure_s"])
    fast = observation_values(AOSPlant.acquisition_observation(science_fast))
    slow = observation_values(AOSPlant.acquisition_observation(science_slow))
    isapprox(slow, exposure_ratio .* fast; rtol=1e-12, atol=1e-12) ||
        error("Gate 2 science detector fan-out changed exposure scaling")
    all(values -> all(isfinite, values) && sum(values) > 0,
        canonical_snapshot) || error(
        "Gate 2 serial plant produced an invalid acquisition product")

    return Dict{String,Any}(
        "declaration_reordering_replay" => true,
        "unique_path_count" => length(canonical_paths),
        "acquisition_count" => length(
            AOSPlant.prepared_acquisitions(canonical.selection)),
        "canonical_path_ids" => path_ids,
        "canonical_acquisition_ids" => collect(map(owner -> String(
                AOSPlant.acquisition_id(owner.definition).name),
            AOSPlant.prepared_acquisitions(canonical.selection))),
        "science_path_reused" => true,
        "science_exposure_ratio" => exposure_ratio,
        "observation_sha256" => collect(map(canonical_snapshot) do values
            bytes2hex(SHA.sha256(reinterpret(UInt8, vec(values))))
        end),
    )
end

function final_observation_summary(operation::SerialPlantOperation)
    acquisitions = AOSPlant.prepared_acquisitions(operation.selection)
    atmosphere = AOSPlant.plant_atmosphere(operation.selection.plant.definition)
    return Dict{String,Any}(
        "model_time_s" => AOS.epoch_time(AOS.current_epoch(atmosphere)),
        "epoch_sequence" => Int(AOS.epoch_sequence(
            AOS.current_epoch(atmosphere))),
        "acquisitions" => collect(map(acquisitions) do owner
            values = Array(observation_values(
                AOSPlant.acquisition_observation(owner)))
            Dict{String,Any}(
                "id" => String(AOSPlant.acquisition_id(owner.definition).name),
                "shape" => collect(size(values)),
                "sum" => sum(values),
                "minimum" => minimum(values),
                "maximum" => maximum(values),
                "sha256" => bytes2hex(SHA.sha256(
                    reinterpret(UInt8, vec(values)))),
            )
        end),
    )
end

export SerialPlantOperation
export prepare_serial_plant_operation
export validate_replay_and_reordering
export final_observation_summary

end
