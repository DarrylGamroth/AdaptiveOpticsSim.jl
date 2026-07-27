abstract type DeviceModelMatrixWFSFamily end
struct DeviceModelMatrixShackHartmann <: DeviceModelMatrixWFSFamily end
struct DeviceModelMatrixPyramid <: DeviceModelMatrixWFSFamily end
struct DeviceModelMatrixBioEdge <: DeviceModelMatrixWFSFamily end

@inline device_model_matrix_wfs_rows() = (
    (
        DeviceModelMatrixShackHartmann(),
        Val(:ngs),
        Val(:monochromatic),
    ),
    (
        DeviceModelMatrixShackHartmann(),
        Val(:lgs),
        Val(:monochromatic),
    ),
    (
        DeviceModelMatrixPyramid(),
        Val(:ngs),
        Val(:spectral),
    ),
    (
        DeviceModelMatrixBioEdge(),
        Val(:ngs),
        Val(:spectral),
    ),
)

struct DeviceModelMatrixWFSPathModel{F<:DeviceModelMatrixWFSFamily}
    family::F
end

Plant.plant_model_definition_style(
    ::Type{<:DeviceModelMatrixWFSPathModel},
) = ColdPlantModelDefinition()

@inline device_model_matrix_family_symbol(
    ::DeviceModelMatrixShackHartmann,
) = :shack_hartmann
@inline device_model_matrix_family_symbol(
    ::DeviceModelMatrixPyramid,
) = :pyramid
@inline device_model_matrix_family_symbol(
    ::DeviceModelMatrixBioEdge,
) = :bioedge

function device_model_matrix_front_end(
    ::DeviceModelMatrixShackHartmann,
    telescope::Telescope,
    source::AdaptiveOpticsSim.AbstractSource,
    ::Type{T},
) where {T<:AbstractFloat}
    sensor = ShackHartmannWFS(
        telescope;
        n_lenslets=2,
        n_pix_subap=2,
        mode=Diffractive(),
        T,
        backend=backend(telescope),
    )
    return ShackHartmannOpticalFrontEnd(sensor.front_end, source)
end

function device_model_matrix_front_end(
    ::DeviceModelMatrixPyramid,
    telescope::Telescope,
    source::AdaptiveOpticsSim.AbstractSource,
    ::Type{T},
) where {T<:AbstractFloat}
    sensor = PyramidWFS(
        telescope;
        pupil_samples=2,
        mode=Diffractive(),
        modulation=T(1),
        modulation_points=4,
        diffraction_padding=2,
        T,
        backend=backend(telescope),
    )
    return PyramidOpticalFrontEnd(sensor, source)
end

function device_model_matrix_front_end(
    ::DeviceModelMatrixBioEdge,
    telescope::Telescope,
    source::AdaptiveOpticsSim.AbstractSource,
    ::Type{T},
) where {T<:AbstractFloat}
    sensor = BioEdgeWFS(
        telescope;
        pupil_samples=2,
        mode=Diffractive(),
        modulation=T(1),
        modulation_points=4,
        diffraction_padding=2,
        T,
        backend=backend(telescope),
    )
    return BioEdgeOpticalFrontEnd(sensor, source)
end

@inline device_model_matrix_rate_map(
    ::DeviceModelMatrixShackHartmann,
    front_end::ShackHartmannOpticalFrontEnd,
    pupil::PupilFunction,
    source::SpectralSource,
) = shack_hartmann_rate_map(front_end, pupil, source)

@inline device_model_matrix_rate_map(
    ::DeviceModelMatrixShackHartmann,
    front_end::ShackHartmannOpticalFrontEnd,
    pupil::PupilFunction,
    ::AdaptiveOpticsSim.AbstractSource,
) = shack_hartmann_rate_map(front_end, pupil)

@inline device_model_matrix_rate_map(
    ::DeviceModelMatrixPyramid,
    front_end::PyramidOpticalFrontEnd,
    pupil::PupilFunction,
    ::AdaptiveOpticsSim.AbstractSource,
) = pyramid_rate_map(front_end, pupil)

@inline device_model_matrix_rate_map(
    ::DeviceModelMatrixBioEdge,
    front_end::BioEdgeOpticalFrontEnd,
    pupil::PupilFunction,
    ::AdaptiveOpticsSim.AbstractSource,
) = bioedge_rate_map(front_end, pupil)

function Plant.prepare_path_executor(
    model::DeviceModelMatrixWFSPathModel,
    definition::OpticalPathDefinition,
    source::AdaptiveOpticsSim.AbstractSource,
    telescope::Telescope,
    atmosphere::AdaptiveOpticsSim.AbstractTimedAtmosphere,
)
    T = eltype(pupil_reflectivity(telescope))
    pupil = PupilFunction(telescope; T, backend=backend(telescope))
    front_end = device_model_matrix_front_end(
        model.family, telescope, source, T)
    output = device_model_matrix_rate_map(
        model.family, front_end, pupil, source)
    plan = prepare_wfs_optical_formation(front_end, pupil, output)
    execution = WFSOpticalPathExecution(plan)
    family = device_model_matrix_family_symbol(model.family)
    return PreparedPathExecutor(
        definition,
        source,
        telescope,
        atmosphere,
        pupil,
        output,
        execution;
        materialization=prepare_pupil_opd_materialization(
            atmosphere,
            telescope,
            source,
            pupil,
        ),
        optical_model=(
            kind=:device_model_matrix_wfs,
            family,
        ),
        propagation_model=family,
        model_revisions=UInt(1),
    )
end

@inline device_model_matrix_source(
    ::Val{:ngs},
    wavelength_m,
    photon_rate,
    coordinates,
    ::Type{T},
) where {T<:AbstractFloat} = Source(
    band=:custom,
    wavelength=T(wavelength_m),
    photon_irradiance=T(photon_rate),
    coordinates=(T(coordinates[1]), T(coordinates[2])),
    T=T,
)

@inline device_model_matrix_source(
    ::Val{:lgs},
    wavelength_m,
    photon_rate,
    coordinates,
    ::Type{T},
) where {T<:AbstractFloat} = LGSSource(
    wavelength=T(wavelength_m),
    photon_irradiance=T(photon_rate),
    coordinates=(T(coordinates[1]), T(coordinates[2])),
    altitude=T(90_000),
    elongation_factor=T(1.4),
    T=T,
)

@inline device_model_matrix_wavelength(
    ::Val{:ngs},
    ::Type{T},
) where {T<:AbstractFloat} = T(0.75e-6)

@inline device_model_matrix_wavelength(
    ::Val{:lgs},
    ::Type{T},
) where {T<:AbstractFloat} = T(589e-9)

@inline device_model_matrix_spectrum(
    source::AdaptiveOpticsSim.AbstractSource,
    ::Val{:monochromatic},
    ::Type{T},
) where {T<:AbstractFloat} = source

function device_model_matrix_spectrum(
    source::AdaptiveOpticsSim.AbstractSource,
    ::Val{:spectral},
    ::Type{T},
) where {T<:AbstractFloat}
    return with_spectrum(
        source,
        SpectralBundle(
            T[T(0.94) * wavelength(source),
                T(1.06) * wavelength(source)],
            T[T(0.45), T(0.55)];
            T,
        ),
    )
end

function device_model_matrix_wfs_fixture(
    family::DeviceModelMatrixWFSFamily;
    backend::AdaptiveOpticsSim.AbstractArrayBackend=CPUBackend(),
    selection::Val=Val(:all),
    direction::Val=Val(:ngs),
    spectral::Val=Val(:monochromatic),
    T::Type{<:AbstractFloat}=Float64,
    r0::Real=T(0.2),
)
    telescope = Telescope(
        resolution=8,
        diameter=T(4),
        central_obstruction=zero(T),
        T=T,
        backend=backend,
    )
    atmosphere = MultiLayerAtmosphere(
        telescope;
        r0=T(r0),
        L0=T(25),
        fractional_cn2=T[0.65, 0.35],
        wind_speed=T[7, 11],
        wind_direction=T[20, 125],
        altitude=T[0, 5_000],
        layer_ids=(:ground, :high),
        T=T,
        backend=backend,
    )
    wavelength_m = device_model_matrix_wavelength(direction, T)
    first_leaf = device_model_matrix_source(
        direction, wavelength_m, T(55), (T(3), T(25)), T)
    second_leaf = device_model_matrix_source(
        direction, wavelength_m, T(45), (T(7), T(145)), T)
    first_source = device_model_matrix_spectrum(first_leaf, spectral, T)
    second_source = device_model_matrix_spectrum(second_leaf, spectral, T)
    witness_source = Source(
        band=:custom,
        wavelength=T(0.8e-6),
        photon_irradiance=T(30),
        coordinates=(zero(T), zero(T)),
        T=T,
    )
    path_definitions = (
        OpticalPathDefinition(
            :wfs_alpha,
            first_source,
            DeviceModelMatrixWFSPathModel(family),
        ),
        OpticalPathDefinition(
            :wfs_beta,
            second_source,
            DeviceModelMatrixWFSPathModel(family),
        ),
        OpticalPathDefinition(
            :readout_witness,
            witness_source,
            DeviceBatchTestPathModel(2),
        ),
    )
    acquisition_definitions = (
        AcquisitionDefinition(
            :readout_witness_camera,
            :readout_witness,
            DeviceBatchTestAcquisitionModel(T(0.08)),
        ),
    )
    definition = PlantDefinition(
        ;
        telescope,
        atmosphere,
        paths=path_definitions,
        acquisitions=acquisition_definitions,
    )
    plant = prepare_plant(definition; run_seed=0x7_500)
    event_definition = PlantEventLoopDefinition((
        OpticalSampleDefinition(
            :wfs_alpha,
            PeriodicSchedule(period_ns=100_000_000),
        ),
        OpticalSampleDefinition(
            :wfs_beta,
            PeriodicSchedule(period_ns=100_000_000),
        ),
        OpticalSampleDefinition(
            :readout_witness,
            PeriodicSchedule(period_ns=100_000_000),
        ),
    ), (
        AcquisitionEventDefinition(
            :readout_witness_camera,
            GlobalShutterAcquisitionDefinition(
                PlantDuration(80_000_000);
                readout_duration=PlantDuration(10_000_000),
                readiness_delay=PlantDuration(5_000_000),
            ),
            PeriodicAcquisitionStart(
                PeriodicSchedule(
                    period_ns=200_000_000,
                    phase_ns=10_000_000,
                ),
            ),
        ),
    ))
    prepared = Plant._prepare_plant_event_loop(
        plant,
        event_definition,
        selection,
    )
    return (
        plant,
        prepared,
        state=PlantEventLoopState(prepared),
        workspace=PlantEventLoopWorkspace(prepared),
        path_ids=(:wfs_alpha, :wfs_beta),
    )
end

@inline function device_model_matrix_product_host(
    product::IntensityMap,
)
    return Array(product.values)
end

function device_model_matrix_product_host(
    product::OpticalProductBundle,
)
    return ntuple(length(product)) do index
        Array(product[index].values)
    end
end

function device_model_matrix_copy_atmosphere_screens!(
    destination::MultiLayerAtmosphere,
    source::MultiLayerAtmosphere,
)
    length(destination.layers) == length(source.layers) ||
        throw(DimensionMismatch(
            "device-model matrix atmosphere layer counts differ",
        ))
    @inbounds for index in eachindex(destination.layers)
        copyto!(
            destination.layers[index].generator.state.opd,
            source.layers[index].generator.state.opd,
        )
    end
    AdaptiveOpticsSim.synchronize_backend!(
        AdaptiveOpticsSim.execution_style(
            first(destination.layers).generator.state.opd,
        ),
    )
    return destination
end
