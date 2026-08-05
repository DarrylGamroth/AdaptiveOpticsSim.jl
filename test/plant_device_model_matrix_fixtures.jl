abstract type DeviceModelMatrixWFSFamily end
struct DeviceModelMatrixShackHartmann <: DeviceModelMatrixWFSFamily end
struct DeviceModelMatrixPyramid <: DeviceModelMatrixWFSFamily end
struct DeviceModelMatrixBiOEdge <: DeviceModelMatrixWFSFamily end
struct DeviceModelMatrixZernike <: DeviceModelMatrixWFSFamily end
struct DeviceModelMatrixCurvature <: DeviceModelMatrixWFSFamily end

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
        DeviceModelMatrixBiOEdge(),
        Val(:ngs),
        Val(:spectral),
    ),
)

struct DeviceModelMatrixWFSPathModel{F<:DeviceModelMatrixWFSFamily}
    family::F
    variant::Int
end

DeviceModelMatrixWFSPathModel(family::DeviceModelMatrixWFSFamily) =
    DeviceModelMatrixWFSPathModel(family, 0)

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
    ::DeviceModelMatrixBiOEdge,
) = :bi_o_edge
@inline device_model_matrix_family_symbol(
    ::DeviceModelMatrixZernike,
) = :zernike
@inline device_model_matrix_family_symbol(
    ::DeviceModelMatrixCurvature,
) = :curvature

function device_model_matrix_front_end(
    ::DeviceModelMatrixShackHartmann,
    telescope::Telescope,
    source::AdaptiveOpticsSim.Optics.AbstractSource,
    ::Type{T},
    variant::Int=0,
) where {T<:AbstractFloat}
    sensor = ShackHartmannWFS(
        telescope;
        n_lenslets=iszero(variant) ? 2 : 4,
        n_pix_subap=2,
        mode=Diffractive(),
        T,
        backend=backend(telescope),
    )
    return shack_hartmann_optics(sensor, source)
end

function device_model_matrix_front_end(
    ::DeviceModelMatrixPyramid,
    telescope::Telescope,
    source::AdaptiveOpticsSim.Optics.AbstractSource,
    ::Type{T},
    ::Int=0,
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
    ::DeviceModelMatrixBiOEdge,
    telescope::Telescope,
    source::AdaptiveOpticsSim.Optics.AbstractSource,
    ::Type{T},
    ::Int=0,
) where {T<:AbstractFloat}
    sensor = BiOEdgeWFS(
        telescope;
        pupil_samples=2,
        mode=Diffractive(),
        modulation=T(1),
        modulation_points=4,
        diffraction_padding=2,
        T,
        backend=backend(telescope),
    )
    return BiOEdgeOpticalFrontEnd(sensor, source)
end

function device_model_matrix_front_end(
    ::DeviceModelMatrixZernike,
    telescope::Telescope,
    source::AdaptiveOpticsSim.Optics.AbstractSource,
    ::Type{T},
    ::Int=0,
) where {T<:AbstractFloat}
    sensor = ZernikeWFS(
        telescope;
        pupil_samples=4,
        T,
        backend=backend(telescope),
    )
    return ZernikeOpticalFrontEnd(sensor, source)
end

function device_model_matrix_front_end(
    ::DeviceModelMatrixCurvature,
    telescope::Telescope,
    source::AdaptiveOpticsSim.Optics.AbstractSource,
    ::Type{T},
    ::Int=0,
) where {T<:AbstractFloat}
    sensor = CurvatureWFS(
        telescope;
        pupil_samples=4,
        T,
        backend=backend(telescope),
    )
    return CurvatureOpticalFrontEnd(sensor, source)
end

@inline device_model_matrix_rate_map(
    ::DeviceModelMatrixShackHartmann,
    optics::WavefrontSensors.ShackHartmannOptics,
    pupil::PupilFunction,
    source::SpectralSource,
) = shack_hartmann_rate_map(optics, pupil, source)

@inline device_model_matrix_rate_map(
    ::DeviceModelMatrixShackHartmann,
    optics::WavefrontSensors.ShackHartmannOptics,
    pupil::PupilFunction,
    ::AdaptiveOpticsSim.Optics.AbstractSource,
) = shack_hartmann_rate_map(optics, pupil)

@inline device_model_matrix_rate_map(
    ::DeviceModelMatrixPyramid,
    front_end::PyramidOpticalFrontEnd,
    pupil::PupilFunction,
    ::AdaptiveOpticsSim.Optics.AbstractSource,
) = pyramid_rate_map(front_end, pupil)

@inline device_model_matrix_rate_map(
    ::DeviceModelMatrixBiOEdge,
    front_end::BiOEdgeOpticalFrontEnd,
    pupil::PupilFunction,
    ::AdaptiveOpticsSim.Optics.AbstractSource,
) = bi_o_edge_rate_map(front_end, pupil)

@inline device_model_matrix_rate_map(
    ::DeviceModelMatrixZernike,
    front_end::ZernikeOpticalFrontEnd,
    pupil::PupilFunction,
    ::AdaptiveOpticsSim.Optics.AbstractSource,
) = zernike_rate_map(front_end, pupil)

@inline device_model_matrix_rate_map(
    ::DeviceModelMatrixCurvature,
    front_end::CurvatureOpticalFrontEnd,
    pupil::PupilFunction,
    ::AdaptiveOpticsSim.Optics.AbstractSource,
) = curvature_rate_maps(front_end, pupil)

function Plant.prepare_path_executor(
    model::DeviceModelMatrixWFSPathModel,
    definition::OpticalPathDefinition,
    source::AdaptiveOpticsSim.Optics.AbstractSource,
    telescope::Telescope,
    atmosphere::AdaptiveOpticsSim.Atmospheres.AbstractTimedAtmosphere,
    context,
)
    T = eltype(pupil_reflectivity(telescope))
    pupil = PupilFunction(telescope; T, backend=backend(telescope))
    front_end = device_model_matrix_front_end(
        model.family, telescope, source, T, model.variant)
    output = device_model_matrix_rate_map(
        model.family, front_end, pupil, source)
    plan = prepare_wfs_optics(front_end, pupil, output)
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
        context=context,
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
    source::AdaptiveOpticsSim.Optics.AbstractSource,
    ::Val{:monochromatic},
    ::Type{T},
) where {T<:AbstractFloat} = source

function device_model_matrix_spectrum(
    source::AdaptiveOpticsSim.Optics.AbstractSource,
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
    backend::AdaptiveOpticsSim.Backends.AbstractArrayBackend=CPUBackend(),
    selection::Val=Val(:all),
    direction::Val=Val(:ngs),
    spectral::Val=Val(:monochromatic),
    include_second::Bool=true,
    second_family::DeviceModelMatrixWFSFamily=family,
    second_direction::Val=direction,
    second_spectral::Val=spectral,
    second_variant::Int=0,
    second_period_ns::Integer=100_000_000,
    second_phase_ns::Integer=0,
    second_origin::PlantTimestamp=zero(PlantTimestamp),
    include_physical_state::Bool=true,
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
        wind_direction_deg=T[20, 125],
        altitude=T[0, 5_000],
        layer_ids=(:ground, :high),
        T=T,
        backend=backend,
    )
    wavelength_m = device_model_matrix_wavelength(direction, T)
    second_wavelength_m = device_model_matrix_wavelength(second_direction, T)
    first_leaf = device_model_matrix_source(
        direction, wavelength_m, T(55), (T(3), T(25)), T)
    second_leaf = device_model_matrix_source(
        second_direction, second_wavelength_m, T(45), (T(7), T(145)), T)
    first_source = device_model_matrix_spectrum(first_leaf, spectral, T)
    second_source =
        device_model_matrix_spectrum(second_leaf, second_spectral, T)
    witness_source = Source(
        band=:custom,
        wavelength=T(0.8e-6),
        photon_irradiance=T(30),
        coordinates=(zero(T), zero(T)),
        T=T,
    )
    first_path = OpticalPathDefinition(
        :wfs_alpha,
        first_source,
        DeviceModelMatrixWFSPathModel(family),
    )
    second_path = OpticalPathDefinition(
        :wfs_beta,
        second_source,
        DeviceModelMatrixWFSPathModel(second_family, second_variant),
    )
    witness_path = OpticalPathDefinition(
        :readout_witness,
        witness_source,
        DeviceBatchTestPathModel(2),
    )
    path_definitions = include_second ?
        (first_path, second_path, witness_path) :
        (first_path, witness_path)
    acquisition_definitions = (
        AcquisitionDefinition(
            :readout_witness_camera,
            :readout_witness,
            DeviceBatchTestAcquisitionModel(T(0.08)),
        ),
    )
    physical = include_physical_state ?
        device_batch_test_physical_definitions(
            telescope,
            backend,
            T;
            selected_path=:wfs_alpha,
        ) : nothing
    definition = PlantDefinition(
        ;
        telescope=AdaptiveOpticsSim.Optics.TelescopeDefinition(
            resolution=8,
            diameter=T(4),
            central_obstruction=zero(T),
            revision=1,
            T=T,
        ),
        atmosphere=MultiLayerAtmosphereDefinition(
            r0=T(r0),
            L0=T(25),
            fractional_cn2=T[0.65, 0.35],
            wind_speed=T[7, 11],
            wind_direction_deg=T[20, 125],
            altitude=T[0, 5_000],
            layer_ids=(:ground, :high),
            T=T,
        ),
        sampled_aberrations=isnothing(physical) ? () :
            physical.sampled_aberrations,
        controllable_optics=isnothing(physical) ? () : (physical.optic,),
        paths=path_definitions,
        acquisitions=acquisition_definitions,
    )
    plant = prepare_plant(
        definition, compute_device(pupil_reflectivity(telescope));
        run_seed=0x7_500,
        command_endpoints=isnothing(physical) ? () :
            (physical.configuration,),
    )
    first_sample = OpticalSampleDefinition(
        :wfs_alpha,
        PeriodicSchedule(period_ns=100_000_000),
    )
    second_sample = OpticalSampleDefinition(
        :wfs_beta,
        PeriodicSchedule(
            period_ns=second_period_ns,
            phase_ns=second_phase_ns,
        ),
        origin=second_origin,
    )
    witness_sample = OpticalSampleDefinition(
        :readout_witness,
        PeriodicSchedule(period_ns=100_000_000),
    )
    sample_definitions = include_second ?
        (first_sample, second_sample, witness_sample) :
        (first_sample, witness_sample)
    event_definition = PlantEventLoopDefinition(sample_definitions, (
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
        path_ids=include_second ?
            (:wfs_alpha, :wfs_beta) : (:wfs_alpha,),
        command_schema=isnothing(physical) ? nothing : physical.schema,
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
    AdaptiveOpticsSim.Backends.synchronize_backend!(
        AdaptiveOpticsSim.Backends.execution_style(
            first(destination.layers).generator.state.opd,
        ),
    )
    return destination
end

abstract type DeviceModelMatrixDetectorRow end

struct DeviceModelMatrixM1CCD <: DeviceModelMatrixDetectorRow end
struct DeviceModelMatrixM2FrameTransferEMCCD <:
    DeviceModelMatrixDetectorRow end
struct DeviceModelMatrixM3GlobalCMOS <: DeviceModelMatrixDetectorRow end
struct DeviceModelMatrixM4GlobalHgCdTe <: DeviceModelMatrixDetectorRow end
struct DeviceModelMatrixM5RollingCMOS <: DeviceModelMatrixDetectorRow end
struct DeviceModelMatrixM6UpTheRampHgCdTe <:
    DeviceModelMatrixDetectorRow end

@inline device_model_matrix_wfs_detector_row(
    ::DeviceModelMatrixShackHartmann,
    ::Val{:ngs},
) = DeviceModelMatrixM1CCD()

@inline device_model_matrix_wfs_detector_row(
    ::DeviceModelMatrixShackHartmann,
    ::Val{:lgs},
) = DeviceModelMatrixM2FrameTransferEMCCD()

@inline device_model_matrix_wfs_detector_row(
    ::DeviceModelMatrixPyramid,
    ::Val{:ngs},
) = DeviceModelMatrixM3GlobalCMOS()

@inline device_model_matrix_wfs_detector_row(
    ::DeviceModelMatrixBiOEdge,
    ::Val{:ngs},
) = DeviceModelMatrixM4GlobalHgCdTe()

@inline device_model_matrix_detector_rows() = (
    DeviceModelMatrixM1CCD(),
    DeviceModelMatrixM2FrameTransferEMCCD(),
    DeviceModelMatrixM3GlobalCMOS(),
    DeviceModelMatrixM4GlobalHgCdTe(),
    DeviceModelMatrixM5RollingCMOS(),
    DeviceModelMatrixM6UpTheRampHgCdTe(),
)

@inline device_model_matrix_detector_row_id(
    ::DeviceModelMatrixM1CCD,
) = :M1
@inline device_model_matrix_detector_row_id(
    ::DeviceModelMatrixM2FrameTransferEMCCD,
) = :M2
@inline device_model_matrix_detector_row_id(
    ::DeviceModelMatrixM3GlobalCMOS,
) = :M3
@inline device_model_matrix_detector_row_id(
    ::DeviceModelMatrixM4GlobalHgCdTe,
) = :M4
@inline device_model_matrix_detector_row_id(
    ::DeviceModelMatrixM5RollingCMOS,
) = :M5
@inline device_model_matrix_detector_row_id(
    ::DeviceModelMatrixM6UpTheRampHgCdTe,
) = :M6

@inline device_model_matrix_detector_sensor_symbol(
    ::DeviceModelMatrixM1CCD,
) = :ccd
@inline device_model_matrix_detector_sensor_symbol(
    ::DeviceModelMatrixM2FrameTransferEMCCD,
) = :emccd
@inline device_model_matrix_detector_sensor_symbol(
    ::Union{DeviceModelMatrixM3GlobalCMOS,DeviceModelMatrixM5RollingCMOS},
) = :cmos
@inline device_model_matrix_detector_sensor_symbol(
    ::Union{DeviceModelMatrixM4GlobalHgCdTe,
        DeviceModelMatrixM6UpTheRampHgCdTe},
) = :hgcdte

@inline device_model_matrix_detector_response_symbol(
    ::DeviceModelMatrixM1CCD,
) = :rectangular_aperture
@inline device_model_matrix_detector_response_symbol(
    ::DeviceModelMatrixM2FrameTransferEMCCD,
) = :gaussian
@inline device_model_matrix_detector_response_symbol(
    ::Union{DeviceModelMatrixM3GlobalCMOS,DeviceModelMatrixM4GlobalHgCdTe},
) = :none
@inline device_model_matrix_detector_response_symbol(
    ::DeviceModelMatrixM5RollingCMOS,
) = :sampled
@inline device_model_matrix_detector_response_symbol(
    ::DeviceModelMatrixM6UpTheRampHgCdTe,
) = :rectangular_aperture

@inline device_model_matrix_detector_timing_symbol(
    ::Union{DeviceModelMatrixM1CCD,DeviceModelMatrixM2FrameTransferEMCCD,
        DeviceModelMatrixM3GlobalCMOS,DeviceModelMatrixM4GlobalHgCdTe,
        DeviceModelMatrixM6UpTheRampHgCdTe},
) = :global_shutter
@inline device_model_matrix_detector_timing_symbol(
    ::DeviceModelMatrixM5RollingCMOS,
) = :rolling_shutter

@inline device_model_matrix_detector_sampling_symbol(
    ::Union{DeviceModelMatrixM1CCD,DeviceModelMatrixM2FrameTransferEMCCD,
        DeviceModelMatrixM3GlobalCMOS,DeviceModelMatrixM4GlobalHgCdTe,
        DeviceModelMatrixM5RollingCMOS},
) = :single_read
@inline device_model_matrix_detector_sampling_symbol(
    ::DeviceModelMatrixM6UpTheRampHgCdTe,
) = :up_the_ramp

@inline device_model_matrix_detector_acquisition_symbol(
    ::DeviceModelMatrixM2FrameTransferEMCCD,
) = :frame_transfer
@inline device_model_matrix_detector_acquisition_symbol(
    ::DeviceModelMatrixDetectorRow,
) = :standard

function device_model_matrix_detector_response(
    ::DeviceModelMatrixM1CCD,
    backend::AdaptiveOpticsSim.Backends.AbstractArrayBackend,
    ::Type{T},
) where {T<:AbstractFloat}
    return RectangularPixelAperture(
        ;
        pitch_x_px=T(2),
        pitch_y_px=T(2),
        fill_factor_x=T(0.6),
        fill_factor_y=T(0.8),
        T,
        backend,
    )
end

function device_model_matrix_detector_response(
    ::DeviceModelMatrixM2FrameTransferEMCCD,
    backend::AdaptiveOpticsSim.Backends.AbstractArrayBackend,
    ::Type{T},
) where {T<:AbstractFloat}
    return GaussianPixelResponse(
        ;
        response_width_px=T(0.75),
        T,
        backend,
    )
end

@inline device_model_matrix_detector_response(
    ::Union{DeviceModelMatrixM3GlobalCMOS,DeviceModelMatrixM4GlobalHgCdTe},
    ::AdaptiveOpticsSim.Backends.AbstractArrayBackend,
    ::Type{<:AbstractFloat},
) = NullFrameResponse()

function device_model_matrix_detector_response(
    ::DeviceModelMatrixM5RollingCMOS,
    backend::AdaptiveOpticsSim.Backends.AbstractArrayBackend,
    ::Type{T},
) where {T<:AbstractFloat}
    return SampledFrameResponse(
        T[0 0.1 0; 0.1 0.6 0.1; 0 0.1 0];
        T,
        backend,
    )
end

function device_model_matrix_detector_response(
    ::DeviceModelMatrixM6UpTheRampHgCdTe,
    backend::AdaptiveOpticsSim.Backends.AbstractArrayBackend,
    ::Type{T},
) where {T<:AbstractFloat}
    return RectangularPixelAperture(
        ;
        pitch_x_px=T(2),
        pitch_y_px=T(2),
        fill_factor_x=T(0.7),
        fill_factor_y=T(0.9),
        T,
        backend,
    )
end

@inline device_model_matrix_detector_sensor(
    ::DeviceModelMatrixM1CCD,
    ::AdaptiveOpticsSim.Backends.AbstractArrayBackend,
    ::Type{T},
) where {T<:AbstractFloat} = CCDSensor(
    sampling_mode=SingleRead(),
    T=T,
)

@inline device_model_matrix_detector_sensor(
    ::DeviceModelMatrixM2FrameTransferEMCCD,
    ::AdaptiveOpticsSim.Backends.AbstractArrayBackend,
    ::Type{T},
) where {T<:AbstractFloat} = EMCCDSensor(
    excess_noise_factor=one(T),
    acquisition_mode=FrameTransferAcquisition(transfer_duration=T(0.1)),
    T=T,
)

@inline device_model_matrix_detector_sensor(
    ::DeviceModelMatrixM3GlobalCMOS,
    backend::AdaptiveOpticsSim.Backends.AbstractArrayBackend,
    ::Type{T},
) where {T<:AbstractFloat} = CMOSSensor(
    timing_model=GlobalShutter(),
    T=T,
    backend=backend,
)

@inline device_model_matrix_detector_sensor(
    ::DeviceModelMatrixM4GlobalHgCdTe,
    ::AdaptiveOpticsSim.Backends.AbstractArrayBackend,
    ::Type{T},
) where {T<:AbstractFloat} = HgCdTeSensor(
    sampling_mode=SingleRead(),
    T=T,
)

@inline device_model_matrix_detector_sensor(
    ::DeviceModelMatrixM5RollingCMOS,
    backend::AdaptiveOpticsSim.Backends.AbstractArrayBackend,
    ::Type{T},
) where {T<:AbstractFloat} = CMOSSensor(
    timing_model=RollingShutter(T(0.1); row_group_size=3),
    T=T,
    backend=backend,
)

@inline device_model_matrix_detector_sensor(
    ::DeviceModelMatrixM6UpTheRampHgCdTe,
    ::AdaptiveOpticsSim.Backends.AbstractArrayBackend,
    ::Type{T},
) where {T<:AbstractFloat} = HgCdTeSensor(
    read_duration=zero(T),
    sampling_mode=UpTheRampSampling(3),
    T=T,
)

function device_model_matrix_detector(
    row::DeviceModelMatrixDetectorRow;
    backend::AdaptiveOpticsSim.Backends.AbstractArrayBackend=CPUBackend(),
    T::Type{<:AbstractFloat}=Float64,
)
    return Detector(
        ;
        exposure_duration=one(T),
        qe=T(0.5),
        gain=one(T),
        noise=NoiseNone(),
        sensor=device_model_matrix_detector_sensor(row, backend, T),
        response_model=device_model_matrix_detector_response(row, backend, T),
        T,
        backend,
    )
end

function device_model_matrix_detector_rate_host(
    ::DeviceModelMatrixDetectorRow,
    ::Type{T},
) where {T<:AbstractFloat}
    values = zeros(T, 9, 9)
    values[5, 5] = T(4)
    return values
end

function device_model_matrix_detector_rate_map(
    row::DeviceModelMatrixDetectorRow;
    backend::AdaptiveOpticsSim.Backends.AbstractArrayBackend=CPUBackend(),
    T::Type{<:AbstractFloat}=Float64,
)
    host = device_model_matrix_detector_rate_host(row, T)
    array_backend = AdaptiveOpticsSim.Backends._resolve_array_backend(backend)
    values = array_backend{T}(undef, size(host)...)
    copyto!(values, host)
    metadata = OpticalPlaneMetadata(
        DetectorPlane(),
        values;
        coordinate_domain=AngularCoordinates(),
        sampling=(one(T), one(T)),
        normalization=PhotonRateNormalization(),
        spatial_measure=CellIntegratedMeasure(),
        coherence=IncoherentIntensityAddition(),
        spectral=MonochromaticChannel(T(0.8e-6)),
    )
    return IntensityMap(metadata, values)
end

@inline device_model_matrix_response_kernel(
    ::NullFrameResponse,
    ::Type{T},
) where {T<:AbstractFloat} = reshape(T[one(T)], 1, 1)

function device_model_matrix_response_kernel(
    response::GaussianPixelResponse,
    ::Type{T},
) where {T<:AbstractFloat}
    kernel = T.(Array(response.kernel))
    return kernel * transpose(kernel)
end

@inline device_model_matrix_response_kernel(
    response::SampledFrameResponse,
    ::Type{T},
) where {T<:AbstractFloat} = T.(Array(response.kernel))

function device_model_matrix_response_kernel(
    response::RectangularPixelAperture,
    ::Type{T},
) where {T<:AbstractFloat}
    kernel_y = T.(Array(response.kernel_y))
    kernel_x = T.(Array(response.kernel_x))
    return kernel_y * transpose(kernel_x)
end

function device_model_matrix_zero_extended_response(
    input::AbstractMatrix{T},
    kernel::AbstractMatrix{T},
) where {T<:AbstractFloat}
    output = zeros(T, size(input))
    radius_i = fld(size(kernel, 1), 2)
    radius_j = fld(size(kernel, 2), 2)
    @inbounds for j in axes(output, 2), i in axes(output, 1)
        value = zero(T)
        for kj in axes(kernel, 2), ki in axes(kernel, 1)
            ii = i - (ki - radius_i - 1)
            jj = j - (kj - radius_j - 1)
            if checkbounds(Bool, input, ii, jj)
                value = muladd(kernel[ki, kj], input[ii, jj], value)
            end
        end
        output[i, j] = value
    end
    return output
end

function device_model_matrix_expected_detector_frame(
    row::DeviceModelMatrixDetectorRow,
    detector::Detector,
    ::Type{T},
) where {T<:AbstractFloat}
    input = device_model_matrix_detector_rate_host(row, T)
    return device_model_matrix_expected_detector_frame(
        row,
        detector,
        input,
        T,
    )
end

function device_model_matrix_expected_detector_frame(
    ::DeviceModelMatrixDetectorRow,
    detector::Detector,
    input::AbstractMatrix,
    ::Type{T},
) where {T<:AbstractFloat}
    kernel = device_model_matrix_response_kernel(
        detector.params.response_model,
        T,
    )
    response = device_model_matrix_zero_extended_response(T.(input), kernel)
    return response .* T(0.5)
end

function device_model_matrix_expected_detector_frame(
    row::DeviceModelMatrixM6UpTheRampHgCdTe,
    detector::Detector,
    ::Type{T},
) where {T<:AbstractFloat}
    input = device_model_matrix_detector_rate_host(row, T)
    kernel = device_model_matrix_response_kernel(
        detector.params.response_model,
        T,
    )
    # Three equally spaced reads see rates 1x then 3x. Their fitted slope is
    # the response to the mean rate, and the integrated product spans 1 s.
    return device_model_matrix_zero_extended_response(input, kernel)
end

function device_model_matrix_expected_detector_mtf(
    detector::Detector,
    spatial_frequency_x::Real,
    spatial_frequency_y::Real,
    ::Type{T},
) where {T<:AbstractFloat}
    kernel = device_model_matrix_response_kernel(
        detector.params.response_model,
        T,
    )
    center_i = fld(size(kernel, 1), 2) + 1
    center_j = fld(size(kernel, 2), 2) + 1
    transfer = zero(Complex{T})
    @inbounds for j in axes(kernel, 2), i in axes(kernel, 1)
        phase = -T(2pi) * (
            T(spatial_frequency_y) * T(i - center_i) +
            T(spatial_frequency_x) * T(j - center_j)
        )
        transfer += kernel[i, j] * cis(phase)
    end
    return abs(transfer) / sum(kernel)
end

@inline function device_model_matrix_response_metadata_signature(metadata)
    return (
        metadata.frame_response,
        metadata.response_width_px,
        metadata.response_application_domain,
        metadata.response_is_separable,
        metadata.response_is_shift_invariant,
        metadata.response_support_rows,
        metadata.response_support_cols,
        metadata.pitch_x_px,
        metadata.pitch_y_px,
        metadata.fill_factor_x,
        metadata.fill_factor_y,
        metadata.aperture_shape,
    )
end

@inline device_model_matrix_detector_start() =
    PlantTimestamp(100_000_000)
@inline device_model_matrix_detector_middle(start::PlantTimestamp) =
    start + PlantDuration(400_000_000)
@inline device_model_matrix_detector_ramp_middle(start::PlantTimestamp) =
    start + PlantDuration(500_000_000)
@inline device_model_matrix_detector_close(start::PlantTimestamp) =
    start + PlantDuration(1_000_000_000)

function device_model_matrix_execute_detector(
    row::Union{DeviceModelMatrixM1CCD,DeviceModelMatrixM3GlobalCMOS,
        DeviceModelMatrixM4GlobalHgCdTe};
    backend::AdaptiveOpticsSim.Backends.AbstractArrayBackend=CPUBackend(),
    T::Type{<:AbstractFloat}=Float64,
    input_map=nothing,
)
    detector = device_model_matrix_detector(row; backend, T)
    map = isnothing(input_map) ?
        device_model_matrix_detector_rate_map(row; backend, T) : input_map
    definition = GlobalShutterAcquisitionDefinition(
        PlantDuration(1_000_000_000);
        readout_duration=PlantDuration(200_000_000),
        readiness_delay=PlantDuration(100_000_000),
    )
    prepared = prepare_global_shutter_acquisition(detector, map, definition)
    metadata_before = detector_export_metadata(detector)
    state = GlobalShutterAcquisitionState(prepared)
    ready_status = detector_acquisition_status(state)
    start = device_model_matrix_detector_start()
    middle = device_model_matrix_detector_middle(start)
    close = device_model_matrix_detector_close(start)
    rng = Xoshiro(0x7_501)

    begin_exposure!(prepared, state, start)
    active_status = detector_acquisition_status(state)
    accumulate_exposure_interval!(prepared, state, start, middle, rng)
    accumulate_exposure_interval!(prepared, state, middle, close, rng)
    close_exposure!(prepared, state, close)
    pending_status = detector_acquisition_status(state)
    readout = close + definition.readout_duration
    output = complete_readout!(prepared, state, readout, rng)
    complete_status = detector_acquisition_status(state)
    readiness = readout + definition.readiness_delay
    mark_acquisition_ready!(prepared, state, readiness)
    final_status = detector_acquisition_status(state)

    return (
        row,
        detector,
        map,
        prepared,
        state,
        output,
        rng,
        metadata_before,
        metadata_after=detector_export_metadata(detector),
        status_trace=(
            ready_status,
            active_status,
            pending_status,
            complete_status,
            final_status,
        ),
        event_times=(; start, close, readout, readiness),
    )
end

function device_model_matrix_execute_detector(
    row::DeviceModelMatrixM2FrameTransferEMCCD;
    backend::AdaptiveOpticsSim.Backends.AbstractArrayBackend=CPUBackend(),
    T::Type{<:AbstractFloat}=Float64,
    input_map=nothing,
)
    detector = device_model_matrix_detector(row; backend, T)
    map = isnothing(input_map) ?
        device_model_matrix_detector_rate_map(row; backend, T) : input_map
    definition = FrameTransferAcquisitionDefinition(
        PlantDuration(1_000_000_000);
        readout_duration=PlantDuration(200_000_000),
    )
    prepared = prepare_frame_transfer_acquisition(detector, map, definition)
    metadata_before = detector_export_metadata(detector)
    state = FrameTransferAcquisitionState(prepared)
    initial_transition = (
        frame_transfer_image_ready(state),
        frame_transfer_storage_empty(state),
        frame_transfer_readout_pending(state),
    )
    start = device_model_matrix_detector_start()
    middle = device_model_matrix_detector_middle(start)
    close = device_model_matrix_detector_close(start)
    rng = Xoshiro(0x7_502)

    begin_exposure!(prepared, state, start)
    active_transition = (
        frame_transfer_image_ready(state),
        frame_transfer_storage_empty(state),
        frame_transfer_readout_pending(state),
    )
    accumulate_exposure_interval!(prepared, state, start, middle, rng)
    accumulate_exposure_interval!(prepared, state, middle, close, rng)
    close_exposure!(prepared, state, close)
    closed_transition = (
        frame_transfer_image_ready(state),
        frame_transfer_storage_empty(state),
        frame_transfer_readout_pending(state),
    )
    transfer = frame_transfer_complete_timestamp(state)
    complete_frame_transfer!(prepared, state, transfer)
    transferred_transition = (
        frame_transfer_image_ready(state),
        frame_transfer_storage_empty(state),
        frame_transfer_readout_pending(state),
    )
    readout = readout_complete_timestamp(state)
    output = complete_readout!(prepared, state, readout, rng)
    final_transition = (
        frame_transfer_image_ready(state),
        frame_transfer_storage_empty(state),
        frame_transfer_readout_pending(state),
    )

    return (
        row,
        detector,
        map,
        prepared,
        state,
        output,
        rng,
        metadata_before,
        metadata_after=detector_export_metadata(detector),
        status_trace=(
            initial_transition,
            active_transition,
            closed_transition,
            transferred_transition,
            final_transition,
        ),
        event_times=(; start, close, transfer, readout),
    )
end

@inline device_model_matrix_detector_facing_product(
    product::IntensityMap,
) = product

@inline device_model_matrix_detector_facing_product(
    bundle::OpticalProductBundle,
) = bundle[1]

function device_model_matrix_execute_detector(
    row::DeviceModelMatrixM5RollingCMOS;
    backend::AdaptiveOpticsSim.Backends.AbstractArrayBackend=CPUBackend(),
    T::Type{<:AbstractFloat}=Float64,
)
    detector = device_model_matrix_detector(row; backend, T)
    map = device_model_matrix_detector_rate_map(row; backend, T)
    definition = RollingShutterAcquisitionDefinition(
        PlantDuration(1_000_000_000);
        readiness_delay=PlantDuration(100_000_000),
    )
    prepared = prepare_rolling_shutter_acquisition(detector, map, definition)
    metadata_before = detector_export_metadata(detector)
    state = RollingShutterAcquisitionState(prepared)
    ready_status = detector_acquisition_status(state)
    start = device_model_matrix_detector_start()
    rng = Xoshiro(0x7_505)

    begin_exposure!(prepared, state, start)
    active_status = detector_acquisition_status(state)
    device_model_matrix_integrate_rolling!(prepared, state, rng)
    pending_status = detector_acquisition_status(state)
    readout = readout_complete_timestamp(state)
    output = complete_readout!(prepared, state, readout, rng)
    complete_status = detector_acquisition_status(state)
    readiness = acquisition_readiness_timestamp(state)
    mark_acquisition_ready!(prepared, state, readiness)
    final_status = detector_acquisition_status(state)

    return (
        row,
        detector,
        map,
        prepared,
        state,
        output,
        rng,
        metadata_before,
        metadata_after=detector_export_metadata(detector),
        status_trace=(
            ready_status,
            active_status,
            pending_status,
            complete_status,
            final_status,
        ),
        event_times=(; start, readout, readiness),
    )
end

function device_model_matrix_execute_detector(
    row::DeviceModelMatrixM6UpTheRampHgCdTe;
    backend::AdaptiveOpticsSim.Backends.AbstractArrayBackend=CPUBackend(),
    T::Type{<:AbstractFloat}=Float64,
)
    detector = device_model_matrix_detector(row; backend, T)
    map = device_model_matrix_detector_rate_map(row; backend, T)
    definition = GlobalShutterAcquisitionDefinition(
        PlantDuration(1_000_000_000);
        readout_duration=PlantDuration(200_000_000),
        readiness_delay=PlantDuration(100_000_000),
    )
    prepared = prepare_global_shutter_acquisition(detector, map, definition)
    metadata_before = detector_export_metadata(detector)
    state = GlobalShutterAcquisitionState(prepared)
    ready_status = detector_acquisition_status(state)
    start = device_model_matrix_detector_start()
    middle = device_model_matrix_detector_ramp_middle(start)
    close = device_model_matrix_detector_close(start)
    rng = Xoshiro(0x7_506)

    begin_exposure!(prepared, state, start)
    active_status = detector_acquisition_status(state)
    take_nondestructive_read!(prepared, state, start, rng)
    accumulate_exposure_interval!(prepared, state, start, middle, rng)
    take_nondestructive_read!(prepared, state, middle, rng)
    map.values .*= T(3)
    accumulate_exposure_interval!(prepared, state, middle, close, rng)
    take_nondestructive_read!(prepared, state, close, rng)
    close_exposure!(prepared, state, close)
    pending_status = detector_acquisition_status(state)
    readout = close + definition.readout_duration
    output = complete_readout!(prepared, state, readout, rng)
    complete_status = detector_acquisition_status(state)
    readiness = readout + definition.readiness_delay
    mark_acquisition_ready!(prepared, state, readiness)
    final_status = detector_acquisition_status(state)

    return (
        row,
        detector,
        map,
        prepared,
        state,
        output,
        rng,
        metadata_before,
        metadata_after=detector_export_metadata(detector),
        status_trace=(
            ready_status,
            active_status,
            pending_status,
            complete_status,
            final_status,
        ),
        event_times=(; start, middle, close, readout, readiness),
    )
end

function device_model_matrix_integrate_rolling!(
    prepared::PreparedRollingShutterAcquisition,
    state::RollingShutterAcquisitionState,
    rng::AbstractRNG,
)
    while true
        next_open = next_rolling_band_open_timestamp(prepared, state)
        next_close = next_rolling_band_close_timestamp(prepared, state)
        next_open === nothing && next_close === nothing && break
        timestamp = if next_open === nothing
            next_close
        elseif next_close === nothing
            next_open
        else
            min(next_open, next_close)
        end
        if integrated_through_timestamp(state) < timestamp &&
                rolling_opened_band_count(state) >
                    rolling_closed_band_count(state)
            accumulate_rolling_exposure_interval!(
                prepared,
                state,
                integrated_through_timestamp(state),
                timestamp,
                rng,
            )
        end
        next_close == timestamp &&
            close_next_rolling_band!(prepared, state, timestamp)
        next_open == timestamp &&
            open_next_rolling_band!(prepared, state, timestamp)
    end
    return nothing
end

function device_model_matrix_repeat_detector!(
    ::Union{DeviceModelMatrixM1CCD,DeviceModelMatrixM3GlobalCMOS,
        DeviceModelMatrixM4GlobalHgCdTe},
    result,
    start::PlantTimestamp,
)
    prepared = result.prepared
    state = result.state
    middle = device_model_matrix_detector_middle(start)
    close = device_model_matrix_detector_close(start)
    begin_exposure!(prepared, state, start)
    accumulate_exposure_interval!(prepared, state, start, middle, result.rng)
    accumulate_exposure_interval!(prepared, state, middle, close, result.rng)
    close_exposure!(prepared, state, close)
    readout = close + prepared.definition.readout_duration
    output = complete_readout!(prepared, state, readout, result.rng)
    mark_acquisition_ready!(
        prepared,
        state,
        readout + prepared.definition.readiness_delay,
    )
    return output
end

function device_model_matrix_repeat_detector!(
    ::DeviceModelMatrixM2FrameTransferEMCCD,
    result,
    start::PlantTimestamp,
)
    prepared = result.prepared
    state = result.state
    middle = device_model_matrix_detector_middle(start)
    close = device_model_matrix_detector_close(start)
    begin_exposure!(prepared, state, start)
    accumulate_exposure_interval!(prepared, state, start, middle, result.rng)
    accumulate_exposure_interval!(prepared, state, middle, close, result.rng)
    close_exposure!(prepared, state, close)
    transfer = frame_transfer_complete_timestamp(state)
    complete_frame_transfer!(prepared, state, transfer)
    return complete_readout!(
        prepared,
        state,
        readout_complete_timestamp(state),
        result.rng,
    )
end

function device_model_matrix_repeat_detector!(
    ::DeviceModelMatrixM5RollingCMOS,
    result,
    start::PlantTimestamp,
)
    prepared = result.prepared
    state = result.state
    begin_exposure!(prepared, state, start)
    device_model_matrix_integrate_rolling!(prepared, state, result.rng)
    output = complete_readout!(
        prepared,
        state,
        readout_complete_timestamp(state),
        result.rng,
    )
    mark_acquisition_ready!(
        prepared,
        state,
        acquisition_readiness_timestamp(state),
    )
    return output
end

function device_model_matrix_repeat_detector!(
    ::DeviceModelMatrixM6UpTheRampHgCdTe,
    result,
    start::PlantTimestamp,
)
    prepared = result.prepared
    state = result.state
    T = eltype(result.detector.products.frame)
    result.map.values .*= inv(T(3))
    middle = device_model_matrix_detector_ramp_middle(start)
    close = device_model_matrix_detector_close(start)
    begin_exposure!(prepared, state, start)
    take_nondestructive_read!(prepared, state, start, result.rng)
    accumulate_exposure_interval!(prepared, state, start, middle, result.rng)
    take_nondestructive_read!(prepared, state, middle, result.rng)
    result.map.values .*= T(3)
    accumulate_exposure_interval!(prepared, state, middle, close, result.rng)
    take_nondestructive_read!(prepared, state, close, result.rng)
    close_exposure!(prepared, state, close)
    readout = close + prepared.definition.readout_duration
    output = complete_readout!(prepared, state, readout, result.rng)
    mark_acquisition_ready!(
        prepared,
        state,
        readout + prepared.definition.readiness_delay,
    )
    return output
end
