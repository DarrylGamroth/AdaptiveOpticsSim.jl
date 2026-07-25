struct ConjugateGeometryPathModel end

struct ConjugateGeometryAcquisitionModel{T<:AbstractFloat}
    exposure_s::T
end

struct ConjugateGeometryOpticModel{T<:AbstractFloat,R}
    surface_value::T
    registration::R
end

struct PreparedConjugateGeometryOptic{
    T<:AbstractFloat,
    A<:AbstractMatrix{T},
    M,
    R,
}
    endpoint::CommandEndpointID
    pattern::A
    metadata::M
    registration::R
end

mutable struct ConjugateGeometryOpticState{T<:AbstractFloat,A}
    visible::T
    surface::A
end

mutable struct ConjugateGeometryOpticWorkspace{T<:AbstractFloat}
    staged::T
end

struct ZeroPupilMaterialization{P}
    destination::P
end

Plant.plant_model_definition_style(
    ::Type{ConjugateGeometryPathModel}) = ColdPlantModelDefinition()
Plant.plant_model_definition_style(
    ::Type{<:ConjugateGeometryAcquisitionModel}) =
    ColdPlantModelDefinition()
Plant.plant_model_definition_style(
    ::Type{<:ConjugateGeometryOpticModel}) =
    ColdPlantModelDefinition()

function Plant.validate_path_materialization_binding(
    materialization::ZeroPupilMaterialization,
    input::PupilFunction,
    ::AbstractAtmosphere,
    ::AbstractSource,
)
    materialization.destination === input || throw(
        PlantPreparationError(:path, :prepared_binding,
            "zero-pupil materialization belongs to another path input"))
    return nothing
end

function Plant.validate_path_materialization(
    materialization::ZeroPupilMaterialization,
    input::PupilFunction,
    ::AbstractAtmosphere,
    ::AtmosphereEpoch,
)
    materialization.destination === input || throw(
        PlantPreparationError(:path, :prepared_binding,
            "zero-pupil materialization belongs to another path input"))
    return input
end

function Plant.materialize_path_input!(
    materialization::ZeroPupilMaterialization,
    input::PupilFunction,
    ::AbstractAtmosphere,
    ::AtmosphereEpoch,
)
    materialization.destination === input || throw(
        PlantPreparationError(:path, :prepared_binding,
            "zero-pupil materialization belongs to another path input"))
    fill!(input.opd, zero(eltype(input.opd)))
    return input
end

function Plant.prepare_path_executor(
    ::ConjugateGeometryPathModel,
    definition::OpticalPathDefinition,
    source::AbstractSource,
    telescope::Telescope,
    atmosphere::AbstractTimedAtmosphere,
)
    T = eltype(pupil_reflectivity(telescope))
    pupil = PupilFunction(telescope; T=T, backend=backend(telescope))
    imaging = prepare_direct_imaging(pupil, source; zero_padding=1)
    return PreparedPathExecutor(
        definition,
        source,
        telescope,
        atmosphere,
        pupil,
        direct_imaging_output(imaging),
        imaging;
        materialization=ZeroPupilMaterialization(pupil),
        optical_model=:conjugate_geometry_test,
        propagation_model=:fraunhofer_fft,
        model_revisions=UInt(1),
    )
end

function Plant.prepare_acquisition_provider(
    model::ConjugateGeometryAcquisitionModel,
    ::AcquisitionDefinition,
    path::PreparedPathExecutor,
)
    require_path_result(path)
    T = eltype(path.result.values)
    detector = Detector(
        integration_time=T(model.exposure_s),
        noise=NoiseNone(),
        qe=one(T),
        gain=one(T),
        response_model=NullFrameResponse(),
        sensor=CMOSSensor(timing_model=GlobalShutter(), T=T),
        T=T,
        backend=path.key.backend,
    )
    execution = FrameAcquisitionExecution(detector, path.result)
    products = AcquisitionProducts(execution.observation;
        metadata=(kind=:conjugate_geometry_test,
            units=:detected_electrons,
            geometry=path.result.metadata))
    return prepare_full_optical_provider(execution, products)
end

function Plant.prepare_controllable_optic(
    model::ConjugateGeometryOpticModel{T},
    definition::ControllableOpticDefinition,
    telescope::Telescope,
    ::AbstractAtmosphere,
) where {T}
    pattern = similar(
        pupil_reflectivity(telescope), T,
        size(pupil_reflectivity(telescope))...)
    fill!(pattern, model.surface_value)
    pupil = PupilFunction(telescope; T=T, backend=backend(telescope))
    metadata = OpticalPlaneMetadata(PupilPlane(), pattern;
        coordinate_domain=MetricCoordinates(),
        sampling=pupil.metadata.sampling,
        origin=pupil.metadata.origin,
        centering=pupil.metadata.centering,
        orientation=pupil.metadata.orientation,
        spectral=AchromaticSpectralCoordinate(),
        normalization=DimensionlessNormalization(),
        spatial_measure=PointSampledMeasure(),
        coherence=NonCombinableProduct())
    return PreparedConjugateGeometryOptic(
        only(command_endpoint_ids(definition)),
        pattern,
        metadata,
        model.registration,
    )
end

function Plant.prepare_controllable_optic_state(
    prepared::PreparedConjugateGeometryOptic{T},
    ::ControllableOpticDefinition,
    endpoint_ids::Tuple,
    initial_commands::Tuple,
) where {T}
    only(endpoint_ids) == prepared.endpoint || throw(
        PlantPreparationError(:controllable_optic, :prepared_binding,
            "conjugate-geometry endpoint binding changed"))
    visible = T(only(initial_commands))
    surface = similar(prepared.pattern)
    @. surface = visible * prepared.pattern
    return ConjugateGeometryOpticState(visible, surface)
end

function Plant.prepare_controllable_optic_workspace(
    ::PreparedConjugateGeometryOptic{T}) where {T}
    return ConjugateGeometryOpticWorkspace(zero(T))
end

function Plant.stage_controllable_optic_command!(
    prepared::PreparedConjugateGeometryOptic{T},
    ::ConjugateGeometryOpticState{T},
    workspace::ConjugateGeometryOpticWorkspace{T},
    endpoint::CommandEndpointID,
    command,
    ::PlantTimestamp,
) where {T}
    endpoint == prepared.endpoint || throw(PlantCommandError(
        :physical_application, :endpoint_mismatch,
        "conjugate-geometry optic received another endpoint"))
    workspace.staged = T(command)
    return nothing
end

function Plant.commit_controllable_optic_command!(
    prepared::PreparedConjugateGeometryOptic{T},
    state::ConjugateGeometryOpticState{T},
    workspace::ConjugateGeometryOpticWorkspace{T},
    ::CommandEndpointID,
    ::PlantTimestamp,
) where {T}
    state.visible = workspace.staged
    @. state.surface = state.visible * prepared.pattern
    return nothing
end

function Plant.prepare_controllable_optic_path_coupling(
    prepared::PreparedConjugateGeometryOptic,
    definition::ControllableOpticDefinition,
    path::PreparedPathExecutor,
)
    return prepare_sampled_pupil_footprint_coupling(
        prepared.metadata,
        prepared.pattern,
        path,
        controllable_optic_placement(definition);
        registration=prepared.registration,
    )
end

function Plant.apply_controllable_optic_surface!(
    input::PupilFunction,
    ::PreparedConjugateGeometryOptic,
    state::ConjugateGeometryOpticState,
    coupling::Union{
        PreparedIdentityPupilFootprintCoupling,
        PreparedPupilFootprintCoupling,
    },
)
    return apply_sampled_pupil_surface!(
        input, state.surface, coupling)
end

function conjugate_geometry_schema(id::Symbol)
    endpoint = Symbol(id, :_command)
    return PlantCommandSchema(
        Float64,
        ();
        id=Symbol(endpoint, :_schema),
        version=1,
        endpoint,
        units=:metre,
        sign_convention=:positive_surface_increases_opd,
        basis=CommandBasis(:sampled_surface_scale, id),
        basis_revision=1,
        semantics=AbsoluteCommand,
        bounds=UniformCommandBounds(-10.0, 10.0),
        value_policy=CommandValuePolicy(),
        sequence_policy=CommandSequencePolicy(),
        effective_time_policy=CommandEffectiveTimePolicy(),
        silence_policy=CommandSilencePolicy(),
    )
end

function conjugate_geometry_atmosphere(telescope::Telescope)
    T = eltype(pupil_reflectivity(telescope))
    return MultiLayerAtmosphere(
        telescope;
        r0=T(0.2),
        L0=T(25),
        fractional_cn2=T[1],
        wind_speed=T[0],
        wind_direction=T[0],
        altitude=T[0],
        layer_ids=(:ground,),
        T=T,
        backend=backend(telescope),
    )
end

function conjugate_geometry_path(
    source::AbstractSource;
    resolution::Int=5,
    diameter::Real=5.0,
    T::Type{<:AbstractFloat}=Float64,
)
    telescope = Telescope(
        resolution=resolution,
        diameter=T(diameter),
        central_obstruction=zero(T),
        T=T,
    )
    atmosphere = conjugate_geometry_atmosphere(telescope)
    definition = PlantDefinition(;
        telescope,
        atmosphere,
        paths=(
            OpticalPathDefinition(
                :geometry,
                source,
                ConjugateGeometryPathModel(),
            ),
        ),
    )
    plant = prepare_plant(definition; run_seed=0x8504)
    return plant, prepared_path(plant, :geometry)
end

function conjugate_surface_metadata(
    prototype::PupilFunction,
    dimensions::NTuple{2,Int};
    sampling=(
        one(eltype(prototype.opd)),
        one(eltype(prototype.opd)),
    ),
    origin=centered_grid_origin(dimensions, sampling),
    orientation::PlaneAxisOrientation=PlaneAxisOrientation(),
)
    T = eltype(prototype.opd)
    sampling_t = (T(sampling[1]), T(sampling[2]))
    origin_t = (T(origin[1]), T(origin[2]))
    values = similar(prototype.opd, T, dimensions)
    metadata = OpticalPlaneMetadata(PupilPlane(), values;
        coordinate_domain=MetricCoordinates(),
        sampling=sampling_t,
        origin=origin_t,
        orientation,
        spectral=AchromaticSpectralCoordinate(),
        normalization=DimensionlessNormalization(),
        spatial_measure=PointSampledMeasure(),
        coherence=NonCombinableProduct(),
        device=prototype.metadata.device)
    return values, metadata
end

function conjugate_logical_coordinates(
    metadata::OpticalPlaneMetadata,
    first::Integer,
    second::Integer,
)
    oriented = (
        metadata.origin[1] + (first - 1) * metadata.sampling[1],
        metadata.origin[2] + (second - 1) * metadata.sampling[2],
    )
    coordinates = (
        metadata.orientation.signs[1] * oriented[1],
        metadata.orientation.signs[2] * oriented[2],
    )
    return metadata.orientation.axes == (:x, :y) ?
        coordinates : (coordinates[2], coordinates[1])
end

function fill_conjugate_linear_surface!(
    surface::AbstractMatrix,
    metadata::OpticalPlaneMetadata;
    x_gain=1.0,
    y_gain=10.0,
)
    @inbounds for second in axes(surface, 2), first in axes(surface, 1)
        x, y = conjugate_logical_coordinates(metadata, first, second)
        surface[first, second] = x_gain * x + y_gain * y
    end
    return surface
end

function conjugate_geometry_coupling(
    path::PreparedPathExecutor,
    surface::AbstractMatrix,
    metadata::OpticalPlaneMetadata,
    placement;
    registration=nothing,
)
    return prepare_sampled_pupil_footprint_coupling(
        metadata,
        surface,
        path,
        placement;
        registration,
    )
end

function apply_conjugate_surface!(
    path::PreparedPathExecutor,
    surface::AbstractMatrix,
    coupling,
)
    pupil = path_input(path)
    fill!(pupil.opd, 0.0)
    return apply_sampled_pupil_surface!(pupil, surface, coupling)
end

@inline function conjugate_surface_application_allocations(
    pupil::PupilFunction,
    surface::AbstractMatrix,
    coupling,
)
    return @allocated apply_sampled_pupil_surface!(
        pupil, surface, coupling)
end

function conjugate_geometry_event_fixture(; split_couplings::Bool=false)
    T = Float64
    telescope = Telescope(
        resolution=5,
        diameter=T(5),
        central_obstruction=zero(T),
        T=T,
    )
    atmosphere = conjugate_geometry_atmosphere(telescope)
    source = LGSSource(
        altitude=T(10_000),
        photon_irradiance=T(100),
        T=T,
    )
    path = OpticalPathDefinition(
        :science,
        source,
        ConjugateGeometryPathModel(),
    )
    acquisition = AcquisitionDefinition(
        :camera,
        :science,
        ConjugateGeometryAcquisitionModel(T(0.02)),
    )
    placement = AtmosphericConjugatePlacement(T(5_000))
    registration = PupilRelayRegistration(T=T)
    second_registration = split_couplings ?
        PupilRelayRegistration(decenter_m=(T(0.25), zero(T)), T=T) :
        registration
    first = ControllableOpticDefinition(
        :first,
        ConjugateGeometryOpticModel(T(1), registration),
        (conjugate_geometry_schema(:first),);
        placement,
        visibility=AllPathVisibility(),
    )
    second = ControllableOpticDefinition(
        :second,
        ConjugateGeometryOpticModel(T(1), second_registration),
        (conjugate_geometry_schema(:second),);
        placement,
        visibility=AllPathVisibility(),
    )
    plant = prepare_plant(
        PlantDefinition(;
            telescope,
            atmosphere,
            controllable_optics=(second, first),
            paths=(path,),
            acquisitions=(acquisition,),
        );
        run_seed=0x8505,
        command_endpoints=(
            CommandEndpointConfiguration(
                :second_command, T(2); capacity=2),
            CommandEndpointConfiguration(
                :first_command, T(1); capacity=2),
        ),
    )
    event_loop = prepare_plant_event_loop(
        plant,
        PlantEventLoopDefinition(
            (
                OpticalSampleDefinition(
                    :science,
                    PeriodicSchedule(period_ns=100_000_000),
                ),
            ),
            (
                AcquisitionEventDefinition(
                    :camera,
                    GlobalShutterAcquisitionDefinition(
                        PlantDuration(20_000_000)),
                    PeriodicAcquisitionStart(
                        PeriodicSchedule(period_ns=1_000_000_000)),
                ),
            ),
        ),
    )
    return plant, event_loop, PlantEventLoopState(event_loop),
        PlantEventLoopWorkspace(event_loop)
end

@testset "Pupil-relay registration contract" begin
    registration = @inferred PupilRelayRegistration(
        magnification=(1.5, 0.75),
        rotation_deg=30.0,
        parity=(-1, 1),
        decenter_m=(0.25, -0.5),
    )
    @test pupil_relay_magnification(registration) == (1.5, 0.75)
    @test pupil_relay_rotation_deg(registration) ≈ 30.0
    @test pupil_relay_rotation_rad(registration) ≈ pi / 6
    @test pupil_relay_parity(registration) == (Int8(-1), Int8(1))
    @test pupil_relay_decenter_m(registration) == (0.25, -0.5)
    @test !Base.ismutable(registration)

    for operation in (
        () -> PupilRelayRegistration(magnification=0.0),
        () -> PupilRelayRegistration(magnification=(1.0, Inf)),
        () -> PupilRelayRegistration(rotation_deg=NaN),
        () -> PupilRelayRegistration(parity=(1, 0)),
        () -> PupilRelayRegistration(parity=(1, 257)),
        () -> PupilRelayRegistration(parity=[1, 1]),
        () -> PupilRelayRegistration(magnification=[1.0, 1.0]),
        () -> PupilRelayRegistration(decenter_m=(0.0, Inf)),
        () -> PupilRelayRegistration(T=AbstractFloat),
    )
        error = try
            operation()
            nothing
        catch caught
            caught
        end
        @test error isa PlantDefinitionError
        @test error.component === :controllable_optic
        @test error.reason === :invalid_pupil_relay_registration
    end

    source32 = Source(
        band=:custom,
        wavelength=Float32(0.8e-6),
        photon_irradiance=Float32(100),
        T=Float32,
    )
    _, path32 = conjugate_geometry_path(source32; T=Float32)
    pupil32 = path_input(path32)
    surface32, metadata32 =
        conjugate_surface_metadata(pupil32, (5, 5))
    fill!(surface32, one(Float32))
    converted = conjugate_geometry_coupling(
        path32,
        surface32,
        metadata32,
        PupilPlanePlacement();
        registration=PupilRelayRegistration(
            magnification=(1.25, 0.75),
            rotation_deg=15.0,
            parity=(-1, 1),
            decenter_m=(0.125, -0.25),
        ),
    )
    @test converted isa PreparedPupilFootprintCoupling
    @test all(value -> value isa Float32, converted.index_transform)
    @test all(value -> value isa Float32, converted.index_offset)

    invalid_registration_error = try
        conjugate_geometry_coupling(
            path32,
            surface32,
            metadata32,
            PupilPlanePlacement();
            registration=:invalid,
        )
        nothing
    catch caught
        caught
    end
    @test invalid_registration_error isa PlantPreparationError
    @test invalid_registration_error.reason ===
        :invalid_pupil_relay_registration

    for overflowing_registration in (
        PupilRelayRegistration(magnification=(1e40, 1.0)),
        PupilRelayRegistration(rotation_deg=1e41),
        PupilRelayRegistration(decenter_m=(1e40, 0.0)),
    )
        error = try
            conjugate_geometry_coupling(
                path32,
                surface32,
                metadata32,
                PupilPlanePlacement();
                registration=overflowing_registration,
            )
            nothing
        catch caught
            caught
        end
        @test error isa PlantPreparationError
        @test error.reason === :invalid_pupil_relay_registration
    end
end

@testset "Analytic NGS and finite-height LGS footprints" begin
    arcsec_per_radian = 180 * 3600 / pi
    cases = (
        (
            source=Source(
                band=:custom,
                wavelength=0.8e-6,
                photon_irradiance=100.0,
            ),
            altitude=2_000.0,
            shift=(0.0, 0.0),
            scale=1.0,
            identity=true,
        ),
        (
            source=Source(
                band=:custom,
                wavelength=0.8e-6,
                coordinates=(
                    sqrt(2) * arcsec_per_radian / 1_000,
                    -45.0,
                ),
                photon_irradiance=100.0,
            ),
            altitude=1_000.0,
            shift=(1.0, -1.0),
            scale=1.0,
            identity=false,
        ),
        (
            source=with_spectrum(
                Source(
                    band=:custom,
                    wavelength=0.8e-6,
                    coordinates=(
                        arcsec_per_radian / 2_000,
                        90.0,
                    ),
                    photon_irradiance=100.0,
                ),
                SpectralBundle(
                    [0.7e-6, 0.9e-6],
                    [0.4, 0.6],
                ),
            ),
            altitude=2_000.0,
            shift=(0.0, 1.0),
            scale=1.0,
            identity=false,
        ),
        (
            source=LGSSource(
                altitude=10_000.0,
                photon_irradiance=100.0,
            ),
            altitude=2_500.0,
            shift=(0.0, 0.0),
            scale=0.75,
            identity=false,
        ),
        (
            source=LGSSource(
                coordinates=(
                    sqrt(0.5) * arcsec_per_radian / 5_000,
                    45.0,
                ),
                altitude=10_000.0,
                photon_irradiance=100.0,
            ),
            altitude=5_000.0,
            shift=(0.5, 0.5),
            scale=0.5,
            identity=false,
        ),
    )

    for case in cases
        _, path = conjugate_geometry_path(case.source)
        pupil = path_input(path)
        surface, metadata =
            conjugate_surface_metadata(pupil, (5, 5))
        fill_conjugate_linear_surface!(surface, metadata)
        coupling = conjugate_geometry_coupling(
            path,
            surface,
            metadata,
            AtmosphericConjugatePlacement(case.altitude),
        )
        @test (coupling isa
            PreparedIdentityPupilFootprintCoupling) == case.identity
        apply_conjugate_surface!(path, surface, coupling)
        expected = similar(pupil.opd)
        @inbounds for second in axes(expected, 2),
            first in axes(expected, 1)
            x, y = conjugate_logical_coordinates(
                pupil.metadata, first, second)
            surface_x = case.scale * x + case.shift[1]
            surface_y = case.scale * y + case.shift[2]
            expected[first, second] =
                -2.0 <= surface_x <= 2.0 &&
                -2.0 <= surface_y <= 2.0 ?
                surface_x + 10surface_y : 0.0
        end
        @test pupil.opd ≈ expected rtol=0.0 atol=64eps(Float64)
    end
end

@testset "Relay registration and sampled-grid orientation" begin
    source = Source(
        band=:custom,
        wavelength=0.8e-6,
        photon_irradiance=100.0,
    )
    _, path = conjugate_geometry_path(source)
    pupil = path_input(path)
    orientation = PlaneAxisOrientation((:y, :x), (-1, 1))
    surface, metadata = conjugate_surface_metadata(
        pupil,
        (9, 9);
        sampling=(1.0, 1.0),
        origin=(-4.0, -4.0),
        orientation,
    )
    fill_conjugate_linear_surface!(
        surface, metadata; x_gain=2.0, y_gain=3.0)
    registration = PupilRelayRegistration(
        magnification=(1.0, 1.0),
        rotation_deg=90.0,
        parity=(-1, 1),
        decenter_m=(0.25, -0.5),
    )
    coupling = conjugate_geometry_coupling(
        path,
        surface,
        metadata,
        PupilPlanePlacement();
        registration,
    )
    @test coupling isa PreparedPupilFootprintCoupling
    apply_conjugate_surface!(path, surface, coupling)
    m11, m12, m21, m22 = registration.transform
    expected = similar(pupil.opd)
    @inbounds for second in axes(expected, 2),
        first in axes(expected, 1)
        x, y = conjugate_logical_coordinates(
            pupil.metadata, first, second)
        surface_x = muladd(
            m12, y, muladd(m11, x, registration.decenter_m[1]))
        surface_y = muladd(
            m22, y, muladd(m21, x, registration.decenter_m[2]))
        expected[first, second] = 2surface_x + 3surface_y
    end
    @test pupil.opd ≈ expected rtol=0.0 atol=128eps(Float64)
end

@testset "Finite support, impulse interpolation, and errors" begin
    source = Source(
        band=:custom,
        wavelength=0.8e-6,
        photon_irradiance=100.0,
    )
    _, path = conjugate_geometry_path(source)
    pupil = path_input(path)
    surface, metadata = conjugate_surface_metadata(pupil, (5, 5))
    fill!(surface, 0.0)
    surface[3, 3] = 1.0
    half_sample = PupilRelayRegistration(
        decenter_m=(0.5, 0.5))
    coupling = conjugate_geometry_coupling(
        path,
        surface,
        metadata,
        PupilPlanePlacement();
        registration=half_sample,
    )
    apply_conjugate_surface!(path, surface, coupling)
    expected = zeros(5, 5)
    expected[2:3, 2:3] .= 0.25
    @test pupil.opd == expected

    outside = conjugate_geometry_coupling(
        path,
        surface,
        metadata,
        PupilPlanePlacement();
        registration=PupilRelayRegistration(decenter_m=(100.0, 0.0)),
    )
    apply_conjugate_surface!(path, surface, outside)
    @test all(iszero, pupil.opd)

    low_lgs = LGSSource(
        altitude=5_000.0,
        photon_irradiance=100.0,
    )
    _, low_path = conjugate_geometry_path(low_lgs)
    low_pupil = path_input(low_path)
    low_surface, low_metadata =
        conjugate_surface_metadata(low_pupil, (5, 5))
    error = try
        conjugate_geometry_coupling(
            low_path,
            low_surface,
            low_metadata,
            AtmosphericConjugatePlacement(5_000.0),
        )
        nothing
    catch caught
        caught
    end
    @test error isa PlantPreparationError
    @test error.reason === :source_below_conjugate

    invalid_values = zeros(5, 5)
    invalid_metadata = OpticalPlaneMetadata(
        PupilPlane(),
        invalid_values;
        coordinate_domain=MetricCoordinates(),
        sampling=(1.0, 1.0),
        orientation=PlaneAxisOrientation((:row, :column)),
    )
    error = try
        conjugate_geometry_coupling(
            path,
            invalid_values,
            invalid_metadata,
            PupilPlanePlacement(),
        )
        nothing
    catch caught
        caught
    end
    @test error isa PlantPreparationError
    @test error.reason === :unsupported_axis_orientation

    directional_source = Asterism([
        source,
        Source(
            band=:custom,
            wavelength=0.8e-6,
            coordinates=(1.0, 0.0),
            photon_irradiance=100.0,
        ),
    ])
    _, directional_path = conjugate_geometry_path(directional_source)
    directional_pupil = path_input(directional_path)
    directional_surface, directional_metadata =
        conjugate_surface_metadata(directional_pupil, (5, 5))
    error = try
        conjugate_geometry_coupling(
            directional_path,
            directional_surface,
            directional_metadata,
            AtmosphericConjugatePlacement(1_000.0),
        )
        nothing
    catch caught
        caught
    end
    @test error isa PlantPreparationError
    @test error.reason === :unsupported_source_geometry
end

@testset "Inferred zero-allocation sampled-surface application" begin
    source = LGSSource(
        altitude=10_000.0,
        photon_irradiance=100.0,
    )
    _, path = conjugate_geometry_path(source)
    pupil = path_input(path)
    surface, metadata = conjugate_surface_metadata(pupil, (5, 5))
    fill_conjugate_linear_surface!(surface, metadata)
    coupling = conjugate_geometry_coupling(
        path,
        surface,
        metadata,
        AtmosphericConjugatePlacement(5_000.0),
    )
    @test @inferred(apply_sampled_pupil_surface!(
        pupil, surface, coupling)) === pupil
    fill!(pupil.opd, 0.0)
    conjugate_surface_application_allocations(
        pupil, surface, coupling)
    fill!(pupil.opd, 0.0)
    @test conjugate_surface_application_allocations(
        pupil, surface, coupling) == 0
end

@testset "Event composition binds compatible co-conjugated geometry" begin
    plant, event_loop, state, workspace =
        conjugate_geometry_event_fixture()
    event_path = only(event_loop.paths)
    @test length(event_path.optic_couplings) == 2
    @test length(event_path.optic_coupling_groups) == 1
    @test all(coupling ->
            coupling isa PreparedPupilFootprintCoupling,
        event_path.optic_couplings)

    split_plant, split_event_loop, split_state, split_workspace =
        conjugate_geometry_event_fixture(split_couplings=true)
    split_event_path = only(split_event_loop.paths)
    @test length(split_event_path.optic_coupling_groups) == 2
    @test map(
        group -> group.binding_count,
        split_event_path.optic_coupling_groups,
    ) == [1, 1]
    @test step_plant_events!(
        split_event_loop, split_state, split_workspace) ==
        PlantTimestamp(0)
    @test all(==(3.0), path_input(
        prepared_path(split_plant, :science)).opd)

    @test step_plant_events!(event_loop, state, workspace) ==
        PlantTimestamp(0)
    pupil = path_input(prepared_path(plant, :science))
    @test all(==(3.0), pupil.opd)

    schema = command_schema(
        first(controllable_optic_definitions(
            getfield(plant, :definition))),
        :second_command,
    )
    admission = admit_plant_command!(
        event_loop,
        state,
        workspace,
        PlantCommand(
            schema,
            1,
            PlantTimestamp(100_000_000),
            4.0,
        ),
        PlantTimestamp(1),
    )
    @test command_admission_status(admission) ==
        CommandAdmittedPending
    @test step_plant_events!(event_loop, state, workspace) ==
        PlantTimestamp(20_000_000)
    @test step_plant_events!(event_loop, state, workspace) ==
        PlantTimestamp(100_000_000)
    @test all(==(5.0), pupil.opd)
end
