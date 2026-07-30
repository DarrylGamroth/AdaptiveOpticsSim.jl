module Gate5OpticalPlacementBenchmark

using AdaptiveOpticsSim
using AdaptiveOpticsSim.Atmospheres
using AdaptiveOpticsSim.Detectors
using AdaptiveOpticsSim.Optics
using AdaptiveOpticsSim.Plant
using SHA

const AOS = AdaptiveOpticsSim
const AOSPlant = AdaptiveOpticsSim.Plant

struct Gate5PathModel end

struct Gate5FrameAcquisitionModel{T<:AbstractFloat}
    exposure_s::T
end

struct ZeroPupilMaterialization{P}
    destination::P
end

AOSPlant.plant_model_definition_style(::Type{Gate5PathModel}) =
    AOSPlant.ColdPlantModelDefinition()
AOSPlant.plant_model_definition_style(
    ::Type{<:Gate5FrameAcquisitionModel},
) = AOSPlant.ColdPlantModelDefinition()

function AOSPlant.validate_path_materialization_binding(
    materialization::ZeroPupilMaterialization,
    input::AOS.Optics.PupilFunction,
    ::AOS.Atmospheres.AbstractAtmosphere,
    ::AOS.Optics.AbstractSource,
)
    materialization.destination === input || throw(
        AOSPlant.PlantPreparationError(
            :path,
            :prepared_binding,
            "Gate 5 zero-pupil materialization belongs to another path",
        ),
    )
    return nothing
end

function AOSPlant.validate_path_materialization(
    materialization::ZeroPupilMaterialization,
    input::AOS.Optics.PupilFunction,
    ::AOS.Atmospheres.AbstractAtmosphere,
    ::AOS.Atmospheres.AtmosphereEpoch,
)
    materialization.destination === input || throw(
        AOSPlant.PlantPreparationError(
            :path,
            :prepared_binding,
            "Gate 5 zero-pupil materialization belongs to another path",
        ),
    )
    return input
end

function AOSPlant.materialize_path_input!(
    materialization::ZeroPupilMaterialization,
    input::AOS.Optics.PupilFunction,
    ::AOS.Atmospheres.AbstractAtmosphere,
    ::AOS.Atmospheres.AtmosphereEpoch,
)
    materialization.destination === input || throw(
        AOSPlant.PlantPreparationError(
            :path,
            :prepared_binding,
            "Gate 5 zero-pupil materialization belongs to another path",
        ),
    )
    fill!(input.opd, zero(eltype(input.opd)))
    return input
end

function AOSPlant.prepare_path_executor(
    ::Gate5PathModel,
    definition::AOSPlant.OpticalPathDefinition,
    source::AOS.Optics.AbstractSource,
    telescope::AOS.Optics.Telescope,
    atmosphere::AOS.Atmospheres.AbstractTimedAtmosphere,
)
    T = eltype(AOS.Optics.pupil_reflectivity(telescope))
    pupil = AOS.Optics.PupilFunction(
        telescope;
        T,
        backend=AOS.backend(telescope),
    )
    imaging = AOS.Optics.prepare_direct_imaging(pupil, source; zero_padding=1)
    return AOSPlant.PreparedPathExecutor(
        definition,
        source,
        telescope,
        atmosphere,
        pupil,
        AOS.Optics.direct_imaging_output(imaging),
        imaging;
        materialization=ZeroPupilMaterialization(pupil),
        optical_model=:gate5_conjugated_optical_placement,
        propagation_model=:fraunhofer_fft,
        model_revisions=UInt(1),
    )
end

function AOSPlant.prepare_acquisition_provider(
    model::Gate5FrameAcquisitionModel,
    ::AOSPlant.AcquisitionDefinition,
    path::AOSPlant.PreparedPathExecutor,
)
    AOSPlant.require_path_result(path)
    T = eltype(path.result.values)
    detector = AOS.Detector(
        integration_time=T(model.exposure_s),
        noise=AOS.NoiseNone(),
        qe=one(T),
        gain=one(T),
        response_model=AOS.NullFrameResponse(),
        sensor=AOS.CMOSSensor(timing_model=AOS.GlobalShutter(), T=T),
        T=T,
        backend=path.key.backend,
    )
    execution = AOSPlant.FrameAcquisitionExecution(detector, path.result)
    products = AOSPlant.AcquisitionProducts(
        execution.observation;
        metadata=(
            kind=:gate5_path_frame,
            units=:detected_electrons,
            geometry=path.result.metadata,
            semantics=:complete_acquisition,
        ),
    )
    return AOSPlant.prepare_full_optical_provider(execution, products)
end

function command_schema(endpoint::Symbol)
    return AOSPlant.PlantCommandSchema(
        Float64,
        (1,);
        id=Symbol(endpoint, :_schema),
        version=1,
        endpoint,
        units=:metre,
        sign_convention=:positive_surface_increases_opd,
        basis=AOSPlant.CommandBasis(:actuator, endpoint),
        basis_revision=1,
        semantics=AOSPlant.AbsoluteCommand,
        bounds=AOSPlant.UniformCommandBounds(-1.0e-3, 1.0e-3),
        value_policy=AOSPlant.CommandValuePolicy(),
        sequence_policy=AOSPlant.CommandSequencePolicy(),
        effective_time_policy=AOSPlant.CommandEffectiveTimePolicy(),
        silence_policy=AOSPlant.CommandSilencePolicy(),
    )
end

function sampled_surface_metadata(
    prototype::AOS.Optics.PupilFunction,
    surface::AbstractMatrix,
)
    return AOS.Optics.OpticalPlaneMetadata(
        AOS.Optics.PupilPlane(),
        surface;
        coordinate_domain=AOS.Optics.MetricCoordinates(),
        sampling=prototype.metadata.sampling,
        orientation=prototype.metadata.orientation,
        spectral=AOS.Optics.AchromaticSpectralCoordinate(),
        normalization=AOS.Optics.DimensionlessNormalization(),
        spatial_measure=AOS.Optics.PointSampledMeasure(),
        coherence=AOS.Optics.NonCombinableProduct(),
    )
end

function gate5_path_id(index::Integer, wfs_count::Integer)
    local_index = index <= wfs_count ? index : index - wfs_count
    prefix = index <= wfs_count ? "wfs" : "science"
    return Symbol(prefix, "_", lpad(string(local_index), 2, '0'))
end

function gate5_source(
    index::Integer,
    wfs_count::Integer,
    path_count::Integer,
    raw::AbstractDict,
)
    radius_arcsec = Float64(raw["source_radius_arcsec"]) *
        (0.5 + 0.5 * mod(index - 1, 3))
    azimuth_deg = 360.0 * (index - 1) / path_count
    coordinates = (radius_arcsec, azimuth_deg)
    if index <= wfs_count && isodd(index)
        return AOS.Optics.LGSSource(
            coordinates=coordinates,
            altitude=Float64(raw["lgs_altitude_m"]),
            photon_irradiance=Float64(raw["photon_irradiance"]),
        )
    end
    wavelength = index <= wfs_count ?
        Float64(raw["wfs_wavelength_m"]) :
        Float64(raw["science_wavelength_m"])
    return AOS.Optics.Source(
        band=:custom,
        wavelength=wavelength,
        coordinates=coordinates,
        photon_irradiance=Float64(raw["photon_irradiance"]),
    )
end

function deformable_mirror_definition(
    id::Symbol,
    placement,
    visibility,
    resolution::Integer,
    command_value::Real,
    registration::AOSPlant.PupilRelayRegistration,
)
    endpoint = Symbol(id, :_command)
    topology = AOS.SampledActuatorTopology(zeros(Float64, 2, 1))
    influence = AOS.DenseInfluenceMatrix(
        fill(1.0, Int(resolution)^2, 1),
    )
    model = AOSPlant.DeformableMirrorModel(
        topology=topology,
        influence_model=influence,
        pupil_relay_registration=registration,
    )
    definition = AOSPlant.ControllableOpticDefinition(
        id,
        model,
        (command_schema(endpoint),);
        placement,
        visibility,
    )
    configuration = AOSPlant.CommandEndpointConfiguration(
        endpoint,
        [Float64(command_value)];
        capacity=4,
    )
    return definition, configuration
end

function _gate5_common_optics(
    resolution::Integer,
    raw::AbstractDict,
    registration::AOSPlant.PupilRelayRegistration,
)
    values = Float64.(raw["common_dm_commands_m"])
    altitudes = Float64.(raw["common_dm_altitudes_m"])
    length(values) == 3 ||
        error("Gate 5 workload requires three common DM command values")
    length(altitudes) == 3 ||
        error("Gate 5 workload requires three common DM altitudes")
    iszero(first(altitudes)) ||
        error("Gate 5 first common DM altitude must be the pupil at 0 m")
    placements = (
        AOSPlant.PupilPlanePlacement(),
        AOSPlant.AtmosphericConjugatePlacement(altitudes[2]),
        AOSPlant.AtmosphericConjugatePlacement(altitudes[3]),
    )
    ids = (:common_ground, :common_mid, :common_high)
    return map(ids, placements, values) do id, placement, value
        deformable_mirror_definition(
            id,
            placement,
            AOSPlant.AllPathVisibility(),
            resolution,
            value,
            registration,
        )
    end
end

function gate5_plant_definition(
    raw::AbstractDict;
    path_count::Integer=Int(raw["path_count"]),
    resolution::Integer=Int(raw["resolution"]),
    reverse_declarations::Bool=false,
)
    path_count >= 4 && iseven(path_count) ||
        error("Gate 5 path count must be an even integer of at least four")
    resolution >= 5 ||
        error("Gate 5 resolution must be at least five")
    wfs_count = path_count ÷ 2
    science_count = path_count - wfs_count
    telescope = AOS.Optics.Telescope(
        resolution=resolution,
        diameter=Float64(raw["diameter_m"]),
        central_obstruction=Float64(raw["central_obstruction"]),
    )
    atmosphere = AOS.Atmospheres.MultiLayerAtmosphere(
        telescope;
        r0=Float64(raw["r0_m"]),
        L0=Float64(raw["outer_scale_m"]),
        fractional_cn2=[1.0],
        wind_speed=[0.0],
        wind_direction=[0.0],
        altitude=[0.0],
        layer_ids=(:ground,),
    )
    path_ids = ntuple(
        index -> gate5_path_id(index, wfs_count),
        path_count,
    )
    sources = ntuple(
        index -> gate5_source(index, wfs_count, path_count, raw),
        path_count,
    )
    paths = map(path_ids, sources) do id, source
        AOSPlant.OpticalPathDefinition(id, source, Gate5PathModel())
    end
    science_ids = path_ids[(wfs_count + 1):end]
    registration = AOSPlant.PupilRelayRegistration(
        magnification=Float64(raw["pupil_relay_magnification"]),
    )

    common_pairs = _gate5_common_optics(resolution, raw, registration)
    common_definitions = map(first, common_pairs)
    common_configurations = map(last, common_pairs)
    moao_base = Float64(raw["moao_command_base_m"])
    moao_step = Float64(raw["moao_command_step_m"])
    moao_pairs = map(enumerate(science_ids)) do (index, path_id)
        value = moao_base + index * moao_step
        deformable_mirror_definition(
            Symbol(path_id, :_moao),
            AOSPlant.AtmosphericConjugatePlacement(
                Float64(raw["moao_altitude_m"]),
            ),
            AOSPlant.SelectedPathVisibility(path_id),
            resolution,
            value,
            registration,
        )
    end
    moao_definitions = map(first, moao_pairs)
    moao_configurations = map(last, moao_pairs)
    optics = (common_definitions..., moao_definitions...)
    configurations =
        (common_configurations..., moao_configurations...)

    prototype = AOS.Optics.PupilFunction(telescope)
    common_static_value = Float64(raw["common_static_opd_m"])
    common_static_opd = fill(common_static_value, size(prototype.opd))
    common_static = AOSPlant.SampledAberrationDefinition(
        :common_static,
        AOS.Optics.OPDMap(common_static_opd),
        sampled_surface_metadata(prototype, common_static_opd);
        placement=AOSPlant.PupilPlanePlacement(),
        visibility=AOSPlant.AllPathVisibility(),
        application=AOS.Optics.DMReplace(),
    )
    ncpa_base = Float64(raw["ncpa_base_m"])
    ncpa_step = Float64(raw["ncpa_step_m"])
    ncpa_definitions = map(enumerate(science_ids)) do (index, path_id)
        value = ncpa_base + index * ncpa_step
        opd = fill(value, size(prototype.opd))
        AOSPlant.SampledAberrationDefinition(
            Symbol(path_id, :_ncpa),
            AOS.Optics.NCPA(opd),
            sampled_surface_metadata(prototype, opd);
            placement=AOSPlant.PupilPlanePlacement(),
            visibility=AOSPlant.SelectedPathVisibility(path_id),
            application=AOS.Optics.DMAdditive(),
            registration,
        )
    end
    aberrations = (common_static, ncpa_definitions...)

    common_sum = common_static_value + sum(Float64.(
        raw["common_dm_commands_m"],
    ))
    expected = map(enumerate(path_ids)) do (index, path_id)
        if index <= wfs_count
            (; id=path_id, opd=common_sum, role=:wfs)
        else
            science_index = index - wfs_count
            moao = moao_base + science_index * moao_step
            ncpa = ncpa_base + science_index * ncpa_step
            (; id=path_id, opd=common_sum + moao + ncpa,
                role=:science)
        end
    end
    acquisitions = map(path_ids) do path_id
        AOSPlant.AcquisitionDefinition(
            Symbol(path_id, :_camera),
            path_id,
            Gate5FrameAcquisitionModel(
                Float64(raw["acquisition_duration_ns"]) / 1.0e9,
            ),
        )
    end

    if reverse_declarations
        optics = reverse(optics)
        configurations = reverse(configurations)
        aberrations = reverse(aberrations)
        paths = reverse(paths)
    end
    definition = AOSPlant.PlantDefinition(
        telescope=telescope,
        atmosphere=atmosphere,
        controllable_optics=optics,
        sampled_aberrations=aberrations,
        paths=paths,
        acquisitions=acquisitions,
    )
    return definition, configurations, expected, science_count
end

function gate5_event_definition(
    paths::AbstractVector,
    raw::AbstractDict;
    reverse_declarations::Bool=false,
)
    schedule = AOSPlant.PeriodicSchedule(
        period_ns=Int64(raw["optical_sample_period_ns"]),
    )
    samples = map(paths) do path
        AOSPlant.OpticalSampleDefinition(
            AOSPlant.path_id(path),
            schedule,
        )
    end
    reverse_declarations && (samples = reverse(samples))
    acquisition_period_ns = Int64(raw["acquisition_period_ns"])
    acquisition_duration_ns = Int64(raw["acquisition_duration_ns"])
    acquisitions = map(paths) do path
        AOSPlant.AcquisitionEventDefinition(
            Symbol(AOSPlant.path_id(path).name, :_camera),
            AOSPlant.GlobalShutterAcquisitionDefinition(
                AOSPlant.PlantDuration(acquisition_duration_ns),
            ),
            AOSPlant.PeriodicAcquisitionStart(
                AOSPlant.PeriodicSchedule(
                    period_ns=acquisition_period_ns,
                ),
            ),
        )
    end
    return AOSPlant.PlantEventLoopDefinition(samples, acquisitions)
end

mutable struct Gate5OpticalOperation{P,E,S,W}
    plant::P
    event_loop::E
    state::S
    workspace::W
    cycles::UInt64
end

function prepare_gate5_operation(
    raw::AbstractDict;
    path_count::Integer=Int(raw["path_count"]),
    resolution::Integer=Int(raw["resolution"]),
    reverse_declarations::Bool=false,
)
    definition, configurations, expected, science_count =
        gate5_plant_definition(
            raw;
            path_count,
            resolution,
            reverse_declarations,
        )
    plant = AOSPlant.prepare_plant(
        definition;
        run_seed=UInt64(raw["run_seed"]),
        command_endpoints=configurations,
    )
    event_definition = gate5_event_definition(
        AOSPlant.path_definitions(definition),
        raw;
        reverse_declarations,
    )
    event_loop = AOSPlant.prepare_plant_event_loop(
        plant,
        event_definition,
    )
    operation = Gate5OpticalOperation(
        plant,
        event_loop,
        AOSPlant.PlantEventLoopState(event_loop),
        AOSPlant.PlantEventLoopWorkspace(event_loop),
        UInt64(0),
    )
    return operation, expected, science_count
end

@inline function (operation::Gate5OpticalOperation)()
    timestamp = AOSPlant.step_plant_events!(
        operation.event_loop,
        operation.state,
        operation.workspace,
    )
    operation.cycles += UInt64(1)
    return timestamp
end

function run_gate5_cycles!(
    operation::Gate5OpticalOperation,
    cycles::Integer,
)
    cycles >= 0 || error("Gate 5 cycle count must be nonnegative")
    timestamp = AOSPlant.scheduler_timestamp(operation.state.scheduler)
    @inbounds for _ in 1:cycles
        timestamp = operation()
    end
    return timestamp
end

function path_opd(operation::Gate5OpticalOperation, id::Symbol)
    path = AOSPlant.prepared_path(operation.plant, id)
    return AOSPlant.path_input(path).opd
end

function path_output(operation::Gate5OpticalOperation, id::Symbol)
    path = AOSPlant.prepared_path(operation.plant, id)
    return AOSPlant.path_result(path).values
end

function array_sha256(array::AbstractArray)
    return bytes2hex(SHA.sha256(reinterpret(UInt8, vec(array))))
end

function path_hashes(
    operation::Gate5OpticalOperation,
    expected,
)
    return Dict(
        String(item.id) => Dict(
            "opd_sha256" => array_sha256(path_opd(operation, item.id)),
            "result_sha256" =>
                array_sha256(path_output(operation, item.id)),
        ) for item in expected
    )
end

function numerical_oracle(
    operation::Gate5OpticalOperation,
    expected,
)
    maximum_error = 0.0
    wfs_maximum_error = 0.0
    science_maximum_error = 0.0
    for item in expected
        error = maximum(abs, path_opd(operation, item.id) .- item.opd)
        maximum_error = max(maximum_error, error)
        if item.role === :wfs
            wfs_maximum_error = max(wfs_maximum_error, error)
        else
            science_maximum_error = max(science_maximum_error, error)
        end
    end
    return (
        maximum_error=maximum_error,
        wfs_maximum_error=wfs_maximum_error,
        science_maximum_error=science_maximum_error,
    )
end

function validate_integrated_oracle(raw::AbstractDict)
    canonical, expected, science_count = prepare_gate5_operation(raw)
    reordered, reordered_expected, reordered_science_count =
        prepare_gate5_operation(raw; reverse_declarations=true)
    science_count == reordered_science_count ||
        error("Gate 5 declaration reordering changed science-path count")
    canonical_timestamp = canonical()
    reordered_timestamp = reordered()
    canonical_timestamp == reordered_timestamp ||
        error("Gate 5 declaration reordering changed event time")
    canonical_error = numerical_oracle(canonical, expected)
    reordered_error = numerical_oracle(reordered, reordered_expected)
    order_error = 0.0
    output_order_error = 0.0
    for item in expected
        order_error = max(order_error, maximum(
            abs,
            path_opd(canonical, item.id) .-
                path_opd(reordered, item.id),
        ))
        output_order_error = max(output_order_error, maximum(
            abs,
            path_output(canonical, item.id) .-
                path_output(reordered, item.id),
        ))
    end
    opd_tolerance = Float64(raw["numerical_atol"])
    output_order_tolerance = Float64(raw["output_order_atol"])
    maximum_opd_error = max(
        canonical_error.maximum_error,
        reordered_error.maximum_error,
        order_error,
    )
    maximum_opd_error <= opd_tolerance ||
        error("Gate 5 integrated numerical oracle exceeded tolerance")
    output_order_error <= output_order_tolerance ||
        error("Gate 5 propagated declaration-order oracle exceeded tolerance")
    return Dict{String,Any}(
        "path_count" => length(expected),
        "science_path_count" => science_count,
        "wfs_path_count" => length(expected) - science_count,
        "maximum_opd_error_m" => max(
            canonical_error.maximum_error,
            reordered_error.maximum_error,
        ),
        "wfs_maximum_opd_error_m" =>
            canonical_error.wfs_maximum_error,
        "science_maximum_opd_error_m" =>
            canonical_error.science_maximum_error,
        "declaration_order_opd_error_m" => order_error,
        "declaration_order_output_error" => output_order_error,
        "tolerance" => opd_tolerance,
        "output_order_tolerance" => output_order_tolerance,
        "passed" => true,
        "path_hashes" => path_hashes(canonical, expected),
    )
end

function validate_finite_support(raw::AbstractDict)
    telescope = AOS.Optics.Telescope(
        resolution=5,
        diameter=5.0,
        central_obstruction=0.0,
    )
    atmosphere = AOS.Atmospheres.MultiLayerAtmosphere(
        telescope;
        r0=0.2,
        L0=25.0,
        fractional_cn2=[1.0],
        wind_speed=[0.0],
        wind_direction=[0.0],
        altitude=[0.0],
        layer_ids=(:ground,),
    )
    source = AOS.Optics.Source(
        band=:custom,
        wavelength=Float64(raw["science_wavelength_m"]),
        photon_irradiance=Float64(raw["photon_irradiance"]),
    )
    definition = AOSPlant.OpticalPathDefinition(
        :finite_support,
        source,
        Gate5PathModel(),
    )
    path = AOSPlant.prepare_path_executor(
        Gate5PathModel(),
        definition,
        source,
        telescope,
        atmosphere,
    )
    pupil = AOSPlant.path_input(path)
    surface = fill(2.0, 3, 3)
    metadata = sampled_surface_metadata(pupil, surface)
    coupling = AOSPlant.prepare_sampled_pupil_footprint_coupling(
        metadata,
        surface,
        path,
        AOSPlant.PupilPlanePlacement(),
    )
    fill!(pupil.opd, 7.0)
    AOSPlant.apply_sampled_pupil_surface!(
        pupil,
        surface,
        coupling,
        AOS.Optics.DMReplace(),
    )
    expected = zeros(5, 5)
    expected[2:4, 2:4] .= 2.0
    maximum_error = maximum(abs, pupil.opd .- expected)
    tolerance = Float64(raw["finite_support_atol"])
    maximum_error <= tolerance ||
        error("Gate 5 finite-support oracle exceeded tolerance")
    return Dict{String,Any}(
        "maximum_error_m" => maximum_error,
        "tolerance" => tolerance,
        "supported_samples" => count(!iszero, pupil.opd),
        "zeroed_samples" => count(iszero, pupil.opd),
        "passed" => true,
    )
end

function storage_signature(operation::Gate5OpticalOperation)
    path_storage = map(AOSPlant.prepared_paths(operation.plant)) do path
        (
            objectid(AOSPlant.path_input(path).opd),
            objectid(AOSPlant.path_result(path).values),
        )
    end
    optic_storage = map(
        operation.state.controllable_optics,
        operation.workspace.controllable_optics,
    ) do state, workspace
        (
            objectid(AOS.Optics.surface_opd(state.active)),
            objectid(AOS.Optics.surface_opd(workspace.staged)),
        )
    end
    return (
        path_storage,
        optic_storage,
        objectid(operation.state.scheduler.cursors),
        objectid(operation.state.command_endpoints),
        objectid(operation.state.command_applications),
        objectid(operation.state.controllable_optics),
        objectid(operation.state.path_sampled),
        objectid(operation.workspace.scheduler.due_slots),
        objectid(operation.workspace.command_endpoints),
        objectid(operation.workspace.controllable_optics),
        objectid(operation.workspace.command_dispositions),
        objectid(operation.workspace.due_paths),
        Base.summarysize(operation.event_loop),
        Base.summarysize(operation.state),
        Base.summarysize(operation.workspace),
    )
end

function validate_fixed_storage(
    raw::AbstractDict,
    warmup_cycles::Integer,
    long_run_cycles::Integer,
)
    operation, expected, _ = prepare_gate5_operation(raw)
    run_gate5_cycles!(operation, warmup_cycles)
    before = storage_signature(operation)
    hashes_before = path_hashes(operation, expected)
    run_gate5_cycles!(operation, long_run_cycles)
    after = storage_signature(operation)
    hashes_after = path_hashes(operation, expected)
    before == after ||
        error("Gate 5 prepared storage changed during the long run")
    hashes_before == hashes_after ||
        error("Gate 5 deterministic optical products changed during replay")
    return Dict{String,Any}(
        "warmup_cycles" => warmup_cycles,
        "long_run_cycles" => long_run_cycles,
        "storage_signature_stable" => true,
        "deterministic_product_hashes_stable" => true,
        "final_cycle" => Int64(operation.cycles),
        "prepared_event_loop_bytes" => Base.summarysize(
            operation.event_loop,
        ),
        "event_state_bytes" => Base.summarysize(operation.state),
        "event_workspace_bytes" => Base.summarysize(operation.workspace),
    )
end

function topology_characterization(
    raw::AbstractDict,
    path_counts,
    resolution::Integer,
)
    results = Vector{Dict{String,Any}}(undef, length(path_counts))
    for (index, path_count_value) in enumerate(path_counts)
        path_count = Int(path_count_value)
        GC.gc()
        timed = @timed prepare_gate5_operation(
            raw;
            path_count,
            resolution,
        )
        operation, _, science_count = timed.value
        controllable_bindings =
            AOSPlant.prepared_controllable_optic_path_bindings(
                operation.plant,
            )
        sampled_bindings =
            AOSPlant.prepared_sampled_aberration_path_bindings(
                operation.plant,
            )
        common_optic_count = 3
        expected_controllable_bindings =
            common_optic_count * path_count + science_count
        expected_sampled_bindings = path_count + science_count
        actual_controllable_bindings =
            AOSPlant.prepared_controllable_optic_binding_count(
                controllable_bindings,
            )
        actual_sampled_bindings =
            AOSPlant.prepared_sampled_aberration_binding_count(
                sampled_bindings,
            )
        actual_controllable_bindings == expected_controllable_bindings ||
            error("Gate 5 topology produced an unexpected optic binding count")
        actual_sampled_bindings == expected_sampled_bindings ||
            error("Gate 5 topology produced an unexpected aberration binding count")
        operation()
        results[index] = Dict{String,Any}(
            "path_count" => path_count,
            "wfs_path_count" => path_count - science_count,
            "science_path_count" => science_count,
            "controllable_optic_count" =>
                length(AOSPlant.prepared_controllable_optics(
                    operation.plant,
                )),
            "sampled_aberration_count" =>
                length(AOSPlant.prepared_sampled_aberrations(
                    operation.plant,
                )),
            "controllable_binding_count" =>
                actual_controllable_bindings,
            "expected_controllable_binding_count" =>
                expected_controllable_bindings,
            "sampled_aberration_binding_count" =>
                actual_sampled_bindings,
            "expected_sampled_aberration_binding_count" =>
                expected_sampled_bindings,
            "plane_group_count" =>
                AOSPlant.prepared_controllable_optic_plane_group_count(
                    controllable_bindings,
                ),
            "preparation_elapsed_ns" =>
                round(Int64, timed.time * 1.0e9),
            "preparation_allocated_bytes" => Int64(timed.bytes),
            "preparation_gc_ns" =>
                round(Int64, timed.gctime * 1.0e9),
            "prepared_summary_bytes" => Base.summarysize((
                operation.plant,
                operation.event_loop,
                operation.state,
                operation.workspace,
            )),
            "first_cycle_completed" => true,
            "counts_passed" => true,
        )
    end
    return results
end

end # module
