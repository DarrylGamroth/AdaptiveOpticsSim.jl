#
# Exact-target validation for one prepared Plant graph
#
# Preparation is a cold boundary. Every model extension must explicitly
# validate the target-local numerical storage it contributes; unknown prepared
# implementations fail closed. Host configuration, fixed registries, scalar
# values, and the RNG subsystem's documented host staging remain excluded.
#

@noinline function _throw_wrong_plant_target(
    target::AbstractComputeDevice,
    label::AbstractString,
    actual::AbstractComputeDevice,
)
    _throw_compute_device_error(
        :validate_prepared_plant_target,
        :wrong_device,
        target,
        "$label occupies $(actual)",
    )
end

@inline function _require_exact_plant_array_target(
    storage::AbstractArray,
    target::AbstractComputeDevice,
    label::AbstractString,
)
    actual = compute_device(storage)
    actual == target || _throw_wrong_plant_target(target, label, actual)
    return storage
end

@inline function _require_exact_plant_metadata_target(
    metadata::Union{
        OpticalPlaneMetadata,
        WFSObservationMetadata,
        WFSMeasurementMetadata,
    },
    target::AbstractComputeDevice,
    label::AbstractString,
)
    actual = metadata.device
    actual == target ||
        _throw_wrong_plant_target(target, "$label metadata", actual)
    return metadata
end

function _require_exact_plant_product_target(
    product::PupilFunction,
    target::AbstractComputeDevice,
    label::AbstractString,
)
    _validate_path_input(product)
    _require_exact_plant_metadata_target(product.metadata, target, label)
    _require_exact_plant_array_target(
        product.support, target, "$label support")
    _require_exact_plant_array_target(
        product.amplitude, target, "$label amplitude")
    _require_exact_plant_array_target(product.opd, target, "$label OPD")
    return product
end

function _require_exact_plant_product_target(
    product::Union{ElectricField,IntensityMap},
    target::AbstractComputeDevice,
    label::AbstractString,
)
    _validate_path_input(product)
    _require_exact_plant_metadata_target(product.metadata, target, label)
    _require_exact_plant_array_target(
        product.values, target, "$label storage")
    return product
end

function _require_exact_plant_product_target(
    product::WFSObservation,
    target::AbstractComputeDevice,
    label::AbstractString,
)
    validate_wfs_observation(product)
    _require_exact_plant_wfs_storage_target(
        product.metadata, product.storage, target, label)
    return product
end

function _require_exact_plant_product_target(
    product::WFSMeasurement,
    target::AbstractComputeDevice,
    label::AbstractString,
)
    validate_wfs_measurement(product)
    _require_exact_plant_wfs_storage_target(
        product.metadata, product.storage, target, label)
    return product
end

@inline function _require_exact_plant_wfs_storage_target(
    metadata,
    storage::AbstractArray,
    target::AbstractComputeDevice,
    label::AbstractString,
)
    _require_exact_plant_metadata_target(metadata, target, label)
    _require_exact_plant_array_target(storage, target, "$label storage")
    return storage
end

@inline function _require_exact_plant_wfs_storage_target(
    metadata,
    storage::Base.RefValue,
    ::AbstractComputeDevice,
    ::AbstractString,
)
    # The ordinary WFS validator has already bound this metadata and scalar
    # storage to the canonical host domain.
    return storage
end

@inline _require_exact_plant_product_target(
    ::Nothing, ::AbstractComputeDevice, ::AbstractString) = nothing

@inline function _require_exact_plant_product_target(
    value::Base.RefValue,
    ::AbstractComputeDevice,
    ::AbstractString,
)
    acquisition_product_contract(value)
    return value
end

@inline _require_exact_plant_product_target(
    ::Tuple{}, ::AbstractComputeDevice, ::AbstractString) = nothing

@inline function _require_exact_plant_product_target(
    products::Tuple,
    target::AbstractComputeDevice,
    label::AbstractString,
)
    _require_exact_plant_product_target(first(products), target, label)
    _require_exact_plant_product_target(Base.tail(products), target, label)
    return products
end

function _require_exact_plant_product_target(
    products::_FixedOpticalProductVector,
    target::AbstractComputeDevice,
    label::AbstractString,
)
    @inbounds for product in products
        _require_exact_plant_product_target(product, target, label)
    end
    return products
end

@inline function _require_exact_plant_product_target(
    bundle::OpticalProductBundle,
    target::AbstractComputeDevice,
    label::AbstractString,
)
    _require_exact_plant_product_target(bundle.products, target, label)
    return bundle
end

function _require_exact_plant_product_target(
    storage::AbstractArray,
    target::AbstractComputeDevice,
    label::AbstractString,
)
    return _require_exact_plant_array_target(storage, target, label)
end

function _require_exact_plant_product_target(
    product,
    ::AbstractComputeDevice,
    label::AbstractString,
)
    throw(PlantPreparationError(
        :plant,
        :unsupported_target_validation,
        "$label type $(typeof(product)) has no exact-target validator",
    ))
end

function _require_exact_acquisition_products_target(
    products::AcquisitionProducts,
    target::AbstractComputeDevice,
)
    _require_exact_plant_product_target(
        products.observation, target, "acquisition observation")
    _require_exact_plant_product_target(
        products.measurement, target, "acquisition measurement")
    return products
end

"""Qualified extension seam for target validation of a prepared path execution."""
function validate_path_execution_target(execution, ::AbstractComputeDevice)
    throw(PlantPreparationError(
        :path,
        :unsupported_target_validation,
        "prepared path execution $(typeof(execution)) does not validate " *
        "its exact-target storage",
    ))
end

function validate_path_execution_target(
    execution::Union{
        PreparedDirectImaging,
        PreparedIncoherentDirectImaging,
        PreparedBundledDirectImaging,
    },
    target::AbstractComputeDevice,
)
    _require_exact_direct_imaging_target(execution, target)
    return execution
end

function validate_path_execution_target(
    execution::WFSOpticalPathExecution,
    target::AbstractComputeDevice,
)
    validate_wfs_target(execution.plan, target)
    return execution
end

"""Qualified extension seam for exact-target path materialization validation."""
function validate_path_materialization_target(
    materialization, input, atmosphere, ::AbstractComputeDevice)
    throw(PlantPreparationError(
        :path,
        :unsupported_target_validation,
        "prepared path materialization $(typeof(materialization)) does not " *
        "validate its exact-target storage",
    ))
end

@inline function validate_path_materialization_target(
    materialization::AtmosphereIndependentPath,
    input,
    atmosphere,
    ::AbstractComputeDevice,
)
    validate_path_materialization_binding(
        materialization, input, atmosphere, nothing)
    return materialization
end

function validate_path_materialization_target(
    materialization::PreparedPupilOPDMaterialization,
    input::PupilFunction,
    atmosphere::AbstractTimedAtmosphere,
    target::AbstractComputeDevice,
)
    validate_path_materialization_binding(
        materialization, input, atmosphere, materialization.source)
    renderer = materialization.renderer
    _require_exact_plant_metadata_target(
        renderer.output_metadata, target, "atmosphere renderer output")
    _require_exact_plant_array_target(
        renderer.pupil, target, "atmosphere renderer pupil support")
    # Per-layer shift and cone-scale vectors are immutable host geometry used
    # to launch scalar-direction rendering. They are configuration, not
    # accelerator data-plane storage.
    return materialization
end

"""Qualified extension seam for target validation of an illumination evaluator."""
function validate_illumination_evaluator_target(
    evaluator, ::AbstractComputeDevice)
    throw(PlantPreparationError(
        :illumination,
        :unsupported_target_validation,
        "prepared illumination evaluator $(typeof(evaluator)) does not " *
        "validate its exact-target storage",
    ))
end

function validate_illumination_evaluator_target(
    evaluator::PreparedUniformIntensityIllumination,
    target::AbstractComputeDevice,
)
    evaluator.device == target || _throw_wrong_plant_target(
        target, "uniform illumination evaluator", evaluator.device)
    return evaluator
end

function validate_path_materialization_target(
    entry::PreparedIlluminationEntry,
    input,
    atmosphere,
    target::AbstractComputeDevice,
)
    validate_illumination_entry_binding(entry, input)
    _require_exact_plant_product_target(
        illumination_destination(entry), target, "illumination destination")
    validate_illumination_evaluator_target(
        illumination_evaluator(entry), target)
    return entry
end

"""Qualified extension seam for exact-target controllable-optic preparation."""
function validate_controllable_optic_target(
    implementation, ::AbstractComputeDevice)
    throw(PlantPreparationError(
        :controllable_optic,
        :unsupported_target_validation,
        "prepared controllable optic $(typeof(implementation)) does not " *
        "validate its exact-target storage",
    ))
end

@inline _require_exact_dm_influence_target(
    ::Union{GaussianInfluenceWidth,GaussianMechanicalCoupling},
    ::AbstractComputeDevice,
) = nothing

function _require_exact_dm_influence_target(
    influence::Union{DenseInfluenceMatrix,MeasuredInfluenceFunctions},
    target::AbstractComputeDevice,
)
    _require_exact_plant_array_target(
        influence.modes, target, "deformable-mirror influence matrix")
    return influence
end

function _require_exact_dm_influence_target(
    influence::GaussianInfluenceOperator,
    target::AbstractComputeDevice,
)
    _require_exact_plant_array_target(
        influence.pupil_backend, target,
        "deformable-mirror Gaussian pupil support")
    _require_exact_plant_array_target(
        influence.coordinates_backend, target,
        "deformable-mirror Gaussian actuator coordinates")
    # `pupil_host` and `coordinates_host` are deliberate CPU mirrors used by
    # scalar indexing and inspection; the backend copies are authoritative for
    # accelerator execution.
    return influence
end

function _require_exact_dm_influence_target(
    modes::AbstractMatrix,
    target::AbstractComputeDevice,
)
    _require_exact_plant_array_target(
        modes, target, "deformable-mirror influence storage")
    return modes
end

function _require_exact_dm_influence_target(
    influence::AbstractDMInfluenceModel,
    ::AbstractComputeDevice,
)
    throw(PlantPreparationError(
        :controllable_optic,
        :unsupported_target_validation,
        "deformable-mirror influence model $(typeof(influence)) has no " *
        "exact-target validator",
    ))
end

@inline _require_exact_dm_actuator_target(
    ::Union{LinearStaticActuators,ClippedActuators},
    ::AbstractComputeDevice,
) = nothing

function _require_exact_dm_actuator_target(
    model::ActuatorHealthMap,
    target::AbstractComputeDevice,
)
    _require_exact_plant_array_target(
        model.gains, target, "deformable-mirror actuator health gains")
    return model
end

@inline _require_exact_dm_actuator_stages_target(
    ::Tuple{}, ::AbstractComputeDevice) = nothing

@inline function _require_exact_dm_actuator_stages_target(
    stages::Tuple,
    target::AbstractComputeDevice,
)
    _require_exact_dm_actuator_target(first(stages), target)
    _require_exact_dm_actuator_stages_target(Base.tail(stages), target)
    return stages
end

function _require_exact_dm_actuator_target(
    model::CompositeDMActuatorModel,
    target::AbstractComputeDevice,
)
    _require_exact_dm_actuator_stages_target(model.stages, target)
    return model
end

function _require_exact_dm_actuator_target(
    model::AbstractDMActuatorModel,
    ::AbstractComputeDevice,
)
    throw(PlantPreparationError(
        :controllable_optic,
        :unsupported_target_validation,
        "deformable-mirror actuator model $(typeof(model)) has no " *
        "exact-target validator",
    ))
end

@inline _require_exact_optional_plant_array_target(
    ::Nothing, ::AbstractComputeDevice, ::AbstractString) = nothing

@inline function _require_exact_optional_plant_array_target(
    storage::AbstractArray,
    target::AbstractComputeDevice,
    label::AbstractString,
)
    return _require_exact_plant_array_target(storage, target, label)
end

function validate_controllable_optic_target(
    implementation::_PreparedPlantDeformableMirror,
    target::AbstractComputeDevice,
)
    _require_exact_plant_metadata_target(
        implementation.surface_metadata, target,
        "deformable-mirror sampled surface")
    _require_exact_dm_influence_target(implementation.modes, target)
    _require_exact_optional_plant_array_target(
        implementation.separable_x, target,
        "deformable-mirror separable x factor")
    _require_exact_optional_plant_array_target(
        implementation.separable_y_t, target,
        "deformable-mirror separable y factor")
    _require_exact_dm_influence_target(
        implementation.params.influence_model, target)
    _require_exact_dm_actuator_target(
        implementation.params.actuator_model, target)
    # Prepared topology coordinates, validity masks, indices, and metadata are
    # immutable host configuration. Runtime Gaussian and dense influence
    # storage above is the target-local numerical representation.
    return implementation
end

@inline function validate_controllable_optic_target(
    implementation::PreparedCircularPyramidModulator,
    ::AbstractComputeDevice,
)
    # The prepared modulator contains only scalar waveform configuration and
    # endpoint identities; its mutable state and workspace are likewise
    # scalar host control-plane values.
    return implementation
end

"""Qualified extension seam for exact-target controllable-optic state."""
function validate_controllable_optic_state_target(
    implementation, state, ::AbstractComputeDevice)
    throw(PlantPreparationError(
        :controllable_optic,
        :unsupported_target_validation,
        "runtime state $(typeof(state)) for controllable optic " *
        "$(typeof(implementation)) has no exact-target validator",
    ))
end

"""Qualified extension seam for exact-target controllable-optic workspace."""
function validate_controllable_optic_workspace_target(
    implementation, workspace, ::AbstractComputeDevice)
    throw(PlantPreparationError(
        :controllable_optic,
        :unsupported_target_validation,
        "runtime workspace $(typeof(workspace)) for controllable optic " *
        "$(typeof(implementation)) has no exact-target validator",
    ))
end

@inline function _require_target_local_command_values_target(
    ::_TargetLocalScalarEffectiveCommandValues,
    ::AbstractComputeDevice,
)
    return nothing
end

function _require_target_local_command_values_target(
    values::_TargetLocalArrayEffectiveCommandValues,
    target::AbstractComputeDevice,
)
    _require_exact_plant_array_target(
        values.active, target, "target-local active effective command")
    _require_exact_plant_array_target(
        values.staging, target, "target-local staged effective command")
    values.active === values.staging && throw(PlantPreparationError(
        :command_replica,
        :aliased_staging,
        "target-local active and staging command storage must be distinct",
    ))
    return nothing
end

function validate_target_local_command_endpoint_target(
    endpoint::PreparedTargetLocalCommandEndpoint,
    state::TargetLocalCommandEndpointState,
    target::AbstractComputeDevice,
)
    endpoint.binding === state.binding || throw(PlantPreparationError(
        :command_replica,
        :foreign_endpoint_state,
        "target-local command state belongs to another prepared endpoint",
    ))
    compute_device(endpoint) == target || _throw_wrong_plant_target(
        target,
        "target-local command endpoint",
        compute_device(endpoint),
    )
    _require_target_local_command_values_target(state.values, target)
    return state
end

@inline function _require_exact_dm_runtime_modes_target(
    modes::AbstractArray,
    target::AbstractComputeDevice,
    label::AbstractString,
)
    return _require_exact_plant_array_target(modes, target, label)
end

function _require_exact_dm_runtime_modes_target(
    operator::GaussianInfluenceOperator,
    target::AbstractComputeDevice,
    label::AbstractString,
)
    _require_exact_plant_array_target(
        operator.pupil_backend, target, "$label pupil storage")
    _require_exact_plant_array_target(
        operator.coordinates_backend, target, "$label coordinate storage")
    # The corresponding host arrays are deliberate immutable geometry and
    # scalar-indexing storage for the lazy AbstractMatrix interface.
    return operator
end

function _require_exact_deformable_mirror_runtime_target(
    mirror::DeformableMirror,
    target::AbstractComputeDevice,
    label::AbstractString,
)
    state = mirror.state
    for (field_label, storage) in (
        ("OPD", state.opd),
        ("OPD vector", state.opd_vec),
        ("command coefficients", state.coefs),
        ("actuator coefficients", state.actuator_coefs),
    )
        _require_exact_plant_array_target(
            storage, target, "$label $field_label")
    end
    _require_exact_dm_runtime_modes_target(
        state.modes, target, "$label influence modes")
    _require_exact_optional_plant_array_target(
        state.coefs_grid, target, "$label coefficient grid")
    _require_exact_optional_plant_array_target(
        state.separable_x, target, "$label separable x factor")
    _require_exact_optional_plant_array_target(
        state.separable_y_t, target, "$label separable y factor")
    _require_exact_optional_plant_array_target(
        state.separable_tmp, target, "$label separable workspace")
    _require_exact_dm_influence_target(
        mirror.params.influence_model, target)
    _require_exact_dm_actuator_target(mirror.params.actuator_model, target)
    return mirror
end

function validate_controllable_optic_state_target(
    implementation::_PreparedPlantDeformableMirror,
    state::_PlantDeformableMirrorState,
    target::AbstractComputeDevice,
)
    _require_exact_deformable_mirror_runtime_target(
        state.active, target, "active deformable mirror")
    return state
end

function validate_controllable_optic_workspace_target(
    implementation::_PreparedPlantDeformableMirror,
    workspace::_PlantDeformableMirrorWorkspace,
    target::AbstractComputeDevice,
)
    _require_exact_deformable_mirror_runtime_target(
        workspace.staged, target, "staged deformable mirror")
    return workspace
end

@inline function validate_controllable_optic_state_target(
    ::PreparedCircularPyramidModulator,
    state::CircularPyramidModulatorState,
    ::AbstractComputeDevice,
)
    return state
end

@inline function validate_controllable_optic_workspace_target(
    ::PreparedCircularPyramidModulator,
    workspace::CircularPyramidModulatorWorkspace,
    ::AbstractComputeDevice,
)
    return workspace
end

"""Qualified extension seam for target validation of a pupil-surface coupling."""
function validate_pupil_surface_coupling_target(
    coupling::AbstractPupilSurfacePathCoupling,
    ::AbstractComputeDevice,
)
    throw(PlantPreparationError(
        :pupil_surface_coupling,
        :unsupported_target_validation,
        "prepared coupling $(typeof(coupling)) does not validate its " *
        "exact-target storage",
    ))
end

"""
Qualified extension seam for exact-target autonomous-optic coupling
validation.
"""
function validate_autonomous_optic_coupling_target(
    coupling, ::AbstractComputeDevice)
    throw(PlantPreparationError(
        :autonomous_periodic_optic,
        :unsupported_target_validation,
        "prepared autonomous-optic coupling $(typeof(coupling)) does not " *
        "validate its exact-target storage",
    ))
end

function validate_autonomous_optic_coupling_target(
    coupling::PreparedCycleAveragedPyramidModulation,
    target::AbstractComputeDevice,
)
    _require_exact_plant_array_target(
        coupling.modulation.phases, target,
        "cycle-averaged Pyramid modulation phases")
    # Quadrature weights and the circular policy are immutable host
    # configuration; the phase cube is the target-authoritative data plane.
    return coupling
end

@inline validate_pupil_surface_coupling_target(
    coupling::_NoPupilSurfacePathCoupling,
    ::AbstractComputeDevice,
) = coupling

function validate_pupil_surface_coupling_target(
    coupling::PreparedDirectPupilSurfaceCoupling,
    target::AbstractComputeDevice,
)
    _require_exact_plant_product_target(
        coupling.destination, target, "pupil-surface coupling destination")
    return coupling
end

function validate_pupil_surface_coupling_target(
    coupling::Union{
        PreparedIdentityPupilFootprintCoupling,
        PreparedPupilFootprintCoupling,
    },
    target::AbstractComputeDevice,
)
    _require_exact_plant_product_target(
        coupling.destination, target, "pupil-footprint coupling destination")
    _require_exact_plant_metadata_target(
        coupling.surface_metadata, target, "pupil-footprint surface")
    return coupling
end

function _require_exact_command_endpoint_target(
    binding::_PreparedPlantCommandEndpoint,
    target::AbstractComputeDevice,
)
    endpoint = binding.endpoint
    schema = command_schema(endpoint)
    return _require_exact_command_endpoint_target(
        binding, endpoint, schema, target)
end

@inline _require_host_scalar_command_value(
    value::Real, ::AbstractString) = value

function _require_host_scalar_command_value(value, label::AbstractString)
    throw(PlantPreparationError(
        :command_endpoint,
        Symbol(label),
        "scalar command endpoint $label must be a host real scalar; got " *
        "$(typeof(value))",
    ))
end

@inline _require_optional_host_scalar_command_value(
    ::Nothing, ::AbstractString) = nothing
@inline _require_optional_host_scalar_command_value(
    value::Real, ::AbstractString) = value

function _require_optional_host_scalar_command_value(
    value, label::AbstractString)
    throw(PlantPreparationError(
        :command_endpoint,
        Symbol(label),
        "scalar command endpoint $label must be a host real scalar or " *
        "nothing; got $(typeof(value))",
    ))
end

function _require_exact_command_endpoint_target(
    binding::_PreparedPlantCommandEndpoint,
    endpoint::PreparedCommandEndpoint{<:PlantCommandSchema{<:Any,0}},
    schema::PlantCommandSchema,
    ::AbstractComputeDevice,
)
    typeof(backend(endpoint)) === CPUBackend || throw(
        PlantPreparationError(
            :command_endpoint,
            :backend,
            "scalar command endpoint storage must remain host-resident",
        ))
    _require_host_scalar_command_value(
        initial_effective_command(binding), "initial_command")
    _require_optional_host_scalar_command_value(
        safe_effective_command(binding), "safe_command")
    return binding
end

function _require_exact_command_endpoint_target(
    binding::_PreparedPlantCommandEndpoint,
    endpoint::PreparedCommandEndpoint{<:PlantCommandSchema{<:Any,N}},
    schema::PlantCommandSchema,
    target::AbstractComputeDevice,
) where {N}
    typeof(backend(endpoint)) === typeof(compute_device_backend(target)) ||
        throw(PlantPreparationError(
            :command_endpoint,
            :backend,
            "array command endpoint backend does not derive from the exact " *
            "plant target",
        ))
    initial = initial_effective_command(binding)
    size(initial) == command_dimensions(schema) || throw(
        PlantPreparationError(
            :command_endpoint,
            :initial_command,
            "array command endpoint initial value has the wrong shape",
        ))
    _require_exact_plant_array_target(
        initial, target, "array command endpoint initial value")
    safe = safe_effective_command(binding)
    if safe !== nothing
        size(safe) == command_dimensions(schema) || throw(
            PlantPreparationError(
                :command_endpoint,
                :safe_command,
                "array command endpoint safe value has the wrong shape",
            ))
        _require_exact_plant_array_target(
            safe, target, "array command endpoint safe value")
    end
    return binding
end

function _require_exact_sampled_aberration_target(
    aberration::PreparedSampledAberration,
    target::AbstractComputeDevice,
)
    _validate_prepared_sampled_aberration_storage(
        aberration.metadata, aberration.opd, target)
    return aberration
end

function _require_exact_path_target(
    path::PreparedPathExecutor,
    telescope::AbstractTelescope,
    atmosphere::AbstractTimedAtmosphere,
    target::AbstractComputeDevice,
    context,
)
    path.telescope === telescope || throw(PlantPreparationError(
        :path, :prepared_binding,
        "prepared path does not retain the plant telescope"))
    path.atmosphere === atmosphere || throw(PlantPreparationError(
        :path, :prepared_binding,
        "prepared path does not retain the authoritative plant atmosphere"))
    path.context === context || throw(PlantPreparationError(
        :path, :prepared_binding,
        "prepared path does not retain the authoritative plant execution context"))
    getfield(path.key, :device) == target || _throw_wrong_plant_target(
        target, "path-result key", getfield(path.key, :device))
    typeof(getfield(path.key, :backend)) ===
        typeof(compute_device_backend(target)) || throw(
        PlantPreparationError(
            :path,
            :backend,
            "prepared path backend does not derive from the exact plant target",
        ))
    _require_exact_plant_product_target(
        path.input, target, "prepared path input")
    _require_exact_plant_product_target(
        path.result, target, "prepared path result")
    validate_path_materialization_target(
        path.materialization, path.input, atmosphere, target)
    validate_path_execution_binding(path.execution, path.input, path.result)
    validate_path_execution_target(path.execution, target)
    return path
end

"""Qualified extension seam for target validation of acquisition execution."""
function validate_acquisition_execution_target(
    execution, ::AbstractComputeDevice)
    throw(PlantPreparationError(
        :acquisition,
        :unsupported_target_validation,
        "prepared acquisition execution $(typeof(execution)) does not " *
        "validate its exact-target storage",
    ))
end

function validate_acquisition_execution_target(
    execution::FrameAcquisitionExecution,
    target::AbstractComputeDevice,
)
    _require_exact_detector_acquisition_target(
        execution.detector, execution.plan, target)
    _require_exact_plant_array_target(
        execution.observation, target, "frame acquisition observation")
    return execution
end

function validate_acquisition_execution_target(
    execution::WFSAcquisitionExecution,
    target::AbstractComputeDevice,
)
    validate_wfs_target(execution.acquisition, target)
    _require_exact_plant_product_target(
        execution.observation, target, "WFS acquisition observation")
    if execution.estimator !== nothing
        validate_wfs_target(execution.estimator, target)
        _require_exact_plant_product_target(
            execution.measurement, target, "WFS acquisition measurement")
    end
    return execution
end

"""Qualified extension seam for exact-target acquisition-provider validation."""
function validate_acquisition_provider_target(
    implementation, path_result, products, ::AbstractComputeDevice)
    throw(PlantPreparationError(
        :acquisition,
        :unsupported_target_validation,
        "prepared acquisition provider $(typeof(implementation)) does not " *
        "validate its exact-target storage",
    ))
end

function validate_acquisition_provider_target(
    implementation::PreparedFullOpticalProvider,
    path_result,
    products::AcquisitionProducts,
    target::AbstractComputeDevice,
)
    validate_acquisition_execution_target(implementation.execution, target)
    return implementation
end

@inline function validate_acquisition_provider_target(
    implementation::PreparedUnchangedSyntheticProvider,
    path_result,
    products::AcquisitionProducts,
    ::AbstractComputeDevice,
)
    return implementation
end

function validate_acquisition_provider_target(
    implementation::PreparedCopiedSyntheticProvider,
    path_result,
    products::AcquisitionProducts,
    target::AbstractComputeDevice,
)
    _require_exact_acquisition_products_target(
        implementation.source, target)
    return implementation
end

function validate_acquisition_provider_target(
    implementation::PreparedCyclicReplayProvider,
    path_result,
    products::AcquisitionProducts,
    target::AbstractComputeDevice,
)
    @inbounds for source in implementation.corpus
        _require_exact_acquisition_products_target(source, target)
    end
    return implementation
end

@inline _require_exact_reduced_order_responses_target(
    ::Tuple{}, ::AbstractComputeDevice) = nothing

@inline function _require_exact_reduced_order_responses_target(
    responses::Tuple,
    target::AbstractComputeDevice,
)
    response = first(responses)
    _require_exact_plant_array_target(
        response.operator, target, "reduced-order command response")
    _require_exact_reduced_order_responses_target(
        Base.tail(responses), target)
    return responses
end

function validate_acquisition_provider_target(
    implementation::PreparedLinearReducedOrderProvider,
    path_result,
    products::AcquisitionProducts,
    target::AbstractComputeDevice,
)
    implementation.device == target || _throw_wrong_plant_target(
        target, "reduced-order provider", implementation.device)
    disturbance = implementation.disturbance
    for (label, storage) in (
        ("reduced-order disturbance offsets", disturbance.offsets),
        ("reduced-order disturbance amplitudes", disturbance.amplitudes),
        ("reduced-order disturbance frequencies", disturbance.frequencies_hz),
        ("reduced-order disturbance phases", disturbance.phases_rad),
        ("reduced-order disturbance state", implementation.state.modes),
        ("reduced-order path projection", implementation.path_projection),
        ("reduced-order sensor operator", implementation.sensor_operator),
        ("reduced-order residual", implementation.residual),
        ("reduced-order command workspace", implementation.command_workspace),
    )
        _require_exact_plant_array_target(storage, target, label)
    end
    _require_exact_reduced_order_responses_target(
        implementation.command_responses, target)
    return implementation
end

function _require_exact_acquisition_target(
    owner::PreparedAcquisitionOwner,
    target::AbstractComputeDevice,
    context,
)
    owner.context === context || throw(PlantPreparationError(
        :acquisition, :prepared_binding,
        "prepared acquisition does not retain the authoritative plant execution context"))
    getfield(owner.path_key, :device) == target || _throw_wrong_plant_target(
        target, "acquisition path-result key",
        getfield(owner.path_key, :device))
    _require_exact_plant_product_target(
        owner.path_result, target, "acquisition path result")
    provider = owner.provider
    validate_acquisition_provider_binding(provider, owner.path_result)
    _require_exact_acquisition_products_target(provider.products, target)
    validate_acquisition_provider_target(
        provider.implementation, owner.path_result, provider.products, target)
    return owner
end

function _require_exact_coupling_bindings_target(
    bindings::PreparedControllableOpticPathBindings,
    ::AbstractComputeDevice,
)
    # This table contains only host-side path/group/optic ordinals. Exact
    # path-local controllable-optic couplings are prepared by the event loop.
    return bindings
end

function _require_exact_coupling_bindings_target(
    bindings::PreparedSampledAberrationPathBindings,
    target::AbstractComputeDevice,
)
    @inbounds for coupling in bindings.couplings
        validate_pupil_surface_coupling_target(coupling, target)
    end
    return bindings
end

function _require_exact_prepared_plant_target(
    plant::PreparedPlant,
    target::AbstractComputeDevice=plant.target,
)
    plant.target == target || _throw_wrong_plant_target(
        target, "prepared plant", plant.target)
    context_target =
        _prepared_device_execution_compute_device(plant.context)
    context_target == target || _throw_wrong_plant_target(
        target, "prepared plant execution context", context_target)
    validate_telescope_target(plant.telescope, target)
    validate_timed_atmosphere_target(plant.atmosphere, target)

    @inbounds for binding in plant.command_endpoints
        _require_exact_command_endpoint_target(binding, target)
    end
    @inbounds for optic in plant.controllable_optics
        validate_controllable_optic_target(optic.implementation, target)
    end
    @inbounds for aberration in plant.sampled_aberrations
        _require_exact_sampled_aberration_target(aberration, target)
    end
    @inbounds for path in plant.paths
        _require_exact_path_target(
            path, plant.telescope, plant.atmosphere, target, plant.context)
    end
    @inbounds for acquisition in plant.acquisitions
        _require_exact_acquisition_target(
            acquisition, target, plant.context)
    end
    _require_exact_coupling_bindings_target(
        plant.controllable_optic_path_bindings, target)
    _require_exact_coupling_bindings_target(
        plant.sampled_aberration_path_bindings, target)
    # Stable RNG owners, their host RNG state, and fixed Memory registries are
    # deliberate host control-plane storage and contain no target-local arrays.
    return plant
end

function _require_exact_device_path_batch_implementation_target(
    implementation::_PreparedDirectImagingDevicePathBatch,
    atmosphere::AbstractTimedAtmosphere,
    target::AbstractComputeDevice,
)
    implementation.key.device == target || _throw_wrong_plant_target(
        target, "direct-imaging device-batch key", implementation.key.device)
    _prepared_device_execution_compute_device(implementation.context) ==
        target || _throw_wrong_plant_target(
        target,
        "direct-imaging device-batch context",
        _prepared_device_execution_compute_device(implementation.context),
    )
    _validate_atmosphere_direction_batch_binding(
        implementation.atmosphere_batch, atmosphere)
    _require_exact_plant_metadata_target(
        implementation.atmosphere_batch.params.output_metadata,
        target,
        "direct-imaging atmosphere batch",
    )
    _require_exact_direct_imaging_target(implementation.optical_batch, target)
    @inbounds for output in implementation.atmosphere_outputs
        _require_exact_plant_array_target(
            output, target, "direct-imaging atmosphere batch output")
    end
    @inbounds for input in implementation.path_inputs
        _require_exact_plant_product_target(
            input, target, "direct-imaging batched path input")
    end
    @inbounds for result in implementation.path_results
        _require_exact_plant_product_target(
            result, target, "direct-imaging batched path result")
    end
    return implementation
end

function _require_exact_device_path_batch_implementation_target(
    implementation::_PreparedWFSDevicePathBatch,
    atmosphere::AbstractTimedAtmosphere,
    target::AbstractComputeDevice,
)
    implementation.key.device == target || _throw_wrong_plant_target(
        target, "WFS device-batch key", implementation.key.device)
    _prepared_device_execution_compute_device(implementation.context) ==
        target || _throw_wrong_plant_target(
        target,
        "WFS device-batch context",
        _prepared_device_execution_compute_device(implementation.context),
    )
    _validate_atmosphere_direction_batch_binding(
        implementation.atmosphere_batch, atmosphere)
    _require_exact_plant_metadata_target(
        implementation.atmosphere_batch.params.output_metadata,
        target,
        "WFS atmosphere batch",
    )
    @inbounds for output in implementation.atmosphere_outputs
        _require_exact_plant_array_target(
            output, target, "WFS atmosphere batch output")
    end
    @inbounds for input in implementation.path_inputs
        _require_exact_plant_product_target(
            input, target, "WFS batched path input")
    end
    @inbounds for result in implementation.path_results
        _require_exact_plant_product_target(
            result, target, "WFS batched path result")
    end
    return implementation
end

function _require_exact_device_path_batch_owner_target(
    owner::PreparedDevicePathBatchOwner,
    atmosphere::AbstractTimedAtmosphere,
    target::AbstractComputeDevice,
)
    owner.device == target || _throw_wrong_plant_target(
        target, "device path-batch owner", owner.device)
    _require_exact_device_path_batch_implementation_target(
        owner.implementation, atmosphere, target)
    return owner
end

function _require_exact_detector_lifecycle_target(
    lifecycle::Union{
        PreparedGlobalShutterAcquisition,
        PreparedRollingShutterAcquisition,
        PreparedFrameTransferAcquisition,
    },
    target::AbstractComputeDevice,
)
    _require_exact_detector_acquisition_target(
        lifecycle.detector, lifecycle.plan, target)
    _require_exact_detector_readout_products_target(
        lifecycle.readout_products, target)
    return lifecycle
end

@inline function _require_exact_acquisition_lifecycle_target(
    lifecycle::Union{
        PreparedGlobalShutterAcquisition,
        PreparedRollingShutterAcquisition,
    },
    target::AbstractComputeDevice,
)
    return _require_exact_detector_lifecycle_target(lifecycle, target)
end

function _require_exact_acquisition_lifecycle_target(
    lifecycle::PreparedFrameTransferAcquisition,
    target::AbstractComputeDevice,
)
    _require_exact_detector_lifecycle_target(lifecycle, target)
    _require_exact_plant_array_target(
        lifecycle.storage_frame, target, "frame-transfer storage frame")
    return lifecycle
end

function _require_exact_acquisition_lifecycle_target(
    lifecycle::PreparedDirectMeasurementAcquisition,
    target::AbstractComputeDevice,
)
    _require_exact_plant_product_target(
        lifecycle.measurement, target, "direct-measurement lifecycle output")
    _require_exact_plant_array_target(
        lifecycle.instantaneous_sample,
        target,
        "direct-measurement instantaneous sample",
    )
    _require_exact_plant_array_target(
        lifecycle.integrated_sample,
        target,
        "direct-measurement integrated sample",
    )
    return lifecycle
end

function _require_exact_acquisition_lifecycle_target(
    lifecycle::AbstractPreparedAcquisitionLifecycle,
    ::AbstractComputeDevice,
)
    throw(PlantPreparationError(
        :acquisition_lifecycle,
        :unsupported_target_validation,
        "prepared acquisition lifecycle $(typeof(lifecycle)) has no " *
        "exact-target validator",
    ))
end

@inline function _require_exact_command_payload_target(
    ::_ScalarCommandPayloadSlots,
    ::AbstractComputeDevice,
    ::AbstractString,
)
    return nothing
end

function _require_exact_command_payload_target(
    payloads::_ArrayCommandPayloadSlots,
    target::AbstractComputeDevice,
    label::AbstractString,
)
    @inbounds for values in payloads.values
        _require_exact_plant_array_target(values, target, label)
    end
    return payloads
end

@inline function _require_exact_command_application_target(
    ::_ScalarCommandApplicationValues,
    ::AbstractComputeDevice,
    ::AbstractString,
)
    return nothing
end

function _require_exact_command_application_target(
    values::_ArrayCommandApplicationValues,
    target::AbstractComputeDevice,
    label::AbstractString,
)
    _require_exact_plant_array_target(
        values.effective, target, "$label effective value")
    _require_exact_plant_array_target(
        values.staging, target, "$label staging value")
    _require_exact_optional_plant_array_target(
        values.safe, target, "$label safe value")
    return values
end

function _require_exact_plant_event_loop_state_target(
    prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState,
    target::AbstractComputeDevice=prepared.target,
)
    _require_plant_event_loop_binding(prepared, state)
    @inbounds for index in eachindex(prepared.command_endpoints)
        endpoint = prepared.command_endpoints[index].binding.endpoint
        endpoint_state = state.command_endpoints[index]
        _require_command_endpoint_binding(endpoint, endpoint_state)
        _require_exact_command_payload_target(
            endpoint_state.payloads,
            target,
            "command endpoint payload slot",
        )
        application = state.command_applications[index]
        _require_command_application_binding(
            endpoint, endpoint_state, application)
        _require_exact_command_application_target(
            application.values,
            target,
            "command application",
        )
    end
    @inbounds for index in eachindex(prepared.optics)
        validate_controllable_optic_state_target(
            prepared.optics[index].implementation,
            state.controllable_optics[index],
            target,
        )
    end
    # Scheduler, acquisition lifecycle, trigger, sequence, and timestamp state
    # contain host control-plane values only.
    return state
end

function _require_exact_plant_event_loop_workspace_target(
    prepared::PreparedPlantEventLoop,
    workspace::PlantEventLoopWorkspace,
    target::AbstractComputeDevice=prepared.target,
)
    _require_plant_event_loop_binding(prepared, workspace)
    @inbounds for index in eachindex(prepared.optics)
        validate_controllable_optic_workspace_target(
            prepared.optics[index].implementation,
            workspace.controllable_optics[index],
            target,
        )
    end
    # Scheduler, command disposition, transaction, due-path, batch-claim, and
    # trigger workspaces are fixed host control-plane storage.
    return workspace
end

function _require_exact_prepared_event_loop_target(
    prepared::PreparedPlantEventLoop,
    target::AbstractComputeDevice=prepared.target,
)
    prepared.target == target || _throw_wrong_plant_target(
        target, "prepared plant event loop", prepared.target)
    context_target =
        _prepared_device_execution_compute_device(prepared.context)
    context_target == target || _throw_wrong_plant_target(
        target, "prepared plant event-loop context", context_target)
    validate_timed_atmosphere_target(prepared.atmosphere, target)

    @inbounds for binding in prepared.command_endpoints
        _require_exact_command_endpoint_target(binding.binding, target)
    end
    @inbounds for optic in prepared.optics
        validate_controllable_optic_target(optic.implementation, target)
    end
    @inbounds for aberration in prepared.sampled_aberrations
        _require_exact_sampled_aberration_target(aberration, target)
    end
    @inbounds for group in prepared.path_groups
        group.requirements.compute_device == target ||
            _throw_wrong_plant_target(
                target,
                "event-loop path execution group",
                group.requirements.compute_device,
            )
        _require_exact_path_target(
            group.path, group.path.telescope, prepared.atmosphere, target,
            prepared.context)
        @inbounds for coupling in group.optic_couplings
            validate_pupil_surface_coupling_target(coupling, target)
        end
    end
    @inbounds for acquisition in prepared.acquisitions
        _require_exact_acquisition_lifecycle_target(
            acquisition.lifecycle, target)
        _require_exact_acquisition_products_target(
            acquisition.products, target)
        _require_exact_plant_product_target(
            acquisition.product, target, "event-loop acquisition product")
        if acquisition.sample_provider !== nothing
            provider = acquisition.sample_provider.provider
            validate_acquisition_provider_target(
                provider, nothing, acquisition.products, target)
        end
    end
    @inbounds for owner in prepared.device_path_batch_owners
        _require_exact_device_path_batch_owner_target(
            owner, prepared.atmosphere, target)
    end
    @inbounds for optic in prepared.autonomous_optics
        validate_controllable_optic_target(optic.implementation, target)
        validate_autonomous_optic_coupling_target(optic.coupling, target)
    end
    # Scheduler storage, topology ordinals, fixed registries, actions, trigger
    # state, and RNG owners are deliberate host control-plane data.
    return prepared
end
