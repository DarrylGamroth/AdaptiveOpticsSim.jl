#
# Target-local preparation resources
#
# These records are preparation-only leaves for a caller-resolved target
# partition.  They deliberately retain neither an atmosphere/materialization
# owner nor a lifecycle, schedule, or execution method.  A later partition
# composition layer remains responsible for binding authoritative atmosphere
# evolution and runtime ownership.
#

mutable struct _PreparedTargetLocalPathResourcesToken end
const _PREPARED_TARGET_LOCAL_PATH_RESOURCES_TOKEN =
    _PreparedTargetLocalPathResourcesToken()

"""
    PreparedTargetLocalPathResources

Immutable binding record for target-local optical-path resources prepared from
one exact definition, frozen source, telescope, and execution context.
`input`, `result`, and `execution` retain co-located static data and reusable
mutable workspace. The record contains no atmosphere authority,
materialization policy, timeline, schedule, or partition executor.
"""
struct PreparedTargetLocalPathResources{D,S,T,X,I,R,E,K<:PathResultKey}
    definition::D
    source::S
    telescope::T
    context::X
    rng_token::RNGOwnerToken
    input::I
    result::R
    execution::E
    key::K

    function PreparedTargetLocalPathResources(
        token::_PreparedTargetLocalPathResourcesToken,
        definition::D,
        source::S,
        telescope::T,
        context::X,
        input::I,
        result::R,
        execution::E,
        key::K,
    ) where {D,S,T,X,I,R,E,K<:PathResultKey}
        token === _PREPARED_TARGET_LOCAL_PATH_RESOURCES_TOKEN || throw(
            ArgumentError("invalid internal target-local path token"))
        return new{D,S,T,X,I,R,E,K}(
            definition,
            source,
            telescope,
            context,
            RNGOwnerToken(),
            input,
            result,
            execution,
            key,
        )
    end
end

@inline _rng_owner_binding_token(resources::PreparedTargetLocalPathResources) =
    resources.rng_token
@inline path_id(resources::PreparedTargetLocalPathResources) =
    path_id(resources.definition)

function PreparedTargetLocalPathResources(
    definition::OpticalPathDefinition,
    source::AbstractSource,
    telescope::AbstractTelescope,
    input,
    result,
    execution;
    context,
    optical_model,
    sampling_contract::AbstractOpticalSamplingContract=
        InstantaneousOpticalSample(),
    propagation_model,
    model_revisions=(),
)
    _validate_path_input(input)
    output_plane = _path_output_contract(result)
    first_result = _first_path_result(result)
    selector = backend(first_result)
    device = compute_device(first_result.values)
    _prepared_device_execution_compute_device(context) == device || throw(
        PlantPreparationError(
            :path,
            :execution_context_target,
            "prepared target-local path execution context does not match its exact compute device",
        ),
    )
    _require_path_result_domain(result, selector, device)
    _require_path_input_domain(input, selector, device)
    typeof(backend(telescope)) === typeof(selector) || throw(
        PlantPreparationError(
            :path,
            :backend,
            "prepared target-local path and telescope backends differ",
        ),
    )
    compute_device(pupil_reflectivity(telescope)) == device || throw(
        PlantPreparationError(
            :path,
            :device,
            "prepared target-local path and telescope occupy different compute devices",
        ),
    )
    revision = aperture_revision(telescope)
    _require_path_input_revisions(input, revision)
    validate_telescope_target(telescope, device)
    _require_exact_plant_product_target(
        input, device, "prepared target-local path input")
    _require_exact_plant_product_target(
        result, device, "prepared target-local path result")
    validate_path_execution_binding(execution, input, result)
    validate_path_execution_target(execution, device)
    revisions = (telescope=revision, model=model_revisions)
    key = PathResultKey(
        path_source_geometry_key(source),
        path_source_spectral_key(source),
        path_source_radiometry_key(source),
        optical_model,
        sampling_contract,
        propagation_model,
        output_plane,
        revisions,
        selector,
        device,
    )
    return PreparedTargetLocalPathResources(
        _PREPARED_TARGET_LOCAL_PATH_RESOURCES_TOKEN,
        definition,
        source,
        telescope,
        context,
        input,
        result,
        execution,
        key,
    )
end

@inline path_input(resources::PreparedTargetLocalPathResources) =
    resources.input
@inline path_result(resources::PreparedTargetLocalPathResources) =
    resources.result
@inline path_result_key(resources::PreparedTargetLocalPathResources) =
    resources.key

@inline _prepare_sampled_aberration_path_coupling(
    aberration::PreparedSampledAberration,
    path::PreparedTargetLocalPathResources,
) = _prepare_sampled_aberration_path_coupling_impl(aberration, path)

function prepare_controllable_optic_path_coupling(
    optic::PreparedTargetLocalControllableOptic,
    path::PreparedTargetLocalPathResources,
)
    path_target = getfield(path.key, :device)
    compute_device(optic) == path_target || throw(PlantPreparationError(
        :controllable_optic,
        :wrong_device,
        "target-local controllable optic $(controllable_optic_id(optic)) " *
        "occupies $(compute_device(optic)); path " *
        "$(path_id(path)) occupies $path_target",
    ))
    return prepare_controllable_optic_path_coupling(
        optic.implementation, optic.definition, path)
end

function require_path_result(path::PreparedTargetLocalPathResources;
    source_geometry=getfield(path.key, :source_geometry),
    spectral_sampling=getfield(path.key, :spectral_sampling),
    radiometry=getfield(path.key, :radiometry),
    optical_model=getfield(path.key, :optical_model),
    sampling_contract=getfield(path.key, :sampling_contract),
    propagation_model=getfield(path.key, :propagation_model),
    output_plane=getfield(path.key, :output_plane),
    revisions=getfield(path.key, :revisions),
    backend::AbstractArrayBackend=getfield(path.key, :backend),
    device::AbstractComputeDevice=getfield(path.key, :device))
    required = PathResultKey(source_geometry, spectral_sampling, radiometry,
        optical_model, sampling_contract, propagation_model, output_plane,
        revisions, backend, device)
    _require_path_result_key(path.key, required)
    return path
end

function _require_prepared_target_local_path_resources(
    resources::PreparedTargetLocalPathResources,
    definition::OpticalPathDefinition,
    source::AbstractSource,
    telescope::AbstractTelescope,
    context,
)
    resources.definition === definition || throw(PlantPreparationError(
        :path,
        :prepared_binding,
        "prepared target-local path does not retain its exact definition",
    ))
    resources.source === source || throw(PlantPreparationError(
        :path,
        :prepared_binding,
        "prepared target-local path does not retain its run-owned frozen source",
    ))
    resources.telescope === telescope || throw(PlantPreparationError(
        :path,
        :prepared_binding,
        "prepared target-local path does not retain its telescope",
    ))
    resources.context === context || throw(PlantPreparationError(
        :path,
        :prepared_binding,
        "prepared target-local path does not retain its execution context",
    ))
    _prepared_device_execution_compute_device(resources.context) ==
        getfield(resources.key, :device) || throw(PlantPreparationError(
        :path,
        :execution_context_target,
        "prepared target-local path execution context does not match its exact compute device",
    ))
    _validate_path_input(resources.input)
    target = getfield(resources.key, :device)
    validate_telescope_target(resources.telescope, target)
    _require_exact_plant_product_target(
        resources.input, target, "prepared target-local path input")
    _require_exact_plant_product_target(
        resources.result, target, "prepared target-local path result")
    _require_path_input_revisions(resources.input,
        getfield(getfield(resources.key, :revisions), :telescope))
    validate_path_execution_binding(resources.execution, resources.input,
        resources.result)
    validate_path_execution_target(resources.execution, target)
    return resources
end

function _require_prepared_target_local_path_resources(
    resources,
    ::OpticalPathDefinition,
    ::AbstractSource,
    ::AbstractTelescope,
    context,
)
    throw(PlantPreparationError(
        :path,
        :invalid_preparation,
        "target-local path model preparation must return PreparedTargetLocalPathResources; got $(typeof(resources))",
    ))
end

"""
    prepare_target_local_path_resources(definition, telescope, context)

Prepare one path's target-local optical resources.  Model extensions implement
the five-argument method below; this wrapper freezes the source and verifies
that the returned resources retain the exact caller-owned bindings.
"""
function prepare_target_local_path_resources(
    definition::OpticalPathDefinition,
    telescope::AbstractTelescope,
    context,
)
    source = freeze_source(path_source(definition))
    resources = _with_completed_prepared_device_execution_context(context) do
        prepare_target_local_path_resources(
            path_model(definition), definition, source, telescope, context)
    end
    return _require_prepared_target_local_path_resources(
        resources, definition, source, telescope, context)
end

"""Qualified cold extension seam for target-local optical-path resources."""
function prepare_target_local_path_resources(
    model,
    definition::OpticalPathDefinition,
    source::AbstractSource,
    telescope::AbstractTelescope,
    context,
)
    throw(PlantPreparationError(
        :path,
        :unsupported_model,
        "path model $(typeof(model)) does not implement prepare_target_local_path_resources",
    ))
end

mutable struct _PreparedTargetLocalAcquisitionResourcesToken end
const _PREPARED_TARGET_LOCAL_ACQUISITION_RESOURCES_TOKEN =
    _PreparedTargetLocalAcquisitionResourcesToken()

"""
    PreparedTargetLocalAcquisitionResources

Immutable binding record for target-local acquisition resources that borrow
one exact path result and retain one exact prepared provider with its mutable
detector or estimator state. It owns no event lifecycle, schedule, transfer,
or partition executor.
"""
struct PreparedTargetLocalAcquisitionResources{D,K<:PathResultKey,R,X,P}
    definition::D
    path_key::K
    path_result::R
    context::X
    rng_token::RNGOwnerToken
    provider::P

    function PreparedTargetLocalAcquisitionResources(
        token::_PreparedTargetLocalAcquisitionResourcesToken,
        definition::D,
        path_key::K,
        path_result::R,
        context::X,
        provider::P,
    ) where {D,K<:PathResultKey,R,X,P}
        token === _PREPARED_TARGET_LOCAL_ACQUISITION_RESOURCES_TOKEN || throw(
            ArgumentError("invalid internal target-local acquisition token"))
        return new{D,K,R,X,P}(
            definition,
            path_key,
            path_result,
            context,
            RNGOwnerToken(),
            provider,
        )
    end
end

@inline _rng_owner_binding_token(
    resources::PreparedTargetLocalAcquisitionResources,
) = resources.rng_token

@inline acquisition_provider(
    resources::PreparedTargetLocalAcquisitionResources,
) = resources.provider
@inline acquisition_provider_style(
    resources::PreparedTargetLocalAcquisitionResources,
) = acquisition_provider_style(resources.provider)
@inline acquisition_provider_payload_work(
    resources::PreparedTargetLocalAcquisitionResources,
) = acquisition_provider_payload_work(resources.provider)
@inline acquisition_product_contract(
    resources::PreparedTargetLocalAcquisitionResources,
) = acquisition_product_contract(resources.provider)
@inline acquisition_products(
    resources::PreparedTargetLocalAcquisitionResources,
) = acquisition_products(resources.provider)
@inline acquisition_observation(
    resources::PreparedTargetLocalAcquisitionResources,
) = acquisition_observation(resources.provider)
@inline acquisition_measurement(
    resources::PreparedTargetLocalAcquisitionResources,
) = acquisition_measurement(resources.provider)
@inline acquisition_product_metadata(
    resources::PreparedTargetLocalAcquisitionResources,
) = acquisition_product_metadata(resources.provider)

function PreparedTargetLocalAcquisitionResources(
    definition::AcquisitionDefinition,
    path::PreparedTargetLocalPathResources,
    provider::PreparedAcquisitionProvider,
)
    acquisition_path_id(definition) == path_id(path) || throw(
        PlantPreparationError(
            :acquisition,
            :unknown_path,
            "acquisition $(definition.id) does not reference prepared target-local path $(path_id(path))",
        ),
    )
    target = getfield(path.key, :device)
    _prepared_device_execution_compute_device(path.context) == target ||
        throw(PlantPreparationError(
            :acquisition,
            :execution_context_target,
            "prepared target-local acquisition context does not match its path-result device",
        ))
    _require_exact_plant_product_target(
        path.result, target, "target-local acquisition path result")
    validate_acquisition_provider_binding(provider, path.result)
    _require_exact_acquisition_products_target(provider.products, target)
    validate_acquisition_provider_target(
        provider.implementation, path.result, provider.products, target)
    return PreparedTargetLocalAcquisitionResources(
        _PREPARED_TARGET_LOCAL_ACQUISITION_RESOURCES_TOKEN,
        definition,
        path.key,
        path.result,
        path.context,
        provider,
    )
end

function _require_prepared_target_local_acquisition_resources(
    resources::PreparedTargetLocalAcquisitionResources,
    definition::AcquisitionDefinition,
    path::PreparedTargetLocalPathResources,
)
    resources.definition === definition || throw(PlantPreparationError(
        :acquisition,
        :prepared_binding,
        "prepared target-local acquisition does not retain its exact definition",
    ))
    resources.path_key === path.key || throw(PlantPreparationError(
        :acquisition,
        :prepared_binding,
        "prepared target-local acquisition does not retain the exact path-result key",
    ))
    resources.path_result === path.result || throw(PlantPreparationError(
        :acquisition,
        :prepared_binding,
        "prepared target-local acquisition does not retain the exact path result",
    ))
    resources.context === path.context || throw(PlantPreparationError(
        :acquisition,
        :prepared_binding,
        "prepared target-local acquisition does not retain the exact path execution context",
    ))
    target = getfield(resources.path_key, :device)
    _prepared_device_execution_compute_device(resources.context) == target ||
        throw(PlantPreparationError(
            :acquisition,
            :execution_context_target,
            "prepared target-local acquisition context does not match its path-result device",
        ))
    _require_exact_plant_product_target(
        resources.path_result, target, "target-local acquisition path result")
    validate_acquisition_provider_binding(resources.provider,
        resources.path_result)
    _require_exact_acquisition_products_target(
        resources.provider.products, target)
    validate_acquisition_provider_target(
        resources.provider.implementation,
        resources.path_result,
        resources.provider.products,
        target,
    )
    return resources
end

function _require_prepared_target_local_acquisition_resources(
    resources,
    ::AcquisitionDefinition,
    ::PreparedTargetLocalPathResources,
)
    throw(PlantPreparationError(
        :acquisition,
        :invalid_preparation,
        "target-local acquisition model preparation must return PreparedTargetLocalAcquisitionResources; got $(typeof(resources))",
    ))
end

"""
    prepare_target_local_acquisition_resources(definition, path_resources)

Prepare one acquisition's target-local provider resources from an already
prepared target-local path.  The model extension owns provider construction;
this wrapper enforces exact path and target co-location.
"""
function prepare_target_local_acquisition_resources(
    definition::AcquisitionDefinition,
    path_resources::PreparedTargetLocalPathResources,
)
    resources = _with_completed_prepared_device_execution_context(
        path_resources.context) do
        provider = prepare_target_local_acquisition_provider(
            acquisition_model(definition), definition, path_resources)
        prepared = _require_prepared_target_local_acquisition_provider(
            provider, definition, path_resources)
        PreparedTargetLocalAcquisitionResources(
            definition, path_resources, prepared)
    end
    return _require_prepared_target_local_acquisition_resources(
        resources, definition, path_resources)
end

function _require_prepared_target_local_acquisition_provider(
    provider::PreparedAcquisitionProvider,
    ::AcquisitionDefinition,
    path_resources::PreparedTargetLocalPathResources,
)
    validate_acquisition_provider_binding(provider, path_resources.result)
    return provider
end

function _require_prepared_target_local_acquisition_provider(
    provider,
    ::AcquisitionDefinition,
    ::PreparedTargetLocalPathResources,
)
    throw(PlantPreparationError(
        :acquisition,
        :invalid_preparation,
        "target-local acquisition model preparation must return PreparedAcquisitionProvider; got $(typeof(provider))",
    ))
end

"""Qualified cold extension seam for one target-local acquisition provider."""
function prepare_target_local_acquisition_provider(
    model,
    definition::AcquisitionDefinition,
    path_resources::PreparedTargetLocalPathResources,
)
    throw(PlantPreparationError(
        :acquisition,
        :unsupported_model,
        "acquisition model $(typeof(model)) does not implement prepare_target_local_acquisition_provider",
    ))
end
