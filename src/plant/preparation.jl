#
# Prepared plant ownership
#
# Preparation is intentionally model-extensible rather than a universal
# optical graph. Cold model definitions dispatch to concrete preparation
# methods; the resulting owners bind exact products and prepared stage plans.
#

@inline function _with_completed_prepared_device_execution_context(
    f::F, context) where {F}
    return _with_prepared_device_execution_context(context) do
        result = f()
        _synchronize_prepared_device_execution_context!(context)
        result
    end
end

"""Semantic contract for when an optical path result samples the plant."""
abstract type AbstractOpticalSamplingContract end

"""The path result is a photon-arrival-rate sample at one plant instant."""
struct InstantaneousOpticalSample <: AbstractOpticalSamplingContract end

struct _PathResultKeyToken end
const _PATH_RESULT_KEY_TOKEN = _PathResultKeyToken()

"""
    PathResultKey

Run-owned compatibility description for one prepared optical result. The key
records physical source geometry and radiometry, spectral sampling, optical
and propagation model identities, the sampling contract, output-plane
semantics, revisions, backend, and compute device. Acquisitions compare all
fields during preparation and retain the exact key and result binding.
"""
struct PathResultKey{G,S,R,M,C,P,O,V,B<:AbstractArrayBackend,
    D<:AbstractComputeDevice}
    source_geometry::G
    spectral_sampling::S
    radiometry::R
    optical_model::M
    sampling_contract::C
    propagation_model::P
    output_plane::O
    revisions::V
    backend::B
    device::D

    function PathResultKey(::_PathResultKeyToken, source_geometry,
        spectral_sampling, radiometry, optical_model, sampling_contract,
        propagation_model, output_plane, revisions,
        backend::B, device::D) where {
        B<:AbstractArrayBackend,D<:AbstractComputeDevice,
    }
        geometry_snapshot = deepcopy(source_geometry)
        spectral_snapshot = deepcopy(spectral_sampling)
        radiometry_snapshot = deepcopy(radiometry)
        optical_snapshot = deepcopy(optical_model)
        sampling_snapshot = deepcopy(sampling_contract)
        propagation_snapshot = deepcopy(propagation_model)
        output_snapshot = deepcopy(output_plane)
        revision_snapshot = deepcopy(revisions)
        return new{
            typeof(geometry_snapshot),typeof(spectral_snapshot),
            typeof(radiometry_snapshot),typeof(optical_snapshot),
            typeof(sampling_snapshot),typeof(propagation_snapshot),
            typeof(output_snapshot),typeof(revision_snapshot),B,D,
        }(geometry_snapshot, spectral_snapshot, radiometry_snapshot,
            optical_snapshot, sampling_snapshot, propagation_snapshot,
            output_snapshot, revision_snapshot, backend, device)
    end
end

function PathResultKey(source_geometry, spectral_sampling, radiometry,
    optical_model, sampling_contract, propagation_model, output_plane,
    revisions, backend::AbstractArrayBackend, device::AbstractComputeDevice)
    return PathResultKey(_PATH_RESULT_KEY_TOKEN, source_geometry,
        spectral_sampling, radiometry, optical_model, sampling_contract,
        propagation_model, output_plane, revisions, backend, device)
end

const _PATH_RESULT_KEY_SNAPSHOT_FIELDS = (
    :source_geometry,
    :spectral_sampling,
    :radiometry,
    :optical_model,
    :sampling_contract,
    :propagation_model,
    :output_plane,
    :revisions,
)

"""Return defensive copies for value descriptions retained by the key."""
function Base.getproperty(key::PathResultKey, name::Symbol)
    name in _PATH_RESULT_KEY_SNAPSHOT_FIELDS &&
        return deepcopy(getfield(key, name))
    return getfield(key, name)
end

function Base.isequal(left::PathResultKey, right::PathResultKey)
    return isequal(getfield(left, :source_geometry),
            getfield(right, :source_geometry)) &&
        isequal(getfield(left, :spectral_sampling),
            getfield(right, :spectral_sampling)) &&
        isequal(getfield(left, :radiometry), getfield(right, :radiometry)) &&
        isequal(getfield(left, :optical_model),
            getfield(right, :optical_model)) &&
        isequal(getfield(left, :sampling_contract),
            getfield(right, :sampling_contract)) &&
        isequal(getfield(left, :propagation_model),
            getfield(right, :propagation_model)) &&
        isequal(getfield(left, :output_plane),
            getfield(right, :output_plane)) &&
        isequal(getfield(left, :revisions), getfield(right, :revisions)) &&
        typeof(getfield(left, :backend)) ===
            typeof(getfield(right, :backend)) &&
        isequal(getfield(left, :device), getfield(right, :device))
end

Base.:(==)(left::PathResultKey, right::PathResultKey) =
    isequal(left, right)

function Base.hash(key::PathResultKey, seed::UInt)
    return hash((getfield(key, :source_geometry),
        getfield(key, :spectral_sampling), getfield(key, :radiometry),
        getfield(key, :optical_model), getfield(key, :sampling_contract),
        getfield(key, :propagation_model), getfield(key, :output_plane),
        getfield(key, :revisions), typeof(getfield(key, :backend)),
        getfield(key, :device)), seed)
end

@inline function _leaf_source_geometry_key(src::Source)
    return (
        kind=Source,
        direction_arcsec=coordinates_xy_arcsec(src),
        height_m=source_height_m(src),
    )
end

@inline function _leaf_source_geometry_key(src::LGSSource)
    params = src.params
    profile = isnothing(params.sodium_layer_profile) ? nothing :
        copy(params.sodium_layer_profile)
    return (
        kind=LGSSource,
        direction_arcsec=coordinates_xy_arcsec(src),
        height_m=source_height_m(src),
        laser_coordinates_m=params.laser_coordinates,
        elongation_factor=params.elongation_factor,
        uplink_fwhm_arcsec=params.fwhm_spot_up,
        sodium_layer_profile=profile,
    )
end

@inline path_source_geometry_key(src::Union{Source,LGSSource}) =
    _leaf_source_geometry_key(src)
@inline path_source_geometry_key(src::SpectralSource) =
    _leaf_source_geometry_key(src.source)

function path_source_geometry_key(src::Asterism)
    return map(path_source_geometry_key, src.sources)
end

function path_source_geometry_key(src::ExtendedSource)
    return map(path_source_geometry_key, src.quadrature.sources)
end

function path_source_geometry_key(src::AbstractSource)
    throw(PlantPreparationError(:path, :source_geometry,
        "source type $(typeof(src)) must implement a prepared path geometry key"))
end

@inline function _leaf_source_spectral_key(src::Union{Source,LGSSource})
    return ((wavelength_m=wavelength(src), weight=one(wavelength(src))),)
end

@inline path_source_spectral_key(src::Union{Source,LGSSource}) =
    _leaf_source_spectral_key(src)

function path_source_spectral_key(src::SpectralSource)
    return map(src.bundle.samples) do sample
        (wavelength_m=sample.wavelength, weight=sample.weight)
    end
end

function path_source_spectral_key(src::Asterism)
    return map(path_source_spectral_key, src.sources)
end

@inline path_source_spectral_key(src::ExtendedSource) =
    _leaf_source_spectral_key(src.source)

function path_source_spectral_key(src::AbstractSource)
    throw(PlantPreparationError(:path, :spectral_sampling,
        "source type $(typeof(src)) must implement a prepared spectral key"))
end

@inline function _leaf_source_radiometry_key(src::Union{Source,LGSSource})
    return (
        policy=typeof(source_radiometry(src)),
        value=source_radiometric_value(src),
    )
end

@inline path_source_radiometry_key(src::Union{Source,LGSSource}) =
    _leaf_source_radiometry_key(src)
@inline path_source_radiometry_key(src::SpectralSource) =
    _leaf_source_radiometry_key(src.source)

function path_source_radiometry_key(src::Asterism)
    return map(path_source_radiometry_key, src.sources)
end

function path_source_radiometry_key(src::ExtendedSource)
    return map(path_source_radiometry_key, src.quadrature.sources)
end

function path_source_radiometry_key(src::AbstractSource)
    throw(PlantPreparationError(:path, :radiometry,
        "source type $(typeof(src)) must implement a prepared radiometry key"))
end

@inline _require_path_result_plane(::FocalPlane) = nothing
@inline _require_path_result_plane(::DetectorPlane) = nothing

function _require_path_result_plane(kind::AbstractOpticalPlaneKind)
    throw(PlantPreparationError(:path, :output_plane,
        "acquisition-facing path results must be on a focal or detector plane; got $(typeof(kind))"))
end

@inline _require_path_result_rate(::PhotonRateNormalization) = nothing

function _require_path_result_rate(normalization::AbstractOpticalNormalization)
    throw(PlantPreparationError(:path, :radiometry,
        "acquisition-facing path results must use photon-rate normalization; got $(typeof(normalization))"))
end

@inline _require_path_result_measure(::SpatialDensityMeasure) = nothing
@inline _require_path_result_measure(::CellIntegratedMeasure) = nothing

function _require_path_result_measure(measure::AbstractSpatialMeasure)
    throw(PlantPreparationError(:path, :radiometry,
        "acquisition-facing path results must use spatial-density or cell-integrated samples; got $(typeof(measure))"))
end

@inline _require_path_result_coherence(::IncoherentIntensityAddition) =
    nothing

function _require_path_result_coherence(coherence::AbstractCombinationPolicy)
    throw(PlantPreparationError(:path, :radiometry,
        "acquisition-facing path results must declare incoherent intensity addition; got $(typeof(coherence))"))
end

@inline _require_path_result_spectral(::MonochromaticChannel) = nothing
@inline _require_path_result_spectral(::IntegratedSpectralChannel) = nothing

function _require_path_result_spectral(spectral::AbstractSpectralCoordinate)
    throw(PlantPreparationError(:path, :spectral_sampling,
        "acquisition-facing path results require a declared spectral channel; got $(typeof(spectral))"))
end

function _path_output_contract(product::IntensityMap)
    metadata = validate_plane_storage(product.metadata, product.values;
        label="prepared optical-path result")
    _require_path_result_plane(metadata.kind)
    _require_path_result_rate(metadata.normalization)
    _require_path_result_measure(metadata.spatial_measure)
    _require_path_result_coherence(metadata.coherence)
    _require_path_result_spectral(metadata.spectral)
    return (
        kind=metadata.kind,
        coordinate_domain=metadata.coordinate_domain,
        dimensions=metadata.dimensions,
        sampling=metadata.sampling,
        origin=metadata.origin,
        centering=metadata.centering,
        orientation=metadata.orientation,
        spectral=metadata.spectral,
        numeric_type=metadata.numeric_type,
        normalization=metadata.normalization,
        spatial_measure=metadata.spatial_measure,
        coherence=metadata.coherence,
    )
end

@inline _path_output_contract(bundle::OpticalProductBundle) =
    _path_output_contract(bundle.products)
@inline _path_output_contract(products::Tuple) =
    map(_path_output_contract, products)

function _path_output_contract(
    products::FixedSizeVector{<:AbstractOpticalProduct})
    isempty(products) && throw(PlantPreparationError(:path, :output_plane,
        "prepared optical-path result bundle must not be empty"))
    return map(_path_output_contract, products)
end

function _path_output_contract(result)
    throw(PlantPreparationError(:path, :output_plane,
        "prepared optical-path result must be an IntensityMap or concrete tuple/bundle of IntensityMap values; got $(typeof(result))"))
end

@inline _first_path_result(result::IntensityMap) = result
@inline _first_path_result(bundle::OpticalProductBundle) =
    _first_path_result(bundle.products)
@inline _first_path_result(products::Tuple{<:Any,Vararg}) =
    _first_path_result(first(products))

function _first_path_result(
    products::FixedSizeVector{<:AbstractOpticalProduct})
    isempty(products) && throw(PlantPreparationError(:path, :output_plane,
        "prepared optical-path result bundle must not be empty"))
    return _first_path_result(first(products))
end

function _first_path_result(::Tuple{})
    throw(PlantPreparationError(:path, :output_plane,
        "prepared optical-path result tuple must not be empty"))
end

@inline _require_path_result_domain(::Tuple{}, ::AbstractArrayBackend,
    ::AbstractComputeDevice) = nothing

@inline function _require_path_result_domain(product::IntensityMap,
    selector::AbstractArrayBackend, device::AbstractComputeDevice)
    typeof(backend(product)) === typeof(selector) || throw(
        PlantPreparationError(:path, :backend,
            "prepared path result leaves must use one array backend"))
    compute_device(product.values) == device || throw(PlantPreparationError(
        :path, :device,
        "prepared path result leaves must occupy one compute device"))
    return nothing
end

@inline _require_path_result_domain(bundle::OpticalProductBundle,
    selector::AbstractArrayBackend, device::AbstractComputeDevice) =
    _require_path_result_domain(bundle.products, selector, device)

@inline function _require_path_result_domain(products::Tuple,
    selector::AbstractArrayBackend, device::AbstractComputeDevice)
    _require_path_result_domain(first(products), selector, device)
    return _require_path_result_domain(Base.tail(products), selector, device)
end

function _require_path_result_domain(
    products::FixedSizeVector{<:AbstractOpticalProduct},
    selector::AbstractArrayBackend, device::AbstractComputeDevice)
    @inbounds for product in products
        _require_path_result_domain(product, selector, device)
    end
    return nothing
end

function _validate_path_input(input::PupilFunction)
    typeof(input.metadata.kind) === PupilPlane || throw(
        PlantPreparationError(:path, :input_plane,
            "prepared path PupilFunction input must be on a pupil plane"))
    validate_plane_storage(input.metadata, input.amplitude;
        label="prepared path pupil amplitude")
    validate_plane_storage(input.metadata, input.opd;
        label="prepared path pupil OPD")
    return input
end

function _validate_path_input(input::ElectricField)
    validate_plane_storage(input.metadata, input.values;
        label="prepared path electric field")
    return input
end

function _validate_path_input(input::IntensityMap)
    validate_plane_storage(input.metadata, input.values;
        label="prepared path intensity map")
    return input
end

@inline _validate_path_input(::Tuple{}) = throw(PlantPreparationError(
    :path, :input_plane, "prepared path input tuple must not be empty"))

@inline function _validate_path_input(inputs::Tuple)
    _validate_path_input(first(inputs))
    _validate_remaining_path_inputs(Base.tail(inputs))
    return inputs
end

@inline _validate_remaining_path_inputs(::Tuple{}) = nothing

@inline function _validate_remaining_path_inputs(inputs::Tuple)
    _validate_path_input(first(inputs))
    return _validate_remaining_path_inputs(Base.tail(inputs))
end

function _validate_path_input(input)
    throw(PlantPreparationError(:path, :input_plane,
        "prepared path input must be a PupilFunction, ElectricField, IntensityMap, or concrete tuple of them; got $(typeof(input))"))
end

@inline _path_input_storage(input::PupilFunction) = input.opd
@inline _path_input_storage(input::ElectricField) = input.values
@inline _path_input_storage(input::IntensityMap) = input.values

@inline function _require_path_input_domain(
    input::Union{PupilFunction,ElectricField,IntensityMap},
    selector::AbstractArrayBackend,
    device::AbstractComputeDevice,
)
    typeof(backend(input)) === typeof(selector) || throw(
        PlantPreparationError(:path, :backend,
            "prepared path input and result backends differ"))
    compute_device(_path_input_storage(input)) == device || throw(
        PlantPreparationError(:path, :device,
            "prepared path input and result occupy different compute devices"))
    return nothing
end

@inline _require_path_input_domain(::Tuple{}, ::AbstractArrayBackend,
    ::AbstractComputeDevice) = nothing

@inline function _require_path_input_domain(inputs::Tuple,
    selector::AbstractArrayBackend, device::AbstractComputeDevice)
    _require_path_input_domain(first(inputs), selector, device)
    return _require_path_input_domain(Base.tail(inputs), selector, device)
end

@inline _require_path_input_revision(::ElectricField, ::UInt) = nothing
@inline _require_path_input_revision(::IntensityMap, ::UInt) = nothing

function _require_path_input_revision(input::PupilFunction, revision::UInt)
    aperture_revision(input) == revision || throw(PlantPreparationError(
        :path, :revision,
        "prepared path pupil aperture revision does not match its telescope"))
    return nothing
end

@inline _require_path_input_revisions(::Tuple{}, ::UInt) = nothing

@inline function _require_path_input_revisions(inputs::Tuple, revision::UInt)
    _require_path_input_revision(first(inputs), revision)
    return _require_path_input_revisions(Base.tail(inputs), revision)
end

@inline _require_path_input_revisions(input, revision::UInt) =
    _require_path_input_revision(input, revision)

"""Qualified extension seam for cold path execution-binding validation."""
function validate_path_execution_binding(execution, input, result)
    throw(PlantPreparationError(:path, :unsupported_binding_validation,
        "prepared path execution type $(typeof(execution)) does not validate its input/result binding"))
end

"""
    AtmosphereIndependentPath()

Explicit materialization policy for a path whose prepared input does not read
the plant atmosphere. This is an assertion by the optical-model extension, not
a fallback for an unsupported atmosphere-dependent input.
"""
struct AtmosphereIndependentPath end

"""
    PreparedPupilOPDMaterialization

Prepared binding from one timed atmosphere direction to the exact path-local
`PupilFunction` whose OPD is materialized before optical-path execution.
Construct it with `prepare_pupil_opd_materialization`.
"""
mutable struct _PreparedPupilOPDMaterializationToken end
const _PREPARED_PUPIL_OPD_MATERIALIZATION_TOKEN =
    _PreparedPupilOPDMaterializationToken()

struct PreparedPupilOPDMaterialization{R,P,S,T}
    renderer::R
    destination::P
    source::S
    telescope::T

    function PreparedPupilOPDMaterialization(
        token::_PreparedPupilOPDMaterializationToken,
        renderer::R,
        destination::P,
        source::S,
        telescope::T,
    ) where {R,P,S,T}
        token === _PREPARED_PUPIL_OPD_MATERIALIZATION_TOKEN || throw(
            ArgumentError("invalid internal pupil-OPD materialization token"))
        return new{R,P,S,T}(renderer, destination, source, telescope)
    end
end

"""
    prepare_pupil_opd_materialization(atmosphere, telescope, source, pupil)

Prepare current-epoch atmospheric OPD materialization for one frozen path
source and exact caller-owned pupil input.
"""
function prepare_pupil_opd_materialization(
    atmosphere::AbstractTimedAtmosphere,
    telescope::Telescope,
    source::AbstractSource,
    pupil::PupilFunction,
)
    renderer = prepare_atmosphere_renderer(atmosphere, telescope, source;
        T=eltype(pupil.opd))
    _validate_atmosphere_renderer_binding(renderer, atmosphere)
    _validate_atmosphere_destination(pupil, renderer)
    return PreparedPupilOPDMaterialization(
        _PREPARED_PUPIL_OPD_MATERIALIZATION_TOKEN,
        renderer,
        pupil,
        source,
        telescope,
    )
end

"""Qualified extension seam for prepared materialization-binding validation."""
function validate_path_materialization_binding(materialization, input,
    atmosphere, source)
    throw(PlantPreparationError(:path,
        :unsupported_materialization_binding,
        "prepared path materialization type $(typeof(materialization)) does not validate its atmosphere/input binding"))
end

@inline validate_path_materialization_binding(
    ::AtmosphereIndependentPath, input, atmosphere, source) = nothing

@inline _validate_path_materialization_telescope(
    materialization, telescope) = nothing

function _validate_path_materialization_telescope(
    materialization::PreparedPupilOPDMaterialization,
    telescope::AbstractTelescope,
)
    materialization.telescope === telescope || throw(PlantPreparationError(
        :path,
        :prepared_binding,
        "pupil-OPD materialization does not retain the exact telescope",
    ))
    return nothing
end

function _validate_pupil_opd_renderer_source_binding(
    renderer::AtmosphereDirectionRenderer,
    source::AbstractSource,
)
    T = renderer.output_metadata.numeric_type
    isequal(
        source_geometry_signature(renderer.source, T),
        source_geometry_signature(source, T),
    ) || throw(PlantPreparationError(
        :path,
        :prepared_binding,
        "pupil-OPD renderer does not retain the requested source geometry",
    ))
    return nothing
end

function _validate_pupil_opd_renderer_source_binding(
    renderer,
    ::AbstractSource,
)
    throw(PlantPreparationError(
        :path,
        :unsupported_materialization_binding,
        "pupil-OPD materialization renderer $(typeof(renderer)) does not " *
        "validate its source geometry",
    ))
end

function validate_path_materialization_binding(
    materialization::PreparedPupilOPDMaterialization,
    input::PupilFunction,
    atmosphere::AbstractTimedAtmosphere,
    source::AbstractSource,
)
    materialization.destination === input || throw(PlantPreparationError(
        :path, :prepared_binding,
        "pupil-OPD materialization does not retain the exact path input"))
    materialization.source === source || throw(PlantPreparationError(
        :path, :prepared_binding,
        "pupil-OPD materialization does not retain the exact frozen source"))
    aperture_revision(input) == aperture_revision(materialization.telescope) ||
        throw(PlantPreparationError(
            :path,
            :revision,
            "pupil-OPD materialization telescope revision does not match its destination",
        ))
    _validate_pupil_opd_renderer_source_binding(
        materialization.renderer, source)
    _validate_atmosphere_renderer_binding(materialization.renderer,
        atmosphere)
    _validate_atmosphere_destination(input, materialization.renderer)
    return nothing
end

"""Qualified extension seam for mutation-free current-epoch preflight."""
function validate_path_materialization(materialization, input, atmosphere,
    epoch)
    throw(PlantPreparationError(:path, :unsupported_materialization,
        "prepared path materialization type $(typeof(materialization)) does not validate current-epoch materialization"))
end

@inline validate_path_materialization(
    ::AtmosphereIndependentPath, input, atmosphere, epoch) = input

function validate_path_materialization(
    materialization::PreparedPupilOPDMaterialization,
    input::PupilFunction,
    atmosphere::AbstractTimedAtmosphere,
    epoch::AtmosphereEpoch,
)
    validate_path_materialization_binding(materialization, input,
        atmosphere, materialization.source)
    validate_atmosphere_rendering(input, materialization.renderer,
        atmosphere, epoch)
    return input
end

"""Qualified extension seam for materializing one validated path input."""
function materialize_path_input!(materialization, input, atmosphere, epoch)
    throw(PlantPreparationError(:path, :unsupported_materialization,
        "prepared path materialization type $(typeof(materialization)) does not implement materialize_path_input!"))
end

@inline materialize_path_input!(
    ::AtmosphereIndependentPath, input, atmosphere, epoch) = input

@inline function materialize_path_input!(materialization, input,
    atmosphere, epoch, ::AbstractRNG)
    return materialize_path_input!(materialization, input, atmosphere,
        epoch)
end

"""Qualified materialization seam with one exact prepared path RNG group."""
@inline function materialize_path_input_rngs!(materialization, input,
    atmosphere, epoch, ::PreparedOwnerRNGs)
    return materialize_path_input!(materialization, input, atmosphere,
        epoch)
end

function materialize_path_input!(
    materialization::PreparedPupilOPDMaterialization,
    input::PupilFunction,
    atmosphere::AbstractTimedAtmosphere,
    epoch::AtmosphereEpoch,
)
    return render_atmosphere!(input, materialization.renderer, atmosphere,
        epoch)
end

struct _PreparedPathExecutorToken end
const _PREPARED_PATH_EXECUTOR_TOKEN = _PreparedPathExecutorToken()

"""
    PreparedPathExecutor

Concrete single-writer owner for one prepared optical path. `input` and
`result` are explicit path-local products, `materialization` binds the path's
atmosphere-dependent input operation, and `execution` owns or references the
concrete prepared propagation/front-end workspace used to update `result`.
An input may begin at a declared pupil-function, electric-field, or intensity
entry boundary; it need not represent the entrance pupil.
"""
struct PreparedPathExecutor{S,T,A,X,I,R,M,E,K<:PathResultKey}
    definition_slot::UInt32
    id::OpticalPathID
    source::S
    telescope::T
    atmosphere::A
    context::X
    rng_token::RNGOwnerToken
    input::I
    result::R
    materialization::M
    execution::E
    key::K

    function PreparedPathExecutor(::_PreparedPathExecutorToken,
        definition_slot::UInt32, id::OpticalPathID, source::S,
        telescope::T, atmosphere::A, context::X, rng_token::RNGOwnerToken,
        input::I, result::R, materialization::M, execution::E,
        key::K) where {
        S,T,A,X,I,R,M,E,K<:PathResultKey,
    }
        return new{S,T,A,X,I,R,M,E,K}(definition_slot, id, source,
            telescope, atmosphere, context, rng_token, input, result,
            materialization, execution, key)
    end
end

@inline _rng_owner_binding_token(path::PreparedPathExecutor) =
    path.rng_token
@inline path_id(path::PreparedPathExecutor) = path.id

function PreparedPathExecutor(definition::OpticalPathDefinition,
    source::AbstractSource, telescope::AbstractTelescope,
    atmosphere::AbstractAtmosphere, input, result, execution;
    context,
    materialization,
    optical_model,
    sampling_contract::AbstractOpticalSamplingContract=
        InstantaneousOpticalSample(),
    propagation_model,
    model_revisions=())
    _validate_path_input(input)
    output_plane = _path_output_contract(result)
    first_result = _first_path_result(result)
    selector = backend(first_result)
    device = compute_device(first_result.values)
    _prepared_device_execution_compute_device(context) == device || throw(
        PlantPreparationError(:path, :execution_context_target,
            "prepared path execution context does not match its exact compute device"))
    _require_path_result_domain(result, selector, device)
    _require_path_input_domain(input, selector, device)
    typeof(backend(telescope)) === typeof(selector) || throw(
        PlantPreparationError(:path, :backend,
            "prepared path and telescope backends differ"))
    compute_device(pupil_reflectivity(telescope)) == device || throw(
        PlantPreparationError(:path, :device,
            "prepared path and telescope occupy different compute devices"))
    revision = aperture_revision(telescope)
    _require_path_input_revisions(input, revision)
    validate_path_materialization_binding(materialization, input,
        atmosphere, source)
    _validate_path_materialization_telescope(materialization, telescope)
    validate_path_execution_binding(execution, input, result)
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
    return PreparedPathExecutor(_PREPARED_PATH_EXECUTOR_TOKEN, UInt32(0),
        path_id(definition), source, telescope, atmosphere, context,
        RNGOwnerToken(), input, result, materialization, execution, key)
end

@inline path_input(path::PreparedPathExecutor) = path.input
@inline path_result(path::PreparedPathExecutor) = path.result
@inline path_result_key(path::PreparedPathExecutor) = path.key
@inline path_materialization(path::PreparedPathExecutor) =
    path.materialization

function _require_current_path_revision(path::PreparedPathExecutor)
    revision = getfield(getfield(path.key, :revisions), :telescope)
    aperture_revision(path.telescope) == revision || throw(
        PlantPreparationError(:path, :revision,
            "telescope aperture changed after path preparation"))
    _require_path_input_revisions(path.input, revision)
    return nothing
end

function _require_current_path_binding(path::PreparedPathExecutor)
    _require_current_path_revision(path)
    validate_path_materialization_binding(path.materialization, path.input,
        path.atmosphere, path.source)
    _validate_path_materialization_telescope(
        path.materialization, path.telescope)
    validate_path_execution_binding(path.execution, path.input, path.result)
    return nothing
end

@inline _require_timed_path_atmosphere(
    atmosphere::AbstractTimedAtmosphere) = atmosphere

function _require_timed_path_atmosphere(atmosphere::AbstractAtmosphere)
    throw(PlantPreparationError(:path, :untimed_atmosphere,
        "current-epoch path materialization requires a timed atmosphere"))
end

function _materialize_path_input_in_context!(path::PreparedPathExecutor,
    epoch::AtmosphereEpoch)
    _require_current_path_binding(path)
    atmosphere = _require_timed_path_atmosphere(path.atmosphere)
    _validate_epoch_identity(atmosphere_identity(atmosphere), atmosphere,
        epoch)
    validate_path_materialization(path.materialization, path.input,
        atmosphere, epoch)
    return materialize_path_input!(path.materialization, path.input,
        atmosphere, epoch)
end

"""
    materialize_path_input!(path, epoch)

Materialize one path's atmosphere-dependent input from the explicitly supplied
current epoch. This operation never advances the atmosphere.
"""
function materialize_path_input!(path::PreparedPathExecutor,
    epoch::AtmosphereEpoch)
    result = _with_completed_prepared_device_execution_context(
        path.context) do
        _materialize_path_input_in_context!(path, epoch)
    end
    return result::typeof(path.input)
end

function _materialize_path_input_in_context!(path::PreparedPathExecutor,
    epoch::AtmosphereEpoch, rng::AbstractRNG)
    _require_current_path_binding(path)
    atmosphere = _require_timed_path_atmosphere(path.atmosphere)
    _validate_epoch_identity(atmosphere_identity(atmosphere), atmosphere,
        epoch)
    validate_path_materialization(path.materialization, path.input,
        atmosphere, epoch)
    return materialize_path_input!(path.materialization, path.input,
        atmosphere, epoch, rng)
end

"""Materialize one path input at an explicit epoch with a caller-owned RNG."""
function materialize_path_input!(path::PreparedPathExecutor,
    epoch::AtmosphereEpoch, rng::AbstractRNG)
    result = _with_completed_prepared_device_execution_context(
        path.context) do
        _materialize_path_input_in_context!(path, epoch, rng)
    end
    return result::typeof(path.input)
end

function _execute_path_in_context!(path::PreparedPathExecutor)
    _require_current_path_binding(path)
    return execute_path!(path.result, path.input, path.execution)
end

"""Execute one already prepared path without selecting or advancing time."""
function execute_path!(path::PreparedPathExecutor)
    result = _with_completed_prepared_device_execution_context(
        path.context) do
        _execute_path_in_context!(path)
    end
    return result::typeof(path.result)
end

function _execute_path_in_context!(path::PreparedPathExecutor,
    rngs::PreparedOwnerRNGs)
    _require_current_path_binding(path)
    _require_rng_owner_binding(rngs, path)
    return execute_path_rngs!(path.result, path.input, path.execution,
        rngs)
end

"""Execute one prepared path with its exact prepared RNG-owner group."""
function execute_path!(path::PreparedPathExecutor,
    rngs::PreparedOwnerRNGs)
    result = _with_completed_prepared_device_execution_context(
        path.context) do
        _execute_path_in_context!(path, rngs)
    end
    return result::typeof(path.result)
end

"""Qualified extension seam for path execution with prepared RNG owners."""
function execute_path_rngs!(result, input, execution,
    rngs::PreparedOwnerRNGs)
    return execute_path!(result, input, execution,
        rng_stream_state(rngs, Val(:provider)))
end

@inline execute_path!(result, input, execution, ::AbstractRNG) =
    execute_path!(result, input, execution)

function execute_path!(result, input, execution)
    throw(PlantPreparationError(:path, :unsupported_execution,
        "prepared path execution type $(typeof(execution)) does not implement execute_path!"))
end

@inline function _require_direct_path_input(
    prepared::PreparedDirectImaging, input)
    prepared.input === input || throw(PlantPreparationError(:path,
        :prepared_binding,
        "direct-imaging input does not match its prepared path"))
    return nothing
end

function _require_direct_path_input(prepared::Union{
    PreparedIncoherentDirectImaging,PreparedBundledDirectImaging}, input)
    @inbounds for component in prepared.components
        component.input === input || throw(PlantPreparationError(:path,
            :prepared_binding,
            "direct-imaging component input does not match its prepared path"))
    end
    return nothing
end

function execute_path!(result, input, prepared::Union{
    PreparedDirectImaging,PreparedIncoherentDirectImaging,
    PreparedBundledDirectImaging})
    validate_path_execution_binding(prepared, input, result)
    return form_direct_image!(prepared)
end

function validate_path_execution_binding(prepared::Union{
    PreparedDirectImaging,PreparedIncoherentDirectImaging,
    PreparedBundledDirectImaging}, input, result)
    direct_imaging_output(prepared) === result || throw(
        PlantPreparationError(:path, :prepared_binding,
            "direct-imaging output does not match its prepared path"))
    _require_direct_path_input(prepared, input)
    return nothing
end

"""Adapter from a Gate 0 prepared WFS optical plan to a plant path."""
struct WFSOpticalPathExecution{P}
    plan::P
end

@inline function validate_path_execution_binding(
    execution::WFSOpticalPathExecution, input, result)
    return validate_wfs_optics_binding(result, input,
        execution.plan)
end

@inline function execute_path!(result, input,
    execution::WFSOpticalPathExecution)
    return form_wfs_optical_products!(result, input, execution.plan)
end

struct _FrameAcquisitionExecutionToken end
const _FRAME_ACQUISITION_EXECUTION_TOKEN = _FrameAcquisitionExecutionToken()

"""Prepared detector capture and copy into a distinct caller-owned frame."""
struct FrameAcquisitionExecution{P,F}
    acquisition::P
    observation::F

    function FrameAcquisitionExecution(::_FrameAcquisitionExecutionToken,
        acquisition::P, observation::F) where {P,F}
        return new{P,F}(acquisition, observation)
    end
end

function FrameAcquisitionExecution(detector::Detector,
    optical_result::IntensityMap, observation::AbstractArray)
    _require_frame_acquisition_observation_ownership(detector,
        optical_result, observation)
    candidate = _prepare_detached_detector_acquisition(detector,
        optical_result)
    return _commit_frame_acquisition_execution!(detector, candidate,
        observation)
end

function FrameAcquisitionExecution(detector::Detector,
    optical_result::IntensityMap)
    candidate = _prepare_detached_detector_acquisition(detector,
        optical_result)
    observation = similar(output_frame(
        detector_acquisition_detector(candidate)))
    fill!(observation, zero(eltype(observation)))
    return _commit_frame_acquisition_execution!(detector, candidate,
        observation)
end

function _commit_frame_acquisition_execution!(detector::Detector,
    candidate::PreparedDetectorAcquisition,
    observation::AbstractArray)
    _require_frame_acquisition_observation(
        detector_acquisition_detector(candidate),
        detector_acquisition_input(candidate), observation)
    acquisition = _rebind_prepared_detector_acquisition(detector, candidate)
    execution = FrameAcquisitionExecution(_FRAME_ACQUISITION_EXECUTION_TOKEN,
        acquisition, observation)
    _commit_prepared_detector_acquisition!(acquisition)
    return execution
end

function _require_frame_acquisition_observation_ownership(
    detector::Detector, optical_result::IntensityMap,
    observation::AbstractArray)
    (_detector_storage_mightalias(detector, observation) ||
        Base.mightalias(observation, optical_result.values)) && throw(
        PlantPreparationError(:acquisition, :ownership,
            "caller-owned acquisition observation must not alias detector storage or its optical input"))
    return nothing
end

function _require_frame_acquisition_observation(detector::Detector,
    optical_result::IntensityMap, observation::AbstractArray)
    frame = output_frame(detector)
    size(observation) == size(frame) || throw(PlantPreparationError(
        :acquisition, :shape,
        "acquisition observation shape must match detector output"))
    eltype(observation) === eltype(frame) || throw(PlantPreparationError(
        :acquisition, :numeric_type,
        "acquisition observation element type must match detector output"))
    typeof(backend(observation)) === typeof(backend(frame)) || throw(
        PlantPreparationError(:acquisition, :backend,
            "acquisition observation and detector output backends differ"))
    compute_device(observation) == compute_device(frame) || throw(
        PlantPreparationError(:acquisition, :device,
            "acquisition observation and detector output occupy different devices"))
    _require_frame_acquisition_observation_ownership(detector,
        optical_result, observation)
    return nothing
end

struct _WFSAcquisitionExecutionToken end
const _WFS_ACQUISITION_EXECUTION_TOKEN = _WFSAcquisitionExecutionToken()

"""Prepared Gate 0 WFS acquisition and optional estimator composition."""
struct WFSAcquisitionExecution{A,E,O,M}
    acquisition::A
    estimator::E
    observation::O
    measurement::M

    function WFSAcquisitionExecution(::_WFSAcquisitionExecutionToken,
        acquisition::A, estimator::E, observation::O,
        measurement::M) where {A,E,O,M}
        return new{A,E,O,M}(acquisition, estimator, observation, measurement)
    end
end

@inline WFSAcquisitionExecution(acquisition, observation) =
    WFSAcquisitionExecution(_WFS_ACQUISITION_EXECUTION_TOKEN, acquisition,
        nothing, observation, nothing)

function WFSAcquisitionExecution(acquisition, estimator,
    observation, measurement::WFSMeasurement)
    return WFSAcquisitionExecution(_WFS_ACQUISITION_EXECUTION_TOKEN,
        acquisition, estimator, observation, measurement)
end

"""Qualified extension seam for cold acquisition execution validation."""
function validate_acquisition_execution_binding(execution, path_result,
    products::AcquisitionProducts)
    throw(PlantPreparationError(:acquisition,
        :unsupported_binding_validation,
        "prepared acquisition execution type $(typeof(execution)) does not validate its path/product binding"))
end

function validate_acquisition_execution_binding(
    execution::FrameAcquisitionExecution,
    path_result::IntensityMap,
    products::AcquisitionProducts{<:AbstractArray,Nothing})
    products.observation === execution.observation || throw(
        PlantPreparationError(:acquisition, :prepared_binding,
            "acquisition observation does not match its prepared storage"))
    detector = detector_acquisition_detector(execution.acquisition)
    input = detector_acquisition_input(execution.acquisition)
    _require_frame_acquisition_observation(detector, input,
        execution.observation)
    path_result.metadata === input.metadata &&
        path_result.values === input.values || throw(
        PlantPreparationError(:acquisition, :prepared_binding,
            "acquisition path result does not match its prepared detector input"))
    try
        _require_prepared_acquisition(execution.acquisition)
    catch error
        error isa Union{InvalidConfiguration,DimensionMismatchError} ||
            rethrow()
        throw(PlantPreparationError(:acquisition, :prepared_binding,
            "acquisition detector ownership changed after preparation"))
    end
    return nothing
end

@inline _require_acquired_wfs_estimator(::AcquiredObservationPath) = nothing

function _require_acquired_wfs_estimator(::DirectMeasurementPath)
    throw(PlantPreparationError(:acquisition, :estimator,
        "WFS acquisition composition requires an acquired-observation estimator"))
end

function _require_acquired_wfs_estimator(path)
    throw(PlantPreparationError(:acquisition, :estimator,
        "WFS estimator declared unsupported measurement path $(typeof(path))"))
end

function validate_acquisition_execution_binding(
    execution::WFSAcquisitionExecution{A,Nothing,O,Nothing}, path_result,
    products::AcquisitionProducts{O,Nothing}) where {A,O}
    products.observation === execution.observation || throw(
        PlantPreparationError(:acquisition, :prepared_binding,
            "WFS observation does not match its prepared acquisition"))
    validate_wfs_acquisition_binding(products.observation, path_result,
        execution.acquisition)
    return nothing
end

function validate_acquisition_execution_binding(
    execution::WFSAcquisitionExecution{A,E,O,M}, path_result,
    products::AcquisitionProducts{O,M}) where {A,E,O,M<:WFSMeasurement}
    products.observation === execution.observation &&
        products.measurement === execution.measurement || throw(
        PlantPreparationError(:acquisition, :prepared_binding,
            "WFS products do not match their prepared acquisition"))
    validate_wfs_acquisition_binding(products.observation, path_result,
        execution.acquisition)
    _require_acquired_wfs_estimator(wfs_measurement_path(execution.estimator))
    validate_wfs_estimation_binding(products.measurement,
        products.observation, execution.estimator)
    return nothing
end

function _validate_acquisition_provider_binding(
    provider::PreparedFullOpticalProvider, path_result,
    products::AcquisitionProducts, contract::AcquisitionProductContract)
    return validate_acquisition_execution_binding(provider.execution,
        path_result, products)
end

function execute_acquisition_provider!(products::AcquisitionProducts,
    path_result, provider::PreparedFullOpticalProvider,
    rngs::PreparedOwnerRNGs)
    return execute_acquisition_rngs!(products, path_result,
        provider.execution, rngs)
end

function execute_acquisition_provider!(products::AcquisitionProducts,
    path_result, provider::PreparedFullOpticalProvider,
    rng::AbstractRNG)
    return execute_acquisition!(products, path_result,
        provider.execution, rng)
end

struct _PreparedAcquisitionOwnerToken end
const _PREPARED_ACQUISITION_OWNER_TOKEN = _PreparedAcquisitionOwnerToken()

"""
    PreparedAcquisitionOwner

Concrete single-writer owner for one run-immutable acquisition provider. It
borrows one exact prepared path result as read-only input and owns separate
caller-visible observation/measurement products plus the provider's detector,
WFS, reduced-order, payload, or replay state.
"""
struct PreparedAcquisitionOwner{K<:PathResultKey,R,X,P}
    definition_slot::UInt32
    id::AcquisitionID
    path_id::OpticalPathID
    path_key::K
    path_result::R
    context::X
    rng_token::RNGOwnerToken
    provider::P

    function PreparedAcquisitionOwner(::_PreparedAcquisitionOwnerToken,
        definition_slot::UInt32, id::AcquisitionID,
        path_id::OpticalPathID,
        path_key::K, path_result::R, context::X,
        rng_token::RNGOwnerToken, provider::P) where {
        K<:PathResultKey,R,X,P,
    }
        return new{K,R,X,P}(definition_slot, id, path_id, path_key,
            path_result, context, rng_token, provider)
    end
end

@inline _rng_owner_binding_token(owner::PreparedAcquisitionOwner) =
    owner.rng_token
@inline acquisition_id(owner::PreparedAcquisitionOwner) = owner.id
@inline acquisition_path_id(owner::PreparedAcquisitionOwner) = owner.path_id

@inline acquisition_provider(owner::PreparedAcquisitionOwner) = owner.provider
@inline acquisition_provider_style(owner::PreparedAcquisitionOwner) =
    acquisition_provider_style(owner.provider)
@inline acquisition_provider_payload_work(owner::PreparedAcquisitionOwner) =
    acquisition_provider_payload_work(owner.provider)
@inline acquisition_product_contract(owner::PreparedAcquisitionOwner) =
    acquisition_product_contract(owner.provider)
@inline acquisition_products(owner::PreparedAcquisitionOwner) =
    owner.provider.products
@inline acquisition_observation(owner::PreparedAcquisitionOwner) =
    owner.provider.products.observation
@inline acquisition_measurement(owner::PreparedAcquisitionOwner) =
    owner.provider.products.measurement
@inline acquisition_product_metadata(owner::PreparedAcquisitionOwner) =
    owner.provider.products.metadata

function _path_key_mismatch_reason(actual::PathResultKey,
    expected::PathResultKey)
    !isequal(getfield(actual, :source_geometry),
        getfield(expected, :source_geometry)) &&
        return :source_geometry
    !isequal(getfield(actual, :spectral_sampling),
        getfield(expected, :spectral_sampling)) &&
        return :spectral_sampling
    !isequal(getfield(actual, :radiometry),
        getfield(expected, :radiometry)) && return :radiometry
    !isequal(getfield(actual, :optical_model),
        getfield(expected, :optical_model)) &&
        return :optical_model
    !isequal(getfield(actual, :sampling_contract),
        getfield(expected, :sampling_contract)) &&
        return :sampling_contract
    !isequal(getfield(actual, :propagation_model),
        getfield(expected, :propagation_model)) &&
        return :propagation_model
    !isequal(getfield(actual, :output_plane),
        getfield(expected, :output_plane)) &&
        return :output_plane
    !isequal(getfield(actual, :revisions), getfield(expected, :revisions)) &&
        return :revision
    typeof(getfield(actual, :backend)) ===
        typeof(getfield(expected, :backend)) || return :backend
    !isequal(getfield(actual, :device), getfield(expected, :device)) &&
        return :device
    return :compatibility
end

function _require_path_result_key(actual::PathResultKey,
    expected::PathResultKey)
    isequal(actual, expected) && return actual
    reason = _path_key_mismatch_reason(actual, expected)
    throw(PlantPreparationError(:acquisition, reason,
        "acquisition requires a path result incompatible in $(reason)"))
end

"""
    require_path_result(path; ...)

Validate an acquisition's cold compatibility requirements against a prepared
path before constructing or mutating acquisition destinations. Omitted fields
retain the prepared path's value. This is a preparation-time extension seam,
not a hot-path lookup.
"""
function require_path_result(path::PreparedPathExecutor;
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

function PreparedAcquisitionOwner(definition::AcquisitionDefinition,
    path::PreparedPathExecutor, provider::PreparedAcquisitionProvider)
    acquisition_path_id(definition) == path_id(path) || throw(
        PlantPreparationError(:acquisition, :unknown_path,
            "acquisition $(definition.id) does not reference prepared path $(path_id(path))"))
    validate_acquisition_provider_binding(provider, path.result)
    return PreparedAcquisitionOwner(_PREPARED_ACQUISITION_OWNER_TOKEN,
        UInt32(0), acquisition_id(definition),
        acquisition_path_id(definition), path.key, path.result,
        path.context, RNGOwnerToken(), provider)
end

@inline function _execute_acquisition_in_context!(
    owner::PreparedAcquisitionOwner,
    rng::AbstractRNG)
    return execute_acquisition_provider!(owner.provider, owner.path_result,
        rng)
end

"""Execute one acquisition provider into its caller-owned products."""
function execute_acquisition!(owner::PreparedAcquisitionOwner,
    rng::AbstractRNG)
    result = _with_completed_prepared_device_execution_context(
        owner.context) do
        _execute_acquisition_in_context!(owner, rng)
    end
    return result::typeof(acquisition_products(owner))
end

function _execute_acquisition_in_context!(owner::PreparedAcquisitionOwner,
    rngs::PreparedOwnerRNGs)
    _require_rng_owner_binding(rngs, owner)
    return execute_acquisition_provider!(owner.provider, owner.path_result,
        rngs)
end

"""Execute one acquisition using its exact prepared RNG-owner group."""
function execute_acquisition!(owner::PreparedAcquisitionOwner,
    rngs::PreparedOwnerRNGs)
    result = _with_completed_prepared_device_execution_context(
        owner.context) do
        _execute_acquisition_in_context!(owner, rngs)
    end
    return result::typeof(acquisition_products(owner))
end

"""Qualified extension seam for acquisition execution with prepared RNG owners."""
function execute_acquisition_rngs!(products, path_result, execution,
    rngs::PreparedOwnerRNGs)
    return execute_acquisition!(products, path_result, execution,
        rng_stream_state(rngs, Val(:detector)))
end

function execute_acquisition!(products, path_result, execution, rng)
    throw(PlantPreparationError(:acquisition, :unsupported_execution,
        "prepared acquisition execution type $(typeof(execution)) does not implement execute_acquisition!"))
end

function execute_acquisition!(products::AcquisitionProducts{<:AbstractArray,
    Nothing}, path_result::IntensityMap,
    execution::FrameAcquisitionExecution, rng::AbstractRNG)
    validate_acquisition_execution_binding(execution, path_result, products)
    frame = capture!(execution.acquisition, rng)
    copyto!(products.observation, frame)
    return products
end

function execute_acquisition!(products::AcquisitionProducts{O,Nothing},
    path_result,
    execution::WFSAcquisitionExecution{A,Nothing,O,Nothing}, rng) where {A,O}
    validate_acquisition_execution_binding(execution, path_result, products)
    acquire_wfs_observation!(products.observation, path_result,
        execution.acquisition, rng)
    return products
end

function execute_acquisition!(products::AcquisitionProducts{O,M},
    path_result,
    execution::WFSAcquisitionExecution{A,E,O,M}, rng) where {
    A,E,O,M<:WFSMeasurement,
}
    validate_acquisition_execution_binding(execution, path_result, products)
    acquire_wfs_observation!(products.observation, path_result,
        execution.acquisition, rng)
    estimate_wfs_measurement!(products.measurement, products.observation,
        execution.estimator)
    return products
end

"""Compact definition-order location in exact path-family storage."""
struct _PreparedPathSlot
    family_slot::UInt32
    member_slot::UInt32
end

"""Fixed-capacity prepared paths with one exact executable owner type."""
struct _PreparedPathFamily{
    O<:PreparedPathExecutor,
    V<:FixedSizeVector{O},
}
    values::V
end

"""
Definition-order view over exact prepared path execution families.

The family tuple type depends only on the distinct normalized executable
owner types. Cold exact-type collection and identity keys make family-group
order process-canonical and independent of definition order;
owner cardinality and definition order remain value-level data in `slots`,
with each owner retaining its stable identity.
"""
struct _PreparedPathRegistry{
    D<:_FixedPlantRegistry{OpticalPathDefinition},
    G<:Tuple,
    S<:FixedSizeVector{_PreparedPathSlot},
}
    definitions::D
    groups::G
    slots::S
end

Base.size(registry::_PreparedPathRegistry) = (length(registry.slots),)
Base.axes(registry::_PreparedPathRegistry) = axes(registry.slots)
Base.length(registry::_PreparedPathRegistry) = length(registry.slots)
Base.keys(registry::_PreparedPathRegistry) = eachindex(registry.slots)
Base.eachindex(registry::_PreparedPathRegistry) = eachindex(registry.slots)
Base.firstindex(registry::_PreparedPathRegistry) =
    firstindex(registry.slots)
Base.lastindex(registry::_PreparedPathRegistry) = lastindex(registry.slots)

@noinline function _prepared_path_slot_error(family::Int, member::Int)
    throw(BoundsError((family_slot=family, member_slot=member)))
end

@inline function _prepared_path_family_value(
    ::Tuple{}, family::Int, member::Int)
    return _prepared_path_slot_error(family, member)
end

@inline function _prepared_path_family_value(
    groups::Tuple, family::Int, member::Int)
    family == 1 && return @inbounds groups[1].values[member]
    return _prepared_path_family_value(
        Base.tail(groups), family - 1, member)
end

@inline function _prepared_path_value(
    registry::_PreparedPathRegistry, slot::_PreparedPathSlot)
    return _prepared_path_family_value(registry.groups,
        Int(slot.family_slot), Int(slot.member_slot))
end

@inline function Base.getindex(registry::_PreparedPathRegistry, index::Int)
    checkbounds(registry.slots, index)
    return _prepared_path_value(registry, @inbounds(registry.slots[index]))
end

@inline function Base.iterate(registry::_PreparedPathRegistry, state::Int=1)
    state > length(registry) && return nothing
    return (@inbounds registry[state], state + 1)
end

@inline function _prepared_path_definition(
    registry::_PreparedPathRegistry, path::PreparedPathExecutor)
    slot = Int(path.definition_slot)
    checkbounds(registry.definitions, slot)
    definition = @inbounds registry.definitions[slot]
    path_id(definition) == path_id(path) || throw(PlantPreparationError(
        :path, :prepared_binding,
        "prepared path definition slot does not match its stable identity"))
    return definition
end

struct _PreparedPathFamilyType{O<:PreparedPathExecutor} end

@inline _prepared_path_family_matches(
    ::_PreparedPathFamilyType{O}, ::Type{O}) where {
    O<:PreparedPathExecutor,
} = true

@inline _prepared_path_family_matches(
    ::_PreparedPathFamilyType, ::Type{<:PreparedPathExecutor}) = false

function _prepared_path_owner_types(paths)
    Base.@nospecialize paths
    owner_types = DataType[]
    @inbounds for path in paths
        owner_type = typeof(path)::DataType
        present = false
        for existing in owner_types
            if existing === owner_type
                present = true
                break
            end
        end
        present && continue
        push!(owner_types, owner_type)
    end
    sort!(owner_types; by=_prepared_owner_family_order_key)
    return _require_unique_prepared_path_family_order_keys(owner_types)
end

@inline _prepared_owner_family_order_key(owner_type::Type) =
    (hash(owner_type, UInt(0)), objectid(owner_type))

function _require_unique_prepared_path_family_order_keys(owner_types)
    length(owner_types) <= 1 && return owner_types
    previous_type = owner_types[1]
    previous_key = _prepared_owner_family_order_key(previous_type)
    @inbounds for index in 2:length(owner_types)
        owner_type = owner_types[index]
        key = _prepared_owner_family_order_key(owner_type)
        key == previous_key && owner_type !== previous_type && throw(
            PlantPreparationError(:path, :family_type_key_collision,
                "distinct prepared path owner types have the same process-local family key"))
        previous_type = owner_type
        previous_key = key
    end
    return owner_types
end

function _prepared_path_family_types(paths)
    owner_types = _prepared_path_owner_types(paths)
    families = ()
    for owner_type in owner_types
        family_type = Core.apply_type(_PreparedPathFamilyType, owner_type)
        families = (families..., family_type())
    end
    return families
end

function _prepare_path_family(
    ::_PreparedPathFamilyType{O}, paths) where {O<:PreparedPathExecutor}
    Base.@nospecialize paths
    count = 0
    @inbounds for path in paths
        typeof(path) === O && (count += 1)
    end
    values = Vector{O}(undef, count)
    next = 1
    @inbounds for path in paths
        typeof(path) === O || continue
        values[next] = path
        next += 1
    end
    fixed = FixedSizeVectorDefault{O}(values)
    return _PreparedPathFamily{O,typeof(fixed)}(fixed)
end

@inline _prepare_path_families(::Tuple{}, paths) = ()

function _prepare_path_families(families::Tuple, paths)
    Base.@nospecialize paths
    first = _prepare_path_family(families[1], paths)
    rest = _prepare_path_families(Base.tail(families), paths)
    return (first, rest...)
end

@inline function _prepared_path_family_index(
    ::Tuple{}, ::Type{<:PreparedPathExecutor}, ::Int=1)
    return 0
end

@inline function _prepared_path_family_index(
    families::Tuple, ::Type{O}, index::Int=1,
) where {O<:PreparedPathExecutor}
    _prepared_path_family_matches(families[1], O) && return index
    return _prepared_path_family_index(
        Base.tail(families), O, index + 1)
end

function _prepare_path_slots(families::Tuple, paths)
    Base.@nospecialize paths
    length(families) <= typemax(UInt32) || throw(PlantPreparationError(
        :path, :family_capacity,
        "prepared path family count exceeds UInt32 capacity"))
    counts = zeros(UInt32, length(families))
    slots = Vector{_PreparedPathSlot}(undef, length(paths))
    @inbounds for index in eachindex(paths)
        family = _prepared_path_family_index(families, typeof(paths[index]))
        iszero(family) && throw(PlantPreparationError(
            :path, :missing_family,
            "prepared path has no exact executable family"))
        counts[family] == typemax(UInt32) && throw(PlantPreparationError(
            :path, :family_capacity,
            "prepared path family exceeds UInt32 capacity"))
        counts[family] += UInt32(1)
        slots[index] = _PreparedPathSlot(UInt32(family), counts[family])
    end
    return FixedSizeVectorDefault{_PreparedPathSlot}(slots)
end

function _require_prepared_path_slots(groups::Tuple, slots)
    @inbounds for slot in slots
        family = Int(slot.family_slot)
        1 <= family <= length(groups) ||
            _prepared_path_slot_error(family, Int(slot.member_slot))
        checkbounds(groups[family].values, Int(slot.member_slot))
    end
    return slots
end

function _prepare_path_registry(definitions, paths)
    Base.@nospecialize paths
    family_types = _prepared_path_family_types(paths)
    groups = _prepare_path_families(family_types, paths)
    slots = _prepare_path_slots(family_types, paths)
    _require_prepared_path_slots(groups, slots)
    return _PreparedPathRegistry(definitions, groups, slots)
end

"""Compact definition-order location in exact acquisition-family storage."""
struct _PreparedAcquisitionSlot
    family_slot::UInt32
    member_slot::UInt32
end

"""Fixed-capacity acquisitions with one exact executable owner type."""
struct _PreparedAcquisitionFamily{
    O<:PreparedAcquisitionOwner,
    V<:FixedSizeVector{O},
}
    values::V
end

"""
Definition-order view over exact prepared acquisition families.

Cold exact-type collection and identity keys make family-group order
process-canonical and independent of definition order. Definition-order
traversal is retained by `slots`, and each owner retains its stable identity.
"""
struct _PreparedAcquisitionRegistry{
    D<:_FixedPlantRegistry{AcquisitionDefinition},
    G<:Tuple,
    S<:FixedSizeVector{_PreparedAcquisitionSlot},
}
    definitions::D
    groups::G
    slots::S
end

Base.size(registry::_PreparedAcquisitionRegistry) =
    (length(registry.slots),)
Base.axes(registry::_PreparedAcquisitionRegistry) = axes(registry.slots)
Base.length(registry::_PreparedAcquisitionRegistry) = length(registry.slots)
Base.keys(registry::_PreparedAcquisitionRegistry) = eachindex(registry.slots)
Base.eachindex(registry::_PreparedAcquisitionRegistry) =
    eachindex(registry.slots)
Base.firstindex(registry::_PreparedAcquisitionRegistry) =
    firstindex(registry.slots)
Base.lastindex(registry::_PreparedAcquisitionRegistry) =
    lastindex(registry.slots)

@noinline function _prepared_acquisition_slot_error(
    family::Int, member::Int)
    throw(BoundsError((family_slot=family, member_slot=member)))
end

@inline function _prepared_acquisition_family_value(
    ::Tuple{}, family::Int, member::Int)
    return _prepared_acquisition_slot_error(family, member)
end

@inline function _prepared_acquisition_family_value(
    groups::Tuple, family::Int, member::Int)
    family == 1 && return @inbounds groups[1].values[member]
    return _prepared_acquisition_family_value(
        Base.tail(groups), family - 1, member)
end

@inline function _prepared_acquisition_value(
    registry::_PreparedAcquisitionRegistry, slot::_PreparedAcquisitionSlot)
    return _prepared_acquisition_family_value(registry.groups,
        Int(slot.family_slot), Int(slot.member_slot))
end

@inline function Base.getindex(
    registry::_PreparedAcquisitionRegistry, index::Int)
    checkbounds(registry.slots, index)
    return _prepared_acquisition_value(
        registry, @inbounds(registry.slots[index]))
end


@inline function Base.iterate(
    registry::_PreparedAcquisitionRegistry, state::Int=1)
    state > length(registry) && return nothing
    return (@inbounds registry[state], state + 1)
end

@inline function _prepared_acquisition_definition(
    registry::_PreparedAcquisitionRegistry,
    owner::PreparedAcquisitionOwner)
    slot = Int(owner.definition_slot)
    checkbounds(registry.definitions, slot)
    definition = @inbounds registry.definitions[slot]
    acquisition_id(definition) == acquisition_id(owner) || throw(
        PlantPreparationError(:acquisition, :prepared_binding,
            "prepared acquisition definition slot does not match its stable identity"))
    acquisition_path_id(definition) == acquisition_path_id(owner) || throw(
        PlantPreparationError(:acquisition, :prepared_binding,
            "prepared acquisition definition slot does not match its path identity"))
    return definition
end

struct _PreparedAcquisitionFamilyType{O<:PreparedAcquisitionOwner} end

@inline _prepared_acquisition_family_matches(
    ::_PreparedAcquisitionFamilyType{O}, ::Type{O}) where {
    O<:PreparedAcquisitionOwner,
} = true

@inline _prepared_acquisition_family_matches(
    ::_PreparedAcquisitionFamilyType,
    ::Type{<:PreparedAcquisitionOwner}) = false

function _prepared_acquisition_owner_types(acquisitions)
    Base.@nospecialize acquisitions
    owner_types = DataType[]
    @inbounds for owner in acquisitions
        owner_type = typeof(owner)::DataType
        present = false
        for existing in owner_types
            if existing === owner_type
                present = true
                break
            end
        end
        present && continue
        push!(owner_types, owner_type)
    end
    sort!(owner_types; by=_prepared_owner_family_order_key)
    return _require_unique_prepared_acquisition_family_order_keys(owner_types)
end

function _require_unique_prepared_acquisition_family_order_keys(owner_types)
    length(owner_types) <= 1 && return owner_types
    previous_type = owner_types[1]
    previous_key = _prepared_owner_family_order_key(previous_type)
    @inbounds for index in 2:length(owner_types)
        owner_type = owner_types[index]
        key = _prepared_owner_family_order_key(owner_type)
        key == previous_key && owner_type !== previous_type && throw(
            PlantPreparationError(:acquisition,
                :family_type_key_collision,
                "distinct prepared acquisition owner types have the same process-local family key"))
        previous_type = owner_type
        previous_key = key
    end
    return owner_types
end

function _prepared_acquisition_family_types(acquisitions)
    owner_types = _prepared_acquisition_owner_types(acquisitions)
    families = ()
    for owner_type in owner_types
        family_type =
            Core.apply_type(_PreparedAcquisitionFamilyType, owner_type)
        families = (families..., family_type())
    end
    return families
end

function _prepare_acquisition_family(
    ::_PreparedAcquisitionFamilyType{O}, acquisitions,
) where {O<:PreparedAcquisitionOwner}
    Base.@nospecialize acquisitions
    count = 0
    @inbounds for owner in acquisitions
        typeof(owner) === O && (count += 1)
    end
    values = Vector{O}(undef, count)
    next = 1
    @inbounds for owner in acquisitions
        typeof(owner) === O || continue
        values[next] = owner
        next += 1
    end
    fixed = FixedSizeVectorDefault{O}(values)
    return _PreparedAcquisitionFamily{O,typeof(fixed)}(fixed)
end

@inline _prepare_acquisition_families(::Tuple{}, acquisitions) = ()

function _prepare_acquisition_families(families::Tuple, acquisitions)
    Base.@nospecialize acquisitions
    first = _prepare_acquisition_family(families[1], acquisitions)
    rest = _prepare_acquisition_families(
        Base.tail(families), acquisitions)
    return (first, rest...)
end

@inline function _prepared_acquisition_family_index(
    ::Tuple{}, ::Type{<:PreparedAcquisitionOwner}, ::Int=1)
    return 0
end

@inline function _prepared_acquisition_family_index(
    families::Tuple, ::Type{O}, index::Int=1,
) where {O<:PreparedAcquisitionOwner}
    _prepared_acquisition_family_matches(families[1], O) && return index
    return _prepared_acquisition_family_index(
        Base.tail(families), O, index + 1)
end

function _prepare_acquisition_slots(families::Tuple, acquisitions)
    Base.@nospecialize acquisitions
    length(families) <= typemax(UInt32) || throw(PlantPreparationError(
        :acquisition, :family_capacity,
        "prepared acquisition family count exceeds UInt32 capacity"))
    counts = zeros(UInt32, length(families))
    slots = Vector{_PreparedAcquisitionSlot}(undef, length(acquisitions))
    @inbounds for index in eachindex(acquisitions)
        family = _prepared_acquisition_family_index(
            families, typeof(acquisitions[index]))
        iszero(family) && throw(PlantPreparationError(
            :acquisition, :missing_family,
            "prepared acquisition has no exact executable family"))
        counts[family] == typemax(UInt32) && throw(PlantPreparationError(
            :acquisition, :family_capacity,
            "prepared acquisition family exceeds UInt32 capacity"))
        counts[family] += UInt32(1)
        slots[index] = _PreparedAcquisitionSlot(
            UInt32(family), counts[family])
    end
    return FixedSizeVectorDefault{_PreparedAcquisitionSlot}(slots)
end

function _require_prepared_acquisition_slots(groups::Tuple, slots)
    @inbounds for slot in slots
        family = Int(slot.family_slot)
        1 <= family <= length(groups) ||
            _prepared_acquisition_slot_error(
                family, Int(slot.member_slot))
        checkbounds(groups[family].values, Int(slot.member_slot))
    end
    return slots
end

function _prepare_acquisition_registry(definitions, acquisitions)
    Base.@nospecialize acquisitions
    family_types =
        _prepared_acquisition_family_types(acquisitions)
    groups = _prepare_acquisition_families(family_types, acquisitions)
    slots = _prepare_acquisition_slots(family_types, acquisitions)
    _require_prepared_acquisition_slots(groups, slots)
    return _PreparedAcquisitionRegistry(definitions, groups, slots)
end

struct _PreparedPlantToken end
const _PREPARED_PLANT_TOKEN = _PreparedPlantToken()

"""
Prepared, schedule-free plant with concrete owners in fixed-size homogeneous
registries, bounded path-to-optic bindings, and RNG streams.
"""
struct PreparedPlant{
    D<:PlantDefinition,
    C<:AbstractComputeDevice,
    X,
    T<:AbstractTelescope,
    A<:AbstractTimedAtmosphere,
    O<:_PreparedControllableOpticRegistry,
    S,
    B,
    P<:_PreparedPathRegistry,
    Q<:_PreparedAcquisitionRegistry,
    R<:PreparedPlantRNGs,
}
    definition::D
    target::C
    context::X
    telescope::T
    atmosphere::A
    controllable_optics::O
    controllable_optic_path_bindings::PreparedControllableOpticPathBindings
    sampled_aberrations::S
    sampled_aberration_path_bindings::B
    command_endpoints::Memory{_PreparedPlantCommandEndpoint}
    paths::P
    acquisitions::Q
    rngs::R

    function PreparedPlant(::_PreparedPlantToken, definition::D,
        target::C,
        context::X,
        telescope::T,
        atmosphere::A,
        controllable_optics::O,
        bindings::PreparedControllableOpticPathBindings,
        sampled_aberrations::S,
        sampled_bindings::B,
        command_endpoints::Memory{_PreparedPlantCommandEndpoint},
        paths::P,
        acquisitions::Q,
        rngs::R,
    ) where {
        D<:PlantDefinition,
        C<:AbstractComputeDevice,
        X,
        T<:AbstractTelescope,
        A<:AbstractTimedAtmosphere,
        O<:_PreparedControllableOpticRegistry,
        S,
        B,
        P<:_PreparedPathRegistry,
        Q<:_PreparedAcquisitionRegistry,
        R<:PreparedPlantRNGs,
    }
        _prepared_device_execution_compute_device(context) == target ||
            throw(PlantPreparationError(:plant, :execution_context_target,
                "prepared execution context does not match the exact plant target"))
        validate_telescope_target(telescope, target)
        validate_timed_atmosphere_target(atmosphere, target)
        plant = new{D,C,X,T,A,O,S,B,P,Q,R}(
            definition,
            target,
            context,
            telescope,
            atmosphere,
            controllable_optics,
            bindings,
            sampled_aberrations,
            sampled_bindings,
            command_endpoints,
            paths,
            acquisitions,
            rngs,
        )
        return _require_exact_prepared_plant_target(plant, target)
    end
end

@inline plant_definition(plant::PreparedPlant) = plant.definition
@inline compute_device(plant::PreparedPlant) = plant.target
@inline prepared_telescope(plant::PreparedPlant) = plant.telescope
@inline prepared_atmosphere(plant::PreparedPlant) = plant.atmosphere
@inline prepared_controllable_optics(plant::PreparedPlant) =
    plant.controllable_optics
@inline prepared_controllable_optic_path_bindings(plant::PreparedPlant) =
    plant.controllable_optic_path_bindings
@inline prepared_sampled_aberrations(plant::PreparedPlant) =
    plant.sampled_aberrations
@inline prepared_sampled_aberration_path_bindings(plant::PreparedPlant) =
    plant.sampled_aberration_path_bindings
@inline prepared_command_endpoints(plant::PreparedPlant) =
    map(binding -> binding.endpoint, plant.command_endpoints)
@inline prepared_paths(plant::PreparedPlant) = plant.paths
@inline prepared_acquisitions(plant::PreparedPlant) = plant.acquisitions

"""Return structured seed, derivation, owner, and derived-stream metadata."""
@inline rng_replay_metadata(plant::PreparedPlant) =
    _rng_replay_metadata(plant.rngs)

function prepared_path(plant::PreparedPlant, id)
    resolved = _as_optical_path_id(id)
    for path in plant.paths
        path_id(path) == resolved && return path
    end
    throw(PlantPreparationError(:path, :unknown_id,
        "prepared plant has no optical path $resolved"))
end

function prepared_acquisition(plant::PreparedPlant, id)
    resolved = _as_acquisition_id(id)
    for acquisition in plant.acquisitions
        acquisition_id(acquisition) == resolved &&
            return acquisition
    end
    throw(PlantPreparationError(:acquisition, :unknown_id,
        "prepared plant has no acquisition $resolved"))
end

function prepared_controllable_optic(plant::PreparedPlant, id)
    resolved = _as_controllable_optic_id(id)
    optics = plant.controllable_optics
    for index in eachindex(optics)
        definition = _prepared_controllable_optic_definition(optics, index)
        controllable_optic_id(definition) == resolved && return optics[index]
    end
    throw(PlantPreparationError(:controllable_optic, :unknown_id,
        "prepared plant has no controllable optic $resolved"))
end

function prepared_sampled_aberration(plant::PreparedPlant, id)
    resolved = _as_sampled_aberration_id(id)
    for aberration in plant.sampled_aberrations
        sampled_aberration_id(aberration) == resolved && return aberration
    end
    throw(PlantPreparationError(:sampled_aberration, :unknown_id,
        "prepared plant has no sampled aberration $resolved"))
end

function prepared_command_endpoint(plant::PreparedPlant, id)
    resolved = _as_command_endpoint_id(id)
    for binding in plant.command_endpoints
        command_endpoint_id(binding) == resolved && return binding.endpoint
    end
    throw(PlantPreparationError(:command_endpoint, :unknown_id,
        "prepared plant has no command endpoint $resolved"))
end

function prepare_path_executor(definition::OpticalPathDefinition,
    telescope::AbstractTelescope, atmosphere::AbstractAtmosphere)
    context = _prepare_device_execution_context(
        pupil_reflectivity(telescope))
    return _with_completed_prepared_device_execution_context(context) do
        prepare_path_executor(
            definition, telescope, atmosphere, context)
    end
end

function prepare_path_executor(definition::OpticalPathDefinition,
    telescope::AbstractTelescope, atmosphere::AbstractAtmosphere,
    context)
    source = freeze_source(path_source(definition))
    prepared = prepare_path_executor(path_model(definition), definition,
        source, telescope, atmosphere, context)
    return _require_prepared_path_executor(prepared, definition, source,
        telescope, atmosphere, context)
end

function _require_prepared_path_executor(prepared::PreparedPathExecutor,
    definition::OpticalPathDefinition, source::AbstractSource,
    telescope::AbstractTelescope, atmosphere::AbstractAtmosphere,
    context)
    path_id(prepared) == path_id(definition) || throw(PlantPreparationError(
        :path, :prepared_binding,
        "prepared path does not retain its declared stable identity"))
    prepared.source === source || throw(PlantPreparationError(:path,
        :prepared_binding,
        "prepared path does not retain its run-owned frozen source"))
    prepared.telescope === telescope || throw(PlantPreparationError(:path,
        :prepared_binding,
        "prepared path does not retain its plant telescope"))
    prepared.atmosphere === atmosphere || throw(PlantPreparationError(:path,
        :prepared_binding,
        "prepared path does not retain its plant atmosphere"))
    prepared.context === context || throw(PlantPreparationError(:path,
        :prepared_binding,
        "prepared path does not retain its plant execution context"))
    validate_path_materialization_binding(prepared.materialization,
        prepared.input, prepared.atmosphere, prepared.source)
    _validate_path_materialization_telescope(
        prepared.materialization, prepared.telescope)
    validate_path_execution_binding(prepared.execution, prepared.input,
        prepared.result)
    return prepared
end

function _seal_prepared_path_executor(
    prepared::PreparedPathExecutor,
    definition::OpticalPathDefinition,
    definition_slot::UInt32,
)
    iszero(definition_slot) && throw(PlantPreparationError(
        :path, :definition_slot,
        "prepared path definition slot must be positive"))
    iszero(prepared.definition_slot) || throw(PlantPreparationError(
        :path, :definition_slot,
        "path model preparation returned an already sealed owner"))
    path_id(prepared) == path_id(definition) || throw(
        PlantPreparationError(:path, :prepared_binding,
            "prepared path stable identity changed before sealing"))
    return PreparedPathExecutor(
        _PREPARED_PATH_EXECUTOR_TOKEN,
        definition_slot,
        path_id(prepared),
        prepared.source,
        prepared.telescope,
        prepared.atmosphere,
        prepared.context,
        prepared.rng_token,
        prepared.input,
        prepared.result,
        prepared.materialization,
        prepared.execution,
        prepared.key,
    )
end


function _require_prepared_path_executor(prepared,
    ::OpticalPathDefinition, ::AbstractSource, ::AbstractTelescope,
    ::AbstractAtmosphere, context)
    throw(PlantPreparationError(
        :path, :invalid_preparation,
        "path model preparation must return PreparedPathExecutor; got $(typeof(prepared))"))
end

function prepare_path_executor(model, definition::OpticalPathDefinition,
    source::AbstractSource, telescope::AbstractTelescope,
    atmosphere::AbstractAtmosphere, context)
    throw(PlantPreparationError(:path, :unsupported_model,
        "path model $(typeof(model)) does not implement prepare_path_executor"))
end

function _prepare_acquisition_owner_in_context(
    definition::AcquisitionDefinition,
    path::PreparedPathExecutor)
    provider = prepare_acquisition_provider(acquisition_model(definition),
        definition, path)
    prepared = _require_prepared_acquisition_provider(provider, definition,
        path)
    owner = PreparedAcquisitionOwner(definition, path, prepared)
    return _require_prepared_acquisition_owner(owner, definition, path)
end

function prepare_acquisition_owner(definition::AcquisitionDefinition,
    path::PreparedPathExecutor)
    result = _with_completed_prepared_device_execution_context(
        path.context) do
        _prepare_acquisition_owner_in_context(definition, path)
    end
    return result::PreparedAcquisitionOwner
end

function _require_prepared_acquisition_owner(
    prepared::PreparedAcquisitionOwner,
    definition::AcquisitionDefinition, path::PreparedPathExecutor)
    acquisition_id(prepared) == acquisition_id(definition) || throw(
        PlantPreparationError(
        :acquisition, :prepared_binding,
        "prepared acquisition does not retain its declared stable identity"))
    acquisition_path_id(prepared) == acquisition_path_id(definition) ||
        throw(PlantPreparationError(:acquisition, :prepared_binding,
            "prepared acquisition does not retain its declared path identity"))
    prepared.path_key === path.key || throw(PlantPreparationError(
        :acquisition, :prepared_binding,
        "prepared acquisition does not retain the exact path-result key"))
    prepared.path_result === path.result || throw(PlantPreparationError(
        :acquisition, :prepared_binding,
        "prepared acquisition does not retain the exact path result"))
    prepared.context === path.context || throw(PlantPreparationError(
        :acquisition, :prepared_binding,
        "prepared acquisition does not retain the exact path execution context"))
    validate_acquisition_provider_binding(prepared.provider,
        prepared.path_result)
    return prepared
end

function _seal_prepared_acquisition_owner(
    prepared::PreparedAcquisitionOwner,
    definition::AcquisitionDefinition,
    definition_slot::UInt32,
)
    iszero(definition_slot) && throw(PlantPreparationError(
        :acquisition, :definition_slot,
        "prepared acquisition definition slot must be positive"))
    iszero(prepared.definition_slot) || throw(PlantPreparationError(
        :acquisition, :definition_slot,
        "acquisition preparation returned an already sealed owner"))
    acquisition_id(prepared) == acquisition_id(definition) || throw(
        PlantPreparationError(:acquisition, :prepared_binding,
            "prepared acquisition stable identity changed before sealing"))
    acquisition_path_id(prepared) == acquisition_path_id(definition) ||
        throw(PlantPreparationError(:acquisition, :prepared_binding,
            "prepared acquisition path identity changed before sealing"))
    return PreparedAcquisitionOwner(
        _PREPARED_ACQUISITION_OWNER_TOKEN,
        definition_slot,
        acquisition_id(prepared),
        acquisition_path_id(prepared),
        prepared.path_key,
        prepared.path_result,
        prepared.context,
        prepared.rng_token,
        prepared.provider,
    )
end

function _require_prepared_acquisition_owner(prepared,
    ::AcquisitionDefinition, ::PreparedPathExecutor)
    throw(PlantPreparationError(
        :acquisition, :invalid_preparation,
        "acquisition owner construction must return PreparedAcquisitionOwner; got $(typeof(prepared))"))
end

function _require_prepared_acquisition_provider(
    prepared::PreparedAcquisitionProvider,
    ::AcquisitionDefinition, path::PreparedPathExecutor)
    validate_acquisition_provider_binding(prepared, path.result)
    return prepared
end

function _require_prepared_acquisition_provider(prepared,
    ::AcquisitionDefinition, ::PreparedPathExecutor)
    throw(PlantPreparationError(:acquisition, :invalid_preparation,
        "acquisition model preparation must return PreparedAcquisitionProvider; got $(typeof(prepared))"))
end

"""Qualified cold preparation seam for one acquisition product provider."""
function prepare_acquisition_provider(model,
    definition::AcquisitionDefinition, path::PreparedPathExecutor)
    throw(PlantPreparationError(:acquisition, :unsupported_model,
        "acquisition model $(typeof(model)) does not implement prepare_acquisition_provider"))
end

function _prepare_path_executors(definitions::AbstractVector, telescope,
    atmosphere, context)
    length(definitions) <= typemax(UInt32) || throw(PlantPreparationError(
        :path, :capacity,
        "prepared path count exceeds UInt32 capacity"))
    paths = PreparedPathExecutor[]
    sizehint!(paths, length(definitions))
    @inbounds for index in eachindex(definitions)
        definition = definitions[index]
        prepared = prepare_path_executor(
            definition, telescope, atmosphere, context)
        push!(paths, _seal_prepared_path_executor(
            prepared, definition, UInt32(index)))
    end
    return _prepare_path_registry(definitions, paths)
end

@inline function _prepared_path_for_acquisition(definition,
    paths::_PreparedPathRegistry)
    id = acquisition_path_id(definition)
    for path in paths
        path_id(path) == id && return path
    end
    throw(PlantPreparationError(:acquisition, :unknown_path,
        "acquisition $(definition.id) references an unprepared path $id"))
end

function _prepare_acquisition_owners(definitions::AbstractVector, paths)
    length(definitions) <= typemax(UInt32) || throw(PlantPreparationError(
        :acquisition, :capacity,
        "prepared acquisition count exceeds UInt32 capacity"))
    acquisitions = PreparedAcquisitionOwner[]
    sizehint!(acquisitions, length(definitions))
    @inbounds for index in eachindex(definitions)
        definition = definitions[index]
        path = _prepared_path_for_acquisition(definition, paths)
        prepared = _prepare_acquisition_owner_in_context(definition, path)
        push!(acquisitions, _seal_prepared_acquisition_owner(
            prepared, definition, UInt32(index)))
    end
    return _prepare_acquisition_registry(definitions, acquisitions)
end

function _path_rng_owner_bindings(paths::_PreparedPathRegistry)
    groups = Vector{Tuple}(undef, length(paths))
    @inbounds for index in eachindex(paths)
        path = paths[index]
        groups[index] = _rng_owner_bindings(path, :path,
            path_id(path).name,
            _path_rng_owner_roles(path.execution, path.materialization))
    end
    return groups
end

function _acquisition_rng_owner_bindings(
    acquisitions::_PreparedAcquisitionRegistry)
    groups = Vector{Tuple}(undef, length(acquisitions))
    @inbounds for index in eachindex(acquisitions)
        acquisition = acquisitions[index]
        groups[index] = _rng_owner_bindings(acquisition, :acquisition,
            acquisition_id(acquisition).name,
            _acquisition_rng_owner_roles(
                acquisition.provider.implementation))
    end
    return groups
end

function _append_rng_binding_groups!(destination, groups)
    for group in groups
        append!(destination, group)
    end
    return destination
end

function _prepare_owner_rng_family(owner_group, bindings,
    run_seed::UInt64, version::RNGDerivationVersion)
    owners = owner_group.values
    isempty(owners) && throw(PlantPreparationError(:rng, :owner_topology,
        "prepared execution families must not be empty"))
    first_owner = first(owners)
    first_rngs = _prepare_owner_rngs(
        bindings[Int(first_owner.definition_slot)], run_seed, version)
    R = typeof(first_rngs)
    values = Vector{R}(undef, length(owners))
    values[1] = first_rngs
    @inbounds for index in 2:length(owners)
        owner = owners[index]
        prepared = _prepare_owner_rngs(
            bindings[Int(owner.definition_slot)], run_seed, version)
        typeof(prepared) === R || throw(PlantPreparationError(
            :rng, :owner_roles,
            "one exact execution family declared inconsistent RNG roles"))
        values[index] = prepared
    end
    fixed = FixedSizeVectorDefault{R}(values)
    return _PreparedOwnerRNGFamily{R,typeof(fixed)}(fixed)
end

@inline _prepare_owner_rng_families(::Tuple{}, bindings,
    run_seed::UInt64, version::RNGDerivationVersion) = ()

function _prepare_owner_rng_families(owner_groups::Tuple, bindings,
    run_seed::UInt64, version::RNGDerivationVersion)
    first = _prepare_owner_rng_family(
        owner_groups[1], bindings, run_seed, version)
    rest = _prepare_owner_rng_families(
        Base.tail(owner_groups), bindings, run_seed, version)
    return (first, rest...)
end

function _prepare_owner_rng_registry(owner_registry, bindings,
    run_seed::UInt64, version::RNGDerivationVersion)
    length(owner_registry) == length(bindings) || throw(
        PlantPreparationError(:rng, :owner_topology,
            "prepared RNG bindings do not match execution-owner cardinality"))
    groups = _prepare_owner_rng_families(
        owner_registry.groups, bindings, run_seed, version)
    length(groups) == length(owner_registry.groups) || throw(
        PlantPreparationError(:rng, :owner_topology,
            "prepared RNG families do not match execution-owner families"))
    return _PreparedOwnerRNGRegistry(groups, owner_registry.slots)
end

function _prepare_plant_rngs(atmosphere::AbstractTimedAtmosphere,
    paths::_PreparedPathRegistry,
    acquisitions::_PreparedAcquisitionRegistry, run_seed::UInt64,
    version::RNGDerivationVersion)
    atmosphere_bindings = _atmosphere_rng_owner_bindings(atmosphere)
    path_bindings = _path_rng_owner_bindings(paths)
    acquisition_bindings = _acquisition_rng_owner_bindings(acquisitions)
    all_bindings = RNGOwnerBinding[]
    sizehint!(all_bindings, length(atmosphere_bindings) +
        sum(length, path_bindings; init=0) +
        sum(length, acquisition_bindings; init=0))
    append!(all_bindings, atmosphere_bindings)
    _append_rng_binding_groups!(all_bindings, path_bindings)
    _append_rng_binding_groups!(all_bindings, acquisition_bindings)
    _require_unique_rng_owner_identities(all_bindings)
    _require_unique_rng_stream_seeds(all_bindings, run_seed, version)
    atmosphere_rngs = _prepare_atmosphere_rngs(atmosphere,
        atmosphere_bindings, run_seed, version)
    path_rngs = _prepare_owner_rng_registry(
        paths, path_bindings, run_seed, version)
    acquisition_rngs = _prepare_owner_rng_registry(
        acquisitions, acquisition_bindings, run_seed, version)
    return PreparedPlantRNGs(run_seed, version, atmosphere_rngs,
        path_rngs, acquisition_rngs)
end

"""
    prepare_plant(definition, target; run_seed, rng_derivation_version,
        command_endpoints=())

Prepare all declared controllable optics, command endpoints, paths, and
acquisitions without scheduling or executing them. Model-specific construction
dispatches on cold model-definition types. Preparation may allocate and perform
fallible backend/device/revision validation. It also derives exact stateful RNG
streams from the required run seed, derivation version, and stable owner
identities. Whole-plant owners are retained in fixed-size homogeneous
registries; dispatch into each owner recovers its concrete small pipeline.
Each declared command endpoint requires one separate
`CommandEndpointConfiguration`; endpoint ordinals and optic order derive from
stable identities rather than declaration position. Required optic placement
and visibility declarations resolve to canonical bounded per-path bindings and
co-placed groups.
"""
function prepare_plant(definition::PlantDefinition,
    target::AbstractComputeDevice;
    run_seed,
    rng_derivation_version=_DEFAULT_RNG_DERIVATION_VERSION,
    command_endpoints=())
    seed = _prepare_run_seed(run_seed)
    version = _prepare_rng_derivation_version(rng_derivation_version)
    context = _prepare_device_execution_context(target)
    _prepared_device_execution_compute_device(context) == target ||
        throw(PlantPreparationError(:plant, :execution_context_target,
            "prepared execution context does not match the requested exact target"))
    return _with_completed_prepared_device_execution_context(context) do
        telescope = validate_telescope_target(
            prepare_telescope(telescope_definition(definition), target),
            target,
        )
        atmosphere = validate_timed_atmosphere_target(
            prepare_timed_atmosphere(
                atmosphere_definition(definition), telescope, target),
            target,
        )
        endpoint_configurations =
            _sorted_command_endpoint_configurations(command_endpoints)
        optic_definitions =
            _canonical_controllable_optic_definitions(definition)
        prepared_endpoints = _prepare_plant_command_endpoints(definition,
            endpoint_configurations, optic_definitions, target)
        optics = _prepare_controllable_optics(definition, optic_definitions,
            prepared_endpoints, telescope, atmosphere)
        sampled_aberrations =
            _prepare_sampled_aberrations(definition, target)
        paths = _prepare_path_executors(path_definitions(definition),
            telescope, atmosphere, context)
        bindings = _prepare_controllable_optic_path_bindings(optics, paths)
        sampled_bindings = _prepare_sampled_aberration_path_bindings(
            sampled_aberrations, paths)
        acquisitions = _prepare_acquisition_owners(
            acquisition_definitions(definition), paths)
        rngs = _prepare_plant_rngs(
            atmosphere, paths, acquisitions, seed, version)
        return PreparedPlant(
            _PREPARED_PLANT_TOKEN,
            definition,
            target,
            context,
            telescope,
            atmosphere,
            optics,
            bindings,
            sampled_aberrations,
            sampled_bindings,
            prepared_endpoints,
            paths,
            acquisitions,
            rngs,
        )
    end
end

struct _PreparedAcquisitionSelectionToken end
const _PREPARED_ACQUISITION_SELECTION_TOKEN =
    _PreparedAcquisitionSelectionToken()

"""Canonical descriptor view into one plant-owned path registry."""
struct _PreparedPathSelection{
    P<:PreparedPlant,
    S<:FixedSizeVector{_PreparedPathSlot},
}
    plant::P
    slots::S
end

Base.size(selection::_PreparedPathSelection) = (length(selection.slots),)
Base.axes(selection::_PreparedPathSelection) = axes(selection.slots)
Base.length(selection::_PreparedPathSelection) = length(selection.slots)
Base.keys(selection::_PreparedPathSelection) = eachindex(selection.slots)
Base.eachindex(selection::_PreparedPathSelection) =
    eachindex(selection.slots)
Base.firstindex(selection::_PreparedPathSelection) =
    firstindex(selection.slots)
Base.lastindex(selection::_PreparedPathSelection) =
    lastindex(selection.slots)
Base.:(==)(left::_PreparedPathSelection, right::_PreparedPathSelection) =
    left.plant === right.plant && left.slots == right.slots

@inline function Base.getindex(
    selection::_PreparedPathSelection, index::Int)
    checkbounds(selection.slots, index)
    registry = getfield(selection.plant, :paths)
    return _prepared_path_value(
        registry, @inbounds(selection.slots[index]))
end

@inline function Base.iterate(
    selection::_PreparedPathSelection, state::Int=1)
    state > length(selection) && return nothing
    return (@inbounds selection[state], state + 1)
end

"""Canonical descriptor view into one plant-owned acquisition registry."""
struct _PreparedAcquisitionOwnerSelection{
    P<:PreparedPlant,
    S<:FixedSizeVector{_PreparedAcquisitionSlot},
}
    plant::P
    slots::S
end

Base.size(selection::_PreparedAcquisitionOwnerSelection) =
    (length(selection.slots),)
Base.axes(selection::_PreparedAcquisitionOwnerSelection) =
    axes(selection.slots)
Base.length(selection::_PreparedAcquisitionOwnerSelection) =
    length(selection.slots)
Base.keys(selection::_PreparedAcquisitionOwnerSelection) =
    eachindex(selection.slots)
Base.eachindex(selection::_PreparedAcquisitionOwnerSelection) =
    eachindex(selection.slots)
Base.firstindex(selection::_PreparedAcquisitionOwnerSelection) =
    firstindex(selection.slots)
Base.lastindex(selection::_PreparedAcquisitionOwnerSelection) =
    lastindex(selection.slots)
Base.:(==)(left::_PreparedAcquisitionOwnerSelection,
    right::_PreparedAcquisitionOwnerSelection) =
    left.plant === right.plant && left.slots == right.slots

@inline function Base.getindex(
    selection::_PreparedAcquisitionOwnerSelection, index::Int)
    checkbounds(selection.slots, index)
    registry = getfield(selection.plant, :acquisitions)
    return _prepared_acquisition_value(
        registry, @inbounds(selection.slots[index]))
end
@inline function Base.iterate(
    selection::_PreparedAcquisitionOwnerSelection, state::Int=1)
    state > length(selection) && return nothing
    return (@inbounds selection[state], state + 1)
end

"""
    PreparedAcquisitionSelection

Cold, canonical selection of acquisition owners and the unique optical paths
required by their full-optical providers. Reduced-order and synthetic/replay
providers do not select an otherwise unused path. The selected paths and
acquisitions use stable-ID order regardless of declaration or caller-selection
order. This value owns no schedule or independent RNG state; fixed
family/member descriptors resolve the exact selected plant-owned RNG groups.
"""
struct PreparedAcquisitionSelection{
    P<:PreparedPlant,
    V,
    W,
    B,
}
    plant::P
    paths::V
    acquisitions::W
    sampled_aberration_path_plans::B

    function PreparedAcquisitionSelection(
        ::_PreparedAcquisitionSelectionToken,
        plant::P,
        paths::V,
        acquisitions::W,
        sampled_aberration_path_plans::B,
    ) where {
        P<:PreparedPlant,
        V,
        W,
        B,
    }
        return new{P,V,W,B}(
            plant,
            paths,
            acquisitions,
            sampled_aberration_path_plans,
        )
    end
end

@inline prepared_paths(selection::PreparedAcquisitionSelection) =
    selection.paths
@inline prepared_acquisitions(selection::PreparedAcquisitionSelection) =
    selection.acquisitions

@inline _selected_acquisition_id(id::AcquisitionID) = id
@inline _selected_acquisition_id(name::Symbol) = AcquisitionID(name)

function _selected_acquisition_id(value)
    throw(PlantPreparationError(:acquisition, :invalid_selection,
        "selected acquisition identity must be a Symbol or AcquisitionID; got $(typeof(value))"))
end

function _requested_acquisition_ids(ids)
    isempty(ids) && throw(PlantPreparationError(:acquisition,
        :empty_selection, "an acquisition selection must not be empty"))
    requested = Set{AcquisitionID}()
    for value in ids
        id = _selected_acquisition_id(value)
        id in requested && throw(PlantPreparationError(:acquisition,
            :duplicate_selection,
            "acquisition selection contains duplicate identity $id"))
        push!(requested, id)
    end
    return requested
end

function _canonical_selected_acquisitions(plant::PreparedPlant, requested)
    selected = Int[]
    for index in eachindex(plant.acquisitions)
        owner = plant.acquisitions[index]
        acquisition_id(owner) in requested && push!(selected, index)
    end
    length(selected) == length(requested) || begin
        for id in requested
            prepared_acquisition(plant, id)
        end
        throw(PlantPreparationError(:acquisition, :unknown_id,
            "acquisition selection contains an unknown identity"))
    end
    sort!(selected; by=index ->
        acquisition_id(plant.acquisitions[index]).name)
    slots = Vector{_PreparedAcquisitionSlot}(undef, length(selected))
    @inbounds for index in eachindex(selected)
        slots[index] = plant.acquisitions.slots[selected[index]]
    end
    fixed = FixedSizeVectorDefault{_PreparedAcquisitionSlot}(slots)
    _require_prepared_acquisition_slots(
        plant.acquisitions.groups, fixed)
    return _PreparedAcquisitionOwnerSelection(plant, fixed)
end

function _canonical_selected_paths(plant::PreparedPlant,
    acquisitions::_PreparedAcquisitionOwnerSelection)
    requested = Set{OpticalPathID}()
    for acquisition in acquisitions
        _request_provider_path!(requested, acquisition,
            acquisition_provider_style(acquisition))
    end
    selected = Int[]
    for index in eachindex(plant.paths)
        path_id(plant.paths[index]) in requested && push!(selected, index)
    end
    length(selected) == length(requested) || throw(PlantPreparationError(
        :path, :prepared_binding,
        "selected acquisitions do not resolve to their prepared paths"))
    sort!(selected; by=index -> path_id(plant.paths[index]).name)
    slots = Vector{_PreparedPathSlot}(undef, length(selected))
    @inbounds for index in eachindex(selected)
        slots[index] = plant.paths.slots[selected[index]]
    end
    fixed = FixedSizeVectorDefault{_PreparedPathSlot}(slots)
    _require_prepared_path_slots(plant.paths.groups, fixed)
    return _PreparedPathSelection(plant, fixed)
end

@inline function _request_provider_path!(requested,
    acquisition::PreparedAcquisitionOwner, ::FullOpticalProviderStyle)
    push!(requested, acquisition_path_id(acquisition))
    return requested
end

@inline _request_provider_path!(requested,
    ::PreparedAcquisitionOwner,
    ::CommandResponsiveReducedOrderProviderStyle) = requested

@inline _request_provider_path!(requested,
    ::PreparedAcquisitionOwner, ::SyntheticReplayProviderStyle) = requested

function _prepared_path_rngs(plant::PreparedPlant,
    path::PreparedPathExecutor)
    index = Int(path.definition_slot)
    1 <= index <= length(plant.paths) &&
        plant.paths[index] === path && return plant.rngs.paths[index]
    throw(PlantPreparationError(:rng, :prepared_binding,
        "selected path has no exact prepared RNG owner"))
end

function _prepared_acquisition_rngs(plant::PreparedPlant,
    acquisition::PreparedAcquisitionOwner)
    index = Int(acquisition.definition_slot)
    1 <= index <= length(plant.acquisitions) &&
        plant.acquisitions[index] === acquisition &&
        return plant.rngs.acquisitions[index]
    throw(PlantPreparationError(:rng, :prepared_binding,
        "selected acquisition has no exact prepared RNG owner"))
end

struct _PreparedSampledAberrationPlanSlot
    family_slot::UInt32
    member_slot::UInt32
end

struct _PreparedSampledAberrationPlanFamily{
    P,
    V<:FixedSizeVector{P},
}
    values::V
end

struct _PreparedSampledAberrationPlanRegistry{
    G<:Tuple,
    S<:FixedSizeVector{_PreparedSampledAberrationPlanSlot},
}
    groups::G
    slots::S
end


@noinline function _prepared_sampled_aberration_plan_slot_error(
    family::Int, member::Int)
    throw(PlantPreparationError(:sampled_aberration, :plan_topology,
        "sampled-aberration plan family/member slot is out of bounds: ($family, $member)"))
end

Base.length(registry::_PreparedSampledAberrationPlanRegistry) =
    length(registry.slots)
Base.eachindex(registry::_PreparedSampledAberrationPlanRegistry) =
    eachindex(registry.slots)

struct _PreparedSampledAberrationPlanFamilyType{
    P,
} end

@inline _sampled_aberration_plan_family_matches(
    ::_PreparedSampledAberrationPlanFamilyType{P}, ::Type{P}) where {
    P,
} = true

@inline _sampled_aberration_plan_family_matches(
    ::_PreparedSampledAberrationPlanFamilyType,
    ::Type) = false

function _sampled_aberration_plan_family_types(plans)
    plan_types = DataType[]
    @inbounds for plan in plans
        plan_type = typeof(plan)::DataType
        plan_type in plan_types || push!(plan_types, plan_type)
    end
    sort!(plan_types; by=_prepared_owner_family_order_key)
    families = ()
    for plan_type in plan_types
        family_type = Core.apply_type(
            _PreparedSampledAberrationPlanFamilyType, plan_type)
        families = (families..., family_type())
    end
    return families
end

function _prepare_sampled_aberration_plan_family(
    ::_PreparedSampledAberrationPlanFamilyType{P}, plans,
) where {P}
    family_count = count(plan -> typeof(plan) === P, plans)
    values = Vector{P}(undef, family_count)
    next = 1
    @inbounds for plan in plans
        typeof(plan) === P || continue
        values[next] = plan
        next += 1
    end
    fixed = FixedSizeVectorDefault{P}(values)
    return _PreparedSampledAberrationPlanFamily{P,typeof(fixed)}(fixed)
end

@inline _prepare_sampled_aberration_plan_families(::Tuple{}, plans) = ()

function _prepare_sampled_aberration_plan_families(
    family_types::Tuple, plans)
    first = _prepare_sampled_aberration_plan_family(family_types[1], plans)
    rest = _prepare_sampled_aberration_plan_families(
        Base.tail(family_types), plans)
    return (first, rest...)
end

@inline function _sampled_aberration_plan_family_index(
    ::Tuple{}, ::Type, ::Int=1)
    return 0
end

@inline function _sampled_aberration_plan_family_index(
    families::Tuple, ::Type{P}, index::Int=1,
) where {P}
    _sampled_aberration_plan_family_matches(families[1], P) && return index
    return _sampled_aberration_plan_family_index(
        Base.tail(families), P, index + 1)
end

function _prepare_sampled_aberration_plan_slots(family_types, plans)
    counts = zeros(UInt32, length(family_types))
    slots = Vector{_PreparedSampledAberrationPlanSlot}(undef, length(plans))
    @inbounds for index in eachindex(plans)
        family = _sampled_aberration_plan_family_index(
            family_types, typeof(plans[index]))
        iszero(family) && throw(PlantPreparationError(
            :sampled_aberration, :missing_plan_family,
            "sampled-aberration path plan has no exact family"))
        counts[family] == typemax(UInt32) && throw(PlantPreparationError(
            :sampled_aberration, :family_capacity,
            "sampled-aberration path-plan family exceeds UInt32 capacity"))
        counts[family] += UInt32(1)
        slots[index] = _PreparedSampledAberrationPlanSlot(
            UInt32(family), counts[family])
    end
    return FixedSizeVectorDefault{_PreparedSampledAberrationPlanSlot}(slots)
end


function _require_sampled_aberration_plan_slots(groups, slots)
    @inbounds for slot in slots
        family = Int(slot.family_slot)
        1 <= family <= length(groups) ||
            _prepared_sampled_aberration_plan_slot_error(
                family, Int(slot.member_slot))
        checkbounds(groups[family].values, Int(slot.member_slot))
    end
    return slots
end

function _selected_sampled_aberration_path_plans(
    plant::PreparedPlant,
    paths::_PreparedPathSelection,
)
    bindings = plant.sampled_aberration_path_bindings
    plans = Vector{_PreparedSampledAberrationPathPlan}(undef, length(paths))
    @inbounds for index in eachindex(paths)
        plans[index] = _prepare_sampled_aberration_path_plan(
            plant.sampled_aberrations,
            bindings,
            path_id(paths[index]),
        )
    end
    family_types = _sampled_aberration_plan_family_types(plans)
    groups = _prepare_sampled_aberration_plan_families(family_types, plans)
    slots = _prepare_sampled_aberration_plan_slots(family_types, plans)
    _require_sampled_aberration_plan_slots(groups, slots)
    return _PreparedSampledAberrationPlanRegistry(groups, slots)
end

function _prepare_acquisition_selection(plant::PreparedPlant, ids)
    requested = _requested_acquisition_ids(ids)
    acquisitions = _canonical_selected_acquisitions(plant, requested)
    paths = _canonical_selected_paths(plant, acquisitions)
    sampled_aberration_path_plans =
        _selected_sampled_aberration_path_plans(plant, paths)
    for path in paths
        _prepared_path_rngs(plant, path)
    end
    for acquisition in acquisitions
        _prepared_acquisition_rngs(plant, acquisition)
    end
    return PreparedAcquisitionSelection(
        _PREPARED_ACQUISITION_SELECTION_TOKEN, plant, paths, acquisitions,
        sampled_aberration_path_plans)
end

"""
    prepare_acquisition_selection(plant, ids)

Resolve a nonempty tuple or vector of acquisition identities to a canonical,
deduplicated path set and exact acquisition owners. Duplicate and unknown
identities are rejected during this cold operation.
"""
prepare_acquisition_selection(plant::PreparedPlant, ids::Tuple) =
    _prepare_acquisition_selection(plant, ids)

prepare_acquisition_selection(plant::PreparedPlant, ids::AbstractVector) =
    _prepare_acquisition_selection(plant, ids)

prepare_acquisition_selection(plant::PreparedPlant,
    id::Union{Symbol,AcquisitionID}) =
    _prepare_acquisition_selection(plant, (id,))

function prepare_acquisition_selection(plant::PreparedPlant, ids)
    throw(PlantPreparationError(:acquisition, :invalid_selection,
        "acquisition selection must be one identity, a Tuple, or an AbstractVector; got $(typeof(ids))"))
end

Base.@noinline function _require_selected_path_binding(
    path::PreparedPathExecutor, atmosphere, context)
    path.atmosphere === atmosphere || throw(PlantPreparationError(:path,
        :prepared_binding,
        "selected path does not retain the prepared plant atmosphere"))
    path.context === context || throw(PlantPreparationError(:path,
        :prepared_binding,
        "selected path does not retain the prepared plant execution context"))
    _require_current_path_binding(path)
    return nothing
end

@inline function _selected_path_family_operation!(operation,
    ::Tuple{}, family::Int, member::Int, arguments::Tuple)
    return _prepared_path_slot_error(family, member)
end

@inline function _selected_path_family_operation!(operation,
    groups::Tuple, family::Int, member::Int, arguments::Tuple)
    family == 1 && return operation(
        @inbounds(groups[1].values[member]), arguments)
    return _selected_path_family_operation!(operation, Base.tail(groups),
        family - 1, member, arguments)
end

@inline function _selected_path_binding!(path, arguments::Tuple)
    atmosphere, context = arguments
    return _require_selected_path_binding(path, atmosphere, context)
end

@inline function _selected_path_materialization!(path, arguments::Tuple)
    atmosphere, epoch = arguments
    return _validate_selected_materialization(path, atmosphere, epoch)
end

@inline function _selected_path_operation!(operation,
    selection::_PreparedPathSelection, index::Int, arguments::Tuple)
    slot = @inbounds selection.slots[index]
    groups = getfield(selection.plant, :paths).groups
    return _selected_path_family_operation!(operation, groups,
        Int(slot.family_slot), Int(slot.member_slot), arguments)
end

function _require_selected_path_bindings(
    paths::_PreparedPathSelection, atmosphere, context)
    @inbounds for index in eachindex(paths)
        _selected_path_operation!(_selected_path_binding!, paths, index,
            (atmosphere, context))
    end
    return nothing
end

Base.@noinline function _require_selected_acquisition_binding(
    acquisition::PreparedAcquisitionOwner, context)
    acquisition.context === context || throw(PlantPreparationError(
        :acquisition, :prepared_binding,
        "selected acquisition does not retain the prepared plant execution context"))
    validate_acquisition_provider_binding(acquisition.provider,
        acquisition.path_result)
    return nothing
end

@inline function _selected_acquisition_family_operation!(operation,
    ::Tuple{}, family::Int, member::Int, arguments::Tuple)
    return _prepared_acquisition_slot_error(family, member)
end

@inline function _selected_acquisition_family_operation!(operation,
    groups::Tuple, family::Int, member::Int, arguments::Tuple)
    family == 1 && return operation(
        @inbounds(groups[1].values[member]), arguments)
    return _selected_acquisition_family_operation!(operation,
        Base.tail(groups), family - 1, member, arguments)
end

@inline function _selected_acquisition_binding!(owner, arguments::Tuple)
    context = arguments[1]
    return _require_selected_acquisition_binding(owner, context)
end

@inline function _selected_acquisition_operation!(operation,
    selection::_PreparedAcquisitionOwnerSelection, index::Int,
    arguments::Tuple)
    slot = @inbounds selection.slots[index]
    groups = getfield(selection.plant, :acquisitions).groups
    return _selected_acquisition_family_operation!(operation, groups,
        Int(slot.family_slot), Int(slot.member_slot), arguments)
end

function _require_selected_acquisition_bindings(
    acquisitions::_PreparedAcquisitionOwnerSelection, context)
    @inbounds for index in eachindex(acquisitions)
        _selected_acquisition_operation!(_selected_acquisition_binding!,
            acquisitions, index, (context,))
    end
    return nothing
end

function _require_selection_bindings(selection::PreparedAcquisitionSelection)
    atmosphere = prepared_atmosphere(selection.plant)
    context = selection.plant.context
    _require_selected_path_bindings(selection.paths, atmosphere, context)
    _require_selected_acquisition_bindings(selection.acquisitions, context)
    validate_atmosphere_rng_binding(selection.plant.rngs.atmosphere,
        atmosphere)
    _require_rng_registry_binding(selection.plant.rngs.paths,
        selection.plant.paths, :path)
    _require_rng_registry_binding(selection.plant.rngs.acquisitions,
        selection.plant.acquisitions, :acquisition)
    _require_selected_path_rng_owner_bindings(selection.plant,
        selection.paths)
    _require_selected_acquisition_rng_owner_bindings(selection.plant,
        selection.acquisitions)
    return atmosphere
end

@inline _require_rng_registry_group_cardinality(
    ::Tuple{}, ::Tuple{}, ::Symbol) = nothing

Base.@noinline function _require_rng_registry_group_cardinality(
    ::Tuple{}, ::Tuple, category::Symbol)
    throw(PlantPreparationError(:rng, :owner_topology,
        "prepared $category RNG families exceed owner families"))
end

Base.@noinline function _require_rng_registry_group_cardinality(
    ::Tuple, ::Tuple{}, category::Symbol)
    throw(PlantPreparationError(:rng, :owner_topology,
        "prepared $category RNG families do not cover all owner families"))
end

@inline function _require_rng_registry_group_cardinality(
    rng_groups::Tuple, owner_groups::Tuple, category::Symbol)
    length(rng_groups[1].values) == length(owner_groups[1].values) ||
        throw(PlantPreparationError(:rng, :owner_topology,
            "prepared $category RNG family cardinality does not match its owner family"))
    return _require_rng_registry_group_cardinality(
        Base.tail(rng_groups), Base.tail(owner_groups), category)
end

function _require_rng_registry_binding(rngs::_PreparedOwnerRNGRegistry,
    owners, category::Symbol)
    rngs.slots === owners.slots || throw(PlantPreparationError(
        :rng, :owner_topology,
        "prepared $category RNG descriptors do not retain their exact owner registry"))
    length(rngs.groups) == length(owners.groups) || throw(
        PlantPreparationError(:rng, :owner_topology,
            "prepared $category RNG families do not match owner families"))
    _require_rng_registry_group_cardinality(
        rngs.groups, owners.groups, category)
    return rngs
end

function _require_selected_path_rng_owner_bindings(plant, owners)
    groups = getfield(owners.plant, :paths).groups
    rng_groups = plant.rngs.paths.groups
    @inbounds for index in eachindex(owners)
        slot = owners.slots[index]
        _selected_path_rng_binding_family!(groups, rng_groups,
            Int(slot.family_slot), Int(slot.member_slot))
    end
    return nothing
end

function _require_selected_acquisition_rng_owner_bindings(
    plant, owners)
    groups = getfield(owners.plant, :acquisitions).groups
    rng_groups = plant.rngs.acquisitions.groups
    @inbounds for index in eachindex(owners)
        slot = owners.slots[index]
        _selected_acquisition_rng_binding_family!(groups, rng_groups,
            Int(slot.family_slot), Int(slot.member_slot))
    end
    return nothing
end

@inline function _selected_path_rng_binding_family!(
    ::Tuple{}, ::Tuple{}, family::Int, member::Int)
    return _prepared_path_slot_error(family, member)
end

Base.@noinline function _selected_path_rng_binding_family!(
    groups::Tuple, rng_groups::Tuple, family::Int, member::Int)
    family == 1 && return _require_selected_rng_owner_binding(
        @inbounds(rng_groups[1].values[member]),
        @inbounds(groups[1].values[member]))
    return _selected_path_rng_binding_family!(Base.tail(groups),
        Base.tail(rng_groups), family - 1, member)
end

@inline function _selected_acquisition_rng_binding_family!(
    ::Tuple{}, ::Tuple{}, family::Int, member::Int)
    return _prepared_acquisition_slot_error(family, member)
end

Base.@noinline function _selected_acquisition_rng_binding_family!(
    groups::Tuple, rng_groups::Tuple, family::Int, member::Int)
    family == 1 && return _require_selected_rng_owner_binding(
        @inbounds(rng_groups[1].values[member]),
        @inbounds(groups[1].values[member]))
    return _selected_acquisition_rng_binding_family!(Base.tail(groups),
        Base.tail(rng_groups), family - 1, member)
end

Base.@noinline function _require_selected_rng_owner_binding(
    rngs::PreparedOwnerRNGs, owner)
    _require_rng_owner_binding(rngs, owner)
    return nothing
end

Base.@noinline function _validate_selected_materialization(
    path::PreparedPathExecutor, atmosphere, epoch)
    validate_path_materialization(path.materialization, path.input,
        atmosphere, epoch)
    return nothing
end

function _validate_selected_materializations(paths::_PreparedPathSelection,
    atmosphere, epoch)
    @inbounds for index in eachindex(paths)
        _selected_path_operation!(_selected_path_materialization!, paths,
            index, (atmosphere, epoch))
    end
    return nothing
end

Base.@noinline function _materialize_selected_path!(
    path::PreparedPathExecutor, rngs::PreparedOwnerRNGs, atmosphere, epoch)
    materialize_path_input_rngs!(path.materialization, path.input,
        atmosphere, epoch, rngs)
    return nothing
end

@inline function _selected_path_materialize_family!(
    ::Tuple{}, ::Tuple{}, family::Int, member::Int,
    atmosphere, epoch)
    return _prepared_path_slot_error(family, member)
end

Base.@noinline function _selected_path_materialize_family!(
    groups::Tuple, rng_groups::Tuple, family::Int, member::Int,
    atmosphere, epoch)
    family == 1 && return _materialize_selected_path!(
        @inbounds(groups[1].values[member]),
        @inbounds(rng_groups[1].values[member]), atmosphere, epoch)
    return _selected_path_materialize_family!(Base.tail(groups),
        Base.tail(rng_groups), family - 1, member, atmosphere, epoch)
end

function _materialize_selected_paths!(plant, paths::_PreparedPathSelection,
    atmosphere, epoch)
    groups = getfield(paths.plant, :paths).groups
    rng_groups = plant.rngs.paths.groups
    @inbounds for index in eachindex(paths)
        slot = paths.slots[index]
        _selected_path_materialize_family!(groups, rng_groups,
            Int(slot.family_slot), Int(slot.member_slot),
            atmosphere, epoch)
    end
    return nothing
end

Base.@noinline function _apply_selected_sampled_aberration!(
    path::PreparedPathExecutor, plan)
    _apply_sampled_aberration_path_plan!(path.input, plan)
    return nothing
end

@inline function _selected_path_sampled_aberration_family!(
    ::Tuple{}, family::Int, member::Int,
    plan)
    return _prepared_path_slot_error(family, member)
end

Base.@noinline function _selected_path_sampled_aberration_family!(
    groups::Tuple, family::Int, member::Int,
    plan)
    family == 1 && return _apply_selected_sampled_aberration!(
        @inbounds(groups[1].values[member]), plan)
    return _selected_path_sampled_aberration_family!(Base.tail(groups),
        family - 1, member, plan)
end

@inline function _selected_sampled_aberration_plan_family!(
    path_groups::Tuple, ::Tuple{}, family::Int, member::Int,
    path_family::Int, path_member::Int)
    return _prepared_sampled_aberration_plan_slot_error(family, member)
end

Base.@noinline function _selected_sampled_aberration_plan_family!(
    path_groups::Tuple, plan_groups::Tuple, family::Int, member::Int,
    path_family::Int, path_member::Int)
    family == 1 && return _selected_path_sampled_aberration_family!(
        path_groups, path_family, path_member,
        @inbounds(plan_groups[1].values[member]))
    return _selected_sampled_aberration_plan_family!(path_groups,
        Base.tail(plan_groups), family - 1, member,
        path_family, path_member)
end

function _apply_selected_sampled_aberrations!(
    paths::_PreparedPathSelection,
    plans::_PreparedSampledAberrationPlanRegistry,
)
    length(paths) == length(plans) || throw(PlantPreparationError(
        :sampled_aberration, :binding_topology,
        "selected sampled-aberration plans do not match selected paths"))
    groups = getfield(paths.plant, :paths).groups
    plan_groups = plans.groups
    @inbounds for index in eachindex(paths)
        path_slot = paths.slots[index]
        plan_slot = plans.slots[index]
        _selected_sampled_aberration_plan_family!(groups, plan_groups,
            Int(plan_slot.family_slot), Int(plan_slot.member_slot),
            Int(path_slot.family_slot), Int(path_slot.member_slot))
    end
    return nothing
end

Base.@noinline function _execute_selected_path!(
    path::PreparedPathExecutor, rngs::PreparedOwnerRNGs)
    _execute_path_in_context!(path, rngs)
    return nothing
end

@inline function _selected_path_execute_family!(
    ::Tuple{}, ::Tuple{}, family::Int, member::Int)
    return _prepared_path_slot_error(family, member)
end

Base.@noinline function _selected_path_execute_family!(
    groups::Tuple, rng_groups::Tuple, family::Int, member::Int)
    family == 1 && return _execute_selected_path!(
        @inbounds(groups[1].values[member]),
        @inbounds(rng_groups[1].values[member]))
    return _selected_path_execute_family!(Base.tail(groups),
        Base.tail(rng_groups), family - 1, member)
end

function _execute_selected_paths!(plant, paths::_PreparedPathSelection)
    groups = getfield(paths.plant, :paths).groups
    rng_groups = plant.rngs.paths.groups
    @inbounds for index in eachindex(paths)
        slot = paths.slots[index]
        _selected_path_execute_family!(groups, rng_groups,
            Int(slot.family_slot), Int(slot.member_slot))
    end
    return nothing
end

Base.@noinline function _execute_selected_acquisition!(
    acquisition::PreparedAcquisitionOwner, rngs::PreparedOwnerRNGs)
    _execute_acquisition_in_context!(acquisition, rngs)
    return nothing
end

@inline function _selected_acquisition_execute_family!(
    ::Tuple{}, ::Tuple{}, family::Int, member::Int)
    return _prepared_acquisition_slot_error(family, member)
end

Base.@noinline function _selected_acquisition_execute_family!(
    groups::Tuple, rng_groups::Tuple, family::Int, member::Int)
    family == 1 && return _execute_selected_acquisition!(
        @inbounds(groups[1].values[member]),
        @inbounds(rng_groups[1].values[member]))
    return _selected_acquisition_execute_family!(Base.tail(groups),
        Base.tail(rng_groups), family - 1, member)
end

function _execute_selected_acquisitions!(
    plant, acquisitions::_PreparedAcquisitionOwnerSelection)
    groups = getfield(acquisitions.plant, :acquisitions).groups
    rng_groups = plant.rngs.acquisitions.groups
    @inbounds for index in eachindex(acquisitions)
        slot = acquisitions.slots[index]
        _selected_acquisition_execute_family!(groups, rng_groups,
            Int(slot.family_slot), Int(slot.member_slot))
    end
    return nothing
end

@inline _require_selection_atmosphere(
    atmosphere::AbstractTimedAtmosphere) = atmosphere

function _require_selection_atmosphere(atmosphere::AbstractAtmosphere)
    throw(PlantPreparationError(:plant, :untimed_atmosphere,
        "selected current-epoch execution requires a timed plant atmosphere"))
end

function _validate_selection_epoch!(selection::PreparedAcquisitionSelection,
    atmosphere::AbstractTimedAtmosphere, epoch::AtmosphereEpoch)
    _validate_epoch_identity(atmosphere_identity(atmosphere), atmosphere,
        epoch)
    _validate_selected_materializations(selection.paths, atmosphere, epoch)
    return epoch
end

function _execute_selected_epoch!(selection::PreparedAcquisitionSelection,
    atmosphere::AbstractTimedAtmosphere, epoch::AtmosphereEpoch)
    _materialize_selected_paths!(
        selection.plant, selection.paths, atmosphere, epoch)
    _apply_selected_sampled_aberrations!(
        selection.paths,
        selection.sampled_aberration_path_plans,
    )
    _execute_selected_paths!(selection.plant, selection.paths)
    _execute_selected_acquisitions!(
        selection.plant, selection.acquisitions)
    return selection
end

function _execute_acquisition_selection_in_context!(
    selection::PreparedAcquisitionSelection,
    epoch::AtmosphereEpoch)
    atmosphere = _require_selection_atmosphere(
        _require_selection_bindings(selection))
    _validate_selection_epoch!(selection, atmosphere, epoch)
    return _execute_selected_epoch!(selection, atmosphere, epoch)
end

"""
    execute_acquisition_selection!(selection, epoch)

Materialize every unique selected path from one explicit current atmosphere
epoch, then form each path exactly once and execute the canonical acquisition
registry with the exact owner-derived RNG streams prepared by `prepare_plant`.
"""
function execute_acquisition_selection!(
    selection::PreparedAcquisitionSelection,
    epoch::AtmosphereEpoch)
    result = _with_completed_prepared_device_execution_context(
        selection.plant.context) do
        _execute_acquisition_selection_in_context!(selection, epoch)
    end
    return result::typeof(selection)
end

function _execute_acquisition_selection_at_in_context!(
    selection::PreparedAcquisitionSelection,
    model_time::Real)
    atmosphere = _require_selection_atmosphere(
        _require_selection_bindings(selection))
    atmosphere_rng = _prepared_atmosphere_rng(atmosphere,
        selection.plant.rngs.atmosphere)
    epoch = advance_to!(atmosphere, model_time, atmosphere_rng)
    _validate_selection_epoch!(selection, atmosphere, epoch)
    return _execute_selected_epoch!(selection, atmosphere, epoch)
end

"""
    execute_acquisition_selection_at!(selection, model_time)

Preflight all selected owner bindings, advance the one shared atmosphere to an
explicit absolute model time, and execute the selection from the newly current
epoch. Equal-time execution reuses the current publication. Atmosphere and
acquisition randomness comes from exact streams owned by the prepared plant.
"""
function execute_acquisition_selection_at!(
    selection::PreparedAcquisitionSelection,
    model_time::Real)
    result = _with_completed_prepared_device_execution_context(
        selection.plant.context) do
        _execute_acquisition_selection_at_in_context!(selection, model_time)
    end
    return result::typeof(selection)
end
