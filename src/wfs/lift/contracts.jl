abstract type AbstractLiFTObservationMapping end

"""No deterministic spatial mapping between focal-plane rate and observation."""
struct LiFTIdentityMapping <: AbstractLiFTObservationMapping end

"""
    LiFTFrameMapping(; response=NullFrameResponse(), sampling=1, binning=1)

Deterministic spatial preprocessing shared by a LiFT forward model and an
acquisition path. `response` is applied first on the optical grid, followed by
cell-summing `sampling` and `binning`. QE, exposure, noise, gain, readout
windowing, and cadence are deliberately not part of this mapping.
"""
struct LiFTFrameMapping{R<:AbstractFrameResponse} <: AbstractLiFTObservationMapping
    response::R
    sampling::Int
    binning::Int
end

function LiFTFrameMapping(; response::AbstractFrameResponse=NullFrameResponse(),
    sampling::Int=1, binning::Int=1)
    sampling >= 1 || throw(InvalidConfiguration(
        "LiFT frame-mapping sampling must be >= 1"))
    binning >= 1 || throw(InvalidConfiguration(
        "LiFT frame-mapping binning must be >= 1"))
    _require_prepared_response_sampling(response, sampling)
    validate_frame_response_model(response)
    return LiFTFrameMapping{typeof(response)}(response, sampling, binning)
end

abstract type AbstractLiFTObservationDomain end

"""
Photon-arrival-rate observations. The noise-equivalent exposure is acquisition
metadata for calibration weighting; it does not scale the values.
"""
struct LiFTPhotonRate{T<:AbstractFloat} <: AbstractLiFTObservationDomain
    noise_equivalent_exposure_s::T
    quantum_efficiency::T
end

function LiFTPhotonRate(; noise_equivalent_exposure_s::Real=1.0,
    quantum_efficiency::Real=1.0)
    exposure, qe = promote(float(noise_equivalent_exposure_s),
        float(quantum_efficiency))
    _require_lift_exposure_qe(exposure, qe, "LiFT photon-rate observation")
    return LiFTPhotonRate{typeof(exposure)}(exposure, qe)
end

"""Expected detected counts formed using an explicit exposure and QE."""
struct LiFTExpectedCounts{T<:AbstractFloat} <: AbstractLiFTObservationDomain
    exposure_duration_s::T
    quantum_efficiency::T
end

function LiFTExpectedCounts(exposure_duration_s::Real;
    quantum_efficiency::Real=1.0)
    exposure, qe = promote(float(exposure_duration_s), float(quantum_efficiency))
    _require_lift_exposure_qe(exposure, qe, "LiFT expected-count observation")
    return LiFTExpectedCounts{typeof(exposure)}(exposure, qe)
end

"""
Dimensionless relative intensity with an explicit photon-rate value per native
unit. The noise-equivalent exposure is acquisition metadata for calibration.
"""
struct LiFTNormalizedIntensity{T<:AbstractFloat} <: AbstractLiFTObservationDomain
    photon_rate_per_unit::T
    noise_equivalent_exposure_s::T
    quantum_efficiency::T
end

function LiFTNormalizedIntensity(photon_rate_per_unit::Real;
    noise_equivalent_exposure_s::Real=1.0,
    quantum_efficiency::Real=1.0)
    scale, exposure, qe = promote(float(photon_rate_per_unit),
        float(noise_equivalent_exposure_s), float(quantum_efficiency))
    isfinite(scale) && scale > zero(scale) || throw(InvalidConfiguration(
        "LiFT normalized-intensity photon_rate_per_unit must be finite and > 0"))
    _require_lift_exposure_qe(exposure, qe,
        "LiFT normalized-intensity observation")
    return LiFTNormalizedIntensity{typeof(scale)}(scale, exposure, qe)
end

@inline function _require_lift_exposure_qe(exposure::T, qe::T,
    label::AbstractString) where {T<:AbstractFloat}
    isfinite(exposure) && exposure > zero(T) || throw(InvalidConfiguration(
        "$label exposure must be finite and > 0"))
    isfinite(qe) && zero(T) < qe <= one(T) || throw(InvalidConfiguration(
        "$label quantum efficiency must be finite and lie in (0, 1]"))
    return nothing
end

@inline lift_observation_to_rate_scale(::LiFTPhotonRate,
    ::Type{T}) where {T<:AbstractFloat} = one(T)
@inline lift_observation_to_rate_scale(domain::LiFTExpectedCounts,
    ::Type{T}) where {T<:AbstractFloat} =
    inv(T(domain.exposure_duration_s) * T(domain.quantum_efficiency))
@inline lift_observation_to_rate_scale(domain::LiFTNormalizedIntensity,
    ::Type{T}) where {T<:AbstractFloat} = T(domain.photon_rate_per_unit)

struct LiFTObservationContract{M<:OpticalPlaneMetadata,S}
    rate_metadata::M
    preprocessing_signature::S
end

struct LiFTObservationMetadata{T<:AbstractFloat,
    C<:LiFTObservationContract,D<:AbstractLiFTObservationDomain,E,
    B<:AbstractArrayBackend,PD<:AbstractComputeDevice}
    contract::C
    domain::D
    readout_noise_std::T
    numeric_type::Type{E}
    backend::B
    device::PD
end

"""Caller-owned acquired data plus its explicit LiFT observation contract."""
struct LiFTObservation{M<:LiFTObservationMetadata,A<:AbstractMatrix}
    metadata::M
    values::A
end

struct LiFTDenseObjectKernel{T<:AbstractFloat,A<:AbstractMatrix{T}}
    kernel::A
    inv_norm::T
end

struct LiFTSeparableObjectKernel{T<:AbstractFloat,V<:AbstractVector{T}}
    row::V
    col::V
    inv_norm::T
end

struct LiFTForwardPlan{T<:AbstractFloat,
    PM<:AbstractMatrix{Bool},PA<:AbstractMatrix{T},B<:AbstractArray{T,3},
    D<:AbstractMatrix{T},K,M<:AbstractLiFTObservationMapping,
    C<:LiFTObservationContract}
    pupil_mask::PM
    pupil_amplitude::PA
    basis::B
    diversity_opd::D
    wavelength_m::T
    photon_irradiance::T
    pupil_cell_area_m2::T
    focal_resolution::Int
    zero_padding::Int
    object_kernel::K
    mapping::M
    observation_contract::C
end

struct LiFTForwardWorkspace{W<:Workspace,B<:AbstractMatrix,
    C<:AbstractMatrix,RB,SB,OB,CB}
    propagation::W
    optical_rate_buffer::B
    amplitude_buffer::B
    field_scratch::B
    focal_buffer::C
    mode_buffer::C
    conjugate_field_buffer::C
    response_buffer::RB
    response_scratch::RB
    sampling_buffer::SB
    mapped_rate_buffer::OB
    convolution_buffer::CB
    convolution_scratch::CB
end

"""Exact single-writer owner for one LiFT forward-model input and output."""
struct PreparedLiFTForward{M<:LiFTForwardPlan,
    W<:LiFTForwardWorkspace,I<:AbstractMatrix,O<:IntensityMap,B,D}
    plan::M
    workspace::W
    input::I
    output::O
    backend::B
    device::D
end

@inline _lift_mapping_factors(::LiFTIdentityMapping) = (1, 1)
@inline _lift_mapping_factors(mapping::LiFTFrameMapping) =
    (mapping.sampling, mapping.binning)

@inline _lift_mapping_signature(::LiFTIdentityMapping) = (:identity,)

function _lift_array_signature(array::AbstractArray)
    host = Array(array)
    signature = hash(size(host), UInt(0))
    @inbounds for value in host
        signature = hash(value, signature)
    end
    return signature
end

@inline _lift_response_signature(::NullFrameResponse) = (:none,)
@inline _lift_response_signature(model::GaussianPixelResponse) =
    (:gaussian, model.response_width_px, size(model.kernel),
        _lift_array_signature(model.kernel))
@inline _lift_response_signature(model::SampledFrameResponse) =
    (:sampled, size(model.kernel), _lift_array_signature(model.kernel))
@inline _lift_response_signature(model::RectangularPixelAperture) =
    (:rectangular_aperture, model.pitch_x_px, model.pitch_y_px,
        model.fill_factor_x, model.fill_factor_y, size(model.kernel_x),
        size(model.kernel_y), _lift_array_signature(model.kernel_x),
        _lift_array_signature(model.kernel_y))

@inline function _lift_mapping_signature(mapping::LiFTFrameMapping)
    return (:frame, mapping.sampling, mapping.binning,
        _lift_response_signature(mapping.response))
end

@inline _copy_lift_response(::NullFrameResponse) = NullFrameResponse()

function _copy_lift_response(response::GaussianPixelResponse{T}) where {T}
    return GaussianPixelResponse{T,typeof(response.kernel)}(
        response.response_width_px, response.kernel)
end

function _copy_lift_response(response::SampledFrameResponse{T}) where {T}
    return SampledFrameResponse{T,typeof(response.kernel)}(response.kernel)
end

function _copy_lift_response(response::RectangularPixelAperture{T}) where {T}
    return RectangularPixelAperture{T,typeof(response.kernel_x),
        typeof(response.kernel_y)}(response.pitch_x_px, response.pitch_y_px,
        response.fill_factor_x, response.fill_factor_y, response.kernel_x,
        response.kernel_y)
end

@inline _prepare_lift_mapping(mapping::LiFTIdentityMapping) = mapping

function _prepare_lift_mapping(mapping::LiFTFrameMapping)
    return LiFTFrameMapping(response=_copy_lift_response(mapping.response),
        sampling=mapping.sampling, binning=mapping.binning)
end

@inline _lift_output_dimensions(focal_resolution::Int,
    ::LiFTIdentityMapping) = (focal_resolution, focal_resolution)

function _lift_output_dimensions(focal_resolution::Int,
    mapping::LiFTFrameMapping)
    divisor = mapping.sampling * mapping.binning
    focal_resolution % divisor == 0 || throw(DimensionMismatchError(
        "LiFT focal resolution must be divisible by sampling * binning"))
    resolution = div(focal_resolution, divisor)
    return (resolution, resolution)
end

@inline _require_lift_response_backend(::NullFrameResponse,
    ::AbstractMatrix) = nothing

function _require_lift_response_array(array::AbstractArray,
    template::AbstractMatrix, ::Type{T}) where {T<:AbstractFloat}
    eltype(array) === T || throw(InvalidConfiguration(
        "LiFT response and forward model must use the same numeric type"))
    typeof(backend(array)) === typeof(backend(template)) || throw(
        InvalidConfiguration(
            "LiFT response and forward model must use the same array backend"))
    compute_device(array) == compute_device(template) || throw(
        InvalidConfiguration(
            "LiFT response and forward model must occupy the same compute device"))
    return nothing
end

@inline function _require_lift_response_backend(model::GaussianPixelResponse,
    template::AbstractMatrix{T}) where {T<:AbstractFloat}
    return _require_lift_response_array(model.kernel, template, T)
end

@inline function _require_lift_response_backend(model::SampledFrameResponse,
    template::AbstractMatrix{T}) where {T<:AbstractFloat}
    return _require_lift_response_array(model.kernel, template, T)
end

@inline function _require_lift_response_backend(model::RectangularPixelAperture,
    template::AbstractMatrix{T}) where {T<:AbstractFloat}
    _require_lift_response_array(model.kernel_x, template, T)
    return _require_lift_response_array(model.kernel_y, template, T)
end

@inline _require_lift_mapping_backend(::LiFTIdentityMapping,
    ::AbstractMatrix) = nothing
@inline _require_lift_mapping_backend(mapping::LiFTFrameMapping,
    template::AbstractMatrix) =
    _require_lift_response_backend(mapping.response, template)

function _copy_lift_array(template::AbstractArray, input::AbstractArray,
    ::Type{T}) where {T<:AbstractFloat}
    output = similar(template, T, size(input)...)
    copyto!(output, input)
    return output
end

function _prepare_lift_object_kernel(object_kernel, template::AbstractMatrix,
    ::Type{T}) where {T<:AbstractFloat}
    object_kernel === nothing && return nothing
    ndims(object_kernel) == 2 || throw(DimensionMismatchError(
        "LiFT object kernel must be a matrix"))
    host_kernel = T.(Array(object_kernel))
    all(isfinite, host_kernel) || throw(InvalidConfiguration(
        "LiFT object kernel values must be finite"))
    all(>=(zero(T)), host_kernel) || throw(InvalidConfiguration(
        "LiFT object kernel values must be nonnegative"))
    sum(host_kernel) > zero(T) || throw(InvalidConfiguration(
        "LiFT object kernel must have positive total intensity"))
    kernel = similar(template, T, size(host_kernel)...)
    copyto!(kernel, host_kernel)
    return _lift_object_kernel(kernel)
end

function _allocate_lift_mapping_buffers(template::AbstractMatrix{T},
    focal_resolution::Int, ::LiFTIdentityMapping) where {T<:AbstractFloat}
    return (nothing, nothing, nothing, nothing)
end

function _allocate_lift_mapping_buffers(template::AbstractMatrix{T},
    focal_resolution::Int, mapping::LiFTFrameMapping) where {T<:AbstractFloat}
    sampled_resolution = div(focal_resolution, mapping.sampling)
    output_resolution = div(sampled_resolution, mapping.binning)
    response = similar(template, T, focal_resolution, focal_resolution)
    response_scratch = similar(response)
    sampled = similar(template, T, sampled_resolution, sampled_resolution)
    mapped = similar(template, T, output_resolution, output_resolution)
    return (response, response_scratch, sampled, mapped)
end

function _allocate_lift_forward_workspace(template::AbstractMatrix{T},
    pupil_resolution::Int, focal_resolution::Int, zero_padding::Int,
    object_kernel, mapping::AbstractLiFTObservationMapping) where {T<:AbstractFloat}
    oversampling = lift_oversampling(zero_padding)
    propagation = Workspace(template,
        lift_pad_size(pupil_resolution, zero_padding); T=T)
    optical_rate = similar(template, T, focal_resolution, focal_resolution)
    pupil_amplitude = similar(template, T, pupil_resolution, pupil_resolution)
    focal_size = focal_resolution * oversampling
    field_scratch = similar(optical_rate, T, focal_size, focal_size)
    focal = similar(optical_rate, Complex{T}, focal_size, focal_size)
    mode = similar(focal)
    conjugate_field = similar(focal)
    response, response_scratch, sampled, mapped =
        _allocate_lift_mapping_buffers(template, focal_resolution, mapping)
    convolution = object_kernel === nothing ? nothing : similar(optical_rate)
    convolution_scratch = object_kernel === nothing ? nothing : similar(optical_rate)
    return LiFTForwardWorkspace(propagation, optical_rate,
        pupil_amplitude, field_scratch, focal, mode, conjugate_field,
        response, response_scratch, sampled, mapped, convolution,
        convolution_scratch)
end

function _require_lift_forward_input(plan::LiFTForwardPlan,
    input::AbstractMatrix)
    size(input) == size(plan.diversity_opd) || throw(DimensionMismatchError(
        "LiFT forward OPD must match the prepared pupil dimensions"))
    eltype(input) === eltype(plan.diversity_opd) || throw(InvalidConfiguration(
        "LiFT forward OPD must use the prepared numeric type"))
    typeof(backend(input)) === typeof(backend(plan.diversity_opd)) || throw(
        InvalidConfiguration("LiFT forward OPD must use the prepared array backend"))
    compute_device(input) == compute_device(plan.diversity_opd) || throw(
        InvalidConfiguration("LiFT forward OPD must occupy the prepared compute device"))
    return input
end

function _prepare_lift_forward(plan::LiFTForwardPlan,
    input::AbstractMatrix)
    _require_lift_forward_input(plan, input)
    workspace = _allocate_lift_forward_workspace(plan.pupil_amplitude,
        size(plan.pupil_amplitude, 1), plan.focal_resolution,
        plan.zero_padding, plan.object_kernel, plan.mapping)
    output_values = similar(plan.pupil_amplitude,
        plan.observation_contract.rate_metadata.numeric_type,
        plan.observation_contract.rate_metadata.dimensions...)
    output = IntensityMap(plan.observation_contract.rate_metadata,
        output_values)
    forward = PreparedLiFTForward(plan, workspace, input, output,
        backend(input), compute_device(input))
    return _require_lift_forward_owner(forward)
end

function _prepare_lift_forward(plan::LiFTForwardPlan)
    workspace = _allocate_lift_forward_workspace(plan.pupil_amplitude,
        size(plan.pupil_amplitude, 1), plan.focal_resolution,
        plan.zero_padding, plan.object_kernel, plan.mapping)
    input = similar(plan.pupil_amplitude)
    copyto!(input, plan.diversity_opd)
    output_values = similar(plan.pupil_amplitude,
        plan.observation_contract.rate_metadata.numeric_type,
        plan.observation_contract.rate_metadata.dimensions...)
    output = IntensityMap(plan.observation_contract.rate_metadata,
        output_values)
    forward = PreparedLiFTForward(plan, workspace, input, output,
        backend(input), compute_device(input))
    return _require_lift_forward_owner(forward)
end

@inline function _lift_propagation_workspace_arrays(
    workspace::LiFTForwardWorkspace,
)
    propagation = workspace.propagation
    return (propagation.pupil_field, propagation.fft_buffer,
        propagation.psf_buffer)
end

@inline function _lift_forward_workspace_arrays(workspace::LiFTForwardWorkspace)
    return (workspace.optical_rate_buffer, workspace.amplitude_buffer,
        workspace.field_scratch, workspace.focal_buffer,
        workspace.mode_buffer, workspace.conjugate_field_buffer,
        workspace.response_buffer, workspace.response_scratch,
        workspace.sampling_buffer, workspace.mapped_rate_buffer,
        workspace.convolution_buffer, workspace.convolution_scratch)
end

function _require_lift_workspace_array(array::AbstractArray,
    ::Type{E}, dimensions::Tuple, plan::LiFTForwardPlan,
    label::AbstractString) where {E<:Number}
    size(array) == dimensions || throw(DimensionMismatchError(
        "$label dimensions do not match the LiFT forward plan"))
    eltype(array) === E || throw(InvalidConfiguration(
        "$label numeric type does not match the LiFT forward plan"))
    typeof(backend(array)) === typeof(backend(plan.pupil_amplitude)) || throw(
        InvalidConfiguration(
            "$label backend does not match the LiFT forward plan"))
    compute_device(array) == compute_device(plan.pupil_amplitude) || throw(
        InvalidConfiguration(
            "$label device does not match the LiFT forward plan"))
    return array
end

@inline function _require_lift_absent_workspace(value,
    label::AbstractString)
    value === nothing || throw(InvalidConfiguration(
        "$label must be absent for this LiFT forward plan"))
    return nothing
end

function _require_lift_mapping_workspace(::LiFTIdentityMapping,
    workspace::LiFTForwardWorkspace, ::LiFTForwardPlan,
    ::Type{<:AbstractFloat})
    _require_lift_absent_workspace(workspace.response_buffer,
        "LiFT response workspace")
    _require_lift_absent_workspace(workspace.response_scratch,
        "LiFT response scratch")
    _require_lift_absent_workspace(workspace.sampling_buffer,
        "LiFT sampling workspace")
    _require_lift_absent_workspace(workspace.mapped_rate_buffer,
        "LiFT mapped-rate workspace")
    return workspace
end

function _require_lift_mapping_workspace(mapping::LiFTFrameMapping,
    workspace::LiFTForwardWorkspace, plan::LiFTForwardPlan,
    ::Type{T}) where {T<:AbstractFloat}
    focal_dimensions = (plan.focal_resolution, plan.focal_resolution)
    sampled_resolution = div(plan.focal_resolution, mapping.sampling)
    sampled_dimensions = (sampled_resolution, sampled_resolution)
    output_dimensions = plan.observation_contract.rate_metadata.dimensions
    _require_lift_workspace_array(workspace.response_buffer, T,
        focal_dimensions, plan, "LiFT response workspace")
    _require_lift_workspace_array(workspace.response_scratch, T,
        focal_dimensions, plan, "LiFT response scratch")
    _require_lift_workspace_array(workspace.sampling_buffer, T,
        sampled_dimensions, plan, "LiFT sampling workspace")
    _require_lift_workspace_array(workspace.mapped_rate_buffer, T,
        output_dimensions, plan, "LiFT mapped-rate workspace")
    return workspace
end

function _require_lift_convolution_workspace(::Nothing,
    workspace::LiFTForwardWorkspace, ::LiFTForwardPlan,
    ::Type{<:AbstractFloat})
    _require_lift_absent_workspace(workspace.convolution_buffer,
        "LiFT convolution workspace")
    _require_lift_absent_workspace(workspace.convolution_scratch,
        "LiFT convolution scratch")
    return workspace
end

function _require_lift_convolution_workspace(::Union{
        LiFTDenseObjectKernel,LiFTSeparableObjectKernel},
    workspace::LiFTForwardWorkspace, plan::LiFTForwardPlan,
    ::Type{T}) where {T<:AbstractFloat}
    dimensions = (plan.focal_resolution, plan.focal_resolution)
    _require_lift_workspace_array(workspace.convolution_buffer, T,
        dimensions, plan, "LiFT convolution workspace")
    _require_lift_workspace_array(workspace.convolution_scratch, T,
        dimensions, plan, "LiFT convolution scratch")
    return workspace
end

function _require_lift_forward_workspace(plan::LiFTForwardPlan,
    workspace::LiFTForwardWorkspace)
    T = eltype(plan.pupil_amplitude)
    pupil_dimensions = size(plan.pupil_amplitude)
    padded_resolution = lift_pad_size(size(plan.pupil_amplitude, 1),
        plan.zero_padding)
    padded_dimensions = (padded_resolution, padded_resolution)
    focal_dimensions = (plan.focal_resolution, plan.focal_resolution)
    field_resolution = plan.focal_resolution *
        lift_oversampling(plan.zero_padding)
    field_dimensions = (field_resolution, field_resolution)
    propagation = workspace.propagation
    _require_lift_workspace_array(propagation.pupil_field, Complex{T},
        padded_dimensions, plan, "LiFT propagation pupil-field workspace")
    _require_lift_workspace_array(propagation.fft_buffer, Complex{T},
        padded_dimensions, plan, "LiFT propagation FFT workspace")
    _require_lift_workspace_array(propagation.psf_buffer, T,
        padded_dimensions, plan, "LiFT propagation intensity workspace")
    _require_lift_workspace_array(workspace.optical_rate_buffer, T,
        focal_dimensions, plan, "LiFT optical-rate workspace")
    _require_lift_workspace_array(workspace.amplitude_buffer, T,
        pupil_dimensions, plan, "LiFT amplitude workspace")
    _require_lift_workspace_array(workspace.field_scratch, T,
        field_dimensions, plan, "LiFT field-intensity scratch")
    _require_lift_workspace_array(workspace.focal_buffer, Complex{T},
        field_dimensions, plan, "LiFT focal-field workspace")
    _require_lift_workspace_array(workspace.mode_buffer, Complex{T},
        field_dimensions, plan, "LiFT mode-field workspace")
    _require_lift_workspace_array(workspace.conjugate_field_buffer,
        Complex{T}, field_dimensions, plan,
        "LiFT conjugate-field workspace")
    _require_lift_mapping_workspace(plan.mapping, workspace, plan, T)
    _require_lift_convolution_workspace(plan.object_kernel, workspace,
        plan, T)
    return workspace
end

@inline _lift_mightalias_any(::AbstractArray, ::Tuple{}) = false
@inline function _lift_mightalias_any(value::AbstractArray, values::Tuple)
    return _wfs_storage_mightalias(value, first(values)) ||
        _lift_mightalias_any(value, Base.tail(values))
end
@inline _lift_mightalias_any(::Any, ::Tuple) = false

@inline function _lift_mightalias_forward_workspace(
    value::AbstractArray,
    workspace::LiFTForwardWorkspace,
)
    return _lift_mightalias_any(
        value,
        _lift_propagation_workspace_arrays(workspace),
    ) || _lift_mightalias_any(value, _lift_forward_workspace_arrays(workspace))
end

function _require_lift_forward_owner(forward::PreparedLiFTForward)
    _require_lift_forward_input(forward.plan, forward.input)
    _require_lift_forward_workspace(forward.plan, forward.workspace)
    forward.output.metadata == forward.plan.observation_contract.rate_metadata ||
        throw(InvalidConfiguration(
            "LiFT forward output metadata does not match its prepared plan"))
    size(forward.output.values) ==
        forward.plan.observation_contract.rate_metadata.dimensions || throw(
        DimensionMismatchError(
            "LiFT forward output does not match its prepared dimensions"))
    eltype(forward.output.values) ===
        forward.plan.observation_contract.rate_metadata.numeric_type || throw(
        InvalidConfiguration(
            "LiFT forward output does not use its prepared numeric type"))
    typeof(backend(forward.output.values)) === typeof(forward.backend) || throw(
        InvalidConfiguration("LiFT forward output backend binding changed"))
    compute_device(forward.output.values) == forward.device || throw(
        InvalidConfiguration(
            "LiFT forward output compute-device binding changed"))
    typeof(forward.backend) === typeof(backend(forward.input)) || throw(
        InvalidConfiguration("LiFT forward backend binding changed"))
    forward.device == compute_device(forward.input) || throw(
        InvalidConfiguration("LiFT forward compute-device binding changed"))
    storages = (_lift_forward_plan_arrays(forward.plan)..., forward.input,
        forward.output.values,
        _lift_forward_workspace_arrays(forward.workspace)...)
    propagation_storages =
        _lift_propagation_workspace_arrays(forward.workspace)
    (_lift_any_alias(storages) ||
        _lift_any_alias(propagation_storages) ||
        _lift_any_cross_alias(propagation_storages, storages)) && throw(
        InvalidConfiguration(
            "LiFT forward plan, input, output, and workspace must not alias"))
    return forward
end

"""
    prepare_lift_forward_model(telescope, source, basis, opd; diversity_opd, ...)

Prepare the monochromatic LiFT focal-plane forward model independently of any
detector acquisition. `basis[:, :, k]` is a dimensionless modal OPD shape;
LiFT coefficients and `diversity_opd` are in metres, so their assembled sum is
an OPD map in metres. The telescope aperture and diversity are frozen into
backend-resident arrays; no mutable telescope or source object is retained.
"""
function prepare_lift_forward_model(tel::Telescope,
    src::Union{Source,LGSSource}, basis::AbstractArray{T,3},
    input::AbstractMatrix;
    diversity_opd::AbstractMatrix,
    focal_resolution::Int=0, zero_padding::Int=1, object_kernel=nothing,
    mapping::AbstractLiFTObservationMapping=LiFTIdentityMapping()) where {T<:AbstractFloat}
    _require_physical_photon_irradiance(src, "LiFT forward model")
    zero_padding >= 1 || throw(InvalidConfiguration(
        "LiFT zero_padding must be >= 1"))
    resolution = tel.params.resolution
    focal_resolution = focal_resolution <= 0 ? resolution * zero_padding :
        focal_resolution
    focal_resolution >= 1 || throw(InvalidConfiguration(
        "LiFT focal_resolution must be >= 1"))
    focal_resolution * lift_oversampling(zero_padding) <=
        lift_pad_size(resolution, zero_padding) || throw(DimensionMismatchError(
            "LiFT focal resolution exceeds the prepared padded focal field"))
    size(basis, 1) == resolution && size(basis, 2) == resolution || throw(
        DimensionMismatchError(
            "LiFT basis pupil dimensions must match telescope resolution"))
    size(diversity_opd) == (resolution, resolution) || throw(
        DimensionMismatchError(
            "LiFT diversity OPD must match telescope resolution"))
    prototype = pupil_reflectivity(tel)
    eltype(prototype) === T || throw(InvalidConfiguration(
        "LiFT telescope and basis must use the same numeric type"))
    typeof(backend(basis)) === typeof(backend(tel)) || throw(
        InvalidConfiguration(
            "LiFT telescope and basis must use the same array backend"))
    compute_device(basis) == compute_device(prototype) || throw(
        InvalidConfiguration(
            "LiFT telescope and basis must occupy the same compute device"))
    typeof(backend(diversity_opd)) === typeof(backend(tel)) || throw(
        InvalidConfiguration(
            "LiFT diversity and telescope must use the same array backend"))
    compute_device(diversity_opd) == compute_device(prototype) || throw(
        InvalidConfiguration(
            "LiFT diversity and telescope must occupy the same compute device"))
    size(input) == (resolution, resolution) || throw(DimensionMismatchError(
        "LiFT forward OPD must match telescope resolution"))
    eltype(input) === T || throw(InvalidConfiguration(
        "LiFT forward OPD must use the basis numeric type"))
    typeof(backend(input)) === typeof(backend(tel)) || throw(
        InvalidConfiguration(
            "LiFT forward OPD and telescope must use the same array backend"))
    compute_device(input) == compute_device(prototype) || throw(
        InvalidConfiguration(
            "LiFT forward OPD and telescope must occupy the same compute device"))
    _require_lift_mapping_backend(mapping, prototype)
    prepared_mapping = _prepare_lift_mapping(mapping)

    pupil = copy(pupil_mask(tel))
    reflectivity = pupil_reflectivity(tel)
    amplitude = similar(prototype, T, resolution, resolution)
    @. amplitude = sqrt(reflectivity)
    owned_basis = copy(basis)
    diversity = _copy_lift_array(prototype, diversity_opd, T)
    kernel = _prepare_lift_object_kernel(object_kernel, prototype, T)
    output_dimensions = _lift_output_dimensions(focal_resolution,
        prepared_mapping)
    output_values = similar(prototype, T, output_dimensions...)
    workspace = _allocate_lift_forward_workspace(prototype, resolution,
        focal_resolution, zero_padding, kernel, prepared_mapping)
    contract = LiFTObservationContract(
        OpticalPlaneMetadata(FocalPlane(), output_values;
            coordinate_domain=AngularCoordinates(),
            sampling=(T(wavelength(src) /
                (tel.params.diameter * zero_padding) *
                prod(_lift_mapping_factors(prepared_mapping))),
                T(wavelength(src) /
                (tel.params.diameter * zero_padding) *
                prod(_lift_mapping_factors(prepared_mapping)))),
            spectral=MonochromaticChannel(T(wavelength(src))),
            normalization=PhotonRateNormalization(),
            spatial_measure=CellIntegratedMeasure(),
            coherence=IncoherentIntensityAddition()),
        _lift_mapping_signature(prepared_mapping))
    plan = LiFTForwardPlan(pupil, amplitude, owned_basis, diversity,
        T(wavelength(src)), T(photon_irradiance(src)),
        T((tel.params.diameter / resolution)^2), focal_resolution,
        zero_padding, kernel, prepared_mapping, contract)
    output = IntensityMap(contract.rate_metadata, output_values)
    forward = PreparedLiFTForward(plan, workspace, input, output,
        backend(input), compute_device(input))
    return _require_lift_forward_owner(forward)
end

"""Return the shareable run-immutable plan from a prepared LiFT forward owner."""
@inline lift_forward_plan(forward::PreparedLiFTForward) = forward.plan

"""Return the replaceable workspace from a prepared LiFT forward owner."""
@inline lift_forward_workspace(forward::PreparedLiFTForward) =
    forward.workspace

"""Return the immutable observation-compatibility contract for `forward`."""
@inline lift_observation_contract(forward::PreparedLiFTForward) =
    forward.plan.observation_contract
"""Return the caller-visible photon-arrival-rate output owned by `forward`."""
@inline lift_forward_output(forward::PreparedLiFTForward) = forward.output

function LiFTObservation(contract::LiFTObservationContract,
    values::AbstractMatrix{E};
    domain::AbstractLiFTObservationDomain=LiFTPhotonRate(),
    readout_noise_std::Real=0, validate_values::Bool=true) where {E<:Real}
    rate_metadata = contract.rate_metadata
    size(values) == rate_metadata.dimensions || throw(DimensionMismatchError(
        "LiFT observation dimensions do not match its prepared contract"))
    typeof(backend(values)) === typeof(rate_metadata.backend) || throw(
        InvalidConfiguration(
            "LiFT observation and forward model must use the same array backend"))
    device = compute_device(values)
    device == rate_metadata.device || throw(InvalidConfiguration(
        "LiFT observation and forward model must occupy the same compute device"))
    sigma = rate_metadata.numeric_type(readout_noise_std)
    isfinite(sigma) && sigma >= zero(sigma) || throw(InvalidConfiguration(
        "LiFT observation readout noise must be finite and nonnegative"))
    validate_values && _require_finite_nonnegative_intensity(values)
    metadata = LiFTObservationMetadata(contract, domain, sigma, E,
        backend(values), device)
    return LiFTObservation(metadata, values)
end

LiFTObservation(forward::PreparedLiFTForward,
    values::AbstractMatrix; kwargs...) =
    LiFTObservation(lift_observation_contract(forward), values; kwargs...)


@inline _lift_object_plan_arrays(::Nothing) = ()
@inline _lift_object_plan_arrays(kernel::LiFTDenseObjectKernel) =
    (kernel.kernel,)
@inline _lift_object_plan_arrays(kernel::LiFTSeparableObjectKernel) =
    (kernel.row, kernel.col)

@inline _lift_response_plan_arrays(::NullFrameResponse) = ()
@inline _lift_response_plan_arrays(response::GaussianPixelResponse) =
    (response.kernel,)
@inline _lift_response_plan_arrays(response::SampledFrameResponse) =
    (response.kernel,)
@inline _lift_response_plan_arrays(response::RectangularPixelAperture) =
    (response.kernel_x, response.kernel_y)

@inline _lift_mapping_plan_arrays(::LiFTIdentityMapping) = ()
@inline _lift_mapping_plan_arrays(mapping::LiFTFrameMapping) =
    _lift_response_plan_arrays(mapping.response)

@inline function _lift_forward_plan_arrays(plan::LiFTForwardPlan)
    return (plan.pupil_mask, plan.pupil_amplitude, plan.basis,
        plan.diversity_opd, _lift_object_plan_arrays(plan.object_kernel)...,
        _lift_mapping_plan_arrays(plan.mapping)...)
end

@inline _lift_any_alias(::Tuple{}) = false
@inline function _lift_any_alias(values::Tuple)
    remaining = Base.tail(values)
    return _lift_mightalias_any(first(values), remaining) ||
        _lift_any_alias(remaining)
end

@inline _lift_any_cross_alias(::Tuple{}, ::Tuple) = false
@inline function _lift_any_cross_alias(values::Tuple, other_values::Tuple)
    return _lift_mightalias_any(first(values), other_values) ||
        _lift_any_cross_alias(Base.tail(values), other_values)
end
