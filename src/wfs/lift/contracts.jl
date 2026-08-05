abstract type LiFTJacobianMethod end
struct LiFTAnalyticJacobian <: LiFTJacobianMethod end
struct LiFTNumericalJacobian <: LiFTJacobianMethod end

abstract type LiFTSolveMode end
struct LiFTSolveAuto <: LiFTSolveMode end
struct LiFTSolveQR <: LiFTSolveMode end
struct LiFTSolveNormalEquations <: LiFTSolveMode end

abstract type LiFTDampingMode end
struct LiFTDampingNone <: LiFTDampingMode end
struct LiFTLevenbergMarquardt{T<:AbstractFloat} <: LiFTDampingMode
    lambda0::T
    growth::T
    condition_rtol::T
end

struct LiFTAdaptiveLevenbergMarquardt{T<:AbstractFloat} <: LiFTDampingMode
    lambda0::T
    growth::T
    shrink::T
    min_lambda::T
    condition_rtol::T
end

function LiFTLevenbergMarquardt(lambda0::Real, growth::Real, condition_rtol::Real)
    lambda, promoted_growth, promoted_rtol = promote(
        float(lambda0), float(growth), float(condition_rtol))
    return LiFTLevenbergMarquardt{typeof(lambda)}(
        lambda, promoted_growth, promoted_rtol)
end

function LiFTAdaptiveLevenbergMarquardt(lambda0::Real, growth::Real, shrink::Real,
    min_lambda::Real, condition_rtol::Real)
    lambda, promoted_growth, promoted_shrink, promoted_min_lambda, promoted_rtol = promote(
        float(lambda0), float(growth), float(shrink), float(min_lambda), float(condition_rtol))
    return LiFTAdaptiveLevenbergMarquardt{typeof(lambda)}(
        lambda, promoted_growth, promoted_shrink, promoted_min_lambda, promoted_rtol)
end

LiFTLevenbergMarquardt(; lambda0::Real=1e-6, growth::Real=10.0, condition_rtol::Real=sqrt(eps(Float64))) =
    LiFTLevenbergMarquardt(float(lambda0), float(growth), float(condition_rtol))

LiFTAdaptiveLevenbergMarquardt(; lambda0::Real=1e-6, growth::Real=10.0, shrink::Real=2.0,
    min_lambda::Real=1e-10, condition_rtol::Real=sqrt(eps(Float64))) =
    LiFTAdaptiveLevenbergMarquardt(float(lambda0), float(growth), float(shrink), float(min_lambda),
        float(condition_rtol))

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
Photon-arrival-rate observations. The noise-equivalent exposure is used only
when model-based Poisson weighting is requested; it does not scale the values.
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
unit. The noise-equivalent exposure is used only for model-based weighting.
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

@inline lift_shot_variance_rate_scale(domain::LiFTPhotonRate,
    ::Type{T}) where {T<:AbstractFloat} =
    inv(T(domain.noise_equivalent_exposure_s) * T(domain.quantum_efficiency))
@inline lift_shot_variance_rate_scale(domain::LiFTExpectedCounts,
    ::Type{T}) where {T<:AbstractFloat} =
    inv(T(domain.exposure_duration_s) * T(domain.quantum_efficiency))
@inline lift_shot_variance_rate_scale(domain::LiFTNormalizedIntensity,
    ::Type{T}) where {T<:AbstractFloat} =
    inv(T(domain.noise_equivalent_exposure_s) * T(domain.quantum_efficiency))

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

abstract type LiFTWeightingMode end
struct LiFTInitialModelWeighting <: LiFTWeightingMode end
struct LiFTIterativeModelWeighting <: LiFTWeightingMode end
struct LiFTReadNoiseWeighting <: LiFTWeightingMode end
struct LiFTVarianceMapWeighting{M<:AbstractMatrix} <: LiFTWeightingMode
    variance::M
end

abstract type LiFTFluxNormalization end
struct LiFTTotalFluxNormalization <: LiFTFluxNormalization end
struct LiFTPeakIntensityNormalization <: LiFTFluxNormalization end
struct LiFTFixedFlux <: LiFTFluxNormalization end

"""Cold configuration for one LiFT phase-retrieval estimator."""
struct LiFT{J<:LiFTJacobianMethod,S<:LiFTSolveMode,D<:LiFTDampingMode,
    I,W<:LiFTWeightingMode,N<:LiFTFluxNormalization}
    iterations::Int
    jacobian_method::J
    solve_mode::S
    damping::D
    mode_ids::I
    weighting::W
    flux_normalization::N
    check_convergence::Bool
end

function LiFT(; iterations::Int=5,
    jacobian_method::LiFTJacobianMethod=LiFTAnalyticJacobian(),
    solve_mode::LiFTSolveMode=LiFTSolveAuto(),
    damping::LiFTDampingMode=LiFTDampingNone(), mode_ids=nothing,
    weighting::LiFTWeightingMode=LiFTReadNoiseWeighting(),
    flux_normalization::LiFTFluxNormalization=LiFTTotalFluxNormalization(),
    check_convergence::Bool=true)
    iterations >= 1 || throw(InvalidConfiguration(
        "LiFT iterations must be >= 1"))
    prepared_mode_ids = if mode_ids === nothing
        nothing
    else
        values = FixedSizeVectorDefault{Int}(collect(Int, mode_ids))
        isempty(values) && throw(InvalidConfiguration(
            "LiFT mode_ids must not be empty"))
        all(>(0), values) || throw(InvalidConfiguration(
            "LiFT mode_ids must be positive"))
        allunique(values) || throw(InvalidConfiguration(
            "LiFT mode_ids must be unique"))
        values
    end
    return LiFT{typeof(jacobian_method),typeof(solve_mode),typeof(damping),
        typeof(prepared_mode_ids),typeof(weighting),typeof(flux_normalization)}(
        iterations, jacobian_method, solve_mode, damping, prepared_mode_ids,
        weighting, flux_normalization, check_convergence)
end

struct LiFTEstimationPlan{J<:LiFTJacobianMethod,S<:LiFTSolveMode,
    D<:LiFTDampingMode,I<:FixedSizeVector,W<:LiFTWeightingMode,
    N<:LiFTFluxNormalization,C<:LiFTObservationContract}
    iterations::Int
    jacobian_method::J
    solve_mode::S
    damping::D
    mode_ids::I
    weighting::W
    flux_normalization::N
    check_convergence::Bool
    observation_contract::C
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

struct LiFTDiagnostics{T<:AbstractFloat}
    residual_norm::T
    weighted_residual_norm::T
    update_norm::T
    condition_ratio::T
    regularization::T
    used_qr::Bool
    used_fallback::Bool
end

mutable struct LiFTIterationWorkspace{T<:AbstractFloat}
    residual_norm::T
    weighted_residual_norm::T
    update_norm::T
    condition_ratio::T
    regularization::T
    used_qr::Bool
    used_fallback::Bool
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
    output_work_buffer::B
    convolution_buffer::CB
    convolution_scratch::CB
    opd_work_buffer::B
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

struct LiFTEstimationWorkspace{T<:AbstractFloat,
    B<:AbstractMatrix{T},
    V<:AbstractVector{T},
    I<:AbstractVector{Int}}
    observation_rate_buffer::B
    residual_buffer::V
    weight_buffer::V
    H_buffer::B
    normal_buffer::B
    factor_buffer::B
    rhs_buffer::V
    full_coefficients_buffer::V
    mode_id_buffer::I
    iteration::LiFTIterationWorkspace{T}
end

struct PreparedLiFTEstimator{F<:PreparedLiFTForward,
    P<:LiFTEstimationPlan,S<:LiFTEstimationWorkspace,
    O<:LiFTObservation,C<:AbstractVector,I,B,D}
    forward::F
    plan::P
    workspace::S
    observation::O
    coefficients::C
    initial_coefficients::I
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

@inline _lift_output_values(::LiFTIdentityMapping,
    workspace::LiFTForwardWorkspace) = workspace.optical_rate_buffer
@inline _lift_output_values(::LiFTFrameMapping,
    workspace::LiFTForwardWorkspace) = workspace.mapped_rate_buffer

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
    output_dimensions = _lift_output_dimensions(focal_resolution, mapping)
    output_work = similar(template, T, output_dimensions...)
    convolution = object_kernel === nothing ? nothing : similar(optical_rate)
    convolution_scratch = object_kernel === nothing ? nothing : similar(optical_rate)
    opd_work = similar(pupil_amplitude)
    return LiFTForwardWorkspace(propagation, optical_rate,
        pupil_amplitude, field_scratch, focal, mode, conjugate_field,
        response, response_scratch, sampled, mapped, output_work, convolution,
        convolution_scratch, opd_work)
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
    output_values = _lift_output_values(plan.mapping, workspace)
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
    output_values = _lift_output_values(plan.mapping, workspace)
    output = IntensityMap(plan.observation_contract.rate_metadata,
        output_values)
    forward = PreparedLiFTForward(plan, workspace, input, output,
        backend(input), compute_device(input))
    return _require_lift_forward_owner(forward)
end

@inline function _lift_forward_workspace_arrays(workspace::LiFTForwardWorkspace)
    return (workspace.optical_rate_buffer, workspace.amplitude_buffer,
        workspace.field_scratch, workspace.focal_buffer,
        workspace.mode_buffer, workspace.conjugate_field_buffer,
        workspace.response_buffer, workspace.response_scratch,
        workspace.sampling_buffer, workspace.mapped_rate_buffer,
        workspace.output_work_buffer, workspace.convolution_buffer,
        workspace.convolution_scratch, workspace.opd_work_buffer)
end

@inline _lift_mightalias_any(::AbstractArray, ::Tuple{}) = false
@inline function _lift_mightalias_any(value::AbstractArray, values::Tuple)
    return _wfs_storage_mightalias(value, first(values)) ||
        _lift_mightalias_any(value, Base.tail(values))
end

function _require_lift_forward_owner(forward::PreparedLiFTForward)
    _require_lift_forward_input(forward.plan, forward.input)
    forward.output.metadata == forward.plan.observation_contract.rate_metadata ||
        throw(InvalidConfiguration(
            "LiFT forward output metadata does not match its prepared plan"))
    forward.output.values === _lift_output_values(
        forward.plan.mapping, forward.workspace) || throw(
        InvalidConfiguration(
            "LiFT forward output does not match its prepared workspace"))
    typeof(forward.backend) === typeof(backend(forward.input)) || throw(
        InvalidConfiguration("LiFT forward backend binding changed"))
    forward.device == compute_device(forward.input) || throw(
        InvalidConfiguration("LiFT forward compute-device binding changed"))
    _lift_mightalias_any(forward.input,
        _lift_forward_workspace_arrays(forward.workspace)) && throw(
        InvalidConfiguration(
            "LiFT forward input must not alias its output or workspace"))
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
    _lift_output_dimensions(focal_resolution, prepared_mapping)

    pupil = copy(pupil_mask(tel))
    reflectivity = pupil_reflectivity(tel)
    amplitude = similar(prototype, T, resolution, resolution)
    @. amplitude = sqrt(reflectivity)
    owned_basis = copy(basis)
    diversity = _copy_lift_array(prototype, diversity_opd, T)
    kernel = _prepare_lift_object_kernel(object_kernel, prototype, T)
    workspace = _allocate_lift_forward_workspace(prototype, resolution,
        focal_resolution, zero_padding, kernel, prepared_mapping)
    output_values = _lift_output_values(prepared_mapping, workspace)
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

"""
    prepare_lift_estimator(definition, forward, observation, coefficients;
        initial_coefficients=nothing)

Prepare an exact single-writer LiFT estimator. `coefficients` is the selected
caller-owned result. When supplied, `initial_coefficients` contains one value
for every mode in the forward plan and is copied into iteration scratch at the
start of each reconstruction.
"""
function prepare_lift_estimator(definition::LiFT,
    forward::PreparedLiFTForward, observation::LiFTObservation,
    coefficients::AbstractVector; initial_coefficients=nothing)
    _require_lift_forward_owner(forward)
    plan = forward.plan
    output = intensity_values(forward.output)
    T = eltype(output)
    prepared_mode_ids = definition.mode_ids === nothing ?
        FixedSizeVectorDefault{Int}(axes(plan.basis, 3)) :
        FixedSizeVectorDefault{Int}(definition.mode_ids)
    all(mode_id -> 1 <= mode_id <= size(plan.basis, 3),
        prepared_mode_ids) || throw(DimensionMismatchError(
            "LiFT mode ids must index the prepared modal basis"))
    observation.metadata.contract == plan.observation_contract || throw(
        InvalidConfiguration(
            "LiFT observation geometry, wavelength, or preprocessing does not match the prepared forward plan"))
    mode_count = length(prepared_mode_ids)
    length(coefficients) == mode_count || throw(DimensionMismatchError(
        "LiFT coefficient product must match the prepared modal subset"))
    eltype(coefficients) === T || throw(InvalidConfiguration(
        "LiFT coefficient product must use the prepared numeric type"))
    typeof(backend(coefficients)) === typeof(backend(output)) || throw(
        InvalidConfiguration(
            "LiFT coefficient product must use the prepared array backend"))
    compute_device(coefficients) == compute_device(output) || throw(
        InvalidConfiguration(
            "LiFT coefficient product must occupy the prepared compute device"))
    _wfs_storage_mightalias(observation.values, coefficients) && throw(
        InvalidConfiguration(
            "LiFT observation and coefficient product must not alias"))
    _require_lift_initial_coefficients(initial_coefficients, output,
        size(plan.basis, 3))
    _require_lift_weighting(definition.weighting, observation.metadata, T)
    prepared_weighting = _prepare_lift_weighting(definition.weighting, output)

    estimator_forward = _prepare_lift_forward(plan)
    observation_rate = similar(output)
    residual = similar(output, T, length(output))
    weights = similar(residual)
    H = similar(output, T, length(output), mode_count)
    normal = similar(output, T, mode_count, mode_count)
    factor = similar(normal)
    rhs = similar(residual, mode_count)
    full_coefficients = similar(residual, T, size(plan.basis, 3))
    mode_ids_device = similar(rhs, Int, mode_count)
    copyto!(mode_ids_device, collect(prepared_mode_ids))
    iteration = LiFTIterationWorkspace(T(NaN), T(NaN), T(NaN), T(NaN), zero(T),
        false, false)
    workspace = LiFTEstimationWorkspace(observation_rate, residual, weights,
        H, normal, factor, rhs, full_coefficients, mode_ids_device, iteration)
    estimation_plan = LiFTEstimationPlan(definition.iterations,
        definition.jacobian_method, definition.solve_mode, definition.damping,
        prepared_mode_ids, prepared_weighting, definition.flux_normalization,
        definition.check_convergence, plan.observation_contract)
    estimator = PreparedLiFTEstimator(estimator_forward, estimation_plan,
        workspace, observation, coefficients, initial_coefficients,
        backend(output),
        compute_device(output))
    _require_lift_estimation_aliases(estimator)
    return _require_lift_estimator(estimator)
end

@inline _prepare_lift_weighting(mode::LiFTWeightingMode,
    ::AbstractMatrix) = mode

function _prepare_lift_weighting(mode::LiFTVarianceMapWeighting,
    output::AbstractMatrix{T}) where {T<:AbstractFloat}
    return LiFTVarianceMapWeighting(_copy_lift_array(output, mode.variance, T))
end

@inline _require_lift_initial_coefficients(::Nothing,
    ::AbstractMatrix, ::Int) = nothing

function _require_lift_initial_coefficients(initial::AbstractVector,
    output::AbstractMatrix{T}, mode_count::Int) where {T<:AbstractFloat}
    length(initial) == mode_count || throw(DimensionMismatchError(
        "LiFT initial coefficients must cover the complete prepared basis"))
    eltype(initial) === T || throw(InvalidConfiguration(
        "LiFT initial coefficients must use the prepared numeric type"))
    typeof(backend(initial)) === typeof(backend(output)) || throw(
        InvalidConfiguration(
            "LiFT initial coefficients must use the prepared array backend"))
    compute_device(initial) == compute_device(output) || throw(
        InvalidConfiguration(
            "LiFT initial coefficients must occupy the prepared compute device"))
    return nothing
end

_require_lift_initial_coefficients(::Any, ::AbstractMatrix, ::Int) = throw(
    InvalidConfiguration(
        "LiFT initial coefficients must be an abstract vector or nothing"))

@inline function _lift_estimation_workspace_arrays(
    workspace::LiFTEstimationWorkspace)
    return (workspace.observation_rate_buffer, workspace.residual_buffer,
        workspace.weight_buffer, workspace.H_buffer,
        workspace.normal_buffer, workspace.factor_buffer,
        workspace.rhs_buffer, workspace.full_coefficients_buffer,
        workspace.mode_id_buffer)
end

function _require_lift_estimation_aliases(lift::PreparedLiFTEstimator)
    forward_arrays = _lift_forward_workspace_arrays(lift.forward.workspace)
    estimation_arrays = _lift_estimation_workspace_arrays(lift.workspace)
    resources = (forward_arrays..., estimation_arrays...)
    (_lift_mightalias_any(lift.observation.values, resources) ||
        _lift_mightalias_any(lift.coefficients, resources) ||
        _lift_initial_mightalias_any(lift.initial_coefficients, resources) ||
        _wfs_storage_mightalias(
            lift.observation.values, lift.coefficients) ||
        _wfs_storage_mightalias(
            lift.observation.values, lift.initial_coefficients) ||
        _wfs_storage_mightalias(
            lift.coefficients, lift.initial_coefficients)) && throw(
        InvalidConfiguration(
            "LiFT observation, coefficient product, and workspaces must not alias"))
    return lift
end


@inline _lift_initial_mightalias_any(::Nothing, ::Tuple) = false
@inline _lift_initial_mightalias_any(initial::AbstractArray,
    resources::Tuple) = _lift_mightalias_any(initial, resources)
@inline _lift_initial_mightalias_any(::Any, ::Tuple) = false

"""Return the diagnostics from the most recent LiFT reconstruction."""
function diagnostics(lift::PreparedLiFTEstimator)
    iteration = lift.workspace.iteration
    return LiFTDiagnostics(iteration.residual_norm,
        iteration.weighted_residual_norm, iteration.update_norm,
        iteration.condition_ratio, iteration.regularization,
        iteration.used_qr, iteration.used_fallback)
end

"""Return the run-immutable plan from a prepared LiFT estimator."""
@inline lift_estimation_plan(lift::PreparedLiFTEstimator) = lift.plan

"""Return the replaceable workspace from a prepared LiFT estimator."""
@inline lift_estimation_workspace(lift::PreparedLiFTEstimator) =
    lift.workspace

@inline _lift_model(lift::PreparedLiFTEstimator) = lift.forward.plan
@inline _lift_workspace(lift::PreparedLiFTEstimator) = lift.forward.workspace
