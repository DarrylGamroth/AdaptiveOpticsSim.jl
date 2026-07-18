#
# LiFT phase retrieval
#
# LiFT fits modal coefficients by matching a separately prepared focal-plane
# forward model to a caller-owned observation. Optical formation never owns or
# triggers a detector; acquisition timing, QE, and stochastic readout remain at
# the detector boundary.
#
# Forward model:
# 1. combine modal coefficients into an OPD map
# 2. add the configured diversity term
# 3. propagate to the focal plane and form intensity
# 4. optionally convolve with an object kernel
#
# Inverse model:
# - `LiFTAnalytic` builds the Jacobian from focal-plane field derivatives
# - `LiFTNumerical` builds the Jacobian by centered finite differences
#
# The update step is then solved by QR or normal equations with explicit
# damping/fallback logic.
#
@kernel function lift_scatter_update_kernel!(coeffs, delta, mode_ids, n_modes::Int)
    i = @index(Global, Linear)
    if i <= n_modes
        @inbounds coeffs[mode_ids[i]] += delta[i]
    end
end

@kernel function lift_gather_kernel!(out, coeffs, mode_ids, n_modes::Int)
    i = @index(Global, Linear)
    if i <= n_modes
        @inbounds out[i] = coeffs[mode_ids[i]]
    end
end

@kernel function lift_sqrt_weights_kernel!(weights, n::Int)
    i = @index(Global, Linear)
    if i <= n
        @inbounds weights[i] = sqrt(weights[i])
    end
end

@kernel function lift_affine_basis_mode_kernel!(dest, base, basis,
    scale, mode_offset::Int, n::Int)
    i = @index(Global, Linear)
    if i <= n
        @inbounds dest[i] = base[i] + scale * basis[i + mode_offset]
    end
end

@kernel function lift_scaled_basis_mode_kernel!(dest, amplitude, basis,
    scale, mode_offset::Int, n::Int)
    i = @index(Global, Linear)
    if i <= n
        @inbounds dest[i] = amplitude[i] * scale * basis[i + mode_offset]
    end
end

@kernel function lift_copy_column_kernel!(dest, column::Int, src, n::Int)
    i = @index(Global, Linear)
    if i <= n
        @inbounds dest[i, column] = src[i]
    end
end

@kernel function lift_residual_kernel!(dest, observation, model, n::Int)
    i = @index(Global, Linear)
    if i <= n
        @inbounds dest[i] = observation[i] - model[i]
    end
end

@kernel function lift_row_weights_kernel!(matrix, weights,
    n_rows::Int, n_cols::Int)
    i, j = @index(Global, NTuple)
    if i <= n_rows && j <= n_cols
        @inbounds matrix[i, j] *= weights[i]
    end
end

@kernel function lift_inverse_variance_kernel!(dest, values,
    scale, offset, floor_value, n::Int)
    i = @index(Global, Linear)
    if i <= n
        @inbounds dest[i] = inv(max(values[i] * scale + offset,
            floor_value))
    end
end

@kernel function lift_add_diagonal_kernel!(matrix, value, n::Int)
    i = @index(Global, Linear)
    if i <= n
        @inbounds matrix[i, i] += value
    end
end

@kernel function lift_dense_convolution_kernel!(dest, src, kernel,
    inv_norm, n::Int, m::Int, kh::Int, kw::Int)
    i, j = @index(Global, NTuple)
    if i <= n && j <= m
        T = eltype(dest)
        acc = zero(T)
        cx = div(kh, 2)
        cy = div(kw, 2)
        @inbounds for ki in 1:kh, kj in 1:kw
            ii = symm_index(i + ki - cx - 1, n)
            jj = symm_index(j + kj - cy - 1, m)
            acc += src[ii, jj] * kernel[ki, kj]
        end
        @inbounds dest[i, j] = acc * inv_norm
    end
end

@kernel function lift_row_convolution_kernel!(dest, src, kernel,
    n::Int, m::Int, nk::Int)
    i, j = @index(Global, NTuple)
    if i <= n && j <= m
        T = eltype(dest)
        acc = zero(T)
        center = div(nk, 2)
        @inbounds for k in 1:nk
            ii = symm_index(i + k - center - 1, n)
            acc += src[ii, j] * kernel[k]
        end
        @inbounds dest[i, j] = acc
    end
end

@kernel function lift_column_convolution_kernel!(dest, src, kernel,
    inv_norm, n::Int, m::Int, nk::Int)
    i, j = @index(Global, NTuple)
    if i <= n && j <= m
        T = eltype(dest)
        acc = zero(T)
        center = div(nk, 2)
        @inbounds for k in 1:nk
            jj = symm_index(j + k - center - 1, m)
            acc += src[i, jj] * kernel[k]
        end
        @inbounds dest[i, j] = acc * inv_norm
    end
end

abstract type LiFTMode end
struct LiFTAnalytic <: LiFTMode end
struct LiFTNumerical <: LiFTMode end

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
    exposure_time_s::T
    quantum_efficiency::T
end

function LiFTExpectedCounts(exposure_time_s::Real;
    quantum_efficiency::Real=1.0)
    exposure, qe = promote(float(exposure_time_s), float(quantum_efficiency))
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
    inv(T(domain.exposure_time_s) * T(domain.quantum_efficiency))
@inline lift_observation_to_rate_scale(domain::LiFTNormalizedIntensity,
    ::Type{T}) where {T<:AbstractFloat} = T(domain.photon_rate_per_unit)

@inline lift_shot_variance_rate_scale(domain::LiFTPhotonRate,
    ::Type{T}) where {T<:AbstractFloat} =
    inv(T(domain.noise_equivalent_exposure_s) * T(domain.quantum_efficiency))
@inline lift_shot_variance_rate_scale(domain::LiFTExpectedCounts,
    ::Type{T}) where {T<:AbstractFloat} =
    inv(T(domain.exposure_time_s) * T(domain.quantum_efficiency))
@inline lift_shot_variance_rate_scale(domain::LiFTNormalizedIntensity,
    ::Type{T}) where {T<:AbstractFloat} =
    inv(T(domain.noise_equivalent_exposure_s) * T(domain.quantum_efficiency))

struct LiFTObservationContract{M<:OpticalPlaneMetadata,S}
    rate_metadata::M
    preprocessing_signature::S
end

struct LiFTObservationMetadata{T<:AbstractFloat,
    C<:LiFTObservationContract,D<:AbstractLiFTObservationDomain,E,
    B<:AbstractArrayBackend,PD<:AbstractPlaneDevice}
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

struct LiFTParams{S<:LiFTSolveMode,D<:LiFTDampingMode,I<:Tuple}
    iterations::Int
    solve_mode::S
    damping::D
    mode_ids::I
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

mutable struct LiFTDiagnostics{T<:AbstractFloat}
    residual_norm::T
    weighted_residual_norm::T
    update_norm::T
    condition_ratio::T
    regularization::T
    used_qr::Bool
    used_fallback::Bool
end

struct LiFTForwardModel{T<:AbstractFloat,
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
    opd_buffer::B
    opd_work_buffer::B
end

"""Run-immutable LiFT optical definition plus a single-writer workspace."""
struct PreparedLiFTForwardModel{M<:LiFTForwardModel,
    W<:LiFTForwardWorkspace,O<:IntensityMap}
    model::M
    workspace::W
    output::O
end

struct LiFTState{T<:AbstractFloat,
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
    mode_id_buffer::I
    diagnostics::LiFTDiagnostics{T}
end

struct LiFT{M<:LiFTMode,F<:PreparedLiFTForwardModel,
    P<:LiFTParams,S<:LiFTState}
    forward::F
    params::P
    state::S
end

abstract type LiFTWeightingMode end
abstract type LiFTWeightingStatic <: LiFTWeightingMode end
abstract type LiFTWeightingDynamic <: LiFTWeightingMode end

struct LiFTWeightModel <: LiFTWeightingDynamic end
struct LiFTWeightIterative <: LiFTWeightingDynamic end
struct LiFTWeightNone <: LiFTWeightingStatic end
struct LiFTWeightMatrix{M<:AbstractMatrix} <: LiFTWeightingStatic
    R_n::M
end

weight_mode(mode::LiFTWeightingMode) = mode
weight_mode(::Nothing) = LiFTWeightNone()
weight_mode(R_n::AbstractMatrix) = LiFTWeightMatrix(R_n)
function weight_mode(R_n::Symbol)
    if R_n === :model
        return LiFTWeightModel()
    elseif R_n === :iterative
        return LiFTWeightIterative()
    end
    throw(InvalidConfiguration("R_n must be :model, :iterative, a matrix, or nothing"))
end
weight_mode(::Any) = throw(InvalidConfiguration("R_n must be :model, :iterative, a matrix, or nothing"))

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
    plane_device(array) == plane_device(template) || throw(
        InvalidConfiguration(
            "LiFT response and forward model must occupy the same physical device"))
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

"""
    prepare_lift_forward_model(telescope, source, basis; diversity_opd, ...)

Prepare the monochromatic LiFT focal-plane forward model independently of any
detector acquisition. `basis[:, :, k]` is a dimensionless modal OPD shape;
LiFT coefficients and `diversity_opd` are in metres, so their assembled sum is
an OPD map in metres. The telescope aperture and diversity are frozen into
backend-resident arrays; no mutable telescope or source object is retained.
"""
function prepare_lift_forward_model(tel::Telescope,
    src::Union{Source,LGSSource}, basis::AbstractArray{T,3};
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
    eltype(opd_map(tel)) === T || throw(InvalidConfiguration(
        "LiFT telescope and basis must use the same numeric type"))
    typeof(backend(basis)) === typeof(backend(tel)) || throw(
        InvalidConfiguration(
            "LiFT telescope and basis must use the same array backend"))
    plane_device(basis) == plane_device(opd_map(tel)) || throw(
        InvalidConfiguration(
            "LiFT telescope and basis must occupy the same physical device"))
    typeof(backend(diversity_opd)) === typeof(backend(tel)) || throw(
        InvalidConfiguration(
            "LiFT diversity and telescope must use the same array backend"))
    plane_device(diversity_opd) == plane_device(opd_map(tel)) || throw(
        InvalidConfiguration(
            "LiFT diversity and telescope must occupy the same physical device"))
    _lift_output_dimensions(focal_resolution, mapping)
    _require_lift_mapping_backend(mapping, opd_map(tel))

    pupil = copy(pupil_mask(tel))
    reflectivity = pupil_reflectivity(tel)
    amplitude = similar(opd_map(tel), T, resolution, resolution)
    @. amplitude = sqrt(reflectivity)
    owned_basis = copy(basis)
    diversity = _copy_lift_array(opd_map(tel), diversity_opd, T)
    kernel = _prepare_lift_object_kernel(object_kernel, opd_map(tel), T)
    oversampling = lift_oversampling(zero_padding)
    propagation = Workspace(opd_map(tel),
        lift_pad_size(resolution, zero_padding); T=T)
    optical_rate = similar(opd_map(tel), T, focal_resolution,
        focal_resolution)
    pupil_amplitude = similar(opd_map(tel), T, resolution, resolution)
    focal_size = focal_resolution * oversampling
    field_scratch = similar(optical_rate, T, focal_size, focal_size)
    focal = similar(optical_rate, Complex{T}, focal_size, focal_size)
    mode = similar(focal)
    conjugate_field = similar(focal)
    response, response_scratch, sampled, mapped =
        _allocate_lift_mapping_buffers(opd_map(tel), focal_resolution, mapping)
    output_dimensions = _lift_output_dimensions(focal_resolution, mapping)
    output_work = similar(opd_map(tel), T, output_dimensions...)
    convolution = kernel === nothing ? nothing : similar(optical_rate)
    convolution_scratch = kernel === nothing ? nothing : similar(optical_rate)
    opd = similar(pupil_amplitude)
    opd_work = similar(pupil_amplitude)

    workspace = LiFTForwardWorkspace(propagation, optical_rate,
        pupil_amplitude, field_scratch, focal, mode, conjugate_field,
        response, response_scratch, sampled, mapped, output_work, convolution,
        convolution_scratch, opd, opd_work)
    output_values = _lift_output_values(mapping, workspace)
    contract = LiFTObservationContract(
        OpticalPlaneMetadata(FocalPlane(), output_values;
            coordinate_domain=AngularCoordinates(),
            sampling=(T(wavelength(src) /
                (tel.params.diameter * zero_padding) *
                prod(_lift_mapping_factors(mapping))),
                T(wavelength(src) /
                (tel.params.diameter * zero_padding) *
                prod(_lift_mapping_factors(mapping)))),
            spectral=MonochromaticChannel(T(wavelength(src))),
            normalization=PhotonRateNormalization(),
            spatial_measure=CellIntegratedMeasure(),
            coherence=IncoherentIntensityAddition()),
        _lift_mapping_signature(mapping))
    model = LiFTForwardModel(pupil, amplitude, owned_basis, diversity,
        T(wavelength(src)), T(photon_irradiance(src)),
        T((tel.params.diameter / resolution)^2), focal_resolution,
        zero_padding, kernel, mapping, contract)
    output = IntensityMap(contract.rate_metadata, output_values)
    return PreparedLiFTForwardModel(model, workspace, output)
end

"""Return the immutable observation-compatibility contract for `forward`."""
@inline lift_observation_contract(forward::PreparedLiFTForwardModel) =
    forward.model.observation_contract
"""Return the caller-visible photon-arrival-rate output owned by `forward`."""
@inline lift_forward_output(forward::PreparedLiFTForwardModel) = forward.output

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
    device = plane_device(values)
    device == rate_metadata.device || throw(InvalidConfiguration(
        "LiFT observation and forward model must occupy the same physical device"))
    sigma = rate_metadata.numeric_type(readout_noise_std)
    isfinite(sigma) && sigma >= zero(sigma) || throw(InvalidConfiguration(
        "LiFT observation readout noise must be finite and nonnegative"))
    validate_values && _require_finite_nonnegative_intensity(values)
    metadata = LiFTObservationMetadata(contract, domain, sigma, E,
        backend(values), device)
    return LiFTObservation(metadata, values)
end

LiFTObservation(forward::PreparedLiFTForwardModel,
    values::AbstractMatrix; kwargs...) =
    LiFTObservation(lift_observation_contract(forward), values; kwargs...)

"""Construct an iterative LiFT estimator over a prepared forward model."""
function LiFT(forward::PreparedLiFTForwardModel; iterations::Int=5,
    mode_ids=axes(forward.model.basis, 3),
    numerical::Bool=false, solve_mode::LiFTSolveMode=LiFTSolveAuto(),
    damping::LiFTDampingMode=LiFTDampingNone())
    iterations >= 1 || throw(InvalidConfiguration(
        "LiFT iterations must be >= 1"))
    output = intensity_values(forward.output)
    T = eltype(output)
    prepared_mode_ids = Tuple(Int(mode_id) for mode_id in mode_ids)
    isempty(prepared_mode_ids) && throw(InvalidConfiguration(
        "LiFT mode_ids must not be empty"))
    all(mode_id -> 1 <= mode_id <= size(forward.model.basis, 3),
        prepared_mode_ids) || throw(DimensionMismatchError(
            "LiFT mode ids must index the prepared modal basis"))
    allunique(prepared_mode_ids) || throw(InvalidConfiguration(
        "LiFT mode_ids must be unique"))
    mode_count = length(prepared_mode_ids)
    observation_rate = similar(output)
    residual = similar(output, T, length(output))
    weights = similar(residual)
    H = similar(output, T, length(output), mode_count)
    normal = similar(output, T, mode_count, mode_count)
    factor = similar(normal)
    rhs = similar(residual, mode_count)
    mode_ids = similar(rhs, Int, mode_count)
    copyto!(mode_ids, collect(prepared_mode_ids))
    diag = LiFTDiagnostics(T(NaN), T(NaN), T(NaN), T(NaN), zero(T),
        false, false)
    state = LiFTState(observation_rate, residual, weights, H, normal,
        factor, rhs, mode_ids, diag)
    params = LiFTParams(iterations, solve_mode, damping, prepared_mode_ids)
    mode = numerical ? LiFTNumerical() : LiFTAnalytic()
    return LiFT{typeof(mode),typeof(forward),typeof(params),typeof(state)}(
        forward, params, state)
end

"""Return the diagnostics from the most recent LiFT reconstruction."""
diagnostics(lift::LiFT) = lift.state.diagnostics
@inline _lift_model(lift::LiFT) = lift.forward.model
@inline _lift_workspace(lift::LiFT) = lift.forward.workspace

"""
    prepare_opd!(lift, coeffs)

Assemble the current model OPD from the modal basis coefficients plus the fixed
diversity OPD.
"""
@inline function prepare_opd!(lift::LiFT, coeffs::AbstractVector)
    model = _lift_model(lift)
    workspace = _lift_workspace(lift)
    combine_basis!(workspace.opd_buffer, model.basis, coeffs,
        model.pupil_mask)
    @. workspace.opd_buffer += model.diversity_opd
    return workspace.opd_buffer
end

@inline function lift_affine_basis_mode!(dest::AbstractMatrix{T},
    base::AbstractMatrix{T}, basis::AbstractArray{T,3}, mode_id::Int,
    scale::T) where {T<:AbstractFloat}
    return lift_affine_basis_mode!(execution_style(dest), dest, base,
        basis, mode_id, scale)
end

function lift_affine_basis_mode!(::ScalarCPUStyle, dest::AbstractMatrix{T},
    base::AbstractMatrix{T}, basis::AbstractArray{T,3}, mode_id::Int,
    scale::T) where {T<:AbstractFloat}
    mode_offset = (mode_id - 1) * length(dest)
    @inbounds @simd for i in eachindex(dest, base)
        dest[i] = base[i] + scale * basis[i + mode_offset]
    end
    return dest
end

function lift_affine_basis_mode!(style::AcceleratorStyle,
    dest::AbstractMatrix{T}, base::AbstractMatrix{T},
    basis::AbstractArray{T,3}, mode_id::Int,
    scale::T) where {T<:AbstractFloat}
    n = length(dest)
    mode_offset = (mode_id - 1) * n
    launch_kernel!(style, lift_affine_basis_mode_kernel!, dest, base,
        basis, scale, mode_offset, n; ndrange=n)
    return dest
end

@inline function lift_scaled_basis_mode!(dest::AbstractMatrix{T},
    amplitude::AbstractMatrix{T}, basis::AbstractArray{T,3},
    mode_id::Int, scale::T) where {T<:AbstractFloat}
    return lift_scaled_basis_mode!(execution_style(dest), dest, amplitude,
        basis, mode_id, scale)
end

function lift_scaled_basis_mode!(::ScalarCPUStyle,
    dest::AbstractMatrix{T}, amplitude::AbstractMatrix{T},
    basis::AbstractArray{T,3}, mode_id::Int,
    scale::T) where {T<:AbstractFloat}
    mode_offset = (mode_id - 1) * length(dest)
    @inbounds @simd for i in eachindex(dest, amplitude)
        dest[i] = amplitude[i] * scale * basis[i + mode_offset]
    end
    return dest
end

function lift_scaled_basis_mode!(style::AcceleratorStyle,
    dest::AbstractMatrix{T}, amplitude::AbstractMatrix{T},
    basis::AbstractArray{T,3}, mode_id::Int,
    scale::T) where {T<:AbstractFloat}
    n = length(dest)
    mode_offset = (mode_id - 1) * n
    launch_kernel!(style, lift_scaled_basis_mode_kernel!, dest, amplitude,
        basis, scale, mode_offset, n; ndrange=n)
    return dest
end

@inline function lift_copy_column!(dest::AbstractMatrix,
    column::Int, src::AbstractMatrix)
    return lift_copy_column!(execution_style(dest), dest, column, src)
end

function lift_copy_column!(::ScalarCPUStyle, dest::AbstractMatrix{T},
    column::Int, src::AbstractMatrix{T}) where {T}
    @inbounds @simd for i in eachindex(src)
        dest[i, column] = src[i]
    end
    return dest
end

function lift_copy_column!(style::AcceleratorStyle,
    dest::AbstractMatrix, column::Int, src::AbstractMatrix)
    launch_kernel!(style, lift_copy_column_kernel!, dest, column, src,
        length(src); ndrange=length(src))
    return dest
end

@inline function lift_residual!(dest::AbstractVector,
    observation::AbstractMatrix, model::AbstractMatrix)
    return lift_residual!(execution_style(dest), dest, observation, model)
end

function lift_residual!(::ScalarCPUStyle, dest::AbstractVector{T},
    observation::AbstractMatrix{T}, model::AbstractMatrix{T}) where {T}
    @inbounds @simd for i in eachindex(dest, observation, model)
        dest[i] = observation[i] - model[i]
    end
    return dest
end

function lift_residual!(style::AcceleratorStyle, dest::AbstractVector,
    observation::AbstractMatrix, model::AbstractMatrix)
    launch_kernel!(style, lift_residual_kernel!, dest, observation, model,
        length(dest); ndrange=length(dest))
    return dest
end

@inline _apply_lift_mapping!(::LiFTIdentityMapping,
    ::LiFTForwardWorkspace, optical_rate::AbstractMatrix) = optical_rate

@inline function _apply_lift_response!(::NullFrameResponse,
    workspace::LiFTForwardWorkspace, optical_rate::AbstractMatrix)
    copyto!(workspace.response_buffer, optical_rate)
    return workspace.response_buffer
end

@inline function _apply_lift_response!(response::AbstractFrameResponse,
    workspace::LiFTForwardWorkspace, optical_rate::AbstractMatrix)
    copyto!(workspace.response_buffer, optical_rate)
    apply_response!(execution_style(workspace.response_buffer), response,
        workspace.response_buffer, workspace.response_scratch)
    return workspace.response_buffer
end

function _apply_lift_mapping!(mapping::LiFTFrameMapping,
    workspace::LiFTForwardWorkspace, optical_rate::AbstractMatrix)
    response_rate = _apply_lift_response!(mapping.response, workspace,
        optical_rate)
    if mapping.sampling > 1
        bin2d!(workspace.sampling_buffer, response_rate,
            mapping.sampling)
    else
        copyto!(workspace.sampling_buffer, response_rate)
    end
    if mapping.binning > 1
        bin2d!(workspace.mapped_rate_buffer, workspace.sampling_buffer,
            mapping.binning)
    else
        copyto!(workspace.mapped_rate_buffer, workspace.sampling_buffer)
    end
    return workspace.mapped_rate_buffer
end

"""
    lift_interaction_matrix!(H, lift, coefficients; rate_scale=1)

Fill the LiFT Jacobian with derivatives of the prepared observation-domain
photon-arrival rate with respect to the requested modal coefficients.

- `LiFTNumerical` uses centered finite differences
- `LiFTAnalytic` uses the field-derivative formulation
"""
function lift_interaction_matrix!(H::AbstractMatrix,
    lift::LiFT{LiFTNumerical}, coefficients::AbstractVector;
    rate_scale::Real=1.0)
    model = _lift_model(lift)
    workspace = _lift_workspace(lift)
    T = eltype(workspace.optical_rate_buffer)
    mode_ids = lift.params.mode_ids
    n_modes = length(mode_ids)
    output_size = length(lift.forward.output.values)
    if size(H, 1) < output_size || size(H, 2) < n_modes
        throw(InvalidConfiguration("H buffer size does not match LiFT dimensions"))
    end
    delta = T(1e-9)

    initial_opd = prepare_opd!(lift, coefficients)
    opd_work = workspace.opd_work_buffer
    @inbounds for (idx, mode_id) in enumerate(mode_ids)
        lift_affine_basis_mode!(opd_work, initial_opd, model.basis,
            mode_id, delta)
        rate_plus = _lift_rate_values_from_opd!(lift.forward, opd_work;
            rate_scale=rate_scale)
        copyto!(workspace.output_work_buffer, rate_plus)
        lift_affine_basis_mode!(opd_work, initial_opd, model.basis,
            mode_id, -delta)
        rate_minus = _lift_rate_values_from_opd!(lift.forward, opd_work;
            rate_scale=rate_scale)
        @. rate_minus = (workspace.output_work_buffer - rate_minus) /
            (2 * delta)
        lift_copy_column!(H, idx, rate_minus)
    end
    return H
end

function lift_interaction_matrix!(H::AbstractMatrix,
    lift::LiFT{LiFTAnalytic}, coefficients::AbstractVector;
    rate_scale::Real=1.0)
    model = _lift_model(lift)
    workspace = _lift_workspace(lift)
    T = eltype(workspace.optical_rate_buffer)
    mode_ids = lift.params.mode_ids
    n_modes = length(mode_ids)
    output_size = length(lift.forward.output.values)
    if size(H, 1) < output_size || size(H, 2) < n_modes
        throw(InvalidConfiguration("H buffer size does not match LiFT dimensions"))
    end

    initial_opd = prepare_opd!(lift, coefficients)

    amplitude_scale = sqrt(T(model.photon_irradiance * rate_scale *
        model.pupil_cell_area_m2))
    @. workspace.amplitude_buffer = model.pupil_amplitude * amplitude_scale

    oversampling = focal_field_from_opd!(workspace.focal_buffer,
        lift.forward, workspace.amplitude_buffer, initial_opd)
    conjugate_field!(workspace.conjugate_field_buffer,
        workspace.focal_buffer)
    conjugated_field = workspace.conjugate_field_buffer

    wavenumber = T(2 * pi) / model.wavelength_m
    @inbounds for (idx, mode_id) in enumerate(mode_ids)
        lift_scaled_basis_mode!(workspace.amplitude_buffer,
            model.pupil_amplitude, model.basis, mode_id, amplitude_scale)
        focal_field_from_opd!(workspace.mode_buffer, lift.forward,
            workspace.amplitude_buffer, initial_opd)
        field_derivative!(workspace.optical_rate_buffer,
            workspace.mode_buffer, conjugated_field, oversampling,
            2 * wavenumber, workspace.field_scratch)
        maybe_object_convolve!(lift.forward,
            workspace.optical_rate_buffer)
        mapped_derivative = _apply_lift_mapping!(model.mapping, workspace,
            workspace.optical_rate_buffer)
        lift_copy_column!(H, idx, mapped_derivative)
    end
    return H
end

function lift_interaction_matrix(lift::LiFT,
    coefficients::AbstractVector; rate_scale::Real=1.0)
    output = lift.forward.output.values
    H = similar(output, eltype(output), length(output),
        length(lift.params.mode_ids))
    return lift_interaction_matrix!(H, lift, coefficients;
        rate_scale=rate_scale)
end

function _require_lift_observation(lift::LiFT,
    observation::LiFTObservation)
    observation.metadata.contract == lift_observation_contract(lift.forward) ||
        throw(InvalidConfiguration(
            "LiFT observation geometry, wavelength, or preprocessing does not match the prepared forward model"))
    return observation
end

@inline _require_lift_weighting(::LiFTWeightingMode,
    ::LiFTObservationMetadata, ::Type{<:AbstractFloat}) = nothing

function _require_lift_weighting(mode::LiFTWeightMatrix,
    metadata::LiFTObservationMetadata,
    ::Type{T}) where {T<:AbstractFloat}
    size(mode.R_n) == metadata.contract.rate_metadata.dimensions || throw(
        DimensionMismatchError(
            "LiFT observation covariance must match observation dimensions"))
    eltype(mode.R_n) === T || throw(InvalidConfiguration(
        "LiFT observation covariance must use the estimator numeric type"))
    typeof(backend(mode.R_n)) === typeof(metadata.backend) || throw(
        InvalidConfiguration(
            "LiFT observation covariance must use the observation array backend"))
    plane_device(mode.R_n) == metadata.device || throw(InvalidConfiguration(
        "LiFT observation covariance must occupy the observation physical device"))
    return nothing
end

function _require_lift_reconstruction(lift::LiFT,
    observation::LiFTObservation, coeffs::AbstractVector{T},
    mode::LiFTWeightingMode, optimize_norm::Symbol) where {T<:AbstractFloat}
    _require_lift_observation(lift, observation)
    model_values = intensity_values(lift.forward.output)
    eltype(model_values) === T || throw(InvalidConfiguration(
        "LiFT coefficient buffer must use the prepared numeric type"))
    typeof(backend(coeffs)) === typeof(backend(model_values)) || throw(
        InvalidConfiguration(
            "LiFT coefficient buffer must use the prepared array backend"))
    plane_device(coeffs) == plane_device(model_values) || throw(
        InvalidConfiguration(
            "LiFT coefficient buffer must occupy the prepared physical device"))
    optimize_norm in (:sum, :max, :none) || throw(InvalidConfiguration(
        "LiFT optimize_norm must be :sum, :max, or :none"))
    length(coeffs) >= size(_lift_model(lift).basis, 3) || throw(
        DimensionMismatchError(
            "LiFT coefficient buffer must cover the prepared modal basis"))
    _require_lift_weighting(mode, observation.metadata, T)
    return nothing
end

@inline _require_lift_initial_coefficients(::AbstractVector,
    ::Nothing) = nothing

function _require_lift_initial_coefficients(dest::AbstractVector{T},
    initial::AbstractVector) where {T<:AbstractFloat}
    length(initial) == length(dest) || throw(DimensionMismatchError(
        "LiFT initial coefficients must match the full coefficient buffer"))
    eltype(initial) === T || throw(InvalidConfiguration(
        "LiFT initial coefficients must use the prepared numeric type"))
    typeof(backend(initial)) === typeof(backend(dest)) || throw(
        InvalidConfiguration(
            "LiFT initial coefficients must use the coefficient array backend"))
    plane_device(initial) == plane_device(dest) || throw(InvalidConfiguration(
        "LiFT initial coefficients must occupy the coefficient physical device"))
    return nothing
end

_require_lift_initial_coefficients(::AbstractVector, initial) = throw(
    InvalidConfiguration(
        "LiFT initial coefficients must be an abstract vector or nothing"))

function _copy_lift_observation_rate!(dest::AbstractMatrix{T},
    observation::LiFTObservation) where {T<:AbstractFloat}
    scale = lift_observation_to_rate_scale(observation.metadata.domain, T)
    @. dest = observation.values * scale
    return dest
end

function reconstruct(lift::LiFT, observation::LiFTObservation;
    coeffs0=nothing, R_n=nothing, optimize_norm::Symbol=:sum, check_convergence::Bool=true)
    T = eltype(lift.state.observation_rate_buffer)
    coeffs = similar(lift.state.rhs_buffer, T,
        size(_lift_model(lift).basis, 3))
    reconstruct!(coeffs, lift, observation; coeffs0=coeffs0, R_n=R_n,
        optimize_norm=optimize_norm, check_convergence=check_convergence)
    mode_ids_buf = lift.state.mode_id_buffer
    out = similar(lift.state.rhs_buffer, T, length(mode_ids_buf))
    gather_selected_coefficients!(execution_style(coeffs), out, coeffs,
        mode_ids_buf, length(mode_ids_buf))
    return out
end

"""
    reconstruct!(coeffs, lift, observation, weighting; ...)

Run the LiFT iterative reconstruction in-place.

Each iteration evaluates the current rate model, builds/weights the Jacobian,
solves for a modal update, and scatters that update back into the full
coefficient vector.
"""
function reconstruct!(coeffs::AbstractVector{T}, lift::LiFT,
    observation::LiFTObservation, mode::LiFTWeightingMode;
    optimize_norm::Symbol=:sum,
    check_convergence::Bool=true) where {T<:AbstractFloat}
    _require_lift_reconstruction(lift, observation, coeffs, mode,
        optimize_norm)
    n_modes = length(lift.params.mode_ids)
    observation_rate = _copy_lift_observation_rate!(
        lift.state.observation_rate_buffer, observation)

    residual = lift.state.residual_buffer
    sqrtw = lift.state.weight_buffer
    H = lift.state.H_buffer
    normal = lift.state.normal_buffer
    factor = lift.state.factor_buffer
    rhs = lift.state.rhs_buffer
    mode_ids_buf = lift.state.mode_id_buffer
    diag = lift.state.diagnostics
    style = execution_style(H)
    effective_mode = effective_solve_mode(style, lift.params.solve_mode)
    λ_state = initial_damping_state(lift.params.damping, T)
    prev_weighted_residual_norm = T(Inf)
    init_weights!(sqrtw, mode, observation_rate, observation.metadata)
    for iter in 1:lift.params.iterations
        current_opd = prepare_opd!(lift, coeffs)
        model_rate = _lift_rate_values_from_opd!(lift.forward, current_opd)
        scale = one(T)
        model_photon_rate = zero(T)
        if optimize_norm == :sum
            denom = backend_sum_value(model_rate)
            if denom > 0
                scale = backend_sum_value(observation_rate) / denom
            end
            model_photon_rate = abs(denom * scale)
        elseif optimize_norm == :max
            denom = backend_maximum_value(model_rate)
            if denom > 0
                scale = backend_maximum_value(observation_rate) / denom
            end
        end
        if scale != one(T)
            model_rate .*= scale
        end
        if optimize_norm != :sum
            model_photon_rate = abs(backend_sum_value(model_rate))
        end
        objective_scale = max(model_photon_rate, one(T))
        lift_residual!(residual, observation_rate, model_rate)
        diag.residual_norm = norm(residual)
        update_weights!(sqrtw, mode, model_rate, observation.metadata)
        # Generate the Jacobian in the normalized objective scale so a
        # physically large photon rate cannot overflow before H'H is formed.
        lift_interaction_matrix!(H, lift, coeffs;
            rate_scale=scale / objective_scale)

        apply_row_weights!(H, sqrtw, n_modes)
        apply_vec_weights!(residual, sqrtw)
        diag.weighted_residual_norm = norm(residual)
        # The same normalization on the residual leaves the unregularized
        # least-squares update unchanged. Damping is evaluated on this
        # normalized objective rather than on the raw photon-rate scale.
        residual .*= inv(objective_scale)
        mul!(normal, adjoint(H), H)
        mul!(rhs, adjoint(H), residual)
        λ_state = update_damping_state(lift.params.damping, λ_state, prev_weighted_residual_norm,
            diag.weighted_residual_norm, normal)
        damping = effective_damping(lift.params.damping, λ_state)
        delta = solve_lift_system!(diag, residual, rhs, H, normal, factor, effective_mode, damping)
        diag.update_norm = delta_norm(delta, n_modes, effective_mode)
        update_coefficients!(coeffs, delta, mode_ids_buf, effective_mode, style)
        prev_weighted_residual_norm = diag.weighted_residual_norm
        if check_convergence && diag.update_norm / max(norm(coeffs), eps(T)) < 1e-3
            break
        end
    end

    return coeffs
end

function reconstruct!(coeffs::AbstractVector, lift::LiFT,
    observation::LiFTObservation;
    coeffs0=nothing, R_n=nothing, optimize_norm::Symbol=:sum, check_convergence::Bool=true)
    mode = weight_mode(R_n)
    _require_lift_reconstruction(lift, observation, coeffs, mode,
        optimize_norm)
    _require_lift_initial_coefficients(coeffs, coeffs0)
    if coeffs0 === nothing
        fill!(coeffs, zero(eltype(coeffs)))
    else
        copyto!(coeffs, coeffs0)
    end
    return reconstruct!(coeffs, lift, observation, mode;
        optimize_norm=optimize_norm, check_convergence=check_convergence)
end

@inline function apply_vec_weights!(vec::AbstractVector{T}, weights::AbstractVector{T}) where {T<:AbstractFloat}
    @. vec *= weights
    return vec
end

@inline function apply_row_weights!(mat::AbstractMatrix{T},
    weights::AbstractVector{T}, n_cols::Int) where {T<:AbstractFloat}
    return apply_row_weights!(execution_style(mat), mat, weights, n_cols)
end

function apply_row_weights!(::ScalarCPUStyle, mat::AbstractMatrix{T},
    weights::AbstractVector{T}, n_cols::Int) where {T<:AbstractFloat}
    n_rows = size(mat, 1)
    @inbounds for j in 1:n_cols
        @simd for i in 1:n_rows
            mat[i, j] *= weights[i]
        end
    end
    return mat
end

function apply_row_weights!(style::AcceleratorStyle,
    mat::AbstractMatrix, weights::AbstractVector, n_cols::Int)
    n_rows = size(mat, 1)
    launch_kernel!(style, lift_row_weights_kernel!, mat, weights,
        n_rows, n_cols; ndrange=(n_rows, n_cols))
    return mat
end

@inline effective_solve_mode(::ScalarCPUStyle, ::LiFTSolveAuto) = LiFTSolveQR()
@inline effective_solve_mode(::AcceleratorStyle, ::LiFTSolveAuto) = LiFTSolveNormalEquations()
@inline effective_solve_mode(::ExecutionStyle, mode::LiFTSolveMode) = mode

function solve_lift_fallback!(diag::LiFTDiagnostics{T}, rhs::AbstractVector{T},
    H::AbstractMatrix{T}, residual::AbstractVector{T}, damping::LiFTDampingMode) where {T<:AbstractFloat}
    λ = fallback_damping_lambda(damping, T, H)
    F = svd(H; full=false)
    s = F.S
    work = similar(rhs)
    mul!(work, transpose(F.U), residual)
    @inbounds for i in eachindex(s)
        denom = s[i]^2 + λ
        work[i] = iszero(denom) ? zero(T) : (s[i] * work[i]) / denom
    end
    mul!(rhs, F.V, work)
    diag.regularization = λ
    diag.used_fallback = true
    return rhs
end

fallback_damping_lambda(::LiFTDampingMode, ::Type{T}, H::AbstractMatrix{T}) where {T<:AbstractFloat} = zero(T)
fallback_damping_lambda(damping::LiFTLevenbergMarquardt, ::Type{T}, H::AbstractMatrix{T}) where {T<:AbstractFloat} =
    damping_lambda(damping, transpose(H) * H)

@inline function add_lift_diagonal!(matrix::AbstractMatrix{T},
    value::T) where {T<:AbstractFloat}
    return add_lift_diagonal!(execution_style(matrix), matrix, value)
end

function add_lift_diagonal!(::ScalarCPUStyle,
    matrix::AbstractMatrix{T}, value::T) where {T<:AbstractFloat}
    @inbounds for i in 1:min(size(matrix, 1), size(matrix, 2))
        matrix[i, i] += value
    end
    return matrix
end

function add_lift_diagonal!(style::AcceleratorStyle,
    matrix::AbstractMatrix{T}, value::T) where {T<:AbstractFloat}
    n = min(size(matrix, 1), size(matrix, 2))
    launch_kernel!(style, lift_add_diagonal_kernel!, matrix, value, n;
        ndrange=n)
    return matrix
end

"""
    solve_normal_system!(diag, rhs, factor, normal, H, residual, damping)

Solve the LiFT normal equations for the current modal update.

The preferred path is Cholesky on the normal matrix. If the factorization is
ill-conditioned, the implementation adds diagonal loading and eventually falls
back to the SVD-based solve.
"""
function solve_normal_system!(diag::LiFTDiagnostics{T}, rhs::AbstractVector{T}, factor::AbstractMatrix{T},
    normal::AbstractMatrix{T}, H::AbstractMatrix{T}, residual::AbstractVector{T},
    ::LiFTDampingNone) where {T<:AbstractFloat}
    copyto!(factor, normal)
    chol = cholesky!(Hermitian(factor), check=false)
    λ = zero(T)
    if !issuccess(chol)
        λ = regularization_load(normal)
        add_lift_diagonal!(factor, λ)
        chol = cholesky!(Hermitian(factor), check=false)
        if !issuccess(chol)
            λ *= T(10)
            add_lift_diagonal!(factor, λ)
            chol = cholesky!(Hermitian(factor), check=false)
            if !issuccess(chol)
                return solve_lift_fallback!(diag, rhs, H, residual, LiFTLevenbergMarquardt(lambda0=λ))
            end
        end
    end
    ldiv!(chol, rhs)
    diag.regularization = λ
    return rhs
end

function solve_normal_system!(diag::LiFTDiagnostics{T}, rhs::AbstractVector{T}, factor::AbstractMatrix{T},
    normal::AbstractMatrix{T}, H::AbstractMatrix{T}, residual::AbstractVector{T},
    damping::LiFTLevenbergMarquardt) where {T<:AbstractFloat}
    copyto!(factor, normal)
    λ = damping_lambda(damping, normal)
    if λ > zero(T)
        add_lift_diagonal!(factor, λ)
    end
    chol = cholesky!(Hermitian(factor), check=false)
    while !issuccess(chol)
        λ = max(λ * T(damping.growth), regularization_load(normal))
        copyto!(factor, normal)
        add_lift_diagonal!(factor, λ)
        chol = cholesky!(Hermitian(factor), check=false)
        if λ > T(1e12)
            return solve_lift_fallback!(diag, rhs, H, residual, damping)
        end
    end
    ldiv!(chol, rhs)
    diag.regularization = λ
    return rhs
end

function solve_lift_system!(diag::LiFTDiagnostics{T}, residual::AbstractVector{T}, rhs::AbstractVector{T},
    H::AbstractMatrix{T}, normal::AbstractMatrix{T}, factor::AbstractMatrix{T},
    effective_mode::LiFTSolveQR, damping::LiFTDampingMode) where {T<:AbstractFloat}
    diag.used_qr = true
    diag.used_fallback = false
    diag.regularization = zero(T)
    qr_factor = qr!(H)
    cond_ratio = qr_condition_ratio(qr_factor, size(normal, 1))
    diag.condition_ratio = cond_ratio
    if cond_ratio > damping_condition_limit(T, damping)
        diag.used_qr = false
        diag.used_fallback = true
        solve_normal_system!(diag, rhs, factor, normal, H, residual, damping)
        return rhs
    end
    try
        ldiv!(qr_factor, residual)
        return residual
    catch err
        return handle_lift_qr_error!(err, diag, rhs, factor, normal, H, residual, damping)
    end
end

function solve_lift_system!(diag::LiFTDiagnostics{T}, residual::AbstractVector{T}, rhs::AbstractVector{T},
    H::AbstractMatrix{T}, normal::AbstractMatrix{T}, factor::AbstractMatrix{T},
    effective_mode::LiFTSolveMode, damping::LiFTDampingMode) where {T<:AbstractFloat}
    diag.used_qr = false
    diag.used_fallback = false
    diag.regularization = zero(T)
    diag.condition_ratio = normal_condition_ratio(normal)
    solve_normal_system!(diag, rhs, factor, normal, H, residual, damping)
    return rhs
end

function handle_lift_qr_error!(::SingularException, diag::LiFTDiagnostics, rhs::AbstractVector,
    factor::AbstractMatrix, normal::AbstractMatrix, H::AbstractMatrix, residual::AbstractVector,
    damping::LiFTDampingMode)
    diag.used_qr = false
    diag.used_fallback = true
    solve_normal_system!(diag, rhs, factor, normal, H, residual, damping)
    return rhs
end

function handle_lift_qr_error!(err, diag::LiFTDiagnostics, rhs::AbstractVector,
    factor::AbstractMatrix, normal::AbstractMatrix, H::AbstractMatrix, residual::AbstractVector,
    damping::LiFTDampingMode)
    throw(err)
end

@inline delta_component(delta::AbstractVector, j::Int, ::LiFTSolveMode) = delta[j]

@inline function update_coefficients!(coeffs::AbstractVector, delta::AbstractVector, mode_ids::AbstractVector,
    mode::LiFTSolveMode, ::ScalarCPUStyle)
    for (j, mode_id) in enumerate(mode_ids)
        coeffs[mode_id] += delta_component(delta, j, mode)
    end
    return coeffs
end

@inline function update_coefficients!(coeffs::AbstractVector, delta::AbstractVector, mode_ids::AbstractVector,
    ::LiFTSolveMode, style::AcceleratorStyle)
    launch_kernel!(style, lift_scatter_update_kernel!, coeffs, delta, mode_ids, length(mode_ids); ndrange=length(mode_ids))
    return coeffs
end

@inline function gather_selected_coefficients!(::ScalarCPUStyle,
    out::AbstractVector, coeffs::AbstractVector,
    mode_ids::AbstractVector, n_modes::Int)
    @inbounds @simd for i in 1:n_modes
        out[i] = coeffs[mode_ids[i]]
    end
    return out
end

@inline function gather_selected_coefficients!(style::AcceleratorStyle, out::AbstractVector, coeffs::AbstractVector,
    mode_ids::AbstractVector, n_modes::Int)
    launch_kernel!(style, lift_gather_kernel!, out, coeffs, mode_ids, n_modes; ndrange=n_modes)
    return out
end

@inline function delta_norm(delta::AbstractVector{T}, n_modes::Int, ::LiFTSolveQR) where {T<:AbstractFloat}
    return norm(@view delta[1:n_modes])
end

@inline function delta_norm(delta::AbstractVector{T}, n_modes::Int, ::LiFTSolveNormalEquations) where {T<:AbstractFloat}
    return norm(delta)
end

function regularization_load(normal::AbstractMatrix{T}) where {T<:AbstractFloat}
    n = min(size(normal, 1), size(normal, 2))
    diagvals = @view normal[diagind(normal)]
    diag_sum = sum(abs, diagvals)
    mean_diag = n == 0 ? zero(T) : diag_sum / T(n)
    return max(sqrt(eps(T)) * max(mean_diag, one(T)), eps(T))
end

function damping_lambda(damping::LiFTLevenbergMarquardt, normal::AbstractMatrix{T}) where {T<:AbstractFloat}
    base = regularization_load(normal)
    return max(T(damping.lambda0) * max(base, one(T)), base)
end

function damping_lambda(damping::LiFTAdaptiveLevenbergMarquardt, normal::AbstractMatrix{T}) where {T<:AbstractFloat}
    base = regularization_load(normal)
    return max(T(damping.lambda0) * max(base, one(T)), T(damping.min_lambda), base)
end

@inline damping_condition_limit(::Type{T}, ::LiFTDampingNone) where {T<:AbstractFloat} = inv(sqrt(eps(T)))
@inline damping_condition_limit(::Type{T}, damping::LiFTLevenbergMarquardt) where {T<:AbstractFloat} =
    inv(max(T(damping.condition_rtol), eps(T)))
@inline damping_condition_limit(::Type{T}, damping::LiFTAdaptiveLevenbergMarquardt) where {T<:AbstractFloat} =
    inv(max(T(damping.condition_rtol), eps(T)))

@inline initial_damping_state(::LiFTDampingNone, ::Type{T}) where {T<:AbstractFloat} = zero(T)
@inline initial_damping_state(damping::LiFTLevenbergMarquardt, ::Type{T}) where {T<:AbstractFloat} = T(damping.lambda0)
@inline initial_damping_state(damping::LiFTAdaptiveLevenbergMarquardt, ::Type{T}) where {T<:AbstractFloat} = T(damping.lambda0)

@inline update_damping_state(::LiFTDampingNone, λ::T, ::T, ::T, ::AbstractMatrix{T}) where {T<:AbstractFloat} = λ
@inline update_damping_state(::LiFTLevenbergMarquardt, λ::T, ::T, ::T, ::AbstractMatrix{T}) where {T<:AbstractFloat} = λ
function update_damping_state(damping::LiFTAdaptiveLevenbergMarquardt, λ::T, prev_residual::T,
    current_residual::T, ::AbstractMatrix{T}) where {T<:AbstractFloat}
    min_lambda = T(damping.min_lambda)
    if !isfinite(prev_residual)
        return max(T(damping.lambda0), min_lambda)
    elseif current_residual > prev_residual * (one(T) + sqrt(eps(T)))
        return max(λ * T(damping.growth), min_lambda)
    end
    return max(λ / T(damping.shrink), min_lambda)
end

@inline effective_damping(::LiFTDampingNone, ::T) where {T<:AbstractFloat} = LiFTDampingNone()
@inline effective_damping(damping::LiFTLevenbergMarquardt, ::T) where {T<:AbstractFloat} = damping
function effective_damping(damping::LiFTAdaptiveLevenbergMarquardt, λ::T) where {T<:AbstractFloat}
    return LiFTLevenbergMarquardt(lambda0=max(λ, T(damping.min_lambda)),
        growth=damping.growth, condition_rtol=damping.condition_rtol)
end

function qr_condition_ratio(qr_factor, n_modes::Int)
    T = real(eltype(qr_factor.factors))
    maxabs = zero(T)
    minabs = typemax(T)
    @inbounds for i in 1:n_modes
        d = abs(qr_factor.factors[i, i])
        maxabs = max(maxabs, d)
        minabs = min(minabs, d)
    end
    maxabs == zero(T) && return T(Inf)
    return maxabs / max(minabs, eps(T))
end

@inline normal_condition_ratio(normal::AbstractMatrix{T}) where {T<:AbstractFloat} =
    normal_condition_ratio(reduction_execution_plan(normal), normal)

function normal_condition_ratio(::DirectReductionPlan,
    normal::Array{T,2}) where {T<:AbstractFloat}
    maxabs = zero(T)
    minabs = typemax(T)
    @inbounds for i in axes(normal, 1)
        value = abs(normal[i, i])
        maxabs = max(maxabs, value)
        minabs = min(minabs, value)
    end
    maxabs == zero(T) && return T(Inf)
    return maxabs / max(minabs, eps(T))
end


function normal_condition_ratio(::DirectReductionPlan,
    normal::AbstractMatrix{T}) where {T<:AbstractFloat}
    diagonal = @view normal[diagind(normal)]
    maxabs = maximum(abs, diagonal)
    minabs = minimum(abs, diagonal)
    maxabs == zero(T) && return T(Inf)
    return maxabs / max(minabs, eps(T))
end

function normal_condition_ratio(::HostMirrorReductionPlan, normal::AbstractMatrix{T}) where {T<:AbstractFloat}
    host_parent = Array(reduction_parent_source(normal))
    host_normal = reduction_host_view(host_parent, normal)
    return normal_condition_ratio(DirectReductionPlan(), host_normal)
end

@inline function init_weights!(sqrtw::AbstractVector{T}, ::LiFTWeightingDynamic,
    ::AbstractMatrix{T}, ::LiFTObservationMetadata) where {T<:AbstractFloat}
    return sqrtw
end

@inline function init_weights!(sqrtw::AbstractVector{T}, mode::LiFTWeightingStatic,
    observation_rate::AbstractMatrix{T},
    metadata::LiFTObservationMetadata) where {T<:AbstractFloat}
    weight_vector!(sqrtw, observation_rate, mode, metadata)
    sqrt_weights!(execution_style(sqrtw), sqrtw)
    return sqrtw
end

@inline function update_weights!(sqrtw::AbstractVector{T}, ::LiFTWeightingStatic,
    ::AbstractMatrix{T}, ::LiFTObservationMetadata) where {T<:AbstractFloat}
    return sqrtw
end

@inline function update_weights!(sqrtw::AbstractVector{T}, mode::LiFTWeightingDynamic,
    model_rate::AbstractMatrix{T},
    metadata::LiFTObservationMetadata) where {T<:AbstractFloat}
    weight_vector!(sqrtw, model_rate, mode, metadata)
    sqrt_weights!(execution_style(sqrtw), sqrtw)
    return sqrtw
end


@inline function sqrt_weights!(::ScalarCPUStyle, weights::AbstractVector)
    @inbounds @simd for i in eachindex(weights)
        weights[i] = sqrt(weights[i])
    end
    return weights
end

@inline function sqrt_weights!(style::AcceleratorStyle, weights::AbstractVector)
    launch_kernel!(style, lift_sqrt_weights_kernel!, weights, length(weights);
        ndrange=length(weights))
    return weights
end

@inline function _lift_readout_variance_rate(
    metadata::LiFTObservationMetadata, ::Type{T}) where {T<:AbstractFloat}
    scale = lift_observation_to_rate_scale(metadata.domain, T)
    sigma_rate = T(metadata.readout_noise_std) * scale
    return sigma_rate * sigma_rate
end

@inline function lift_inverse_variance!(dest::AbstractVector{T},
    values::AbstractArray{T}, scale::T, offset::T) where {T<:AbstractFloat}
    return lift_inverse_variance!(execution_style(dest), dest, values,
        scale, offset)
end

function lift_inverse_variance!(::ScalarCPUStyle,
    dest::AbstractVector{T}, values::AbstractArray{T}, scale::T,
    offset::T) where {T<:AbstractFloat}
    floor_value = eps(T)
    @inbounds @simd for i in eachindex(dest, values)
        dest[i] = inv(max(values[i] * scale + offset, floor_value))
    end
    return dest
end

function lift_inverse_variance!(style::AcceleratorStyle,
    dest::AbstractVector{T}, values::AbstractArray{T}, scale::T,
    offset::T) where {T<:AbstractFloat}
    launch_kernel!(style, lift_inverse_variance_kernel!, dest, values,
        scale, offset, eps(T), length(dest); ndrange=length(dest))
    return dest
end

function weight_vector!(out::AbstractVector{T},
    model_rate::AbstractMatrix{T}, ::LiFTWeightModel,
    metadata::LiFTObservationMetadata) where {T<:AbstractFloat}
    readout_variance = _lift_readout_variance_rate(metadata, T)
    shot_scale = lift_shot_variance_rate_scale(metadata.domain, T)
    return lift_inverse_variance!(out, model_rate, shot_scale,
        readout_variance)
end

function weight_vector!(out::AbstractVector{T},
    model_rate::AbstractMatrix{T}, ::LiFTWeightIterative,
    metadata::LiFTObservationMetadata) where {T<:AbstractFloat}
    return weight_vector!(out, model_rate, LiFTWeightModel(), metadata)
end

function weight_vector!(out::AbstractVector{T}, ::AbstractMatrix{T}, ::LiFTWeightNone,
    metadata::LiFTObservationMetadata) where {T<:AbstractFloat}
    readout_variance = _lift_readout_variance_rate(metadata, T)
    fill!(out, one(T) / max(readout_variance, eps(T)))
    return out
end

function weight_vector!(out::AbstractVector{T}, ::AbstractMatrix{T}, mode::LiFTWeightMatrix,
    metadata::LiFTObservationMetadata) where {T<:AbstractFloat}
    variance_scale = lift_observation_to_rate_scale(metadata.domain, T)^2
    return lift_inverse_variance!(out, mode.R_n, variance_scale, zero(T))
end

function weight_vector!(out::AbstractVector{T}, model_rate::AbstractMatrix{T},
    R_n, metadata::LiFTObservationMetadata) where {T<:AbstractFloat}
    return weight_vector!(out, model_rate, weight_mode(R_n), metadata)
end

function weight_vector(model_rate::AbstractMatrix{T}, R_n,
    metadata::LiFTObservationMetadata) where {T<:AbstractFloat}
    out = similar(model_rate, T, length(model_rate))
    return weight_vector!(out, model_rate, weight_mode(R_n), metadata)
end

"""
    evaluate_lift_forward!(forward, opd; rate_scale=1)

Evaluate the prepared LiFT forward model from an OPD map in metres.

The returned `IntensityMap` contains cell-integrated photon-arrival rate on the
prepared observation grid after object convolution and deterministic spatial
preprocessing. No exposure, QE, noise, or readout operation is applied.
"""
function _lift_rate_values_from_opd!(forward::PreparedLiFTForwardModel,
    opd::AbstractMatrix; rate_scale::Real=1.0)
    model = forward.model
    workspace = forward.workspace
    size(opd) == size(model.diversity_opd) || throw(DimensionMismatchError(
        "LiFT forward OPD must match the prepared pupil dimensions"))
    T = eltype(workspace.optical_rate_buffer)
    scale = T(rate_scale)
    isfinite(scale) && scale >= zero(T) || throw(InvalidConfiguration(
        "LiFT forward rate_scale must be finite and nonnegative"))
    amplitude_scale = sqrt(model.photon_irradiance * scale *
        model.pupil_cell_area_m2)
    @. workspace.amplitude_buffer = model.pupil_amplitude * amplitude_scale
    oversampling = focal_field_from_opd!(workspace.focal_buffer, forward,
        workspace.amplitude_buffer, opd)
    field_intensity!(workspace.optical_rate_buffer, workspace.focal_buffer,
        oversampling, workspace.field_scratch)
    maybe_object_convolve!(forward, workspace.optical_rate_buffer)
    return _apply_lift_mapping!(model.mapping, workspace,
        workspace.optical_rate_buffer)
end

function evaluate_lift_forward!(forward::PreparedLiFTForwardModel,
    opd::AbstractMatrix; rate_scale::Real=1.0)
    values = _lift_rate_values_from_opd!(forward, opd;
        rate_scale=rate_scale)
    values === forward.output.values || throw(InvalidConfiguration(
        "LiFT forward workspace output does not match its prepared contract"))
    return forward.output
end

"""
    predict_lift_observation!(dest, forward, opd, domain)

Evaluate `forward` at `opd` and write values in the explicitly selected native
observation `domain`. The caller owns `dest`; no detector is invoked.
"""
function predict_lift_observation!(dest::AbstractMatrix,
    forward::PreparedLiFTForwardModel, opd::AbstractMatrix,
    domain::AbstractLiFTObservationDomain)
    size(dest) == size(forward.output.values) || throw(DimensionMismatchError(
        "LiFT prediction destination must match the observation dimensions"))
    typeof(backend(dest)) === typeof(backend(forward.output.values)) || throw(
        InvalidConfiguration(
            "LiFT prediction destination must use the prepared array backend"))
    plane_device(dest) == plane_device(forward.output.values) || throw(
        InvalidConfiguration(
            "LiFT prediction destination must occupy the prepared physical device"))
    rate = _lift_rate_values_from_opd!(forward, opd)
    T = eltype(rate)
    native_scale = inv(lift_observation_to_rate_scale(domain, T))
    @. dest = rate * native_scale
    return dest
end

function center_crop!(dest::AbstractMatrix, src::AbstractMatrix)
    if size(dest) == size(src)
        copyto!(dest, src)
        return dest
    end
    n = size(dest, 1)
    cx = div(size(src, 1) - n, 2)
    cy = div(size(src, 2) - n, 2)
    @views copyto!(dest, src[cx+1:cx+n, cy+1:cy+n])
    return dest
end

@inline function maybe_object_convolve!(
    forward::PreparedLiFTForwardModel, matrix::AbstractMatrix)
    return _maybe_object_convolve!(forward.model.object_kernel,
        forward.workspace, matrix)
end

@inline _maybe_object_convolve!(::Nothing, ::LiFTForwardWorkspace,
    matrix::AbstractMatrix) = matrix

function _maybe_object_convolve!(kernel::LiFTDenseObjectKernel,
    workspace::LiFTForwardWorkspace, matrix::AbstractMatrix)
    conv2d_same!(workspace.convolution_buffer, matrix, kernel.kernel,
        kernel.inv_norm)
    copyto!(matrix, workspace.convolution_buffer)
    return matrix
end

function _maybe_object_convolve!(kernel::LiFTSeparableObjectKernel,
    workspace::LiFTForwardWorkspace, matrix::AbstractMatrix)
    conv2d_same_separable!(workspace.convolution_buffer,
        workspace.convolution_scratch, matrix, kernel.row, kernel.col,
        kernel.inv_norm)
    copyto!(matrix, workspace.convolution_buffer)
    return matrix
end

function _lift_object_kernel(kernel::AbstractMatrix{T}) where {T<:AbstractFloat}
    row, col = _separable_kernel_factors(kernel)
    if row === nothing
        norm = sum(Array(kernel))
        inv_norm = iszero(norm) ? one(T) : inv(T(norm))
        return LiFTDenseObjectKernel{T,typeof(kernel)}(kernel, inv_norm)
    end
    row_norm = sum(Array(row))
    col_norm = sum(Array(col))
    inv_norm = (iszero(row_norm) || iszero(col_norm)) ? one(T) :
        inv(T(row_norm * col_norm))
    return LiFTSeparableObjectKernel{T,typeof(row)}(row, col, inv_norm)
end

function _separable_kernel_factors(kernel::AbstractMatrix{T}) where {T<:AbstractFloat}
    F = svd(Matrix(kernel); full=false)
    isempty(F.S) && return nothing, nothing
    σ1 = F.S[1]
    σ2 = length(F.S) >= 2 ? F.S[2] : zero(T)
    tol = sqrt(eps(T)) * max(one(T), σ1)
    σ2 <= tol || return nothing, nothing
    scale = sqrt(σ1)
    host_row = T.(F.U[:, 1] .* scale)
    host_col = T.(F.V[:, 1] .* scale)
    row = similar(kernel, T, length(host_row))
    col = similar(kernel, T, length(host_col))
    copyto!(row, host_row)
    copyto!(col, host_col)
    return row, col
end

function conv2d_same!(dest::AbstractMatrix{T}, src::AbstractMatrix{T},
    kernel::AbstractMatrix) where {T<:AbstractFloat}
    norm = backend_sum_value(kernel)
    inv_norm = iszero(norm) ? one(T) : inv(T(norm))
    return conv2d_same!(dest, src, kernel, inv_norm)
end

@inline function conv2d_same!(dest::AbstractMatrix{T},
    src::AbstractMatrix{T}, kernel::AbstractMatrix,
    inv_norm::T) where {T<:AbstractFloat}
    return conv2d_same!(execution_style(dest), dest, src, kernel, inv_norm)
end

function conv2d_same!(::ScalarCPUStyle, dest::AbstractMatrix{T},
    src::AbstractMatrix{T}, kernel::AbstractMatrix,
    inv_norm::T) where {T<:AbstractFloat}
    n, m = size(src)
    kh, kw = size(kernel)
    cx = div(kh, 2)
    cy = div(kw, 2)
    fill!(dest, zero(T))
    @inbounds for ki in 1:kh, kj in 1:kw
        offset_i = ki - cx - 1
        offset_j = kj - cy - 1
        weight = kernel[ki, kj]
        for j in 1:m
            jj = symm_index(j + offset_j, m)
            @simd for i in 1:n
                ii = symm_index(i + offset_i, n)
                dest[i, j] += src[ii, jj] * weight
            end
        end
    end
    @inbounds @simd for index in eachindex(dest)
        dest[index] *= inv_norm
    end
    return dest
end

function conv2d_same!(style::AcceleratorStyle,
    dest::AbstractMatrix{T}, src::AbstractMatrix{T},
    kernel::AbstractMatrix, inv_norm::T) where {T<:AbstractFloat}
    n, m = size(src)
    kh, kw = size(kernel)
    launch_kernel!(style, lift_dense_convolution_kernel!, dest, src,
        kernel, inv_norm, n, m, kh, kw; ndrange=(n, m))
    return dest
end

function conv2d_same_separable!(dest::AbstractMatrix{T},
    tmp::AbstractMatrix{T}, src::AbstractMatrix{T},
    row_kernel::AbstractVector{T},
    col_kernel::AbstractVector{T}) where {T<:AbstractFloat}
    row_norm = backend_sum_value(row_kernel)
    col_norm = backend_sum_value(col_kernel)
    inv_norm = (iszero(row_norm) || iszero(col_norm)) ? one(T) :
        inv(T(row_norm * col_norm))
    return conv2d_same_separable!(dest, tmp, src, row_kernel,
        col_kernel, inv_norm)
end

@inline function conv2d_same_separable!(dest::AbstractMatrix{T},
    tmp::AbstractMatrix{T}, src::AbstractMatrix{T},
    row_kernel::AbstractVector{T}, col_kernel::AbstractVector{T},
    inv_norm::T) where {T<:AbstractFloat}
    return conv2d_same_separable!(execution_style(dest), dest, tmp, src,
        row_kernel, col_kernel, inv_norm)
end

function conv2d_same_separable!(::ScalarCPUStyle,
    dest::AbstractMatrix{T}, tmp::AbstractMatrix{T},
    src::AbstractMatrix{T}, row_kernel::AbstractVector{T},
    col_kernel::AbstractVector{T}, inv_norm::T) where {T<:AbstractFloat}
    n, m = size(src)
    kr = length(row_kernel)
    kc = length(col_kernel)
    cx = div(kr, 2)
    cy = div(kc, 2)
    @inbounds for j in 1:m
        for i in 1:n
            acc = zero(T)
            for ki in 1:kr
                ii = symm_index(i + ki - cx - 1, n)
                acc += src[ii, j] * row_kernel[ki]
            end
            tmp[i, j] = acc
        end
    end
    fill!(dest, zero(T))
    @inbounds for kj in 1:kc
        offset_j = kj - cy - 1
        weight = col_kernel[kj]
        for j in 1:m
            jj = symm_index(j + offset_j, m)
            @simd for i in 1:n
                dest[i, j] += tmp[i, jj] * weight
            end
        end
    end
    @inbounds @simd for index in eachindex(dest)
        dest[index] *= inv_norm
    end
    return dest
end

function conv2d_same_separable!(style::AcceleratorStyle,
    dest::AbstractMatrix{T}, tmp::AbstractMatrix{T},
    src::AbstractMatrix{T}, row_kernel::AbstractVector{T},
    col_kernel::AbstractVector{T}, inv_norm::T) where {T<:AbstractFloat}
    n, m = size(src)
    phase = begin_kernel_phase(style)
    queue_kernel!(phase, lift_row_convolution_kernel!, tmp, src,
        row_kernel, n, m, length(row_kernel); ndrange=(n, m))
    queue_kernel!(phase, lift_column_convolution_kernel!, dest, tmp,
        col_kernel, inv_norm, n, m, length(col_kernel); ndrange=(n, m))
    finish_kernel_phase!(phase)
    return dest
end

@inline function symm_index(i::Int, n::Int)
    n == 1 && return 1
    while i < 1 || i > n
        if i < 1
            i = 2 - i
        else
            i = 2 * n - i
        end
    end
    return i
end

@inline function lift_oversampling(zero_padding::Int)
    zero_padding < 1 && throw(InvalidConfiguration("LiFT zero_padding must be >= 1"))
    return zero_padding < 2 ? cld(2, zero_padding) : 1
end

@inline function lift_pad_size(resolution::Int, zero_padding::Int)
    oversampling = lift_oversampling(zero_padding)
    nominal = zero_padding * oversampling * resolution
    pad_width = cld(nominal - resolution, 2)
    return resolution + 2 * pad_width
end

function focal_field_from_opd!(dest::AbstractMatrix{Complex{T}},
    forward::PreparedLiFTForwardModel,
    amplitude::AbstractMatrix{T}, opd::AbstractMatrix) where {T<:AbstractFloat}
    model = forward.model
    n = size(model.pupil_mask, 1)
    oversampling = lift_oversampling(model.zero_padding)
    n_pad = lift_pad_size(n, model.zero_padding)
    image_size = model.focal_resolution * oversampling
    ws = forward.workspace.propagation
    ensure_psf_buffers!(ws, n_pad)
    if size(dest) != (image_size, image_size)
        throw(DimensionMismatchError("LiFT focal field buffer size must match oversampled image size"))
    end

    opd_to_cycles = T(2) / model.wavelength_m
    ox = cld(n_pad - n, 2)
    oy = cld(n_pad - n, 2)
    fill!(ws.pupil_field, zero(eltype(ws.pupil_field)))
    @views @. ws.pupil_field[ox+1:ox+n, oy+1:oy+n] = amplitude * cispi(opd_to_cycles * opd)
    if iseven(model.focal_resolution)
        phase_shift = -T(pi) * (T(n_pad) + one(T)) / T(n_pad)
        apply_centering_phase!(execution_style(ws.pupil_field), ws.pupil_field, phase_shift)
    end

    copyto!(ws.fft_buffer, ws.pupil_field)
    execute_fft_plan!(ws.fft_buffer, ws.fft_plan)
    ws.fft_buffer ./= T(n_pad)

    shift_pix = if n_pad % 2 == image_size % 2
        0
    elseif iseven(n_pad)
        1
    else
        -1
    end
    start = Int(ceil(n_pad / 2)) - div(image_size, 2) + (1 - (n_pad % 2))
    stop = Int(ceil(n_pad / 2)) + div(image_size, 2) + shift_pix
    @views copyto!(dest, ws.fft_buffer[start:stop, start:stop])
    return oversampling
end

function field_intensity!(dest::AbstractMatrix{T}, field::AbstractMatrix{Complex{T}}, oversampling::Int,
    scratch::AbstractMatrix{T}) where {T<:AbstractFloat}
    if oversampling == 1
        @. dest = abs2(field)
        return dest
    end
    n_out, m_out = size(dest)
    if size(field) != (n_out * oversampling, m_out * oversampling)
        throw(DimensionMismatchError("LiFT field size does not match oversampled rate dimensions"))
    end
    if size(scratch) != size(field)
        throw(DimensionMismatchError("LiFT scratch buffer size must match oversampled field size"))
    end
    @. scratch = abs2(field)
    bin2d!(dest, scratch, oversampling)
    return dest
end

function field_derivative!(dest::AbstractMatrix{T}, buf::AbstractMatrix{Complex{T}},
    Pd::AbstractMatrix{Complex{T}}, oversampling::Int, scale::T, scratch::AbstractMatrix{T}) where {T<:AbstractFloat}
    if size(buf) != size(Pd)
        throw(DimensionMismatchError("LiFT focal fields must have matching sizes"))
    end
    if oversampling == 1
        @. dest = scale * real(im * buf * Pd)
        return dest
    end
    n_out, m_out = size(dest)
    if size(buf) != (n_out * oversampling, m_out * oversampling)
        throw(DimensionMismatchError("LiFT derivative field size does not match oversampled image size"))
    end
    if size(scratch) != size(buf)
        throw(DimensionMismatchError("LiFT scratch buffer size must match oversampled derivative field size"))
    end
    @. scratch = real(im * buf * Pd)
    bin2d!(dest, scratch, oversampling)
    @. dest = scale * dest
    return dest
end

function conjugate_field!(dest::AbstractMatrix{Complex{T}}, src::AbstractMatrix{Complex{T}}) where {T<:AbstractFloat}
    @. dest = conj(src)
    return dest
end
