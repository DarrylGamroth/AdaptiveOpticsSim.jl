@kernel function lift_basis_expansion_kernel!(opd, basis, coeffs, pupil, n_modes::Int)
    I = @index(Global, Cartesian)
    i, j = Tuple(I)
    if i <= size(opd, 1) && j <= size(opd, 2)
        value = zero(eltype(opd))
        @inbounds for k in 1:n_modes
            value += coeffs[k] * basis[i, j, k]
        end
        @inbounds opd[i, j] = ifelse(pupil[i, j], value, zero(value))
    end
end

@inline function lift_basis_expansion!(opd::AbstractMatrix{T},
    basis::AbstractArray{T,3}, coeffs::AbstractVector{T},
    pupil::AbstractMatrix{Bool}) where {T<:AbstractFloat}
    return lift_basis_expansion!(execution_style(opd), opd, basis, coeffs,
        pupil)
end

function lift_basis_expansion!(::ScalarCPUStyle, opd::AbstractMatrix{T},
    basis::AbstractArray{T,3}, coeffs::AbstractVector{T},
    pupil::AbstractMatrix{Bool}) where {T<:AbstractFloat}
    n_modes = min(size(basis, 3), length(coeffs))
    fill!(opd, zero(T))
    @inbounds for k in 1:n_modes
        @views @. opd += coeffs[k] * basis[:, :, k]
    end
    @. opd *= pupil
    return opd
end

function lift_basis_expansion!(style::AcceleratorStyle,
    opd::AbstractMatrix{T}, basis::AbstractArray{T,3},
    coeffs::AbstractVector{T}, pupil::AbstractMatrix{Bool}) where {T<:AbstractFloat}
    n_modes = min(size(basis, 3), length(coeffs))
    if iszero(n_modes)
        fill!(opd, zero(T))
        return opd
    end
    launch_kernel!(style, lift_basis_expansion_kernel!, opd, basis, coeffs,
        pupil, n_modes; ndrange=size(opd))
    return opd
end

"""
    prepare_opd!(lift, coeffs)

Assemble the current model OPD from the modal basis coefficients plus the fixed
diversity OPD.
"""
@inline function prepare_opd!(lift::PreparedLiFTEstimator, coeffs::AbstractVector)
    model = _lift_model(lift)
    input = lift.forward.input
    lift_basis_expansion!(input, model.basis, coeffs,
        model.pupil_mask)
    @. input += model.diversity_opd
    return input
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

"""
    lift_interaction_matrix!(H, lift, coefficients; rate_scale=1)

Fill the LiFT Jacobian with derivatives of the prepared observation-domain
photon-arrival rate with respect to the requested modal coefficients.

- `LiFTNumericalJacobian` uses centered finite differences
- `LiFTAnalyticJacobian` uses the field-derivative formulation
"""
function lift_interaction_matrix!(H::AbstractMatrix,
    lift::PreparedLiFTEstimator, coefficients::AbstractVector;
    rate_scale::Real=1.0)
    _require_lift_interaction_inputs(H, lift, coefficients, rate_scale)
    return _lift_interaction_matrix!(lift.plan.jacobian_method, H, lift,
        coefficients; rate_scale)
end

function _require_lift_interaction_inputs(H::AbstractMatrix,
    lift::PreparedLiFTEstimator, coefficients::AbstractVector,
    rate_scale::Real)
    _require_lift_estimator(lift)
    output = lift.forward.output.values
    T = eltype(output)
    expected_size = (length(output), length(lift.plan.mode_ids))
    size(H) == expected_size || throw(DimensionMismatchError(
        "LiFT Jacobian destination must have size $expected_size"))
    eltype(H) === T || throw(InvalidConfiguration(
        "LiFT Jacobian destination must use the prepared numeric type"))
    typeof(backend(H)) === typeof(lift.backend) || throw(
        InvalidConfiguration(
            "LiFT Jacobian destination must use the prepared array backend"))
    compute_device(H) == lift.device || throw(InvalidConfiguration(
        "LiFT Jacobian destination must occupy the prepared compute device"))
    length(coefficients) == size(lift.forward.plan.basis, 3) || throw(
        DimensionMismatchError(
            "LiFT Jacobian coefficients must cover the complete prepared basis"))
    eltype(coefficients) === T || throw(InvalidConfiguration(
        "LiFT Jacobian coefficients must use the prepared numeric type"))
    typeof(backend(coefficients)) === typeof(lift.backend) || throw(
        InvalidConfiguration(
            "LiFT Jacobian coefficients must use the prepared array backend"))
    compute_device(coefficients) == lift.device || throw(
        InvalidConfiguration(
            "LiFT Jacobian coefficients must occupy the prepared compute device"))
    _wfs_storage_mightalias(H, coefficients) && throw(InvalidConfiguration(
        "LiFT Jacobian destination and coefficient input must not alias"))
    scale = T(rate_scale)
    isfinite(scale) && scale >= zero(T) || throw(InvalidConfiguration(
        "LiFT Jacobian rate_scale must be finite and nonnegative"))
    return nothing
end

function _lift_interaction_matrix!(::LiFTNumericalJacobian,
    H::AbstractMatrix, lift::PreparedLiFTEstimator,
    coefficients::AbstractVector;
    rate_scale::Real=1.0)
    model = _lift_model(lift)
    workspace = _lift_workspace(lift)
    T = eltype(workspace.optical_rate_buffer)
    mode_ids = lift.plan.mode_ids
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

function _lift_interaction_matrix!(::LiFTAnalyticJacobian,
    H::AbstractMatrix, lift::PreparedLiFTEstimator,
    coefficients::AbstractVector;
    rate_scale::Real=1.0)
    model = _lift_model(lift)
    workspace = _lift_workspace(lift)
    T = eltype(workspace.optical_rate_buffer)
    mode_ids = lift.plan.mode_ids
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

function lift_interaction_matrix(lift::PreparedLiFTEstimator,
    coefficients::AbstractVector; rate_scale::Real=1.0)
    output = lift.forward.output.values
    H = similar(output, eltype(output), length(output),
        length(lift.plan.mode_ids))
    return lift_interaction_matrix!(H, lift, coefficients;
        rate_scale=rate_scale)
end

@inline _require_lift_weighting(::LiFTWeightingMode,
    ::LiFTObservationMetadata, ::Type{<:AbstractFloat}) = nothing

function _require_lift_weighting(mode::LiFTVarianceMapWeighting,
    metadata::LiFTObservationMetadata,
    ::Type{T}) where {T<:AbstractFloat}
    size(mode.variance) == metadata.contract.rate_metadata.dimensions || throw(
        DimensionMismatchError(
            "LiFT variance map must match observation dimensions"))
    eltype(mode.variance) === T || throw(InvalidConfiguration(
        "LiFT variance map must use the estimator numeric type"))
    typeof(backend(mode.variance)) === typeof(metadata.backend) || throw(
        InvalidConfiguration(
            "LiFT variance map must use the observation array backend"))
    compute_device(mode.variance) == metadata.device || throw(InvalidConfiguration(
        "LiFT variance map must occupy the observation compute device"))
    return nothing
end

function _copy_lift_observation_rate!(dest::AbstractMatrix{T},
    observation::LiFTObservation) where {T<:AbstractFloat}
    scale = lift_observation_to_rate_scale(observation.metadata.domain, T)
    @. dest = observation.values * scale
    return dest
end

function scale_lift_model!(::LiFTTotalRateMatching,
    model_rate::AbstractMatrix{T}, observation_rate::AbstractMatrix{T}) where {T}
    denominator = backend_sum_value(model_rate)
    scale = denominator > zero(T) ?
        backend_sum_value(observation_rate) / denominator : one(T)
    scale == one(T) || (model_rate .*= scale)
    return scale, abs(denominator * scale)
end

function scale_lift_model!(::LiFTPeakRateMatching,
    model_rate::AbstractMatrix{T}, observation_rate::AbstractMatrix{T}) where {T}
    denominator = backend_maximum_value(model_rate)
    scale = denominator > zero(T) ?
        backend_maximum_value(observation_rate) / denominator : one(T)
    scale == one(T) || (model_rate .*= scale)
    return scale, abs(backend_sum_value(model_rate))
end

function scale_lift_model!(::LiFTPhysicalRatePreservation,
    model_rate::AbstractMatrix{T}, ::AbstractMatrix{T}) where {T}
    return one(T), abs(backend_sum_value(model_rate))
end

function reconstruct(lift::PreparedLiFTEstimator)
    reconstruct!(lift)
    return copy(lift.coefficients)
end

@inline function _initialize_lift_coefficients!(dest::AbstractVector,
    ::Nothing)
    fill!(dest, zero(eltype(dest)))
    return dest
end

@inline function _initialize_lift_coefficients!(dest::AbstractVector,
    initial::AbstractVector)
    copyto!(dest, initial)
    return dest
end

"""
    reconstruct!(lift)

Run the PreparedLiFTEstimator iterative reconstruction in-place.

Each iteration evaluates the current rate model, builds/weights the Jacobian,
solves for a modal update, and scatters that update back into the full
coefficient vector.
"""
function reconstruct!(lift::PreparedLiFTEstimator)
    _require_lift_estimator(lift)
    coeffs = lift.workspace.full_coefficients_buffer
    T = eltype(coeffs)
    _initialize_lift_coefficients!(coeffs, lift.initial_coefficients)
    mode = lift.plan.weighting
    observation = lift.observation
    n_modes = length(lift.plan.mode_ids)
    observation_rate = _copy_lift_observation_rate!(
        lift.workspace.observation_rate_buffer, observation)

    residual = lift.workspace.residual_buffer
    sqrtw = lift.workspace.weight_buffer
    H = lift.workspace.H_buffer
    normal = lift.workspace.normal_buffer
    factor = lift.workspace.factor_buffer
    rhs = lift.workspace.rhs_buffer
    mode_ids_buf = lift.plan.mode_ids_device
    diag = lift.workspace.diagnostics
    style = execution_style(H)
    effective_mode = effective_solve_mode(style, lift.plan.solve_mode)
    λ_state = initial_damping_state(lift.plan.damping, T)
    prev_weighted_residual_norm = T(Inf)
    init_weights!(sqrtw, mode, observation_rate, observation.metadata)
    for iter in 1:lift.plan.iterations
        current_opd = prepare_opd!(lift, coeffs)
        model_rate = _lift_rate_values_from_opd!(lift.forward, current_opd)
        scale, model_photon_rate = scale_lift_model!(
            lift.plan.model_scaling, model_rate, observation_rate)
        objective_scale = max(model_photon_rate, one(T))
        lift_residual!(residual, observation_rate, model_rate)
        diag.residual_norm = norm(residual)
        update_weights!(sqrtw, mode, model_rate, observation.metadata, iter)
        # Generate the Jacobian in the normalized objective scale so a
        # physically large photon rate cannot overflow before H'H is formed.
        _lift_interaction_matrix!(lift.plan.jacobian_method, H, lift, coeffs;
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
        λ_state = update_damping_state(lift.plan.damping, λ_state, prev_weighted_residual_norm,
            diag.weighted_residual_norm, normal)
        damping = effective_damping(lift.plan.damping, λ_state)
        delta = solve_lift_system!(diag, residual, rhs, H, normal, factor, effective_mode, damping)
        diag.update_norm = delta_norm(delta, n_modes, effective_mode)
        update_coefficients!(coeffs, delta, mode_ids_buf, effective_mode, style)
        prev_weighted_residual_norm = diag.weighted_residual_norm
        if lift.plan.check_convergence &&
                diag.update_norm / max(norm(coeffs), eps(T)) < 1e-3
            break
        end
    end
    gather_selected_coefficients!(style, lift.coefficients, coeffs,
        mode_ids_buf, n_modes)
    return lift.coefficients
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

function solve_lift_fallback!(diag::LiFTDiagnosticsWorkspace{T}, rhs::AbstractVector{T},
    H::AbstractMatrix{T}, residual::AbstractVector{T}, damping::LiFTDampingMode) where {T<:AbstractFloat}
    λ = fallback_damping_lambda(damping, T, H)
    F = svd(H; full=false)
    s = F.S
    work = similar(rhs)
    mul!(work, transpose(F.U), residual)
    @. work = ifelse(iszero(s^2 + λ), zero(T),
        (s * work) / (s^2 + λ))
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
function solve_normal_system!(diag::LiFTDiagnosticsWorkspace{T}, rhs::AbstractVector{T}, factor::AbstractMatrix{T},
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

function solve_normal_system!(diag::LiFTDiagnosticsWorkspace{T}, rhs::AbstractVector{T}, factor::AbstractMatrix{T},
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

function solve_lift_system!(diag::LiFTDiagnosticsWorkspace{T}, residual::AbstractVector{T}, rhs::AbstractVector{T},
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

function solve_lift_system!(diag::LiFTDiagnosticsWorkspace{T}, residual::AbstractVector{T}, rhs::AbstractVector{T},
    H::AbstractMatrix{T}, normal::AbstractMatrix{T}, factor::AbstractMatrix{T},
    effective_mode::LiFTSolveMode, damping::LiFTDampingMode) where {T<:AbstractFloat}
    diag.used_qr = false
    diag.used_fallback = false
    diag.regularization = zero(T)
    diag.condition_ratio = normal_condition_ratio(normal)
    solve_normal_system!(diag, rhs, factor, normal, H, residual, damping)
    return rhs
end

function handle_lift_qr_error!(::SingularException, diag::LiFTDiagnosticsWorkspace, rhs::AbstractVector,
    factor::AbstractMatrix, normal::AbstractMatrix, H::AbstractMatrix, residual::AbstractVector,
    damping::LiFTDampingMode)
    diag.used_qr = false
    diag.used_fallback = true
    solve_normal_system!(diag, rhs, factor, normal, H, residual, damping)
    return rhs
end

function handle_lift_qr_error!(err, diag::LiFTDiagnosticsWorkspace, rhs::AbstractVector,
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

@inline regularization_load(normal::AbstractMatrix{T}) where {T<:AbstractFloat} =
    regularization_load(reduction_execution_strategy(normal), normal)

function regularization_load(
    ::DirectReductionStrategy,
    normal::AbstractMatrix{T},
) where {T<:AbstractFloat}
    n = min(size(normal, 1), size(normal, 2))
    diagvals = @view normal[diagind(normal)]
    diag_sum = sum(abs, diagvals)
    mean_diag = n == 0 ? zero(T) : diag_sum / T(n)
    return max(sqrt(eps(T)) * max(mean_diag, one(T)), eps(T))
end

function regularization_load(
    ::HostMirrorReductionStrategy,
    normal::AbstractMatrix{T},
) where {T<:AbstractFloat}
    host_parent = Array(reduction_parent_source(normal))
    host_normal = reduction_host_view(host_parent, normal)
    return regularization_load(DirectReductionStrategy(), host_normal)
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
    normal_condition_ratio(reduction_execution_strategy(normal), normal)

function normal_condition_ratio(::DirectReductionStrategy,
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


function normal_condition_ratio(::DirectReductionStrategy,
    normal::AbstractMatrix{T}) where {T<:AbstractFloat}
    diagonal = @view normal[diagind(normal)]
    maxabs = maximum(abs, diagonal)
    minabs = minimum(abs, diagonal)
    maxabs == zero(T) && return T(Inf)
    return maxabs / max(minabs, eps(T))
end

function normal_condition_ratio(::HostMirrorReductionStrategy, normal::AbstractMatrix{T}) where {T<:AbstractFloat}
    host_parent = Array(reduction_parent_source(normal))
    host_normal = reduction_host_view(host_parent, normal)
    return normal_condition_ratio(DirectReductionStrategy(), host_normal)
end

@inline function init_weights!(sqrtw::AbstractVector{T},
    ::Union{LiFTInitialModelWeighting,LiFTIterativeModelWeighting},
    ::AbstractMatrix{T}, ::LiFTObservationMetadata) where {T<:AbstractFloat}
    return sqrtw
end

@inline function init_weights!(sqrtw::AbstractVector{T},
    mode::Union{LiFTReadNoiseWeighting,LiFTVarianceMapWeighting},
    observation_rate::AbstractMatrix{T},
    metadata::LiFTObservationMetadata) where {T<:AbstractFloat}
    weight_vector!(sqrtw, observation_rate, mode, metadata)
    sqrt_weights!(execution_style(sqrtw), sqrtw)
    return sqrtw
end

@inline function update_weights!(sqrtw::AbstractVector{T},
    ::Union{LiFTReadNoiseWeighting,LiFTVarianceMapWeighting},
    ::AbstractMatrix{T}, ::LiFTObservationMetadata,
    ::Int) where {T<:AbstractFloat}
    return sqrtw
end

@inline function update_weights!(sqrtw::AbstractVector{T},
    mode::LiFTInitialModelWeighting,
    model_rate::AbstractMatrix{T}, metadata::LiFTObservationMetadata,
    iteration::Int) where {T<:AbstractFloat}
    isone(iteration) || return sqrtw
    weight_vector!(sqrtw, model_rate, mode, metadata)
    sqrt_weights!(execution_style(sqrtw), sqrtw)
    return sqrtw
end

@inline function update_weights!(sqrtw::AbstractVector{T},
    mode::LiFTIterativeModelWeighting,
    model_rate::AbstractMatrix{T},
    metadata::LiFTObservationMetadata, ::Int) where {T<:AbstractFloat}
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
    model_rate::AbstractMatrix{T}, ::LiFTInitialModelWeighting,
    metadata::LiFTObservationMetadata) where {T<:AbstractFloat}
    readout_variance = _lift_readout_variance_rate(metadata, T)
    shot_scale = lift_shot_variance_rate_scale(metadata.domain, T)
    return lift_inverse_variance!(out, model_rate, shot_scale,
        readout_variance)
end

function weight_vector!(out::AbstractVector{T},
    model_rate::AbstractMatrix{T}, ::LiFTIterativeModelWeighting,
    metadata::LiFTObservationMetadata) where {T<:AbstractFloat}
    return weight_vector!(out, model_rate, LiFTInitialModelWeighting(), metadata)
end

function weight_vector!(out::AbstractVector{T}, ::AbstractMatrix{T},
    ::LiFTReadNoiseWeighting,
    metadata::LiFTObservationMetadata) where {T<:AbstractFloat}
    readout_variance = _lift_readout_variance_rate(metadata, T)
    fill!(out, one(T) / max(readout_variance, eps(T)))
    return out
end

function weight_vector!(out::AbstractVector{T}, ::AbstractMatrix{T},
    mode::LiFTVarianceMapWeighting,
    metadata::LiFTObservationMetadata) where {T<:AbstractFloat}
    variance_scale = lift_observation_to_rate_scale(metadata.domain, T)^2
    return lift_inverse_variance!(out, mode.variance, variance_scale, zero(T))
end
