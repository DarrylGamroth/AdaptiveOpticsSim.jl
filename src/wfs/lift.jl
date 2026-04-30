#
# LiFT phase retrieval
#
# LiFT fits modal coefficients by matching a parameterized PSF model to an
# observed detector image.
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

LiFTLevenbergMarquardt(; lambda0::Real=1e-6, growth::Real=10.0, condition_rtol::Real=sqrt(eps(Float64))) =
    LiFTLevenbergMarquardt(float(lambda0), float(growth), float(condition_rtol))

LiFTAdaptiveLevenbergMarquardt(; lambda0::Real=1e-6, growth::Real=10.0, shrink::Real=2.0,
    min_lambda::Real=1e-10, condition_rtol::Real=sqrt(eps(Float64))) =
    LiFTAdaptiveLevenbergMarquardt(float(lambda0), float(growth), float(shrink), float(min_lambda),
        float(condition_rtol))

struct LiFTParams{T<:AbstractFloat,A<:AbstractMatrix{T},K,S<:LiFTSolveMode,D<:LiFTDampingMode}
    diversity_opd::A
    iterations::Int
    img_resolution::Int
    zero_padding::Int
    object_kernel::K
    solve_mode::S
    damping::D
end

struct LiFTDenseObjectKernel{T<:AbstractFloat,A<:AbstractMatrix{T}}
    kernel::A
end

struct LiFTSeparableObjectKernel{T<:AbstractFloat,V<:AbstractVector{T}}
    row::V
    col::V
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

struct LiFTState{T<:AbstractFloat,
    W<:Workspace,
    B<:AbstractMatrix{T},
    C<:AbstractMatrix{Complex{T}},
    V<:AbstractVector{T},
    I<:AbstractVector{Int}}
    workspace::W
    psf_buffer::B
    amp_buffer::B
    field_scratch::B
    focal_buffer::C
    mode_buffer::C
    pd_buffer::C
    conv_buffer::B
    conv_aux_buffer::B
    opd_buffer::B
    residual_buffer::V
    weight_buffer::V
    H_buffer::B
    normal_buffer::B
    factor_buffer::B
    rhs_buffer::V
    mode_id_buffer::I
    diagnostics::LiFTDiagnostics{T}
end

struct LiFT{M<:LiFTMode,P<:LiFTParams,S<:LiFTState,B<:AbstractArray{<:AbstractFloat,3},SRC<:AbstractSource,D<:AbstractDetector}
    tel::Telescope
    src::SRC
    det::D
    basis::B
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

"""
    LiFT(tel, src, basis, det; ...)

Construct a LiFT phase-retrieval model.

`basis[:, :, k]` defines the modal OPD basis being estimated. The reconstruction
iteratively matches those coefficients to the detector-plane PSF, optionally
including diversity and object convolution.
"""
function LiFT(tel::Telescope, src::AbstractSource, basis::AbstractArray, det::AbstractDetector;
    diversity_opd::AbstractMatrix, iterations::Int=5, img_resolution::Int=0,
    numerical::Bool=false, ang_pixel_arcsec=nothing, object_kernel=nothing,
    solve_mode::LiFTSolveMode=LiFTSolveAuto(), damping::LiFTDampingMode=LiFTDampingNone())
    T = eltype(tel.state.opd)

    if size(basis, 1) != tel.params.resolution || size(basis, 2) != tel.params.resolution
        throw(InvalidConfiguration("basis resolution must match telescope resolution"))
    end
    if size(diversity_opd) != size(tel.state.opd)
        throw(InvalidConfiguration("diversity_opd must match telescope resolution"))
    end

    zero_padding = det.params.psf_sampling
    if ang_pixel_arcsec !== nothing
        scale = psf_pixel_scale_arcsec(tel, src, 1)
        zero_padding = max(1, round(Int, scale / ang_pixel_arcsec))
        @info "LiFT using angular sampling override", zero_padding
    end
    if img_resolution <= 0
        img_resolution = tel.params.resolution * zero_padding
    end
    kernel = object_kernel === nothing ? nothing : _lift_object_kernel(T.(object_kernel))
    params = LiFTParams(float.(diversity_opd), iterations, img_resolution, zero_padding, kernel, solve_mode, damping)
    oversampling = lift_oversampling(zero_padding)
    ws = Workspace(tel.state.opd, lift_pad_size(tel.params.resolution, zero_padding); T=T)
    psf_buffer = similar(tel.state.opd, T, img_resolution, img_resolution)
    amp_buffer = similar(tel.state.opd, T, tel.params.resolution, tel.params.resolution)
    focal_size = img_resolution * oversampling
    field_scratch = similar(psf_buffer, T, focal_size, focal_size)
    focal_buffer = similar(psf_buffer, Complex{eltype(psf_buffer)}, focal_size, focal_size)
    mode_buffer = similar(focal_buffer)
    pd_buffer = similar(focal_buffer)
    conv_buffer = similar(psf_buffer)
    conv_aux_buffer = similar(psf_buffer)
    opd_buffer = similar(amp_buffer)
    residual_buffer = similar(psf_buffer, eltype(psf_buffer), img_resolution * img_resolution)
    weight_buffer = similar(residual_buffer)
    H_buffer = similar(psf_buffer, eltype(psf_buffer), img_resolution * img_resolution, size(basis, 3))
    normal_buffer = similar(psf_buffer, eltype(psf_buffer), size(basis, 3), size(basis, 3))
    factor_buffer = similar(normal_buffer)
    rhs_buffer = similar(residual_buffer, size(basis, 3))
    mode_id_buffer = similar(rhs_buffer, Int, size(basis, 3))
    diagnostics = LiFTDiagnostics(T(NaN), T(NaN), T(NaN), T(NaN), zero(T), false, false)
    state = LiFTState(ws, psf_buffer, amp_buffer, field_scratch, focal_buffer, mode_buffer, pd_buffer, conv_buffer,
        conv_aux_buffer,
        opd_buffer, residual_buffer, weight_buffer, H_buffer, normal_buffer, factor_buffer, rhs_buffer, mode_id_buffer,
        diagnostics)
    mode = numerical ? LiFTNumerical() : LiFTAnalytic()
    return LiFT{typeof(mode), typeof(params), typeof(state), typeof(basis), typeof(src), typeof(det)}(
        tel,
        src,
        det,
        basis,
        params,
        state,
    )
end

diagnostics(lift::LiFT) = lift.state.diagnostics

"""
    prepare_opd!(lift, coeffs)

Assemble the current model OPD from the modal basis coefficients plus the fixed
diversity OPD.
"""
@inline function prepare_opd!(lift::LiFT, coeffs::AbstractVector)
    combine_basis!(lift.state.opd_buffer, lift.basis, coeffs, lift.tel.state.pupil)
    @. lift.state.opd_buffer += lift.params.diversity_opd
    return lift.state.opd_buffer
end

"""
    lift_interaction_matrix!(H, lift, coefficients, mode_ids; flux_norm=1)

Fill the LiFT Jacobian matrix with derivatives of detector intensity with
respect to the requested modal coefficients.

- `LiFTNumerical` uses centered finite differences
- `LiFTAnalytic` uses the field-derivative formulation
"""
function lift_interaction_matrix!(H::AbstractMatrix, lift::LiFT{LiFTNumerical}, coefficients::AbstractVector,
    mode_ids::AbstractVector; flux_norm::Real=1.0)
    T = eltype(lift.state.psf_buffer)
    n_modes = length(mode_ids)
    psf_size = lift.params.img_resolution
    if size(H, 1) < psf_size * psf_size || size(H, 2) < n_modes
        throw(InvalidConfiguration("H buffer size does not match LiFT dimensions"))
    end
    delta = T(1e-9)

    initial_opd = prepare_opd!(lift, coefficients)
    opd_work = lift.state.amp_buffer
    @inbounds for (idx, mode_id) in enumerate(mode_ids)
        @views mode = lift.basis[:, :, mode_id]
        @. opd_work = initial_opd + delta * mode
        psf_p = psf_from_opd!(lift, opd_work; flux_norm=flux_norm)
        copyto!(lift.state.conv_buffer, psf_p)
        @. opd_work = initial_opd - delta * mode
        psf_m = psf_from_opd!(lift, opd_work; flux_norm=flux_norm)
        @. lift.state.psf_buffer = (lift.state.conv_buffer - psf_m) / (2 * delta)
        maybe_object_convolve!(lift, lift.state.psf_buffer)
        @views H[:, idx] .= reshape(lift.state.psf_buffer, :)
    end
    return H
end

function lift_interaction_matrix!(H::AbstractMatrix, lift::LiFT{LiFTAnalytic}, coefficients::AbstractVector,
    mode_ids::AbstractVector; flux_norm::Real=1.0)
    T = eltype(lift.state.psf_buffer)
    n_modes = length(mode_ids)
    psf_size = lift.params.img_resolution
    if size(H, 1) < psf_size * psf_size || size(H, 2) < n_modes
        throw(InvalidConfiguration("H buffer size does not match LiFT dimensions"))
    end

    initial_opd = prepare_opd!(lift, coefficients)

    amp_scale = sqrt(T(photon_flux(lift.src) * flux_norm *
        lift.tel.params.sampling_time * (lift.tel.params.diameter / lift.tel.params.resolution)^2))
    @. lift.state.amp_buffer = ifelse(lift.tel.state.pupil, amp_scale, zero(T))

    oversampling = focal_field_from_opd!(lift.state.focal_buffer, lift, lift.state.amp_buffer, initial_opd)
    conjugate_field!(lift.state.pd_buffer, lift.state.focal_buffer)
    Pd = lift.state.pd_buffer

    k = T(2 * pi) / wavelength(lift.src)
    @inbounds for (idx, mode_id) in enumerate(mode_ids)
        @views mode = lift.basis[:, :, mode_id]
        @. lift.state.amp_buffer = ifelse(lift.tel.state.pupil, amp_scale * mode, zero(T))
        focal_field_from_opd!(lift.state.mode_buffer, lift, lift.state.amp_buffer, initial_opd)
        field_derivative!(lift.state.psf_buffer, lift.state.mode_buffer, Pd, oversampling, 2 * k,
            lift.state.field_scratch)
        maybe_object_convolve!(lift, lift.state.psf_buffer)
        @views H[:, idx] .= reshape(lift.state.psf_buffer, :)
    end
    return H
end

function lift_interaction_matrix(lift::LiFT{LiFTNumerical}, coefficients::AbstractVector,
    mode_ids::AbstractVector; flux_norm::Real=1.0)
    T = eltype(lift.state.psf_buffer)
    psf_size = lift.params.img_resolution
    H = zeros(T, psf_size * psf_size, length(mode_ids))
    return lift_interaction_matrix!(H, lift, coefficients, mode_ids; flux_norm=flux_norm)
end

function lift_interaction_matrix(lift::LiFT{LiFTAnalytic}, coefficients::AbstractVector,
    mode_ids::AbstractVector; flux_norm::Real=1.0)
    T = eltype(lift.state.psf_buffer)
    psf_size = lift.params.img_resolution
    H = zeros(T, psf_size * psf_size, length(mode_ids))
    return lift_interaction_matrix!(H, lift, coefficients, mode_ids; flux_norm=flux_norm)
end

function reconstruct(lift::LiFT, psf_in::AbstractMatrix, mode_ids::AbstractVector;
    coeffs0=nothing, R_n=nothing, optimize_norm::Symbol=:sum, check_convergence::Bool=true)
    T = eltype(lift.state.psf_buffer)
    coeffs = similar(lift.state.rhs_buffer, T, maximum(mode_ids))
    if coeffs0 === nothing
        fill!(coeffs, zero(T))
    else
        copyto!(coeffs, coeffs0)
    end
    reconstruct!(coeffs, lift, psf_in, mode_ids, weight_mode(R_n);
        optimize_norm=optimize_norm, check_convergence=check_convergence)
    mode_ids_buf = @view lift.state.mode_id_buffer[1:length(mode_ids)]
    copyto!(mode_ids_buf, mode_ids)
    out = similar(lift.state.rhs_buffer, T, length(mode_ids))
    gather_selected_coefficients!(execution_style(coeffs), out, coeffs, mode_ids_buf, length(mode_ids))
    return out
end

"""
    reconstruct!(coeffs, lift, psf_in, mode_ids, weighting; ...)

Run the LiFT iterative reconstruction in-place.

Each iteration evaluates the current PSF model, builds/weights the Jacobian,
solves for a modal update, and scatters that update back into the full
coefficient vector.
"""
function reconstruct!(coeffs::AbstractVector{T}, lift::LiFT, psf_in::AbstractMatrix{T}, mode_ids::AbstractVector,
    mode::LiFTWeightingMode; optimize_norm::Symbol=:sum, check_convergence::Bool=true) where {T<:AbstractFloat}
    maximum(mode_ids) <= length(coeffs) ||
        throw(DimensionMismatchError("coefficient buffer length must cover requested mode ids"))
    n_modes = length(mode_ids)

    residual = lift.state.residual_buffer
    sqrtw = lift.state.weight_buffer
    H = @view lift.state.H_buffer[:, 1:n_modes]
    normal = @view lift.state.normal_buffer[1:n_modes, 1:n_modes]
    factor = @view lift.state.factor_buffer[1:n_modes, 1:n_modes]
    rhs = @view lift.state.rhs_buffer[1:n_modes]
    mode_ids_buf = @view lift.state.mode_id_buffer[1:n_modes]
    diag = lift.state.diagnostics
    style = execution_style(H)
    effective_mode = effective_solve_mode(style, lift.params.solve_mode)
    λ_state = initial_damping_state(lift.params.damping, T)
    prev_weighted_residual_norm = T(Inf)
    copyto!(mode_ids_buf, mode_ids)
    init_weights!(sqrtw, mode, psf_in, lift.det)
    for iter in 1:lift.params.iterations
        current_opd = prepare_opd!(lift, coeffs)
        model_psf = psf_from_opd!(lift, current_opd)
        scale = one(T)
        if optimize_norm == :sum
            denom = sum(model_psf)
            if denom > 0
                scale = sum(psf_in) / denom
            end
        elseif optimize_norm == :max
            denom = maximum(model_psf)
            if denom > 0
                scale = maximum(psf_in) / denom
            end
        end
        if scale != one(T)
            model_psf .*= scale
        end
        copyto!(residual, vec(psf_in))
        residual .-= vec(model_psf)
        diag.residual_norm = norm(residual)
        lift_interaction_matrix!(H, lift, coeffs, mode_ids; flux_norm=scale)

        update_weights!(sqrtw, mode, model_psf, lift.det)
        apply_row_weights!(H, sqrtw, n_modes)
        apply_vec_weights!(residual, sqrtw)
        diag.weighted_residual_norm = norm(residual)
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

function reconstruct!(coeffs::AbstractVector, lift::LiFT, psf_in::AbstractMatrix, mode_ids::AbstractVector;
    coeffs0=nothing, R_n=nothing, optimize_norm::Symbol=:sum, check_convergence::Bool=true)
    if coeffs0 === nothing
        fill!(coeffs, zero(eltype(coeffs)))
    else
        copyto!(coeffs, coeffs0)
    end
    return reconstruct!(coeffs, lift, psf_in, mode_ids, weight_mode(R_n);
        optimize_norm=optimize_norm, check_convergence=check_convergence)
end

@inline function apply_vec_weights!(vec::AbstractVector{T}, weights::AbstractVector{T}) where {T<:AbstractFloat}
    @. vec *= weights
    return vec
end

@inline function apply_row_weights!(mat::AbstractMatrix{T}, weights::AbstractVector{T}, n_cols::Int) where {T<:AbstractFloat}
    row_weights = reshape(weights, :, 1)
    @views @. mat[:, 1:n_cols] *= row_weights
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
        denom = s[i]^2 + λ^2
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
        @views factor[diagind(factor)] .+= λ
        chol = cholesky!(Hermitian(factor), check=false)
        if !issuccess(chol)
            λ *= T(10)
            @views factor[diagind(factor)] .+= λ
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
        @views factor[diagind(factor)] .+= λ
    end
    chol = cholesky!(Hermitian(factor), check=false)
    while !issuccess(chol)
        λ = max(λ * T(damping.growth), regularization_load(normal))
        copyto!(factor, normal)
        @views factor[diagind(factor)] .+= λ
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
    synchronize_backend!(style)
    return coeffs
end

@inline gather_selected_coefficients!(::ScalarCPUStyle, out::AbstractVector, coeffs::AbstractVector,
    mode_ids::AbstractVector, n_modes::Int) = copyto!(out, coeffs[mode_ids])

@inline function gather_selected_coefficients!(style::AcceleratorStyle, out::AbstractVector, coeffs::AbstractVector,
    mode_ids::AbstractVector, n_modes::Int)
    launch_kernel!(style, lift_gather_kernel!, out, coeffs, mode_ids, n_modes; ndrange=n_modes)
    synchronize_backend!(style)
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
    current_residual::T, normal::AbstractMatrix{T}) where {T<:AbstractFloat}
    base = regularization_load(normal)
    if !isfinite(prev_residual)
        return max(T(damping.lambda0), T(damping.min_lambda), base)
    elseif current_residual > prev_residual * (one(T) + sqrt(eps(T)))
        return max(λ * T(damping.growth), T(damping.min_lambda), base)
    end
    return max(λ / T(damping.shrink), T(damping.min_lambda), base)
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

function normal_condition_ratio(normal::AbstractMatrix{T}) where {T<:AbstractFloat}
    diagvals = abs.(@view normal[diagind(normal)])
    maxabs = maximum(diagvals)
    minabs = minimum(diagvals)
    maxabs == zero(T) && return T(Inf)
    return maxabs / max(minabs, eps(T))
end

@inline function init_weights!(sqrtw::AbstractVector{T}, ::LiFTWeightingDynamic,
    ::AbstractMatrix{T}, ::AbstractDetector) where {T<:AbstractFloat}
    return sqrtw
end

@inline function init_weights!(sqrtw::AbstractVector{T}, mode::LiFTWeightingStatic,
    psf_in::AbstractMatrix{T}, det::AbstractDetector) where {T<:AbstractFloat}
    weight_vector!(sqrtw, psf_in, mode, det)
    map!(sqrt, sqrtw, sqrtw)
    return sqrtw
end

@inline function update_weights!(sqrtw::AbstractVector{T}, ::LiFTWeightingStatic,
    ::AbstractMatrix{T}, ::AbstractDetector) where {T<:AbstractFloat}
    return sqrtw
end

@inline function update_weights!(sqrtw::AbstractVector{T}, mode::LiFTWeightingDynamic,
    model_psf::AbstractMatrix{T}, det::AbstractDetector) where {T<:AbstractFloat}
    weight_vector!(sqrtw, model_psf, mode, det)
    map!(sqrt, sqrtw, sqrtw)
    return sqrtw
end

function weight_vector!(out::AbstractVector{T}, psf::AbstractMatrix{T}, ::LiFTWeightModel,
    det::AbstractDetector) where {T<:AbstractFloat}
    σ = readout_noise(det)
    σ2 = T(σ * σ)
    vals = vec(psf)
    @. out = inv(max(vals + σ2, eps(T)))
    return out
end

function weight_vector!(out::AbstractVector{T}, psf::AbstractMatrix{T}, ::LiFTWeightIterative,
    det::AbstractDetector) where {T<:AbstractFloat}
    return weight_vector!(out, psf, LiFTWeightModel(), det)
end

function weight_vector!(out::AbstractVector{T}, ::AbstractMatrix{T}, ::LiFTWeightNone,
    det::AbstractDetector) where {T<:AbstractFloat}
    σ = readout_noise(det)
    σ2 = T(σ * σ)
    fill!(out, T(1) / max(σ2, eps(T)))
    return out
end

function weight_vector!(out::AbstractVector{T}, ::AbstractMatrix{T}, mode::LiFTWeightMatrix,
    ::AbstractDetector) where {T<:AbstractFloat}
    vals = vec(mode.R_n)
    @. out = inv(max(vals, eps(T)))
    return out
end

function weight_vector!(out::AbstractVector{T}, psf::AbstractMatrix{T}, R_n,
    det::AbstractDetector) where {T<:AbstractFloat}
    return weight_vector!(out, psf, weight_mode(R_n), det)
end

function weight_vector(psf::AbstractMatrix{T}, R_n, det::AbstractDetector) where {T<:AbstractFloat}
    out = similar(psf, T, length(psf))
    return weight_vector!(out, psf, weight_mode(R_n), det)
end

"""
    psf_from_opd!(lift, opd; flux_norm=1)

Evaluate the LiFT forward model PSF from an OPD map.

This includes coherent propagation, intensity formation, detector-size
resampling, and optional object convolution.
"""
function psf_from_opd!(lift::LiFT, opd::AbstractMatrix; flux_norm::Real=1.0)
    amp_scale = sqrt(eltype(lift.state.psf_buffer)(photon_flux(lift.src) * flux_norm *
        lift.tel.params.sampling_time * (lift.tel.params.diameter / lift.tel.params.resolution)^2))
    @. lift.state.amp_buffer = ifelse(lift.tel.state.pupil, amp_scale, zero(eltype(lift.state.amp_buffer)))
    oversampling = focal_field_from_opd!(lift.state.focal_buffer, lift, lift.state.amp_buffer, opd)
    field_intensity!(lift.state.psf_buffer, lift.state.focal_buffer, oversampling, lift.state.field_scratch)
    maybe_object_convolve!(lift, lift.state.psf_buffer)
    return lift.state.psf_buffer
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

function readout_noise(det::Detector{<:NoiseNone})
    return 0.0
end

function readout_noise(det::Detector{<:NoisePhoton})
    return 0.0
end

function readout_noise(det::Detector{<:NoiseReadout})
    return det.noise.sigma
end

function readout_noise(det::Detector{<:NoisePhotonReadout})
    return det.noise.sigma
end

function maybe_object_convolve!(lift::LiFT{<:LiFTMode,<:LiFTParams{<:AbstractFloat,<:AbstractMatrix,Nothing}}, mat::AbstractMatrix)
    return mat
end

function maybe_object_convolve!(lift::LiFT{<:LiFTMode,<:LiFTParams{<:AbstractFloat,<:AbstractMatrix,<:LiFTDenseObjectKernel}}, mat::AbstractMatrix)
    conv2d_same!(lift.state.conv_buffer, mat, lift.params.object_kernel.kernel)
    copyto!(mat, lift.state.conv_buffer)
    return mat
end

function maybe_object_convolve!(lift::LiFT{<:LiFTMode,<:LiFTParams{<:AbstractFloat,<:AbstractMatrix,<:LiFTSeparableObjectKernel}}, mat::AbstractMatrix)
    conv2d_same_separable!(lift.state.conv_buffer, lift.state.conv_aux_buffer, mat,
        lift.params.object_kernel.row, lift.params.object_kernel.col)
    copyto!(mat, lift.state.conv_buffer)
    return mat
end

function _lift_object_kernel(kernel::AbstractMatrix{T}) where {T<:AbstractFloat}
    row, col = _separable_kernel_factors(kernel)
    if row === nothing
        return LiFTDenseObjectKernel{T,typeof(kernel)}(kernel)
    end
    return LiFTSeparableObjectKernel{T,typeof(row)}(row, col)
end

function _separable_kernel_factors(kernel::AbstractMatrix{T}) where {T<:AbstractFloat}
    F = svd(Matrix(kernel); full=false)
    isempty(F.S) && return nothing, nothing
    σ1 = F.S[1]
    σ2 = length(F.S) >= 2 ? F.S[2] : zero(T)
    tol = sqrt(eps(T)) * max(one(T), σ1)
    σ2 <= tol || return nothing, nothing
    scale = sqrt(σ1)
    row = T.(F.U[:, 1] .* scale)
    col = T.(F.V[:, 1] .* scale)
    return row, col
end

function conv2d_same!(dest::AbstractMatrix{T}, src::AbstractMatrix{T}, kernel::AbstractMatrix) where {T<:AbstractFloat}
    n, m = size(src)
    kh, kw = size(kernel)
    cx = div(kh, 2)
    cy = div(kw, 2)
    norm = sum(kernel)
    inv_norm = norm == 0 ? one(T) : T(1) / T(norm)
    @inbounds for i in 1:n, j in 1:m
        acc = zero(T)
        for ki in 1:kh, kj in 1:kw
            ii = symm_index(i + ki - cx - 1, n)
            jj = symm_index(j + kj - cy - 1, m)
            acc += src[ii, jj] * kernel[ki, kj]
        end
        dest[i, j] = acc * inv_norm
    end
    return dest
end

function conv2d_same_separable!(dest::AbstractMatrix{T}, tmp::AbstractMatrix{T}, src::AbstractMatrix{T},
    row_kernel::AbstractVector{T}, col_kernel::AbstractVector{T}) where {T<:AbstractFloat}
    n, m = size(src)
    kr = length(row_kernel)
    kc = length(col_kernel)
    cx = div(kr, 2)
    cy = div(kc, 2)
    row_norm = sum(row_kernel)
    col_norm = sum(col_kernel)
    inv_norm = (iszero(row_norm) || iszero(col_norm)) ? one(T) : inv(row_norm * col_norm)
    @inbounds for i in 1:n, j in 1:m
        acc = zero(T)
        for ki in 1:kr
            ii = symm_index(i + ki - cx - 1, n)
            acc += src[ii, j] * row_kernel[ki]
        end
        tmp[i, j] = acc
    end
    @inbounds for i in 1:n, j in 1:m
        acc = zero(T)
        for kj in 1:kc
            jj = symm_index(j + kj - cy - 1, m)
            acc += tmp[i, jj] * col_kernel[kj]
        end
        dest[i, j] = acc * inv_norm
    end
    return dest
end

@inline function symm_index(i::Int, n::Int)
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

function focal_field_from_opd!(dest::AbstractMatrix{Complex{T}}, lift::LiFT,
    amplitude::AbstractMatrix{T}, opd::AbstractMatrix) where {T<:AbstractFloat}
    n = lift.tel.params.resolution
    oversampling = lift_oversampling(lift.params.zero_padding)
    n_pad = lift_pad_size(n, lift.params.zero_padding)
    img_size = lift.params.img_resolution * oversampling
    ws = lift.state.workspace
    ensure_psf_buffers!(ws, n_pad)
    if size(dest) != (img_size, img_size)
        throw(DimensionMismatchError("LiFT focal field buffer size must match oversampled image size"))
    end

    opd_to_cycles = T(2) / wavelength(lift.src)
    ox = cld(n_pad - n, 2)
    oy = cld(n_pad - n, 2)
    fill!(ws.pupil_field, zero(eltype(ws.pupil_field)))
    @views @. ws.pupil_field[ox+1:ox+n, oy+1:oy+n] = amplitude * cispi(opd_to_cycles * opd)
    if iseven(lift.params.img_resolution)
        phase_shift = -T(pi) * (T(n_pad) + one(T)) / T(n_pad)
        apply_centering_phase!(execution_style(ws.pupil_field), ws.pupil_field, phase_shift)
    end

    copyto!(ws.fft_buffer, ws.pupil_field)
    execute_fft_plan!(ws.fft_buffer, ws.fft_plan)
    ws.fft_buffer ./= T(n_pad)

    shift_pix = if n_pad % 2 == img_size % 2
        0
    elseif iseven(n_pad)
        1
    else
        -1
    end
    start = Int(ceil(n_pad / 2)) - div(img_size, 2) + (1 - (n_pad % 2))
    stop = Int(ceil(n_pad / 2)) + div(img_size, 2) + shift_pix
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
        throw(DimensionMismatchError("LiFT field size does not match oversampled PSF dimensions"))
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
