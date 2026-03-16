using LinearAlgebra
using Logging

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

LiFTLevenbergMarquardt(; lambda0::Real=1e-6, growth::Real=10.0, condition_rtol::Real=sqrt(eps(Float64))) =
    LiFTLevenbergMarquardt(float(lambda0), float(growth), float(condition_rtol))

struct LiFTParams{T<:AbstractFloat,A<:AbstractMatrix{T},K,S<:LiFTSolveMode,D<:LiFTDampingMode}
    diversity_opd::A
    iterations::Int
    img_resolution::Int
    zero_padding::Int
    object_kernel::K
    solve_mode::S
    damping::D
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

mutable struct LiFTState{T<:AbstractFloat,
    W<:Workspace,
    B<:AbstractMatrix{T},
    C<:AbstractMatrix{Complex{T}},
    V<:AbstractVector{T}}
    workspace::W
    psf_buffer::B
    amp_buffer::B
    focal_buffer::C
    mode_buffer::C
    pd_buffer::C
    conv_buffer::B
    opd_buffer::B
    residual_buffer::V
    weight_buffer::V
    H_buffer::B
    normal_buffer::B
    factor_buffer::B
    rhs_buffer::V
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
    kernel = object_kernel === nothing ? nothing : T.(object_kernel)
    params = LiFTParams(float.(diversity_opd), iterations, img_resolution, zero_padding, kernel, solve_mode, damping)
    oversampling = lift_oversampling(zero_padding)
    ws = Workspace(tel.state.opd, lift_pad_size(tel.params.resolution, zero_padding); T=T)
    psf_buffer = similar(tel.state.opd, T, img_resolution, img_resolution)
    amp_buffer = similar(tel.state.opd, T, tel.params.resolution, tel.params.resolution)
    focal_buffer = similar(psf_buffer, Complex{eltype(psf_buffer)}, img_resolution * oversampling, img_resolution * oversampling)
    mode_buffer = similar(focal_buffer)
    pd_buffer = similar(focal_buffer)
    conv_buffer = similar(psf_buffer)
    opd_buffer = similar(amp_buffer)
    residual_buffer = similar(psf_buffer, eltype(psf_buffer), img_resolution * img_resolution)
    weight_buffer = similar(residual_buffer)
    H_buffer = similar(psf_buffer, eltype(psf_buffer), img_resolution * img_resolution, size(basis, 3))
    normal_buffer = similar(psf_buffer, eltype(psf_buffer), size(basis, 3), size(basis, 3))
    factor_buffer = similar(normal_buffer)
    rhs_buffer = similar(residual_buffer, size(basis, 3))
    diagnostics = LiFTDiagnostics(T(NaN), T(NaN), T(NaN), T(NaN), zero(T), false, false)
    state = LiFTState(ws, psf_buffer, amp_buffer, focal_buffer, mode_buffer, pd_buffer, conv_buffer,
        opd_buffer, residual_buffer, weight_buffer, H_buffer, normal_buffer, factor_buffer, rhs_buffer, diagnostics)
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

@inline function prepare_opd!(lift::LiFT, coeffs::AbstractVector)
    combine_basis!(lift.state.opd_buffer, lift.basis, coeffs, lift.tel.state.pupil)
    @. lift.state.opd_buffer += lift.params.diversity_opd
    return lift.state.opd_buffer
end

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
        field_derivative!(lift.state.psf_buffer, lift.state.mode_buffer, Pd, oversampling, 2 * k)
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
    coeffs = coeffs0 === nothing ? zeros(T, maximum(mode_ids)) : copy(coeffs0)
    reconstruct!(coeffs, lift, psf_in, mode_ids, weight_mode(R_n);
        optimize_norm=optimize_norm, check_convergence=check_convergence)
    return coeffs[mode_ids]
end

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
    diag = lift.state.diagnostics
    style = execution_style(H)
    effective_mode = effective_solve_mode(style, lift.params.solve_mode)
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
        @inbounds for idx in eachindex(model_psf)
            residual[idx] = psf_in[idx] - model_psf[idx]
        end
        diag.residual_norm = norm(residual)
        lift_interaction_matrix!(H, lift, coeffs, mode_ids; flux_norm=scale)

        update_weights!(sqrtw, mode, model_psf, lift.det)
        apply_row_weights!(H, sqrtw, n_modes)
        apply_vec_weights!(residual, sqrtw)
        diag.weighted_residual_norm = norm(residual)
        mul!(normal, adjoint(H), H)
        mul!(rhs, adjoint(H), residual)
        delta = solve_lift_system!(diag, residual, rhs, H, normal, factor, effective_mode, lift.params.damping)
        diag.update_norm = delta_norm(delta, n_modes, effective_mode)
        for (j, mode_id) in enumerate(mode_ids)
            coeffs[mode_id] += delta_component(delta, j, effective_mode)
        end
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
    @inbounds for idx in eachindex(vec, weights)
        vec[idx] *= weights[idx]
    end
    return vec
end

@inline function apply_row_weights!(mat::AbstractMatrix{T}, weights::AbstractVector{T}, n_cols::Int) where {T<:AbstractFloat}
    @inbounds for col in 1:n_cols
        for row in eachindex(weights)
            mat[row, col] *= weights[row]
        end
    end
    return mat
end

@inline effective_solve_mode(::ScalarCPUStyle, ::LiFTSolveAuto) = LiFTSolveQR()
@inline effective_solve_mode(::AcceleratorStyle, ::LiFTSolveAuto) = LiFTSolveNormalEquations()
@inline effective_solve_mode(::ExecutionStyle, mode::LiFTSolveMode) = mode

function solve_normal_system!(diag::LiFTDiagnostics{T}, rhs::AbstractVector{T}, factor::AbstractMatrix{T},
    normal::AbstractMatrix{T}, ::LiFTDampingNone) where {T<:AbstractFloat}
    copyto!(factor, normal)
    chol = cholesky!(Hermitian(factor), check=false)
    λ = zero(T)
    if !issuccess(chol)
        λ = regularization_load(normal)
        @inbounds for i in axes(factor, 1)
            factor[i, i] += λ
        end
        chol = cholesky!(Hermitian(factor), check=false)
        if !issuccess(chol)
            λ *= T(10)
            @inbounds for i in axes(factor, 1)
                factor[i, i] += λ
            end
            chol = cholesky!(Hermitian(factor), check=false)
            if !issuccess(chol)
                rhs .= pinv(normal) * rhs
                diag.regularization = λ
                diag.used_fallback = true
                return rhs
            end
        end
    end
    ldiv!(chol, rhs)
    diag.regularization = λ
    return rhs
end

function solve_normal_system!(diag::LiFTDiagnostics{T}, rhs::AbstractVector{T}, factor::AbstractMatrix{T},
    normal::AbstractMatrix{T}, damping::LiFTLevenbergMarquardt) where {T<:AbstractFloat}
    copyto!(factor, normal)
    λ = damping_lambda(damping, normal)
    if λ > zero(T)
        @inbounds for i in axes(factor, 1)
            factor[i, i] += λ
        end
    end
    chol = cholesky!(Hermitian(factor), check=false)
    while !issuccess(chol)
        λ = max(λ * T(damping.growth), regularization_load(normal))
        copyto!(factor, normal)
        @inbounds for i in axes(factor, 1)
            factor[i, i] += λ
        end
        chol = cholesky!(Hermitian(factor), check=false)
        if λ > T(1e12)
            rhs .= pinv(normal) * rhs
            diag.regularization = λ
            diag.used_fallback = true
            return rhs
        end
    end
    ldiv!(chol, rhs)
    diag.regularization = λ
    return rhs
end

function solve_lift_system!(diag::LiFTDiagnostics{T}, residual::AbstractVector{T}, rhs::AbstractVector{T},
    H::AbstractMatrix{T}, normal::AbstractMatrix{T}, factor::AbstractMatrix{T},
    effective_mode::LiFTSolveMode, damping::LiFTDampingMode) where {T<:AbstractFloat}
    diag.used_qr = effective_mode isa LiFTSolveQR
    diag.used_fallback = false
    diag.regularization = zero(T)
    if effective_mode isa LiFTSolveQR
        qr_factor = qr!(H)
        cond_ratio = qr_condition_ratio(qr_factor, size(normal, 1))
        diag.condition_ratio = cond_ratio
        if cond_ratio > damping_condition_limit(T, damping)
            diag.used_qr = false
            diag.used_fallback = true
            solve_normal_system!(diag, rhs, factor, normal, damping)
            return rhs
        end
        try
            ldiv!(qr_factor, residual)
            return residual
        catch err
            if err isa SingularException
                diag.used_qr = false
                diag.used_fallback = true
                solve_normal_system!(diag, rhs, factor, normal, damping)
                return rhs
            end
            rethrow()
        end
    end
    diag.condition_ratio = normal_condition_ratio(normal)
    solve_normal_system!(diag, rhs, factor, normal, damping)
    return rhs
end

@inline delta_component(delta::AbstractVector, j::Int, ::LiFTSolveMode) = delta[j]

@inline function delta_norm(delta::AbstractVector{T}, n_modes::Int, ::LiFTSolveQR) where {T<:AbstractFloat}
    return norm(@view delta[1:n_modes])
end

@inline function delta_norm(delta::AbstractVector{T}, n_modes::Int, ::LiFTSolveNormalEquations) where {T<:AbstractFloat}
    return norm(delta)
end

function regularization_load(normal::AbstractMatrix{T}) where {T<:AbstractFloat}
    diag_sum = zero(T)
    n = min(size(normal, 1), size(normal, 2))
    @inbounds for i in 1:n
        diag_sum += abs(normal[i, i])
    end
    mean_diag = n == 0 ? zero(T) : diag_sum / T(n)
    return max(sqrt(eps(T)) * max(mean_diag, one(T)), eps(T))
end

function damping_lambda(damping::LiFTLevenbergMarquardt, normal::AbstractMatrix{T}) where {T<:AbstractFloat}
    base = regularization_load(normal)
    return max(T(damping.lambda0) * max(base, one(T)), base)
end

@inline damping_condition_limit(::Type{T}, ::LiFTDampingNone) where {T<:AbstractFloat} = inv(sqrt(eps(T)))
@inline damping_condition_limit(::Type{T}, damping::LiFTLevenbergMarquardt) where {T<:AbstractFloat} =
    inv(max(T(damping.condition_rtol), eps(T)))

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
    maxabs = zero(T)
    minabs = typemax(T)
    n = min(size(normal, 1), size(normal, 2))
    @inbounds for i in 1:n
        d = abs(normal[i, i])
        maxabs = max(maxabs, d)
        minabs = min(minabs, d)
    end
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
    @inbounds for idx in eachindex(out)
        denom = vals[idx] + σ2
        out[idx] = inv(max(denom, eps(T)))
    end
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
    @inbounds for idx in eachindex(out)
        out[idx] = inv(max(vals[idx], eps(T)))
    end
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

function psf_from_opd!(lift::LiFT, opd::AbstractMatrix; flux_norm::Real=1.0)
    amp_scale = sqrt(eltype(lift.state.psf_buffer)(photon_flux(lift.src) * flux_norm *
        lift.tel.params.sampling_time * (lift.tel.params.diameter / lift.tel.params.resolution)^2))
    @. lift.state.amp_buffer = ifelse(lift.tel.state.pupil, amp_scale, zero(eltype(lift.state.amp_buffer)))
    oversampling = focal_field_from_opd!(lift.state.focal_buffer, lift, lift.state.amp_buffer, opd)
    field_intensity!(lift.state.psf_buffer, lift.state.focal_buffer, oversampling)
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

function maybe_object_convolve!(lift::LiFT{<:LiFTMode,<:LiFTParams{<:AbstractFloat,<:AbstractMatrix,<:AbstractMatrix}}, mat::AbstractMatrix)
    conv2d_same!(lift.state.conv_buffer, mat, lift.params.object_kernel)
    copyto!(mat, lift.state.conv_buffer)
    return mat
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
        @inbounds for j in 1:n_pad, i in 1:n_pad
            ws.pupil_field[i, j] *= cis(phase_shift * (i + j - 2))
        end
    end

    copyto!(ws.fft_buffer, ws.pupil_field)
    mul!(ws.fft_buffer, ws.fft_plan, ws.fft_buffer)
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

function field_intensity!(dest::AbstractMatrix{T}, field::AbstractMatrix{Complex{T}}, oversampling::Int) where {T<:AbstractFloat}
    if oversampling == 1
        @. dest = abs2(field)
        return dest
    end
    n_out, m_out = size(dest)
    if size(field) != (n_out * oversampling, m_out * oversampling)
        throw(DimensionMismatchError("LiFT field size does not match oversampled PSF dimensions"))
    end
    @inbounds for j in 1:m_out, i in 1:n_out
        acc = zero(T)
        xs = (i - 1) * oversampling + 1
        ys = (j - 1) * oversampling + 1
        for jj in ys:ys+oversampling-1, ii in xs:xs+oversampling-1
            acc += abs2(field[ii, jj])
        end
        dest[i, j] = acc
    end
    return dest
end

function field_derivative!(dest::AbstractMatrix{T}, buf::AbstractMatrix{Complex{T}},
    Pd::AbstractMatrix{Complex{T}}, oversampling::Int, scale::T) where {T<:AbstractFloat}
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
    @inbounds for j in 1:m_out, i in 1:n_out
        acc = zero(T)
        xs = (i - 1) * oversampling + 1
        ys = (j - 1) * oversampling + 1
        for jj in ys:ys+oversampling-1, ii in xs:xs+oversampling-1
            acc += real(im * buf[ii, jj] * Pd[ii, jj])
        end
        dest[i, j] = scale * acc
    end
    return dest
end

function conjugate_field!(dest::AbstractMatrix{Complex{T}}, src::AbstractMatrix{Complex{T}}) where {T<:AbstractFloat}
    @inbounds for idx in eachindex(dest, src)
        dest[idx] = conj(src[idx])
    end
    return dest
end
