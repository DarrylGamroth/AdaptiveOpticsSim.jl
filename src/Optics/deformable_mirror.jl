using LinearAlgebra

#
# Deformable mirror model
#
# The DM surface is represented as a linear superposition of influence
# functions:
#
#   opd(x, y) = sum_k modes(x, y, k) * coefs[k]
#
# For the general path this is materialized as a dense matrix-vector product.
# For the common Gaussian/no-rotation/no-anamorphosis case we also build a
# separable approximation
#
#   opd(x, y) = X * C * Y'
#
# where `C` is the actuator coefficient grid. That formulation keeps the same
# mathematical model for this influence family while substantially reducing the
# runtime cost of repeated DM application.
#
@kernel function dm_mode_kernel!(mode, pupil, x_m, y_m, cx, cy, scale, sigma2, n::Int)
    i, j = @index(Global, NTuple)
    if i <= n && j <= n
        x = (i - cx) / scale
        y = (j - cy) / scale
        r2 = (x - x_m)^2 + (y - y_m)^2
        val = exp(-r2 / (2 * sigma2))
        @inbounds mode[(j - 1) * n + i] = ifelse(pupil[i, j], val, zero(eltype(mode)))
    end
end

@kernel function dm_separable_tmp_kernel!(tmp, xbasis, coefs_grid, n_act::Int)
    i, a = @index(Global, NTuple)
    if i <= size(tmp, 1) && a <= size(tmp, 2)
        acc = zero(eltype(tmp))
        @inbounds for b in 1:n_act
            acc = muladd(xbasis[i, b], coefs_grid[b, a], acc)
        end
        tmp[i, a] = acc
    end
end

@kernel function dm_separable_finalize_kernel!(opd, tmp, ybasis_t, pupil, n_act::Int)
    i, j = @index(Global, NTuple)
    if i <= size(opd, 1) && j <= size(opd, 2)
        acc = zero(eltype(opd))
        @inbounds for a in 1:n_act
            acc = muladd(tmp[i, a], ybasis_t[a, j], acc)
        end
        opd[i, j] = ifelse(pupil[i, j], acc, zero(eltype(opd)))
    end
end

abstract type AbstractDMInfluenceModel end

struct GaussianInfluenceWidth{T<:AbstractFloat} <: AbstractDMInfluenceModel
    width::T
end

struct GaussianMechanicalCoupling{T<:AbstractFloat} <: AbstractDMInfluenceModel
    coupling::T
end

struct DenseInfluenceMatrix{T<:AbstractFloat,M<:AbstractMatrix{T}} <: AbstractDMInfluenceModel
    modes::M
end

GaussianInfluenceWidth(width::Real) = GaussianInfluenceWidth(float(width))
GaussianMechanicalCoupling(coupling::Real) = GaussianMechanicalCoupling(float(coupling))
DenseInfluenceMatrix(modes::AbstractMatrix{<:Real}) = DenseInfluenceMatrix(Matrix{float(eltype(modes))}(modes))

struct DeformableMirrorParams{T<:AbstractFloat,I<:AbstractDMInfluenceModel}
    n_act::Int
    influence_model::I
    misregistration::Misregistration{T}
end

mutable struct DeformableMirrorState{
    T<:AbstractFloat,
    A<:AbstractMatrix{T},
    M<:AbstractMatrix{T},
    V<:AbstractVector{T},
    C,
}
    opd::A
    opd_vec::V
    modes::M
    coefs::V
    coefs_grid::C
    separable_x::Union{Nothing,A}
    separable_y_t::Union{Nothing,A}
    separable_tmp::Union{Nothing,A}
end

struct DeformableMirror{P<:DeformableMirrorParams,S<:DeformableMirrorState,B<:AbstractArrayBackend} <: AbstractDeformableMirror
    params::P
    state::S
end

@inline backend(::DeformableMirror{<:Any,<:Any,B}) where {B} = B()

@inline command_storage(dm::DeformableMirror) = dm.state.coefs
@inline command_layout(dm::DeformableMirror) = RuntimeCommandLayout(:dm => length(dm.state.coefs))
@inline influence_model(dm::DeformableMirror) = dm.params.influence_model

@inline dm_normalized_pitch(n_act::Integer) = 2.0 / (n_act - 1)

function mechanical_coupling(n_act::Integer, influence_width::Real)
    n_act > 1 || throw(InvalidConfiguration("mechanical coupling requires n_act > 1"))
    influence_width > 0 || throw(InvalidConfiguration("influence_width must be positive"))
    pitch = dm_normalized_pitch(n_act)
    return exp(-(pitch^2) / (2 * influence_width^2))
end

function influence_width_from_mechanical_coupling(n_act::Integer, coupling::Real)
    n_act > 1 || throw(InvalidConfiguration("mechanical coupling requires n_act > 1"))
    zero(coupling) < coupling < one(coupling) ||
        throw(InvalidConfiguration("mechanical_coupling must lie in (0, 1)"))
    pitch = dm_normalized_pitch(n_act)
    return pitch / sqrt(-2 * log(float(coupling)))
end

@inline influence_width(model::GaussianInfluenceWidth, ::Integer) = model.width
@inline influence_width(model::GaussianMechanicalCoupling, n_act::Integer) =
    typeof(model.coupling)(influence_width_from_mechanical_coupling(n_act, model.coupling))
influence_width(::DenseInfluenceMatrix, ::Integer) =
    throw(UnsupportedAlgorithm("influence_width is only defined for Gaussian DM influence models"))
@inline influence_width(dm::DeformableMirror) = influence_width(influence_model(dm), dm.params.n_act)

@inline mechanical_coupling(model::GaussianInfluenceWidth, n_act::Integer) =
    typeof(model.width)(mechanical_coupling(n_act, model.width))
@inline mechanical_coupling(model::GaussianMechanicalCoupling, ::Integer) = model.coupling
mechanical_coupling(::DenseInfluenceMatrix, ::Integer) =
    throw(UnsupportedAlgorithm("mechanical_coupling is only defined for Gaussian DM influence models"))
@inline mechanical_coupling(dm::DeformableMirror) = mechanical_coupling(influence_model(dm), dm.params.n_act)

@inline _is_identity_misregistration(mis::Misregistration) =
    iszero(mis.shift_x) &&
    iszero(mis.shift_y) &&
    iszero(mis.rotation_rad) &&
    iszero(mis.anamorphosis_angle_rad) &&
    mis.tangential_scaling == one(mis.tangential_scaling) &&
    mis.radial_scaling == one(mis.radial_scaling)

function _convert_dm_influence_model(model::GaussianInfluenceWidth, ::Type{T}, backend, n::Integer, n_act::Integer, misregistration::Misregistration) where {T<:AbstractFloat}
    return GaussianInfluenceWidth{T}(T(model.width))
end

function _convert_dm_influence_model(model::GaussianMechanicalCoupling, ::Type{T}, backend, n::Integer, n_act::Integer, misregistration::Misregistration) where {T<:AbstractFloat}
    return GaussianMechanicalCoupling{T}(T(model.coupling))
end

function _convert_dm_influence_model(model::DenseInfluenceMatrix, ::Type{T}, backend, n::Integer, n_act::Integer, misregistration::Misregistration) where {T<:AbstractFloat}
    _is_identity_misregistration(misregistration) ||
        throw(InvalidConfiguration("DenseInfluenceMatrix already encodes the sampled DM surface and requires identity misregistration"))
    expected = (n * n, n_act * n_act)
    size(model.modes) == expected ||
        throw(DimensionMismatchError("DenseInfluenceMatrix size $(size(model.modes)) does not match expected $(expected)"))
    modes = backend{T}(undef, expected...)
    copyto!(modes, model.modes)
    return DenseInfluenceMatrix{T,typeof(modes)}(modes)
end

function _resolve_dm_influence_model(; n::Integer, n_act::Integer,
    influence_width::Union{Nothing,Real}, mechanical_coupling::Union{Nothing,Real},
    influence_model::Union{Nothing,AbstractDMInfluenceModel}, T::Type{<:AbstractFloat},
    backend, misregistration::Misregistration)

    if influence_model !== nothing
        influence_width === nothing ||
            throw(InvalidConfiguration("specify either influence_model or influence_width, not both"))
        mechanical_coupling === nothing ||
            throw(InvalidConfiguration("specify either influence_model or mechanical_coupling, not both"))
        return _convert_dm_influence_model(influence_model, T, backend, n, n_act, misregistration)
    end
    if mechanical_coupling !== nothing
        influence_width === nothing ||
            throw(InvalidConfiguration("specify either influence_width or mechanical_coupling, not both"))
        return GaussianMechanicalCoupling{T}(T(mechanical_coupling))
    end
    resolved_width = influence_width === nothing ? T(0.2) : T(influence_width)
    return GaussianInfluenceWidth{T}(resolved_width)
end

"""
    DeformableMirror(tel; ...)

Build a DM whose surface is a weighted sum of Gaussian influence functions on
the telescope pupil grid.

The stored dense operator is the flattened matrix `modes`, mapping actuator
coefficients to OPD samples. When the misregistration keeps the Gaussian basis
separable, the constructor also prepares a faster `X * C * Y'` runtime path.
"""
function DeformableMirror(tel::Telescope; n_act::Int,
    influence_width::Union{Nothing,Real}=nothing, mechanical_coupling::Union{Nothing,Real}=nothing,
    influence_model::Union{Nothing,AbstractDMInfluenceModel}=nothing,
    T::Type{<:AbstractFloat}=Float64, misregistration::Misregistration=Misregistration(T=T), backend::AbstractArrayBackend=backend(tel))
    selector = require_same_backend(tel, _resolve_backend_selector(backend))
    backend = _resolve_array_backend(selector)
    n = tel.params.resolution
    resolved_model = _resolve_dm_influence_model(
        n=n,
        n_act=n_act,
        influence_width=influence_width,
        mechanical_coupling=mechanical_coupling,
        influence_model=influence_model,
        T=T,
        backend=backend,
        misregistration=misregistration,
    )
    params = DeformableMirrorParams{T,typeof(resolved_model)}(n_act, resolved_model, misregistration)
    opd = backend{T}(undef, n, n)
    fill!(opd, zero(T))
    opd_vec = reshape(opd, :)
    modes = backend{T}(undef, n * n, n_act * n_act)
    coefs = backend{T}(undef, n_act * n_act)
    coefs_grid = reshape(coefs, n_act, n_act)
    fill!(coefs, zero(T))
    state = DeformableMirrorState{T, typeof(opd), typeof(modes), typeof(opd_vec), typeof(coefs_grid)}(
        opd, opd_vec, modes, coefs, coefs_grid, nothing, nothing, nothing)
    dm = DeformableMirror{typeof(params), typeof(state), typeof(selector)}(params, state)
    build_influence_functions!(dm, tel)
    build_separable_influence!(dm, tel)
    return dm
end

@inline function supports_separable_influence(mis::Misregistration)
    return iszero(mis.rotation_rad) &&
           iszero(mis.anamorphosis_angle_rad) &&
           mis.tangential_scaling == one(mis.tangential_scaling) &&
           mis.radial_scaling == one(mis.radial_scaling)
end

@inline supports_separable_influence(::AbstractDMInfluenceModel, ::Misregistration) = false
@inline supports_separable_influence(::GaussianInfluenceWidth, mis::Misregistration) = supports_separable_influence(mis)
@inline supports_separable_influence(::GaussianMechanicalCoupling, mis::Misregistration) = supports_separable_influence(mis)

@inline function _clear_separable_influence!(dm::DeformableMirror)
    dm.state.separable_x = nothing
    dm.state.separable_y_t = nothing
    dm.state.separable_tmp = nothing
    return dm
end

"""
    build_separable_influence!(dm, tel)

Build the factored Gaussian influence basis used by the fast separable DM
application path.

This optimization is only valid when the DM misregistration preserves the
separable tensor-product structure of the Gaussian influence model.
"""
function build_separable_influence!(dm::DeformableMirror, tel::Telescope)
    return build_separable_influence!(dm, tel, influence_model(dm))
end

function build_separable_influence!(dm::DeformableMirror, tel::Telescope, model::AbstractDMInfluenceModel)
    supports_separable_influence(model, dm.params.misregistration) || return _clear_separable_influence!(dm)
    n = tel.params.resolution
    n_act = dm.params.n_act
    T = eltype(dm.state.opd)
    sigma2 = T(influence_width(model, n_act))^2
    xs = range(T(-1), T(1); length=n_act)
    ys = range(T(-1), T(1); length=n_act)
    cx = T((n + 1) / 2)
    cy = T((n + 1) / 2)
    scale = T(n / 2)
    xbasis_host = Matrix{T}(undef, n, n_act)
    ybasis_host = Matrix{T}(undef, n_act, n)
    @inbounds for (ai, x0) in enumerate(xs)
        x_m, _ = apply_misregistration(dm.params.misregistration, x0, zero(T))
        for i in 1:n
            x = (T(i) - cx) / scale
            xbasis_host[i, ai] = exp(-((x - T(x_m))^2) / (2 * sigma2))
        end
    end
    @inbounds for (aj, y0) in enumerate(ys)
        _, y_m = apply_misregistration(dm.params.misregistration, zero(T), y0)
        for j in 1:n
            y = (T(j) - cy) / scale
            ybasis_host[aj, j] = exp(-((y - T(y_m))^2) / (2 * sigma2))
        end
    end
    xbasis = similar(dm.state.opd, n, n_act)
    ybasis_t = similar(dm.state.opd, n_act, n)
    tmp = similar(dm.state.opd, n, n_act)
    copyto!(xbasis, xbasis_host)
    copyto!(ybasis_t, ybasis_host)
    fill!(tmp, zero(T))
    dm.state.separable_x = xbasis
    dm.state.separable_y_t = ybasis_t
    dm.state.separable_tmp = tmp
    return dm
end

function build_influence_functions!(dm::DeformableMirror, tel::Telescope)
    return build_influence_functions!(execution_style(dm.state.modes), dm, tel, influence_model(dm))
end

"""
    build_influence_functions!(dm, tel)

Materialize the dense actuator-to-OPD operator.

Each column corresponds to one actuator's influence function sampled on the
telescope pupil grid. This is the reference linear model used for calibration
and for the dense DM application fallback.
"""
function build_influence_functions!(::ScalarCPUStyle, dm::DeformableMirror, tel::Telescope, model::AbstractDMInfluenceModel)
    Base.require_one_based_indexing(tel.state.pupil, dm.state.modes)
    n = tel.params.resolution
    n_act = dm.params.n_act
    sigma = influence_width(model, n_act)
    xs = range(-1.0, 1.0; length=n_act)
    ys = range(-1.0, 1.0; length=n_act)

    cx = (n + 1) / 2
    cy = (n + 1) / 2
    scale = n / 2

    idx = 1
    @inbounds for x0 in xs, y0 in ys
        x_m, y_m = apply_misregistration(dm.params.misregistration, x0, y0)
        for i in 1:n, j in 1:n
            x = (i - cx) / scale
            y = (j - cy) / scale
            r2 = (x - x_m)^2 + (y - y_m)^2
            dm.state.modes[(j - 1) * n + i, idx] = exp(-r2 / (2 * sigma^2)) * tel.state.pupil[i, j]
        end
        idx += 1
    end
    return dm
end

function build_influence_functions!(style::AcceleratorStyle, dm::DeformableMirror, tel::Telescope, model::AbstractDMInfluenceModel)
    Base.require_one_based_indexing(tel.state.pupil, dm.state.modes)
    n = tel.params.resolution
    n_act = dm.params.n_act
    T = eltype(dm.state.modes)
    sigma2 = T(influence_width(model, n_act))^2
    xs = range(T(-1), T(1); length=n_act)
    ys = range(T(-1), T(1); length=n_act)

    cx = T((n + 1) / 2)
    cy = T((n + 1) / 2)
    scale = T(n / 2)

    idx = 1
    for x0 in xs, y0 in ys
        x_m, y_m = apply_misregistration(dm.params.misregistration, x0, y0)
        launch_kernel!(style, dm_mode_kernel!, @view(dm.state.modes[:, idx]), tel.state.pupil,
            T(x_m), T(y_m), cx, cy, scale, sigma2, n; ndrange=(n, n))
        idx += 1
    end
    synchronize_backend!(style)
    return dm
end

function build_influence_functions!(::ScalarCPUStyle, dm::DeformableMirror, tel::Telescope, model::DenseInfluenceMatrix)
    Base.require_one_based_indexing(dm.state.modes)
    size(model.modes) == size(dm.state.modes) ||
        throw(DimensionMismatchError("DenseInfluenceMatrix size $(size(model.modes)) does not match DM state size $(size(dm.state.modes))"))
    copyto!(dm.state.modes, model.modes)
    return dm
end

function build_influence_functions!(::AcceleratorStyle, dm::DeformableMirror, tel::Telescope, model::DenseInfluenceMatrix)
    size(model.modes) == size(dm.state.modes) ||
        throw(DimensionMismatchError("DenseInfluenceMatrix size $(size(model.modes)) does not match DM state size $(size(dm.state.modes))"))
    copyto!(dm.state.modes, model.modes)
    return dm
end

function apply!(dm::DeformableMirror, tel::Telescope, ::DMAdditive)
    apply_opd!(dm, tel)
    tel.state.opd .+= dm.state.opd
    return tel
end

function apply!(dm::DeformableMirror, tel::Telescope, ::DMReplace)
    apply_opd!(dm, tel)
    tel.state.opd .= dm.state.opd
    return tel
end

@inline function _apply_opd_separable!(::ScalarCPUStyle, dm::DeformableMirror, tel::Telescope)
    xbasis = dm.state.separable_x::typeof(dm.state.opd)
    ybasis_t = dm.state.separable_y_t::typeof(dm.state.opd)
    tmp = dm.state.separable_tmp::typeof(dm.state.opd)
    mul!(tmp, xbasis, dm.state.coefs_grid)
    mul!(dm.state.opd, tmp, ybasis_t)
    dm.state.opd .*= tel.state.pupil
    return dm.state.opd
end

@inline function _apply_opd_separable!(style::AcceleratorStyle, dm::DeformableMirror, tel::Telescope)
    xbasis = dm.state.separable_x::typeof(dm.state.opd)
    ybasis_t = dm.state.separable_y_t::typeof(dm.state.opd)
    tmp = dm.state.separable_tmp::typeof(dm.state.opd)
    n_act = dm.params.n_act
    launch_kernel_async!(style, dm_separable_tmp_kernel!,
        tmp, xbasis, dm.state.coefs_grid, n_act; ndrange=size(tmp))
    launch_kernel!(style, dm_separable_finalize_kernel!,
        dm.state.opd, tmp, ybasis_t, tel.state.pupil, n_act; ndrange=size(dm.state.opd))
    return dm.state.opd
end

@inline function apply_dense!(dm::DeformableMirror, tel::Telescope, ::DMAdditive)
    mul!(dm.state.opd_vec, dm.state.modes, dm.state.coefs)
    tel.state.opd .+= dm.state.opd
    return tel
end

@inline function apply_dense!(dm::DeformableMirror, tel::Telescope, ::DMReplace)
    mul!(dm.state.opd_vec, dm.state.modes, dm.state.coefs)
    tel.state.opd .= dm.state.opd
    return tel
end

@inline function apply_opd!(dm::DeformableMirror, tel::Telescope)
    if !isnothing(dm.state.separable_x)
        return apply_opd_separable!(dm, tel)
    end
    mul!(dm.state.opd_vec, dm.state.modes, dm.state.coefs)
    return dm.state.opd
end

"""
    apply_opd!(dm, tel)

Assemble the DM OPD from the current actuator coefficients without mutating the
telescope OPD.

This chooses the separable `X * C * Y'` path when available and otherwise falls
back to the dense matrix-vector application `modes * coefs`.
"""
@inline function apply_opd_separable!(dm::DeformableMirror, tel::Telescope)
    return _apply_opd_separable!(execution_style(dm.state.opd), dm, tel)
end
