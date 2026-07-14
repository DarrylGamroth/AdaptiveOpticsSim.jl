#
# Reconstructor operators
#
# The core control law is always a linear map from WFS slopes to DM commands.
# This file exposes two closely related forms:
#
# - `ModalReconstructor`: c = g * R * s
# - `MappedReconstructor`: u = B * (g * R * s)
#
# where:
# - `s` is the measured slope vector
# - `R` is the inverse/pseudoinverse of the interaction matrix
# - `g` is a scalar loop gain
# - `B` is an explicit command-basis map, typically M2C
#
# Keeping the two-stage mapped form explicit matters for both numerical
# clarity and performance. It avoids baking large dense `B * R` products into
# one operator when the runtime wants separate modal and command-basis stages.
#
"""
    AbstractReconstructorOperator

Abstract slopes-to-command operator family used by closed-loop runtime control.

The maintained contract is:

- `reconstruct!(out, recon, slopes)`
- `reconstruct(recon, slopes)`

Concrete implementations may internally reconstruct in modal space, command
space, or a mapped two-stage basis, but they should present the same external
slopes-to-command interface.
"""
abstract type AbstractReconstructorOperator end

"""
    runtime_reconstructor_storage(reconstructor)

Return the tuple of array-backed state used by `reconstruct!` at runtime, or
`nothing` when the reconstructor has not declared its runtime storage.

Device-resident execution plans use this trait to reject control operators
that would otherwise introduce an implicit host/device boundary. Custom
device-resident reconstructors should extend this function with their control
matrix and any hot-path workspaces. Backend-agnostic reconstructors return an
empty tuple.
"""
runtime_reconstructor_storage(::Any) = nothing
runtime_reconstructor_ownership_roots(::Any) = ()

function inverse_policy(::AbstractReconstructorOperator)
    throw(InvalidConfiguration("inverse_policy is not defined for this reconstructor family"))
end

function singular_values(::AbstractReconstructorOperator)
    throw(InvalidConfiguration("singular_values is not defined for this reconstructor family"))
end

function condition_number(::AbstractReconstructorOperator)
    throw(InvalidConfiguration("condition_number is not defined for this reconstructor family"))
end

function effective_rank(::AbstractReconstructorOperator)
    throw(InvalidConfiguration("effective_rank is not defined for this reconstructor family"))
end

"""
    NullReconstructor()

Explicit no-op control operator for external-command or HIL-facing runtimes.

`NullReconstructor` exists so a runtime boundary can be assembled without
constructing a fake slopes-to-command map. It is valid with `sense!` plus
`set_command!`, but `step!` / `reconstruct!` will reject it because there is no
internal reconstruction stage to run.
"""
struct NullReconstructor <: AbstractReconstructorOperator end

@inline runtime_reconstructor_storage(::NullReconstructor) = ()

function reconstruct!(out::AbstractVector, ::NullReconstructor, slopes::AbstractVector)
    throw(InvalidConfiguration("NullReconstructor does not define an internal slopes-to-command update; use set_command! with sense! for external-control runtimes"))
end

function reconstruct(::NullReconstructor, slopes::AbstractVector)
    throw(InvalidConfiguration("NullReconstructor does not define an internal slopes-to-command update; use set_command! with sense! for external-control runtimes"))
end

"""
    ModalReconstructor

Linear slopes-to-modes operator built from an interaction matrix inverse.

Mathematically this applies

`c = g * R * s`

where `s` is the slope vector, `R` is the selected inverse of the interaction
matrix, and `g` is the scalar gain stored in the reconstructor.
"""
struct ModalReconstructor{T<:AbstractFloat,M<:AbstractMatrix{T},P<:InversePolicy,V<:AbstractVector{T}} <: AbstractReconstructorOperator
    reconstructor::M
    gain::T
    policy::P
    singular_values::V
    cond::T
    effective_rank::Int
end

@inline runtime_reconstructor_storage(recon::ModalReconstructor) =
    (recon.reconstructor,)

@inline inverse_policy(recon::ModalReconstructor) = recon.policy
@inline singular_values(recon::ModalReconstructor) = recon.singular_values
@inline condition_number(recon::ModalReconstructor) = recon.cond
@inline effective_rank(recon::ModalReconstructor) = recon.effective_rank

function ModalReconstructor(imat::InteractionMatrix; gain::Real=1.0,
    policy::InversePolicy=default_modal_inverse_policy(eltype(imat.matrix)),
    build_backend::BuildBackend=default_runtime_calibration_build_backend(imat.matrix))
    T = eltype(imat.matrix)
    recon, stats = inverse_operator(build_backend, imat.matrix, policy)
    recon_native = materialize_runtime_build_result(build_backend, imat.matrix, recon)
    singular_values = materialize_runtime_build_result(build_backend, similar(imat.matrix, T, 0), stats.singular_values)
    return ModalReconstructor{T, typeof(recon_native), typeof(policy), typeof(singular_values)}(
        recon_native,
        T(gain),
        policy,
        singular_values,
        stats.cond,
        stats.effective_rank,
    )
end

"""
    reconstruct!(out, recon::ModalReconstructor, slopes)

Apply the modal control law `out = gain * R * slopes` in-place.

This is the hot-path runtime form used by closed-loop control and keeps the
output buffer under caller ownership.
"""
function reconstruct!(out::AbstractVector, recon::ModalReconstructor, slopes::AbstractVector)
    mul!(out, recon.reconstructor, slopes)
    out .*= recon.gain
    return out
end

function reconstruct(recon::ModalReconstructor, slopes::AbstractVector)
    out = similar(slopes, size(recon.reconstructor, 1))
    return reconstruct!(out, recon, slopes)
end

"""
    FactorizedReconstructor(imat; gain=1, policy=..., max_rank=nothing,
                            build_backend=...)

SVD-factor slopes-to-command operator that applies

`c = gain * Vᵣ * Diagonal(wᵣ) * Uᵣ' * s`

without materializing the dense inverse matrix. `max_rank` provides an
explicit runtime storage/compute bound; omitted trailing zero-weight factors
are always discarded. The factor matrices and workspace are materialized on
the interaction matrix backend after any requested CPU build.
"""
struct FactorizedReconstructor{
    T<:AbstractFloat,
    U<:AbstractMatrix{T},
    V<:AbstractMatrix{T},
    W<:AbstractVector{T},
    WS<:AbstractVector{T},
    P<:InversePolicy,
    SV<:AbstractVector{T},
} <: AbstractReconstructorOperator
    left_modes::U
    command_modes::V
    inverse_weights::W
    workspace::WS
    gain::T
    policy::P
    singular_values::SV
    cond::T
    effective_rank::Int
    factor_rank::Int
    n_trunc::Int
end

@inline runtime_reconstructor_storage(recon::FactorizedReconstructor) = (
    recon.left_modes,
    recon.command_modes,
    recon.inverse_weights,
    recon.workspace,
)
@inline runtime_reconstructor_ownership_roots(recon::FactorizedReconstructor) =
    (recon.workspace,)
@inline inverse_policy(recon::FactorizedReconstructor) = recon.policy
@inline singular_values(recon::FactorizedReconstructor) = recon.singular_values
@inline condition_number(recon::FactorizedReconstructor) = recon.cond
@inline effective_rank(recon::FactorizedReconstructor) = recon.effective_rank
@inline factorized_rank(recon::FactorizedReconstructor) = recon.factor_rank
@inline truncation_count(recon::FactorizedReconstructor) = recon.n_trunc

function _factorized_reconstructor_rank(inverse_weights::AbstractVector,
    max_rank::Union{Integer,Nothing})
    available = something(findlast(!iszero, inverse_weights), 0)
    isnothing(max_rank) && return available
    max_rank >= 0 || throw(InvalidConfiguration("max_rank must be >= 0"))
    return min(Int(max_rank), available)
end

function _compact_factor_matrix(build_backend::BuildBackend,
    ref::AbstractMatrix, factor::AbstractMatrix, rank::Int)
    compact = copy(@view(factor[:, 1:rank]))
    return materialize_runtime_build_result(build_backend, ref, compact)
end

function FactorizedReconstructor(imat::InteractionMatrix{T};
    gain::Real=1.0,
    policy::InversePolicy=default_modal_inverse_policy(T),
    max_rank::Union{Integer,Nothing}=nothing,
    build_backend::BuildBackend=
        default_runtime_calibration_build_backend(imat.matrix),
) where {T<:AbstractFloat}
    F, inverse_weights_host, stats = inverse_factorization(
        build_backend,
        imat.matrix,
        policy,
    )
    rank = _factorized_reconstructor_rank(inverse_weights_host, max_rank)
    left_modes = _compact_factor_matrix(build_backend, imat.matrix, F.U,
        rank)
    command_modes = _compact_factor_matrix(build_backend, imat.matrix, F.V,
        rank)
    weight_ref = similar(imat.matrix, T, 0)
    inverse_weights = materialize_runtime_build_result(
        build_backend,
        weight_ref,
        copy(@view(inverse_weights_host[1:rank])),
    )
    singular_values = materialize_runtime_build_result(
        build_backend,
        weight_ref,
        stats.singular_values,
    )
    workspace = similar(inverse_weights)
    fill!(workspace, zero(T))
    effective = min(stats.effective_rank, rank)
    cond = if effective == 0
        T(Inf)
    elseif effective < stats.effective_rank
        stats.singular_values[begin] / stats.singular_values[effective]
    else
        stats.cond
    end
    return FactorizedReconstructor{
        T,
        typeof(left_modes),
        typeof(command_modes),
        typeof(inverse_weights),
        typeof(workspace),
        typeof(policy),
        typeof(singular_values),
    }(
        left_modes,
        command_modes,
        inverse_weights,
        workspace,
        T(gain),
        policy,
        singular_values,
        cond,
        effective,
        rank,
        length(stats.singular_values) - rank,
    )
end

function reconstruct!(out::AbstractVector, recon::FactorizedReconstructor,
    slopes::AbstractVector)
    mul!(recon.workspace, adjoint(recon.left_modes), slopes)
    recon.workspace .= recon.gain .* recon.inverse_weights .* recon.workspace
    mul!(out, recon.command_modes, recon.workspace)
    return out
end

function reconstruct(recon::FactorizedReconstructor, slopes::AbstractVector)
    out = similar(slopes, size(recon.command_modes, 1))
    return reconstruct!(out, recon, slopes)
end

"""
    MappedReconstructor

Two-stage slopes-to-commands operator with an explicit command-basis map.

This represents

`u = B * (g * R * s)`

where `R` is the modal inverse operator and `B` maps modal coefficients to the
final command basis. This is the natural runtime form for systems that
reconstruct in one basis but actuate in another.
"""
struct MappedReconstructor{
    T<:AbstractFloat,
    M<:AbstractMatrix{T},
    B<:AbstractMatrix{T},
    P<:InversePolicy,
    V<:AbstractVector{T},
    W<:AbstractVector{T},
} <: AbstractReconstructorOperator
    reconstructor::M
    command_basis::B
    modal_workspace::W
    gain::T
    policy::P
    singular_values::V
    cond::T
    effective_rank::Int
    n_control_modes::Int
end


@inline runtime_reconstructor_storage(recon::MappedReconstructor) =
    (recon.reconstructor, recon.command_basis, recon.modal_workspace)
@inline runtime_reconstructor_ownership_roots(recon::MappedReconstructor) =
    (recon.modal_workspace,)

@inline inverse_policy(recon::MappedReconstructor) = recon.policy
@inline singular_values(recon::MappedReconstructor) = recon.singular_values
@inline condition_number(recon::MappedReconstructor) = recon.cond
@inline effective_rank(recon::MappedReconstructor) = recon.effective_rank

function MappedReconstructor(command_basis::AbstractMatrix{T}, imat::InteractionMatrix{T};
    gain::Real=1.0,
    policy::InversePolicy=default_modal_inverse_policy(T),
    inverse_build_backend::BuildBackend=default_runtime_calibration_build_backend(imat.matrix),
    materialize_backend::BuildBackend=inverse_build_backend,
    ref::AbstractMatrix{T}=imat.matrix) where {T<:AbstractFloat}
    recon, stats = inverse_operator(inverse_build_backend, imat.matrix, policy)
    recon_native = materialize_build(materialize_backend, ref, recon)
    command_basis_native = materialize_build(materialize_backend, ref, command_basis)
    singular_values = materialize_build(materialize_backend, similar(ref, T, 0), stats.singular_values)
    modal_workspace = similar(ref, T, size(command_basis, 2))
    fill!(modal_workspace, zero(T))
    return MappedReconstructor{
        T,
        typeof(recon_native),
        typeof(command_basis_native),
        typeof(policy),
        typeof(singular_values),
        typeof(modal_workspace),
    }(
        recon_native,
        command_basis_native,
        modal_workspace,
        T(gain),
        policy,
        singular_values,
        stats.cond,
        stats.effective_rank,
        size(command_basis, 2),
    )
end

"""
    reconstruct!(out, recon::MappedReconstructor, slopes)

Apply the two-stage mapped control law

`out = command_basis * (gain * reconstructor * slopes)`

using the reconstructor's preallocated modal workspace.
"""
function reconstruct!(out::AbstractVector, recon::MappedReconstructor, slopes::AbstractVector)
    mul!(recon.modal_workspace, recon.reconstructor, slopes)
    recon.modal_workspace .*= recon.gain
    mul!(out, recon.command_basis, recon.modal_workspace)
    return out
end

function reconstruct(recon::MappedReconstructor, slopes::AbstractVector)
    out = similar(recon.command_basis, size(recon.command_basis, 1))
    return reconstruct!(out, recon, slopes)
end

@inline reconstructor_output_length(::NullReconstructor) = nothing
@inline reconstructor_output_length(recon::ModalReconstructor) =
    size(recon.reconstructor, 1)
@inline reconstructor_output_length(recon::FactorizedReconstructor) =
    size(recon.command_modes, 1)
@inline reconstructor_output_length(recon::MappedReconstructor) =
    size(recon.command_basis, 1)

@inline function _controlled_reconstructor_storage(reconstructor, controller,
    workspace)
    reconstructor_storage = runtime_reconstructor_storage(reconstructor)
    controller_storage = runtime_controller_storage(controller)
    isnothing(reconstructor_storage) && return nothing
    isnothing(controller_storage) && return nothing
    return (reconstructor_storage..., controller_storage..., workspace)
end

"""
    ControlledReconstructor(reconstructor, controller; dt)

Compose a slopes-to-command reconstructor with a stateful command controller.
The intermediate command vector is preallocated, so repeated runtime
`reconstruct!` calls remain allocation-free. Device-resident validation sees
the storage owned by both operators as one residency contract.
"""
struct ControlledReconstructor{
    R<:AbstractReconstructorOperator,
    C<:AbstractController,
    W<:AbstractVector,
    T<:AbstractFloat,
} <: AbstractReconstructorOperator
    reconstructor::R
    controller::C
    workspace::W
    dt::T
end

function ControlledReconstructor(reconstructor::AbstractReconstructorOperator,
    controller::AbstractController; dt::Real)
    output_length = reconstructor_output_length(reconstructor)
    isnothing(output_length) &&
        throw(InvalidConfiguration("ControlledReconstructor requires a reconstructor with a declared output length"))
    output = controller_output(controller)
    length(output) == output_length ||
        throw(DimensionMismatchError("controller state length must match reconstructor output length"))
    T = eltype(output)
    dt_value = T(dt)
    dt_value > zero(T) ||
        throw(InvalidConfiguration("controlled reconstructor dt must be greater than zero at controller precision"))
    storage = runtime_reconstructor_storage(reconstructor)
    controller_storage = runtime_controller_storage(controller)
    if !isnothing(storage) && !isnothing(controller_storage)
        require_same_backend(storage..., controller_storage...)
    end
    workspace = similar(output)
    fill!(workspace, zero(eltype(workspace)))
    return ControlledReconstructor{
        typeof(reconstructor),
        typeof(controller),
        typeof(workspace),
        T,
    }(
        reconstructor,
        controller,
        workspace,
        dt_value,
    )
end

@inline runtime_reconstructor_storage(recon::ControlledReconstructor) =
    _controlled_reconstructor_storage(recon.reconstructor, recon.controller,
        recon.workspace)
@inline runtime_reconstructor_ownership_roots(recon::ControlledReconstructor) =
    (
        runtime_reconstructor_ownership_roots(recon.reconstructor)...,
        runtime_controller_ownership_roots(recon.controller)...,
        recon.workspace,
    )
@inline inverse_policy(recon::ControlledReconstructor) =
    inverse_policy(recon.reconstructor)
@inline singular_values(recon::ControlledReconstructor) =
    singular_values(recon.reconstructor)
@inline condition_number(recon::ControlledReconstructor) =
    condition_number(recon.reconstructor)
@inline effective_rank(recon::ControlledReconstructor) =
    effective_rank(recon.reconstructor)
@inline controller_output(recon::ControlledReconstructor) =
    controller_output(recon.controller)
@inline supports_controller_reset(recon::ControlledReconstructor) =
    supports_controller_reset(recon.controller)

function reset_controller!(recon::ControlledReconstructor)
    reset_controller!(recon.controller)
    fill!(recon.workspace, zero(eltype(recon.workspace)))
    return recon
end

function reconstruct!(out::AbstractVector, recon::ControlledReconstructor,
    slopes::AbstractVector)
    reconstruct!(recon.workspace, recon.reconstructor, slopes)
    controlled = update!(recon.controller, recon.workspace, recon.dt)
    copyto!(out, controlled)
    return out
end

function reconstruct(recon::ControlledReconstructor, slopes::AbstractVector)
    out = similar(controller_output(recon))
    return reconstruct!(out, recon, slopes)
end
