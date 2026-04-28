using LinearAlgebra

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
