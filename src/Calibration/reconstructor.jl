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
    ModalReconstructor

Linear slopes-to-modes operator built from an interaction matrix inverse.

Mathematically this applies

`c = g * R * s`

where `s` is the slope vector, `R` is the selected inverse of the interaction
matrix, and `g` is the scalar gain stored in the reconstructor.
"""
struct ModalReconstructor{T<:AbstractFloat,M<:AbstractMatrix{T},P<:InversePolicy,V<:AbstractVector{T}}
    reconstructor::M
    gain::T
    policy::P
    singular_values::V
    cond::T
    effective_rank::Int
end

function ModalReconstructor(imat::InteractionMatrix; gain::Real=1.0,
    policy::InversePolicy=default_modal_inverse_policy(eltype(imat.matrix)),
    build_backend::BuildBackend=default_runtime_calibration_build_backend(imat.matrix))
    T = eltype(imat.matrix)
    recon, stats = inverse_operator(build_backend, imat.matrix, policy)
    recon_native = build_backend isa CPUBuildBackend ?
        materialize_build(NativeBuildBackend(), imat.matrix, recon) : recon
    singular_values = build_backend isa CPUBuildBackend ?
        materialize_build(NativeBuildBackend(), similar(imat.matrix, T, 0), stats.singular_values) :
        stats.singular_values
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
}
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
