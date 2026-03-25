using LinearAlgebra

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

function reconstruct!(out::AbstractVector, recon::ModalReconstructor, slopes::AbstractVector)
    mul!(out, recon.reconstructor, slopes)
    out .*= recon.gain
    return out
end

function reconstruct(recon::ModalReconstructor, slopes::AbstractVector)
    out = similar(slopes, size(recon.reconstructor, 1))
    return reconstruct!(out, recon, slopes)
end

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
