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
