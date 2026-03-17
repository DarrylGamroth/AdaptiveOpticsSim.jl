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
    policy::InversePolicy=TSVDInverse())
    T = eltype(imat.matrix)
    recon, stats = inverse_operator(imat.matrix, policy)
    return ModalReconstructor{T, typeof(recon), typeof(policy), typeof(stats.singular_values)}(
        recon,
        T(gain),
        policy,
        stats.singular_values,
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
