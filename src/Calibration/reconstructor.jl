using LinearAlgebra

struct ModalReconstructor{T<:AbstractFloat,M<:AbstractMatrix{T}}
    reconstructor::M
    gain::T
end

function ModalReconstructor(imat::InteractionMatrix; gain::Real=1.0)
    T = eltype(imat.matrix)
    recon = pinv(imat.matrix)
    return ModalReconstructor{T, typeof(recon)}(recon, T(gain))
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
