abstract type InversePolicy end
abstract type BuildBackend end

struct NativeBuildBackend <: BuildBackend end
struct CPUBuildBackend <: BuildBackend end
struct GPUArrayBuildBackend{B} <: BuildBackend end

struct ExactPseudoInverse <: InversePolicy end

struct TSVDInverse{T<:AbstractFloat} <: InversePolicy
    rtol::T
    atol::T
    n_trunc::Int
end

struct TikhonovInverse{T<:AbstractFloat} <: InversePolicy
    lambda::T
    rtol::T
    atol::T
end

function TSVDInverse(; rtol::Real=eps(Float64), atol::Real=0.0, n_trunc::Integer=0)
    T = promote_type(typeof(float(rtol)), typeof(float(atol)))
    return TSVDInverse(T(rtol), T(atol), Int(n_trunc))
end

function TikhonovInverse(lambda::Real; rtol::Real=eps(Float64), atol::Real=0.0)
    T = promote_type(typeof(float(lambda)), typeof(float(rtol)), typeof(float(atol)))
    return TikhonovInverse(T(lambda), T(rtol), T(atol))
end

default_modal_inverse_policy(::Type{T}) where {T<:AbstractFloat} = TSVDInverse(rtol=sqrt(eps(T)))
default_calibration_inverse_policy(::Type{T}) where {T<:AbstractFloat} = TSVDInverse(rtol=sqrt(eps(T)))
default_projector_inverse_policy(::Type{T}) where {T<:AbstractFloat} = TSVDInverse(rtol=sqrt(eps(T)))
default_build_backend(::AbstractArray) = NativeBuildBackend()
function default_runtime_calibration_build_backend(A::AbstractArray)
    gpu_backend_name(typeof(A)) === nothing && return NativeBuildBackend()
    return CPUBuildBackend()
end
GPUArrayBuildBackend(::Type{B}) where {B} = GPUArrayBuildBackend{B}()

struct InverseStats{T<:AbstractFloat,V<:AbstractVector{T}}
    singular_values::V
    cond::T
    effective_rank::Int
    n_trunc::Int
end

function _inverse_cutoff(s::AbstractVector{T}, rtol::T, atol::T) where {T<:AbstractFloat}
    isempty(s) && return zero(T)
    return max(atol, rtol * s[begin])
end

prepare_build_matrix(::NativeBuildBackend, A::AbstractMatrix) = A
prepare_build_matrix(::CPUBuildBackend, A::AbstractMatrix) = Matrix(A)
prepare_build_matrix(backend::GPUArrayBuildBackend, A::AbstractMatrix) = materialize_build(backend, A)

materialize_build(::NativeBuildBackend, A::AbstractMatrix) = A
materialize_build(::CPUBuildBackend, A::AbstractMatrix) = Matrix(A)
materialize_build(::NativeBuildBackend, A::AbstractVector) = A
materialize_build(::CPUBuildBackend, A::AbstractVector) = Vector(A)

function _backend_array(::Type{B}, ::Type{T}, dims::Vararg{Int,N}) where {B,T,N}
    return backend_fill(B, zero(T), dims...)
end

function _backend_array(::Type{B}, ::Type{Bool}, dims::Vararg{Int,N}) where {B,N}
    return backend_fill(B, false, dims...)
end

@inline function _copy_build_data!(out, A::AbstractArray)
    copyto!(out, A)
    return out
end

@inline function _copy_build_data!(out, A::SubArray{T,N,<:Array}) where {T,N}
    copyto!(out, Array{T,N}(A))
    return out
end

@inline function _copy_build_data!(out, A::Transpose{T,<:Array{T,2}}) where {T}
    copyto!(out, Matrix{T}(A))
    return out
end

@inline function _copy_build_data!(out, A::Adjoint{T,<:Array{T,2}}) where {T}
    copyto!(out, Matrix{T}(A))
    return out
end

function materialize_build(::GPUArrayBuildBackend{B}, A::AbstractMatrix{T}) where {B,T}
    out = _backend_array(B, T, size(A)...)
    _copy_build_data!(out, A)
    return out
end

function materialize_build(::GPUArrayBuildBackend{B}, A::SparseMatrixCSC{T}) where {B,T}
    out = _backend_array(B, T, size(A)...)
    copyto!(out, Matrix(A))
    return out
end

function materialize_build(::GPUArrayBuildBackend{B}, A::BitMatrix) where {B}
    out = _backend_array(B, Bool, size(A)...)
    copyto!(out, Matrix{Bool}(A))
    return out
end

function materialize_build(::GPUArrayBuildBackend{B}, A::AbstractVector{T}) where {B,T}
    out = _backend_array(B, T, length(A))
    _copy_build_data!(out, A)
    return out
end

function materialize_build(::GPUArrayBuildBackend{B}, A::BitVector) where {B}
    out = _backend_array(B, Bool, length(A))
    copyto!(out, Vector{Bool}(A))
    return out
end

function materialize_build(::NativeBuildBackend, ref::AbstractMatrix, data::AbstractMatrix)
    out = similar(ref, eltype(data), size(data)...)
    copyto!(out, data)
    return out
end

materialize_build(::CPUBuildBackend, ::AbstractMatrix, data::AbstractMatrix) = Matrix(data)

function materialize_build(::GPUArrayBuildBackend{B}, ref::AbstractMatrix, data::AbstractMatrix) where {B}
    out = _backend_array(B, eltype(data), size(data)...)
    _copy_build_data!(out, data)
    return out
end

function materialize_build(::GPUArrayBuildBackend{B}, ref::AbstractMatrix, data::SparseMatrixCSC{T}) where {B,T}
    out = _backend_array(B, T, size(data)...)
    copyto!(out, Matrix(data))
    return out
end

function materialize_build(::GPUArrayBuildBackend{B}, ref::AbstractMatrix, data::BitMatrix) where {B}
    out = _backend_array(B, Bool, size(data)...)
    copyto!(out, Matrix{Bool}(data))
    return out
end

function materialize_build(::NativeBuildBackend, ref::AbstractVector{T}, data::AbstractVector{T}) where {T}
    out = similar(ref, T, length(data))
    copyto!(out, data)
    return out
end

materialize_build(::CPUBuildBackend, ::AbstractVector{T}, data::AbstractVector{T}) where {T} = Vector{T}(data)

function materialize_build(::GPUArrayBuildBackend{B}, ref::AbstractVector{T}, data::AbstractVector{T}) where {B,T}
    out = _backend_array(B, T, length(data))
    _copy_build_data!(out, data)
    return out
end

function materialize_build(::GPUArrayBuildBackend{B}, ref::AbstractVector, data::BitVector) where {B}
    out = _backend_array(B, Bool, length(data))
    copyto!(out, Vector{Bool}(data))
    return out
end

singular_values_host(s::AbstractVector{T}) where {T} = Vector{T}(Array(s))

@inline use_host_build_algebra(s::AbstractArray) = use_host_build_algebra(execution_style(s))

function finalize_inverse_operator(backend::BuildBackend, ref::AbstractMatrix, F, inv_s_host::AbstractVector)
    if use_host_build_algebra(F.S)
        return materialize_build(backend, ref, F.V * Diagonal(inv_s_host) * F.U')
    end
    inv_s = materialize_build(backend, F.S, inv_s_host)
    return materialize_build(backend, ref, F.V * Diagonal(inv_s) * F.U')
end

inverse_operator(A::AbstractMatrix{T}, policy::InversePolicy) where {T<:AbstractFloat} =
    inverse_operator(default_build_backend(A), A, policy)

function inverse_operator(::BuildBackend, A::AbstractMatrix{T}, policy::InversePolicy) where {T<:AbstractFloat}
    throw(UnsupportedAlgorithm("inverse_operator is not implemented for $(typeof(policy))"))
end

function inverse_operator(backend::BuildBackend, A::AbstractMatrix{T}, ::ExactPseudoInverse) where {T<:AbstractFloat}
    F = svd(prepare_build_matrix(backend, A); full=false)
    s_host = singular_values_host(F.S)
    inv_s_host = similar(s_host)
    @inbounds for i in eachindex(s_host)
        inv_s_host[i] = iszero(s_host[i]) ? zero(T) : inv(s_host[i])
    end
    M = finalize_inverse_operator(backend, A, F, inv_s_host)
    effective_rank = count(!iszero, s_host)
    cond = effective_rank == 0 ? T(Inf) : s_host[begin] / s_host[effective_rank]
    return M, InverseStats(s_host, cond, effective_rank, 0)
end

function inverse_operator(backend::BuildBackend, A::AbstractMatrix{T}, policy::TSVDInverse) where {T<:AbstractFloat}
    policy.n_trunc >= 0 || throw(InvalidConfiguration("TSVD n_trunc must be >= 0"))
    F = svd(prepare_build_matrix(backend, A); full=false)
    s_host = singular_values_host(F.S)
    isempty(s_host) && return materialize_build(backend, similar(A, T, size(A, 2), size(A, 1))),
        InverseStats(s_host, T(Inf), 0, 0)
    cutoff = _inverse_cutoff(s_host, T(policy.rtol), T(policy.atol))
    rank_by_tol = count(>(cutoff), s_host)
    effective_rank = max(rank_by_tol - policy.n_trunc, 0)
    inv_s_host = similar(s_host)
    fill!(inv_s_host, zero(T))
    @inbounds for i in 1:effective_rank
        inv_s_host[i] = inv(s_host[i])
    end
    M = finalize_inverse_operator(backend, A, F, inv_s_host)
    cond = effective_rank == 0 ? T(Inf) : s_host[begin] / s_host[effective_rank]
    return M, InverseStats(s_host, cond, effective_rank, length(s_host) - effective_rank)
end

function inverse_operator(backend::BuildBackend, A::AbstractMatrix{T}, policy::TikhonovInverse) where {T<:AbstractFloat}
    policy.lambda >= 0 || throw(InvalidConfiguration("Tikhonov lambda must be >= 0"))
    F = svd(prepare_build_matrix(backend, A); full=false)
    s_host = singular_values_host(F.S)
    inv_s_host = similar(s_host)
    λ2 = T(policy.lambda)^2
    @inbounds for i in eachindex(s_host)
        inv_s_host[i] = s_host[i] / (s_host[i]^2 + λ2)
    end
    cutoff = _inverse_cutoff(s_host, T(policy.rtol), T(policy.atol))
    effective_rank = count(>(cutoff), s_host)
    denom = max(isempty(s_host) ? zero(T) : s_host[end], T(policy.lambda))
    cond = (isempty(s_host) || iszero(denom)) ? T(Inf) : s_host[begin] / denom
    M = finalize_inverse_operator(backend, A, F, inv_s_host)
    return M, InverseStats(s_host, cond, effective_rank, 0)
end
