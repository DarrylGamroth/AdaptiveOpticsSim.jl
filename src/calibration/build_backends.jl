# Cold numerical-product backend materialization.
abstract type BuildBackend end

struct NativeBuildBackend <: BuildBackend end
struct CPUBuildBackend <: BuildBackend end
struct GPUArrayBuildBackend{B} <: BuildBackend end

default_build_backend(::AbstractArray) = NativeBuildBackend()
function default_runtime_calibration_build_backend(A::AbstractArray)
    gpu_backend_name(typeof(A)) === nothing && return NativeBuildBackend()
    return CPUBuildBackend()
end
GPUArrayBuildBackend(::Type{B}) where {B} = GPUArrayBuildBackend{B}()

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

materialize_runtime_build_result(::CPUBuildBackend, ref, data) =
    materialize_build(NativeBuildBackend(), ref, data)
materialize_runtime_build_result(backend::BuildBackend, ref, data) =
    materialize_build(backend, ref, data)
