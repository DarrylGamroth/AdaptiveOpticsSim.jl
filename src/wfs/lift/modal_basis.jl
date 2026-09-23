@kernel function lift_basis_expansion_kernel!(opd, basis, coeffs, pupil, n_modes::Int)
    I = @index(Global, Cartesian)
    i, j = Tuple(I)
    if i <= size(opd, 1) && j <= size(opd, 2)
        value = zero(eltype(opd))
        @inbounds for k in 1:n_modes
            value += coeffs[k] * basis[i, j, k]
        end
        @inbounds opd[i, j] = ifelse(pupil[i, j], value, zero(value))
    end
end

@inline function lift_basis_expansion!(opd::AbstractMatrix{T},
    basis::AbstractArray{T,3}, coeffs::AbstractVector{T},
    pupil::AbstractMatrix{Bool}) where {T<:AbstractFloat}
    return lift_basis_expansion!(execution_style(opd), opd, basis, coeffs,
        pupil)
end

function lift_basis_expansion!(::ScalarCPUStyle, opd::AbstractMatrix{T},
    basis::AbstractArray{T,3}, coeffs::AbstractVector{T},
    pupil::AbstractMatrix{Bool}) where {T<:AbstractFloat}
    n_modes = min(size(basis, 3), length(coeffs))
    fill!(opd, zero(T))
    @inbounds for k in 1:n_modes
        @views @. opd += coeffs[k] * basis[:, :, k]
    end
    @. opd *= pupil
    return opd
end

function lift_basis_expansion!(style::AcceleratorStyle,
    opd::AbstractMatrix{T}, basis::AbstractArray{T,3},
    coeffs::AbstractVector{T}, pupil::AbstractMatrix{Bool}) where {T<:AbstractFloat}
    n_modes = min(size(basis, 3), length(coeffs))
    if iszero(n_modes)
        fill!(opd, zero(T))
        return opd
    end
    launch_kernel!(style, lift_basis_expansion_kernel!, opd, basis, coeffs,
        pupil, n_modes; ndrange=size(opd))
    return opd
end

@inline function lift_scaled_basis_mode!(dest::AbstractMatrix{T},
    amplitude::AbstractMatrix{T}, basis::AbstractArray{T,3},
    mode_id::Int, scale::T) where {T<:AbstractFloat}
    return lift_scaled_basis_mode!(execution_style(dest), dest, amplitude,
        basis, mode_id, scale)
end

function lift_scaled_basis_mode!(::ScalarCPUStyle,
    dest::AbstractMatrix{T}, amplitude::AbstractMatrix{T},
    basis::AbstractArray{T,3}, mode_id::Int,
    scale::T) where {T<:AbstractFloat}
    mode_offset = (mode_id - 1) * length(dest)
    @inbounds @simd for i in eachindex(dest, amplitude)
        dest[i] = amplitude[i] * scale * basis[i + mode_offset]
    end
    return dest
end

function lift_scaled_basis_mode!(style::AcceleratorStyle,
    dest::AbstractMatrix{T}, amplitude::AbstractMatrix{T},
    basis::AbstractArray{T,3}, mode_id::Int,
    scale::T) where {T<:AbstractFloat}
    n = length(dest)
    mode_offset = (mode_id - 1) * n
    launch_kernel!(style, lift_scaled_basis_mode_kernel!, dest, amplitude,
        basis, scale, mode_offset, n; ndrange=n)
    return dest
end

@inline function lift_copy_column!(dest::AbstractMatrix,
    column::Int, src::AbstractMatrix)
    return lift_copy_column!(execution_style(dest), dest, column, src)
end

function lift_copy_column!(::ScalarCPUStyle, dest::AbstractMatrix{T},
    column::Int, src::AbstractMatrix{T}) where {T}
    @inbounds @simd for i in eachindex(src)
        dest[i, column] = src[i]
    end
    return dest
end

function lift_copy_column!(style::AcceleratorStyle,
    dest::AbstractMatrix, column::Int, src::AbstractMatrix)
    launch_kernel!(style, lift_copy_column_kernel!, dest, column, src,
        length(src); ndrange=length(src))
    return dest
end
