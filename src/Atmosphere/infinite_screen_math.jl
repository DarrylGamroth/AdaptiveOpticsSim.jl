"""
Geometry for one-sided infinite-screen boundary injection.

The coordinate matrices are geometry coordinates, not Julia array indices.
Interior stencil coordinates live on the positive side of the boundary and the
boundary coordinates lie exactly one pixel outside the maintained screen.
"""
struct InfiniteBoundaryStencil{T<:AbstractFloat,I<:AbstractMatrix{Int},P<:AbstractMatrix{T}}
    stencil_coords::I
    boundary_coords::I
    stencil_positions::P
    boundary_positions::P
    orientation::Symbol
    side::Symbol
end

"""
Conditional Gaussian operator for one boundary update.

`predictor` is the linear conditional mean operator `A`.
`residual_factor` is a square-root factor `B` such that `B * ξ` reproduces the
conditional residual covariance for standard normal `ξ`.
"""
struct InfiniteBoundaryOperator{T<:AbstractFloat,M<:AbstractMatrix{T},V<:AbstractVector{T}}
    predictor::M
    residual_factor::M
    cov_zz::M
    cov_xx::M
    cov_xz::M
    cov_zx::M
    residual_covariance::M
    singular_values::V
    condition_ratio::T
    orientation::Symbol
    side::Symbol
end

struct InfiniteBoundaryModel{
    S<:InfiniteBoundaryStencil,
    O<:InfiniteBoundaryOperator,
}
    stencil::S
    operator::O
end

function validate_infinite_stencil(stencil_size::Int, screen_resolution::Int,
    orientation::Symbol, side::Symbol, tail_stride::Int)
    stencil_size >= 3 || throw(InvalidConfiguration("stencil_size must be >= 3"))
    isodd(stencil_size) || throw(InvalidConfiguration("stencil_size must be odd"))
    stencil_size > screen_resolution || throw(InvalidConfiguration("stencil_size must exceed screen_resolution"))
    orientation in (:column, :row) || throw(InvalidConfiguration("orientation must be :column or :row"))
    side in (:positive, :negative) || throw(InvalidConfiguration("side must be :positive or :negative"))
    tail_stride >= 1 || throw(InvalidConfiguration("tail_stride must be >= 1"))
    return nothing
end

function _append_unique_coord!(coords::Vector{NTuple{2,Int}}, coord::NTuple{2,Int})
    coord in coords || push!(coords, coord)
    return coords
end

function _column_stencil_coords(screen_resolution::Int, stencil_size::Int; tail_stride::Int=screen_resolution)
    coords = NTuple{2,Int}[]
    max_n = floor(Int, log2(stencil_size))
    for n in 0:max_n
        col = n == 0 ? 1 : Int(2^(n - 1) + 1)
        n_points = Int(2^(max_n - n) + 1)
        rows = round.(Int, range(1, stencil_size; length=n_points))
        for row in rows
            _append_unique_coord!(coords, (col, row))
        end
    end
    center_row = (stencil_size + 1) ÷ 2
    for col in tail_stride:tail_stride:(stencil_size - 1)
        _append_unique_coord!(coords, (col, center_row))
    end
    return reduce(vcat, (reshape(collect(coord), 1, 2) for coord in coords))
end

function infinite_boundary_stencil(screen_resolution::Int, pixel_scale::Real;
    stencil_size::Int=default_infinite_stencil_size(screen_resolution),
    orientation::Symbol=:column,
    side::Symbol=:positive,
    tail_stride::Int=screen_resolution,
    T::Type{<:AbstractFloat}=Float64)
    validate_infinite_stencil(stencil_size, screen_resolution, orientation, side, tail_stride)
    base_coords = _column_stencil_coords(screen_resolution, stencil_size; tail_stride=tail_stride)

    stencil_coords = copy(base_coords)
    boundary_coords = hcat(fill(0, stencil_size), collect(1:stencil_size))
    if side === :negative
        stencil_coords[:, 1] .= stencil_size .+ 1 .- stencil_coords[:, 1]
        boundary_coords[:, 1] .= stencil_size + 1
    end
    if orientation === :row
        stencil_coords = stencil_coords[:, [2, 1]]
        boundary_coords = boundary_coords[:, [2, 1]]
    end

    stencil_positions = T.(stencil_coords) .* T(pixel_scale)
    boundary_positions = T.(boundary_coords) .* T(pixel_scale)
    return InfiniteBoundaryStencil(
        stencil_coords,
        boundary_coords,
        stencil_positions,
        boundary_positions,
        orientation,
        side,
    )
end

function pairwise_separations(positions_a::AbstractMatrix{T}, positions_b::AbstractMatrix{T}) where {T<:AbstractFloat}
    size(positions_a, 2) == 2 || throw(DimensionMismatchError("positions_a must have shape (N, 2)"))
    size(positions_b, 2) == 2 || throw(DimensionMismatchError("positions_b must have shape (M, 2)"))
    out = Matrix{T}(undef, size(positions_a, 1), size(positions_b, 1))
    @inbounds for i in axes(positions_a, 1), j in axes(positions_b, 1)
        dx = positions_a[i, 1] - positions_b[j, 1]
        dy = positions_a[i, 2] - positions_b[j, 2]
        out[i, j] = sqrt(dx * dx + dy * dy)
    end
    return out
end

function boundary_injection_covariances(stencil::InfiniteBoundaryStencil, r0::Real, L0::Real)
    cov_zz = phase_covariance(pairwise_separations(stencil.stencil_positions, stencil.stencil_positions), r0, L0)
    cov_xx = phase_covariance(pairwise_separations(stencil.boundary_positions, stencil.boundary_positions), r0, L0)
    cov_xz = phase_covariance(pairwise_separations(stencil.boundary_positions, stencil.stencil_positions), r0, L0)
    cov_zx = phase_covariance(pairwise_separations(stencil.stencil_positions, stencil.boundary_positions), r0, L0)
    return (; cov_zz, cov_xx, cov_xz, cov_zx)
end

function _condition_ratio(values::AbstractVector{T}) where {T<:AbstractFloat}
    σmax = maximum(values)
    σmin = minimum(values)
    σmin > zero(T) || return T(Inf)
    return σmax / σmin
end

function boundary_injection_operator(stencil::InfiniteBoundaryStencil, r0::Real, L0::Real;
    conditioning_tol::Real=1e12,
    eigen_floor::Real=1e-10)
    T = eltype(stencil.stencil_positions)
    cov = boundary_injection_covariances(stencil, r0, L0)
    zz_singular_values = T.(svdvals(Symmetric(cov.cov_zz)))
    cond_ratio = _condition_ratio(zz_singular_values)
    isfinite(cond_ratio) || throw(NumericalConditionError("infinite-screen stencil covariance is singular"))
    cond_ratio <= T(conditioning_tol) || throw(NumericalConditionError("infinite-screen stencil covariance is ill-conditioned (condition ratio $(cond_ratio))"))

    predictor = cov.cov_xz / cov.cov_zz
    residual_covariance = (cov.cov_xx - predictor * cov.cov_zx)
    residual_covariance = (residual_covariance + transpose(residual_covariance)) ./ T(2)
    eig = eigen(Symmetric(residual_covariance))
    λmax = maximum(abs, eig.values)
    λmin = minimum(eig.values)
    λmin >= -T(eigen_floor) * max(λmax, one(T)) ||
        throw(NumericalConditionError("infinite-screen residual covariance is not positive semidefinite"))
    clamped_values = max.(eig.values, zero(T))
    residual_factor = eig.vectors * Diagonal(sqrt.(clamped_values))
    return InfiniteBoundaryOperator(
        Matrix{T}(predictor),
        Matrix{T}(residual_factor),
        Matrix{T}(cov.cov_zz),
        Matrix{T}(cov.cov_xx),
        Matrix{T}(cov.cov_xz),
        Matrix{T}(cov.cov_zx),
        Matrix{T}(residual_covariance),
        Vector{T}(clamped_values),
        cond_ratio,
        stencil.orientation,
        stencil.side,
    )
end

function boundary_injection_model(screen_resolution::Int, pixel_scale::Real, r0::Real, L0::Real;
    stencil_size::Int=default_infinite_stencil_size(screen_resolution),
    orientation::Symbol=:column,
    side::Symbol=:positive,
    tail_stride::Int=screen_resolution,
    conditioning_tol::Real=1e12,
    eigen_floor::Real=1e-10,
    T::Type{<:AbstractFloat}=Float64)
    stencil = infinite_boundary_stencil(screen_resolution, pixel_scale;
        stencil_size=stencil_size,
        orientation=orientation,
        side=side,
        tail_stride=tail_stride,
        T=T)
    op = boundary_injection_operator(stencil, r0, L0;
        conditioning_tol=conditioning_tol,
        eigen_floor=eigen_floor)
    return InfiniteBoundaryModel(stencil, op)
end

function sample_boundary_line(op::InfiniteBoundaryOperator{T}, stencil_data::AbstractVector{T},
    rng::AbstractRNG=Random.default_rng()) where {T<:AbstractFloat}
    length(stencil_data) == size(op.predictor, 2) ||
        throw(DimensionMismatchError("stencil_data length must match the predictor input dimension"))
    return op.predictor * stencil_data + op.residual_factor * randn(rng, T, size(op.residual_factor, 2))
end

function sample_boundary_line!(out::AbstractVector{T}, op::InfiniteBoundaryOperator{T},
    stencil_data::AbstractVector{T}, noise::AbstractVector{T},
    rng::AbstractRNG=Random.default_rng()) where {T<:AbstractFloat}
    length(stencil_data) == size(op.predictor, 2) ||
        throw(DimensionMismatchError("stencil_data length must match the predictor input dimension"))
    length(out) == size(op.predictor, 1) ||
        throw(DimensionMismatchError("output length must match the boundary dimension"))
    length(noise) == size(op.residual_factor, 2) ||
        throw(DimensionMismatchError("noise length must match the residual-factor input dimension"))
    mul!(out, op.predictor, stencil_data)
    randn_backend!(rng, noise)
    mul!(out, op.residual_factor, noise, one(T), one(T))
    synchronize_backend!(execution_style(out))
    return out
end
