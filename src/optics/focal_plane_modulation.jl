#
# Prepared focal-plane modulation
#
# Pyramid and BioEdge sensors may use the same pupil-plane tip/tilt path while
# retaining independent focal-plane masks and propagation workspaces. The
# weights below are optical quadrature weights: they average instantaneous
# intensities over one modulation cycle and never represent elapsed exposure.
#

abstract type AbstractFocalPlaneModulation end

"""A stationary focal-plane sensor with one unit-weight optical sample."""
struct NoModulation <: AbstractFocalPlaneModulation end

@kernel function circular_modulation_phases_kernel!(
    phases, radius, phase_offset, coordinate_start, coordinate_step,
    point_count::Int, resolution::Int)
    axis_1, axis_2, point = @index(Global, NTuple)
    if axis_1 <= resolution && axis_2 <= resolution && point <= point_count
        angle = phase_offset +
            typeof(radius)(2pi) * typeof(radius)(point - 1) /
                typeof(radius)(point_count)
        sine, cosine = sincos(angle)
        offset_x = radius * cosine
        offset_y = radius * sine
        coordinate_1 = coordinate_start +
            typeof(radius)(axis_1 - 1) * coordinate_step
        coordinate_2 = coordinate_start +
            typeof(radius)(axis_2 - 1) * coordinate_step
        @inbounds phases[axis_1, axis_2, point] = cis(
            offset_x * coordinate_2 + offset_y * coordinate_1)
    end
end

"""
    CircularModulation(radius; samples, phase_offset=0)

Uniform circular modulation in units of `lambda / D`. `samples` is the number
of equally weighted points over one cycle and `phase_offset` is in radians.
"""
struct CircularModulation{T<:AbstractFloat} <: AbstractFocalPlaneModulation
    radius::T
    samples::Int
    phase_offset::T

    function CircularModulation(radius::T, samples::Int,
        phase_offset::T) where {T<:AbstractFloat}
        isfinite(radius) && radius >= zero(T) || throw(InvalidConfiguration(
            "circular modulation radius must be finite and nonnegative"))
        samples >= 1 || throw(InvalidConfiguration(
            "circular modulation samples must be >= 1"))
        isfinite(phase_offset) || throw(InvalidConfiguration(
            "circular modulation phase offset must be finite"))
        return new{T}(radius, samples, phase_offset)
    end
end

function CircularModulation(radius::Real; samples::Int,
    phase_offset::Real=0, T::Type{<:AbstractFloat}=Float64)
    return CircularModulation(T(radius), samples, T(phase_offset))
end

"""
    SampledModulation(points; weights=nothing)

An arbitrary focal-plane modulation path. Each point is an `(x, y)` offset in
units of `lambda / D`. Weights are copied, normalized, finite, nonnegative
optical quadrature weights; at least one must be positive.
"""
struct SampledModulation{T<:AbstractFloat} <: AbstractFocalPlaneModulation
    points::Vector{NTuple{2,T}}
    weights::Vector{T}
end

function SampledModulation(points; weights=nothing,
    T::Type{<:AbstractFloat}=Float64)
    isempty(points) && throw(InvalidConfiguration(
        "sampled modulation must contain at least one point"))
    owned_points = Vector{NTuple{2,T}}(undef, length(points))
    @inbounds for index in eachindex(points)
        point = points[index]
        length(point) >= 2 || throw(InvalidConfiguration(
            "each sampled modulation point must contain x and y offsets"))
        x = T(point[1])
        y = T(point[2])
        isfinite(x) && isfinite(y) || throw(InvalidConfiguration(
            "sampled modulation offsets must be finite"))
        owned_points[index] = (x, y)
    end
    raw_weights = if weights === nothing
        inferred = Vector{T}(undef, length(points))
        @inbounds for index in eachindex(points)
            point = points[index]
            inferred[index] = length(point) >= 3 ? T(point[3]) : one(T)
        end
        inferred
    else
        length(weights) == length(points) || throw(InvalidConfiguration(
            "sampled modulation weights must match the point count"))
        T.(collect(weights))
    end
    weight_sum = zero(T)
    @inbounds for weight in raw_weights
        isfinite(weight) && weight >= zero(T) || throw(InvalidConfiguration(
            "sampled modulation weights must be finite and nonnegative"))
        weight_sum += weight
    end
    weight_sum > zero(T) || throw(InvalidConfiguration(
        "sampled modulation requires at least one positive weight"))
    raw_weights ./= weight_sum
    return SampledModulation{T}(owned_points, raw_weights)
end

"""Backend-bound modulation phases and host-resident amplitude weights."""
struct PreparedFocalPlaneModulation{P<:AbstractFocalPlaneModulation,A,V}
    policy::P
    phases::A
    amplitude_weights::V
end

@inline modulation_point_count(::NoModulation) = 1
@inline modulation_point_count(policy::CircularModulation) = policy.samples
@inline modulation_point_count(policy::SampledModulation) = length(policy.points)
@inline modulation_point_count(prepared::PreparedFocalPlaneModulation) =
    size(prepared.phases, 3)

@inline modulation_offset(::NoModulation, ::Int, ::Type{T}) where {T} =
    (zero(T), zero(T))

@inline function modulation_offset(policy::CircularModulation, index::Int,
    ::Type{T}) where {T}
    angle = T(policy.phase_offset) + T(2pi * (index - 1) / policy.samples)
    sine, cosine = sincos(angle)
    radius = T(policy.radius)
    return radius * cosine, radius * sine
end

@inline function modulation_offset(policy::SampledModulation, index::Int,
    ::Type{T}) where {T}
    point = @inbounds policy.points[index]
    return T(point[1]), T(point[2])
end

@inline modulation_weight(::NoModulation, ::Int, ::Type{T}) where {T} = one(T)
@inline modulation_weight(policy::CircularModulation, ::Int,
    ::Type{T}) where {T} = inv(T(policy.samples))
@inline modulation_weight(policy::SampledModulation, index::Int,
    ::Type{T}) where {T} = T(@inbounds policy.weights[index])

function prepare_focal_plane_modulation(policy::AbstractFocalPlaneModulation,
    resolution::Int, backend_storage::AbstractArray, ::Type{T}) where {
    T<:AbstractFloat,
}
    resolution >= 1 || throw(InvalidConfiguration(
        "modulation pupil resolution must be positive"))
    point_count = modulation_point_count(policy)
    phases = similar(backend_storage, Complex{T}, resolution, resolution,
        point_count)
    host_phases = Array{Complex{T}}(undef, resolution, resolution, point_count)
    amplitude_weights = Vector{T}(undef, point_count)
    coordinates = range(-T(pi), T(pi); length=resolution)
    @inbounds for point_index in 1:point_count
        offset_x, offset_y = modulation_offset(policy, point_index, T)
        amplitude_weights[point_index] = sqrt(modulation_weight(policy,
            point_index, T))
        for axis_2 in 1:resolution, axis_1 in 1:resolution
            phase = offset_x * coordinates[axis_2] +
                offset_y * coordinates[axis_1]
            host_phases[axis_1, axis_2, point_index] = cis(phase)
        end
    end
    copyto!(phases, host_phases)
    return PreparedFocalPlaneModulation(policy, phases, amplitude_weights)
end

@inline function _circular_modulation_coordinate_step(
    ::Type{T}, resolution::Int) where {T<:AbstractFloat}
    resolution == 1 && return zero(T)
    return T(2pi) / T(resolution - 1)
end

function _update_cycle_averaged_circular_modulation!(
    ::ScalarCPUStyle, prepared::PreparedFocalPlaneModulation,
    radius::T, phase_offset::T) where {T<:AbstractFloat}
    phases = prepared.phases
    resolution = size(phases, 1)
    point_count = size(phases, 3)
    coordinates = range(-T(pi), T(pi); length=resolution)
    @inbounds for point in axes(phases, 3)
        angle = phase_offset + T(2pi * (point - 1) / point_count)
        sine, cosine = sincos(angle)
        offset_x = radius * cosine
        offset_y = radius * sine
        for axis_2 in axes(phases, 2), axis_1 in axes(phases, 1)
            phases[axis_1, axis_2, point] = cis(
                offset_x * coordinates[axis_2] +
                offset_y * coordinates[axis_1])
        end
    end
    return prepared
end

function _update_cycle_averaged_circular_modulation!(
    style::AcceleratorStyle, prepared::PreparedFocalPlaneModulation,
    radius::T, phase_offset::T) where {T<:AbstractFloat}
    phases = prepared.phases
    resolution = size(phases, 1)
    point_count = size(phases, 3)
    coordinate_start = -T(pi)
    coordinate_step = _circular_modulation_coordinate_step(T, resolution)
    launch_kernel!(style, circular_modulation_phases_kernel!, phases,
        radius, phase_offset, coordinate_start, coordinate_step,
        point_count, resolution; ndrange=size(phases))
    return prepared
end

"""
    update_cycle_averaged_circular_modulation!(
        prepared, radius; enabled=true)

Regenerate the already allocated circular focal-plane quadrature at `radius`
without changing its point count, normalized weights, backend, or numerical
quadrature origin. A disabled waveform is centered by setting every prepared
quadrature point to zero offset. This operation does not integrate detector
time or evaluate trigger-relative waveform phase.
"""
function update_cycle_averaged_circular_modulation!(
    prepared::PreparedFocalPlaneModulation{<:CircularModulation},
    radius::T; enabled::Bool=true) where {T<:AbstractFloat}
    isfinite(radius) && radius >= zero(T) || throw(InvalidConfiguration(
        "cycle-averaged circular modulation radius must be finite and nonnegative"))
    eltype(prepared.amplitude_weights) === T || throw(InvalidConfiguration(
        "cycle-averaged circular modulation radius precision must match its prepared weights"))
    size(prepared.phases, 1) == size(prepared.phases, 2) || throw(
        InvalidConfiguration(
            "cycle-averaged circular modulation phases must use a square pupil grid"))
    size(prepared.phases, 3) == prepared.policy.samples || throw(
        InvalidConfiguration(
            "cycle-averaged circular modulation point count changed after preparation"))
    effective_radius = enabled ? radius : zero(T)
    return _update_cycle_averaged_circular_modulation!(
        execution_style(prepared.phases), prepared, effective_radius,
        T(prepared.policy.phase_offset))
end
