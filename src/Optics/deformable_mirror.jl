using LinearAlgebra

#
# Deformable mirror model
#
# The DM surface is represented as a linear superposition of sampled influence
# functions:
#
#   opd(x, y) = sum_k modes(x, y, k) * actuator_coefs[k]
#
# The implementation now separates:
# - actuator topology
# - static influence basis
# - actuator behavior
#
# The common Gaussian/no-rotation/no-anamorphosis/grid case still uses the
# faster separable
#
#   opd(x, y) = X * C * Y'
#
# path, while richer sampled-basis models fall back to the dense matrix-vector
# application.
#
@kernel function dm_mode_kernel!(mode, pupil, x_m, y_m, cx, cy, scale, sigma2, n::Int)
    i, j = @index(Global, NTuple)
    if i <= n && j <= n
        x = (i - cx) / scale
        y = (j - cy) / scale
        r2 = (x - x_m)^2 + (y - y_m)^2
        val = exp(-r2 / (2 * sigma2))
        @inbounds mode[(j - 1) * n + i] = ifelse(pupil[i, j], val, zero(eltype(mode)))
    end
end

@kernel function dm_separable_tmp_kernel!(tmp, xbasis, coefs_grid, n_act::Int)
    i, a = @index(Global, NTuple)
    if i <= size(tmp, 1) && a <= size(tmp, 2)
        acc = zero(eltype(tmp))
        @inbounds for b in 1:n_act
            acc = muladd(xbasis[i, b], coefs_grid[b, a], acc)
        end
        tmp[i, a] = acc
    end
end

@kernel function dm_separable_finalize_kernel!(opd, tmp, ybasis_t, pupil, n_act::Int)
    i, j = @index(Global, NTuple)
    if i <= size(opd, 1) && j <= size(opd, 2)
        acc = zero(eltype(opd))
        @inbounds for a in 1:n_act
            acc = muladd(tmp[i, a], ybasis_t[a, j], acc)
        end
        opd[i, j] = ifelse(pupil[i, j], acc, zero(eltype(opd)))
    end
end

abstract type AbstractDMTopology end
abstract type AbstractDMInfluenceModel end
abstract type AbstractDMActuatorModel end

struct ActuatorGridTopology{
    T<:AbstractFloat,
    M1<:AbstractMatrix{T},
    M2<:AbstractMatrix{T},
    V1<:AbstractVector{Bool},
    V2<:AbstractVector{Int},
    MD,
} <: AbstractDMTopology
    n_act::Int
    coords::M1
    active_coords::M2
    valid_actuators::V1
    active_indices::V2
    metadata::MD
end

struct SampledActuatorTopology{
    T<:AbstractFloat,
    M1<:AbstractMatrix{T},
    M2<:AbstractMatrix{T},
    V1<:AbstractVector{Bool},
    V2<:AbstractVector{Int},
    MD,
} <: AbstractDMTopology
    coords::M1
    active_coords::M2
    valid_actuators::V1
    active_indices::V2
    metadata::MD
end

struct GaussianInfluenceWidth{T<:AbstractFloat} <: AbstractDMInfluenceModel
    width::T
end

struct GaussianMechanicalCoupling{T<:AbstractFloat} <: AbstractDMInfluenceModel
    coupling::T
end

struct DenseInfluenceMatrix{T<:AbstractFloat,M<:AbstractMatrix{T}} <: AbstractDMInfluenceModel
    modes::M
end

struct MeasuredInfluenceFunctions{T<:AbstractFloat,M<:AbstractMatrix{T},MD} <: AbstractDMInfluenceModel
    modes::M
    metadata::MD
end

struct LinearStaticActuators <: AbstractDMActuatorModel end

struct ClippedActuators{T<:AbstractFloat} <: AbstractDMActuatorModel
    lo::T
    hi::T
end

struct ActuatorHealthMap{T<:AbstractFloat,V<:AbstractVector{T}} <: AbstractDMActuatorModel
    gains::V
end

struct CompositeDMActuatorModel{A<:Tuple} <: AbstractDMActuatorModel
    stages::A
end

GaussianInfluenceWidth(width::Real) = GaussianInfluenceWidth(float(width))
GaussianMechanicalCoupling(coupling::Real) = GaussianMechanicalCoupling(float(coupling))
DenseInfluenceMatrix(modes::AbstractMatrix{<:Real}) = DenseInfluenceMatrix(Matrix{float(eltype(modes))}(modes))
MeasuredInfluenceFunctions(modes::AbstractMatrix{<:Real}; metadata=(;)) =
    MeasuredInfluenceFunctions(Matrix{float(eltype(modes))}(modes), metadata)
ClippedActuators(lo::Real, hi::Real) = ClippedActuators(float(lo), float(hi))
ActuatorHealthMap(gains::AbstractVector{<:Real}) = ActuatorHealthMap(collect(float.(gains)))
CompositeDMActuatorModel(stages::AbstractDMActuatorModel...) = CompositeDMActuatorModel(stages)

function ActuatorGridTopology(n_act::Integer; valid_actuators::Union{Nothing,AbstractVector{Bool}}=nothing,
    metadata=(;), T::Type{<:AbstractFloat}=Float64)
    n_act > 0 || throw(InvalidConfiguration("ActuatorGridTopology n_act must be > 0"))
    total = Int(n_act) * Int(n_act)
    mask = valid_actuators === nothing ? fill(true, total) : collect(valid_actuators)
    length(mask) == total ||
        throw(DimensionMismatchError("ActuatorGridTopology valid_actuators length $(length(mask)) does not match expected $(total)"))
    any(mask) || throw(InvalidConfiguration("ActuatorGridTopology must contain at least one active actuator"))
    xs = collect(range(T(-1), T(1); length=Int(n_act)))
    ys = collect(range(T(-1), T(1); length=Int(n_act)))
    coords = Matrix{T}(undef, 2, total)
    idx = 1
    @inbounds for x0 in xs, y0 in ys
        coords[1, idx] = x0
        coords[2, idx] = y0
        idx += 1
    end
    active_idx = findall(identity, mask)
    active_coords = coords[:, active_idx]
    return ActuatorGridTopology{T,typeof(coords),typeof(active_coords),typeof(mask),typeof(active_idx),typeof(metadata)}(
        Int(n_act), coords, active_coords, mask, active_idx, metadata)
end

function SampledActuatorTopology(coords::AbstractMatrix{<:Real};
    valid_actuators::Union{Nothing,AbstractVector{Bool}}=nothing, metadata=(;),
    T::Type{<:AbstractFloat}=Float64)
    size(coords, 1) == 2 ||
        throw(DimensionMismatchError("SampledActuatorTopology coords must have size (2, n), got $(size(coords))"))
    total = size(coords, 2)
    total > 0 || throw(InvalidConfiguration("SampledActuatorTopology must contain at least one actuator"))
    mask = valid_actuators === nothing ? fill(true, total) : collect(valid_actuators)
    length(mask) == total ||
        throw(DimensionMismatchError("SampledActuatorTopology valid_actuators length $(length(mask)) does not match expected $(total)"))
    any(mask) || throw(InvalidConfiguration("SampledActuatorTopology must contain at least one active actuator"))
    coords_t = Matrix{T}(coords)
    active_idx = findall(identity, mask)
    active_coords = coords_t[:, active_idx]
    return SampledActuatorTopology{T,typeof(coords_t),typeof(active_coords),typeof(mask),typeof(active_idx),typeof(metadata)}(
        coords_t, active_coords, mask, active_idx, metadata)
end

struct DeformableMirrorParams{
    T<:AbstractFloat,
    I<:AbstractDMInfluenceModel,
    TP<:AbstractDMTopology,
    AM<:AbstractDMActuatorModel,
}
    n_act::Int
    topology::TP
    influence_model::I
    actuator_model::AM
    misregistration::Misregistration{T}
end

mutable struct DeformableMirrorState{
    T<:AbstractFloat,
    A<:AbstractMatrix{T},
    M<:AbstractMatrix{T},
    V<:AbstractVector{T},
    C,
}
    opd::A
    opd_vec::V
    modes::M
    coefs::V
    actuator_coefs::V
    coefs_grid::C
    separable_x::Union{Nothing,A}
    separable_y_t::Union{Nothing,A}
    separable_tmp::Union{Nothing,A}
end

struct DeformableMirror{P<:DeformableMirrorParams,S<:DeformableMirrorState,B<:AbstractArrayBackend} <: AbstractDeformableMirror
    params::P
    state::S
end

@inline backend(::DeformableMirror{<:Any,<:Any,B}) where {B} = B()

@inline command_storage(dm::DeformableMirror) = dm.state.coefs
@inline surface_opd(dm::DeformableMirror) = dm.state.opd
@inline command_layout(dm::DeformableMirror) = RuntimeCommandLayout(:dm => length(dm.state.coefs))
@inline influence_model(dm::DeformableMirror) = dm.params.influence_model
@inline topology(dm::DeformableMirror) = dm.params.topology
@inline actuator_model(dm::DeformableMirror) = dm.params.actuator_model
@inline n_actuators(dm::DeformableMirror) = length(dm.state.coefs)

@inline topology_axis_count(::AbstractDMTopology) = 0
@inline topology_axis_count(topology::ActuatorGridTopology) = topology.n_act
@inline topology_command_count(topology::AbstractDMTopology) = size(actuator_coordinates(topology), 2)
@inline actuator_coordinates(topology::ActuatorGridTopology) = topology.active_coords
@inline actuator_coordinates(topology::SampledActuatorTopology) = topology.active_coords
@inline valid_actuator_mask(topology::ActuatorGridTopology) = topology.valid_actuators
@inline valid_actuator_mask(topology::SampledActuatorTopology) = topology.valid_actuators
@inline active_actuator_indices(topology::ActuatorGridTopology) = topology.active_indices
@inline active_actuator_indices(topology::SampledActuatorTopology) = topology.active_indices
@inline topology_metadata(topology::ActuatorGridTopology) = topology.metadata
@inline topology_metadata(topology::SampledActuatorTopology) = topology.metadata
@inline actuator_coordinates(dm::DeformableMirror) = actuator_coordinates(topology(dm))
@inline valid_actuator_mask(dm::DeformableMirror) = valid_actuator_mask(topology(dm))
@inline active_actuator_indices(dm::DeformableMirror) = active_actuator_indices(topology(dm))
@inline topology_metadata(dm::DeformableMirror) = topology_metadata(topology(dm))

@inline supports_separable_topology(::AbstractDMTopology) = false
@inline supports_separable_topology(topology::ActuatorGridTopology) = all(valid_actuator_mask(topology))

@inline dm_normalized_pitch(n_act::Integer) = 2.0 / (n_act - 1)
@inline dm_normalized_pitch(topology::ActuatorGridTopology) = dm_normalized_pitch(topology.n_act)

function mechanical_coupling(n_act::Integer, influence_width::Real)
    n_act > 1 || throw(InvalidConfiguration("mechanical coupling requires n_act > 1"))
    influence_width > 0 || throw(InvalidConfiguration("influence_width must be positive"))
    pitch = dm_normalized_pitch(n_act)
    return exp(-(pitch^2) / (2 * influence_width^2))
end

function influence_width_from_mechanical_coupling(n_act::Integer, coupling::Real)
    n_act > 1 || throw(InvalidConfiguration("mechanical coupling requires n_act > 1"))
    zero(coupling) < coupling < one(coupling) ||
        throw(InvalidConfiguration("mechanical_coupling must lie in (0, 1)"))
    pitch = dm_normalized_pitch(n_act)
    return pitch / sqrt(-2 * log(float(coupling)))
end

@inline influence_width(model::GaussianInfluenceWidth, ::AbstractDMTopology) = model.width
@inline influence_width(model::GaussianInfluenceWidth, ::Integer) = model.width
@inline function influence_width(model::GaussianMechanicalCoupling, topology::ActuatorGridTopology)
    typeof(model.coupling)(influence_width_from_mechanical_coupling(topology.n_act, model.coupling))
end
@inline influence_width(model::GaussianMechanicalCoupling, n_act::Integer) =
    typeof(model.coupling)(influence_width_from_mechanical_coupling(n_act, model.coupling))
function influence_width(::GaussianMechanicalCoupling, ::AbstractDMTopology)
    throw(UnsupportedAlgorithm("GaussianMechanicalCoupling requires a grid-backed DM topology"))
end
influence_width(::DenseInfluenceMatrix, ::AbstractDMTopology) =
    throw(UnsupportedAlgorithm("influence_width is only defined for Gaussian DM influence models"))
influence_width(::MeasuredInfluenceFunctions, ::AbstractDMTopology) =
    throw(UnsupportedAlgorithm("influence_width is only defined for Gaussian DM influence models"))
@inline influence_width(dm::DeformableMirror) = influence_width(influence_model(dm), topology(dm))

@inline mechanical_coupling(model::GaussianInfluenceWidth, topology::ActuatorGridTopology) =
    typeof(model.width)(mechanical_coupling(topology.n_act, model.width))
@inline mechanical_coupling(model::GaussianInfluenceWidth, n_act::Integer) =
    typeof(model.width)(mechanical_coupling(n_act, model.width))
@inline mechanical_coupling(model::GaussianMechanicalCoupling, ::ActuatorGridTopology) = model.coupling
@inline mechanical_coupling(model::GaussianMechanicalCoupling, ::Integer) = model.coupling
function mechanical_coupling(::GaussianMechanicalCoupling, ::AbstractDMTopology)
    throw(UnsupportedAlgorithm("GaussianMechanicalCoupling requires a grid-backed DM topology"))
end
mechanical_coupling(::DenseInfluenceMatrix, ::AbstractDMTopology) =
    throw(UnsupportedAlgorithm("mechanical_coupling is only defined for Gaussian DM influence models"))
mechanical_coupling(::MeasuredInfluenceFunctions, ::AbstractDMTopology) =
    throw(UnsupportedAlgorithm("mechanical_coupling is only defined for Gaussian DM influence models"))
@inline mechanical_coupling(dm::DeformableMirror) = mechanical_coupling(influence_model(dm), topology(dm))

function validate_dm_actuator_model(::LinearStaticActuators)
    return LinearStaticActuators()
end

function validate_dm_actuator_model(model::ClippedActuators)
    model.lo <= model.hi ||
        throw(InvalidConfiguration("ClippedActuators requires lo <= hi"))
    return model
end

function validate_dm_actuator_model(model::ActuatorHealthMap)
    minimum(model.gains) >= zero(eltype(model.gains)) ||
        throw(InvalidConfiguration("ActuatorHealthMap gains must be >= 0"))
    return model
end

function validate_dm_actuator_model(model::CompositeDMActuatorModel)
    length(model.stages) > 0 ||
        throw(InvalidConfiguration("CompositeDMActuatorModel requires at least one stage"))
    map(validate_dm_actuator_model, model.stages)
    return model
end

function _convert_dm_actuator_model(::LinearStaticActuators, ::Type{T}, backend, n_commands::Integer) where {T<:AbstractFloat}
    return LinearStaticActuators()
end

function _convert_dm_actuator_model(model::ClippedActuators, ::Type{T}, backend, n_commands::Integer) where {T<:AbstractFloat}
    return ClippedActuators{T}(T(model.lo), T(model.hi))
end

function _convert_dm_actuator_model(model::ActuatorHealthMap, ::Type{T}, backend, n_commands::Integer) where {T<:AbstractFloat}
    length(model.gains) == n_commands ||
        throw(DimensionMismatchError("ActuatorHealthMap length $(length(model.gains)) does not match expected $(n_commands)"))
    gains = backend{T}(undef, n_commands)
    copyto!(gains, T.(Array(model.gains)))
    return ActuatorHealthMap{T,typeof(gains)}(gains)
end

function _convert_dm_actuator_model(model::CompositeDMActuatorModel, ::Type{T}, backend, n_commands::Integer) where {T<:AbstractFloat}
    stages = map(stage -> _convert_dm_actuator_model(stage, T, backend, n_commands), model.stages)
    return CompositeDMActuatorModel{typeof(stages)}(stages)
end

function _resolve_dm_actuator_model(model::Union{Nothing,AbstractDMActuatorModel},
    ::Type{T}, backend, n_commands::Integer) where {T<:AbstractFloat}
    resolved = model === nothing ? LinearStaticActuators() : model
    return validate_dm_actuator_model(_convert_dm_actuator_model(resolved, T, backend, n_commands))
end

@inline _is_identity_misregistration(mis::Misregistration) =
    iszero(mis.shift_x) &&
    iszero(mis.shift_y) &&
    iszero(mis.rotation_rad) &&
    iszero(mis.anamorphosis_angle_rad) &&
    mis.tangential_scaling == one(mis.tangential_scaling) &&
    mis.radial_scaling == one(mis.radial_scaling)

function _resolve_dm_topology(; n_act::Union{Nothing,Integer}, topology::Union{Nothing,AbstractDMTopology},
    T::Type{<:AbstractFloat})
    if topology !== nothing
        n_act === nothing || topology_axis_count(topology) == 0 || topology_axis_count(topology) == Int(n_act) ||
            throw(InvalidConfiguration("n_act does not match the supplied DM topology"))
        return topology
    end
    n_act === nothing &&
        throw(InvalidConfiguration("specify either n_act or topology when constructing a DeformableMirror"))
    return ActuatorGridTopology(Int(n_act); T=T)
end

function _convert_dm_influence_model(model::GaussianInfluenceWidth, ::Type{T}, backend,
    n::Integer, topology::AbstractDMTopology, misregistration::Misregistration) where {T<:AbstractFloat}
    return GaussianInfluenceWidth{T}(T(model.width))
end

function _convert_dm_influence_model(model::GaussianMechanicalCoupling, ::Type{T}, backend,
    n::Integer, topology::AbstractDMTopology, misregistration::Misregistration) where {T<:AbstractFloat}
    return GaussianMechanicalCoupling{T}(T(model.coupling))
end

function _convert_dm_influence_model(model::DenseInfluenceMatrix, ::Type{T}, backend,
    n::Integer, topology::AbstractDMTopology, misregistration::Misregistration) where {T<:AbstractFloat}
    _is_identity_misregistration(misregistration) ||
        throw(InvalidConfiguration("DenseInfluenceMatrix already encodes the sampled DM surface and requires identity misregistration"))
    expected = (n * n, topology_command_count(topology))
    size(model.modes) == expected ||
        throw(DimensionMismatchError("DenseInfluenceMatrix size $(size(model.modes)) does not match expected $(expected)"))
    modes = backend{T}(undef, expected...)
    copyto!(modes, T.(Array(model.modes)))
    return DenseInfluenceMatrix{T,typeof(modes)}(modes)
end

function _convert_dm_influence_model(model::MeasuredInfluenceFunctions, ::Type{T}, backend,
    n::Integer, topology::AbstractDMTopology, misregistration::Misregistration) where {T<:AbstractFloat}
    _is_identity_misregistration(misregistration) ||
        throw(InvalidConfiguration("MeasuredInfluenceFunctions already encode the sampled DM surface and require identity misregistration"))
    expected = (n * n, topology_command_count(topology))
    size(model.modes) == expected ||
        throw(DimensionMismatchError("MeasuredInfluenceFunctions size $(size(model.modes)) does not match expected $(expected)"))
    modes = backend{T}(undef, expected...)
    copyto!(modes, T.(Array(model.modes)))
    return MeasuredInfluenceFunctions{T,typeof(modes),typeof(model.metadata)}(modes, model.metadata)
end

function _resolve_dm_influence_model(; n::Integer, topology::AbstractDMTopology,
    influence_width::Union{Nothing,Real}, mechanical_coupling::Union{Nothing,Real},
    influence_model::Union{Nothing,AbstractDMInfluenceModel}, T::Type{<:AbstractFloat},
    backend, misregistration::Misregistration)

    if influence_model !== nothing
        influence_width === nothing ||
            throw(InvalidConfiguration("specify either influence_model or influence_width, not both"))
        mechanical_coupling === nothing ||
            throw(InvalidConfiguration("specify either influence_model or mechanical_coupling, not both"))
        return _convert_dm_influence_model(influence_model, T, backend, n, topology, misregistration)
    end
    if mechanical_coupling !== nothing
        influence_width === nothing ||
            throw(InvalidConfiguration("specify either influence_width or mechanical_coupling, not both"))
        return GaussianMechanicalCoupling{T}(T(mechanical_coupling))
    end
    resolved_width = influence_width === nothing ? T(0.2) : T(influence_width)
    return GaussianInfluenceWidth{T}(resolved_width)
end

"""
    DeformableMirror(tel; ...)

Build a DM whose surface is a weighted sum of sampled influence functions on
the telescope pupil grid.

The maintained concise constructor remains the grid-based Gaussian path:

```julia
DeformableMirror(tel; n_act=16, influence_width=0.3)
```

For richer models, supply a topology and an explicit influence model.
"""
function DeformableMirror(tel::Telescope; n_act::Union{Nothing,Int}=nothing,
    topology::Union{Nothing,AbstractDMTopology}=nothing,
    influence_width::Union{Nothing,Real}=nothing, mechanical_coupling::Union{Nothing,Real}=nothing,
    influence_model::Union{Nothing,AbstractDMInfluenceModel}=nothing,
    actuator_model::Union{Nothing,AbstractDMActuatorModel}=nothing,
    T::Type{<:AbstractFloat}=Float64, misregistration::Misregistration=Misregistration(T=T),
    backend::AbstractArrayBackend=backend(tel))
    selector = require_same_backend(tel, _resolve_backend_selector(backend))
    backend = _resolve_array_backend(selector)
    n = tel.params.resolution
    resolved_topology = _resolve_dm_topology(n_act=n_act, topology=topology, T=T)
    n_commands = topology_command_count(resolved_topology)
    resolved_model = _resolve_dm_influence_model(
        n=n,
        topology=resolved_topology,
        influence_width=influence_width,
        mechanical_coupling=mechanical_coupling,
        influence_model=influence_model,
        T=T,
        backend=backend,
        misregistration=misregistration,
    )
    resolved_actuator_model = _resolve_dm_actuator_model(actuator_model, T, backend, n_commands)
    params = DeformableMirrorParams{T,typeof(resolved_model),typeof(resolved_topology),typeof(resolved_actuator_model)}(
        topology_axis_count(resolved_topology), resolved_topology, resolved_model, resolved_actuator_model, misregistration)
    opd = backend{T}(undef, n, n)
    fill!(opd, zero(T))
    opd_vec = reshape(opd, :)
    modes = backend{T}(undef, n * n, n_commands)
    coefs = backend{T}(undef, n_commands)
    actuator_coefs = similar(coefs)
    fill!(coefs, zero(T))
    fill!(actuator_coefs, zero(T))
    coefs_grid = supports_separable_topology(resolved_topology) ?
        reshape(actuator_coefs, resolved_topology.n_act, resolved_topology.n_act) :
        nothing
    state = DeformableMirrorState{T,typeof(opd),typeof(modes),typeof(opd_vec),typeof(coefs_grid)}(
        opd, opd_vec, modes, coefs, actuator_coefs, coefs_grid, nothing, nothing, nothing)
    dm = DeformableMirror{typeof(params),typeof(state),typeof(selector)}(params, state)
    build_influence_functions!(dm, tel)
    build_separable_influence!(dm, tel)
    return dm
end

@inline function supports_separable_influence(mis::Misregistration)
    return iszero(mis.rotation_rad) &&
           iszero(mis.anamorphosis_angle_rad) &&
           mis.tangential_scaling == one(mis.tangential_scaling) &&
           mis.radial_scaling == one(mis.radial_scaling)
end

@inline supports_separable_influence(::AbstractDMInfluenceModel, ::AbstractDMTopology, ::Misregistration) = false
@inline supports_separable_influence(::GaussianInfluenceWidth, topology::ActuatorGridTopology, mis::Misregistration) =
    supports_separable_topology(topology) && supports_separable_influence(mis)
@inline supports_separable_influence(::GaussianMechanicalCoupling, topology::ActuatorGridTopology, mis::Misregistration) =
    supports_separable_topology(topology) && supports_separable_influence(mis)

@inline supports_dm_misregistration_identification(::AbstractDMInfluenceModel, ::AbstractDMTopology) = false
@inline supports_dm_misregistration_identification(::GaussianInfluenceWidth, ::ActuatorGridTopology) = true
@inline supports_dm_misregistration_identification(::GaussianMechanicalCoupling, ::ActuatorGridTopology) = true

@inline function _clear_separable_influence!(dm::DeformableMirror)
    dm.state.separable_x = nothing
    dm.state.separable_y_t = nothing
    dm.state.separable_tmp = nothing
    return dm
end

@inline function prepare_actuator_commands!(dm::DeformableMirror)
    copyto!(dm.state.actuator_coefs, dm.state.coefs)
    apply_actuator_model!(dm.state.actuator_coefs, actuator_model(dm))
    return dm.state.actuator_coefs
end

@inline apply_actuator_model!(buffer, ::LinearStaticActuators) = buffer

@inline function apply_actuator_model!(buffer, model::ClippedActuators)
    @. buffer = clamp(buffer, model.lo, model.hi)
    return buffer
end

@inline function apply_actuator_model!(buffer, model::ActuatorHealthMap)
    @. buffer = buffer * model.gains
    return buffer
end

@inline function apply_actuator_model!(buffer, model::CompositeDMActuatorModel)
    for stage in model.stages
        apply_actuator_model!(buffer, stage)
    end
    return buffer
end

"""
    build_separable_influence!(dm, tel)

Build the factored Gaussian influence basis used by the fast separable DM
application path.

This optimization is only valid when the topology is a fully active regular
grid and the misregistration preserves the separable Gaussian structure.
"""
function build_separable_influence!(dm::DeformableMirror, tel::Telescope)
    return build_separable_influence!(dm, tel, influence_model(dm), topology(dm))
end

function build_separable_influence!(dm::DeformableMirror, tel::Telescope,
    model::AbstractDMInfluenceModel, topology::AbstractDMTopology)
    supports_separable_influence(model, topology, dm.params.misregistration) || return _clear_separable_influence!(dm)
    n = tel.params.resolution
    n_act = topology_axis_count(topology)
    T = eltype(dm.state.opd)
    sigma2 = T(influence_width(model, topology))^2
    xs = range(T(-1), T(1); length=n_act)
    ys = range(T(-1), T(1); length=n_act)
    cx = T((n + 1) / 2)
    cy = T((n + 1) / 2)
    scale = T(n / 2)
    xbasis_host = Matrix{T}(undef, n, n_act)
    ybasis_host = Matrix{T}(undef, n_act, n)
    @inbounds for (ai, x0) in enumerate(xs)
        x_m, _ = apply_misregistration(dm.params.misregistration, x0, zero(T))
        for i in 1:n
            x = (T(i) - cx) / scale
            xbasis_host[i, ai] = exp(-((x - T(x_m))^2) / (2 * sigma2))
        end
    end
    @inbounds for (aj, y0) in enumerate(ys)
        _, y_m = apply_misregistration(dm.params.misregistration, zero(T), y0)
        for j in 1:n
            y = (T(j) - cy) / scale
            ybasis_host[aj, j] = exp(-((y - T(y_m))^2) / (2 * sigma2))
        end
    end
    xbasis = similar(dm.state.opd, n, n_act)
    ybasis_t = similar(dm.state.opd, n_act, n)
    tmp = similar(dm.state.opd, n, n_act)
    copyto!(xbasis, xbasis_host)
    copyto!(ybasis_t, ybasis_host)
    fill!(tmp, zero(T))
    dm.state.separable_x = xbasis
    dm.state.separable_y_t = ybasis_t
    dm.state.separable_tmp = tmp
    return dm
end

function build_influence_functions!(dm::DeformableMirror, tel::Telescope)
    return build_influence_functions!(execution_style(dm.state.modes), dm, tel, influence_model(dm), topology(dm))
end

"""
    build_influence_functions!(dm, tel)

Materialize the dense actuator-to-OPD operator.

Each column corresponds to one active actuator's influence function sampled on
the telescope pupil grid. This is the reference linear model used for
calibration and for the dense DM application fallback.
"""
function build_influence_functions!(::ScalarCPUStyle, dm::DeformableMirror, tel::Telescope,
    model::GaussianInfluenceWidth, topology::AbstractDMTopology)
    Base.require_one_based_indexing(tel.state.pupil, dm.state.modes)
    n = tel.params.resolution
    sigma = influence_width(model, topology)
    coords = actuator_coordinates(topology)
    cx = (n + 1) / 2
    cy = (n + 1) / 2
    scale = n / 2
    n_commands = topology_command_count(topology)
    @inbounds for idx in 1:n_commands
        x0 = coords[1, idx]
        y0 = coords[2, idx]
        x_m, y_m = apply_misregistration(dm.params.misregistration, x0, y0)
        for i in 1:n, j in 1:n
            x = (i - cx) / scale
            y = (j - cy) / scale
            r2 = (x - x_m)^2 + (y - y_m)^2
            dm.state.modes[(j - 1) * n + i, idx] = exp(-r2 / (2 * sigma^2)) * tel.state.pupil[i, j]
        end
    end
    return dm
end

function build_influence_functions!(style::AcceleratorStyle, dm::DeformableMirror, tel::Telescope,
    model::GaussianInfluenceWidth, topology::AbstractDMTopology)
    Base.require_one_based_indexing(tel.state.pupil, dm.state.modes)
    n = tel.params.resolution
    T = eltype(dm.state.modes)
    sigma2 = T(influence_width(model, topology))^2
    coords = actuator_coordinates(topology)
    cx = T((n + 1) / 2)
    cy = T((n + 1) / 2)
    scale = T(n / 2)
    n_commands = topology_command_count(topology)
    @inbounds for idx in 1:n_commands
        x0 = T(coords[1, idx])
        y0 = T(coords[2, idx])
        x_m, y_m = apply_misregistration(dm.params.misregistration, x0, y0)
        launch_kernel!(style, dm_mode_kernel!, @view(dm.state.modes[:, idx]), tel.state.pupil,
            T(x_m), T(y_m), cx, cy, scale, sigma2, n; ndrange=(n, n))
    end
    synchronize_backend!(style)
    return dm
end

function build_influence_functions!(::ScalarCPUStyle, dm::DeformableMirror, tel::Telescope,
    model::GaussianMechanicalCoupling, topology::AbstractDMTopology)
    return build_influence_functions!(ScalarCPUStyle(), dm, tel,
        GaussianInfluenceWidth(influence_width(model, topology)), topology)
end

function build_influence_functions!(style::AcceleratorStyle, dm::DeformableMirror, tel::Telescope,
    model::GaussianMechanicalCoupling, topology::AbstractDMTopology)
    return build_influence_functions!(style, dm, tel,
        GaussianInfluenceWidth(influence_width(model, topology)), topology)
end

function build_influence_functions!(::ScalarCPUStyle, dm::DeformableMirror, tel::Telescope,
    model::DenseInfluenceMatrix, topology::AbstractDMTopology)
    Base.require_one_based_indexing(dm.state.modes)
    size(model.modes) == size(dm.state.modes) ||
        throw(DimensionMismatchError("DenseInfluenceMatrix size $(size(model.modes)) does not match DM state size $(size(dm.state.modes))"))
    copyto!(dm.state.modes, model.modes)
    return dm
end

function build_influence_functions!(::AcceleratorStyle, dm::DeformableMirror, tel::Telescope,
    model::DenseInfluenceMatrix, topology::AbstractDMTopology)
    size(model.modes) == size(dm.state.modes) ||
        throw(DimensionMismatchError("DenseInfluenceMatrix size $(size(model.modes)) does not match DM state size $(size(dm.state.modes))"))
    copyto!(dm.state.modes, model.modes)
    return dm
end

function build_influence_functions!(::ScalarCPUStyle, dm::DeformableMirror, tel::Telescope,
    model::MeasuredInfluenceFunctions, topology::AbstractDMTopology)
    Base.require_one_based_indexing(dm.state.modes)
    size(model.modes) == size(dm.state.modes) ||
        throw(DimensionMismatchError("MeasuredInfluenceFunctions size $(size(model.modes)) does not match DM state size $(size(dm.state.modes))"))
    copyto!(dm.state.modes, model.modes)
    return dm
end

function build_influence_functions!(::AcceleratorStyle, dm::DeformableMirror, tel::Telescope,
    model::MeasuredInfluenceFunctions, topology::AbstractDMTopology)
    size(model.modes) == size(dm.state.modes) ||
        throw(DimensionMismatchError("MeasuredInfluenceFunctions size $(size(model.modes)) does not match DM state size $(size(dm.state.modes))"))
    copyto!(dm.state.modes, model.modes)
    return dm
end

function apply!(dm::DeformableMirror, tel::Telescope, ::DMAdditive)
    apply_opd!(dm, tel)
    tel.state.opd .+= dm.state.opd
    return tel
end

function apply!(dm::DeformableMirror, tel::Telescope, ::DMReplace)
    apply_opd!(dm, tel)
    tel.state.opd .= dm.state.opd
    return tel
end

@inline function _apply_opd_separable!(::ScalarCPUStyle, dm::DeformableMirror, tel::Telescope)
    xbasis = dm.state.separable_x::typeof(dm.state.opd)
    ybasis_t = dm.state.separable_y_t::typeof(dm.state.opd)
    tmp = dm.state.separable_tmp::typeof(dm.state.opd)
    coefs_grid = dm.state.coefs_grid
    coefs_grid === nothing && throw(InvalidConfiguration("separable DM application requires a grid-backed command buffer"))
    mul!(tmp, xbasis, coefs_grid)
    mul!(dm.state.opd, tmp, ybasis_t)
    dm.state.opd .*= tel.state.pupil
    return dm.state.opd
end

@inline function _apply_opd_separable!(style::AcceleratorStyle, dm::DeformableMirror, tel::Telescope)
    xbasis = dm.state.separable_x::typeof(dm.state.opd)
    ybasis_t = dm.state.separable_y_t::typeof(dm.state.opd)
    tmp = dm.state.separable_tmp::typeof(dm.state.opd)
    coefs_grid = dm.state.coefs_grid
    coefs_grid === nothing && throw(InvalidConfiguration("separable DM application requires a grid-backed command buffer"))
    n_act = topology_axis_count(topology(dm))
    launch_kernel_async!(style, dm_separable_tmp_kernel!,
        tmp, xbasis, coefs_grid, n_act; ndrange=size(tmp))
    launch_kernel!(style, dm_separable_finalize_kernel!,
        dm.state.opd, tmp, ybasis_t, tel.state.pupil, n_act; ndrange=size(dm.state.opd))
    return dm.state.opd
end

@inline function apply_dense!(dm::DeformableMirror, tel::Telescope, ::DMAdditive)
    prepare_actuator_commands!(dm)
    mul!(dm.state.opd_vec, dm.state.modes, dm.state.actuator_coefs)
    tel.state.opd .+= dm.state.opd
    return tel
end

@inline function apply_dense!(dm::DeformableMirror, tel::Telescope, ::DMReplace)
    prepare_actuator_commands!(dm)
    mul!(dm.state.opd_vec, dm.state.modes, dm.state.actuator_coefs)
    tel.state.opd .= dm.state.opd
    return tel
end

@inline function apply_opd!(dm::DeformableMirror, tel::Telescope)
    prepare_actuator_commands!(dm)
    if !isnothing(dm.state.separable_x)
        return apply_opd_separable!(dm, tel)
    end
    mul!(dm.state.opd_vec, dm.state.modes, dm.state.actuator_coefs)
    return dm.state.opd
end

"""
    apply_opd!(dm, tel)

Assemble the DM OPD from the current actuator coefficients without mutating the
telescope OPD.

This chooses the separable `X * C * Y'` path when available and otherwise
falls back to the dense matrix-vector application `modes * actuator_coefs`.
"""
@inline function apply_opd_separable!(dm::DeformableMirror, tel::Telescope)
    return _apply_opd_separable!(execution_style(dm.state.opd), dm, tel)
end
