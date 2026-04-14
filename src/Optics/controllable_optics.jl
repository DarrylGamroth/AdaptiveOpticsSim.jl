using LinearAlgebra

command_storage(::AbstractControllableOptic) =
    throw(InvalidConfiguration("command_storage is not implemented for this controllable optic"))

@inline n_control_dofs(optic::AbstractControllableOptic) = length(command_storage(optic))
@inline command_layout(optic::AbstractControllableOptic) = RuntimeCommandLayout(:optic => n_control_dofs(optic))
@inline controllable_surface_labels(optic::AbstractControllableOptic) = command_segment_labels(command_layout(optic))
@inline supports_segmented_command(optic::AbstractControllableOptic) = length(command_segments(command_layout(optic))) > 1

function _set_command_segments!(setter, total_length::Int, layout::RuntimeCommandLayout, command::NamedTuple;
    require_all::Bool=true)
    layout.total_length == total_length ||
        throw(InvalidConfiguration("command layout total length must match the target command length"))
    labels = Tuple(command_segment_labels(layout))
    provided = Tuple(keys(command))
    if require_all
        provided == labels ||
            throw(InvalidConfiguration("structured command keys $(provided) must match command layout labels $(labels)"))
    else
        all(label -> label in labels, provided) ||
            throw(InvalidConfiguration("structured command keys $(provided) must be a subset of command layout labels $(labels)"))
    end
    @inbounds for seg in command_segments(layout)
        if hasproperty(command, seg.label)
            segment = getproperty(command, seg.label)
            length(segment) == seg.length ||
                throw(DimensionMismatchError("command segment $(seg.label) must have length $(seg.length)"))
            setter(seg, segment)
        elseif require_all
            throw(InvalidConfiguration("structured command is missing segment $(seg.label)"))
        end
    end
    return command
end

@inline function set_command!(optic::AbstractControllableOptic, command::AbstractVector)
    length(command_storage(optic)) == length(command) ||
        throw(DimensionMismatchError("command length must match controllable optic command length"))
    copyto!(command_storage(optic), command)
    return command_storage(optic)
end

@inline function _stage_command!(optic::AbstractControllableOptic, command::AbstractVector{T}, sign::T) where {T<:AbstractFloat}
    apply_command!(command_storage(optic), command, sign)
    return command_storage(optic)
end

@inline function set_command!(optic::AbstractControllableOptic, command::NamedTuple)
    _set_command_segments!(length(command_storage(optic)), command_layout(optic), command) do seg, segment
        copyto!(@view(command_storage(optic)[command_segment_range(seg)]), segment)
    end
    return command_storage(optic)
end

@inline function update_command!(optic::AbstractControllableOptic, command::NamedTuple)
    _set_command_segments!(length(command_storage(optic)), command_layout(optic), command; require_all=false) do seg, segment
        copyto!(@view(command_storage(optic)[command_segment_range(seg)]), segment)
    end
    return command_storage(optic)
end

@inline function _apply_selected!(optic::AbstractControllableOptic, tel::Telescope, mode)
    apply!(optic, tel, mode)
    return tel
end

@inline function _apply_selected!(optic::AbstractControllableOptic, tel::Telescope, mode,
    labels::Tuple{Vararg{Symbol}})
    isempty(labels) || apply!(optic, tel, mode)
    return tel
end

@inline _apply_selected!(optic::AbstractControllableOptic, tel::Telescope, mode,
    label::Symbol) = _apply_selected!(optic, tel, mode, (label,))

function _backend_copy_matrix(host::AbstractMatrix{T}, backend::AbstractArrayBackend, ::Type{T}) where {T<:AbstractFloat}
    backend_array = _resolve_array_backend(backend)
    out = backend_array{T}(undef, size(host)...)
    copyto!(out, host)
    return out
end

function _backend_copy_matrix(host::AbstractMatrix{T}, backend::Type{<:AbstractArray}, ::Type{T}) where {T<:AbstractFloat}
    out = backend{T}(undef, size(host)...)
    copyto!(out, host)
    return out
end

function _modal_command_layout(spec::Union{Nothing,Symbol}, n::Int)
    isnothing(spec) && return RuntimeCommandLayout(:optic => n)
    return RuntimeCommandLayout(spec => n)
end

function _modal_mode_matrix(tel::Telescope, definitions::Tuple, ::Type{T}, backend) where {T<:AbstractFloat}
    selector = require_same_backend(tel, _resolve_backend_selector(backend))
    backend = _resolve_array_backend(selector)
    n = tel.params.resolution
    pupil = Array(tel.state.pupil)
    # Match the existing reference-harness and OOPAO bundle convention:
    # use an n+1 half-open normalized grid and drop the final endpoint.
    xs = collect(range(T(-1), T(1); length=n + 1))[1:n]
    ys = collect(range(T(-1), T(1); length=n + 1))[1:n]
    host = Matrix{T}(undef, n * n, length(definitions))
    for (k, f) in pairs(definitions)
        for j in 1:n, i in 1:n
            host[(j - 1) * n + i, k] = pupil[i, j] ? f(xs[i], ys[j]) : zero(T)
        end
    end
    return _backend_copy_matrix(host, backend, T)
end

function _normalize_pupil_mode!(mode::AbstractMatrix{T}, pupil::AbstractMatrix{Bool}) where {T<:AbstractFloat}
    vals = T[]
    n = size(pupil, 1)
    for j in 1:size(pupil, 2), i in 1:n
        if pupil[i, j]
            push!(vals, mode[(j - 1) * n + i, 1])
        end
    end
    isempty(vals) && return mode
    μ = mean(vals)
    σ = std(vals; corrected=false)
    σ > eps(T) || (σ = one(T))
    @inbounds for idx in axes(mode, 1)
        mode[idx, 1] = (mode[idx, 1] - μ) / σ
    end
    return mode
end

abstract type AbstractModalOpticBasis end

struct FunctionModalBasis{D} <: AbstractModalOpticBasis
    definitions::D
    normalize_modes::Bool
end

FunctionModalBasis(definitions::Tuple; normalize_modes::Bool=false) =
    FunctionModalBasis{typeof(definitions)}(definitions, normalize_modes)

struct MatrixModalBasis{M<:AbstractMatrix} <: AbstractModalOpticBasis
    matrix::M
    normalize_modes::Bool
end

MatrixModalBasis(matrix::AbstractMatrix; normalize_modes::Bool=false) =
    MatrixModalBasis{typeof(matrix)}(matrix, normalize_modes)

struct ZernikeOpticBasis{I,S<:Real} <: AbstractModalOpticBasis
    mode_ids::I
    scale::S
    normalize_modes::Bool
end

function ZernikeOpticBasis(mode_ids::AbstractVector{<:Integer};
    scale::Real=1.0, normalize_modes::Bool=true)
    isempty(mode_ids) && throw(InvalidConfiguration("mode_ids must contain at least one mode"))
    return ZernikeOpticBasis{typeof(mode_ids),typeof(scale)}(mode_ids, scale, normalize_modes)
end

struct CartesianTiltBasis{S<:Real} <: AbstractModalOpticBasis
    scale::S
end

CartesianTiltBasis(; scale::Real=1.0) = CartesianTiltBasis(scale)

struct QuadraticFocusBasis{S<:Real} <: AbstractModalOpticBasis
    scale::S
end

QuadraticFocusBasis(; scale::Real=1.0) = QuadraticFocusBasis(scale)

mutable struct ModalControllableOpticState{T<:AbstractFloat,A<:AbstractMatrix{T},M<:AbstractMatrix{T},V<:AbstractVector{T}}
    opd::A
    opd_vec::V
    modes::M
    coefs::V
end

struct ModalControllableOpticParams{L}
    labels::L
    normalized_modes::Bool
end

struct ModalControllableOptic{P<:ModalControllableOpticParams,S<:ModalControllableOpticState,CL<:RuntimeCommandLayout,B<:AbstractArrayBackend} <: AbstractControllableOptic
    params::P
    state::S
    layout::CL
end

@inline backend(::ModalControllableOptic{<:Any,<:Any,<:Any,B}) where {B} = B()

function _modal_layout(labels::Symbol, n_modes::Int)
    return RuntimeCommandLayout(labels => n_modes)
end

function _modal_layout(labels::Tuple{Vararg{Symbol}}, n_modes::Int)
    length(labels) == n_modes ||
        throw(DimensionMismatchError("mode labels must have length $(n_modes)"))
    return RuntimeCommandLayout((label => 1 for label in labels)...)
end

function _default_zernike_labels(zernike_modes)
    return Tuple(Symbol(:zernike_, j) for j in zernike_modes)
end

@inline _modal_basis_normalize(basis::FunctionModalBasis) = basis.normalize_modes
@inline _modal_basis_normalize(basis::MatrixModalBasis) = basis.normalize_modes
@inline _modal_basis_normalize(basis::ZernikeOpticBasis) = basis.normalize_modes
@inline _modal_basis_normalize(::CartesianTiltBasis) = false
@inline _modal_basis_normalize(::QuadraticFocusBasis) = true
@inline _default_modal_labels(::FunctionModalBasis) = :modal_optic
@inline _default_modal_labels(::MatrixModalBasis) = :modal_optic
@inline _default_modal_labels(basis::ZernikeOpticBasis) = _default_zernike_labels(basis.mode_ids)
@inline _default_modal_labels(::CartesianTiltBasis) = (:x_tilt, :y_tilt)
@inline _default_modal_labels(::QuadraticFocusBasis) = :focus

function _modal_basis_matrix(tel::Telescope, basis::FunctionModalBasis, ::Type{T},
    selector::AbstractArrayBackend) where {T<:AbstractFloat}
    return Array(_modal_mode_matrix(tel, basis.definitions, T, selector))
end

function _modal_basis_matrix(tel::Telescope, basis::MatrixModalBasis, ::Type{T},
    selector::AbstractArrayBackend) where {T<:AbstractFloat}
    n = tel.params.resolution
    size(basis.matrix, 1) == n * n ||
        throw(DimensionMismatchError("mode matrix must have $(n * n) rows for telescope resolution $(n)"))
    return Matrix{T}(basis.matrix)
end

function _modal_basis_matrix(tel::Telescope, basis::ZernikeOpticBasis, ::Type{T},
    selector::AbstractArrayBackend) where {T<:AbstractFloat}
    max_mode = maximum(basis.mode_ids)
    zb = compute_zernike!(ZernikeBasis(tel, max_mode; T=T), tel)
    n = tel.params.resolution
    host_modes = Matrix{T}(undef, n * n, length(basis.mode_ids))
    for (k, j) in pairs(basis.mode_ids)
        copyto!(view(host_modes, :, k), vec(T(basis.scale) .* zb.modes[:, :, j]))
    end
    return host_modes
end

function _modal_basis_matrix(tel::Telescope, basis::CartesianTiltBasis, ::Type{T},
    selector::AbstractArrayBackend) where {T<:AbstractFloat}
    return Array(_modal_mode_matrix(tel, (
        (x, y) -> T(basis.scale) * x,
        (x, y) -> T(basis.scale) * y,
    ), T, selector))
end

function _modal_basis_matrix(tel::Telescope, basis::QuadraticFocusBasis, ::Type{T},
    selector::AbstractArrayBackend) where {T<:AbstractFloat}
    n = tel.params.resolution
    pupil = Array(tel.state.pupil)
    # Use a cell-centered symmetric grid so the quadratic focus mode does not
    # leak tip/tilt under even-resolution sampling.
    coords = collect(T.((2 .* (1:n) .- (n + 1)) ./ n))
    host = Matrix{T}(undef, n * n, 1)
    @inbounds for j in 1:n, i in 1:n
        x = coords[i]
        y = coords[j]
        host[(j - 1) * n + i, 1] = pupil[i, j] ? T(basis.scale) * (x^2 + y^2) : zero(T)
    end
    return host
end

function _materialize_modal_modes(tel::Telescope, host_modes::AbstractMatrix, ::Type{T},
    selector::AbstractArrayBackend; normalize_modes::Bool=false) where {T<:AbstractFloat}
    array_backend = _resolve_array_backend(selector)
    n = tel.params.resolution
    size(host_modes, 1) == n * n ||
        throw(DimensionMismatchError("mode matrix must have $(n * n) rows for telescope resolution $(n)"))
    host = Matrix{T}(undef, size(host_modes)...)
    copyto!(host, host_modes)
    if normalize_modes
        pupil = Array(tel.state.pupil)
        for k in axes(host, 2)
            _normalize_pupil_mode!(view(host, :, k:k), pupil)
        end
    end
    modes = _backend_copy_matrix(host, selector, T)
    opd = array_backend{T}(undef, n, n)
    fill!(opd, zero(T))
    opd_vec = reshape(opd, :)
    coefs = array_backend{T}(undef, size(host, 2))
    fill!(coefs, zero(T))
    return ModalControllableOpticState{T,typeof(opd),typeof(modes),typeof(coefs)}(opd, opd_vec, modes, coefs)
end

function ModalControllableOptic(tel::Telescope, basis::AbstractModalOpticBasis;
    labels::Union{Nothing,Symbol,Tuple{Vararg{Symbol}}}=nothing,
    T::Type{<:AbstractFloat}=Float64,
    backend::AbstractArrayBackend=backend(tel))
    selector = require_same_backend(tel, _resolve_backend_selector(backend))
    host_modes = _modal_basis_matrix(tel, basis, T, selector)
    resolved_labels = isnothing(labels) ? _default_modal_labels(basis) : labels
    normalize_modes = _modal_basis_normalize(basis)
    layout = _modal_layout(resolved_labels, size(host_modes, 2))
    state = _materialize_modal_modes(tel, host_modes, T, selector; normalize_modes=normalize_modes)
    params = ModalControllableOpticParams(resolved_labels, normalize_modes)
    return ModalControllableOptic{typeof(params),typeof(state),typeof(layout),typeof(selector)}(params, state, layout)
end

function ModalControllableOptic(tel::Telescope, mode_definitions::Tuple;
    labels::Union{Symbol,Tuple{Vararg{Symbol}}}=:modal_optic,
    normalize_modes::Bool=false,
    T::Type{<:AbstractFloat}=Float64,
    backend::AbstractArrayBackend=backend(tel))
    return ModalControllableOptic(tel, FunctionModalBasis(mode_definitions; normalize_modes=normalize_modes);
        labels=labels, T=T, backend=backend)
end

function ModalControllableOptic(tel::Telescope, mode_matrix::AbstractMatrix;
    labels::Union{Symbol,Tuple{Vararg{Symbol}}}=:modal_optic,
    normalize_modes::Bool=false,
    T::Type{<:AbstractFloat}=Float64,
    backend::AbstractArrayBackend=backend(tel))
    return ModalControllableOptic(tel, MatrixModalBasis(mode_matrix; normalize_modes=normalize_modes);
        labels=labels, T=T, backend=backend)
end

function ModalControllableOptic(tel::Telescope;
    zernike_modes::AbstractVector{<:Integer},
    labels::Union{Nothing,Symbol,Tuple{Vararg{Symbol}}}=nothing,
    scale::Real=1.0,
    normalize_modes::Bool=true,
    T::Type{<:AbstractFloat}=Float64,
    backend::AbstractArrayBackend=backend(tel))
    return ModalControllableOptic(tel,
        ZernikeOpticBasis(zernike_modes; scale=scale, normalize_modes=normalize_modes);
        labels=labels, T=T, backend=backend)
end

function TipTiltMirror(tel::Telescope; scale::Real=1.0,
    T::Type{<:AbstractFloat}=Float64, backend::AbstractArrayBackend=backend(tel), label::Symbol=:tiptilt)
    return ModalControllableOptic(tel, CartesianTiltBasis(; scale=scale);
        labels=label, T=T, backend=backend)
end

function FocusStage(tel::Telescope; scale::Real=1.0,
    T::Type{<:AbstractFloat}=Float64, backend::AbstractArrayBackend=backend(tel), label::Symbol=:focus)
    return ModalControllableOptic(tel, QuadraticFocusBasis(; scale=scale);
        labels=label, T=T, backend=backend)
end

@inline command_storage(optic::ModalControllableOptic) = optic.state.coefs
@inline command_layout(optic::ModalControllableOptic) = optic.layout

@inline function apply!(optic::ModalControllableOptic, tel::Telescope, ::DMAdditive)
    mul!(optic.state.opd_vec, optic.state.modes, optic.state.coefs)
    tel.state.opd .+= optic.state.opd
    return tel
end

@inline function apply!(optic::ModalControllableOptic, tel::Telescope, ::DMReplace)
    mul!(optic.state.opd_vec, optic.state.modes, optic.state.coefs)
    tel.state.opd .= optic.state.opd
    return tel
end

struct CompositeControllableOptic{O,CV,CL,CR,B<:AbstractArrayBackend} <: AbstractControllableOptic
    optics::O
    command::CV
    layout::CL
    child_ranges::CR
end

@inline backend(::CompositeControllableOptic{<:Any,<:Any,<:Any,<:Any,B}) where {B} = B()

function _composite_segment_specs(label::Symbol, optic::AbstractControllableOptic)
    layout = command_layout(optic)
    segments = command_segments(layout)
    if length(segments) == 1
        return (label => segments[1].length,)
    end
    return Tuple(Symbol(label, :__, seg.label) => seg.length for seg in segments)
end

function CompositeControllableOptic(entries::Vararg{Pair{<:Symbol,<:AbstractControllableOptic},N}) where {N}
    N > 0 || throw(InvalidConfiguration("CompositeControllableOptic requires at least one child optic"))
    optics = ntuple(i -> entries[i].second, N)
    specs = Tuple(vcat(map(entry -> collect(_composite_segment_specs(entry.first, entry.second)), entries)...))
    layout = RuntimeCommandLayout(specs...)
    first_storage = command_storage(optics[1])
    command = similar(first_storage, layout.total_length)
    fill!(command, zero(eltype(command)))
    child_ranges = ntuple(N) do i
        start = i == 1 ? 1 : 1 + sum(length(command_storage(optics[j])) for j in 1:(i - 1))
        start:(start + length(command_storage(optics[i])) - 1)
    end
    selector = require_same_backend(optics...)
    optic = CompositeControllableOptic{typeof(optics),typeof(command),typeof(layout),typeof(child_ranges),typeof(selector)}(optics, command, layout, child_ranges)
    for i in 1:N
        copyto!(@view(optic.command[child_ranges[i]]), command_storage(optics[i]))
    end
    return optic
end

function CompositeControllableOptic(optics::AbstractControllableOptic...)
    length(optics) > 0 || throw(InvalidConfiguration("CompositeControllableOptic requires at least one child optic"))
    labels = map(optic -> begin
        labels = controllable_surface_labels(optic)
        length(labels) == 1 || throw(InvalidConfiguration("CompositeControllableOptic requires explicit labels when a child optic has multiple command segments"))
        labels[1]
    end, optics)
    entries = ntuple(i -> labels[i] => optics[i], length(optics))
    return CompositeControllableOptic(entries...)
end

@inline command_storage(optic::CompositeControllableOptic) = optic.command
@inline command_layout(optic::CompositeControllableOptic) = optic.layout

@inline function set_command!(optic::CompositeControllableOptic, command::AbstractVector)
    length(optic.command) == length(command) ||
        throw(DimensionMismatchError("command length must match composite optic command length"))
    copyto!(optic.command, command)
    @inbounds for i in eachindex(optic.optics)
        set_command!(optic.optics[i], @view(optic.command[optic.child_ranges[i]]))
    end
    return optic.command
end

@inline function _stage_command!(optic::CompositeControllableOptic, command::AbstractVector{T}, sign::T) where {T<:AbstractFloat}
    apply_command!(optic.command, command, sign)
    @inbounds for i in eachindex(optic.optics)
        set_command!(optic.optics[i], @view(optic.command[optic.child_ranges[i]]))
    end
    return optic.command
end

@inline function set_command!(optic::CompositeControllableOptic, command::NamedTuple)
    _set_command_segments!(length(optic.command), optic.layout, command) do seg, segment
        copyto!(@view(optic.command[command_segment_range(seg)]), segment)
    end
    @inbounds for i in eachindex(optic.optics)
        set_command!(optic.optics[i], @view(optic.command[optic.child_ranges[i]]))
    end
    return optic.command
end

@inline function update_command!(optic::CompositeControllableOptic, command::NamedTuple)
    _set_command_segments!(length(optic.command), optic.layout, command; require_all=false) do seg, segment
        copyto!(@view(optic.command[command_segment_range(seg)]), segment)
    end
    @inbounds for i in eachindex(optic.optics)
        set_command!(optic.optics[i], @view(optic.command[optic.child_ranges[i]]))
    end
    return optic.command
end

@inline function apply!(optic::CompositeControllableOptic, tel::Telescope, ::DMAdditive)
    @inbounds for child in optic.optics
        apply!(child, tel, DMAdditive())
    end
    return tel
end

@inline function apply!(optic::CompositeControllableOptic, tel::Telescope, ::DMReplace)
    fill!(tel.state.opd, zero(eltype(tel.state.opd)))
    @inbounds for child in optic.optics
        apply!(child, tel, DMAdditive())
    end
    return tel
end

@inline function _apply_selected!(optic::CompositeControllableOptic, tel::Telescope, ::DMAdditive,
    labels::Tuple{Vararg{Symbol}})
    label_set = Set(labels)
    @inbounds for child in optic.optics
        any(in(label_set), controllable_surface_labels(child)) && apply!(child, tel, DMAdditive())
    end
    return tel
end

@inline function _apply_selected!(optic::CompositeControllableOptic, tel::Telescope, ::DMReplace,
    labels::Tuple{Vararg{Symbol}})
    fill!(tel.state.opd, zero(eltype(tel.state.opd)))
    label_set = Set(labels)
    @inbounds for child in optic.optics
        any(in(label_set), controllable_surface_labels(child)) && apply!(child, tel, DMAdditive())
    end
    return tel
end
