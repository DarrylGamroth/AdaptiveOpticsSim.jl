command_storage(::AbstractControllableOptic) =
    throw(InvalidConfiguration("command_storage is not implemented for this controllable optic"))

@inline n_control_dofs(optic::AbstractControllableOptic) = length(command_storage(optic))
@inline surface_opd(optic::AbstractControllableOptic) = optic.state.opd

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

@inline function _apply_selected!(optic::AbstractControllableOptic,
    pupil::PupilFunction, mode)
    update_surface!(optic)
    apply_surface!(pupil, optic, mode)
    return pupil
end

@inline function _apply_selected!(optic::AbstractControllableOptic,
    pupil::PupilFunction, mode,
    labels::Tuple{Vararg{Symbol}})
    isempty(labels) || _apply_selected!(optic, pupil, mode)
    return pupil
end

@inline _apply_selected!(optic::AbstractControllableOptic,
    pupil::PupilFunction, mode, label::Symbol) =
    _apply_selected!(optic, pupil, mode, (label,))

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

function _modal_mode_matrix(tel::Telescope, definitions::Tuple, ::Type{T}, backend) where {T<:AbstractFloat}
    selector = require_same_backend(tel, _resolve_backend_selector(backend))
    backend = _resolve_array_backend(selector)
    n = tel.params.resolution
    pupil = Array(pupil_mask(tel))
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

struct ModalControllableOpticState{T<:AbstractFloat,A<:AbstractMatrix{T},M<:AbstractMatrix{T},V<:AbstractVector{T}}
    opd::A
    opd_vec::V
    modes::M
    coefs::V
end

struct ModalControllableOpticParams{L}
    labels::L
    normalized_modes::Bool
end

struct ModalControllableOptic{
    P<:ModalControllableOpticParams,
    S<:ModalControllableOpticState,
    B<:AbstractArrayBackend,
} <: AbstractControllableOptic
    params::P
    state::S
end

@inline backend(::ModalControllableOptic{<:Any,<:Any,B}) where {B} = B()

@inline _validate_modal_labels(labels::Symbol, ::Int) = labels

function _validate_modal_labels(labels::Tuple{Vararg{Symbol}}, n_modes::Int)
    length(labels) == n_modes ||
        throw(DimensionMismatchError("mode labels must have length $(n_modes)"))
    return labels
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
    n = tel.params.resolution
    pupil = Array(pupil_mask(tel))
    cx = (n + 1) / 2
    cy = (n + 1) / 2
    scale = n / 2
    host_modes = Matrix{T}(undef, n * n, length(basis.mode_ids))
    @inbounds for (k, j) in pairs(basis.mode_ids)
        n_mode, m_mode = noll_to_nm(j)
        for col in 1:n, row in 1:n
            idx = (col - 1) * n + row
            if pupil[row, col]
                x = (row - cx) / scale
                y = (col - cy) / scale
                r = sqrt(x^2 + y^2)
                theta = atan(y, x)
                radial = zernike_radial(n_mode, m_mode, r)
                value = if m_mode >= 0
                    radial * cos(m_mode * theta)
                else
                    radial * sin(abs(m_mode) * theta)
                end
                host_modes[idx, k] = T(basis.scale) * T(value)
            else
                host_modes[idx, k] = zero(T)
            end
        end
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
    pupil = Array(pupil_mask(tel))
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
        pupil = Array(pupil_mask(tel))
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
    resolved_labels = _validate_modal_labels(resolved_labels,
        size(host_modes, 2))
    state = _materialize_modal_modes(tel, host_modes, T, selector; normalize_modes=normalize_modes)
    params = ModalControllableOpticParams(resolved_labels, normalize_modes)
    return ModalControllableOptic{
        typeof(params),typeof(state),typeof(selector)}(params, state)
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

@kernel function modal_opd_one_kernel!(opd, modes, coefs, n::Int)
    idx = @index(Global, Linear)
    if idx <= n
        @inbounds opd[idx] = modes[idx, 1] * coefs[1]
    end
end

@kernel function modal_opd_two_kernel!(opd, modes, coefs, n::Int)
    idx = @index(Global, Linear)
    if idx <= n
        @inbounds opd[idx] = modes[idx, 1] * coefs[1] + modes[idx, 2] * coefs[2]
    end
end

@inline function _apply_modal_opd!(::ExecutionStyle, optic::ModalControllableOptic)
    mul!(optic.state.opd_vec, optic.state.modes, optic.state.coefs)
    return optic.state.opd
end

@inline function _apply_modal_opd!(style::AcceleratorStyle, optic::ModalControllableOptic)
    n_modes = size(optic.state.modes, 2)
    n = length(optic.state.opd_vec)
    if n_modes == 1
        launch_kernel!(style, modal_opd_one_kernel!, optic.state.opd_vec,
            optic.state.modes, optic.state.coefs, n; ndrange=n)
    elseif n_modes == 2
        launch_kernel!(style, modal_opd_two_kernel!, optic.state.opd_vec,
            optic.state.modes, optic.state.coefs, n; ndrange=n)
    else
        mul!(optic.state.opd_vec, optic.state.modes, optic.state.coefs)
    end
    return optic.state.opd
end

@inline function _apply_modal_opd!(optic::ModalControllableOptic)
    return _apply_modal_opd!(execution_style(optic.state.opd_vec), optic)
end

function update_surface!(optic::ModalControllableOptic)
    _apply_modal_opd!(optic)
    return optic
end
