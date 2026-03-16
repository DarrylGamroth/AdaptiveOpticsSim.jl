using LinearAlgebra

struct InteractionMatrix{T<:AbstractFloat,M<:AbstractMatrix{T}}
    matrix::M
    amplitude::T
end

@inline function _measure_for_calibration!(wfs::AbstractWFS, tel::Telescope, src::Union{Nothing,AbstractSource})
    if src === nothing
        return measure!(wfs, tel)
    end
    return measure!(wfs, tel, src)
end

@inline function _interaction_matrix_buffer(ref::AbstractArray{T}, n_rows::Int, n_cols::Int) where {T}
    return similar(ref, T, n_rows, n_cols)
end

function interaction_matrix(dm::DeformableMirror, wfs::AbstractWFS, tel::Telescope; amplitude::Real=1.0)
    n_act = length(dm.state.coefs)
    T = eltype(dm.state.coefs)
    mat = nothing

    opd_base = copy(tel.state.opd)
    coefs = dm.state.coefs

    for k in 1:n_act
        fill!(coefs, zero(T))
        @views coefs[k:k] .= T(amplitude)
        apply!(dm, tel, DMReplace())
        _measure_for_calibration!(wfs, tel, nothing)
        if mat === nothing
            mat = _interaction_matrix_buffer(tel.state.opd, length(wfs.state.slopes), n_act)
        end
        copyto!(@view(mat[:, k]), wfs.state.slopes)
    end

    tel.state.opd .= opd_base
    mat === nothing && throw(InvalidConfiguration("interaction matrix requires at least one actuator"))
    return InteractionMatrix{T, typeof(mat)}(mat, T(amplitude))
end

function interaction_matrix(dm::DeformableMirror, wfs::AbstractWFS, tel::Telescope,
    src::AbstractSource; amplitude::Real=1.0)
    n_act = length(dm.state.coefs)
    T = eltype(dm.state.coefs)
    mat = nothing

    opd_base = copy(tel.state.opd)
    coefs = dm.state.coefs

    for k in 1:n_act
        fill!(coefs, zero(T))
        @views coefs[k:k] .= T(amplitude)
        apply!(dm, tel, DMReplace())
        _measure_for_calibration!(wfs, tel, src)
        if mat === nothing
            mat = _interaction_matrix_buffer(tel.state.opd, length(wfs.state.slopes), n_act)
        end
        copyto!(@view(mat[:, k]), wfs.state.slopes)
    end

    tel.state.opd .= opd_base
    mat === nothing && throw(InvalidConfiguration("interaction matrix requires at least one actuator"))
    return InteractionMatrix{T, typeof(mat)}(mat, T(amplitude))
end

function interaction_matrix(dm::DeformableMirror, wfs::AbstractWFS, tel::Telescope,
    commands::AbstractMatrix; amplitude::Real=1.0)
    n_act = length(dm.state.coefs)
    if size(commands, 1) != n_act
        throw(DimensionMismatchError("command matrix row count must match DM coefficients"))
    end
    n_modes = size(commands, 2)
    T = eltype(dm.state.coefs)
    mat = nothing

    opd_base = copy(tel.state.opd)
    coefs = dm.state.coefs

    for k in 1:n_modes
        @views coefs .= T(amplitude) .* commands[:, k]
        apply!(dm, tel, DMReplace())
        _measure_for_calibration!(wfs, tel, nothing)
        if mat === nothing
            mat = _interaction_matrix_buffer(tel.state.opd, length(wfs.state.slopes), n_modes)
        end
        copyto!(@view(mat[:, k]), wfs.state.slopes)
    end

    tel.state.opd .= opd_base
    mat === nothing && throw(InvalidConfiguration("interaction matrix requires at least one command mode"))
    return InteractionMatrix{T, typeof(mat)}(mat, T(amplitude))
end

function interaction_matrix(dm::DeformableMirror, wfs::AbstractWFS, tel::Telescope,
    commands::AbstractMatrix, src::AbstractSource; amplitude::Real=1.0)
    n_act = length(dm.state.coefs)
    if size(commands, 1) != n_act
        throw(DimensionMismatchError("command matrix row count must match DM coefficients"))
    end
    n_modes = size(commands, 2)
    T = eltype(dm.state.coefs)
    mat = nothing

    opd_base = copy(tel.state.opd)
    coefs = dm.state.coefs

    for k in 1:n_modes
        @views coefs .= T(amplitude) .* commands[:, k]
        apply!(dm, tel, DMReplace())
        _measure_for_calibration!(wfs, tel, src)
        if mat === nothing
            mat = _interaction_matrix_buffer(tel.state.opd, length(wfs.state.slopes), n_modes)
        end
        copyto!(@view(mat[:, k]), wfs.state.slopes)
    end

    tel.state.opd .= opd_base
    mat === nothing && throw(InvalidConfiguration("interaction matrix requires at least one command mode"))
    return InteractionMatrix{T, typeof(mat)}(mat, T(amplitude))
end
