#
# Interaction-matrix calibration
#
# The interaction matrix records the local linear response of the WFS slopes to
# DM command perturbations.
#
# Each column is built by:
# 1. applying one actuator or command-basis perturbation to the DM
# 2. propagating the modified telescope phase through the WFS
# 3. recording the measured slope vector
#
# This provides the linear operator used by control matrices, modal
# reconstructors, and interaction-matrix tomography.
#
struct InteractionMatrix{T<:AbstractFloat,M<:AbstractMatrix{T}}
    matrix::M
    amplitude::T
end

@inline forward_operator(imat::InteractionMatrix) = imat.matrix
@inline calibration_amplitude(imat::InteractionMatrix) = imat.amplitude

@inline function _measure_for_calibration!(wfs::AbstractWFS, tel::Telescope, src::Union{Nothing,AbstractSource})
    if src === nothing
        return measure!(wfs, tel)
    end
    return measure!(wfs, tel, src)
end

@inline function _interaction_matrix_buffer(ref::AbstractArray{T}, n_rows::Int, n_cols::Int) where {T}
    return similar(ref, T, n_rows, n_cols)
end

"""
    interaction_matrix(dm, wfs, tel; amplitude=1)
    interaction_matrix(dm, wfs, tel, src; amplitude=1)
    interaction_matrix(dm, wfs, tel, commands; amplitude=1)
    interaction_matrix(dm, wfs, tel, commands, src; amplitude=1)

Build the WFS interaction matrix for either actuator-space pushes or an
explicit command basis.

The returned matrix stores one measured slope vector per commanded
perturbation, scaled by the requested calibration amplitude.
"""
function interaction_matrix(dm::DeformableMirror, wfs::AbstractWFS, tel::Telescope; amplitude::Real=1.0)
    n_act = length(dm.state.coefs)
    T = eltype(dm.state.coefs)
    mat = nothing

    opd_base = copy(tel.state.opd)
    coefs = dm.state.coefs

    for k in 1:n_act
        fill!(coefs, zero(T))
        @views coefs[k:k] .= T(amplitude)
        apply_dense!(dm, tel, DMReplace())
        _measure_for_calibration!(wfs, tel, nothing)
        if mat === nothing
            mat = _interaction_matrix_buffer(tel.state.opd, length(slopes(wfs)), n_act)
        end
        copyto!(@view(mat[:, k]), slopes(wfs))
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
        apply_dense!(dm, tel, DMReplace())
        _measure_for_calibration!(wfs, tel, src)
        if mat === nothing
            mat = _interaction_matrix_buffer(tel.state.opd, length(slopes(wfs)), n_act)
        end
        copyto!(@view(mat[:, k]), slopes(wfs))
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
        apply_dense!(dm, tel, DMReplace())
        _measure_for_calibration!(wfs, tel, nothing)
        if mat === nothing
            mat = _interaction_matrix_buffer(tel.state.opd, length(slopes(wfs)), n_modes)
        end
        copyto!(@view(mat[:, k]), slopes(wfs))
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
        apply_dense!(dm, tel, DMReplace())
        _measure_for_calibration!(wfs, tel, src)
        if mat === nothing
            mat = _interaction_matrix_buffer(tel.state.opd, length(slopes(wfs)), n_modes)
        end
        copyto!(@view(mat[:, k]), slopes(wfs))
    end

    tel.state.opd .= opd_base
    mat === nothing && throw(InvalidConfiguration("interaction matrix requires at least one command mode"))
    return InteractionMatrix{T, typeof(mat)}(mat, T(amplitude))
end
