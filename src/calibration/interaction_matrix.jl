#
# Interaction-matrix calibration
#
# The interaction matrix records the local linear response of the WFS slopes to
# DM command perturbations.
#
# Each column is built by:
# 1. applying one actuator or command-basis perturbation to the DM
# 2. propagating the modified pupil path through the WFS
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

@inline function _measure_for_calibration!(wfs::AbstractWFS,
    pupil::PupilFunction, src::Union{Nothing,AbstractSource})
    if src === nothing
        return measure!(wfs, pupil)
    end
    return measure!(wfs, pupil, src)
end

@inline function _prepare_interaction_matrix_wfs!(wfs::AbstractWFS,
    pupil::PupilFunction, src::Union{Nothing,AbstractSource})
    _measure_for_calibration!(wfs, pupil, src)
    n_rows = length(slopes(wfs))
    n_rows > 0 || throw(InvalidConfiguration(
        "interaction matrix requires a non-empty calibrated WFS signal"))
    return n_rows
end

@inline function _interaction_matrix_buffer(ref::AbstractArray{T}, n_rows::Int, n_cols::Int) where {T}
    return similar(ref, T, n_rows, n_cols)
end

abstract type AbstractCalibrationCommandPlan end

struct ActuatorCalibrationCommands <: AbstractCalibrationCommandPlan
    count::Int
end

struct MatrixCalibrationCommands{M<:AbstractMatrix} <: AbstractCalibrationCommandPlan
    commands::M
end

@inline calibration_command_count(plan::ActuatorCalibrationCommands) = plan.count
@inline calibration_command_count(plan::MatrixCalibrationCommands) =
    size(plan.commands, 2)

@inline function stage_calibration_command!(coefs::AbstractVector{T},
    ::ActuatorCalibrationCommands, k::Int, amplitude::T) where {T}
    fill!(coefs, zero(T))
    @views coefs[k:k] .= amplitude
    return coefs
end

@inline function stage_calibration_command!(coefs::AbstractVector{T},
    plan::MatrixCalibrationCommands, k::Int, amplitude::T) where {T}
    copyto!(coefs, @view(plan.commands[:, k]))
    coefs .*= amplitude
    return coefs
end

function validate_interaction_matrix_output(out::AbstractMatrix,
    dm::DeformableMirror, wfs::AbstractWFS,
    plan::AbstractCalibrationCommandPlan)
    size(out, 1) == length(slopes(wfs)) ||
        throw(DimensionMismatchError("interaction-matrix output row count must match WFS slopes"))
    size(out, 2) == calibration_command_count(plan) ||
        throw(DimensionMismatchError("interaction-matrix output column count must match calibration commands"))
    eltype(out) == eltype(dm.state.coefs) ||
        throw(InvalidConfiguration("interaction-matrix output element type must match DM coefficients"))
    return out
end

function validate_calibration_commands(dm::DeformableMirror,
    commands::AbstractMatrix)
    size(commands, 1) == length(dm.state.coefs) ||
        throw(DimensionMismatchError("command matrix row count must match DM coefficients"))
    size(commands, 2) > 0 ||
        throw(InvalidConfiguration("interaction matrix requires at least one command mode"))
    return MatrixCalibrationCommands(commands)
end

@inline _interaction_source(::Nothing) = nothing
@inline _interaction_source(src::AbstractSource) = src

function _fill_prepared_interaction_matrix!(out::AbstractMatrix{T},
    dm::DeformableMirror, wfs::AbstractWFS, pupil::PupilFunction,
    plan::AbstractCalibrationCommandPlan, src, amplitude::T) where {T<:AbstractFloat}
    validate_interaction_matrix_output(out, dm, wfs, plan)
    coefs = dm.state.coefs
    opd_base = copy(pupil.opd)
    coefs_base = copy(coefs)
    try
        @inbounds for k in axes(out, 2)
            stage_calibration_command!(coefs, plan, k, amplitude)
            update_surface!(dm)
            apply_surface!(pupil, dm, DMReplace())
            _measure_for_calibration!(wfs, pupil, _interaction_source(src))
            copyto!(@view(out[:, k]), slopes(wfs))
        end
    finally
        copyto!(coefs, coefs_base)
        copyto!(pupil.opd, opd_base)
    end
    return out
end

function _fill_interaction_matrix!(out::AbstractMatrix{T},
    dm::DeformableMirror, wfs::AbstractWFS, pupil::PupilFunction,
    plan::AbstractCalibrationCommandPlan, src, amplitude::T) where {T<:AbstractFloat}
    _prepare_interaction_matrix_wfs!(wfs, pupil, src)
    return _fill_prepared_interaction_matrix!(out, dm, wfs, pupil,
        plan, src, amplitude)
end

"""
    interaction_matrix!(out, dm, wfs, pupil; amplitude=1)
    interaction_matrix!(out, dm, wfs, pupil, src; amplitude=1)
    interaction_matrix!(out, dm, wfs, pupil, commands; amplitude=1)
    interaction_matrix!(out, dm, wfs, pupil, commands, src; amplitude=1)

Fill caller-owned interaction-matrix storage and return an `InteractionMatrix`
view of that storage. `out` may be a normal array, backend-native array, or a
disk-backed `AbstractMatrix` supplied by an integration package.
"""
function interaction_matrix!(out::AbstractMatrix{T}, dm::DeformableMirror,
    wfs::AbstractWFS, pupil::PupilFunction;
    amplitude::Real=1.0) where {T<:AbstractFloat}
    plan = ActuatorCalibrationCommands(length(dm.state.coefs))
    _fill_interaction_matrix!(out, dm, wfs, pupil, plan, nothing,
        T(amplitude))
    return InteractionMatrix(out, T(amplitude))
end

function interaction_matrix!(out::AbstractMatrix{T}, dm::DeformableMirror,
    wfs::AbstractWFS, pupil::PupilFunction, src::AbstractSource;
    amplitude::Real=1.0) where {T<:AbstractFloat}
    plan = ActuatorCalibrationCommands(length(dm.state.coefs))
    _fill_interaction_matrix!(out, dm, wfs, pupil, plan, src, T(amplitude))
    return InteractionMatrix(out, T(amplitude))
end

function interaction_matrix!(out::AbstractMatrix{T}, dm::DeformableMirror,
    wfs::AbstractWFS, pupil::PupilFunction, commands::AbstractMatrix;
    amplitude::Real=1.0) where {T<:AbstractFloat}
    plan = validate_calibration_commands(dm, commands)
    _fill_interaction_matrix!(out, dm, wfs, pupil, plan, nothing,
        T(amplitude))
    return InteractionMatrix(out, T(amplitude))
end

function interaction_matrix!(out::AbstractMatrix{T}, dm::DeformableMirror,
    wfs::AbstractWFS, pupil::PupilFunction, commands::AbstractMatrix,
    src::AbstractSource; amplitude::Real=1.0) where {T<:AbstractFloat}
    plan = validate_calibration_commands(dm, commands)
    _fill_interaction_matrix!(out, dm, wfs, pupil, plan, src, T(amplitude))
    return InteractionMatrix(out, T(amplitude))
end

"""
    interaction_matrix(dm, wfs, pupil; amplitude=1)
    interaction_matrix(dm, wfs, pupil, src; amplitude=1)
    interaction_matrix(dm, wfs, pupil, commands; amplitude=1)
    interaction_matrix(dm, wfs, pupil, commands, src; amplitude=1)

Build the WFS interaction matrix for either actuator-space pushes or an
explicit command basis.

The returned matrix stores one measured slope vector per commanded
perturbation, scaled by the requested calibration amplitude.
"""
function interaction_matrix(dm::DeformableMirror, wfs::AbstractWFS,
    pupil::PupilFunction; amplitude::Real=1.0)
    n_act = length(dm.state.coefs)
    T = eltype(dm.state.coefs)
    n_act > 0 || throw(InvalidConfiguration("interaction matrix requires at least one actuator"))
    plan = ActuatorCalibrationCommands(n_act)
    n_rows = _prepare_interaction_matrix_wfs!(wfs, pupil, nothing)
    out = _interaction_matrix_buffer(pupil.opd, n_rows, n_act)
    _fill_prepared_interaction_matrix!(out, dm, wfs, pupil, plan,
        nothing, T(amplitude))
    return InteractionMatrix(out, T(amplitude))
end

function interaction_matrix(dm::DeformableMirror, wfs::AbstractWFS,
    pupil::PupilFunction,
    src::AbstractSource; amplitude::Real=1.0)
    n_act = length(dm.state.coefs)
    T = eltype(dm.state.coefs)
    n_act > 0 || throw(InvalidConfiguration("interaction matrix requires at least one actuator"))
    plan = ActuatorCalibrationCommands(n_act)
    n_rows = _prepare_interaction_matrix_wfs!(wfs, pupil, src)
    out = _interaction_matrix_buffer(pupil.opd, n_rows, n_act)
    _fill_prepared_interaction_matrix!(out, dm, wfs, pupil, plan,
        src, T(amplitude))
    return InteractionMatrix(out, T(amplitude))
end

function interaction_matrix(dm::DeformableMirror, wfs::AbstractWFS,
    pupil::PupilFunction,
    commands::AbstractMatrix; amplitude::Real=1.0)
    plan = validate_calibration_commands(dm, commands)
    T = eltype(dm.state.coefs)
    n_rows = _prepare_interaction_matrix_wfs!(wfs, pupil, nothing)
    out = _interaction_matrix_buffer(pupil.opd, n_rows,
        calibration_command_count(plan))
    _fill_prepared_interaction_matrix!(out, dm, wfs, pupil, plan,
        nothing, T(amplitude))
    return InteractionMatrix(out, T(amplitude))
end

function interaction_matrix(dm::DeformableMirror, wfs::AbstractWFS,
    pupil::PupilFunction,
    commands::AbstractMatrix, src::AbstractSource; amplitude::Real=1.0)
    plan = validate_calibration_commands(dm, commands)
    T = eltype(dm.state.coefs)
    n_rows = _prepare_interaction_matrix_wfs!(wfs, pupil, src)
    out = _interaction_matrix_buffer(pupil.opd, n_rows,
        calibration_command_count(plan))
    _fill_prepared_interaction_matrix!(out, dm, wfs, pupil, plan,
        src, T(amplitude))
    return InteractionMatrix(out, T(amplitude))
end
