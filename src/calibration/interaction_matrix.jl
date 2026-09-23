#
# Interaction-matrix calibration
#
# The interaction matrix records the local linear response of a declared WFS
# measurement to DM command perturbations.
#
# Each column is built by:
# 1. applying one actuator or command-basis perturbation to the DM
# 2. invoking the composed optical, detector, and estimator callback
# 3. recording the declared measurement vector
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

@inline function _calibration_measurement_storage(measurement::WFSMeasurement,
    pupil::PupilFunction)
    WavefrontSensors.validate_wfs_measurement(measurement)
    values = measurement_storage(measurement)
    values isa AbstractVector || throw(InvalidConfiguration(
        "interaction matrix requires an explicitly ordered vector WFS measurement"))
    typeof(backend(values)) === typeof(backend(pupil.opd)) &&
        compute_device(values) == compute_device(pupil.opd) ||
        throw(InvalidConfiguration(
            "WFS measurement and calibration pupil must share a backend and device"))
    return values
end

@inline function _interaction_measurement_length(measurement::WFSMeasurement,
    pupil::PupilFunction)
    values = _calibration_measurement_storage(measurement, pupil)
    n_rows = length(values)
    n_rows > 0 || throw(InvalidConfiguration(
        "interaction matrix requires a non-empty WFS measurement"))
    return n_rows
end

@inline function _interaction_matrix_buffer(ref::AbstractArray{T}, n_rows::Int, n_cols::Int) where {T}
    return similar(ref, T, n_rows, n_cols)
end

"""
    AbstractCalibrationCommandPlan

Internal run-immutable interface for interaction-matrix command sequences.
Implementations define `calibration_command_count` and
`stage_calibration_command!`. Staging overwrites and returns caller-owned
coefficient storage; it must not retain that destination. A plan may be reused
concurrently only when its command data are not mutated and every invocation
has distinct DM, measurement, pupil, and destination owners. Validation rejects shape,
command-count, and output numeric-type incompatibilities before calibration
mutation. Command storage must not alias the destination coefficients, and the
two storage backends must support the staged copy and scaling operations.
"""
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
@inline calibration_command_storage(::ActuatorCalibrationCommands) = nothing
@inline calibration_command_storage(plan::MatrixCalibrationCommands) =
    plan.commands

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
    dm::DeformableMirror, measurement::WFSMeasurement,
    pupil::PupilFunction, plan::AbstractCalibrationCommandPlan)
    values = measurement_storage(measurement)
    size(out, 1) == length(values) ||
        throw(DimensionMismatchError("interaction-matrix output row count must match WFS measurement"))
    size(out, 2) == calibration_command_count(plan) ||
        throw(DimensionMismatchError("interaction-matrix output column count must match calibration commands"))
    eltype(out) == eltype(dm.state.coefs) ||
        throw(InvalidConfiguration("interaction-matrix output element type must match DM coefficients"))
    typeof(backend(dm.state.coefs)) === typeof(backend(values)) &&
        compute_device(dm.state.coefs) == compute_device(values) ||
        throw(InvalidConfiguration(
            "DM coefficients and WFS measurement must share a backend and device"))
    Base.mightalias(out, values) && throw(InvalidConfiguration(
        "interaction-matrix output must not alias WFS measurement storage"))
    commands = calibration_command_storage(plan)
    if commands !== nothing
        (Base.mightalias(out, commands) ||
         Base.mightalias(values, commands)) && throw(InvalidConfiguration(
            "interaction-matrix output or WFS measurement must not alias calibration commands"))
    end
    for storage in (dm.state.coefs, dm.state.actuator_coefs,
                    dm.state.opd, pupil.opd)
        (Base.mightalias(out, storage) ||
         Base.mightalias(values, storage)) && throw(InvalidConfiguration(
            "interaction-matrix output or WFS measurement must not alias DM or pupil storage"))
    end
    return out
end

function validate_calibration_commands(dm::DeformableMirror,
    commands::AbstractMatrix)
    size(commands, 1) == length(dm.state.coefs) ||
        throw(DimensionMismatchError("command matrix row count must match DM coefficients"))
    size(commands, 2) > 0 ||
        throw(InvalidConfiguration("interaction matrix requires at least one command mode"))
    Base.mightalias(commands, dm.state.coefs) && throw(
        InvalidConfiguration(
            "calibration commands must not alias DM coefficient storage"))
    return MatrixCalibrationCommands(commands)
end

function _fill_prepared_interaction_matrix!(out::AbstractMatrix{T},
    dm::DeformableMirror, measurement::WFSMeasurement,
    pupil::PupilFunction, plan::AbstractCalibrationCommandPlan,
    measure_callback, amplitude::T) where {T<:AbstractFloat}
    values = _calibration_measurement_storage(measurement, pupil)
    validate_interaction_matrix_output(out, dm, measurement, pupil, plan)
    coefs = dm.state.coefs
    opd_base = copy(pupil.opd)
    coefs_base = copy(coefs)
    surface_base = copy(dm.state.opd)
    actuator_coefs_base = copy(dm.state.actuator_coefs)
    try
        @inbounds for k in axes(out, 2)
            stage_calibration_command!(coefs, plan, k, amplitude)
            update_surface!(dm)
            apply_surface!(pupil, dm, DMReplace())
            measure_callback(measurement, pupil)
            WavefrontSensors.validate_wfs_measurement(measurement)
            copyto!(@view(out[:, k]), values)
        end
    finally
        copyto!(coefs, coefs_base)
        copyto!(dm.state.actuator_coefs, actuator_coefs_base)
        copyto!(dm.state.opd, surface_base)
        copyto!(pupil.opd, opd_base)
    end
    return out
end

function _fill_interaction_matrix!(out::AbstractMatrix{T},
    dm::DeformableMirror, measurement::WFSMeasurement,
    pupil::PupilFunction, plan::AbstractCalibrationCommandPlan,
    measure_callback, amplitude::T) where {T<:AbstractFloat}
    _interaction_measurement_length(measurement, pupil)
    return _fill_prepared_interaction_matrix!(out, dm, measurement, pupil,
        plan, measure_callback, amplitude)
end

"""
    interaction_matrix!(out, dm, measurement, pupil, measure_callback; amplitude=1)
    interaction_matrix!(out, dm, measurement, pupil, commands, measure_callback; amplitude=1)

Fill caller-owned interaction-matrix storage and return an `InteractionMatrix`
view of that storage. The `WFSMeasurement` storage must be an explicitly
ordered vector; a caller using an image or scalar product must define its
vector ordering and units before calibration. `measure_callback(measurement,
pupil)` writes that exact caller-owned vector; its callable owns any source,
detector, or estimator context. `out` may be a normal array, backend-native
array, or a disk-backed `AbstractMatrix` supplied by an integration package.
Output, measurement, commands, DM storage, and pupil OPD must not alias.
"""
function interaction_matrix!(out::AbstractMatrix{T}, dm::DeformableMirror,
    measurement::WFSMeasurement, pupil::PupilFunction, measure_callback;
    amplitude::Real=1.0) where {T<:AbstractFloat}
    plan = ActuatorCalibrationCommands(length(dm.state.coefs))
    _fill_interaction_matrix!(out, dm, measurement, pupil, plan, measure_callback,
        T(amplitude))
    return InteractionMatrix(out, T(amplitude))
end

function interaction_matrix!(out::AbstractMatrix{T}, dm::DeformableMirror,
    measurement::WFSMeasurement, pupil::PupilFunction,
    commands::AbstractMatrix, measure_callback;
    amplitude::Real=1.0) where {T<:AbstractFloat}
    plan = validate_calibration_commands(dm, commands)
    _fill_interaction_matrix!(out, dm, measurement, pupil, plan,
        measure_callback, T(amplitude))
    return InteractionMatrix(out, T(amplitude))
end

"""
    interaction_matrix(dm, measurement, pupil, measure_callback; amplitude=1)
    interaction_matrix(dm, measurement, pupil, commands, measure_callback; amplitude=1)

Build the declared WFS measurement's interaction matrix for either
actuator-space pushes or an explicit command basis.

The returned matrix stores one measured signal vector per commanded
perturbation, scaled by the requested calibration amplitude.
"""
function interaction_matrix(dm::DeformableMirror,
    measurement::WFSMeasurement, pupil::PupilFunction, measure_callback;
    amplitude::Real=1.0)
    n_act = length(dm.state.coefs)
    T = eltype(dm.state.coefs)
    n_act > 0 || throw(InvalidConfiguration("interaction matrix requires at least one actuator"))
    plan = ActuatorCalibrationCommands(n_act)
    n_rows = _interaction_measurement_length(measurement, pupil)
    out = _interaction_matrix_buffer(pupil.opd, n_rows, n_act)
    _fill_prepared_interaction_matrix!(out, dm, measurement, pupil, plan,
        measure_callback, T(amplitude))
    return InteractionMatrix(out, T(amplitude))
end

function interaction_matrix(dm::DeformableMirror,
    measurement::WFSMeasurement, pupil::PupilFunction,
    commands::AbstractMatrix, measure_callback; amplitude::Real=1.0)
    plan = validate_calibration_commands(dm, commands)
    T = eltype(dm.state.coefs)
    n_rows = _interaction_measurement_length(measurement, pupil)
    out = _interaction_matrix_buffer(pupil.opd, n_rows,
        calibration_command_count(plan))
    _fill_prepared_interaction_matrix!(out, dm, measurement, pupil, plan,
        measure_callback, T(amplitude))
    return InteractionMatrix(out, T(amplitude))
end
