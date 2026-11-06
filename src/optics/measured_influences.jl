function _measured_influence_integer(value::Integer, name)
    1 <= value <= typemax(Int) ||
        throw(InvalidConfiguration("$name must be a bounded positive integer"))
    return Int(value)
end
_measured_influence_integer(value, name) =
    throw(InvalidConfiguration("$name must be an integer"))
_measured_influence_integer(value::Bool, name) =
    throw(InvalidConfiguration("$name must not be Bool"))

function _measured_influence_number(value::Real, name)
    isfinite(value) || throw(InvalidConfiguration("$name must be finite"))
    number = Float64(value)
    isfinite(number) ||
        throw(InvalidConfiguration("$name must be representable as finite Float64"))
    return number
end
_measured_influence_number(value, name) =
    throw(InvalidConfiguration("$name must be real"))
_measured_influence_number(value::Bool, name) =
    throw(InvalidConfiguration("$name must not be Bool"))

function _measured_influence_tuple(value::Union{Tuple,AbstractVector}, count, convert, name)
    length(value) == count ||
        throw(DimensionMismatchError("$name must have $count elements"))
    return Tuple(convert(element, name) for element in value)
end
_measured_influence_tuple(value, count, convert, name) =
    throw(InvalidConfiguration("$name must be a tuple or vector"))

"""
    resample_measured_influences(samples; grid_shape, sample_origin, sample_pitch,
        pupil_resolution, actuator_pitch_pupil_radius, command_indices,
        opd_scale_m, reference_command_m)

Resample decoded measured influence functions onto a square pupil. `samples`
is a one-based, host-accessible Float32 or Float64 cube with logical axes
`(source column, source row, influence plane)`; views are supported. `grid_shape`
is the mirror's `(columns, rows)`. `sample_origin` and positive `sample_pitch`
locate samples in that source grid. `actuator_pitch_pupil_radius` is the positive
actuator pitch divided by pupil radius. Only identity pupil registration is
represented by this separable operation.

The pupil's first array axis increases in x and its second increases in y;
source columns increase rightward and source rows downward. Bilinear sampling
must cover the full pupil square, without extrapolation, aperture masking,
normalization, or an added reflection factor. `command_indices` maps each plane
to a unique physical-command column and must permute `1:size(samples, 3)`.

Multiply sampled values by `opd_scale_m / reference_command_m`, where
`opd_scale_m` is positive metres per stored sample unit and the finite nonzero signed
`reference_command_m` is the measured command in metres. The ratio's sign is
preserved. Return a dense `Matrix{Float32}` of dimensionless OPD response per
command, with column-major pupil samples in each column. Nonfinite Float32
results are rejected. This cold host operation performs no file decoding or
artifact admission.
"""
function resample_measured_influences(
    samples::AbstractArray{<:Union{Float32,Float64},3};
    grid_shape, sample_origin, sample_pitch, pupil_resolution,
    actuator_pitch_pupil_radius, command_indices, opd_scale_m, reference_command_m,
)
    all(axis -> axis == Base.OneTo(length(axis)), axes(samples)) ||
        throw(InvalidConfiguration("measured influence samples must use one-based axes"))
    size(samples, 1) >= 2 && size(samples, 2) >= 2 && size(samples, 3) > 0 ||
        throw(DimensionMismatchError("measured influences require two spatial samples per axis and at least one plane"))
    columns, rows = _measured_influence_tuple(
        grid_shape, 2, _measured_influence_integer, "grid_shape")
    origin = _measured_influence_tuple(
        sample_origin, 2, _measured_influence_number, "sample_origin")
    pitch = _measured_influence_tuple(
        sample_pitch, 2, _measured_influence_number, "sample_pitch")
    all(>(0), pitch) || throw(InvalidConfiguration("sample_pitch must be positive"))
    resolution = _measured_influence_integer(pupil_resolution, "pupil_resolution")
    actuator_pitch = _measured_influence_number(
        actuator_pitch_pupil_radius, "actuator_pitch_pupil_radius")
    actuator_pitch > 0 ||
        throw(InvalidConfiguration("actuator_pitch_pupil_radius must be positive"))
    count = size(samples, 3)
    commands = _measured_influence_tuple(
        command_indices, count, _measured_influence_integer, "command_indices")
    Set(commands) == Set(1:count) ||
        throw(InvalidConfiguration("command_indices must permute the physical commands"))
    opd_scale = _measured_influence_number(opd_scale_m, "opd_scale_m")
    opd_scale > 0 || throw(InvalidConfiguration("opd_scale_m must be positive"))
    reference = _measured_influence_number(reference_command_m, "reference_command_m")
    !iszero(reference) || throw(InvalidConfiguration("reference_command_m must be nonzero"))
    scale = opd_scale / reference
    isfinite(scale) && !iszero(scale) ||
        throw(InvalidConfiguration("influence OPD conversion must be finite and nonzero"))
    resolution <= typemax(Int) ÷ resolution ÷ count ÷ sizeof(Float32) ||
        throw(InvalidConfiguration("measured influence output exceeds machine array bounds"))
    half = resolution / 2
    coordinates = ntuple(2) do axis
        centre = (axis == 1 ? columns : rows) / 2 + 0.5
        direction = axis == 1 ? 1 : -1
        values = [
            (
                (
                    centre +
                    direction * ((i - (resolution + 1) / 2) / half) /
                    actuator_pitch - origin[axis]
                ) / pitch[axis] + 1
            ) for i in 1:resolution
        ]
        all(isfinite, values) && minimum(values) >= 1 &&
            maximum(values) <= size(samples, axis) || throw(InvalidConfiguration(
                "influence sampling domain does not cover the full target pupil square"))
        values
    end
    column_samples, row_samples = coordinates
    modes = Matrix{Float32}(undef, resolution * resolution, count)
    for (plane, command) in pairs(commands)
        interpolant = Interpolations.interpolate(
            view(samples, :, :, plane), Interpolations.BSpline(Interpolations.Linear()))
        for j in 1:resolution, i in 1:resolution
            value = Float32(interpolant(column_samples[i], row_samples[j]) * scale)
            isfinite(value) || throw(InvalidConfiguration(
                "resampled influence response must remain finite Float32"))
            modes[i+(j-1)*resolution, command] = value
        end
    end
    return modes
end

resample_measured_influences(samples; kwargs...) = throw(InvalidConfiguration(
    "measured influence samples must be a Float32 or Float64 three-dimensional array"))
