#
# Misregistration identification
#
# This file acquires the physical meta-sensitivity of interaction matrices to
# deformable-mirror misregistration parameters.
#
# The workflow is:
# 1. build a reference interaction matrix at a chosen zero-point misregistration
# 2. form the derivative of the interaction matrix with respect to each
#    selected parameter
# 3. pass the explicit response matrices to
#    AdaptiveOpticsCalibration.Misregistration for numerical estimation
#
const MISREG_FIELDS = (:shift_x, :shift_y, :rotation_deg, :radial_scaling, :tangential_scaling)
const MISREG_FIELD_UNITS = (
    shift_x = "DM actuator-coordinate convention",
    shift_y = "DM actuator-coordinate convention",
    rotation_deg = "degree",
    anamorphosis_angle_deg = "degree",
    radial_scaling = "dimensionless",
    tangential_scaling = "dimensionless",
)

"""
    MetaSensitivity

Store physical interaction-matrix data for complete-data misregistration
estimation.

`D0` is the reference interaction matrix at the physical zero point. `J` has
column-major `vec(D0)` row order and one column for each name in
`field_order`. `epsilon` stores finite-difference step sizes for validation or
fallback acquisition. `field_units` describes the ordered parameter values;
this physical product does not construct an inverse or estimate offsets.
"""
struct MetaSensitivity{T<:AbstractFloat,
    D<:AbstractMatrix{T},
    J<:AbstractMatrix{T},
    O<:Tuple{Vararg{Symbol}},
    U<:Tuple{Vararg{String}}}
    D0::D
    J::J
    epsilon::Misregistration{T}
    field_order::O
    field_units::U
end

function MetaSensitivity(D0::D, J::M, epsilon::Misregistration{T},
    field_order) where {T<:AbstractFloat,D<:AbstractMatrix{T},M<:AbstractMatrix{T}}
    fields = Tuple(field_order)
    size(J, 1) == length(D0) || throw(DimensionMismatchError(
        "meta-sensitivity row count must equal length(D0)",
    ))
    size(J, 2) == length(fields) || throw(DimensionMismatchError(
        "meta-sensitivity column count must equal field_order length",
    ))
    return MetaSensitivity(D0, J, epsilon, fields,
        _misregistration_field_units(fields))
end

function _misregistration_field_units(fields::Tuple{Vararg{Symbol}})
    return map(fields) do field
        hasproperty(MISREG_FIELD_UNITS, field) || throw(InvalidConfiguration(
            "unsupported misregistration field $(field)",
        ))
        return getproperty(MISREG_FIELD_UNITS, field)
    end
end

"""
    compute_meta_sensitivity_matrix(tel, dm, wfs, basis; ...)

Build the sensitivity of the interaction matrix to selected misregistration
parameters.

CPU grid-backed Gaussian DM misregistration uses ForwardDiff by default. Use
`sensitivity=:finite_difference` for validation or accelerator-backed arrays.
The caller-owned `MetaSensitivity` retains physical D₀, J, finite-difference
epsilon, and ordered parameter metadata. This operation performs no numerical
inversion, cache lookup, serialization, or filesystem I/O.
"""
function compute_meta_sensitivity_matrix(tel::Telescope, dm::DeformableMirror, wfs::AbstractWFS,
    basis::AbstractMatrix; misregistration_zero::Misregistration=Misregistration(T=eltype(pupil_reflectivity(tel))),
    epsilon::Misregistration=Misregistration(shift_x=1e-3, shift_y=1e-3, rotation_deg=1e-3, radial_scaling=1e-3,
        tangential_scaling=1e-3, T=eltype(pupil_reflectivity(tel))),
    n_mis_reg::Int=3, field_order=collect(MISREG_FIELDS),
    sensitivity::Symbol=:ad, source=nothing)

    if sensitivity === :ad
        return _compute_meta_sensitivity_matrix_ad(tel, dm, wfs, basis;
            source=source,
            misregistration_zero=misregistration_zero,
            epsilon=epsilon,
            n_mis_reg=n_mis_reg,
            field_order=field_order)
    elseif sensitivity === :finite_difference || sensitivity === :fd
        return _compute_meta_sensitivity_matrix_fd(tel, dm, wfs, basis;
            source=source,
            misregistration_zero=misregistration_zero,
            epsilon=epsilon,
            n_mis_reg=n_mis_reg,
            field_order=field_order)
    end
    throw(InvalidConfiguration("sensitivity must be :ad or :finite_difference"))
end

function _compute_meta_sensitivity_matrix_fd(tel::Telescope, dm::DeformableMirror, wfs::AbstractWFS,
    basis::AbstractMatrix; source=nothing,
    misregistration_zero::Misregistration=Misregistration(T=eltype(pupil_reflectivity(tel))),
    epsilon::Misregistration=Misregistration(shift_x=1e-3, shift_y=1e-3, rotation_deg=1e-3, radial_scaling=1e-3,
        tangential_scaling=1e-3, T=eltype(pupil_reflectivity(tel))),
    n_mis_reg::Int=3, field_order=collect(MISREG_FIELDS),
    amplitude::Real=1e-9)

    T = eltype(pupil_reflectivity(tel))
    pupil = PupilFunction(tel; T=T)
    fields = collect(field_order)[1:min(n_mis_reg, length(field_order))]
    dm_model = influence_model(dm)

    supports_dm_misregistration_identification(dm_model, topology(dm)) ||
        throw(UnsupportedAlgorithm("misregistration identification is only supported for grid-backed Gaussian DeformableMirror models"))

    dm0 = DeformableMirror(tel; topology=topology(dm), influence_model=dm_model,
        misregistration=misregistration_zero, T=T)
    calib0 = _interaction_matrix_for_sensitivity(dm0, wfs, pupil, basis,
        source, amplitude)

    n_elements = length(calib0.matrix)
    meta = zeros(T, n_elements, length(fields))
    for (idx, field) in enumerate(fields)
        eps_val = misregistration_component(epsilon, field)
        if eps_val == 0
            throw(InvalidConfiguration("epsilon for $(field) must be non-zero"))
        end
        base_val = misregistration_component(misregistration_zero, field)
        mis_p = update_misregistration(misregistration_zero, field,
            base_val + eps_val)
        mis_n = update_misregistration(misregistration_zero, field,
            base_val - eps_val)

        dm_p = DeformableMirror(tel; topology=topology(dm),
            influence_model=dm_model, misregistration=mis_p, T=T)
        dm_n = DeformableMirror(tel; topology=topology(dm),
            influence_model=dm_model, misregistration=mis_n, T=T)
        imat_p = _interaction_matrix_for_sensitivity(dm_p, wfs, pupil,
            basis, source, amplitude)
        imat_n = _interaction_matrix_for_sensitivity(dm_n, wfs, pupil,
            basis, source, amplitude)
        meta[:, idx] .= vec((imat_p.matrix .- imat_n.matrix) ./ (2 * eps_val))
    end

    return MetaSensitivity(calib0.matrix, meta, epsilon, fields)
end

"""
    update_misregistration(mis, field, value)

Return a copy of `mis` with one named misregistration field replaced.
"""
function update_misregistration(mis::Misregistration{T}, field::Symbol, value::Real) where {T<:AbstractFloat}
    return Misregistration(; shift_x=field == :shift_x ? T(value) : mis.shift_x,
        shift_y=field == :shift_y ? T(value) : mis.shift_y,
        rotation_deg=field == :rotation_deg ? T(value) : rotation_deg(mis),
        anamorphosis_angle_deg=field == :anamorphosis_angle_deg ? T(value) : anamorphosis_angle_deg(mis),
        tangential_scaling=field == :tangential_scaling ? T(value) : mis.tangential_scaling,
        radial_scaling=field == :radial_scaling ? T(value) : mis.radial_scaling,
        T=T)
end
