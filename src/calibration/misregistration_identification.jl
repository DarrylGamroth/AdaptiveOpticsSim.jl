#
# Misregistration identification and SPRINT
#
# This file implements a finite-difference meta-sensitivity approach for
# estimating DM or WFS misregistration parameters from interaction matrices.
#
# The workflow is:
# 1. build a reference interaction matrix at a chosen zero-point misregistration
# 2. perturb selected misregistration parameters one at a time by `±epsilon`
# 3. form the derivative of the interaction matrix with respect to each
#    parameter
# 4. invert that meta-sensitivity matrix
# 5. estimate parameter offsets from the difference between a measured
#    interaction matrix and the stored reference
#
# `SPRINT` wraps that linearized estimator and optionally refreshes the
# zero-point so the finite-difference model can be iterated around a new
# operating point.
#
const MISREG_FIELDS = (:shift_x, :shift_y, :rotation_deg, :radial_scaling, :tangential_scaling)

"""
    MetaSensitivity

Store the linearized sensitivity of an interaction matrix to misregistration
parameters.

`meta` is the inverted meta-sensitivity operator, `calib0` is the reference
interaction matrix at the zero-point, `epsilon` stores the finite-difference
step sizes, and `field_order` defines which misregistration parameters were
included and in what order.
"""
struct MetaSensitivity{T<:AbstractFloat,
    M<:ControlMatrix{T},
    C<:ControlMatrix{T},
    V<:AbstractVector{Symbol}}
    meta::M
    calib0::C
    epsilon::Misregistration{T}
    field_order::V
end

"""
    SPRINT

Stateful wrapper around the meta-sensitivity misregistration estimator.

It stores both the linearized sensitivity model and the currently assumed
zero-point / output misregistration so the estimate can be iterated if the
user requests zero-point updates.
"""
mutable struct SPRINT{T<:AbstractFloat}
    meta::MetaSensitivity{T}
    misregistration_zero::Misregistration{T}
    misregistration_out::Misregistration{T}
    cache_path::Union{Nothing,String}
    save_sensitivity::Bool
    recompute_sensitivity::Bool
    wfs_mis_registered::Bool
end

"""
    compute_meta_sensitivity_matrix(tel, dm, wfs, basis; ...)

Build the finite-difference sensitivity of the interaction matrix to selected
misregistration parameters.

For each requested field, the implementation measures interaction matrices at
positive and negative perturbations around the chosen zero point and forms a
centered derivative. The resulting dense meta-sensitivity matrix is then
inverted through `ControlMatrix`.
"""
function compute_meta_sensitivity_matrix(tel::Telescope, dm::DeformableMirror, wfs::AbstractWFS,
    basis::AbstractMatrix; misregistration_zero::Misregistration=Misregistration(T=eltype(tel.state.opd)),
    epsilon::Misregistration=Misregistration(shift_x=1e-3, shift_y=1e-3, rotation_deg=1e-3, radial_scaling=1e-3,
        tangential_scaling=1e-3, T=eltype(tel.state.opd)),
    n_mis_reg::Int=3, field_order::AbstractVector{Symbol}=collect(MISREG_FIELDS),
    cache_path::Union{Nothing,String}=nothing, save_sensitivity::Bool=true, recompute_sensitivity::Bool=false,
    wfs_mis_registered::Bool=false)

    if cache_path !== nothing && !recompute_sensitivity && isfile(cache_path)
        return deserialize(cache_path)
    end

    T = eltype(tel.state.opd)
    fields = collect(field_order)[1:min(n_mis_reg, length(field_order))]
    dm_model = influence_model(dm)

    supports_dm_misregistration_identification(dm_model, topology(dm)) ||
        throw(UnsupportedAlgorithm("misregistration identification is only supported for grid-backed Gaussian DeformableMirror models"))

    dm0 = DeformableMirror(tel; topology=topology(dm), influence_model=dm_model,
        misregistration=misregistration_zero, T=T)
    calib0 = interaction_matrix(dm0, wfs, tel, basis; amplitude=1e-9)
    calib0_control_matrix = ControlMatrix(calib0.matrix)

    n_elements = length(calib0.matrix)
    meta = zeros(T, n_elements, length(fields))
    for (idx, field) in enumerate(fields)
        eps_val = misregistration_component(epsilon, field)
        if eps_val == 0
            throw(InvalidConfiguration("epsilon for $(field) must be non-zero"))
        end
        if wfs_mis_registered
            if field != :shift_x && field != :shift_y
                throw(InvalidConfiguration("wfs_mis_registered supports shift_x/shift_y only"))
            end
            sx = field == :shift_x ? (misregistration_zero.shift_x + eps_val) : misregistration_zero.shift_x
            sy = field == :shift_y ? (misregistration_zero.shift_y + eps_val) : misregistration_zero.shift_y
            apply_shift_wfs!(wfs; sx=sx, sy=sy)
            imat_p = interaction_matrix(dm, wfs, tel, basis; amplitude=1e-9)
            sx = field == :shift_x ? (misregistration_zero.shift_x - eps_val) : misregistration_zero.shift_x
            sy = field == :shift_y ? (misregistration_zero.shift_y - eps_val) : misregistration_zero.shift_y
            apply_shift_wfs!(wfs; sx=sx, sy=sy)
            imat_n = interaction_matrix(dm, wfs, tel, basis; amplitude=1e-9)
            apply_shift_wfs!(wfs; sx=0, sy=0)
        else
            base_val = misregistration_component(misregistration_zero, field)
            mis_p = update_misregistration(misregistration_zero, field, base_val + eps_val)
            mis_n = update_misregistration(misregistration_zero, field, base_val - eps_val)

            dm_p = DeformableMirror(tel; topology=topology(dm), influence_model=dm_model,
                misregistration=mis_p, T=T)
            dm_n = DeformableMirror(tel; topology=topology(dm), influence_model=dm_model,
                misregistration=mis_n, T=T)
            imat_p = interaction_matrix(dm_p, wfs, tel, basis; amplitude=1e-9)
            imat_n = interaction_matrix(dm_n, wfs, tel, basis; amplitude=1e-9)
        end
        meta[:, idx] .= vec((imat_p.matrix .- imat_n.matrix) ./ (2 * eps_val))
    end

    meta_control_matrix = ControlMatrix(meta)
    out = MetaSensitivity(meta_control_matrix, calib0_control_matrix, epsilon, fields)
    if cache_path !== nothing && save_sensitivity
        serialize(cache_path, out)
    end
    return out
end

"""
    estimate_misregistration(meta, calib_in; misregistration_zero, precision=3, gain_estimation=1)

Estimate misregistration offsets from an input interaction matrix.

This compares `calib_in` to the stored zero-point interaction matrix, projects
that difference through the inverted meta-sensitivity operator, and applies the
resulting parameter offsets to `misregistration_zero`.
"""
function estimate_misregistration(meta::MetaSensitivity, calib_in::AbstractMatrix;
    misregistration_zero::Misregistration, precision::Int=3, gain_estimation::Real=1.0)

    calib_in_control_matrix = ControlMatrix(calib_in)
    diff = vec(calib_in_control_matrix.D .- meta.calib0.D)
    if meta.meta.M === nothing
        throw(InvalidConfiguration("meta sensitivity matrix is not inverted"))
    end
    delta = meta.meta.M * diff
    delta .*= gain_estimation
    delta = round.(delta; digits=precision)

    out = misregistration_zero
    for (i, field) in enumerate(meta.field_order)
        out = update_misregistration(out, field, misregistration_component(out, field) + delta[i])
    end
    return out
end

"""
    SPRINT(tel, dm, wfs, basis; ...)

Construct the iterative misregistration estimator around a chosen zero point.

This first computes the meta-sensitivity matrix and then packages it together
with the zero-point and persistence options used by later `estimate!` calls.
"""
function SPRINT(tel::Telescope, dm::DeformableMirror, wfs::AbstractWFS, basis::AbstractMatrix;
    misregistration_zero::Misregistration=Misregistration(T=eltype(tel.state.opd)),
    epsilon::Misregistration=Misregistration(shift_x=1e-3, shift_y=1e-3, rotation_deg=1e-3, radial_scaling=1e-3,
        tangential_scaling=1e-3, T=eltype(tel.state.opd)),
    n_mis_reg::Int=3, field_order::AbstractVector{Symbol}=collect(MISREG_FIELDS),
    cache_path::Union{Nothing,String}=nothing, save_sensitivity::Bool=true, recompute_sensitivity::Bool=false,
    wfs_mis_registered::Bool=false)

    meta = compute_meta_sensitivity_matrix(tel, dm, wfs, basis;
        misregistration_zero=misregistration_zero,
        epsilon=epsilon,
        n_mis_reg=n_mis_reg,
        field_order=field_order,
        cache_path=cache_path,
        save_sensitivity=save_sensitivity,
        recompute_sensitivity=recompute_sensitivity,
        wfs_mis_registered=wfs_mis_registered)
    return SPRINT(meta, misregistration_zero, misregistration_zero, cache_path, save_sensitivity,
        recompute_sensitivity, wfs_mis_registered)
end

"""
    estimate!(sprint, calib_in; ...)

Run the SPRINT estimator in-place.

With `n_update_zero_point == 0`, this applies one linearized estimate around
the stored zero point. With a positive zero-point update count, the estimator
rebuilds the meta-sensitivity matrix around each newly estimated operating
point before the next update.
"""
function estimate!(sprint::SPRINT, calib_in::AbstractMatrix; precision::Int=3, gain_estimation::Real=1.0,
    n_update_zero_point::Int=0, tel::Union{Nothing,Telescope}=nothing,
    dm::Union{Nothing,DeformableMirror}=nothing, wfs::Union{Nothing,AbstractWFS}=nothing,
    basis::Union{Nothing,AbstractMatrix}=nothing)
    if n_update_zero_point > 0
        if tel === nothing || dm === nothing || wfs === nothing || basis === nothing
            throw(InvalidConfiguration("tel, dm, wfs, and basis are required for zero-point updates"))
        end
        sprint.misregistration_out = sprint.misregistration_zero
        for _ in 1:n_update_zero_point
            sprint.misregistration_out = estimate_misregistration(sprint.meta, calib_in;
                misregistration_zero=sprint.misregistration_out,
                precision=precision,
                gain_estimation=gain_estimation)
            sprint.meta = compute_meta_sensitivity_matrix(tel, dm, wfs, basis;
                misregistration_zero=sprint.misregistration_out,
                epsilon=sprint.meta.epsilon,
                n_mis_reg=length(sprint.meta.field_order),
                field_order=sprint.meta.field_order,
                cache_path=sprint.cache_path,
                save_sensitivity=sprint.save_sensitivity,
                recompute_sensitivity=true,
                wfs_mis_registered=sprint.wfs_mis_registered)
        end
    else
        sprint.misregistration_out = estimate_misregistration(sprint.meta, calib_in;
            misregistration_zero=sprint.misregistration_zero,
            precision=precision,
            gain_estimation=gain_estimation)
    end
    return sprint.misregistration_out
end

"""
    update_misregistration(mis, field, value)

Return a copy of `mis` with one named misregistration field replaced.
"""
function update_misregistration(mis::Misregistration{T}, field::Symbol, value::Real) where {T<:AbstractFloat}
    return Misregistration(; shift_x=field == :shift_x ? T(value) : mis.shift_x,
        shift_y=field == :shift_y ? T(value) : mis.shift_y,
        rotation_deg=field == :rotation_deg ? T(value) : rotation_deg(mis),
        anamorphosis_angle=field == :anamorphosis_angle ? T(value) : anamorphosis_angle_deg(mis),
        tangential_scaling=field == :tangential_scaling ? T(value) : mis.tangential_scaling,
        radial_scaling=field == :radial_scaling ? T(value) : mis.radial_scaling,
        T=T)
end
