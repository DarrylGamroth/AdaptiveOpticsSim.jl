using LinearAlgebra

const MISREG_FIELDS = (:shift_x, :shift_y, :rotation_deg, :radial_scaling, :tangential_scaling)

struct MetaSensitivity{T<:AbstractFloat}
    meta::CalibrationVault{T,Matrix{T}}
    calib0::CalibrationVault{T,Matrix{T}}
    epsilon::Misregistration{T}
    field_order::Vector{Symbol}
end

mutable struct SPRINT{T<:AbstractFloat}
    meta::MetaSensitivity{T}
    misregistration_zero::Misregistration{T}
    misregistration_out::Misregistration{T}
end

function compute_meta_sensitivity_matrix(tel::Telescope, dm::DeformableMirror, wfs::AbstractWFS,
    basis::AbstractMatrix; misregistration_zero::Misregistration=Misregistration(T=eltype(tel.state.opd)),
    epsilon::Misregistration=Misregistration(shift_x=1e-3, shift_y=1e-3, rotation_deg=1e-3, radial_scaling=1e-3,
        tangential_scaling=1e-3, T=eltype(tel.state.opd)),
    n_mis_reg::Int=3, field_order::AbstractVector{Symbol}=collect(MISREG_FIELDS))

    T = eltype(tel.state.opd)
    fields = collect(field_order)[1:min(n_mis_reg, length(field_order))]

    dm0 = DeformableMirror(tel; n_act=dm.params.n_act, influence_width=dm.params.influence_width,
        misregistration=misregistration_zero, T=T)
    calib0 = interaction_matrix(dm0, wfs, tel, basis; amplitude=1e-9)
    calib0_vault = CalibrationVault(Matrix(calib0.matrix))

    n_elements = length(calib0.matrix)
    meta = zeros(T, n_elements, length(fields))
    for (idx, field) in enumerate(fields)
        eps_val = getfield(epsilon, field)
        if eps_val == 0
            throw(InvalidConfiguration("epsilon for $(field) must be non-zero"))
        end
        mis_p = update_misregistration(misregistration_zero, field, getfield(misregistration_zero, field) + eps_val)
        mis_n = update_misregistration(misregistration_zero, field, getfield(misregistration_zero, field) - eps_val)

        dm_p = DeformableMirror(tel; n_act=dm.params.n_act, influence_width=dm.params.influence_width,
            misregistration=mis_p, T=T)
        dm_n = DeformableMirror(tel; n_act=dm.params.n_act, influence_width=dm.params.influence_width,
            misregistration=mis_n, T=T)
        imat_p = interaction_matrix(dm_p, wfs, tel, basis; amplitude=1e-9)
        imat_n = interaction_matrix(dm_n, wfs, tel, basis; amplitude=1e-9)
        meta[:, idx] .= vec((imat_p.matrix .- imat_n.matrix) ./ (2 * eps_val))
    end

    meta_vault = CalibrationVault(meta)
    return MetaSensitivity(meta_vault, calib0_vault, epsilon, fields)
end

function estimate_misregistration(meta::MetaSensitivity, calib_in::AbstractMatrix;
    misregistration_zero::Misregistration, precision::Int=3, gain_estimation::Real=1.0)

    calib_in_vault = CalibrationVault(Matrix(calib_in))
    diff = vec(calib_in_vault.D .- meta.calib0.D)
    if meta.meta.M === nothing
        throw(InvalidConfiguration("meta sensitivity matrix is not inverted"))
    end
    delta = meta.meta.M * diff
    delta .*= gain_estimation
    delta = round.(delta; digits=precision)

    out = misregistration_zero
    for (i, field) in enumerate(meta.field_order)
        out = update_misregistration(out, field, getfield(out, field) + delta[i])
    end
    return out
end

function SPRINT(tel::Telescope, dm::DeformableMirror, wfs::AbstractWFS, basis::AbstractMatrix;
    misregistration_zero::Misregistration=Misregistration(T=eltype(tel.state.opd)),
    epsilon::Misregistration=Misregistration(shift_x=1e-3, shift_y=1e-3, rotation_deg=1e-3, radial_scaling=1e-3,
        tangential_scaling=1e-3, T=eltype(tel.state.opd)),
    n_mis_reg::Int=3, field_order::AbstractVector{Symbol}=collect(MISREG_FIELDS))

    meta = compute_meta_sensitivity_matrix(tel, dm, wfs, basis;
        misregistration_zero=misregistration_zero,
        epsilon=epsilon,
        n_mis_reg=n_mis_reg,
        field_order=field_order)
    return SPRINT(meta, misregistration_zero, misregistration_zero)
end

function estimate!(sprint::SPRINT, calib_in::AbstractMatrix; precision::Int=3, gain_estimation::Real=1.0)
    sprint.misregistration_out = estimate_misregistration(sprint.meta, calib_in;
        misregistration_zero=sprint.misregistration_zero,
        precision=precision,
        gain_estimation=gain_estimation)
    return sprint.misregistration_out
end

function update_misregistration(mis::Misregistration{T}, field::Symbol, value::Real) where {T<:AbstractFloat}
    return Misregistration(; shift_x=field == :shift_x ? T(value) : mis.shift_x,
        shift_y=field == :shift_y ? T(value) : mis.shift_y,
        rotation_deg=field == :rotation_deg ? T(value) : mis.rotation_deg,
        anamorphosis_angle=field == :anamorphosis_angle ? T(value) : mis.anamorphosis_angle,
        tangential_scaling=field == :tangential_scaling ? T(value) : mis.tangential_scaling,
        radial_scaling=field == :radial_scaling ? T(value) : mis.radial_scaling,
        T=T)
end
