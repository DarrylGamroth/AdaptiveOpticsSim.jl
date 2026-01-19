using LinearAlgebra

struct InteractionMatrix{T<:AbstractFloat,M<:AbstractMatrix{T}}
    matrix::M
    amplitude::T
end

function interaction_matrix(dm::DeformableMirror, wfs::AbstractWFS, tel::Telescope; amplitude::Real=1.0)
    n_act = length(dm.state.coefs)
    n_slopes = length(wfs.state.slopes)
    T = eltype(dm.state.coefs)
    mat = zeros(T, n_slopes, n_act)

    opd_base = copy(tel.state.opd)
    coefs = dm.state.coefs

    for k in 1:n_act
        fill!(coefs, zero(T))
        coefs[k] = T(amplitude)
        apply!(dm, tel, DMReplace())
        measure!(wfs, tel)
        mat[:, k] .= wfs.state.slopes
    end

    tel.state.opd .= opd_base
    return InteractionMatrix{T, typeof(mat)}(mat, T(amplitude))
end

function interaction_matrix(dm::DeformableMirror, wfs::AbstractWFS, tel::Telescope,
    commands::AbstractMatrix; amplitude::Real=1.0)
    n_act = length(dm.state.coefs)
    if size(commands, 1) != n_act
        throw(DimensionMismatchError("command matrix row count must match DM coefficients"))
    end
    n_slopes = length(wfs.state.slopes)
    n_modes = size(commands, 2)
    T = eltype(dm.state.coefs)
    mat = zeros(T, n_slopes, n_modes)

    opd_base = copy(tel.state.opd)
    coefs = dm.state.coefs

    for k in 1:n_modes
        @views coefs .= T(amplitude) .* commands[:, k]
        apply!(dm, tel, DMReplace())
        measure!(wfs, tel)
        mat[:, k] .= wfs.state.slopes
    end

    tel.state.opd .= opd_base
    return InteractionMatrix{T, typeof(mat)}(mat, T(amplitude))
end
