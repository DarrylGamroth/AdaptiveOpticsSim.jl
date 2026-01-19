using LinearAlgebra

struct DeformableMirrorParams{T<:AbstractFloat}
    n_act::Int
    influence_width::T
    misregistration::Misregistration{T}
end

mutable struct DeformableMirrorState{T<:AbstractFloat,A<:AbstractMatrix{T},M<:AbstractMatrix{T},V<:AbstractVector{T}}
    opd::A
    modes::M
    coefs::V
end

struct DeformableMirror{P<:DeformableMirrorParams,S<:DeformableMirrorState} <: AbstractDeformableMirror
    params::P
    state::S
end

function DeformableMirror(tel::Telescope; n_act::Int, influence_width::Real=0.2,
    T::Type{<:AbstractFloat}=Float64, misregistration::Misregistration=Misregistration(T=T), backend=Array)
    params = DeformableMirrorParams{T}(n_act, T(influence_width), misregistration)
    n = tel.params.resolution
    opd = backend{T}(undef, n, n)
    fill!(opd, zero(T))
    modes = backend{T}(undef, n * n, n_act * n_act)
    coefs = backend{T}(undef, n_act * n_act)
    fill!(coefs, zero(T))
    state = DeformableMirrorState{T, typeof(opd), typeof(modes), typeof(coefs)}(opd, modes, coefs)
    dm = DeformableMirror(params, state)
    build_influence_functions!(dm, tel)
    return dm
end

function build_influence_functions!(dm::DeformableMirror, tel::Telescope)
    pupil = tel.state.pupil
    ax1 = axes(pupil, 1)
    ax2 = axes(pupil, 2)
    n = length(ax1)
    n_act = dm.params.n_act
    sigma = dm.params.influence_width
    xs = range(-1.0, 1.0; length=n_act)
    ys = range(-1.0, 1.0; length=n_act)

    cx = (n + 1) / 2
    cy = (n + 1) / 2
    scale = n / 2
    i0 = first(ax1)
    j0 = first(ax2)

    idx = 1
    @inbounds for x0 in xs, y0 in ys
        x_m, y_m = apply_misregistration(x0, y0, dm.params.misregistration)
        for i in ax1, j in ax2
            ii = i - i0 + 1
            jj = j - j0 + 1
            x = (ii - cx) / scale
            y = (jj - cy) / scale
            r2 = (x - x_m)^2 + (y - y_m)^2
            dm.state.modes[(jj - 1) * n + ii, idx] = exp(-r2 / (2 * sigma^2)) * pupil[i, j]
        end
        idx += 1
    end
    return dm
end

function apply!(dm::DeformableMirror, tel::Telescope, ::DMAdditive)
    opd_vec = reshape(dm.state.opd, :)
    mul!(opd_vec, dm.state.modes, dm.state.coefs)
    tel.state.opd .+= dm.state.opd
    return tel
end

function apply!(dm::DeformableMirror, tel::Telescope, ::DMReplace)
    opd_vec = reshape(dm.state.opd, :)
    mul!(opd_vec, dm.state.modes, dm.state.coefs)
    tel.state.opd .= dm.state.opd
    return tel
end
