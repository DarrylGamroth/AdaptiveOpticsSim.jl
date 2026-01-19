using LinearAlgebra

struct DeformableMirrorParams{T<:AbstractFloat}
    n_act::Int
    influence_width::T
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

function DeformableMirror(tel::Telescope; n_act::Int, influence_width::Real=0.2, T::Type{<:AbstractFloat}=Float64, backend=Array)
    params = DeformableMirrorParams{T}(n_act, T(influence_width))
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
    n = tel.params.resolution
    n_act = dm.params.n_act
    sigma = dm.params.influence_width
    xs = range(-1.0, 1.0; length=n_act)
    ys = range(-1.0, 1.0; length=n_act)

    cx = (n + 1) / 2
    cy = (n + 1) / 2
    scale = n / 2

    idx = 1
    @inbounds for x0 in xs, y0 in ys
        for i in 1:n, j in 1:n
            x = (i - cx) / scale
            y = (j - cy) / scale
            r2 = (x - x0)^2 + (y - y0)^2
            dm.state.modes[(j - 1) * n + i, idx] = exp(-r2 / (2 * sigma^2)) * tel.state.pupil[i, j]
        end
        idx += 1
    end
    return dm
end

function apply!(dm::DeformableMirror, tel::Telescope; additive::Bool=true)
    opd_vec = reshape(dm.state.opd, :)
    mul!(opd_vec, dm.state.modes, dm.state.coefs)
    if additive
        tel.state.opd .+= dm.state.opd
    else
        tel.state.opd .= dm.state.opd
    end
    return tel
end
