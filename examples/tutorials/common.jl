using AdaptiveOptics
using Logging
using Random
using Statistics

function tutorial_rng(seed::Integer=0)
    return MersenneTwister(seed)
end

function base_telescope(; resolution::Int=32, diameter::Real=8.0, sampling_time::Real=1e-3,
    central_obstruction::Real=0.1, fov_arcsec::Real=0.0)
    return Telescope(
        resolution=resolution,
        diameter=diameter,
        sampling_time=sampling_time,
        central_obstruction=central_obstruction,
        fov_arcsec=fov_arcsec,
    )
end

function base_source(; band::Symbol=:I, magnitude::Real=8.0, coordinates::Tuple{<:Real,<:Real}=(0.0, 0.0))
    return Source(band=band, magnitude=magnitude, coordinates=coordinates)
end

function base_atmosphere(tel::Telescope; r0::Real=0.15, L0::Real=25.0)
    return MultiLayerAtmosphere(
        tel;
        r0=r0,
        L0=L0,
        fractional_r0=[1.0],
        wind_speed=[8.0],
        wind_direction=[0.0],
        altitude=[0.0],
    )
end

function apply_demo_ramp!(tel::Telescope; scale_x::Real=0.0, scale_y::Real=0.0, bias::Real=0.0)
    @inbounds for j in axes(tel.state.opd, 2), i in axes(tel.state.opd, 1)
        tel.state.opd[i, j] = bias + scale_x * (i - 1) + scale_y * (j - 1)
    end
    return tel
end

function pupil_rms(opd::AbstractMatrix, pupil::AbstractMatrix{Bool})
    vals = opd[pupil]
    return sqrt(sum(abs2, vals) / length(vals))
end

function combine_modes(basis::AbstractArray{T,3}, coeffs::AbstractVector{<:Real}) where {T}
    opd = zeros(T, size(basis, 1), size(basis, 2))
    n_modes = min(size(basis, 3), length(coeffs))
    @inbounds for k in 1:n_modes
        @views @. opd += T(coeffs[k]) * basis[:, :, k]
    end
    return opd
end

function run_closed_loop_example(make_wfs::Function; n_iter::Int=4, seed::Integer=0,
    resolution::Int=16, n_subap::Int=4, n_act::Int=3, amplitude::Real=1e-9, gain::Real=0.4)
    rng = tutorial_rng(seed)
    tel = base_telescope(resolution=resolution, central_obstruction=0.0)
    src = base_source()
    atm = base_atmosphere(tel)
    dm = DeformableMirror(tel; n_act=n_act, influence_width=0.35)
    wfs = make_wfs(tel, n_subap)
    imat = interaction_matrix(dm, wfs, tel, src; amplitude=amplitude)
    recon = ModalReconstructor(imat; gain=gain)
    cmd = similar(dm.state.coefs)
    residual_before = zeros(Float64, n_iter)
    residual_after = similar(residual_before)

    for k in 1:n_iter
        advance!(atm, tel; rng=rng)
        propagate!(atm, tel)
        residual_before[k] = pupil_rms(tel.state.opd, tel.state.pupil)
        measure!(wfs, tel, src)
        reconstruct!(cmd, recon, wfs.state.slopes)
        dm.state.coefs .= -cmd
        apply!(dm, tel, DMAdditive())
        residual_after[k] = pupil_rms(tel.state.opd, tel.state.pupil)
    end

    psf = copy(compute_psf!(tel, src; zero_padding=2))
    return (
        residual_before=residual_before,
        residual_after=residual_after,
        final_psf=psf,
        final_slopes=copy(wfs.state.slopes),
    )
end
