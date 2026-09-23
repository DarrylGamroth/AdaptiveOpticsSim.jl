using AdaptiveOpticsSim
using AdaptiveOpticsSim.Atmospheres
using AdaptiveOpticsSim.Optics
using AdaptiveOpticsSim.WavefrontSensors
using AdaptiveOpticsSim.Calibration
using AdaptiveOpticsSim.Tomography
import AdaptiveOpticsSim.Optics: filter!
using Logging
using Random
using Statistics

function tutorial_rng(seed::Integer=0)
    return MersenneTwister(seed)
end

function base_telescope(; resolution::Int=32, diameter::Real=8.0,
    central_obstruction::Real=0.1, fov_arcsec::Real=0.0)
    return Telescope(
        resolution=resolution,
        diameter=diameter,
        central_obstruction=central_obstruction,
        fov_arcsec=fov_arcsec,
    )
end

function base_source(; band::Symbol=:I, magnitude::Real=8.0,
    separation_arcsec::Real=0.0, position_angle_deg::Real=0.0)
    return Source(band=band, magnitude=magnitude,
        separation_arcsec=separation_arcsec, position_angle_deg=position_angle_deg)
end

function base_atmosphere(tel::Telescope; r0::Real=0.15,
    reference_wavelength_m::Real=500e-9, L0::Real=25.0)
    return MultiLayerAtmosphere(
        tel;
        r0=r0,
        reference_wavelength_m=reference_wavelength_m,
        L0=L0,
        fractional_cn2=[1.0],
        wind_speed=[8.0],
        wind_direction_deg=[0.0],
        altitude=[0.0],
    )
end

function apply_demo_ramp!(pupil::PupilFunction; scale_x::Real=0.0,
    scale_y::Real=0.0, bias::Real=0.0)
    @inbounds for j in axes(pupil.opd, 2), i in axes(pupil.opd, 1)
        pupil.opd[i, j] = bias + scale_x * (i - 1) + scale_y * (j - 1)
    end
    return pupil
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
