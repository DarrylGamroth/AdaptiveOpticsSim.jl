include(joinpath(@__DIR__, "common.jl"))
using LinearAlgebra

function cartesian_basis(tel::Telescope, n_modes::Int)
    n = tel.params.resolution
    basis = zeros(Float64, n, n, n_modes)
    x = collect(range(-1.0, 1.0; length=n + 1))[1:n]
    y = collect(range(-1.0, 1.0; length=n + 1))[1:n]
    @inbounds for j in 1:n, i in 1:n
        px = x[i]
        py = y[j]
        pupil = tel.state.pupil[i, j] ? 1.0 : 0.0
        if n_modes >= 1
            basis[i, j, 1] = pupil * px
        end
        if n_modes >= 2
            basis[i, j, 2] = pupil * py
        end
        if n_modes >= 3
            basis[i, j, 3] = pupil * px * py
        end
        if n_modes >= 4
            basis[i, j, 4] = pupil * (px * px - py * py)
        end
    end
    return basis
end

function combine_modes!(out::AbstractMatrix{T}, basis::AbstractArray{<:Real,3},
    coeffs::AbstractVector{<:Real}) where {T<:AbstractFloat}
    fill!(out, zero(T))
    n_modes = min(size(basis, 3), length(coeffs))
    @inbounds for k in 1:n_modes
        coeff = T(coeffs[k])
        @views @. out += coeff * basis[:, :, k]
    end
    return out
end

function atmosphere_gsc_trace(
    tel::Telescope,
    ngs::Source,
    sci::Source,
    wfs::PyramidWFS,
    basis::AbstractArray{<:Real,3},
    atm::AbstractAtmosphere;
    gain::Real=0.2,
    frame_delay::Int=2,
    calibration_amplitude::Real=1e-9,
    psf_zero_padding::Int=2,
    og_floor::Real=0.05,
    n_iter::Int=6,
    seed::Integer=7,
)
    rng = tutorial_rng(seed)
    H = zeros(Float64, length(wfs.state.slopes), size(basis, 3))
    for k in 1:size(basis, 3)
        apply_opd!(tel, calibration_amplitude .* view(basis, :, :, k))
        measure!(wfs, tel, ngs)
        H[:, k] .= wfs.state.slopes
    end
    reset_opd!(tel)

    recon = pinv(H)
    control_coeffs = zeros(Float64, size(basis, 3))
    delayed_slopes = zeros(Float64, size(H, 1))
    forcing_opd = zeros(Float64, size(tel.state.opd))
    correction_opd = similar(forcing_opd)
    residual_opd = similar(forcing_opd)

    gsc = GainSensingCamera(wfs.state.pyramid_mask, basis)
    calibration_frame = similar(wfs.state.intensity)
    reset_opd!(tel)
    pyramid_modulation_frame!(calibration_frame, wfs, tel, ngs)
    calibrate!(gsc, calibration_frame)
    frame = similar(calibration_frame)
    og_safe = similar(gsc.og)

    ngs_psf_ref = begin
        reset_opd!(tel)
        copy(compute_psf!(tel, ngs; zero_padding=psf_zero_padding))
    end
    sci_psf_ref = begin
        reset_opd!(tel)
        copy(compute_psf!(tel, sci; zero_padding=psf_zero_padding))
    end

    trace = Matrix{Float64}(undef, n_iter, 7)
    for iter in 1:n_iter
        advance!(atm, tel; rng=rng)
        propagate!(atm, tel)
        forcing_opd .= tel.state.opd
        combine_modes!(correction_opd, basis, control_coeffs)
        @. residual_opd = forcing_opd - correction_opd
        apply_opd!(tel, residual_opd)

        trace[iter, 1] = pupil_rms(forcing_opd, tel.state.pupil) * 1e9
        trace[iter, 2] = pupil_rms(residual_opd, tel.state.pupil) * 1e9
        ngs_psf = compute_psf!(tel, ngs; zero_padding=psf_zero_padding)
        trace[iter, 4] = maximum(ngs_psf) / maximum(ngs_psf_ref)
        slopes = measure!(wfs, tel, ngs)
        trace[iter, 6] = norm(slopes)
        pyramid_modulation_frame!(frame, wfs, tel, ngs)
        og = compute_optical_gains!(gsc, frame)
        @. og_safe = max(abs(og), og_floor)

        trace[iter, 3] = pupil_rms(residual_opd, tel.state.pupil) * 1e9
        sci_psf = compute_psf!(tel, sci; zero_padding=psf_zero_padding)
        trace[iter, 5] = maximum(sci_psf) / maximum(sci_psf_ref)

        if frame_delay == 1
            delayed_slopes .= slopes
        end
        control_coeffs .+= gain .* ((recon * delayed_slopes) ./ og_safe)
        trace[iter, 7] = sum(og_safe) / length(og_safe)
        if frame_delay == 2
            delayed_slopes .= slopes
        end
    end

    return trace
end

function main(; resolution::Int=24, pupil_samples::Int=4)
    tel = base_telescope(resolution=resolution, central_obstruction=0.0)
    src = base_source(band=:R, magnitude=8.0)
    sci = base_source(band=:K, magnitude=8.0, coordinates=(0.5, 0.0))
    wfs = PyramidWFS(tel; pupil_samples=pupil_samples, mode=Diffractive(), threshold=0.5, modulation=3.0,
        normalization=IncidenceFluxNormalization(),
        modulation_points=8, diffraction_padding=2, n_pix_separation=2, n_pix_edge=1)
    basis = cartesian_basis(tel, 4)
    gsc = GainSensingCamera(wfs.state.pyramid_mask, basis)

    calibration_frame = similar(wfs.state.intensity)
    reset_opd!(tel)
    pyramid_modulation_frame!(calibration_frame, wfs, tel, src)
    calibrate!(gsc, calibration_frame)

    coeffs = [20e-9, -12e-9, 8e-9, 0.0]
    apply_opd!(tel, combine_modes(basis, coeffs))
    frame = similar(calibration_frame)
    pyramid_modulation_frame!(frame, wfs, tel, src)
    optical_gains = copy(compute_optical_gains!(gsc, frame))

    forcing_coeffs = [
        2.0e-8 -1.0e-8 0.5e-8 -0.25e-8
        1.5e-8 0.5e-8 -0.75e-8 0.5e-8
        -1.0e-8 1.25e-8 0.5e-8 -0.5e-8
        0.5e-8 -0.75e-8 1.0e-8 0.25e-8
    ]
    H = zeros(Float64, length(wfs.state.slopes), size(basis, 3))
    for k in 1:size(basis, 3)
        apply_opd!(tel, 1e-9 .* view(basis, :, :, k))
        measure!(wfs, tel, src)
        H[:, k] .= wfs.state.slopes
    end
    reset_opd!(tel)
    recon = pinv(H)
    delayed_slopes = zeros(Float64, size(H, 1))
    control_coeffs = zeros(Float64, size(basis, 3))
    trace = Matrix{Float64}(undef, size(forcing_coeffs, 1), 3)
    psf_ref = copy(compute_psf!(tel, src; zero_padding=2))

    for iter in 1:size(forcing_coeffs, 1)
        opd = combine_modes(basis, @view forcing_coeffs[iter, :]) .- combine_modes(basis, control_coeffs)
        apply_opd!(tel, opd)
        trace[iter, 1] = pupil_rms(tel.state.opd, tel.state.pupil) * 1e9
        psf = compute_psf!(tel, src; zero_padding=2)
        trace[iter, 2] = maximum(psf) / maximum(psf_ref)
        slopes = measure!(wfs, tel, src)
        pyramid_modulation_frame!(frame, wfs, tel, src)
        og = compute_optical_gains!(gsc, frame)
        control_coeffs .+= 0.2 .* ((recon * delayed_slopes) ./ max.(abs.(og), 0.05))
        delayed_slopes .= slopes
        trace[iter, 3] = sum(og) / length(og)
    end

    atm = MultiLayerAtmosphere(
        tel;
        r0=0.15,
        L0=25.0,
        fractional_cn2=[0.6, 0.4],
        wind_speed=[10.0, 18.0],
        wind_direction=[0.0, 144.0],
        altitude=[0.0, 5000.0],
    )
    atmosphere_trace = atmosphere_gsc_trace(tel, src, sci, wfs, basis, atm)

    @info "Gain sensing camera tutorial complete" n_modes=length(optical_gains) final_mean_og=trace[end, 3] final_loop_mean_og=atmosphere_trace[end, 7]
    return (
        coeffs=coeffs,
        calibration_frame=calibration_frame,
        frame=frame,
        optical_gains=optical_gains,
        trace=trace,
        atmosphere_trace=atmosphere_trace,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
