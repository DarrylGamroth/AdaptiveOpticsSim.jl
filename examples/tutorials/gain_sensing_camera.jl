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

function main(; resolution::Int=24, n_subap::Int=4)
    tel = base_telescope(resolution=resolution, central_obstruction=0.0)
    src = base_source(band=:R, magnitude=8.0)
    wfs = PyramidWFS(tel; n_subap=n_subap, mode=Diffractive(), threshold=0.5, modulation=3.0,
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

    @info "Gain sensing camera tutorial complete" n_modes=length(optical_gains) final_mean_og=trace[end, 3]
    return (
        coeffs=coeffs,
        calibration_frame=calibration_frame,
        frame=frame,
        optical_gains=optical_gains,
        trace=trace,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
