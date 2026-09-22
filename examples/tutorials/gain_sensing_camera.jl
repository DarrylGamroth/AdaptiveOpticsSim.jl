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
        pupil = pupil_mask(tel)[i, j] ? 1.0 : 0.0
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
    psf_zero_padding::Int=2,
    n_iter::Int=6,
    seed::Integer=7,
    atmosphere_step::Real=1e-3,
)
    rng = tutorial_rng(seed)
    pupil = PupilFunction(tel)

    gsc = GainSensingCamera(wfs, basis)
    reset_opd!(pupil)
    calibration_frame = pyramid_modulation_frame(wfs, pupil, ngs)
    calibrate!(gsc, calibration_frame)
    frame = similar(calibration_frame)

    reset_opd!(pupil)
    ngs_imaging = prepare_direct_imaging(pupil, ngs;
        zero_padding=psf_zero_padding)
    sci_imaging = prepare_direct_imaging(pupil, sci;
        zero_padding=psf_zero_padding)
    ngs_image_ref = copy(intensity_values(form_direct_image!(ngs_imaging)))
    sci_image_ref = copy(intensity_values(form_direct_image!(sci_imaging)))
    atmosphere_renderer = prepare_atmosphere_renderer(atm, tel, ngs)
    atmosphere_output = PupilFunction(tel)

    trace = Matrix{Float64}(undef, n_iter, 7)
    for iter in 1:n_iter
        epoch = advance_by!(atm, atmosphere_step; rng=rng)
        render_atmosphere!(atmosphere_output, atmosphere_renderer, atm, epoch)
        apply_opd!(pupil, atmosphere_output.opd)

        trace[iter, 1] = pupil_rms(
            atmosphere_output.opd, pupil_support(pupil)) * 1e9
        ngs_image = intensity_values(form_direct_image!(ngs_imaging))
        trace[iter, 2] = maximum(ngs_image) / maximum(ngs_image_ref)
        pyramid_modulation_frame!(frame, wfs, pupil, ngs)
        og = compute_optical_gains!(gsc, frame)

        sci_image = intensity_values(form_direct_image!(sci_imaging))
        trace[iter, 3] = maximum(sci_image) / maximum(sci_image_ref)
        trace[iter, 4] = norm(frame)
        trace[iter, 5] = norm(og)
        trace[iter, 6] = minimum(abs, og)
        trace[iter, 7] = sum(abs, og) / length(og)
    end

    return trace
end

function main(; resolution::Int=24, pupil_samples::Int=4)
    tel = base_telescope(resolution=resolution, central_obstruction=0.0)
    src = base_source(band=:R, magnitude=8.0)
    sci = base_source(band=:K, magnitude=8.0, coordinates=(0.5, 0.0))
    wfs = PyramidWFS(tel; pupil_samples=pupil_samples, modulation=3.0,
        modulation_points=8, diffraction_padding=2, n_pix_separation=2, n_pix_edge=1)
    basis = cartesian_basis(tel, 4)
    gsc = GainSensingCamera(wfs, basis)
    pupil = PupilFunction(tel)

    reset_opd!(pupil)
    calibration_frame = pyramid_modulation_frame(wfs, pupil, src)
    calibrate!(gsc, calibration_frame)

    coeffs = [20e-9, -12e-9, 8e-9, 0.0]
    apply_opd!(pupil, combine_modes(basis, coeffs))
    frame = similar(calibration_frame)
    pyramid_modulation_frame!(frame, wfs, pupil, src)
    optical_gains = copy(compute_optical_gains!(gsc, frame))

    forcing_coeffs = [
        2.0e-8 -1.0e-8 0.5e-8 -0.25e-8
        1.5e-8 0.5e-8 -0.75e-8 0.5e-8
        -1.0e-8 1.25e-8 0.5e-8 -0.5e-8
        0.5e-8 -0.75e-8 1.0e-8 0.25e-8
    ]
    trace = Matrix{Float64}(undef, size(forcing_coeffs, 1), 3)
    imaging = prepare_direct_imaging(pupil, src; zero_padding=2)
    image_ref = copy(intensity_values(form_direct_image!(imaging)))

    for iter in 1:size(forcing_coeffs, 1)
        opd = combine_modes(basis, @view forcing_coeffs[iter, :])
        apply_opd!(pupil, opd)
        trace[iter, 1] = pupil_rms(pupil.opd, pupil_support(pupil)) * 1e9
        image = intensity_values(form_direct_image!(imaging))
        trace[iter, 2] = maximum(image) / maximum(image_ref)
        pyramid_modulation_frame!(frame, wfs, pupil, src)
        og = compute_optical_gains!(gsc, frame)
        trace[iter, 3] = sum(abs, og) / length(og)
    end

    atm = MultiLayerAtmosphere(
        tel;
        r0=0.15,
        reference_wavelength_m=500e-9,
        L0=25.0,
        fractional_cn2=[0.6, 0.4],
        wind_speed=[10.0, 18.0],
        wind_direction_deg=[0.0, 144.0],
        altitude=[0.0, 5000.0],
    )
    atmosphere_trace = atmosphere_gsc_trace(tel, src, sci, wfs, basis, atm)

    @info "Gain sensing camera tutorial complete" n_modes=length(optical_gains) final_mean_og=trace[end, 3] final_atmosphere_mean_og=atmosphere_trace[end, 7]
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
