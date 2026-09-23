include(joinpath(@__DIR__, "common.jl"))
using LinearAlgebra
import AdaptiveOpticsCalibration

const OpticalGains = AdaptiveOpticsCalibration.OpticalGains

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

function centered_focal_basis(pupil_basis::AbstractArray{T,3},
    mask::AbstractMatrix{Complex{T}}) where {T<:AbstractFloat}
    side = size(mask, 1)
    size(mask, 2) == side || throw(ArgumentError(
        "Pyramid focal mask must be square for complete-image gain sensing",
    ))
    size(pupil_basis, 1) == size(pupil_basis, 2) || throw(ArgumentError(
        "modal pupil basis must be square for complete-image gain sensing",
    ))
    side >= size(pupil_basis, 1) || throw(ArgumentError(
        "Pyramid focal mask must not be smaller than the modal pupil basis",
    ))
    iseven(side - size(pupil_basis, 1)) || throw(ArgumentError(
        "Pyramid focal mask and modal pupil basis must have aligned centers",
    ))

    basis = zeros(T, side, side, size(pupil_basis, 3))
    offset = div(side - size(pupil_basis, 1), 2)
    @views basis[offset+1:offset+size(pupil_basis, 1),
        offset+1:offset+size(pupil_basis, 2), :] .= pupil_basis
    return basis
end

function estimate_optical_gains!(product, workspace, plan,
    frame::AbstractMatrix{T}) where {T<:AbstractFloat}
    AdaptiveOpticsCalibration.process!(product, workspace, plan,
        OpticalGains.GainSensingInputs(frame))
    return OpticalGains.optical_gains(product)
end

function atmosphere_gsc_trace(
    tel::Telescope,
    ngs::Source,
    sci::Source,
    wfs::PyramidWFS,
    atm::AbstractAtmosphere;
    gain_sensing_plan,
    gain_sensing_product,
    gain_sensing_workspace,
    psf_zero_padding::Int=2,
    n_iter::Int=6,
    seed::Integer=7,
    atmosphere_step::Real=1e-3,
)
    rng = tutorial_rng(seed)
    pupil = PupilFunction(tel)

    reset_opd!(pupil)
    frame = pyramid_modulation_frame(wfs, pupil, ngs)

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
        og = estimate_optical_gains!(gain_sensing_product,
            gain_sensing_workspace, gain_sensing_plan, frame)

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
    pupil_basis = cartesian_basis(tel, 4)
    pupil = PupilFunction(tel)

    reset_opd!(pupil)
    calibration_frame = pyramid_modulation_frame(wfs, pupil, src)
    focal_mask = pyramid_focal_mask(wfs)
    focal_basis = centered_focal_basis(pupil_basis, focal_mask)
    specification = OpticalGains.GainSensingSpecification(
        focal_mask, focal_basis, calibration_frame,
    )
    gain_sensing_plan = AdaptiveOpticsCalibration.prepare(
        OpticalGains.GainSensing(), specification,
    )
    gain_sensing_product = AdaptiveOpticsCalibration.allocate_result(gain_sensing_plan)
    gain_sensing_workspace = AdaptiveOpticsCalibration.allocate_workspace(gain_sensing_plan)

    coeffs = [20e-9, -12e-9, 8e-9, 0.0]
    apply_opd!(pupil, combine_modes(pupil_basis, coeffs))
    frame = similar(calibration_frame)
    pyramid_modulation_frame!(frame, wfs, pupil, src)
    optical_gains = copy(estimate_optical_gains!(gain_sensing_product,
        gain_sensing_workspace, gain_sensing_plan, frame))

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
        opd = combine_modes(pupil_basis, @view forcing_coeffs[iter, :])
        apply_opd!(pupil, opd)
        trace[iter, 1] = pupil_rms(pupil.opd, pupil_support(pupil)) * 1e9
        image = intensity_values(form_direct_image!(imaging))
        trace[iter, 2] = maximum(image) / maximum(image_ref)
        pyramid_modulation_frame!(frame, wfs, pupil, src)
        og = estimate_optical_gains!(gain_sensing_product,
            gain_sensing_workspace, gain_sensing_plan, frame)
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
    atmosphere_trace = atmosphere_gsc_trace(tel, src, sci, wfs, atm;
        gain_sensing_plan=gain_sensing_plan,
        gain_sensing_product=gain_sensing_product,
        gain_sensing_workspace=gain_sensing_workspace,
    )

    @info "Gain-sensing tutorial complete" n_modes=length(optical_gains) final_mean_og=trace[end, 3] final_atmosphere_mean_og=atmosphere_trace[end, 7]
    return (
        coeffs=coeffs,
        calibration_frame=calibration_frame,
        focal_mask=focal_mask,
        focal_basis=focal_basis,
        frame=frame,
        optical_gains=optical_gains,
        trace=trace,
        atmosphere_trace=atmosphere_trace,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
