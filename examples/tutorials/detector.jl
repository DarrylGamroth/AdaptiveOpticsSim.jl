include(joinpath(@__DIR__, "common.jl"))

function main(; resolution::Int=24, zero_padding::Int=2)
    tel = base_telescope(resolution=resolution)
    src = base_source()
    pupil = PupilFunction(tel)
    imaging = prepare_direct_imaging(pupil, src; zero_padding=zero_padding)
    rate_map = form_direct_image!(imaging)

    native = Detector(noise=NoiseNone(), integration_time=1.0, qe=1.0, psf_sampling=1, binning=1)
    sampled = Detector(noise=NoiseNone(), integration_time=1.0, qe=1.0, psf_sampling=2, binning=2)
    native_acquisition = prepare_detector_acquisition(native, rate_map)
    sampled_acquisition = prepare_detector_acquisition(sampled, rate_map)
    rng = tutorial_rng()

    frame_native = copy(capture!(native, rate_map, native_acquisition; rng=rng))
    frame_sampled = copy(capture!(sampled, rate_map, sampled_acquisition; rng=rng))

    @info "Detector tutorial complete" native_size=size(frame_native) sampled_size=size(frame_sampled)
    return (
        photon_rate_image=copy(intensity_values(rate_map)),
        frame_native=frame_native,
        frame_sampled=frame_sampled,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
