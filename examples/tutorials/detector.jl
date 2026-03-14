include(joinpath(@__DIR__, "common.jl"))

function main(; resolution::Int=24, zero_padding::Int=2)
    tel = base_telescope(resolution=resolution)
    src = base_source()
    psf = copy(compute_psf!(tel, src; zero_padding=zero_padding))

    native = Detector(noise=NoiseNone(), integration_time=1.0, qe=1.0, psf_sampling=1, binning=1)
    sampled = Detector(noise=NoiseNone(), integration_time=1.0, qe=1.0, psf_sampling=2, binning=2)

    frame_native = copy(capture!(native, psf))
    frame_sampled = copy(capture!(sampled, psf))

    @info "Detector tutorial complete" native_size=size(frame_native) sampled_size=size(frame_sampled)
    return (
        psf=psf,
        frame_native=frame_native,
        frame_sampled=frame_sampled,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
