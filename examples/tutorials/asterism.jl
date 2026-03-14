include(joinpath(@__DIR__, "common.jl"))

function main(; resolution::Int=24, zero_padding::Int=2)
    tel = base_telescope(resolution=resolution, fov_arcsec=4.0)
    sources = (
        base_source(coordinates=(0.0, 0.0)),
        base_source(coordinates=(1.0, 0.0)),
        base_source(coordinates=(1.0, 120.0)),
        base_source(coordinates=(1.0, 240.0)),
    )
    asterism = Asterism(collect(sources))
    psf = copy(compute_psf!(tel, asterism; zero_padding=zero_padding))

    @info "Asterism tutorial complete" n_sources=length(asterism) psf_size=size(psf)
    return (
        combined_psf=psf,
        per_source_psf=copy(tel.state.psf_stack),
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
