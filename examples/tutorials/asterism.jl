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
    pupil = PupilFunction(tel)
    imaging = prepare_direct_imaging(pupil, asterism; zero_padding=zero_padding)
    form_direct_image!(imaging)
    combined = copy(intensity_values(direct_imaging_output(imaging)))
    component_maps = map(direct_imaging_components(imaging)) do component
        copy(intensity_values(direct_imaging_output(component)))
    end

    @info "Asterism tutorial complete" n_sources=length(asterism) image_size=size(combined)
    return (
        combined_image=combined,
        component_images=component_maps,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
