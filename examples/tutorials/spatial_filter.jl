include(joinpath(@__DIR__, "common.jl"))

function main(; resolution::Int=24)
    tel = base_telescope(resolution=resolution, central_obstruction=0.0)
    src = base_source()
    apply_demo_ramp!(tel; scale_x=4e-9, scale_y=-2e-9)
    filter = SpatialFilter(tel; shape=CircularFilter(), diameter=resolution ÷ 3, zero_padding=2)
    phase, amplitude = filter!(filter, tel, src)

    @info "Spatial-filter tutorial complete" phase_rms=pupil_rms(phase, tel.state.pupil)
    return (
        filtered_phase=copy(phase),
        filtered_amplitude=copy(amplitude),
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
