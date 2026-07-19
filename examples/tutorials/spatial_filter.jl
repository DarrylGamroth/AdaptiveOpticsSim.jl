include(joinpath(@__DIR__, "common.jl"))

function main(; resolution::Int=24)
    tel = base_telescope(resolution=resolution, central_obstruction=0.0)
    src = base_source()
    filter = SpatialFilter(tel; shape=CircularFilter(), diameter=resolution ÷ 3, zero_padding=2)
    wavefront = PupilFunction(tel)
    apply_demo_ramp!(wavefront; scale_x=4e-9, scale_y=-2e-9)
    field = ElectricField(wavefront, src; zero_padding=2,
        normalization=DimensionlessNormalization(),
        spatial_measure=PointSampledMeasure(),
        coherence=CoherentFieldCombination())
    formation = prepare_pupil_field(wavefront, src, field;
        center_even_grid=false, amplitude_scale=1)
    fill_electric_field!(field, wavefront, formation)
    filtered = PupilFunction(tel)
    plan = prepare_spatial_filter(tel, filter, field, filtered)
    workspace = SpatialFilterWorkspace(filter)
    filter!(filtered, field, filter, plan, workspace)
    phase = filtered.opd .* (2pi / wavelength(src))

    @info "Spatial-filter tutorial complete" phase_rms=pupil_rms(phase, pupil_mask(tel))
    return (
        filtered_phase=copy(phase),
        filtered_amplitude=copy(filtered.amplitude),
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
