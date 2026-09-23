include(joinpath(@__DIR__, "common.jl"))

function compact_tomography_setup()
    atmosphere = TomographyAtmosphereParams(
        zenith_angle_deg=0.0,
        layer_altitudes_m=[0.0],
        L0=25.0,
        r0_zenith=0.2,
        fractional_cn2=[1.0],
        reference_wavelength_m=5.0e-7,
        wind_direction_deg=[0.0],
        wind_speed=[10.0],
    )
    asterism = LGSAsterismParams(
        radius_arcsec=7.6,
        wavelength_m=5.89e-7,
        base_height_m=9.0e4,
        n_lgs=1,
    )
    wfs = LGSWFSParams(
        pupil_diameter_m=8.0,
        n_lenslets=2,
        n_px=4,
        field_stop_size_arcsec=2.0,
        valid_lenslet_map=Bool[1 1; 1 1],
        lenslet_grid_rotations_rad=[0.0],
        lenslet_grid_offsets_fraction=zeros(2, 1),
    )
    tomography = TomographyParams(
        n_fit_src=1,
        fov_optimization_arcsec=0.0,
        fit_src_height_m=Inf,
    )
    dm = TomographyDMParams(
        heights_m=[0.0],
        pitch_m=[0.5],
        cross_coupling=0.2,
        n_actuators=[2],
        valid_actuators=Bool[1 1; 1 1],
    )
    return atmosphere, asterism, wfs, tomography, dm
end

function main()
    atmosphere, asterism, wfs, tomography, dm = compact_tomography_setup()
    recon = build_reconstructor(ModelBasedTomography(), atmosphere, asterism, wfs, tomography, dm)
    command_recon = assemble_reconstructor_and_fitting(
        recon,
        dm;
        n_channels=1,
        slope_order=SimulationSlopes(),
        scaling_factor=1.5e7,
    )
    @info "Tomography matrices prepared for an external RTC" n_phase_samples=size(recon.reconstructor, 1) n_actuators=size(command_recon.matrix, 1)
    return (
        phase_reconstructor=recon,
        command_reconstructor=command_recon,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
