@testset "Modal bases and fitting" begin
    tel = Telescope(resolution=16, diameter=8.0, central_obstruction=0.0)
    dm = DeformableMirror(tel; n_act=2, influence_width=0.4)
    basis = @inferred modal_basis(dm, tel; n_modes=2)
    basis_without_projector = @inferred modal_basis(
        dm, tel; n_modes=2, projector=false)
    @test typeof(basis_without_projector) === typeof(basis)
    @test basis_without_projector.projector === nothing
    @test !isdefined(AdaptiveOpticsSim, :KLDMModes)
    @test !isdefined(Calibration, :KLDMModes)
    @test !isdefined(AdaptiveOpticsSim, :KLBasis)
    @test !isdefined(Calibration, :KLBasis)
    @test !isdefined(AdaptiveOpticsSim, :KLHHtPSD)
    @test !isdefined(Calibration, :KLHHtPSD)
    @test !isdefined(AdaptiveOpticsSim, :kl_modal_basis)
    @test !isdefined(Calibration, :kl_modal_basis)
    @test !isdefined(AdaptiveOpticsSim, :fitting_error)
    @test !isdefined(Calibration, :fitting_error)
    @test parentmodule(KarhunenLoeveBasis) === AOCModalBases
    @test size(basis.M2C, 2) == 2
    sampled_influences = Matrix(sampled_influence_matrix(dm))
    direct_plan = AdaptiveOpticsCalibration.prepare(
        AOCModalBases.InfluenceFunctionEigenbasis(),
        AOCModalBases.SampledInfluenceBasisSpecification(
            size(sampled_influences, 1),
            size(sampled_influences, 2),
            2,
            vec(pupil_mask(tel)),
            eltype(sampled_influences),
        ),
    )
    direct_basis = AdaptiveOpticsCalibration.process(direct_plan, sampled_influences)
    direct_m2c = AOCModalBases.modal_to_command(direct_basis)
    direct_sampled_modes = AOCModalBases.sampled_modes(direct_basis)
    for mode in 1:2
        sign = dot(
            @view(basis.M2C[:, mode]),
            @view(direct_m2c[:, mode]),
        ) < 0 ? -1 : 1
        @test @view(basis.M2C[:, mode]) ≈
            sign .* (@view direct_m2c[:, mode])
        @test @view(basis.basis[:, mode]) ≈
            sign .* (@view direct_sampled_modes[:, mode])
    end
    opd = rand(16, 16)
    fitting_plan = AdaptiveOpticsCalibration.prepare(
        AOCModalBases.ModalFitting(),
        AOCModalBases.ModalFittingSpecification(
            size(opd, 1),
            size(opd, 2),
            size(basis.basis, 2),
            eltype(opd),
        ),
    )
    fitting = AdaptiveOpticsCalibration.process(
        fitting_plan,
        AOCModalBases.ModalFittingInputs(opd, basis.projector, basis.basis),
    )
    @test size(AOCModalBases.residual_opd(fitting)) == size(opd)
    @test size(AOCModalBases.fitted_opd(fitting)) == size(opd)
    @test size(AOCModalBases.input_opd(fitting)) == size(opd)

    atm = KolmogorovAtmosphere(tel; r0=0.2,
        reference_wavelength_m=TEST_ATMOSPHERE_REFERENCE_WAVELENGTH_M,
        L0=25.0)
    atmospheric_basis = modal_basis(
        dm,
        tel;
        n_modes=2,
        projector=true,
        method=KarhunenLoeveBasis(),
        atm=atm,
    )
    @test size(atmospheric_basis.M2C, 2) == 2
    @test size(atmospheric_basis.basis, 2) == 2
    sampled_influences = sampled_influence_matrix(dm)
    @test atmospheric_basis.basis ≈
        sampled_influences * atmospheric_basis.M2C atol=2e-12
    support = vec(Array(pupil_mask(tel)))
    @test atmospheric_basis.basis[support, :]' *
          atmospheric_basis.basis[support, :] / count(support) ≈
        Matrix{Float64}(I, 2, 2) atol=2e-12
    @test maximum(abs, atmospheric_basis.basis[.!support, :]) == 0
    expected_atmospheric_projector =
        atmospheric_basis.basis' * Diagonal(Float64.(support)) / count(support)
    @test atmospheric_basis.projector ≈ expected_atmospheric_projector atol=2e-14
    @test maximum(abs, atmospheric_basis.projector[:, .!support]) == 0
    pupil_opd = randn(length(support))
    exterior_changed_opd = copy(pupil_opd)
    exterior_changed_opd[.!support] .= 1e6
    @test atmospheric_basis.projector * exterior_changed_opd ≈
        atmospheric_basis.projector * pupil_opd atol=2e-12

    default_atmospheric_basis = @inferred modal_basis(
        dm,
        tel;
        method=KarhunenLoeveBasis(),
        atm=atm,
    )
    @test size(default_atmospheric_basis.M2C, 2) ==
        size(sampled_influences, 2) - 1

    sampling_m = tel.aperture.sampling_m[1]
    covariance = projected_atmospheric_opd_covariance(
        Matrix(sampled_influences),
        support,
        tel.params.resolution,
        sampling_m,
        atm,
    )
    @test covariance ≈ covariance' atol=2e-28
    @test minimum(eigvals(Symmetric(covariance))) >= -1e-25
    coordinate_scaling = Diagonal(range(0.5, 1.5; length=size(sampled_influences, 2)))
    scaled_covariance = projected_atmospheric_opd_covariance(
        Matrix(sampled_influences) * coordinate_scaling,
        support,
        tel.params.resolution,
        sampling_m,
        atm,
    )
    @test scaled_covariance ≈ coordinate_scaling' * covariance *
        coordinate_scaling rtol=2e-12 atol=2e-28
    outside_changed = Matrix(sampled_influences)
    outside_changed[.!support, :] .= 100
    @test projected_atmospheric_opd_covariance(
        outside_changed,
        support,
        tel.params.resolution,
        sampling_m,
        atm,
    ) ≈ covariance rtol=2e-12 atol=2e-28

    oracle_resolution = 2
    oracle_sampling_m = 0.7
    oracle_influences = [
        1.0 0.2
        -0.3 0.7
        0.4 -0.5
        0.8 0.1
    ]
    oracle_support = trues(4)
    oracle_telescope = Telescope(
        resolution=oracle_resolution,
        diameter=oracle_resolution * oracle_sampling_m,
        central_obstruction=0.0,
    )
    oracle_atmosphere = KolmogorovAtmosphere(
        oracle_telescope;
        r0=0.31,
        reference_wavelength_m=632.8e-9,
        L0=18.0,
    )
    oracle_covariance = projected_atmospheric_opd_covariance(
        oracle_influences,
        oracle_support,
        oracle_resolution,
        oracle_sampling_m,
        oracle_atmosphere,
    )
    spectral_resolution = 2 * oracle_resolution
    frequency_step = 1 / (spectral_resolution * oracle_sampling_m)
    frequencies = [0.0, frequency_step, 2 * frequency_step, -frequency_step]
    direct_covariance = zeros(2, 2)
    for frequency_column in 1:spectral_resolution,
        frequency_row in 1:spectral_resolution
        fx = frequencies[frequency_row]
        fy = frequencies[frequency_column]
        radial_frequency_squared = fx^2 + fy^2
        phase_psd = 0.023 * 0.31^(-5 / 3) *
                    (radial_frequency_squared + 18.0^(-2))^(-11 / 6)
        opd_psd = phase_psd * (632.8e-9 / (2π))^2
        transformed = zeros(ComplexF64, 2)
        for corrector in 1:2, column in 1:oracle_resolution,
            row in 1:oracle_resolution
            sample = (column - 1) * oracle_resolution + row
            phase = -2π * ((frequency_row - 1) * (row - 1) +
                              (frequency_column - 1) * (column - 1)) /
                    spectral_resolution
            transformed[corrector] +=
                oracle_influences[sample, corrector] * cis(phase) / 4
        end
        direct_covariance .+=
            real.(conj.(transformed) * transpose(transformed)) .* opd_psd .* frequency_step^2
    end
    @test oracle_covariance ≈ direct_covariance rtol=3e-14 atol=1e-30
    doubled_wavelength_atmosphere = KolmogorovAtmosphere(
        tel;
        r0=0.2,
        reference_wavelength_m=2 * TEST_ATMOSPHERE_REFERENCE_WAVELENGTH_M,
        L0=25.0,
    )
    doubled_wavelength_covariance = projected_atmospheric_opd_covariance(
        Matrix(sampled_influences),
        support,
        tel.params.resolution,
        sampling_m,
        doubled_wavelength_atmosphere,
    )
    @test doubled_wavelength_covariance ≈ 4 .* covariance rtol=2e-12 atol=2e-28

    telescope_f32 = Telescope(
        resolution=8,
        diameter=8.0f0,
        central_obstruction=0.0f0,
        T=Float32,
    )
    dm_f32 = DeformableMirror(
        telescope_f32;
        n_act=2,
        influence_width=0.4f0,
        T=Float32,
    )
    atmosphere_f32 = KolmogorovAtmosphere(
        telescope_f32;
        r0=0.2f0,
        reference_wavelength_m=Float32(TEST_ATMOSPHERE_REFERENCE_WAVELENGTH_M),
        L0=25.0f0,
    )
    atmospheric_basis_f32 = @inferred modal_basis(
        dm_f32,
        telescope_f32;
        n_modes=2,
        projector=false,
        method=KarhunenLoeveBasis(),
        atm=atmosphere_f32,
    )
    @test eltype(atmospheric_basis_f32.M2C) === Float32
    @test eltype(atmospheric_basis_f32.basis) === Float32
end

@testset "Interaction-matrix WFS output finalization" begin
    T = Float32
    tel = Telescope(
        resolution=480,
        diameter=8.0,
        central_obstruction=0.0,
        T=T,
    )
    src = Source(band=:R, magnitude=3.0, T=T)
    dm = DeformableMirror(tel; n_act=2, influence_width=T(0.2), T=T)
    wfs = PyramidWFS(
        tel;
        pupil_samples=20,
        threshold=T(0.1),
        modulation=zero(T),
        modulation_points=1,
        light_ratio=T(0.1),
        n_pix_separation=4,
        n_pix_edge=2,
        psf_centering=true,
        mode=Diffractive(),
        T=T,
    )
    initial_rows = length(slopes(wfs))
    imat = interaction_matrix(
        dm,
        wfs,
        PupilFunction(tel; T=T),
        src;
        amplitude=T(5e-9),
    )

    @test length(slopes(wfs)) < initial_rows
    @test size(imat.matrix) ==
        (length(slopes(wfs)), length(dm.state.coefs))
    @test all(isfinite, imat.matrix)

    stale_wfs = PyramidWFS(
        tel;
        pupil_samples=20,
        threshold=T(0.1),
        modulation=zero(T),
        modulation_points=1,
        light_ratio=T(0.1),
        n_pix_separation=4,
        n_pix_edge=2,
        psf_centering=true,
        mode=Diffractive(),
        T=T,
    )
    stale_out = zeros(T, initial_rows, length(dm.state.coefs))
    @test_throws DimensionMismatchError interaction_matrix!(
        stale_out,
        dm,
        stale_wfs,
        PupilFunction(tel; T=T),
        src;
        amplitude=T(5e-9),
    )
    @test length(slopes(stale_wfs)) == size(imat.matrix, 1)
end

@testset "Mis-registration identification" begin
    tel = Telescope(resolution=8, diameter=8.0, central_obstruction=0.0)
    dm = DeformableMirror(tel; n_act=2, influence_width=0.4)
    wfs = PyramidWFS(tel; pupil_samples=2)
    basis = modal_basis(dm, tel; n_modes=2)
    fields = collect(Calibration.MISREG_FIELDS)
    meta, meta_fd, meta_ad = mktempdir() do root
        cd(root) do
            @test isempty(readdir())
            local_meta = Calibration.compute_meta_sensitivity_matrix(
                tel, dm, wfs, basis.M2C[:, 1:2]; n_mis_reg=2)
            local_fd = Calibration.compute_meta_sensitivity_matrix(
                tel, dm, wfs, basis.M2C[:, 1:2];
                n_mis_reg=length(fields), field_order=fields,
                sensitivity=:finite_difference)
            local_ad = Calibration.compute_meta_sensitivity_matrix(
                tel, dm, wfs, basis.M2C[:, 1:2];
                n_mis_reg=length(fields), field_order=fields)
            @test isempty(readdir())
            return local_meta, local_fd, local_ad
        end
    end
    est = Calibration.estimate_misregistration(
        meta, meta.calib0.D; misregistration_zero=Misregistration())
    @test est.shift_x ≈ 0.0
    @test est.shift_y ≈ 0.0

    @test meta_ad.field_order == fields
    @test isapprox(meta_ad.calib0.D, meta_fd.calib0.D; rtol=1e-12, atol=1e-12)
    @test isapprox(meta_ad.meta.D, meta_fd.meta.D; rtol=2e-3, atol=1e-9)

    mktempdir() do root
        cache_path = joinpath(root, "meta-sensitivity.bin")
        @test_throws MethodError Calibration.compute_meta_sensitivity_matrix(
            tel, dm, wfs, basis.M2C[:, 1:2]; n_mis_reg=2,
            cache_path=cache_path)
        @test_throws MethodError Calibration.SPRINT(
            tel, dm, wfs, basis.M2C[:, 1:2]; n_mis_reg=2,
            save_sensitivity=false)
        @test_throws MethodError Calibration.SPRINT(
            tel, dm, wfs, basis.M2C[:, 1:2]; n_mis_reg=2,
            recompute_sensitivity=true)
        @test !ispath(cache_path)
        @test isempty(readdir(root))
    end

    sampled_topology = SampledActuatorTopology(actuator_coordinates(dm)[:, 1:2])
    measured_dm = DeformableMirror(tel; topology=sampled_topology,
        influence_model=MeasuredInfluenceFunctions(Array(dm.state.modes[:, 1:2])))
    @test_throws UnsupportedAlgorithm Calibration.compute_meta_sensitivity_matrix(
        tel, measured_dm, wfs, basis.M2C[:, 1:2]; n_mis_reg=2)
    @test_throws UnsupportedAlgorithm Calibration.compute_meta_sensitivity_matrix(
        tel, dm, wfs, basis.M2C[:, 1:2]; n_mis_reg=2, wfs_mis_registered=true)
    wfs_meta = Calibration.compute_meta_sensitivity_matrix(
        tel, dm, wfs, basis.M2C[:, 1:2]; n_mis_reg=2, wfs_mis_registered=true,
        sensitivity=:finite_difference)
    @test wfs_meta.field_order == [:shift_x, :shift_y]
    assert_meta_sensitivity_contract(wfs_meta, 2)
end

@testset "Calibration workflow contracts" begin
    tel = Telescope(resolution=8, diameter=8.0, central_obstruction=0.0)
    pupil = PupilFunction(tel)
    src = Source(band=:I, magnitude=0.0)
    atm = KolmogorovAtmosphere(tel; r0=0.2,
        reference_wavelength_m=TEST_ATMOSPHERE_REFERENCE_WAVELENGTH_M,
        L0=25.0)
    dm = DeformableMirror(tel; n_act=2, influence_width=0.4)
    wfs = PyramidWFS(tel; pupil_samples=2)
    det = Detector(noise=NoiseNone(), exposure_duration=1.0, qe=1.0, binning=1)

    basis = modal_basis(dm, tel; n_modes=2)
    assert_modal_basis_contract(basis, length(dm.state.coefs), 2)

    imat = interaction_matrix(dm, wfs, pupil; amplitude=0.1)
    assert_interaction_matrix_contract(imat, length(slopes(wfs)), length(dm.state.coefs), 0.1)

    imat_basis = interaction_matrix(dm, wfs, pupil, basis.M2C;
        amplitude=0.1)
    assert_interaction_matrix_contract(imat_basis, length(slopes(wfs)), size(basis.M2C, 2), 0.1)

    aliased_commands = reshape(dm.state.coefs, :, 1)
    coefficients_before_alias_rejection = copy(dm.state.coefs)
    @test_throws InvalidConfiguration interaction_matrix(
        dm, wfs, pupil, aliased_commands; amplitude=0.1)
    @test dm.state.coefs == coefficients_before_alias_rejection

    control_matrix = ControlMatrix(imat.matrix)
    assert_control_matrix_contract(control_matrix, imat.matrix)
    noninverted_control_matrix = ControlMatrix(imat.matrix; invert=false)
    assert_control_matrix_contract(noninverted_control_matrix, imat.matrix; inverted=false)
    truncated_control_matrix = with_truncation(control_matrix, 0)
    assert_control_matrix_contract(truncated_control_matrix, imat.matrix)
    @test truncated_control_matrix.n_trunc == 0
    calib = ao_calibration(tel, dm, wfs; n_modes=2, amplitude=0.1)
    assert_ao_calibration_contract(calib, length(dm.state.coefs), 2)

    calib = ao_calibration(tel, dm, wfs; n_modes=2, amplitude=0.1, basis=basis)
    assert_ao_calibration_contract(calib, length(dm.state.coefs), 2)
    @test calib.calibration.D == imat_basis.matrix

    meta = Calibration.compute_meta_sensitivity_matrix(
        tel, dm, wfs, basis.M2C[:, 1:2]; n_mis_reg=2)
    assert_meta_sensitivity_contract(meta, 2)

    sprint = Calibration.SPRINT(
        tel, dm, wfs, basis.M2C[:, 1:2]; n_mis_reg=2)
    @test sprint.meta isa Calibration.MetaSensitivity
    @test !hasfield(typeof(sprint), :cache_path)
    @test !hasfield(typeof(sprint), :save_sensitivity)
    @test !hasfield(typeof(sprint), :recompute_sensitivity)
    est = Calibration.estimate!(sprint, meta.calib0.D)
    @test est isa Misregistration
    mktempdir() do root
        cd(root) do
            refreshed = Calibration.estimate!(sprint, meta.calib0.D;
                n_update_zero_point=1, tel=tel, dm=dm, wfs=wfs,
                basis=basis.M2C[:, 1:2])
            @test refreshed isa Misregistration
            @test sprint.meta isa Calibration.MetaSensitivity
            @test isempty(readdir())
        end
    end

    diversity = fill(eltype(pupil.opd)(1e-9), size(pupil.opd))
    lift_basis = basis_from_m2c(dm, tel, basis.M2C)
    lift_forward = prepare_lift_forward_model(tel, src, lift_basis,
        diversity;
        diversity_opd=diversity, focal_resolution=8)
    psf_in = reference_direct_image(tel, src; zero_padding=1)
    lift_observation = LiFTObservation(lift_forward, copy(psf_in))
    coeffs = zeros(eltype(psf_in), 2)
    lift = prepare_lift_estimator(LiFT(iterations=2, mode_ids=1:2,
            jacobian_method=LiFTNumericalJacobian(),
            check_convergence=false), lift_forward, lift_observation, coeffs)
    WavefrontSensors.reconstruct!(lift)
    @test length(coeffs) == 2
    @test diagnostics(lift).residual_norm >= 0
end
