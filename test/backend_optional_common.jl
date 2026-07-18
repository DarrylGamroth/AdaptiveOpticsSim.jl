backend_package_name(::Type{AdaptiveOpticsSim.CUDABackendTag}) = "CUDA"
backend_package_name(::Type{AdaptiveOpticsSim.AMDGPUBackendTag}) = "AMDGPU"

backend_label(::Type{AdaptiveOpticsSim.CUDABackendTag}) = "CUDA"
backend_label(::Type{AdaptiveOpticsSim.AMDGPUBackendTag}) = "AMDGPU"

backend_full_smoke_env(::Type{AdaptiveOpticsSim.CUDABackendTag}) = "ADAPTIVEOPTICS_TEST_FULL_CUDA"
backend_full_smoke_env(::Type{AdaptiveOpticsSim.AMDGPUBackendTag}) = "ADAPTIVEOPTICS_TEST_FULL_AMDGPU"

if !isdefined(@__MODULE__, :ContractRateModel)
    include(joinpath(@__DIR__, "wfs_stage_contract_fixtures.jl"))
end

backend_selector(::Type{AdaptiveOpticsSim.CUDABackendTag}) = AdaptiveOpticsSim.CUDABackend()
backend_selector(::Type{AdaptiveOpticsSim.AMDGPUBackendTag}) = AdaptiveOpticsSim.AMDGPUBackend()

struct OptionalStaticAtmosphere{A,B<:AbstractArrayBackend} <: AdaptiveOpticsSim.AbstractAtmosphere
    screen::A
end

AdaptiveOpticsSim.backend(::OptionalStaticAtmosphere{<:Any,B}) where {B} = B()
AdaptiveOpticsSim.advance!(atm::OptionalStaticAtmosphere, tel::Telescope, rng::AbstractRNG) = atm
AdaptiveOpticsSim.advance!(atm::OptionalStaticAtmosphere, tel::Telescope; rng::AbstractRNG=Random.default_rng()) = atm

function optional_detector_calibration_signature_allocation_bytes(
    det, seed::UInt)
    AdaptiveOpticsSim.detector_calibration_signature(det, seed)
    return @allocated AdaptiveOpticsSim.detector_calibration_signature(
        det, seed)
end

function AdaptiveOpticsSim.propagate!(atm::OptionalStaticAtmosphere, tel::Telescope)
    copyto!(tel.state.opd, atm.screen)
    return tel
end

function OptionalStaticAtmosphere(tel::Telescope; T::Type{<:AbstractFloat}=Float32, backend::AbstractArrayBackend=backend(tel))
    selector = AdaptiveOpticsSim.require_same_backend(tel, AdaptiveOpticsSim._resolve_backend_selector(backend))
    array_backend = AdaptiveOpticsSim._resolve_array_backend(selector)
    host = zeros(T, tel.params.resolution, tel.params.resolution)
    host .*= Array(pupil_mask(tel))
    screen = array_backend{T}(undef, size(host)...)
    copyto!(screen, host)
    return OptionalStaticAtmosphere{typeof(screen),typeof(selector)}(screen)
end

function run_optional_backend_selector_smoke(::Type{B}, BackendArray) where {B<:AdaptiveOpticsSim.GPUBackendTag}
    selector = backend_selector(B)
    array_backend = AdaptiveOpticsSim._resolve_array_backend(selector)
    T = Float32
    tel = Telescope(resolution=8, diameter=T(1), central_obstruction=T(0),
        T=T, backend=selector)
    dm = DeformableMirror(tel; n_act=2, influence_width=T(0.3), T=T, backend=selector)
    dm_dense = DeformableMirror(tel; n_act=2, influence_model=DenseInfluenceMatrix(Array(dm.state.modes)), T=T, backend=selector)
    zernike_modal = ModalControllableOptic(tel, ZernikeOpticBasis([2, 3]); T=T, backend=selector)
    cartesian_modal = ModalControllableOptic(tel, CartesianTiltBasis(; scale=T(0.1)); T=T, backend=selector)
    wfs = ShackHartmannWFS(tel; n_lenslets=2, mode=Diffractive(), T=T, backend=selector)
    det = Detector(noise=NoiseNone(), integration_time=T(1), qe=T(1), binning=1, T=T, backend=selector)
    sampled_response = SampledFrameResponse(
        T[0 0.1 0; 0.1 0.6 0.1; 0 0.1 0]; T=T,
        backend=selector)
    sampled_qe = AdaptiveOpticsSim.SampledQuantumEfficiency(
        T[0.5e-6, 0.6e-6, 0.7e-6], T[0.2, 0.8, 0.4]; T=T)
    calibration_det = Detector(noise=NoiseNone(), integration_time=T(1),
        qe=sampled_qe, response_model=sampled_response, T=T,
        backend=selector)
    @test tel.state.opd isa BackendArray
    @test dm.state.coefs isa BackendArray
    @test dm.state.modes isa AdaptiveOpticsSim.GaussianInfluenceOperator
    @test typeof(backend(dm.state.modes)) === typeof(selector)
    @test similar(dm.state.modes, T, 2, 2) isa BackendArray
    @test AdaptiveOpticsSim.materialize_influence_matrix(dm) isa BackendArray
    @test dm_dense.state.coefs isa BackendArray
    @test dm_dense.state.modes isa BackendArray
    @test Array(dm_dense.state.modes) ≈ Array(dm.state.modes) atol=0 rtol=0
    @test zernike_modal.state.coefs isa BackendArray
    @test zernike_modal.state.modes isa BackendArray
    @test cartesian_modal.state.coefs isa BackendArray
    @test cartesian_modal.state.modes isa BackendArray
    @test slopes(wfs) isa BackendArray
    @test det.state.frame isa BackendArray
    @test calibration_det.params.response_model.kernel isa BackendArray
    seed = UInt(0x51a7)
    AdaptiveOpticsSim.detector_calibration_signature(calibration_det, seed)
    optional_detector_calibration_signature_allocation_bytes(
        calibration_det, seed)
    @test optional_detector_calibration_signature_allocation_bytes(
        calibration_det, seed) == 0

    gain_map = reshape(T.(range(T(0.6), T(1.2); length=16)), 4, 4)
    bad_mask = falses(4, 4)
    bad_mask[2, 3] = true
    defect_det = Detector(noise=NoiseNone(), sensor=CMOSSensor(T=T,
            backend=selector),
        defect_model=CompositeDetectorDefectModel(
            PixelResponseNonuniformity(gain_map; T=T, backend=selector),
            BadPixelMask(bad_mask; T=T, backend=selector)),
        T=T, backend=selector)
    photon_rate = array_backend(fill(T(2), 4, 4))
    calibration_frame = AdaptiveOpticsSim.detector_calibration_frame!(
        defect_det, photon_rate, one(T))
    expected_frame = T(2) .* gain_map
    expected_frame[2, 3] = zero(T)
    @test Array(calibration_frame) ≈ expected_frame rtol=T(1e-6)
    @test optional_detector_calibration_signature_allocation_bytes(
        defect_det, seed) == 0

    invalid_values = array_backend(T[-1 1; 1 1])
    invalid_metadata = OpticalPlaneMetadata(FocalPlane(), invalid_values;
        coordinate_domain=AngularCoordinates(), sampling=(one(T), one(T)),
        normalization=PhotonRateNormalization(),
        spatial_measure=CellIntegratedMeasure(),
        coherence=IncoherentIntensityAddition())
    invalid_map = IntensityMap(invalid_metadata, invalid_values)
    @test_throws InvalidConfiguration prepare_detector_acquisition(
        det, invalid_map)
    return nothing
end

function run_optional_lgs_convolution_normalization(::Type{B}) where {B<:AdaptiveOpticsSim.GPUBackendTag}
    selector = backend_selector(B)
    array_backend = AdaptiveOpticsSim._resolve_array_backend(selector)
    T = Float32
    n = 8
    expected = reshape(collect(range(T(0.25), T(2); length=n * n)), n, n)
    expected_stack = cat(expected, reverse(expected; dims=2); dims=3)
    intensity_stack = array_backend(copy(expected_stack))
    kernel_fft = array_backend(ones(Complex{T}, n, n, 2))
    fft_stack = array_backend(zeros(Complex{T}, n, n, 2))
    fft_plan = AdaptiveOpticsSim.plan_fft_backend!(fft_stack, (1, 2))
    ifft_plan = AdaptiveOpticsSim.plan_ifft_backend!(fft_stack, (1, 2))

    AdaptiveOpticsSim.apply_lgs_convolution_stack!(
        intensity_stack, kernel_fft, fft_stack, fft_plan, ifft_plan)
    actual = Array(intensity_stack)
    @test actual ≈ expected_stack rtol=T(1e-5) atol=T(1e-6)
    @test vec(sum(actual; dims=(1, 2))) ≈
          vec(sum(expected_stack; dims=(1, 2))) rtol=T(1e-5)
    return nothing
end

function run_optional_sodium_profile_wfs(::Type{B},
    BackendArray) where {B<:AdaptiveOpticsSim.GPUBackendTag}
    selector = backend_selector(B)
    T = Float32
    tel = Telescope(resolution=16, diameter=T(8),
        central_obstruction=zero(T), T=T, backend=selector)

    for family in (:pyramid, :bioedge)
        src = LGSSource(
            na_profile=T[80000 90000 100000; 0.2 0.6 0.2],
            laser_coordinates=(T(1), T(-0.5)),
            fwhm_spot_up=T(0.8),
            photon_irradiance=one(T),
            T=T,
        )
        wfs = family === :pyramid ?
            PyramidWFS(tel; pupil_samples=4, mode=Diffractive(),
                modulation=zero(T), T=T, backend=selector) :
            BioEdgeWFS(tel; pupil_samples=4, mode=Diffractive(),
                modulation=zero(T), T=T, backend=selector)

        AdaptiveOpticsSim.ensure_lgs_kernel!(wfs, tel, src)
        propagation = wfs.front_end.propagation
        @test propagation.lgs_kernel_fft isa BackendArray
        original_tag = propagation.lgs_kernel_tag
        original_kernel = Array(propagation.lgs_kernel_fft)
        @test all(isfinite, original_kernel)
        slopes = measure!(wfs, tel, src)
        AdaptiveOpticsSim.synchronize_backend!(
            AdaptiveOpticsSim.execution_style(slopes))
        @test slopes isa BackendArray
        @test all(isfinite, Array(slopes))

        src.params.na_profile[2, :] .= T[0.8, 0.1, 0.1]
        AdaptiveOpticsSim.ensure_lgs_kernel!(wfs, tel, src)
        @test propagation.lgs_kernel_tag != original_tag
        @test !isapprox(Array(propagation.lgs_kernel_fft), original_kernel;
            rtol=T(1e-5), atol=T(1e-6))
    end
    return nothing
end

function run_optional_zernike_normalization(
    ::Type{B}, BackendArray) where {B<:AdaptiveOpticsSim.GPUBackendTag}
    selector = backend_selector(B)
    T = Float32
    cpu_tel = Telescope(resolution=16, diameter=T(8),
        central_obstruction=zero(T), T=T, backend=CPUBackend())
    gpu_tel = Telescope(resolution=16, diameter=T(8),
        central_obstruction=zero(T), T=T, backend=selector)
    src = Source(band=:custom, wavelength=T(0.75e-6),
        photon_irradiance=T(10), T=T)
    normalization_scale = T(0.375)
    frame_host = reshape(collect(range(T(0.25), T(2); length=16)),
        4, 4)

    for normalization in (MeanValidFluxNormalization(),
            IncidenceFluxNormalization())
        cpu_wfs = ZernikeWFS(cpu_tel; pupil_samples=8, binning=2,
            normalization=normalization, T=T, backend=CPUBackend())
        gpu_wfs = ZernikeWFS(gpu_tel; pupil_samples=8, binning=2,
            normalization=normalization, T=T, backend=selector)
        fill!(cpu_wfs.state.reference_signal_2d, zero(T))
        fill!(gpu_wfs.state.reference_signal_2d, zero(T))
        frame = BackendArray(copy(frame_host))

        expected_normalization = AdaptiveOpticsSim.zernike_normalization(
            normalization, cpu_wfs, cpu_tel, src, frame_host,
            normalization_scale)
        actual_normalization = AdaptiveOpticsSim.zernike_normalization(
            normalization, gpu_wfs, gpu_tel, src, frame,
            normalization_scale)
        @test actual_normalization ≈ expected_normalization rtol=T(2e-5)

        expected = copy(AdaptiveOpticsSim.zernike_signal!(cpu_wfs,
            cpu_tel, frame_host, src, normalization_scale))
        actual = AdaptiveOpticsSim.zernike_signal!(gpu_wfs, gpu_tel,
            frame, src, normalization_scale)
        AdaptiveOpticsSim.synchronize_backend!(
            AdaptiveOpticsSim.execution_style(actual))
        @test gpu_wfs.state.normalization_sum isa BackendArray
        @test actual isa BackendArray
        @test Array(actual) ≈ expected rtol=T(2e-5) atol=T(2e-6)
    end

    zero_src = Source(band=:custom, wavelength=wavelength(src),
        photon_irradiance=zero(T), T=T)
    zero_wfs = ZernikeWFS(gpu_tel; pupil_samples=8, binning=2,
        normalization=IncidenceFluxNormalization(), T=T,
        backend=selector)
    fill!(zero_wfs.state.reference_signal_2d, zero(T))
    zero_slopes = AdaptiveOpticsSim.zernike_signal!(zero_wfs, gpu_tel,
        BackendArray(copy(frame_host)), zero_src, one(T))
    AdaptiveOpticsSim.synchronize_backend!(
        AdaptiveOpticsSim.execution_style(zero_slopes))
    @test all(iszero, Array(zero_slopes))
    @test all(isfinite, Array(zero_slopes))
    return nothing
end

function run_optional_wfs_stage_contracts(
    ::Type{B}, BackendArray) where {B<:AdaptiveOpticsSim.GPUBackendTag}
    selector = backend_selector(B)
    T = Float32
    tel = Telescope(resolution=4, diameter=T(2),
        central_obstruction=zero(T), T=T, backend=selector)
    src = Source(band=:custom, wavelength=T(0.75e-6),
        photon_irradiance=T(4), T=T)
    pupil = PupilFunction(tel; T=T, backend=selector)
    fill!(pupil.opd, T(2e-9))
    field = ElectricField(pupil, src; zero_padding=1, T=T)
    field_formation = prepare_pupil_field(tel, pupil, src, field;
        center_even_grid=false)
    fill_electric_field!(field, pupil, field_formation)

    rate_values = BackendArray(zeros(T, 4, 4))
    rate_metadata = OpticalPlaneMetadata(DetectorPlane(), rate_values;
        coordinate_domain=AngularCoordinates(),
        sampling=(one(T), one(T)),
        spectral=MonochromaticChannel(T(wavelength(src))),
        normalization=PhotonRateNormalization(),
        spatial_measure=CellIntegratedMeasure(),
        coherence=IncoherentIntensityAddition())
    rate = IntensityMap(rate_metadata, rate_values)
    optical_model = ContractRateModel(T(3), T(1e6), rate.metadata)
    optical_plan = prepare_wfs_optical_formation(optical_model, pupil, rate)

    detector = Detector(noise=NoiseNone(), integration_time=T(0.4),
        qe=T(0.5), response_model=NullFrameResponse(), T=T,
        backend=selector)
    binding = contract_detector_binding(detector, rate)
    observation = binding.observation
    acquisition_plan = prepare_wfs_acquisition(
        ContractDetectorAcquisitionModel((binding,)), rate, observation)
    measurement_values = similar(observation.storage)
    measurement = WFSMeasurement(measurement_values;
        units=:detector_signal, kind=:copied_detector_frame)
    estimator_plan = prepare_wfs_estimation(ContractCopyEstimator(
        :detector_signal, :copied_detector_frame),
        observation, measurement)
    rng = Xoshiro(0x574653)

    @test @inferred(form_wfs_optical_products!(rate, pupil,
        optical_plan)) === rate
    @test @inferred(acquire_wfs_observation!(observation, rate,
        acquisition_plan, rng)) === observation
    @test @inferred(estimate_wfs_measurement!(measurement, observation,
        estimator_plan)) === measurement
    AdaptiveOpticsSim.synchronize_backend!(
        AdaptiveOpticsSim.execution_style(measurement.storage))
    @test rate.values isa BackendArray
    @test observation.storage isa BackendArray
    @test measurement.storage isa BackendArray
    @test observation.metadata.device == plane_device(observation.storage)
    @test measurement.metadata.device == plane_device(measurement.storage)
    @test Array(measurement.storage) == Array(observation.storage)
    @test sum(Array(observation.storage)) ≈
        sum(Array(rate.values)) * T(0.4) * T(0.5) rtol=T(2e-6)

    direct_values = similar(rate.values)
    direct_measurement = WFSMeasurement(direct_values;
        units=:field_intensity, kind=:direct_field_measurement)
    direct_plan = prepare_wfs_estimation(ContractDirectCopyEstimator(
        :field_intensity, :direct_field_measurement), field,
        direct_measurement)
    @test wfs_measurement_path(direct_plan) isa DirectMeasurementPath
    @test @inferred(estimate_wfs_measurement!(direct_measurement, field,
        direct_plan)) === direct_measurement
    AdaptiveOpticsSim.synchronize_backend!(
        AdaptiveOpticsSim.execution_style(direct_measurement.storage))
    @test direct_measurement.storage isa BackendArray
    @test Array(direct_measurement.storage) ≈ abs2.(Array(field.values))

    second_values = BackendArray(zeros(T, 4, 4))
    second_metadata = OpticalPlaneMetadata(DetectorPlane(), second_values;
        coordinate_domain=AngularCoordinates(),
        sampling=(T(0.5), T(0.5)),
        spectral=MonochromaticChannel(T(0.9e-6)),
        normalization=PhotonRateNormalization(),
        spatial_measure=CellIntegratedMeasure(),
        coherence=IncoherentIntensityAddition())
    second_rate = IntensityMap(second_metadata, second_values)
    bundle = OpticalProductBundle(rate, second_rate)
    bundle_model = ContractBundleRateModel((
        ContractRateModel(T(2), T(1e6), rate.metadata),
        ContractRateModel(T(5), T(1e6), second_rate.metadata),
    ))
    bundle_plan = prepare_wfs_optical_formation(bundle_model, pupil, bundle)
    form_wfs_optical_products!(bundle, pupil, bundle_plan)
    second_detector = Detector(noise=NoiseNone(), integration_time=T(0.7),
        qe=T(0.25), response_model=NullFrameResponse(), T=T,
        backend=selector)
    first_binding = contract_detector_binding(detector, rate)
    second_binding = contract_detector_binding(second_detector, second_rate)
    observations = (first_binding.observation, second_binding.observation)
    multi_acquisition = prepare_wfs_acquisition(
        ContractDetectorAcquisitionModel((first_binding, second_binding)),
        bundle, observations)
    @test @inferred(acquire_wfs_observation!(observations, bundle,
        multi_acquisition, (Xoshiro(1), Xoshiro(2)))) === observations
    AdaptiveOpticsSim.synchronize_backend!(
        AdaptiveOpticsSim.execution_style(second_binding.observation.storage))
    @test first_binding.observation.storage isa BackendArray
    @test second_binding.observation.storage isa BackendArray
    @test first_binding.observation.metadata.device ==
        second_binding.observation.metadata.device

    packed_values = BackendArray(zeros(T, 8, 4))
    regions = ((1:4, 1:4), (5:8, 1:4))
    packed = WFSObservation(packed_values; units=:detector_signal,
        layout=regions)
    packed_plan = prepare_wfs_acquisition(
        ContractPackedAcquisition(regions, T(0.2)), bundle, packed)
    @test @inferred(acquire_wfs_observation!(packed, bundle, packed_plan,
        Xoshiro(3))) === packed
    AdaptiveOpticsSim.synchronize_backend!(
        AdaptiveOpticsSim.execution_style(packed.storage))
    @test packed.storage isa BackendArray
    packed_host = Array(packed.storage)
    @test packed_host[1:4, :] ≈ Array(rate.values) .* T(0.2)
    @test packed_host[5:8, :] ≈ Array(second_rate.values) .* T(0.2)

    physical_wfs = ShackHartmannWFS(tel; n_lenslets=2,
        n_pix_subap=2, mode=Diffractive(), T=T, backend=selector)
    physical_rate = shack_hartmann_rate_map(physical_wfs, pupil, src)
    physical_front_end = ShackHartmannOpticalFrontEnd(physical_wfs, src)
    @test !hasfield(typeof(physical_front_end), :sensor)
    @test physical_front_end.propagation.fft_stack isa BackendArray
    physical_optical_plan = prepare_wfs_optical_formation(
        physical_front_end, pupil, physical_rate)
    form_wfs_optical_products!(physical_rate, pupil, physical_optical_plan)
    @test physical_rate.values isa BackendArray
    @test all(isfinite, Array(physical_rate.values))
    @test sum(Array(physical_rate.values)) > zero(T)

    physical_field_wfs = ShackHartmannWFS(tel; n_lenslets=2,
        n_pix_subap=2, mode=Diffractive(), T=T, backend=selector)
    physical_field_rate = shack_hartmann_rate_map(physical_field_wfs, field)
    physical_field_plan = prepare_wfs_optical_formation(
        ShackHartmannOpticalFrontEnd(physical_field_wfs), field,
        physical_field_rate)
    form_wfs_optical_products!(physical_field_rate, field,
        physical_field_plan)
    AdaptiveOpticsSim.synchronize_backend!(
        AdaptiveOpticsSim.execution_style(physical_field_rate.values))
    @test physical_field_rate.values isa BackendArray
    @test isapprox(Array(physical_field_rate.values),
        Array(physical_rate.values); rtol=T(2e-5), atol=T(2e-5))

    physical_asterism = Asterism([
        Source(band=:custom, wavelength=wavelength(src),
            photon_irradiance=T(1), T=T),
        Source(band=:custom, wavelength=wavelength(src),
            photon_irradiance=T(3), T=T),
    ])
    physical_asterism_wfs = ShackHartmannWFS(tel; n_lenslets=2,
        n_pix_subap=2, mode=Diffractive(), T=T, backend=selector)
    physical_asterism_rate = shack_hartmann_rate_map(
        physical_asterism_wfs, pupil, physical_asterism)
    physical_asterism_plan = prepare_wfs_optical_formation(
        ShackHartmannOpticalFrontEnd(physical_asterism_wfs,
            physical_asterism), pupil, physical_asterism_rate)
    form_wfs_optical_products!(physical_asterism_rate, pupil,
        physical_asterism_plan)
    AdaptiveOpticsSim.synchronize_backend!(
        AdaptiveOpticsSim.execution_style(physical_asterism_rate.values))
    @test physical_asterism_rate.values isa BackendArray
    @test isapprox(Array(physical_asterism_rate.values),
        Array(physical_rate.values); rtol=T(2e-5), atol=T(2e-5))

    physical_detector = Detector(noise=NoiseNone(), integration_time=T(0.4),
        qe=T(0.5), response_model=NullFrameResponse(), T=T,
        backend=selector)
    physical_observation = WFSObservation(similar(physical_rate.values);
        units=:electron_count, layout=:lenslet_mosaic)
    physical_acquisition = prepare_wfs_acquisition(physical_detector,
        physical_rate, physical_observation)
    acquire_wfs_observation!(physical_observation, physical_rate,
        physical_acquisition, Xoshiro(0x5348))
    @test physical_observation.storage isa BackendArray
    @test sum(Array(physical_observation.storage)) ≈
        sum(Array(physical_rate.values)) * T(0.4) * T(0.5) rtol=T(2e-5)

    cpu_tel = Telescope(resolution=4, diameter=T(2),
        central_obstruction=zero(T), T=T)
    cpu_wfs = ShackHartmannWFS(cpu_tel; n_lenslets=2,
        n_pix_subap=2, mode=Diffractive(), T=T)
    mixed_backend_rate = shack_hartmann_rate_map(cpu_wfs, pupil, src)
    mixed_optical_error = try
        prepare_wfs_optical_formation(
            ShackHartmannOpticalFrontEnd(cpu_wfs, src), pupil,
            mixed_backend_rate)
        nothing
    catch err
        err
    end
    @test mixed_optical_error isa WFSPreparationError
    @test mixed_optical_error.reason === :backend

    cpu_observation = WFSObservation(zeros(T, size(physical_rate.values));
        units=:electron_count, layout=:lenslet_mosaic)
    gpu_measurement = WFSMeasurement(similar(slopes(physical_wfs));
        units=:pixel, kind=:centroid_slopes)
    @test_throws WFSPreparationError prepare_wfs_estimation(physical_wfs,
        cpu_observation, gpu_measurement)
    cpu_measurement = WFSMeasurement(zeros(T, length(slopes(physical_wfs)));
        units=:pixel, kind=:centroid_slopes)
    @test_throws WFSPreparationError prepare_wfs_estimation(physical_wfs,
        physical_observation, cpu_measurement)

    set_subaperture_calibration!(physical_wfs.calibration,
        zeros(T, size(physical_wfs.calibration.reference_signal_2d));
        centroid_response=one(T), wavelength=wavelength(src),
        signature=UInt(0x50485953))
    physical_measurement = WFSMeasurement(similar(slopes(physical_wfs));
        units=:pixel, kind=:centroid_slopes)
    physical_estimator = prepare_wfs_estimation(physical_wfs,
        physical_observation, physical_measurement)
    estimate_wfs_measurement!(physical_measurement, physical_observation,
        physical_estimator)
    AdaptiveOpticsSim.synchronize_backend!(
        AdaptiveOpticsSim.execution_style(physical_measurement.storage))
    @test physical_measurement.storage isa BackendArray
    @test all(isfinite, Array(physical_measurement.storage))

    four_pupil_cpu_tel = Telescope(resolution=4, diameter=T(2),
        central_obstruction=zero(T), T=T, backend=CPUBackend())
    four_pupil_cpu = PupilFunction(four_pupil_cpu_tel; T=T)
    copyto!(four_pupil_cpu.opd, Array(pupil.opd))
    for family in (:pyramid, :bioedge)
        cpu_sensor = family === :pyramid ?
            PyramidWFS(four_pupil_cpu_tel; pupil_samples=2,
                mode=Diffractive(), modulation=0, T=T) :
            BioEdgeWFS(four_pupil_cpu_tel; pupil_samples=2,
                mode=Diffractive(), modulation=0, T=T)
        gpu_sensor = family === :pyramid ?
            PyramidWFS(tel; pupil_samples=2, mode=Diffractive(),
                modulation=0, T=T, backend=selector) :
            BioEdgeWFS(tel; pupil_samples=2, mode=Diffractive(),
                modulation=0, T=T, backend=selector)
        cpu_front_end = family === :pyramid ?
            PyramidOpticalFrontEnd(cpu_sensor, src) :
            BioEdgeOpticalFrontEnd(cpu_sensor, src)
        gpu_front_end = family === :pyramid ?
            PyramidOpticalFrontEnd(gpu_sensor, src) :
            BioEdgeOpticalFrontEnd(gpu_sensor, src)
        cpu_rate = family === :pyramid ?
            pyramid_rate_map(cpu_front_end, four_pupil_cpu) :
            bioedge_rate_map(cpu_front_end, four_pupil_cpu)
        gpu_rate = family === :pyramid ?
            pyramid_rate_map(gpu_front_end, pupil) :
            bioedge_rate_map(gpu_front_end, pupil)
        cpu_plan = prepare_wfs_optical_formation(cpu_front_end,
            four_pupil_cpu, cpu_rate)
        gpu_plan = prepare_wfs_optical_formation(gpu_front_end, pupil,
            gpu_rate)
        form_wfs_optical_products!(cpu_rate, four_pupil_cpu, cpu_plan)
        form_wfs_optical_products!(gpu_rate, pupil, gpu_plan)
        AdaptiveOpticsSim.synchronize_backend!(
            AdaptiveOpticsSim.execution_style(gpu_rate.values))
        @test gpu_rate.values isa BackendArray
        @test isapprox(Array(gpu_rate.values), cpu_rate.values;
            rtol=T(3e-5), atol=T(3e-5))

        field_sensor = family === :pyramid ?
            PyramidWFS(tel; pupil_samples=2, mode=Diffractive(),
                modulation=0, T=T, backend=selector) :
            BioEdgeWFS(tel; pupil_samples=2, mode=Diffractive(),
                modulation=0, T=T, backend=selector)
        field_front_end = family === :pyramid ?
            PyramidOpticalFrontEnd(field_sensor) :
            BioEdgeOpticalFrontEnd(field_sensor)
        field_rate = family === :pyramid ?
            pyramid_rate_map(field_front_end, field) :
            bioedge_rate_map(field_front_end, field)
        field_plan = prepare_wfs_optical_formation(field_front_end, field,
            field_rate)
        form_wfs_optical_products!(field_rate, field, field_plan)
        AdaptiveOpticsSim.synchronize_backend!(
            AdaptiveOpticsSim.execution_style(field_rate.values))
        @test field_rate.values isa BackendArray
        @test isapprox(Array(field_rate.values), Array(gpu_rate.values);
            rtol=T(3e-5), atol=T(3e-5))

        four_pupil_detector = Detector(noise=NoiseNone(),
            integration_time=T(0.4), qe=T(0.5),
            response_model=NullFrameResponse(), T=T, backend=selector)
        four_pupil_observation = WFSObservation(similar(gpu_rate.values);
            units=:electron_count, layout=:four_pupil_mosaic)
        four_pupil_acquisition = prepare_wfs_acquisition(
            four_pupil_detector, gpu_rate, four_pupil_observation)
        acquire_wfs_observation!(four_pupil_observation, gpu_rate,
            four_pupil_acquisition, Xoshiro(0x5042))
        AdaptiveOpticsSim.synchronize_backend!(
            AdaptiveOpticsSim.execution_style(
                four_pupil_observation.storage))
        @test four_pupil_observation.storage isa BackendArray
        @test isapprox(Array(four_pupil_observation.storage),
            Array(gpu_rate.values) .* T(0.2); rtol=T(3e-5), atol=T(3e-5))

        reference = similar(gpu_sensor.estimator.state.reference_signal_2d)
        fill!(reference, zero(T))
        if family === :pyramid
            set_pyramid_calibration!(gpu_sensor, reference;
                wavelength_m=wavelength(src), signature=UInt(0x5042))
        else
            set_bioedge_calibration!(gpu_sensor, reference;
                wavelength_m=wavelength(src), signature=UInt(0x5042))
        end
        four_pupil_measurement = WFSMeasurement(similar(slopes(gpu_sensor));
            units=:dimensionless, kind=:differential_slopes)
        four_pupil_estimator = prepare_wfs_estimation(gpu_sensor,
            four_pupil_observation, four_pupil_measurement)
        estimate_wfs_measurement!(four_pupil_measurement,
            four_pupil_observation, four_pupil_estimator)
        AdaptiveOpticsSim.synchronize_backend!(
            AdaptiveOpticsSim.execution_style(
                four_pupil_measurement.storage))
        @test four_pupil_measurement.storage isa BackendArray
        @test all(isfinite, Array(four_pupil_measurement.storage))

        geometric = family === :pyramid ?
            PyramidWFS(tel; pupil_samples=2, mode=Geometric(), T=T,
                backend=selector) :
            BioEdgeWFS(tel; pupil_samples=2, mode=Geometric(), T=T,
                backend=selector)
        geometric_measurement = WFSMeasurement(similar(slopes(geometric));
            units=:metre, kind=:geometric_slopes)
        geometric_plan = prepare_wfs_estimation(geometric, pupil,
            geometric_measurement)
        estimate_wfs_measurement!(geometric_measurement, pupil,
            geometric_plan)
        AdaptiveOpticsSim.synchronize_backend!(
            AdaptiveOpticsSim.execution_style(
                geometric_measurement.storage))
        @test geometric_measurement.storage isa BackendArray
        @test all(isfinite, Array(geometric_measurement.storage))

        spectral_source = with_spectrum(src,
            SpectralBundle(T[0.7e-6, 0.9e-6], T[0.25, 0.75]; T=T))
        spectral_front_end = family === :pyramid ?
            PyramidOpticalFrontEnd(gpu_sensor, spectral_source) :
            BioEdgeOpticalFrontEnd(gpu_sensor, spectral_source)
        spectral_rates = family === :pyramid ?
            pyramid_rate_map(spectral_front_end, pupil) :
            bioedge_rate_map(spectral_front_end, pupil)
        spectral_plan = prepare_wfs_optical_formation(spectral_front_end,
            pupil, spectral_rates)
        form_wfs_optical_products!(spectral_rates, pupil, spectral_plan)
        @test all(product -> product.values isa BackendArray,
            spectral_rates)

        path_source = Asterism([
            Source(band=:custom, wavelength=wavelength(src),
                photon_irradiance=T(1), T=T),
            Source(band=:custom, wavelength=wavelength(src),
                photon_irradiance=T(3), T=T),
        ])
        second_pupil = PupilFunction(tel; T=T, backend=selector)
        copyto!(second_pupil.opd, pupil.opd)
        path_front_end = family === :pyramid ?
            PyramidOpticalFrontEnd(gpu_sensor, path_source) :
            BioEdgeOpticalFrontEnd(gpu_sensor, path_source)
        path_inputs = (pupil, second_pupil)
        path_rates = family === :pyramid ?
            pyramid_rate_map(path_front_end, path_inputs) :
            bioedge_rate_map(path_front_end, path_inputs)
        path_plan = prepare_wfs_optical_formation(path_front_end,
            path_inputs, path_rates)
        form_wfs_optical_products!(path_rates, path_inputs, path_plan)
        @test all(product -> product.values isa BackendArray, path_rates)
        @test isapprox(sum(Array(path_rates[1].values)) /
            sum(Array(path_rates[2].values)), T(1 / 3); rtol=T(3e-5))

        for lgs in (
            LGSSource(wavelength=wavelength(src),
                photon_irradiance=T(4), elongation_factor=T(1.8), T=T),
            LGSSource(wavelength=wavelength(src),
                photon_irradiance=T(4),
                na_profile=T[80_000 90_000 100_000; 0.2 0.6 0.2],
                laser_coordinates=(T(1), T(-0.5)),
                fwhm_spot_up=T(0.8), T=T),
        )
            cpu_lgs_sensor = family === :pyramid ?
                PyramidWFS(four_pupil_cpu_tel; pupil_samples=2,
                    mode=Diffractive(), modulation=0, T=T) :
                BioEdgeWFS(four_pupil_cpu_tel; pupil_samples=2,
                    mode=Diffractive(), modulation=0, T=T)
            gpu_lgs_sensor = family === :pyramid ?
                PyramidWFS(tel; pupil_samples=2, mode=Diffractive(),
                    modulation=0, T=T, backend=selector) :
                BioEdgeWFS(tel; pupil_samples=2, mode=Diffractive(),
                    modulation=0, T=T, backend=selector)
            cpu_lgs_front_end = family === :pyramid ?
                PyramidOpticalFrontEnd(cpu_lgs_sensor, lgs) :
                BioEdgeOpticalFrontEnd(cpu_lgs_sensor, lgs)
            gpu_lgs_front_end = family === :pyramid ?
                PyramidOpticalFrontEnd(gpu_lgs_sensor, lgs) :
                BioEdgeOpticalFrontEnd(gpu_lgs_sensor, lgs)
            cpu_lgs_rate = family === :pyramid ?
                pyramid_rate_map(cpu_lgs_front_end, four_pupil_cpu) :
                bioedge_rate_map(cpu_lgs_front_end, four_pupil_cpu)
            gpu_lgs_rate = family === :pyramid ?
                pyramid_rate_map(gpu_lgs_front_end, pupil) :
                bioedge_rate_map(gpu_lgs_front_end, pupil)
            cpu_lgs_plan = prepare_wfs_optical_formation(
                cpu_lgs_front_end, four_pupil_cpu, cpu_lgs_rate)
            gpu_lgs_plan = prepare_wfs_optical_formation(
                gpu_lgs_front_end, pupil, gpu_lgs_rate)
            form_wfs_optical_products!(cpu_lgs_rate, four_pupil_cpu,
                cpu_lgs_plan)
            form_wfs_optical_products!(gpu_lgs_rate, pupil, gpu_lgs_plan)
            AdaptiveOpticsSim.synchronize_backend!(
                AdaptiveOpticsSim.execution_style(gpu_lgs_rate.values))
            @test isapprox(Array(gpu_lgs_rate.values), cpu_lgs_rate.values;
                rtol=T(5e-5), atol=T(5e-5))
        end
    end

    cpu_measurement = WFSMeasurement(zeros(T, 4, 4);
        units=:detector_signal, kind=:cpu_copy)
    mismatch = try
        prepare_wfs_estimation(ContractCopyEstimator(
            :detector_signal, :cpu_copy), observation,
            cpu_measurement)
        nothing
    catch err
        err
    end
    @test mismatch isa WFSPreparationError
    @test mismatch.reason === :backend
    return nothing
end

function run_optional_plane_product_checks(tel::Telescope,
    src::AdaptiveOpticsSim.AbstractSource,
    selector::AdaptiveOpticsSim.AbstractArrayBackend, BackendArray,
    ::Type{T}) where {T<:AbstractFloat}
    wavefront = PupilFunction(tel; T=T, backend=selector)
    apply_opd!(wavefront, opd_map(tel))
    field = ElectricField(wavefront, src; zero_padding=2, T=T)
    formation = prepare_pupil_field(tel, wavefront, src, field)
    fill_electric_field!(field, wavefront, formation)
    @test wavefront.amplitude isa BackendArray
    @test wavefront.opd isa BackendArray
    @test field.values isa BackendArray
    @test field.metadata.device == plane_device(field.values)
    field_view = @view field.values[:, :]
    for wrapper in (
        field_view,
        reshape(field_view, 1, length(field_view)),
        transpose(field_view),
        PermutedDimsArray(field_view, (2, 1)),
    )
        @test plane_device(wrapper) == field.metadata.device
    end

    wrapped_intensity_parent = similar(wavefront.opd, T, 4, 4)
    fill!(wrapped_intensity_parent, one(T))
    wrapped_intensity_view = @view wrapped_intensity_parent[:, :]
    for wrapper in (
        wrapped_intensity_view,
        reshape(wrapped_intensity_view, 1, length(wrapped_intensity_view)),
        transpose(wrapped_intensity_view),
        PermutedDimsArray(wrapped_intensity_view, (2, 1)),
    )
        @test typeof(backend(wrapper)) === typeof(selector)
        @test typeof(AdaptiveOpticsSim.array_backend_selector(
            typeof(wrapper))) === typeof(selector)
        wrapper_metadata = OpticalPlaneMetadata(FocalPlane(), wrapper;
            coordinate_domain=AngularCoordinates(),
            sampling=(one(T), one(T)),
            normalization=PhotonRateNormalization(),
            spatial_measure=CellIntegratedMeasure(),
            coherence=IncoherentIntensityAddition())
        wrapper_map = IntensityMap(wrapper_metadata, wrapper)
        @test wrapper_map.values === wrapper
        @test wrapper_map.metadata.device == field.metadata.device
    end

    prepared = prepare_direct_imaging(tel, wavefront, src;
        zero_padding=2)
    formed = form_direct_image!(prepared)
    @test formed.values isa BackendArray
    @test formed.metadata.device == wavefront.metadata.device
    @test sum(Array(formed.values)) ≈
        sum(Array(pupil_photon_rate_map(tel, src))) atol=T(2e-5) rtol=T(2e-5)

    off_axis_src = Source(band=:I, magnitude=zero(T),
        coordinates=(T(0.08), T(90)), T=T)
    off_axis = prepare_direct_imaging(tel, wavefront, off_axis_src;
        zero_padding=2)
    off_axis_map = form_direct_image!(off_axis)
    @test off_axis_map.values isa BackendArray
    @test off_axis_map.metadata.device == formed.metadata.device
    @test off_axis.plan.shift_samples isa NTuple{2,Int}
    @test off_axis.plan.shift_samples != (0, 0)
    @test Array(off_axis_map.values) ≈ circshift(Array(formed.values),
        off_axis.plan.shift_samples) rtol=T(2e-5) atol=T(2e-5)

    spectral_src = with_spectrum(src,
        SpectralBundle(T[0.7e-6, 0.9e-6], T[0.4, 0.6]; T=T))
    spectral = prepare_direct_imaging(tel, wavefront, spectral_src;
        zero_padding=2)
    spectral_products = form_direct_image!(spectral)
    @test spectral_products isa OpticalProductBundle
    @test length(spectral_products) == 2
    @test all(product -> product.values isa BackendArray,
        spectral_products)
    @test all(product -> product.metadata.device == formed.metadata.device,
        spectral_products)
    @test spectral_products[1].metadata.sampling !=
        spectral_products[2].metadata.sampling

    extended_src = with_extended_source(src,
        GaussianDiskSourceModel(sigma_arcsec=T(0.02), n_side=5, T=T))
    extended = prepare_direct_imaging(tel, wavefront,
        extended_source_asterism(extended_src); zero_padding=1)
    @test extended.components isa Vector
    @test extended.products isa Vector
    @test isconcretetype(eltype(extended.components))
    @test isconcretetype(eltype(extended.products))
    @test length(extended.components) == 25
    extended_map = form_direct_image!(extended)
    @test extended_map.values isa BackendArray
    @test all(product -> product.values isa BackendArray,
        extended.products)
    @test all(product -> product.metadata.device == formed.metadata.device,
        extended.products)
    @test sum(Array(extended_map.values)) ≈
        sum(Array(pupil_photon_rate_map(tel, src))) atol=T(2e-5) rtol=T(2e-5)

    host_values = zeros(T, size(formed.values))
    host_metadata = OpticalPlaneMetadata(FocalPlane(), host_values;
        coordinate_domain=AngularCoordinates(),
        sampling=formed.metadata.sampling,
        origin=formed.metadata.origin,
        spectral=formed.metadata.spectral,
        normalization=PhotonRateNormalization(),
        spatial_measure=CellIntegratedMeasure(),
        coherence=IncoherentIntensityAddition())
    host_output = IntensityMap(host_metadata, host_values)
    @test_throws InvalidConfiguration prepare_direct_imaging(tel,
        wavefront, src, prepared.field, host_output)

    sum_output_values = similar(formed.values)
    first_sum_values = similar(formed.values)
    second_sum_values = similar(formed.values)
    fill!(first_sum_values, one(T))
    fill!(second_sum_values, T(2))
    sum_output = IntensityMap(formed.metadata, sum_output_values)
    first_sum_input = IntensityMap(formed.metadata, first_sum_values)
    second_sum_input = IntensityMap(formed.metadata, second_sum_values)
    sum_plan = prepare_incoherent_sum(sum_output, first_sum_input,
        second_sum_input)
    accumulate_intensity!(sum_output,
        (first_sum_input, second_sum_input), sum_plan)
    AdaptiveOpticsSim.synchronize_backend!(
        AdaptiveOpticsSim.execution_style(sum_output.values))
    @test sum_output.values isa BackendArray
    @test Array(sum_output.values) == fill(T(3), size(sum_output.values))

    detector = Detector(integration_time=T(0.5), noise=NoiseNone(),
        qe=T(0.5), response_model=NullFrameResponse(), T=T,
        backend=selector)
    acquisition = prepare_detector_acquisition(detector, prepared.output)
    detector_frame = capture!(detector, prepared.output, acquisition;
        rng=MersenneTwister(301))
    AdaptiveOpticsSim.synchronize_backend!(
        AdaptiveOpticsSim.execution_style(detector_frame))
    @test detector_frame isa BackendArray
    @test sum(Array(detector_frame)) ≈
        sum(Array(prepared.output.values)) * T(0.25) atol=T(2e-5) rtol=T(2e-5)
    identical_detector = Detector(integration_time=T(0.5), noise=NoiseNone(),
        qe=T(0.5), response_model=NullFrameResponse(), T=T,
        backend=selector)
    @test identical_detector.params === detector.params
    @test identical_detector.state !== detector.state
    @test_throws InvalidConfiguration capture!(identical_detector,
        prepared.output, acquisition; rng=MersenneTwister(301))

    long_detector = Detector(integration_time=T(1.0), noise=NoiseNone(),
        qe=T(0.5), response_model=NullFrameResponse(), T=T,
        backend=selector)
    long_acquisition = prepare_detector_acquisition(long_detector,
        prepared.output)
    @test acquisition.input_values === prepared.output.values
    @test long_acquisition.input_values === prepared.output.values
    short_snapshot = copy(Array(detector_frame))
    long_frame = capture!(long_detector, prepared.output,
        long_acquisition; rng=MersenneTwister(307))
    AdaptiveOpticsSim.synchronize_backend!(
        AdaptiveOpticsSim.execution_style(long_frame))
    @test Array(long_frame) ≈ 2 .* short_snapshot atol=T(2e-5) rtol=T(2e-5)

    incremental_detector = Detector(integration_time=T(0.5),
        noise=NoiseNone(), qe=T(0.5), response_model=NullFrameResponse(),
        T=T, backend=selector)
    incremental_acquisition = prepare_detector_acquisition(
        incremental_detector, prepared.output)
    capture!(incremental_detector, prepared.output,
        incremental_acquisition; rng=MersenneTwister(305),
        sample_time=T(0.2))
    @test !readout_ready(incremental_detector)
    incremental_frame = capture!(incremental_detector, prepared.output,
        incremental_acquisition; rng=MersenneTwister(306),
        sample_time=T(0.3))
    AdaptiveOpticsSim.synchronize_backend!(
        AdaptiveOpticsSim.execution_style(incremental_frame))
    @test incremental_frame isa BackendArray
    @test readout_ready(incremental_detector)
    @test sum(Array(incremental_frame)) ≈
        sum(Array(prepared.output.values)) * T(0.25) atol=T(2e-5) rtol=T(2e-5)
    cpu_detector = Detector(integration_time=T(0.5), noise=NoiseNone(),
        qe=T(0.5), response_model=NullFrameResponse(), T=T,
        backend=CPUBackend())
    @test_throws InvalidConfiguration prepare_detector_acquisition(
        cpu_detector, prepared.output)
    @test_throws InvalidConfiguration capture!(cpu_detector,
        prepared.output, acquisition; rng=MersenneTwister(306))

    density_host = zeros(T, 9, 9)
    density_host[3, 5] = T(8)
    density_values = BackendArray(density_host)
    density_metadata = OpticalPlaneMetadata(FocalPlane(), density_values;
        coordinate_domain=AngularCoordinates(), sampling=(T(0.5), T(0.25)),
        spectral=MonochromaticChannel(T(wavelength(src))),
        normalization=AdaptiveOpticsSim.PhotonRateNormalization(),
        spatial_measure=AdaptiveOpticsSim.SpatialDensityMeasure(),
        coherence=AdaptiveOpticsSim.IncoherentIntensityAddition())
    density_map = IntensityMap(density_metadata, density_values)
    response_kernel = T[0 0.1 0; 0.1 0.6 0.1; 0 0.1 0]
    response_detector = Detector(integration_time=T(2), noise=NoiseNone(),
        qe=T(0.5), binning=3,
        response_model=SampledFrameResponse(response_kernel; T=T), T=T,
        backend=selector)
    response_acquisition = prepare_detector_acquisition(response_detector,
        density_map)
    response_frame = capture!(response_detector, density_map,
        response_acquisition; rng=MersenneTwister(302))
    AdaptiveOpticsSim.synchronize_backend!(
        AdaptiveOpticsSim.execution_style(response_frame))
    manual_response = zeros(T, 9, 9)
    manual_response[3, 5] = T(4.8)
    manual_response[2, 5] = T(0.8)
    manual_response[4, 5] = T(0.8)
    manual_response[3, 4] = T(0.8)
    manual_response[3, 6] = T(0.8)
    manual_binned = zeros(T, 3, 3)
    AdaptiveOpticsSim.bin2d!(manual_binned, manual_response, 3)
    @test Array(response_frame) ≈ manual_binned .* T(0.125) atol=T(1e-6) rtol=T(1e-6)

    asymmetric_kernel = T[0 0 0; 0.1 0.2 0.7; 0 0 0]
    edge_host = zeros(T, 9, 9)
    edge_host[5, end] = one(T)
    edge_values = BackendArray(edge_host)
    edge_metadata = OpticalPlaneMetadata(FocalPlane(), edge_values;
        coordinate_domain=AngularCoordinates(), sampling=(one(T), one(T)),
        spectral=MonochromaticChannel(T(wavelength(src))),
        normalization=AdaptiveOpticsSim.PhotonRateNormalization(),
        spatial_measure=AdaptiveOpticsSim.CellIntegratedMeasure(),
        coherence=AdaptiveOpticsSim.IncoherentIntensityAddition())
    edge_map = IntensityMap(edge_metadata, edge_values)
    edge_detector = Detector(integration_time=one(T), noise=NoiseNone(),
        qe=one(T), response_model=SampledFrameResponse(asymmetric_kernel; T=T),
        T=T, backend=selector)
    edge_acquisition = prepare_detector_acquisition(edge_detector, edge_map)
    edge_frame = capture!(edge_detector, edge_map, edge_acquisition;
        rng=MersenneTwister(303))
    AdaptiveOpticsSim.synchronize_backend!(
        AdaptiveOpticsSim.execution_style(edge_frame))
    @test sum(Array(edge_frame)) ≈ T(0.9) atol=T(2e-6) rtol=T(2e-6)
    @test sum(Array(edge_frame)) <= sum(edge_host) + T(2e-6)

    edge_cube_host = zeros(T, 2, 9, 9)
    edge_cube_host[1, 5, end] = one(T)
    edge_cube_host[2, 5, 1] = one(T)
    edge_cube = BackendArray(edge_cube_host)
    edge_scratch = similar(edge_cube)
    edge_stack = AdaptiveOpticsSim.capture_stack!(edge_detector, edge_cube,
        edge_scratch; rng=MersenneTwister(304))
    AdaptiveOpticsSim.synchronize_backend!(
        AdaptiveOpticsSim.execution_style(edge_stack))
    edge_stack_host = Array(edge_stack)
    @test sum(@view(edge_stack_host[1, :, :])) ≈ T(0.9) atol=T(2e-6) rtol=T(2e-6)
    @test sum(@view(edge_stack_host[2, :, :])) ≈ T(0.3) atol=T(2e-6) rtol=T(2e-6)

    fraunhofer = FraunhoferPropagation(field)
    propagated = propagation_output(field, fraunhofer)
    propagate_field!(propagated, field, fraunhofer)
    @test propagated.values isa BackendArray
    @test propagated.metadata.kind isa FocalPlane

    spatial_filter = SpatialFilter(tel; shape=SquareFilter(), diameter=5,
        zero_padding=2, T=T, backend=selector)
    spatial_field = ElectricField(wavefront, src; zero_padding=2, T=T,
        normalization=DimensionlessNormalization(),
        spatial_measure=PointSampledMeasure(),
        coherence=CoherentFieldCombination())
    spatial_formation = prepare_pupil_field(tel, wavefront, src, spatial_field;
        center_even_grid=false, amplitude_scale=1)
    fill_electric_field!(spatial_field, wavefront, spatial_formation)
    filtered = PupilFunction(tel; T=T, backend=selector)
    spatial_plan = prepare_spatial_filter(tel, spatial_filter, spatial_field,
        filtered)
    spatial_workspace = SpatialFilterWorkspace(spatial_filter)
    filter!(filtered, spatial_field, spatial_filter, spatial_plan,
        spatial_workspace)
    @test filtered.opd isa BackendArray
    @test filtered.amplitude isa BackendArray
    @test all(isfinite, Array(filtered.opd))
    @test all(isfinite, Array(filtered.amplitude))

    telescope_opd = copy(Array(opd_map(tel)))
    dm = DeformableMirror(tel; n_act=4, influence_width=T(0.3), T=T,
        backend=selector)
    fill!(dm.state.coefs, T(1e-8))
    update_surface!(dm, tel)
    apply_surface!(wavefront, dm, DMReplace())
    @test Array(wavefront.opd) ≈ Array(surface_opd(dm)) atol=zero(T) rtol=zero(T)
    @test Array(opd_map(tel)) == telescope_opd
    return nothing
end

function build_optional_control_loop_branch(::Type{T}, backend, label::Symbol; sensor::Symbol=:sh, seed::Integer=1) where {T<:AbstractFloat}
    tel = Telescope(resolution=16, diameter=T(8.0),
        central_obstruction=T(0.0), T=T, backend=backend)
    src = Source(band=:I, magnitude=0.0, T=T)
    atm = KolmogorovAtmosphere(tel; r0=T(0.2), L0=T(25.0), T=T, backend=backend)
    dm = DeformableMirror(tel; n_act=4, influence_width=T(0.3), T=T, backend=backend)
    wfs = sensor == :sh ?
        ShackHartmannWFS(tel; n_lenslets=4, mode=Diffractive(), T=T, backend=backend) :
        PyramidWFS(tel; pupil_samples=4, modulation=T(1.0), mode=Diffractive(), T=T, backend=backend)
    det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0), binning=1, T=T, backend=backend)
    sim = AOSimulation(tel, src, atm, dm, wfs)
    imat = interaction_matrix(dm, wfs, tel, src; amplitude=T(0.05))
    recon = ModalReconstructor(imat; gain=T(0.5))
    return ControlLoopBranch(label, sim, recon; wfs_detector=det, rng=MersenneTwister(seed))
end

function run_optional_scalable_reconstructor_checks(::Type{T}, selector,
    BackendArray) where {T<:AbstractFloat}
    rng = MersenneTwister(20260713)
    n_slopes = 12
    n_commands = 8
    D_host = randn(rng, T, n_slopes, n_commands)
    slopes_host = randn(rng, T, n_slopes)
    D = BackendArray(D_host)
    input = BackendArray(slopes_host)
    interaction = InteractionMatrix(D, T(0.1))
    dense = ModalReconstructor(interaction; gain=T(0.5))
    factorized = FactorizedReconstructor(interaction; gain=T(0.5))
    dense_out = BackendArray(zeros(T, n_commands))
    factorized_out = BackendArray(zeros(T, n_commands))
    reconstruct!(dense_out, dense, input)
    reconstruct!(factorized_out, factorized, input)
    AdaptiveOpticsSim.synchronize_backend!(
        AdaptiveOpticsSim.execution_style(factorized_out))
    @test Array(factorized_out) ≈ Array(dense_out) atol=T(2e-4) rtol=T(2e-4)
    @test AdaptiveOpticsSim.factorized_rank(factorized) == n_commands
    @test all(storage -> storage isa BackendArray,
        AdaptiveOpticsSim.runtime_reconstructor_storage(factorized))

    compact = FactorizedReconstructor(interaction; gain=T(0.5), max_rank=3)
    @test AdaptiveOpticsSim.factorized_rank(compact) == 3
    @test sum(length,
        AdaptiveOpticsSim.runtime_reconstructor_storage(compact)) <
        sum(length,
            AdaptiveOpticsSim.runtime_reconstructor_storage(factorized))

    controller = DiscreteIntegratorController(n_commands;
        gain=T(0.7), tau=T(0.01), T=T, backend=selector)
    controlled = ControlledReconstructor(factorized, controller; dt=T(1e-3))
    reconstruct!(factorized_out, controlled, input)
    AdaptiveOpticsSim.synchronize_backend!(
        AdaptiveOpticsSim.execution_style(factorized_out))
    @test all(isfinite, Array(factorized_out))
    @test all(storage -> storage isa BackendArray,
        AdaptiveOpticsSim.runtime_reconstructor_storage(controlled))
    @test AdaptiveOpticsSim.reset_controller!(controlled) === controlled
    @test_throws InvalidConfiguration ControlledReconstructor(
        factorized,
        DiscreteIntegratorController(n_commands; T=T, backend=CPUBackend());
        dt=T(1e-3),
    )
    return nothing
end

@inline _optional_low_order_label(::Val{:tiptilt}) = :tiptilt
@inline _optional_low_order_label(::Val{:steering}) = :steering
@inline _optional_low_order_label(::Val{:focus}) = :focus

function _build_optional_low_order_wfs(tel::Telescope, backend, ::Type{T}, ::Val{:sh}) where {T<:AbstractFloat}
    return ShackHartmannWFS(tel; n_lenslets=4, mode=Diffractive(), T=T, backend=backend)
end

function _build_optional_low_order_wfs(tel::Telescope, backend, ::Type{T}, ::Val{:pyr}) where {T<:AbstractFloat}
    return PyramidWFS(tel; pupil_samples=4, modulation=T(1.0), mode=Diffractive(), T=T, backend=backend)
end

function _build_optional_low_order_wfs(tel::Telescope, backend, ::Type{T}, ::Val{:bio}) where {T<:AbstractFloat}
    return BioEdgeWFS(tel; pupil_samples=4, modulation=T(1.0), mode=Diffractive(), T=T, backend=backend)
end

function _build_optional_composite_optic_case(backend, ::Type{T}, ::Val{:tiptilt}, wfs_case::Val{W}=Val(:sh)) where {T<:AbstractFloat,W}
    tel = Telescope(resolution=16, diameter=T(8.0),
        central_obstruction=T(0.0), T=T, backend=backend)
    src = Source(band=:I, magnitude=0.0, T=T)
    atm = OptionalStaticAtmosphere(tel; T=T, backend=backend)
    tiptilt = TipTiltMirror(tel; scale=T(0.1), T=T, backend=backend, label=:tiptilt)
    dm = DeformableMirror(tel; n_act=4, influence_width=T(0.3), T=T, backend=backend)
    optic = CompositeControllableOptic(:tiptilt => tiptilt, :dm => dm)
    wfs = _build_optional_low_order_wfs(tel, backend, T, wfs_case)
    det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0), binning=1, T=T, backend=backend)
    sim = AOSimulation(tel, src, atm, optic, wfs)
    return build_control_loop_scenario(
        SingleControlLoopConfig(atmosphere_step=T(1e-3),
            name=Symbol(:optional_backend_composite_tiptilt_, W), branch_label=:main,
            outputs=RuntimeOutputRequirements(slopes=true, wfs_pixels=true, science_pixels=false)),
        ControlLoopBranch(:main, sim, NullReconstructor(); wfs_detector=det, rng=MersenneTwister(91)),
    )
end

function _build_optional_composite_optic_case(backend, ::Type{T}, ::Val{:steering}, ::Val{:sh}=Val(:sh)) where {T<:AbstractFloat}
    tel = Telescope(resolution=16, diameter=T(8.0),
        central_obstruction=T(0.0), T=T, backend=backend)
    src = Source(band=:I, magnitude=0.0, T=T)
    atm = OptionalStaticAtmosphere(tel; T=T, backend=backend)
    steering = ModalControllableOptic(tel, CartesianTiltBasis(; scale=T(0.1));
        labels=:steering, T=T, backend=backend)
    dm = DeformableMirror(tel; n_act=4, influence_width=T(0.3), T=T, backend=backend)
    optic = CompositeControllableOptic(:steering => steering, :dm => dm)
    wfs = ShackHartmannWFS(tel; n_lenslets=4, mode=Diffractive(), T=T, backend=backend)
    det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0), binning=1, T=T, backend=backend)
    sim = AOSimulation(tel, src, atm, optic, wfs)
    return build_control_loop_scenario(
        SingleControlLoopConfig(atmosphere_step=T(1e-3),
            name=:optional_backend_composite_steering, branch_label=:main,
            outputs=RuntimeOutputRequirements(slopes=true, wfs_pixels=true, science_pixels=false)),
        ControlLoopBranch(:main, sim, NullReconstructor(); wfs_detector=det, rng=MersenneTwister(91)),
    )
end

function _build_optional_composite_optic_case(backend, ::Type{T}, ::Val{:focus}, ::Val{:sh}=Val(:sh)) where {T<:AbstractFloat}
    tel = Telescope(resolution=16, diameter=T(8.0),
        central_obstruction=T(0.0), T=T, backend=backend)
    src = Source(band=:I, magnitude=0.0, T=T)
    atm = OptionalStaticAtmosphere(tel; T=T, backend=backend)
    focus = FocusStage(tel; scale=T(0.1), T=T, backend=backend, label=:focus)
    dm = DeformableMirror(tel; n_act=4, influence_width=T(0.3), T=T, backend=backend)
    optic = CompositeControllableOptic(:focus => focus, :dm => dm)
    wfs = ShackHartmannWFS(tel; n_lenslets=4, mode=Diffractive(), T=T, backend=backend)
    det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0), binning=1, T=T, backend=backend)
    sim = AOSimulation(tel, src, atm, optic, wfs)
    return build_control_loop_scenario(
        SingleControlLoopConfig(atmosphere_step=T(1e-3),
            name=:optional_backend_composite_focus, branch_label=:main,
            outputs=RuntimeOutputRequirements(slopes=true, wfs_pixels=true, science_pixels=false)),
        ControlLoopBranch(:main, sim, NullReconstructor(); wfs_detector=det, rng=MersenneTwister(91)),
    )
end

function _optional_low_order_commands(::Type{T}, ::Val{:focus}) where {T<:AbstractFloat}
    return fill(T(1.25e-7), 1), fill(T(2.5e-7), 1), fill(T(2e-8), 16)
end

function _optional_low_order_commands(::Type{T}, ::Union{Val{:tiptilt},Val{:steering}}) where {T<:AbstractFloat}
    # These coefficients produce nanometre-scale OPD. Centimetre-scale test
    # commands make diffraction output chaotic under otherwise sub-ulp
    # Float32 backend differences and do not represent an AO operating point.
    return fill(T(1.25e-7), 2), fill(T(2.5e-7), 2), fill(T(2e-8), 16)
end

function _optional_low_order_tolerances(::Val{:focus}, ::Val{:sh}=Val(:sh))
    return (slopes_rtol=5f-3, slopes_atol=6f-3, frame_rtol=6f-3, frame_atol=1f6)
end

function _optional_low_order_tolerances(::Union{Val{:tiptilt},Val{:steering}}, ::Union{Val{:sh},Val{:pyr}})
    return (slopes_rtol=3f-3, slopes_atol=4f-3, frame_rtol=4f-3, frame_atol=1f6)
end

function _optional_low_order_tolerances(::Val{:tiptilt}, ::Val{:bio})
    return (slopes_rtol=1.5f-1, slopes_atol=4f-3, frame_rtol=6f-1, frame_atol=1f6)
end

function _run_optional_composite_optic_case(::Type{B}, case::Val{K}, wfs_case::Val{W}=Val(:sh)) where {B<:AdaptiveOpticsSim.GPUBackendTag,K,W}
    T = Float32
    cpu = _build_optional_composite_optic_case(CPUBackend(), T, case, wfs_case)
    gpu = _build_optional_composite_optic_case(backend_selector(B), T, case, wfs_case)
    prepare!(cpu)
    prepare!(gpu)
    initial_low_order, updated_low_order, dm_cmd = _optional_low_order_commands(T, case)
    label = _optional_low_order_label(case)
    tol = _optional_low_order_tolerances(case, wfs_case)
    set_command!(cpu, NamedTuple{(label, :dm)}((initial_low_order, dm_cmd)))
    set_command!(gpu, NamedTuple{(label, :dm)}((initial_low_order, dm_cmd)))
    sense!(cpu)
    sense!(gpu)
    AdaptiveOpticsSim.synchronize_backend!(AdaptiveOpticsSim.execution_style(command(gpu)))
    AdaptiveOpticsSim.synchronize_backend!(AdaptiveOpticsSim.execution_style(slopes(gpu)))
    AdaptiveOpticsSim.synchronize_backend!(AdaptiveOpticsSim.execution_style(wfs_frame(gpu)))
    @test AdaptiveOpticsSim.command_segment_labels(AdaptiveOpticsSim.command_layout(gpu)) == (label, :dm)
    @test isapprox(Array(command(gpu)), Array(command(cpu)); rtol=1f-6, atol=1f-6)
    @test isapprox(Array(slopes(gpu)), Array(slopes(cpu)); rtol=tol.slopes_rtol, atol=tol.slopes_atol)
    @test isapprox(Array(wfs_frame(gpu)), Array(wfs_frame(cpu)); rtol=tol.frame_rtol, atol=tol.frame_atol)

    update_command!(cpu, NamedTuple{(label,)}((updated_low_order,)))
    update_command!(gpu, NamedTuple{(label,)}((updated_low_order,)))
    sense!(cpu)
    sense!(gpu)
    AdaptiveOpticsSim.synchronize_backend!(AdaptiveOpticsSim.execution_style(command(gpu)))
    AdaptiveOpticsSim.synchronize_backend!(AdaptiveOpticsSim.execution_style(slopes(gpu)))
    AdaptiveOpticsSim.synchronize_backend!(AdaptiveOpticsSim.execution_style(wfs_frame(gpu)))
    @test isapprox(Array(command(gpu)), Array(command(cpu)); rtol=1f-6, atol=1f-6)
    @test isapprox(Array(slopes(gpu)), Array(slopes(cpu)); rtol=tol.slopes_rtol, atol=tol.slopes_atol)
    @test isapprox(Array(wfs_frame(gpu)), Array(wfs_frame(cpu)); rtol=tol.frame_rtol, atol=tol.frame_atol)
    return nothing
end

function run_optional_composite_optic_parity(::Type{B}, BackendArray) where {B<:AdaptiveOpticsSim.GPUBackendTag}
    _run_optional_composite_optic_case(B, Val(:tiptilt), Val(:sh))
    _run_optional_composite_optic_case(B, Val(:tiptilt), Val(:pyr))
    _run_optional_composite_optic_case(B, Val(:tiptilt), Val(:bio))
    _run_optional_composite_optic_case(B, Val(:steering), Val(:sh))
    _run_optional_composite_optic_case(B, Val(:focus), Val(:sh))
    return nothing
end

function run_optional_counting_detector_parity(::Type{B}, BackendArray) where {B<:AdaptiveOpticsSim.GPUBackendTag}
    T = Float32
    selector = backend_selector(B)
    correlation = ChannelCrosstalkModel(T(0.4))
    sensor_kwargs = (
        qe=T(1),
        dark_count_rate=T(0),
        fill_factor=T(1),
        energy_resolution=T(12),
        timing_jitter_s=T(2e-6),
        wavelength_range_m=(T(0.8e-6), T(1.4e-6)),
        correlation_model=correlation,
        T=T,
    )
    cpu = MKIDArrayDetector(
        noise=NoiseNone(),
        sensor=MKIDArraySensor(; sensor_kwargs...),
        output_type=UInt16,
        T=T,
        backend=CPUBackend(),
    )
    gpu = MKIDArrayDetector(
        noise=NoiseNone(),
        sensor=MKIDArraySensor(; sensor_kwargs...),
        output_type=UInt16,
        T=T,
        backend=selector,
    )
    input = zeros(T, 5, 5)
    input[3, 3] = T(10)
    cpu_output = capture!(cpu, input, MersenneTwister(19))
    gpu_output = capture!(gpu, BackendArray(input), MersenneTwister(19))
    AdaptiveOpticsSim.synchronize_backend!(AdaptiveOpticsSim.execution_style(gpu_output))
    @test gpu_output isa BackendArray
    @test isapprox(Array(gpu_output), cpu_output; rtol=1f-6, atol=1f-6)
    @test isapprox(sum(Array(gpu_output)), sum(cpu_output); rtol=1f-6, atol=1f-6)

    inside_band = Source(band=:custom, wavelength=T(1.0e-6), T=T)
    outside_band = Source(band=:custom, wavelength=T(0.55e-6), T=T)
    gpu_inside = capture!(gpu, BackendArray(fill(T(2), 2, 2)), inside_band,
        MersenneTwister(20))
    AdaptiveOpticsSim.synchronize_backend!(AdaptiveOpticsSim.execution_style(gpu_inside))
    inside_host = Array(gpu_inside)
    gpu_outside = capture!(gpu, BackendArray(fill(T(2), 2, 2)), outside_band,
        MersenneTwister(20))
    AdaptiveOpticsSim.synchronize_backend!(AdaptiveOpticsSim.execution_style(gpu_outside))
    outside_host = Array(gpu_outside)
    @test inside_host == fill(T(2), 2, 2)
    @test outside_host == zeros(T, 2, 2)
    return nothing
end

function run_optional_avalanche_detector_parity(::Type{B}, BackendArray) where {B<:AdaptiveOpticsSim.GPUBackendTag}
    T = Float32
    selector = backend_selector(B)
    input = BackendArray(fill(one(T), 128, 128))

    pc_detector = Detector(noise=NoiseNone(), qe=one(T), gain=T(10),
        sensor=EMCCDSensor(operating_mode=PhotonCountingEMMode(
            threshold=T(5), detection_efficiency=T(0.5)), T=T),
        response_model=NullFrameResponse(), T=T, backend=selector)
    pc_output = capture!(pc_detector, input; rng=MersenneTwister(2026))
    AdaptiveOpticsSim.synchronize_backend!(
        AdaptiveOpticsSim.execution_style(pc_output))
    pc_host = Array(pc_output)
    @test all(x -> x == zero(T) || x == one(T), pc_host)
    @test T(0.47) <= sum(pc_host) / length(pc_host) <= T(0.53)

    stochastic_detector = Detector(noise=NoiseNone(), qe=one(T), gain=T(5),
        sensor=EMCCDSensor(excess_noise_factor=T(1.4),
            multiplication_model=AdaptiveOpticsSim.StochasticMultiplicationRegister(T(0.6)),
            T=T), response_model=NullFrameResponse(), T=T, backend=selector)
    stochastic_input = BackendArray(fill(T(50), 128, 128))
    stochastic_output = capture!(stochastic_detector, stochastic_input;
        rng=MersenneTwister(2027))
    AdaptiveOpticsSim.synchronize_backend!(
        AdaptiveOpticsSim.execution_style(stochastic_output))
    stochastic_host = Array(stochastic_output)
    @test all(x -> x >= zero(T), stochastic_host)
    @test isapprox(sum(stochastic_host) / length(stochastic_host), T(250);
        rtol=T(0.03))

    ramp_detector = Detector(integration_time=T(2), noise=NoiseNone(),
        qe=one(T), gain=one(T),
        sensor=HgCdTeAvalancheArraySensor(read_time=zero(T),
            sampling_mode=UpTheRampSampling(5), T=T),
        response_model=NullFrameResponse(), T=T, backend=selector)
    ramp_output = capture!(ramp_detector,
        BackendArray(fill(T(3), 32, 32)); rng=MersenneTwister(2028))
    AdaptiveOpticsSim.synchronize_backend!(
        AdaptiveOpticsSim.execution_style(ramp_output))
    @test Array(ramp_output) == fill(T(6), 32, 32)
    @test Array(detector_ramp_slope(ramp_detector)) == fill(T(3), 32, 32)
    @test maximum(abs, Array(detector_ramp_intercept(ramp_detector))) <= eps(T)
    @test size(detector_ramp_cube(ramp_detector)) == (32, 32, 5)
    @test detector_ramp_times(ramp_detector) == T[0, 0.5, 1, 1.5, 2]

    linear_apd = LinearAPDDetector(topology=APDChannelBank(4),
        integration_time=T(0.5), qe=T(0.5), avalanche_gain=T(4),
        conversion_gain=T(2), noise=NoiseNone(), T=T, backend=selector)
    linear_output = capture!(linear_apd, BackendArray(fill(T(10), 4));
        rng=MersenneTwister(2029))
    AdaptiveOpticsSim.synchronize_backend!(
        AdaptiveOpticsSim.execution_style(linear_output))
    @test linear_output isa BackendArray
    @test Array(linear_output) == fill(T(20), 4)
    return nothing
end

function import_backend_package!(::Type{AdaptiveOpticsSim.CUDABackendTag})
    @eval import CUDA
    return nothing
end

function import_backend_package!(::Type{AdaptiveOpticsSim.AMDGPUBackendTag})
    @eval import AMDGPU
    return nothing
end

function backend_functional(::Type{AdaptiveOpticsSim.CUDABackendTag})
    return Base.invokelatest(getproperty(getfield(Main, :CUDA), :functional))
end

function backend_functional(::Type{AdaptiveOpticsSim.AMDGPUBackendTag})
    return Base.invokelatest(getproperty(getfield(Main, :AMDGPU), :functional))
end

run_optional_backend_plan_checks(::Type{<:AdaptiveOpticsSim.GPUBackendTag}, tel, backend) = nothing

function run_optional_lift_fallback_check(array_backend, ::Type{T}) where {T<:AbstractFloat}
    H_host = T[1 0; 0 2; 1 1]
    residual_host = T[2, -1, 0.5]
    H = array_backend(H_host)
    residual = array_backend(residual_host)
    rhs = array_backend(zeros(T, 2))
    damping = LiFTLevenbergMarquardt(lambda0=T(0.1),
        growth=T(10), condition_rtol=T(1e-3))
    normal = transpose(H_host) * H_host
    λ = AdaptiveOpticsSim.damping_lambda(damping, normal)
    expected = (normal + λ * I) \ (transpose(H_host) * residual_host)
    diag = AdaptiveOpticsSim.LiFTDiagnostics(
        T(NaN), T(NaN), T(NaN), T(NaN), zero(T), false, false)
    AdaptiveOpticsSim.solve_lift_fallback!(
        diag, rhs, H, residual, damping)
    @test Array(rhs) ≈ expected rtol=T(1e-4) atol=T(1e-5)
    @test diag.regularization == λ
    @test diag.used_fallback
    return nothing
end

function run_optional_backend_plan_checks(::Type{AdaptiveOpticsSim.AMDGPUBackendTag}, tel, backend)
    T = Float32
    array_backend = AdaptiveOpticsSim._resolve_array_backend(backend)
    sh = ShackHartmannWFS(tel; n_lenslets=4, mode=Diffractive(), T=T, backend=backend)
    sh_large = ShackHartmannWFS(tel; n_lenslets=8, mode=Diffractive(), T=T, backend=backend)
    pyr = PyramidWFS(tel; pupil_samples=4, modulation=T(1.0), mode=Diffractive(), T=T, backend=backend)
    bio = BioEdgeWFS(tel; pupil_samples=4, modulation=T(1.0), mode=Diffractive(), T=T, backend=backend)
    det = Detector(noise=NoiseReadout(T(1.0)), qe=1.0, sensor=HgCdTeAvalancheArraySensor(T=T), T=T, backend=backend)
    det_capture = Detector(noise=NoiseReadout(T(1.0)), qe=1.0, bits=12, full_well=T(100),
        sensor=CMOSSensor(T=T), T=T, backend=backend)
    src = Source(band=:I, magnitude=0.0, T=T)
    atm = MultiLayerAtmosphere(tel;
        r0=T(0.2),
        L0=T(25.0),
        fractional_cn2=T[1.0],
        wind_speed=T[0.0],
        wind_direction=T[0.0],
        altitude=T[0.0],
        T=T,
        backend=backend,
    )
    geom_prop = AtmosphericFieldPropagation(atm, tel, src;
        model=GeometricAtmosphericPropagation(T=T),
        zero_padding=1,
        T=T)
    fresnel_prop = AtmosphericFieldPropagation(atm, tel, src;
        model=LayeredFresnelAtmosphericPropagation(T=T),
        zero_padding=1,
        T=T)
    @test AdaptiveOpticsSim.grouped_accumulation_plan(AdaptiveOpticsSim.execution_style(pyr.front_end.propagation.intensity), pyr) isa AdaptiveOpticsSim.GroupedStaged2DPlan
    @test AdaptiveOpticsSim.grouped_accumulation_plan(AdaptiveOpticsSim.execution_style(bio.front_end.propagation.intensity), bio) isa AdaptiveOpticsSim.GroupedStaged2DPlan
    @test typeof(AdaptiveOpticsSim.sh_sensing_execution_plan(
        AdaptiveOpticsSim.execution_style(slopes(sh)), sh)) ===
        AdaptiveOpticsSim.ShackHartmannWFSRocmHostStatsPlan
    @test typeof(AdaptiveOpticsSim.sh_sensing_execution_plan(
        AdaptiveOpticsSim.execution_style(slopes(sh_large)), sh_large)) ===
        AdaptiveOpticsSim.ShackHartmannWFSRocmHostStatsPlan
    AdaptiveOpticsSim.prepare_sampling!(sh, tel, src)
    sh_sub = div(tel.params.resolution, AdaptiveOpticsSim.n_lenslets(sh))
    sh_pad = size(sh.optical_workspace.field, 1)
    sh_offset = div(sh_pad - sh_sub, 2)
    safe_intensity = AdaptiveOpticsSim.compute_intensity_safe!(
        AdaptiveOpticsSim.execution_style(sh.optical_workspace.intensity),
        sh, tel, src, 1, 1, sh_sub, sh_sub, sh_offset, sh_offset, sh_sub)
    @test safe_intensity === sh.optical_workspace.intensity
    @test all(isfinite, Array(safe_intensity))
    @test AdaptiveOpticsSim.detector_execution_plan(typeof(AdaptiveOpticsSim.execution_style(det.state.frame)), typeof(det)) isa AdaptiveOpticsSim.DetectorHostMirrorPlan
    capture_psf = array_backend{T}(undef, 4, 4)
    fill!(capture_psf, T(10))
    captured = capture!(det_capture, capture_psf; rng=MersenneTwister(2))
    @test captured isa array_backend
    @test maximum(Array(captured)) <= Float64(exp2(T(12)) - one(T))
    cpu_poisson_det = Detector(noise=NoisePhoton(), integration_time=T(1.0), qe=T(1.0),
        sensor=CMOSSensor(T=T), response_model=NullFrameResponse(), T=T,
        backend=CPUBackend())
    gpu_poisson_det = Detector(noise=NoisePhoton(), integration_time=T(1.0), qe=T(1.0),
        sensor=CMOSSensor(T=T), response_model=NullFrameResponse(), T=T,
        backend=backend)
    poisson_input = fill(T(10), 4, 4)
    cpu_poisson = capture!(cpu_poisson_det, copy(poisson_input); rng=MersenneTwister(23))
    gpu_poisson = capture!(gpu_poisson_det, array_backend(copy(poisson_input));
        rng=MersenneTwister(23))
    @test Array(gpu_poisson) == cpu_poisson
    poisson_method = which(
        AdaptiveOpticsSim._poisson_noise_frame!,
        (
            AdaptiveOpticsSim.DetectorHostMirrorPlan,
            typeof(gpu_poisson_det),
            typeof(MersenneTwister(1)),
            typeof(gpu_poisson_det.state.frame),
        ),
    )
    @test occursin("AdaptiveOpticsSimAMDGPUExt", String(poisson_method.file))
    @test AdaptiveOpticsSim.reduction_execution_plan(pyr.front_end.propagation.intensity) isa AdaptiveOpticsSim.HostMirrorReductionPlan
    @test AdaptiveOpticsSim.backend_sum_value(capture_psf) == T(160)

    phase_freqs = T[-0.2, -0.1, 0.1, 0.2]
    cpu_phase_psd = zeros(T, 4, 4)
    gpu_phase_freqs = array_backend(phase_freqs)
    gpu_phase_psd = array_backend(zeros(T, 4, 4))
    phase_args = (T(0.02), T(4pi^2), T(0.01), T(-11 / 6), zero(T), 4)
    AdaptiveOpticsSim._fill_phase_psd!(AdaptiveOpticsSim.ScalarCPUStyle(),
        cpu_phase_psd, phase_freqs, phase_args...)
    phase_style = AdaptiveOpticsSim.execution_style(gpu_phase_psd)
    AdaptiveOpticsSim._fill_phase_psd!(phase_style, gpu_phase_psd,
        gpu_phase_freqs, phase_args...)
    @test isapprox(Array(gpu_phase_psd), cpu_phase_psd; rtol=1f-6, atol=1f-7)
    phase_method = which(
        AdaptiveOpticsSim._fill_phase_psd!,
        (typeof(phase_style), typeof(gpu_phase_psd), typeof(gpu_phase_freqs),
            T, T, T, T, T, Int),
    )
    @test occursin("AdaptiveOpticsSimAMDGPUExt", String(phase_method.file))

    host_opd = reshape(T.(1:256), 16, 16) .* T(1e-9)
    host_valid = trues(4, 4)
    cpu_geometric_slopes = zeros(T, 32)
    gpu_geometric_slopes = array_backend(zeros(T, 32))
    gpu_opd = array_backend(host_opd)
    gpu_valid = array_backend(host_valid)
    AdaptiveOpticsSim.geometric_slopes!(cpu_geometric_slopes, host_opd, host_valid)
    AdaptiveOpticsSim.geometric_slopes!(gpu_geometric_slopes, gpu_opd, gpu_valid)
    @test isapprox(Array(gpu_geometric_slopes), cpu_geometric_slopes;
        rtol=1f-6, atol=1f-8)
    geometric_method = which(
        AdaptiveOpticsSim._geometric_slopes!,
        (typeof(AdaptiveOpticsSim.execution_style(gpu_geometric_slopes)),
            typeof(gpu_geometric_slopes), typeof(gpu_opd), typeof(gpu_valid),
            Int, Int, Int),
    )
    @test occursin("AdaptiveOpticsSimAMDGPUExt", String(geometric_method.file))
    randn_method = which(
        AdaptiveOpticsSim._randn_frame_noise!,
        (
            AdaptiveOpticsSim.DetectorHostMirrorPlan,
            typeof(det),
            typeof(MersenneTwister(1)),
            typeof(det.state.frame),
        ),
    )
    @test occursin("AdaptiveOpticsSimAMDGPUExt", String(randn_method.file))
    @test AdaptiveOpticsSim.atmospheric_field_execution_plan(
        AdaptiveOpticsSim.execution_style(first(geom_prop.state.slices).field.values),
        geom_prop.params.model,
    ) isa AdaptiveOpticsSim.GeometricFieldAsyncPlan
    @test AdaptiveOpticsSim.atmospheric_field_execution_plan(
        AdaptiveOpticsSim.execution_style(first(fresnel_prop.state.slices).field.values),
        fresnel_prop.params.model,
    ) isa AdaptiveOpticsSim.LayeredFresnelFieldAsyncPlan
    correction_models = (
        ReferencePixelCommonModeCorrection(1, 1),
        ReferenceRowCommonModeCorrection(1),
        ReferenceColumnCommonModeCorrection(1),
        ReferenceOutputCommonModeCorrection(2; edge_rows=1, edge_cols=1),
        CompositeFrameReadoutCorrection((
            ReferenceRowCommonModeCorrection(1),
            ReferenceColumnCommonModeCorrection(1),
        )),
    )
    corr_input = reshape(T.(1:96), 2, 6, 8)
    cpu_quant_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
        bits=8, full_well=T(100), sensor=CMOSSensor(T=T),
        response_model=NullFrameResponse(), T=T, backend=CPUBackend())
    gpu_quant_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
        bits=8, full_well=T(100), sensor=CMOSSensor(T=T),
        response_model=NullFrameResponse(), T=T, backend=backend)
    quant_input = reshape(T.(1:48), 6, 8) .* T(3)
    cpu_quant_frame = capture!(cpu_quant_det, copy(quant_input); rng=MersenneTwister(13))
    gpu_quant_frame = capture!(gpu_quant_det, array_backend(copy(quant_input)); rng=MersenneTwister(13))
    @test gpu_quant_frame isa array_backend
    @test isapprox(Array(gpu_quant_frame), cpu_quant_frame; rtol=1f-5, atol=1f-4)

    for correction_model in correction_models
        cpu_frame_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
            sensor=HgCdTeAvalancheArraySensor(T=T),
            response_model=NullFrameResponse(),
            correction_model=correction_model,
            T=T,
            backend=CPUBackend())
        gpu_frame_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
            sensor=HgCdTeAvalancheArraySensor(T=T),
            response_model=NullFrameResponse(),
            correction_model=correction_model,
            T=T,
            backend=backend)
        frame_input = reshape(T.(1:48), 6, 8)
        cpu_frame = capture!(cpu_frame_det, frame_input; rng=MersenneTwister(12))
        gpu_frame = capture!(gpu_frame_det, array_backend(copy(frame_input)); rng=MersenneTwister(12))
        @test gpu_frame isa array_backend
        @test isapprox(Array(gpu_frame), cpu_frame; rtol=1f-5, atol=1f-4)

        cpu_corr_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
            sensor=HgCdTeAvalancheArraySensor(T=T),
            response_model=NullFrameResponse(),
            correction_model=correction_model,
            T=T,
            backend=CPUBackend())
        gpu_corr_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
            sensor=HgCdTeAvalancheArraySensor(T=T),
            response_model=NullFrameResponse(),
            correction_model=correction_model,
            T=T,
            backend=backend)
        cpu_corr_cube = copy(corr_input)
        gpu_corr_cube = array_backend(copy(corr_input))
        cpu_corr_scratch = similar(cpu_corr_cube)
        gpu_corr_scratch = similar(gpu_corr_cube)
        AdaptiveOpticsSim.capture_stack!(cpu_corr_det, cpu_corr_cube, cpu_corr_scratch; rng=MersenneTwister(11))
        AdaptiveOpticsSim.capture_stack!(gpu_corr_det, gpu_corr_cube, gpu_corr_scratch; rng=MersenneTwister(11))
        @test gpu_corr_cube isa array_backend
        @test isapprox(Array(gpu_corr_cube), cpu_corr_cube; rtol=1f-5, atol=1f-4)
    end

    generalized_correction = CompositeFrameReadoutCorrection((
        ReferenceRowCommonModeCorrection(1),
        ReferenceColumnCommonModeCorrection(1),
    ))
    cpu_gen_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
        bits=8, full_well=T(100), readout_window=AdaptiveOpticsSim.FrameWindow(2:5, 3:7), output_type=UInt16,
        sensor=HgCdTeAvalancheArraySensor(T=T),
        response_model=NullFrameResponse(),
        correction_model=generalized_correction,
        T=T,
        backend=CPUBackend())
    gpu_gen_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
        bits=8, full_well=T(100), readout_window=AdaptiveOpticsSim.FrameWindow(2:5, 3:7), output_type=UInt16,
        sensor=HgCdTeAvalancheArraySensor(T=T),
        response_model=NullFrameResponse(),
        correction_model=generalized_correction,
        T=T,
        backend=backend)
    cpu_gen_input = copy(corr_input)
    gpu_gen_input = array_backend(copy(corr_input))
    cpu_gen_out = Array{UInt16}(undef, 2, 4, 5)
    gpu_gen_out = array_backend{UInt16}(undef, 2, 4, 5)
    AdaptiveOpticsSim.capture_stack!(cpu_gen_det, cpu_gen_out, cpu_gen_input; rng=MersenneTwister(14))
    AdaptiveOpticsSim.capture_stack!(gpu_gen_det, gpu_gen_out, gpu_gen_input; rng=MersenneTwister(14))
    @test gpu_gen_out isa array_backend
    @test Array(gpu_gen_out) == cpu_gen_out

    cpu_windowed_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0), binning=1,
        gain=T(1.0),
        sensor=HgCdTeAvalancheArraySensor(T=T, sampling_mode=CorrelatedDoubleSampling()),
        response_model=NullFrameResponse(),
        readout_window=AdaptiveOpticsSim.FrameWindow(2:3, 2:3),
        T=T,
        backend=CPUBackend())
    gpu_windowed_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0), binning=1,
        gain=T(1.0),
        sensor=HgCdTeAvalancheArraySensor(T=T, sampling_mode=CorrelatedDoubleSampling()),
        response_model=NullFrameResponse(),
        readout_window=AdaptiveOpticsSim.FrameWindow(2:3, 2:3),
        T=T,
        backend=backend)
    windowed_input = fill(T(10), 4, 4)
    cpu_windowed_frame = capture!(cpu_windowed_det, windowed_input; rng=MersenneTwister(15))
    gpu_windowed_frame = capture!(gpu_windowed_det, array_backend(copy(windowed_input)); rng=MersenneTwister(15))
    @test isapprox(Array(gpu_windowed_frame), cpu_windowed_frame; rtol=1f-5, atol=1f-4)
    @test isapprox(Array(AdaptiveOpticsSim.detector_reference_frame(gpu_windowed_det)), AdaptiveOpticsSim.detector_reference_frame(cpu_windowed_det); rtol=1f-5, atol=1f-4)
    @test isapprox(Array(AdaptiveOpticsSim.detector_signal_frame(gpu_windowed_det)), AdaptiveOpticsSim.detector_signal_frame(cpu_windowed_det); rtol=1f-5, atol=1f-4)
    @test isapprox(Array(AdaptiveOpticsSim.detector_combined_frame(gpu_windowed_det)), AdaptiveOpticsSim.detector_combined_frame(cpu_windowed_det); rtol=1f-5, atol=1f-4)
    @test isapprox(Array(AdaptiveOpticsSim.detector_reference_cube(gpu_windowed_det)), AdaptiveOpticsSim.detector_reference_cube(cpu_windowed_det); rtol=1f-5, atol=1f-4)
    @test isapprox(Array(AdaptiveOpticsSim.detector_signal_cube(gpu_windowed_det)), AdaptiveOpticsSim.detector_signal_cube(cpu_windowed_det); rtol=1f-5, atol=1f-4)
    @test isapprox(Array(AdaptiveOpticsSim.detector_read_cube(gpu_windowed_det)), AdaptiveOpticsSim.detector_read_cube(cpu_windowed_det); rtol=1f-5, atol=1f-4)
    @test AdaptiveOpticsSim.detector_read_times(gpu_windowed_det) == AdaptiveOpticsSim.detector_read_times(cpu_windowed_det)

    gsc_mask = array_backend(fill(one(T), 8, 8))
    gsc_basis = array_backend(reshape(T.(1:192), 8, 8, 3) .* T(1e-3))
    gsc_frame = array_backend(reshape(T.(1:64), 8, 8))
    gsc = GainSensingCamera(gsc_mask, gsc_basis; T=T)
    calibrate!(gsc, gsc_frame)
    optical_gains = compute_optical_gains!(gsc, gsc_frame)
    @test optical_gains isa array_backend
    @test all(isfinite, Array(optical_gains))

    lift_tel = Telescope(resolution=8, diameter=T(8),
        central_obstruction=zero(T), T=T, backend=backend)
    lift_src = Source(band=:I, magnitude=zero(T), T=T)
    lift_det = Detector(noise=NoiseNone(), psf_sampling=1, T=T,
        backend=backend)
    lift_basis = array_backend(rand(MersenneTwister(29), T, 8, 8, 3))
    lift_diversity = array_backend(zeros(T, 8, 8))
    lift = LiFT(lift_tel, lift_src, lift_basis, lift_det;
        diversity_opd=lift_diversity, iterations=2, img_resolution=8,
        solve_mode=LiFTSolveAuto())
    lift_direct = prepare_direct_imaging(lift_tel, PupilFunction(lift_tel),
        lift_src; zero_padding=1)
    lift_psf = intensity_values(form_direct_image!(lift_direct))
    lift_coeffs = reconstruct(lift, lift_psf, [1, 2])
    @test lift_coeffs isa array_backend
    @test all(isfinite, Array(lift_coeffs))
    run_optional_lift_fallback_check(array_backend, T)
    return nothing
end

function run_optional_backend_plan_checks(::Type{AdaptiveOpticsSim.CUDABackendTag}, tel, backend)
    T = Float32
    array_backend = AdaptiveOpticsSim._resolve_array_backend(backend)
    sh = ShackHartmannWFS(tel; n_lenslets=4, mode=Diffractive(), T=T, backend=backend)
    pyr = PyramidWFS(tel; pupil_samples=4, modulation=T(1.0), mode=Diffractive(), T=T, backend=backend)
    bio = BioEdgeWFS(tel; pupil_samples=4, modulation=T(1.0), mode=Diffractive(), T=T, backend=backend)
    det = Detector(noise=NoiseReadout(T(1.0)), qe=1.0, sensor=HgCdTeAvalancheArraySensor(T=T), T=T, backend=backend)
    src = Source(band=:I, magnitude=0.0, T=T)
    atm = MultiLayerAtmosphere(tel;
        r0=T(0.2),
        L0=T(25.0),
        fractional_cn2=T[1.0],
        wind_speed=T[0.0],
        wind_direction=T[0.0],
        altitude=T[0.0],
        T=T,
        backend=backend,
    )
    geom_prop = AtmosphericFieldPropagation(atm, tel, src;
        model=GeometricAtmosphericPropagation(T=T),
        zero_padding=1,
        T=T)
    fresnel_prop = AtmosphericFieldPropagation(atm, tel, src;
        model=LayeredFresnelAtmosphericPropagation(T=T),
        zero_padding=1,
        T=T)
    @test AdaptiveOpticsSim.grouped_accumulation_plan(AdaptiveOpticsSim.execution_style(pyr.front_end.propagation.intensity), pyr) isa AdaptiveOpticsSim.GroupedStackReducePlan
    @test AdaptiveOpticsSim.grouped_accumulation_plan(AdaptiveOpticsSim.execution_style(bio.front_end.propagation.intensity), bio) isa AdaptiveOpticsSim.GroupedStackReducePlan
    @test AdaptiveOpticsSim.sh_sensing_execution_plan(AdaptiveOpticsSim.execution_style(slopes(sh)), sh) isa AdaptiveOpticsSim.ShackHartmannWFSBatchedPlan
    @test AdaptiveOpticsSim.detector_execution_plan(typeof(AdaptiveOpticsSim.execution_style(det.state.frame)), typeof(det)) isa AdaptiveOpticsSim.DetectorDirectPlan
    @test AdaptiveOpticsSim.reduction_execution_plan(pyr.front_end.propagation.intensity) isa AdaptiveOpticsSim.DirectReductionPlan
    @test AdaptiveOpticsSim.atmospheric_field_execution_plan(
        AdaptiveOpticsSim.execution_style(first(geom_prop.state.slices).field.values),
        geom_prop.params.model,
    ) isa AdaptiveOpticsSim.GeometricFieldAsyncPlan
    @test AdaptiveOpticsSim.atmospheric_field_execution_plan(
        AdaptiveOpticsSim.execution_style(first(fresnel_prop.state.slices).field.values),
        fresnel_prop.params.model,
    ) isa AdaptiveOpticsSim.LayeredFresnelFieldAsyncPlan
    cpu_tel = Telescope(resolution=16, diameter=8.0f0,
        central_obstruction=0.0f0, T=T, backend=CPUBackend())
    gpu_tel = Telescope(resolution=16, diameter=8.0f0,
        central_obstruction=0.0f0, T=T, backend=backend)
    cpu_src = Source(band=:I, magnitude=0.0, T=T)
    gpu_src = Source(band=:I, magnitude=0.0, T=T)
    cpu_sh = ShackHartmannWFS(cpu_tel; n_lenslets=4, mode=Diffractive(), T=T, backend=CPUBackend())
    gpu_sh = ShackHartmannWFS(gpu_tel; n_lenslets=4, mode=Diffractive(), T=T, backend=backend)
    cpu_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
        sensor=CMOSSensor(T=T), response_model=NullFrameResponse(), T=T, backend=CPUBackend())
    gpu_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
        sensor=CMOSSensor(T=T), response_model=NullFrameResponse(), T=T, backend=backend)
    measure!(cpu_sh, cpu_tel, cpu_src, cpu_det; rng=MersenneTwister(3))
    measure!(gpu_sh, gpu_tel, gpu_src, gpu_det; rng=MersenneTwister(3))
    cpu_export = Array(AdaptiveOpticsSim.sh_exported_spot_cube(cpu_sh))
    gpu_export = Array(AdaptiveOpticsSim.sh_exported_spot_cube(gpu_sh))
    cpu_frame = Array(AdaptiveOpticsSim.wfs_output_frame(cpu_sh, cpu_det))
    gpu_frame = Array(AdaptiveOpticsSim.wfs_output_frame(gpu_sh, gpu_det))
    @test size(gpu_export) == size(cpu_export)
    @test isapprox(gpu_export, cpu_export; rtol=1f-5, atol=1f-4)
    @test size(gpu_frame) == size(cpu_frame)
    @test isapprox(gpu_frame, cpu_frame; rtol=1f-5, atol=1f-4)
    cpu_sh_stats = ShackHartmannWFS(cpu_tel; n_lenslets=4, mode=Diffractive(), T=T, backend=CPUBackend(),
        valid_subaperture_policy=FluxThresholdValidSubapertures(light_ratio=0.5f0))
    gpu_sh_stats = ShackHartmannWFS(gpu_tel; n_lenslets=4, mode=Diffractive(), T=T, backend=backend,
        valid_subaperture_policy=FluxThresholdValidSubapertures(light_ratio=0.5f0))
    measure!(cpu_sh_stats, cpu_tel, cpu_src, cpu_det; rng=MersenneTwister(3))
    measure!(gpu_sh_stats, gpu_tel, gpu_src, gpu_det; rng=MersenneTwister(3))
    cpu_peak = AdaptiveOpticsSim.sh_safe_peak_value(cpu_sh_stats.acquisition.spot_cube)
    cpu_cutoff = AdaptiveOpticsSim.centroid_threshold(cpu_sh_stats) * cpu_peak
    AdaptiveOpticsSim.sh_signal_from_spots!(cpu_sh_stats, cpu_cutoff)
    gpu_peak = AdaptiveOpticsSim.sh_safe_peak_value(gpu_sh_stats.acquisition.spot_cube)
    gpu_cutoff = AdaptiveOpticsSim.centroid_threshold(gpu_sh_stats) * gpu_peak
    AdaptiveOpticsSim.sh_signal_from_spots_device_stats!(
        AdaptiveOpticsSim.execution_style(slopes(gpu_sh_stats)),
        gpu_sh_stats,
        gpu_cutoff,
    )
    @test isapprox(Array(slopes(gpu_sh_stats)), slopes(cpu_sh_stats); rtol=1f-5, atol=1f-4)

    correction_models = (
        ReferencePixelCommonModeCorrection(1, 1),
        ReferenceRowCommonModeCorrection(1),
        ReferenceColumnCommonModeCorrection(1),
        ReferenceOutputCommonModeCorrection(2; edge_rows=1, edge_cols=1),
        CompositeFrameReadoutCorrection((
            ReferenceRowCommonModeCorrection(1),
            ReferenceColumnCommonModeCorrection(1),
        )),
    )
    corr_input = reshape(T.(1:96), 2, 6, 8)
    cpu_quant_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
        bits=8, full_well=T(100), sensor=CMOSSensor(T=T),
        response_model=NullFrameResponse(), T=T, backend=CPUBackend())
    gpu_quant_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
        bits=8, full_well=T(100), sensor=CMOSSensor(T=T),
        response_model=NullFrameResponse(), T=T, backend=backend)
    quant_input = reshape(T.(1:48), 6, 8) .* T(3)
    cpu_quant_frame = capture!(cpu_quant_det, copy(quant_input); rng=MersenneTwister(13))
    gpu_quant_frame = capture!(gpu_quant_det, array_backend(copy(quant_input)); rng=MersenneTwister(13))
    @test gpu_quant_frame isa array_backend
    @test isapprox(Array(gpu_quant_frame), cpu_quant_frame; rtol=1f-5, atol=1f-4)

    for correction_model in correction_models
        cpu_frame_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
            sensor=HgCdTeAvalancheArraySensor(T=T),
            response_model=NullFrameResponse(),
            correction_model=correction_model,
            T=T,
            backend=CPUBackend())
        gpu_frame_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
            sensor=HgCdTeAvalancheArraySensor(T=T),
            response_model=NullFrameResponse(),
            correction_model=correction_model,
            T=T,
            backend=backend)
        frame_input = reshape(T.(1:48), 6, 8)
        cpu_frame = capture!(cpu_frame_det, frame_input; rng=MersenneTwister(12))
        gpu_frame = capture!(gpu_frame_det, array_backend(copy(frame_input)); rng=MersenneTwister(12))
        @test gpu_frame isa array_backend
        @test isapprox(Array(gpu_frame), cpu_frame; rtol=1f-5, atol=1f-4)

        cpu_corr_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
            sensor=HgCdTeAvalancheArraySensor(T=T),
            response_model=NullFrameResponse(),
            correction_model=correction_model,
            T=T,
            backend=CPUBackend())
        gpu_corr_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
            sensor=HgCdTeAvalancheArraySensor(T=T),
            response_model=NullFrameResponse(),
            correction_model=correction_model,
            T=T,
            backend=backend)
        cpu_corr_cube = copy(corr_input)
        gpu_corr_cube = array_backend(copy(corr_input))
        cpu_corr_scratch = similar(cpu_corr_cube)
        gpu_corr_scratch = similar(gpu_corr_cube)
        AdaptiveOpticsSim.capture_stack!(cpu_corr_det, cpu_corr_cube, cpu_corr_scratch; rng=MersenneTwister(11))
        AdaptiveOpticsSim.capture_stack!(gpu_corr_det, gpu_corr_cube, gpu_corr_scratch; rng=MersenneTwister(11))
        @test gpu_corr_cube isa array_backend
        @test isapprox(Array(gpu_corr_cube), cpu_corr_cube; rtol=1f-5, atol=1f-4)
    end

    generalized_correction = CompositeFrameReadoutCorrection((
        ReferenceRowCommonModeCorrection(1),
        ReferenceColumnCommonModeCorrection(1),
    ))
    cpu_gen_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
        bits=8, full_well=T(100), readout_window=AdaptiveOpticsSim.FrameWindow(2:5, 3:7), output_type=UInt16,
        sensor=HgCdTeAvalancheArraySensor(T=T),
        response_model=NullFrameResponse(),
        correction_model=generalized_correction,
        T=T,
        backend=CPUBackend())
    gpu_gen_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
        bits=8, full_well=T(100), readout_window=AdaptiveOpticsSim.FrameWindow(2:5, 3:7), output_type=UInt16,
        sensor=HgCdTeAvalancheArraySensor(T=T),
        response_model=NullFrameResponse(),
        correction_model=generalized_correction,
        T=T,
        backend=backend)
    cpu_gen_input = copy(corr_input)
    gpu_gen_input = array_backend(copy(corr_input))
    cpu_gen_out = Array{UInt16}(undef, 2, 4, 5)
    gpu_gen_out = array_backend{UInt16}(undef, 2, 4, 5)
    AdaptiveOpticsSim.capture_stack!(cpu_gen_det, cpu_gen_out, cpu_gen_input; rng=MersenneTwister(14))
    AdaptiveOpticsSim.capture_stack!(gpu_gen_det, gpu_gen_out, gpu_gen_input; rng=MersenneTwister(14))
    @test gpu_gen_out isa array_backend
    @test Array(gpu_gen_out) == cpu_gen_out

    cpu_windowed_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0), binning=1,
        gain=T(1.0),
        sensor=HgCdTeAvalancheArraySensor(T=T, sampling_mode=CorrelatedDoubleSampling()),
        response_model=NullFrameResponse(),
        readout_window=AdaptiveOpticsSim.FrameWindow(2:3, 2:3),
        T=T,
        backend=CPUBackend())
    gpu_windowed_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0), binning=1,
        gain=T(1.0),
        sensor=HgCdTeAvalancheArraySensor(T=T, sampling_mode=CorrelatedDoubleSampling()),
        response_model=NullFrameResponse(),
        readout_window=AdaptiveOpticsSim.FrameWindow(2:3, 2:3),
        T=T,
        backend=backend)
    windowed_input = fill(T(10), 4, 4)
    cpu_windowed_frame = capture!(cpu_windowed_det, windowed_input; rng=MersenneTwister(15))
    gpu_windowed_frame = capture!(gpu_windowed_det, array_backend(copy(windowed_input)); rng=MersenneTwister(15))
    @test isapprox(Array(gpu_windowed_frame), cpu_windowed_frame; rtol=1f-5, atol=1f-4)
    @test isapprox(Array(AdaptiveOpticsSim.detector_reference_frame(gpu_windowed_det)), AdaptiveOpticsSim.detector_reference_frame(cpu_windowed_det); rtol=1f-5, atol=1f-4)
    @test isapprox(Array(AdaptiveOpticsSim.detector_signal_frame(gpu_windowed_det)), AdaptiveOpticsSim.detector_signal_frame(cpu_windowed_det); rtol=1f-5, atol=1f-4)
    @test isapprox(Array(AdaptiveOpticsSim.detector_combined_frame(gpu_windowed_det)), AdaptiveOpticsSim.detector_combined_frame(cpu_windowed_det); rtol=1f-5, atol=1f-4)
    @test isapprox(Array(AdaptiveOpticsSim.detector_reference_cube(gpu_windowed_det)), AdaptiveOpticsSim.detector_reference_cube(cpu_windowed_det); rtol=1f-5, atol=1f-4)
    @test isapprox(Array(AdaptiveOpticsSim.detector_signal_cube(gpu_windowed_det)), AdaptiveOpticsSim.detector_signal_cube(cpu_windowed_det); rtol=1f-5, atol=1f-4)
    @test isapprox(Array(AdaptiveOpticsSim.detector_read_cube(gpu_windowed_det)), AdaptiveOpticsSim.detector_read_cube(cpu_windowed_det); rtol=1f-5, atol=1f-4)
    @test AdaptiveOpticsSim.detector_read_times(gpu_windowed_det) == AdaptiveOpticsSim.detector_read_times(cpu_windowed_det)

    run_optional_lift_fallback_check(array_backend, T)
    return nothing
end

function run_optional_backend_smoke(::Type{B}) where {B<:AdaptiveOpticsSim.GPUBackendTag}
    pkg = backend_package_name(B)
    pkg_path = Base.find_package(pkg)
    if pkg_path === nothing
        @info "Skipping $(backend_label(B)) smoke: $(pkg).jl is not available in this environment"
        @test true
        return nothing
    end

    import_backend_package!(B)
    if !backend_functional(B)
        @info "Skipping $(backend_label(B)) smoke: backend runtime/device is not functional on this host"
        @test true
        return nothing
    end

    AdaptiveOpticsSim.disable_scalar_backend!(B)
    backend = AdaptiveOpticsSim.gpu_backend_array_type(B)
    @test backend !== nothing
    run_optional_backend_selector_smoke(B, backend)
    run_optional_lgs_convolution_normalization(B)
    run_optional_sodium_profile_wfs(B, backend)
    run_optional_zernike_normalization(B, backend)
    run_optional_wfs_stage_contracts(B, backend)

    if get(ENV, backend_full_smoke_env(B), "0") == "1"
        include(joinpath(dirname(@__DIR__), "scripts", "gpu_smoke_contract.jl"))
        run_smoke_matrix = Base.invokelatest(
            getfield, Main, :run_gpu_smoke_matrix)
        Base.invokelatest(run_smoke_matrix, B)
        @test true
        return nothing
    end

    T = Float32
    rng = MersenneTwister(7)
    BackendArray = backend
    selector = backend_selector(B)

    run_optional_counting_detector_parity(B, BackendArray)
    run_optional_avalanche_detector_parity(B, BackendArray)

    tel = Telescope(resolution=16, diameter=8.0f0,
        central_obstruction=0.0f0, T=T, backend=selector)
    src = Source(band=:I, magnitude=0.0, T=T)
    run_optional_plane_product_checks(tel, src, selector, BackendArray, T)

    atm = MultiLayerAtmosphere(tel;
        r0=T(0.2),
        L0=T(25.0),
        fractional_cn2=T[0.7, 0.3],
        wind_speed=T[8.0, 4.0],
        wind_direction=T[0.0, 90.0],
        altitude=T[0.0, 5000.0],
        T=T,
        backend=selector,
    )
    epoch = advance_by!(atm, T(1e-3); rng=rng)
    renderer = prepare_atmosphere_renderer(atm, tel, src)
    atmosphere_output = PupilFunction(tel; T=T)
    render_atmosphere!(atmosphere_output, renderer, atm, epoch)
    @test atm.layers[1].generator.state.opd isa BackendArray
    @test atmosphere_output.opd isa BackendArray
    copyto!(tel.state.opd, atmosphere_output.opd)

    inf_atm = InfiniteMultiLayerAtmosphere(tel;
        r0=T(0.2),
        L0=T(25.0),
        fractional_cn2=T[0.7, 0.3],
        wind_speed=T[8.0, 4.0],
        wind_direction=T[0.0, 90.0],
        altitude=T[0.0, 5000.0],
        screen_resolution=33,
        stencil_size=35,
        T=T,
        backend=selector,
    )
    infinite_epoch = advance_by!(inf_atm, T(1e-3); rng=rng)
    infinite_renderer = prepare_atmosphere_renderer(inf_atm, tel, src)
    render_atmosphere!(atmosphere_output, infinite_renderer, inf_atm,
        infinite_epoch)
    @test inf_atm.layers[1].screen.state.screen isa BackendArray
    @test atmosphere_output.opd isa BackendArray

    prop = AtmosphericFieldPropagation(atm, tel, src;
        model=GeometricAtmosphericPropagation(T=T),
        zero_padding=2,
        T=T)
    field = propagate_atmosphere_field!(prop, atm, epoch)
    @test field.values isa BackendArray
    intensity = atmospheric_intensity!(prop, atm, epoch)
    @test intensity isa BackendArray

    bundle = SpectralBundle(fill(wavelength(src), 2), T[0.4, 0.6]; T=T)
    poly = with_spectrum(src, bundle)
    sh = ShackHartmannWFS(tel; n_lenslets=4, mode=Diffractive(), T=T, backend=selector)
    slopes = measure!(sh, tel, poly)
    @test slopes isa BackendArray

    spectral_optical_sh = ShackHartmannWFS(tel; n_lenslets=4,
        mode=Diffractive(), T=T, backend=selector)
    AdaptiveOpticsSim.sampled_spots_peak!(spectral_optical_sh, tel, poly)
    spectral_optical_spots = Array(spectral_optical_sh.acquisition.spot_cube)
    spectral_qe = AdaptiveOpticsSim.SampledQuantumEfficiency(
        T[0.9 * wavelength(src), 1.1 * wavelength(src)], T[0.2, 0.8])
    spectral_exposure = T(2.5)
    spectral_detector = Detector(noise=NoiseNone(), qe=spectral_qe,
        integration_time=spectral_exposure, binning=1,
        response_model=NullFrameResponse(), T=T, backend=selector)
    spectral_detector_sh = ShackHartmannWFS(tel; n_lenslets=4,
        mode=Diffractive(), T=T, backend=selector)
    AdaptiveOpticsSim.sampled_spots_peak!(spectral_detector_sh, tel, poly,
        spectral_detector, MersenneTwister(149))
    expected_spectral_scale = spectral_exposure *
        T(AdaptiveOpticsSim.qe_at(spectral_qe, wavelength(src)))
    @test Array(spectral_detector_sh.acquisition.spot_cube) ≈
        spectral_optical_spots .* expected_spectral_scale rtol=5e-5

    distinct = with_spectrum(src, SpectralBundle(
        T[0.9 * wavelength(src), 1.1 * wavelength(src)], T[0.4, 0.6];
        T=T))
    @test_throws InvalidConfiguration measure!(sh, tel, distinct)

    science_src = Source(band=:K, magnitude=1.0, coordinates=(4.0, 90.0), T=T)
    split_dm = DeformableMirror(tel; n_act=4, influence_width=T(0.3), T=T,
        backend=selector)
    split_det = Detector(noise=NoiseNone(), integration_time=one(T), qe=one(T),
        binning=1, T=T, backend=selector)
    split_sim = AOSimulation(tel, src, atm, split_dm, sh;
        science_source=science_src)
    split_runtime = AdaptiveOpticsSim.ClosedLoopRuntime(
        split_sim,
        NullReconstructor();
        atmosphere_step=T(1e-3),
        science_detector=split_det,
        outputs=RuntimeOutputRequirements(
            slopes=true,
            wfs_pixels=false,
            science_pixels=true,
        ),
        science_zero_padding=1,
        rng=MersenneTwister(17),
    )
    @test runtime_execution_plan(split_runtime) isa DeviceResidentExecutionPlan
    @test_throws InvalidConfiguration AdaptiveOpticsSim.ClosedLoopRuntime(
        split_sim,
        NullReconstructor();
        atmosphere_step=T(1e-3),
        science_detector=split_det,
        execution_plan=CPUHILExecutionPlan(),
    )
    sense!(split_runtime)
    @test synchronize_runtime!(split_runtime) === split_runtime
    @test wfs_source(split_runtime) === src
    @test science_source(split_runtime) === science_src
    @test split_runtime.science_path isa AdaptiveOpticsSim.RepropagateScienceOpticalPath
    @test science_frame(split_runtime) isa BackendArray
    @test all(isfinite, Array(science_frame(split_runtime)))

    shared_detector_a = Detector(noise=NoiseNone(), integration_time=one(T),
        qe=one(T), binning=1, T=T, backend=selector)
    shared_detector_b = Detector(noise=NoiseNone(), integration_time=one(T),
        qe=one(T), binning=1, T=T, backend=selector)
    shared_arm = SharedOpticalArm(
        :shared_science,
        science_src;
        science_detectors=(shared_detector_a, shared_detector_b),
        science_zero_padding=1,
    )
    shared_runtime = SharedOpticalRuntime(split_runtime, shared_arm)
    sense!(shared_runtime)
    @test runtime_execution_plan(shared_runtime) isa DeviceResidentExecutionPlan
    @test all(frame -> frame isa BackendArray, science_frames(shared_arm))
    @test Array(science_frames(shared_arm)[1]) ≈
        Array(science_frames(shared_arm)[2]) atol=0 rtol=0
    @test synchronize_runtime!(shared_runtime) === shared_runtime

    run_optional_backend_plan_checks(B, tel, selector)
    run_optional_composite_optic_parity(B, BackendArray)
    run_optional_scalable_reconstructor_checks(T, selector, BackendArray)

    curv = CurvatureWFS(tel; pupil_samples=4, T=T, backend=selector)
    curv_slopes = measure!(curv, tel, src, atm)
    @test curv_slopes isa BackendArray

    single_control_loop = build_control_loop_scenario(
        SingleControlLoopConfig(
            atmosphere_step=T(1e-3),
            name=:optional_backend_single,
            branch_label=:main,
            outputs=RuntimeOutputRequirements(slopes=true, wfs_pixels=true, science_pixels=false),
        ),
        build_optional_control_loop_branch(T, selector, :main; sensor=:sh, seed=21),
    )
    prepare!(single_control_loop)
    step!(single_control_loop)
    @test runtime_execution_plan(single_control_loop) isa DeviceResidentExecutionPlan
    @test synchronize_runtime!(single_control_loop) === single_control_loop
    single_runtime = AdaptiveOpticsSim.simulation_interface(single_control_loop).runtime
    @test all(storage -> typeof(AdaptiveOpticsSim.backend(storage)) === typeof(selector),
        AdaptiveOpticsSim.runtime_reconstructor_storage(single_runtime.reconstructor))
    @test AdaptiveOpticsSim.simulation_interface(single_control_loop) isa AdaptiveOpticsSim.SimulationInterface
    @test wfs_frame(single_control_loop) isa BackendArray
    @test control_loop_branch_labels(single_control_loop) == (:main,)

    grouped_control_loop = build_control_loop_scenario(
        GroupedControlLoopConfig(
            (:hi, :lo);
            atmosphere_step=T(1e-3),
            name=:optional_backend_grouped,
            outputs=GroupedRuntimeOutputRequirements(wfs_frames=true, science_frames=false, wfs_stack=true, science_stack=false),
        ),
        build_optional_control_loop_branch(T, selector, :hi; sensor=:sh, seed=31),
        build_optional_control_loop_branch(T, selector, :lo; sensor=:sh, seed=32),
    )
    prepare!(grouped_control_loop)
    step!(grouped_control_loop)
    @test AdaptiveOpticsSim.simulation_interface(grouped_control_loop) isa AdaptiveOpticsSim.CompositeSimulationInterface
    @test grouped_wfs_stack(grouped_control_loop) isa BackendArray
    @test size(grouped_wfs_stack(grouped_control_loop), ndims(grouped_wfs_stack(grouped_control_loop))) == 2
    @test control_loop_branch_labels(grouped_control_loop) == (:hi, :lo)
    return nothing
end
