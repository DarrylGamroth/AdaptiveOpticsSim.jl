@testset "PE-01 calibration command-plan interface" begin
    actuator_commands = Matrix{Float64}(I, 3, 3)
    assert_calibration_command_plan_conformance(
        Calibration.ActuatorCalibrationCommands(3),
        actuator_commands,
        0.25,
    )

    matrix_commands = [1.0 3.0; 2.0 4.0; 5.0 6.0]
    assert_calibration_command_plan_conformance(
        Calibration.MatrixCalibrationCommands(matrix_commands),
        matrix_commands,
        0.5,
    )
end

@testset "PE-01 four-pupil LGS plan interface" begin
    target = HostComputeDevice()
    intensity = ones(4, 4)
    scratch = similar(intensity)
    fft_buffer = zeros(ComplexF64, size(intensity))
    ifft_buffer = similar(fft_buffer)
    fft_plan = Backends.plan_fft_backend!(fft_buffer)
    ifft_plan = Backends.plan_ifft_backend!(ifft_buffer)

    no_lgs = WavefrontSensors.NoPreparedFourPupilLGS()
    baseline = copy(intensity)
    assert_prepared_four_pupil_lgs_conformance!(no_lgs, intensity,
        scratch, fft_buffer, fft_plan, ifft_buffer, ifft_plan, target)
    @test intensity == baseline

    elongation = WavefrontSensors.PreparedFourPupilElongation(
        [0.25, 0.5, 0.25], 1)
    fill!(intensity, 1.0)
    assert_prepared_four_pupil_lgs_conformance!(elongation, intensity,
        scratch, fft_buffer, fft_plan, ifft_buffer, ifft_plan, target)
    @test intensity == ones(4, 4)

    sodium = WavefrontSensors.PreparedFourPupilSodiumLayerProfile(
        ones(ComplexF64, size(intensity)))
    intensity .= reshape(collect(1.0:16.0), 4, 4)
    baseline = copy(intensity)
    assert_prepared_four_pupil_lgs_conformance!(sodium, intensity,
        scratch, fft_buffer, fft_plan, ifft_buffer, ifft_plan, target)
    @test intensity ≈ baseline atol=1e-12 rtol=1e-12

    unavailable = AcceleratorComputeDevice(CUDABackend(), 0)
    @test WavefrontSensors._require_exact_prepared_four_pupil_lgs_target(
        no_lgs, unavailable) === nothing
    @test_throws WFSPreparationError begin
        WavefrontSensors._require_exact_prepared_four_pupil_lgs_target(
            elongation, unavailable)
    end
    @test_throws WFSPreparationError begin
        WavefrontSensors._require_exact_prepared_four_pupil_lgs_target(
            sodium, unavailable)
    end
end

@testset "PE-01 prepared device-context interface" begin
    target = HostComputeDevice()
    context = Backends._prepare_device_execution_context(target)
    observe_host() = HostComputeDevice()
    @test @inferred(Backends._with_prepared_device_execution_context(
        observe_host, context)) == target
    assert_prepared_device_execution_context_conformance(
        context, target, observe_host, target)

    @test_throws ComputeDeviceError begin
        Backends._prepare_device_execution_context(
            AcceleratorComputeDevice(CUDABackend(), 0))
    end
end
