function detector_qualification_response_kernel(::NullFrameResponse)
    return ones(1, 1)
end

function detector_qualification_response_kernel(
    response::GaussianPixelResponse)
    kernel = Array(response.kernel)
    return kernel * transpose(kernel)
end

function detector_qualification_response_kernel(
    response::SampledFrameResponse)
    return Array(response.kernel)
end

function detector_qualification_response_kernel(
    response::RectangularPixelAperture)
    return Array(response.kernel_y) * transpose(Array(response.kernel_x))
end

function detector_qualification_zero_extended_convolution(input, kernel)
    output = zeros(promote_type(eltype(input), eltype(kernel)), size(input))
    center_i = fld(size(kernel, 1), 2) + 1
    center_j = fld(size(kernel, 2), 2) + 1
    @inbounds for j in axes(output, 2), i in axes(output, 1)
        value = zero(eltype(output))
        for kj in axes(kernel, 2), ki in axes(kernel, 1)
            input_i = i - (ki - center_i)
            input_j = j - (kj - center_j)
            if checkbounds(Bool, input, input_i, input_j)
                value = muladd(kernel[ki, kj],
                    input[input_i, input_j], value)
            end
        end
        output[i, j] = value
    end
    return output
end

function detector_qualification_block_sum(input, factor)
    output = zeros(eltype(input), div(size(input, 1), factor),
        div(size(input, 2), factor))
    @inbounds for j in axes(output, 2), i in axes(output, 1)
        value = zero(eltype(output))
        for input_j in ((j - 1) * factor + 1):(j * factor)
            for input_i in ((i - 1) * factor + 1):(i * factor)
                value += input[input_i, input_j]
            end
        end
        output[i, j] = value
    end
    return output
end

function detector_qualification_transfer_magnitude(kernel, fx, fy)
    center_i = fld(size(kernel, 1), 2) + 1
    center_j = fld(size(kernel, 2), 2) + 1
    response = zero(ComplexF64)
    normalization = zero(Float64)
    @inbounds for j in axes(kernel, 2), i in axes(kernel, 1)
        weight = Float64(kernel[i, j])
        phase = -2pi * (fy * (i - center_i) + fx * (j - center_j))
        response += weight * cis(phase)
        normalization += weight
    end
    return abs(response) / normalization
end

function detector_qualification_samples(detector, input, rng,
    frame_count::Int)
    samples = Vector{Float64}(undef, frame_count * length(input))
    offset = 0
    for _ in 1:frame_count
        frame = capture!(detector, input, rng)
        copyto!(samples, offset + 1, vec(frame), 1, length(frame))
        offset += length(frame)
    end
    return samples
end

function test_detector_poisson_moments(samples, expected_mean;
    sigma_limit=6.0)
    sample_count = length(samples)
    mean_limit = sigma_limit * sqrt(expected_mean / sample_count)
    variance_limit = sigma_limit *
        sqrt((expected_mean + 2 * expected_mean^2) / (sample_count - 1))
    @test abs(mean(samples) - expected_mean) <= mean_limit
    @test abs(var(samples) - expected_mean) <= variance_limit
end

function test_detector_gaussian_moments(samples, expected_sigma;
    sigma_limit=6.0)
    sample_count = length(samples)
    expected_variance = expected_sigma^2
    mean_limit = sigma_limit * expected_sigma / sqrt(sample_count)
    variance_limit = sigma_limit * expected_variance *
        sqrt(2 / (sample_count - 1))
    @test abs(mean(samples)) <= mean_limit
    @test abs(var(samples) - expected_variance) <= variance_limit
end

@testset "Shared frame-detector deterministic qualification" begin
    sampled_kernel = [
        0.00 0.05 0.00
        0.10 0.55 0.20
        0.00 0.10 0.00
    ]
    responses = (
        NullFrameResponse(),
        GaussianPixelResponse(response_width_px=0.65),
        RectangularPixelAperture(pitch_x_px=2.0, pitch_y_px=3.0,
            fill_factor_x=0.7, fill_factor_y=0.8),
        SampledFrameResponse(sampled_kernel),
    )
    impulse = zeros(11, 11)
    impulse[6, 6] = 1.0
    frequencies = ((0.0, 0.0), (0.19, 0.0), (0.0, 0.31),
        (0.19, 0.31))

    for response in responses
        detector = Detector(noise=NoiseNone(), integration_time=1.0, qe=1.0,
            response_model=response)
        observed = copy(capture!(detector, impulse, Xoshiro(2001)))
        kernel = detector_qualification_response_kernel(
            detector.params.response_model)
        expected = detector_qualification_zero_extended_convolution(
            impulse, kernel)
        @test observed ≈ expected atol=2e-15 rtol=2e-15
        for (fx, fy) in frequencies
            expected_mtf = detector_qualification_transfer_magnitude(
                kernel, fx, fy)
            @test detector_mtf(detector, fx, fy) ≈ expected_mtf
        end
    end

    asymmetric = SampledFrameResponse(sampled_kernel)
    right_edge = zeros(11, 11)
    right_edge[6, end] = 1.0
    expected_edge = detector_qualification_zero_extended_convolution(
        right_edge, Array(asymmetric.kernel))
    edge_detector = Detector(noise=NoiseNone(), integration_time=1.0, qe=1.0,
        response_model=asymmetric)
    observed_edge = copy(capture!(edge_detector, right_edge, Xoshiro(2002)))
    @test observed_edge == expected_edge
    @test sum(observed_edge) < sum(right_edge)

    ipc_kernel = [
        0.00 0.02 0.00
        0.03 0.90 0.01
        0.00 0.04 0.00
    ]
    ipc_detector = Detector(noise=NoiseNone(), integration_time=2.0, qe=0.5,
        response_model=NullFrameResponse(),
        charge_coupling_model=InterpixelCapacitance(ipc_kernel))
    collected_charge = impulse
    expected_ipc = detector_qualification_zero_extended_convolution(
        collected_charge, ipc_kernel)
    observed_ipc = copy(capture!(ipc_detector, impulse, Xoshiro(2003)))
    @test observed_ipc ≈ expected_ipc
    @test !supports_detector_mtf(ipc_detector)
    @test detector_mtf(ipc_detector, 0.19, 0.31) == 1.0

    input = reshape(collect(1.0:36.0), 6, 6)
    response = SampledFrameResponse([
        0.00 0.10 0.00
        0.05 0.70 0.10
        0.00 0.05 0.00
    ])
    prnu = reshape(collect(range(0.8, 1.2; length=9)), 3, 3)
    dsnu = reshape(collect(range(0.05, 0.45; length=9)), 3, 3)
    bad = falses(3, 3)
    bad[2, 3] = true
    background_map = reshape(collect(0.0:0.25:2.0), 3, 3)
    ipc = InterpixelCapacitance([
        0.00 0.05 0.00
        0.05 0.80 0.05
        0.00 0.05 0.00
    ])
    detector = Detector(noise=NoiseNone(), integration_time=2.0, qe=0.5,
        binning=2, gain=1.5, full_well=15.0, bits=8,
        output_type=UInt16,
        sensor=InGaAsSensor(),
        response_model=response, charge_coupling_model=ipc,
        defect_model=CompositeDetectorDefectModel(
            PixelResponseNonuniformity(prnu),
            DarkSignalNonuniformity(dsnu),
            BadPixelMask(bad; throughput=0.25)),
        nonlinearity_model=SaturatingFrameNonlinearity(0.02),
        background_map=background_map)

    expected = detector_qualification_zero_extended_convolution(
        input, Array(response.kernel))
    expected = detector_qualification_block_sum(expected, 2)
    expected .*= 0.5 * 2.0
    expected .*= prnu
    expected[bad] .*= 0.25
    expected .+= dsnu .* 2.0
    @. expected = expected / (1 + 0.02 * expected)
    clamp!(expected, 0.0, 15.0)
    expected = detector_qualification_zero_extended_convolution(
        expected, Array(ipc.response.kernel))
    expected .*= 1.5
    expected .*= 255 / 15
    clamp!(expected, 0.0, 255.0)
    expected .-= background_map
    expected_output = round.(UInt16, clamp.(expected,
        typemin(UInt16), typemax(UInt16)))
    observed = capture!(detector, input, Xoshiro(2004))
    @test observed == expected_output

    scalar_qe_detector = Detector(noise=NoiseNone(), integration_time=2.0,
        qe=0.25, response_model=NullFrameResponse())
    scalar_map = detector_test_intensity_map(fill(8.0, 2, 2);
        spectral=MonochromaticChannel(0.55e-6))
    scalar_plan = prepare_detector_acquisition(scalar_qe_detector, scalar_map)
    @test capture!(scalar_plan,
        Xoshiro(2005)) == fill(4.0, 2, 2)

    sampled_qe = SampledQuantumEfficiency(
        [0.50e-6, 0.60e-6, 0.70e-6], [0.2, 0.8, 0.4];
        out_of_band=0.1)
    for (wavelength_m, expected_qe) in (
        (0.45e-6, 0.1),
        (0.55e-6, 0.5),
        (0.65e-6, 0.6),
        (0.75e-6, 0.1),
    )
        qe_detector = Detector(noise=NoiseNone(), integration_time=2.0,
            qe=sampled_qe, response_model=NullFrameResponse())
        qe_map = detector_test_intensity_map(fill(8.0, 2, 2);
            spectral=MonochromaticChannel(wavelength_m))
        qe_plan = prepare_detector_acquisition(qe_detector, qe_map)
        @test capture!(qe_plan,
            Xoshiro(2006)) ≈ fill(16 * expected_qe, 2, 2)
    end

    low_fidelity_detector = Detector(noise=NoiseNone(),
        integration_time=1e-3, qe=1.0,
        response_model=NullFrameResponse())
    low_fidelity_map = detector_test_intensity_map(fill(1000.0, 32, 32))
    low_fidelity_plan = prepare_detector_acquisition(low_fidelity_detector,
        low_fidelity_map)
    low_fidelity_rng = Xoshiro(2007)
    @test @inferred(capture!(low_fidelity_plan, low_fidelity_rng)) isa
        Matrix{Float64}
    @test output_frame(low_fidelity_detector) == ones(32, 32)
    @test readout_products(low_fidelity_detector) isa NoFrameReadoutProducts
    low_fidelity_metadata = detector_export_metadata(low_fidelity_detector)
    @test low_fidelity_metadata.noise == :none
    @test low_fidelity_metadata.frame_response == :none
    @test low_fidelity_metadata.charge_coupling == :none
    @test low_fidelity_metadata.detector_defects == :none
    @test low_fidelity_metadata.nonlinearity_model == :none
    @test_detector_allocation prepared_detector_capture_allocations(
        low_fidelity_plan, low_fidelity_rng) == 0

    rejected_response = Detector(noise=NoiseNone(), psf_sampling=2,
        response_model=SampledFrameResponse([0.0 0.1 0.0;
            0.1 0.6 0.1; 0.0 0.1 0.0]))
    @test_throws InvalidConfiguration prepare_detector_acquisition(
        rejected_response, detector_test_intensity_map(ones(4, 4)))
end

@testset "Shared frame-detector stochastic qualification" begin
    # Each contract uses 16 independent frames of 32x32 independent pixels.
    # The six-standard-error limits are declared before any samples are drawn.
    frame_count = 16
    input = fill(50.0, 32, 32)

    shot_mean = 20.0
    shot_detector = Detector(noise=NoisePhoton(), integration_time=0.5,
        qe=0.8, response_model=NullFrameResponse())
    shot_samples = detector_qualification_samples(shot_detector, input,
        Xoshiro(2101), frame_count)
    test_detector_poisson_moments(shot_samples, shot_mean)

    dark_mean = 24.0
    dark_detector = Detector(noise=NoiseNone(), integration_time=2.0,
        qe=1.0, dark_current=12.0, response_model=NullFrameResponse())
    dark_samples = detector_qualification_samples(dark_detector,
        zeros(32, 32), Xoshiro(2102), frame_count)
    test_detector_poisson_moments(dark_samples, dark_mean)

    read_sigma = 3.0
    read_detector = Detector(noise=NoiseReadout(read_sigma),
        integration_time=1.0, qe=1.0,
        response_model=NullFrameResponse())
    read_samples = detector_qualification_samples(read_detector,
        zeros(32, 32), Xoshiro(2103), frame_count)
    test_detector_gaussian_moments(read_samples, read_sigma)
end
