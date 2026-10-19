import FFTW
import KernelAbstractions

function pixel_mtf_fixture(::Type{T}, shape) where {T}
    fy = FFTW.fftfreq(shape[1])
    fx = FFTW.fftfreq(shape[2])
    return T[sinc(1.03y) * sinc(1.03x) *
        exp(-2pi^2 * 0.09 * (y^2 + x^2)) for y in fy, x in fx]
end

function pixel_mtf_oracle(input, transfer, padding)
    padded = zeros(eltype(input), size(transfer))
    padded[(1:size(input, 1)) .+ padding, (1:size(input, 2)) .+ padding] .= input
    filtered = real.(FFTW.ifft(FFTW.fft(padded) .* transfer))
    return filtered[(1:size(input, 1)) .+ padding, (1:size(input, 2)) .+ padding]
end

pixel_mtf_allocations(owner) = @allocated apply_pixel_mtf!(owner)

@testset "Explicit pixel MTF transfer builder" begin
    width, variance, coupling = 1.02984, 0.0399599, 0.07
    response(shape; kwargs...) = pixel_mtf_transfer(shape;
        box_width_pixels=width, diffusion_variance_pixels2=variance,
        ipc_per_neighbor=coupling, kwargs...)
    for T in (Float32, Float64), shape in ((8, 8), (7, 11), (6, 9), (1, 1))
        transfer = response(shape; T)
        @test transfer isa Matrix{T}
        @test size(transfer) == shape
        @test all(isfinite, transfer)
        @test transfer[1, 1] == one(T)
        fy, fx = FFTW.fftfreq(shape[1]), FFTW.fftfreq(shape[2])
        expected = [sinc(width * y) * sinc(width * x) *
            exp(-2pi^2 * variance * (y^2 + x^2)) *
            (1 - 4coupling + 2coupling * (cos(2pi * y) + cos(2pi * x)))
            for y in fy, x in fx]
        @test transfer ≈ expected rtol=8eps(T)
        @test response(shape; T, box_width_pixels=0,
            diffusion_variance_pixels2=0, ipc_per_neighbor=0) == ones(T, shape)
    end
    # Axis and diagonal values at the even-grid Nyquist frequency.
    ipc = response((4, 6); box_width_pixels=0, diffusion_variance_pixels2=0)
    @test ipc[3, 1] ≈ 1 - 4coupling
    @test ipc[1, 4] ≈ 1 - 4coupling
    @test ipc[3, 4] ≈ 1 - 8coupling
    box = response((4, 6); diffusion_variance_pixels2=0, ipc_per_neighbor=0)
    @test box[1, 4] ≈ sinc(width / 2)
    diffusion = response((4, 6); box_width_pixels=0, ipc_per_neighbor=0)
    @test diffusion[3, 4] ≈ exp(-pi^2 * variance)
    # Signed coupling transfer is preserved instead of taking its magnitude.
    @test response((4, 6); box_width_pixels=0,
        diffusion_variance_pixels2=0, ipc_per_neighbor=0.25)[3, 4] == -1
    broad = response((4, 6); diffusion_variance_pixels2=floatmax(Float64))
    @test broad[1, 1] == 1
    @test all(isfinite, broad)

    input = reshape(collect(1.0:35.0), 5, 7)
    transfer = response(size(input) .+ 4)
    owner = prepare_pixel_mtf(input, transfer; padding=2)
    @test apply_pixel_mtf!(owner) ≈ pixel_mtf_oracle(input, transfer, 2) atol=1e-13
    identity = response(size(input) .+ 4; box_width_pixels=0,
        diffusion_variance_pixels2=0, ipc_per_neighbor=0)
    @test apply_pixel_mtf!(prepare_pixel_mtf(input, identity; padding=2)) ≈ input atol=1e-13

    @testset "Validation" begin
        for shape in ((0, 4), (-1, 4), (true, 4), (4, false), (4.0, 5),
                (4,), (4, 5, 6), [4, 5], (typemax(UInt), 4))
            @test_throws InvalidConfiguration response(shape)
        end
        for invalid in (-1, Inf, -Inf, NaN, 1 + im, "1")
            @test_throws InvalidConfiguration response((4, 6); box_width_pixels=invalid)
            @test_throws InvalidConfiguration response((4, 6); diffusion_variance_pixels2=invalid)
            @test_throws InvalidConfiguration response((4, 6); ipc_per_neighbor=invalid)
        end
        @test_throws InvalidConfiguration response((4, 6); ipc_per_neighbor=0.25001)
        for T in (Float16, Int, Float64(1), "Float64")
            @test_throws InvalidConfiguration response((4, 6); T)
        end
        @test_throws UndefKeywordError pixel_mtf_transfer((4, 6))
    end
end

@testset "Prepared calibrated pixel MTF" begin
    for T in (Float32, Float64), shape in ((8, 8), (7, 11)), pad in (0, 2, 16)
        input = T[sin(0.37i + 0.61j) + 0.2cos(0.73j) for i in 1:shape[1], j in 1:shape[2]]
        original = copy(input)
        padded_shape = shape .+ 2pad
        transfer = pixel_mtf_fixture(T, padded_shape)
        calibration = copy(transfer)
        output = similar(input)
        owner = prepare_pixel_mtf(input, transfer; padding=pad, output)
        @test owner isa Detectors.PreparedPixelMTF
        @test owner.input === input
        @test owner.products.output === output
        @test owner.plan.transfer !== transfer
        @test owner.plan.padded_shape == padded_shape
        @test owner.plan.padding == pad
        @test apply_pixel_mtf!(owner) === output
        tol = 32eps(T)
        @test output ≈ pixel_mtf_oracle(input, calibration, pad) rtol=tol atol=tol
        @test input == original
        fill!(transfer, zero(T))
        @test apply_pixel_mtf!(owner) ≈ pixel_mtf_oracle(input, calibration, pad) rtol=tol atol=tol
        first = copy(output)
        other = T[cos(0.23i - 0.71j) for i in 1:shape[1], j in 1:shape[2]]
        input .= other
        apply_pixel_mtf!(owner)
        second = copy(output)
        input .= T(0.3) .* original .+ T(1.7) .* other
        @test apply_pixel_mtf!(owner) ≈ T(0.3) .* first .+ T(1.7) .* second rtol=tol atol=tol
        apply_pixel_mtf!(owner)
        pixel_mtf_allocations(owner)
        if !coverage_instrumented()
            @test pixel_mtf_allocations(owner) == 0
        end
    end

    @testset "DC and boundaries" begin
        input = fill(2.0, 9, 7)
        mtf = pixel_mtf_fixture(Float64, size(input))
        @test apply_pixel_mtf!(prepare_pixel_mtf(input, mtf; padding=0)) ≈ input atol=1e-14
        # Unit transfer is identity even when the input is surrounded by zeros.
        @test apply_pixel_mtf!(prepare_pixel_mtf(input, ones(41, 39))) ≈ input atol=1e-14
        input .= 0
        input[1, 1] = 1
        padded_mtf = pixel_mtf_fixture(Float64, (41, 39))
        @test apply_pixel_mtf!(prepare_pixel_mtf(input, padded_mtf)) ≈
            pixel_mtf_oracle(input, padded_mtf, 16) atol=1e-14
        # SubArray and transposed one-based input preserve their actual sample order.
        storage = reshape(collect(1.0:120.0), 12, 10)
        for wrapped in (view(storage, 2:8, 3:8), transpose(storage))
            transfer = ones(size(wrapped) .+ 4)
            owner = prepare_pixel_mtf(wrapped, transfer; padding=2)
            @test apply_pixel_mtf!(owner) ≈ wrapped atol=1e-13
        end
    end

    @testset "Validation" begin
        input = zeros(Float64, 4, 5)
        transfer = ones(8, 9)
        @test_throws InvalidConfiguration prepare_pixel_mtf(input, transfer; padding=-1)
        @test_throws InvalidConfiguration prepare_pixel_mtf(input, transfer; padding=2.5)
        @test_throws InvalidConfiguration prepare_pixel_mtf(input, transfer; padding=true)
        @test_throws DimensionMismatchError prepare_pixel_mtf(input, ones(8, 8); padding=2)
        @test_throws DimensionMismatchError prepare_pixel_mtf(input, transfer; padding=2, output=zeros(5, 4))
        @test_throws InvalidConfiguration prepare_pixel_mtf(input, transfer; padding=2, output=input)
        @test_throws InvalidConfiguration prepare_pixel_mtf(input, transfer; padding=2, output=zeros(Float32, 4, 5))
        @test_throws InvalidConfiguration prepare_pixel_mtf(zeros(Int, 4, 5), transfer; padding=2)
        @test_throws InvalidConfiguration prepare_pixel_mtf(input, complex.(transfer); padding=2)
        for invalid in (NaN, Inf, -Inf)
            bad = copy(transfer)
            bad[2, 3] = invalid
            @test_throws InvalidConfiguration prepare_pixel_mtf(input, bad; padding=2)
        end
        bad_dc = copy(transfer)
        bad_dc[1, 1] = 0.9
        @test_throws InvalidConfiguration prepare_pixel_mtf(input, bad_dc; padding=2)
        oversized = copy(transfer)
        oversized[2, 3] = floatmax(Float64)
        @test_throws InvalidConfiguration prepare_pixel_mtf(zeros(Float32, 4, 5), oversized; padding=2)
        @test_throws InvalidConfiguration prepare_pixel_mtf(input, input; padding=0)
        output = ones(size(input))
        @test_throws InvalidConfiguration prepare_pixel_mtf(input, output; padding=0, output)
        @test_throws InvalidConfiguration prepare_pixel_mtf(zeros(0, 5), ones(4, 9); padding=2)
    end

    @testset "Accelerator kernel parity on KA CPU" begin
        input = reshape(collect(1.0:77.0), 7, 11)
        transfer = pixel_mtf_fixture(Float64, (39, 43))
        owner = prepare_pixel_mtf(input, transfer)
        style = Backends.AcceleratorStyle(KernelAbstractions.CPU())
        plan, workspace = owner.plan, owner.workspace
        Detectors._pixel_mtf_pad!(style, workspace.buffer, input, plan)
        KernelAbstractions.synchronize(style.backend)
        Backends.enqueue_fft_plan!(workspace.buffer, workspace.fft_plan)
        Detectors._pixel_mtf_multiply!(style, workspace.buffer, plan.transfer, plan)
        KernelAbstractions.synchronize(style.backend)
        Backends.enqueue_fft_plan!(workspace.buffer, workspace.ifft_plan)
        Detectors._pixel_mtf_crop!(style, owner.products.output, workspace.buffer, plan)
        KernelAbstractions.synchronize(style.backend)
        @test owner.products.output ≈ pixel_mtf_oracle(input, transfer, 16) atol=1e-13
    end
end

pixel_mtf_graph_allocations(graph) = @allocated step_graph!(graph)

@testset "Pixel MTF graph adapter" begin
    for T in (Float32, Float64), shape in ((6, 8), (5, 7)), padding in (0, 2)
        input = T[sin(0.31i + 0.59j) for i in 1:shape[1], j in 1:shape[2]]
        transfer = pixel_mtf_transfer(shape .+ 2padding;
            box_width_pixels=1.03, diffusion_variance_pixels2=0.09,
            ipc_per_neighbor=0.03, T)
        calibration = copy(transfer)
        node = pixel_mtf_node(:response; input_shape=shape, padding, T,
            input_schema="test.photon-rate/1", output_schema="test.photon-rate/1",
            transfer)
        @test node.config.transfer !== transfer
        @test node.config.transfer == calibration
        fill!(transfer, T(NaN))
        graph = prepare_algorithm_graph(algorithm_graph((node,);
            inputs=(graph_input(:signal, :response => :input, input),),
            outputs=(graph_output(:filtered, :response => :output),)))
        owner = AlgorithmGraphs.prepared_graph_node(graph, Val(:response))
        output = graph_output(graph, Val(:filtered))
        @test owner.prepared.input === input
        @test owner.prepared.products.output === output
        @test owner.prepared.plan.transfer !== node.config.transfer
        @test AlgorithmGraphs.graph_node_capture_capability(owner) isa
            AlgorithmGraphs.GraphNodeCaptureUnsupported
        @test step_graph!(graph) === graph
        tolerance = 32eps(T)
        @test output ≈ pixel_mtf_oracle(input, calibration, padding) rtol=tolerance atol=tolerance
        input .*= T(0.37)
        step_graph!(graph)
        @test output ≈ pixel_mtf_oracle(input, calibration, padding) rtol=tolerance atol=tolerance
        pixel_mtf_graph_allocations(graph)
        if !coverage_instrumented()
            @test pixel_mtf_graph_allocations(graph) == 0
        end
        @test reset_graph!(graph) === graph
        @test all(iszero, output)
        @test all(iszero, owner.prepared.workspace.buffer)
        @test owner.prepared.products.output === output
        step_graph!(graph)
        @test output ≈ pixel_mtf_oracle(input, calibration, padding) rtol=tolerance atol=tolerance
    end

    declare(; kwargs...) = pixel_mtf_node(:response; input_shape=(4, 6),
        padding=0, T=Float64, input_schema="test.signal/1",
        output_schema="test.signal/1", transfer=ones(4, 6), kwargs...)
    @test_throws AlgorithmGraphError declare(input_schema="")
    @test_throws AlgorithmGraphError declare(output_schema="")
    @test_throws AlgorithmGraphError declare(padding=true)
    @test_throws AlgorithmGraphError declare(padding=-1)
    @test_throws AlgorithmGraphError declare(T=Float16)
    @test_throws AlgorithmGraphError declare(transfer=ones(4, 5))
    @test_throws AlgorithmGraphError declare(transfer=complex.(ones(4, 6)))
    @test_throws AlgorithmGraphError declare(transfer=fill(NaN, 4, 6))
    @test_throws AlgorithmGraphError declare(transfer=fill(0.9, 4, 6))
end
