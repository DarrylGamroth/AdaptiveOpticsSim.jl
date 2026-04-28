using KernelAbstractions

const KA_CPU_STYLE = AdaptiveOpticsSim.AcceleratorStyle(KernelAbstractions.CPU())
const SCALAR_CPU_STYLE = AdaptiveOpticsSim.ScalarCPUStyle()
const KA_CPU_EXERCISED_KERNELS = Set{Symbol}()

function mark_ka_cpu_kernel!(names::Symbol...)
    union!(KA_CPU_EXERCISED_KERNELS, names)
    return nothing
end

function source_kernel_names()
    src_root = normpath(joinpath(@__DIR__, "..", "src"))
    rx = r"^@kernel function\s+([A-Za-z_]\w*!)"
    names = Set{Symbol}()
    for (root, _, files) in walkdir(src_root)
        for file in files
            endswith(file, ".jl") || continue
            for line in eachline(joinpath(root, file))
                m = match(rx, line)
                m === nothing || push!(names, Symbol(m.captures[1]))
            end
        end
    end
    return names
end

function ka_cpu_close(actual, expected; atol=1e-12, rtol=1e-12)
    return actual == expected || isapprox(actual, expected; atol=atol, rtol=rtol)
end

@testset "KA CPU Matrix" begin
    @testset "Core kernels" begin
        src = reshape(collect(1.0:16.0), 4, 4)

        scalar_fftshift = similar(src)
        ka_fftshift = similar(src)
        AdaptiveOpticsSim._fftshift2d!(SCALAR_CPU_STYLE, scalar_fftshift, src)
        AdaptiveOpticsSim._fftshift2d!(KA_CPU_STYLE, ka_fftshift, src)
        mark_ka_cpu_kernel!(:fftshift2d_kernel!)
        @test ka_fftshift == scalar_fftshift

        scalar_shift = similar(src)
        ka_shift = similar(src)
        AdaptiveOpticsSim._circshift2d!(SCALAR_CPU_STYLE, scalar_shift, src, (1, -1))
        AdaptiveOpticsSim._circshift2d!(KA_CPU_STYLE, ka_shift, src, (1, -1))
        mark_ka_cpu_kernel!(:circshift2d_kernel!)
        @test ka_shift == scalar_shift

        scalar_bin = Matrix{Float64}(undef, 2, 2)
        ka_bin = similar(scalar_bin)
        AdaptiveOpticsSim._bin2d!(SCALAR_CPU_STYLE, scalar_bin, src, 2)
        AdaptiveOpticsSim._bin2d!(KA_CPU_STYLE, ka_bin, src, 2)
        mark_ka_cpu_kernel!(:bin2d_kernel!)
        @test ka_bin == scalar_bin

        scalar_freqs = Vector{Float64}(undef, 8)
        ka_freqs = similar(scalar_freqs)
        AdaptiveOpticsSim._fftfreq!(SCALAR_CPU_STYLE, scalar_freqs, 8, 0.25, 0.0)
        AdaptiveOpticsSim._fftfreq!(KA_CPU_STYLE, ka_freqs, 8, 0.25, 0.0)
        mark_ka_cpu_kernel!(:fftfreq_kernel!)
        @test ka_freqs == scalar_freqs

        rng1 = MersenneTwister(11)
        rng2 = MersenneTwister(11)
        ka_randn_a = zeros(Float64, 8)
        ka_randn_b = zeros(Float64, 8)
        AdaptiveOpticsSim.randn_backend_async!(KA_CPU_STYLE, rng1, ka_randn_a)
        KernelAbstractions.synchronize(KA_CPU_STYLE.backend)
        AdaptiveOpticsSim.randn_backend_async!(KA_CPU_STYLE, rng2, ka_randn_b)
        KernelAbstractions.synchronize(KA_CPU_STYLE.backend)
        mark_ka_cpu_kernel!(:randn_fill_kernel!)
        @test ka_randn_a == ka_randn_b
        @test all(isfinite, ka_randn_a)

        rng3 = MersenneTwister(12)
        rng4 = MersenneTwister(12)
        ka_poisson_a = fill(4.0, 8)
        ka_poisson_b = fill(4.0, 8)
        AdaptiveOpticsSim.poisson_noise_async!(KA_CPU_STYLE, rng3, ka_poisson_a)
        KernelAbstractions.synchronize(KA_CPU_STYLE.backend)
        AdaptiveOpticsSim.poisson_noise_async!(KA_CPU_STYLE, rng4, ka_poisson_b)
        KernelAbstractions.synchronize(KA_CPU_STYLE.backend)
        mark_ka_cpu_kernel!(:poisson_noise_kernel!)
        @test ka_poisson_a == ka_poisson_b
        @test all(x -> x >= 0, ka_poisson_a)
    end

    @testset "Telescope and pupil kernels" begin
        params = TelescopeParams{Float64}(16, 8.0, 1e-3, 0.15, 0.0)
        scalar_pupil = Matrix{Bool}(undef, 16, 16)
        ka_pupil = similar(scalar_pupil)
        AdaptiveOpticsSim._generate_pupil!(SCALAR_CPU_STYLE, scalar_pupil, params)
        AdaptiveOpticsSim._generate_pupil!(KA_CPU_STYLE, ka_pupil, params)
        mark_ka_cpu_kernel!(:radial_mask_kernel!)
        @test ka_pupil == scalar_pupil

        scalar_spiders = copy(scalar_pupil)
        ka_spiders = copy(ka_pupil)
        angles = [0.0, 45.0, 90.0]
        AdaptiveOpticsSim._apply_spiders!(SCALAR_CPU_STYLE, scalar_spiders, angles, 0.1, 0.0, 0.0, 8.5, 8.5, 8.0, 16)
        AdaptiveOpticsSim._apply_spiders!(KA_CPU_STYLE, ka_spiders, angles, 0.1, 0.0, 0.0, 8.5, 8.5, 8.0, 16)
        mark_ka_cpu_kernel!(:spider_mask_kernel!)
        @test ka_spiders == scalar_spiders

        scalar_roi = falses(8, 8)
        ka_roi = similar(scalar_roi)
        roi = AdaptiveOpticsSim.RectangularROI(2:6, 3:7)
        AdaptiveOpticsSim._build_rectangular_roi!(SCALAR_CPU_STYLE, scalar_roi, roi, true, false)
        AdaptiveOpticsSim._build_rectangular_roi!(KA_CPU_STYLE, ka_roi, roi, true, false)
        mark_ka_cpu_kernel!(:rectangular_roi_kernel!)
        @test ka_roi == scalar_roi
    end

    @testset "Wavefront support kernels" begin
        pupil = trues(8, 8)
        pupil[1:2, :] .= false
        threshold = 0.4
        scalar_valid = Matrix{Bool}(undef, 2, 2)
        ka_valid = similar(scalar_valid)
        AdaptiveOpticsSim._set_valid_subapertures!(SCALAR_CPU_STYLE, scalar_valid, pupil, threshold, 4, 2)
        AdaptiveOpticsSim._set_valid_subapertures!(KA_CPU_STYLE, ka_valid, pupil, threshold, 4, 2)
        mark_ka_cpu_kernel!(:subaperture_grid_mask_kernel!)
        @test ka_valid == scalar_valid

        opd = reshape(collect(1.0:64.0), 8, 8)
        scalar_slopes = Vector{Float64}(undef, 8)
        ka_slopes = similar(scalar_slopes)
        AdaptiveOpticsSim._geometric_slopes!(SCALAR_CPU_STYLE, scalar_slopes, opd, scalar_valid, 4, 2, 4)
        AdaptiveOpticsSim._geometric_slopes!(KA_CPU_STYLE, ka_slopes, opd, ka_valid, 4, 2, 4)
        mark_ka_cpu_kernel!(:geometric_slopes_kernel!)
        @test ka_slopes == scalar_slopes

        edge_mask = falses(8, 8)
        edge_mask[:, 1] .= true
        edge_mask[:, end] .= true
        edge_mask[1, :] .= true
        edge_mask[end, :] .= true
        scalar_edge_slopes = Vector{Float64}(undef, 8)
        ka_edge_slopes = similar(scalar_edge_slopes)
        AdaptiveOpticsSim._edge_geometric_slopes!(SCALAR_CPU_STYLE, scalar_edge_slopes, opd, scalar_valid, edge_mask, 4, 2, 4)
        AdaptiveOpticsSim._edge_geometric_slopes!(KA_CPU_STYLE, ka_edge_slopes, opd, ka_valid, edge_mask, 4, 2, 4)
        mark_ka_cpu_kernel!(:edge_geometric_slopes_kernel!)
        @test ka_edge_slopes == scalar_edge_slopes

        spot_cube = reshape(collect(1.0:36.0), 4, 3, 3)
        scalar_image = Matrix{Float64}(undef, 7, 7)
        ka_image = similar(scalar_image)
        AdaptiveOpticsSim._shack_hartmann_detector_image!(
            SCALAR_CPU_STYLE, SCALAR_CPU_STYLE, scalar_image, spot_cube, 2, 3, 3, 1, -1.0)
        AdaptiveOpticsSim._shack_hartmann_detector_image!(
            KA_CPU_STYLE, KA_CPU_STYLE, ka_image, spot_cube, 2, 3, 3, 1, -1.0)
        mark_ka_cpu_kernel!(:shack_hartmann_detector_image_kernel!)
        @test ka_image == scalar_image
    end

    @testset "Grouped and Shack-Hartmann stack kernels" begin
        stack = reshape(collect(1.0:48.0), 4, 4, 3)
        scalar_group = Matrix{Float64}(undef, 4, 4)
        ka_group = similar(scalar_group)
        AdaptiveOpticsSim.reduce_grouped_stack!(SCALAR_CPU_STYLE, scalar_group, stack, 3)
        AdaptiveOpticsSim.reduce_grouped_stack!(KA_CPU_STYLE, ka_group, stack, 3)
        mark_ka_cpu_kernel!(:reduce_grouped_stack_kernel!)
        @test ka_group == scalar_group

        block_stack = reshape(collect(1.0:128.0), 4, 4, 8)
        block_out = zeros(Float64, 4, 4, 4)
        AdaptiveOpticsSim.reduce_grouped_blocks!(KA_CPU_STYLE, block_out, block_stack, 4, 2)
        mark_ka_cpu_kernel!(:reduce_grouped_blocks_kernel!)
        @test block_out == block_stack[:, :, 1:4] .+ block_stack[:, :, 5:8]

        @test size(AdaptiveOpticsSim.grouped_stack_view(stack, 2), 3) == 2
        grouped_tel = Telescope(resolution=8, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
        grouped_wfs = ShackHartmannWFS(grouped_tel; n_lenslets=2, mode=Diffractive(), n_pix_subap=2)
        grouped_sources = [1.0, 2.0, 3.0]
        function fill_grouped_stage!(dest, scale, value)
            fill!(dest, scale * value)
            return dest
        end
        scalar_accum = zeros(Float64, 2, 2)
        scalar_stack = zeros(Float64, 2, 2, length(grouped_sources))
        AdaptiveOpticsSim.accumulate_grouped_sources!(
            AdaptiveOpticsSim.GroupedStackReducePlan(), SCALAR_CPU_STYLE, grouped_wfs,
            scalar_accum, scalar_stack, grouped_sources, fill_grouped_stage!, 2.0)
        @test scalar_accum == fill(12.0, 2, 2)
        ka_accum = zeros(Float64, 2, 2)
        ka_stack = zeros(Float64, 2, 2, length(grouped_sources))
        AdaptiveOpticsSim.accumulate_grouped_sources!(
            AdaptiveOpticsSim.GroupedStackReducePlan(), KA_CPU_STYLE, grouped_wfs,
            ka_accum, ka_stack, grouped_sources, fill_grouped_stage!, 2.0)
        @test ka_accum == scalar_accum
        staged_scalar = zeros(Float64, 2, 2)
        staged_scalar_stack = zeros(Float64, 2, 2, length(grouped_sources))
        AdaptiveOpticsSim.accumulate_grouped_sources!(
            AdaptiveOpticsSim.GroupedStaged2DPlan(), SCALAR_CPU_STYLE, grouped_wfs,
            staged_scalar, staged_scalar_stack, grouped_sources, fill_grouped_stage!, 3.0)
        @test staged_scalar == fill(18.0, 2, 2)
        staged_ka = zeros(Float64, 2, 2)
        staged_ka_stack = zeros(Float64, 2, 2, length(grouped_sources))
        AdaptiveOpticsSim.accumulate_grouped_sources!(
            AdaptiveOpticsSim.GroupedStaged2DPlan(), KA_CPU_STYLE, grouped_wfs,
            staged_ka, staged_ka_stack, grouped_sources, fill_grouped_stage!, 3.0)
        @test staged_ka == staged_scalar
        dispatch_accum = zeros(Float64, 2, 2)
        dispatch_stack = zeros(Float64, 2, 2, length(grouped_sources))
        @test AdaptiveOpticsSim.grouped_accumulation_plan(KA_CPU_STYLE, grouped_wfs) isa
              AdaptiveOpticsSim.GroupedStackReducePlan
        AdaptiveOpticsSim.accumulate_grouped_sources!(
            KA_CPU_STYLE, grouped_wfs, dispatch_accum, dispatch_stack,
            grouped_sources, fill_grouped_stage!, 1.0)
        @test dispatch_accum == fill(6.0, 2, 2)

        mask = Bool[1 0; 1 1]
        values = [1.0 2.0; 3.0 4.0]
        values_view = @view values[1:2, 1:2]
        @test AdaptiveOpticsSim.masked_sum2d(values, mask) == 8.0
        scalar_buffer = zeros(Float64, 1)
        scalar_host = zeros(Float64, 1)
        host_parent = zeros(Float64, size(values)...)
        @test first(AdaptiveOpticsSim.masked_sum2d(KA_CPU_STYLE, values_view, mask, mask,
            scalar_buffer, scalar_host, host_parent)) == 8.0
        @test AdaptiveOpticsSim.backend_maximum_value(KA_CPU_STYLE, values) == 4.0
        signal = collect(1.0:8.0)
        @test AdaptiveOpticsSim.packed_valid_pair_mean(KA_CPU_STYLE, signal, mask) ==
              AdaptiveOpticsSim.packed_valid_pair_mean(AdaptiveOpticsSim.DirectReductionPlan(), signal, mask)
        @test AdaptiveOpticsSim.packed_valid_pair_mean(AdaptiveOpticsSim.HostMirrorReductionPlan(), signal, mask) ==
              AdaptiveOpticsSim.packed_valid_pair_mean(SCALAR_CPU_STYLE, signal, mask)

        n_sub = 2
        sub = 2
        n = 4
        pad = 4
        ox = 1
        oy = 1
        n_spots = n_sub * n_sub
        valid_mask = Bool[1 0; 1 1]
        pupil = trues(n, n)
        opd = zeros(Float64, n, n)
        phasor = ones(ComplexF64, pad, pad)
        fft_stack = fill(99.0 + 0.0im, pad, pad, n_spots)
        AdaptiveOpticsSim.launch_kernel!(KA_CPU_STYLE, AdaptiveOpticsSim.sh_field_stack_kernel!,
            fft_stack, valid_mask, pupil, opd, phasor, 2.0, 0.0, n_sub, sub, ox, oy, n, pad;
            ndrange=size(fft_stack))
        mark_ka_cpu_kernel!(:sh_field_stack_kernel!)
        @test fft_stack[2, 2, 1] == 2.0 + 0.0im
        @test fft_stack[2, 2, 2] == 0.0 + 0.0im
        @test fft_stack[3, 3, 4] == 2.0 + 0.0im
        @test fft_stack[1, 1, 1] == 0.0 + 0.0im

        intensity_stack = zeros(Float64, pad, pad, n_spots)
        AdaptiveOpticsSim.launch_kernel!(KA_CPU_STYLE, AdaptiveOpticsSim.complex_abs2_stack_kernel!,
            intensity_stack, fft_stack, pad, n_spots; ndrange=size(intensity_stack))
        mark_ka_cpu_kernel!(:complex_abs2_stack_kernel!)
        @test intensity_stack[2, 2, 1] == 4.0
        @test intensity_stack[2, 2, 2] == 0.0

        amp_scales = [2.0, 3.0]
        opd_to_cycles = [0.0, 0.0]
        fft_ast = fill(99.0 + 0.0im, pad, pad, 2 * n_spots)
        AdaptiveOpticsSim.launch_kernel!(KA_CPU_STYLE, AdaptiveOpticsSim.sh_field_asterism_stack_kernel!,
            fft_ast, valid_mask, pupil, opd, phasor, amp_scales, opd_to_cycles, n_sub, sub, ox, oy, n, pad,
            n_spots; ndrange=size(fft_ast))
        mark_ka_cpu_kernel!(:sh_field_asterism_stack_kernel!)
        @test fft_ast[2, 2, 1] == 2.0 + 0.0im
        @test fft_ast[2, 2, n_spots + 1] == 3.0 + 0.0im
        @test fft_ast[2, 2, 2] == 0.0 + 0.0im

        sample_input = reshape(collect(1.0:64.0), 4, 4, n_spots)
        sampled = fill(-1.0, n_spots, 2, 2)
        AdaptiveOpticsSim.launch_kernel!(KA_CPU_STYLE, AdaptiveOpticsSim.sh_sample_spot_stack_kernel!,
            sampled, sample_input, valid_mask, 2, n_sub, 2, 2; ndrange=size(sampled))
        mark_ka_cpu_kernel!(:sh_sample_spot_stack_kernel!)
        @test sampled[1, 1, 1] == sum(sample_input[1:2, 1:2, 1])
        @test sampled[2, 1, 1] == 0.0
        @test sampled[4, 2, 2] == sum(sample_input[3:4, 3:4, 4])

        spot_cube = zeros(Float64, n_spots, 3, 3)
        spot_cube[1, 2, 3] = 10.0
        spot_cube[1, 1, 1] = 1.0
        spot_cube[3, 3, 1] = 5.0
        spot_cube[4, 2, 2] = 4.0
        stats = zeros(Float64, 3 * n_spots)
        AdaptiveOpticsSim.launch_kernel!(KA_CPU_STYLE, AdaptiveOpticsSim.sh_spot_centroid_stats_kernel!,
            stats, spot_cube, valid_mask, 0.5, n_sub, 3, 3; ndrange=n_spots)
        mark_ka_cpu_kernel!(:sh_spot_centroid_stats_kernel!)
        @test stats[1] == 10.0
        @test stats[2] == 1.0
        @test stats[3] == 2.0
        @test stats[4:6] == zeros(3)

        accum = ones(Float64, 3 * n_spots)
        AdaptiveOpticsSim.launch_kernel!(KA_CPU_STYLE, AdaptiveOpticsSim.accumulate_spot_stats_kernel!,
            accum, stats, n_spots; ndrange=3 * n_spots)
        mark_ka_cpu_kernel!(:accumulate_spot_stats_kernel!)
        @test accum == stats .+ 1.0

        slopes = fill(-1.0, 2 * n_spots)
        reference = zeros(Float64, 2 * n_spots)
        AdaptiveOpticsSim.launch_kernel!(KA_CPU_STYLE, AdaptiveOpticsSim.sh_finalize_asterism_slopes_kernel!,
            slopes, accum, reference, valid_mask, 2.0, 2, n_sub, n_spots; ndrange=n_spots)
        mark_ka_cpu_kernel!(:sh_finalize_asterism_slopes_kernel!)
        @test slopes[1] == (accum[3] / 2) / 2
        @test slopes[n_spots + 1] == (accum[2] / 2) / 2
        @test slopes[2] == 0.0
        @test slopes[n_spots + 2] == 0.0

        centroid_slopes = fill(-1.0, 2 * n_spots)
        centroid_cube = copy(spot_cube)
        AdaptiveOpticsSim.launch_kernel!(KA_CPU_STYLE, AdaptiveOpticsSim.sh_spot_centroid_kernel!,
            centroid_slopes, centroid_cube, valid_mask, 2.0, n_sub, n_spots, 3, 3; ndrange=n_spots)
        mark_ka_cpu_kernel!(:sh_spot_centroid_kernel!)
        @test centroid_slopes[1] == 2.0
        @test centroid_slopes[n_spots + 1] == 1.0
        @test centroid_slopes[2] == 0.0
        @test centroid_cube[1, 1, 1] == 0.0

        ref_scaled_slopes = fill(-1.0, 2 * n_spots)
        reference .= 0.5
        AdaptiveOpticsSim.launch_kernel!(KA_CPU_STYLE, AdaptiveOpticsSim.sh_spot_centroid_reference_scale_kernel!,
            ref_scaled_slopes, copy(spot_cube), reference, valid_mask, 2.0, 2.0, n_sub, n_spots, 3, 3;
            ndrange=n_spots)
        mark_ka_cpu_kernel!(:sh_spot_centroid_reference_scale_kernel!)
        @test ref_scaled_slopes[1] == (2.0 - 0.5) / 2
        @test ref_scaled_slopes[n_spots + 1] == (1.0 - 0.5) / 2
        @test ref_scaled_slopes[2] == 0.0

        cutoff_stats = zeros(Float64, 3 * n_spots)
        AdaptiveOpticsSim.launch_kernel!(KA_CPU_STYLE, AdaptiveOpticsSim.sh_spot_cutoff_stats_kernel!,
            cutoff_stats, spot_cube, valid_mask, 2.0, n_sub, 3, 3; ndrange=n_spots)
        mark_ka_cpu_kernel!(:sh_spot_cutoff_stats_kernel!)
        @test cutoff_stats[1:3] == [10.0, 10.0, 20.0]

        simple_slopes = fill(-1.0, 2 * n_spots)
        AdaptiveOpticsSim.launch_kernel!(KA_CPU_STYLE, AdaptiveOpticsSim.sh_finalize_spot_slopes_kernel!,
            simple_slopes, cutoff_stats, valid_mask, n_sub, n_spots; ndrange=n_spots)
        mark_ka_cpu_kernel!(:sh_finalize_spot_slopes_kernel!)
        @test simple_slopes[1] == 2.0
        @test simple_slopes[n_spots + 1] == 1.0

        scaled_slopes = fill(-1.0, 2 * n_spots)
        AdaptiveOpticsSim.launch_kernel!(KA_CPU_STYLE, AdaptiveOpticsSim.sh_finalize_spot_slopes_reference_scale_kernel!,
            scaled_slopes, cutoff_stats, reference, valid_mask, 2.0, n_sub, n_spots; ndrange=n_spots)
        mark_ka_cpu_kernel!(:sh_finalize_spot_slopes_reference_scale_kernel!)
        @test scaled_slopes[1] == (2.0 - 0.5) / 2
        @test scaled_slopes[n_spots + 1] == (1.0 - 0.5) / 2

        invalid_cube = ones(Float64, n_spots, 2, 2)
        AdaptiveOpticsSim.launch_kernel!(KA_CPU_STYLE, AdaptiveOpticsSim.zero_invalid_spots_kernel!,
            invalid_cube, valid_mask, n_sub, 2, 2; ndrange=size(invalid_cube))
        mark_ka_cpu_kernel!(:zero_invalid_spots_kernel!)
        @test all(iszero, invalid_cube[2, :, :])
        @test all(==(1.0), invalid_cube[1, :, :])

        invalid_slopes = ones(Float64, 2 * n_spots)
        AdaptiveOpticsSim.launch_kernel!(KA_CPU_STYLE, AdaptiveOpticsSim.zero_invalid_sh_slopes_kernel!,
            invalid_slopes, valid_mask, n_sub, n_spots; ndrange=n_spots)
        mark_ka_cpu_kernel!(:zero_invalid_sh_slopes_kernel!)
        @test invalid_slopes[2] == 0.0
        @test invalid_slopes[n_spots + 2] == 0.0
        @test invalid_slopes[1] == 1.0
    end

    @testset "Atmosphere kernels" begin
        freqs = collect(range(0.0, 1.0; length=8))
        scalar_psd = Matrix{Float64}(undef, 8, 8)
        ka_psd = similar(scalar_psd)
        AdaptiveOpticsSim.update_psd!(SCALAR_CPU_STYLE, scalar_psd, freqs, 0.02, 4pi^2, 0.01, -11 / 6, 8)
        AdaptiveOpticsSim.update_psd!(KA_CPU_STYLE, ka_psd, freqs, 0.02, 4pi^2, 0.01, -11 / 6, 8)
        mark_ka_cpu_kernel!(:kolmogorov_psd_kernel!)
        @test ka_cpu_close(ka_psd, scalar_psd)

        scalar_phase_psd = similar(scalar_psd)
        ka_phase_psd = similar(scalar_psd)
        AdaptiveOpticsSim._fill_phase_psd!(SCALAR_CPU_STYLE, scalar_phase_psd, freqs, 0.02, 4pi^2, 0.01, -11 / 6, 0.0, 8)
        AdaptiveOpticsSim._fill_phase_psd!(KA_CPU_STYLE, ka_phase_psd, freqs, 0.02, 4pi^2, 0.01, -11 / 6, 0.0, 8)
        mark_ka_cpu_kernel!(:phase_screen_psd_kernel!)
        @test ka_cpu_close(ka_phase_psd, scalar_phase_psd)

        scalar_spectrum = AdaptiveOpticsSim.phase_spectrum(freqs, 0.15, 25.0)
        ka_spectrum_out = similar(freqs)
        T = eltype(freqs)
        coeff = T(0.023) * T(0.15)^(-T(5) / T(3))
        AdaptiveOpticsSim._phase_spectrum!(KA_CPU_STYLE, ka_spectrum_out, freqs, coeff, T(2pi)^2, inv(T(25.0))^2, -T(11) / T(6))
        mark_ka_cpu_kernel!(:phase_spectrum_kernel!)
        @test ka_cpu_close(ka_spectrum_out, scalar_spectrum)

        rho = reshape(collect(range(0.0, 1.0; length=16)), 4, 4)
        scalar_cov = AdaptiveOpticsSim.phase_covariance(rho, 0.15, 25.0)
        ka_cov = AdaptiveOpticsSim._phase_covariance(KA_CPU_STYLE, rho, 0.15, 25.0)
        mark_ka_cpu_kernel!(:phase_covariance_kernel!)
        @test ka_cpu_close(ka_cov, scalar_cov; rtol=1e-8, atol=1e-8)

        screen = reshape(collect(1.0:100.0), 10, 10)
        scalar_extract = Matrix{Float64}(undef, 4, 4)
        ka_extract = similar(scalar_extract)
        AdaptiveOpticsSim.extract_shifted_screen!(scalar_extract, screen, 0.25, -0.5, 0.75, 1.0)
        AdaptiveOpticsSim._extract_shifted_screen!(KA_CPU_STYLE, ka_extract, screen, 3.75, 4.5, 1.0, 0.75, 4, 10)
        mark_ka_cpu_kernel!(:moving_layer_extract_kernel!)
        @test ka_cpu_close(ka_extract, scalar_extract)

        scalar_accum = fill(2.0, 4, 4)
        ka_accum = fill(2.0, 4, 4)
        AdaptiveOpticsSim.accumulate_shifted_screen!(scalar_accum, screen, 0.25, -0.5, 0.75, 1.0)
        AdaptiveOpticsSim._accumulate_shifted_screen!(KA_CPU_STYLE, ka_accum, screen, 3.75, 4.5, 1.0, 0.75, 4, 10)
        mark_ka_cpu_kernel!(:moving_layer_accumulate_kernel!)
        @test ka_cpu_close(ka_accum, scalar_accum)

        scalar_sub = zeros(Float64, 8, 8)
        ka_sub = zeros(Float64, 8, 8)
        AdaptiveOpticsSim._add_subharmonics!(SCALAR_CPU_STYLE, scalar_sub, 0.15, 25.0, 0.2, 1e-10;
            rng=MersenneTwister(2), n_levels=1, radius=1)
        AdaptiveOpticsSim._add_subharmonics!(KA_CPU_STYLE, ka_sub, 0.15, 25.0, 0.2, 1e-10;
            rng=MersenneTwister(2), n_levels=1, radius=1)
        mark_ka_cpu_kernel!(:add_subharmonics_kernel!)
        @test ka_cpu_close(ka_sub, scalar_sub; rtol=1e-10, atol=1e-10)
    end

    @testset "Optics kernels" begin
        tel = Telescope(resolution=8, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
        src = Source(band=:I, magnitude=0.0)
        scalar_field = Matrix{ComplexF64}(undef, 16, 16)
        ka_field = similar(scalar_field)
        AdaptiveOpticsSim._fill_telescope_field!(SCALAR_CPU_STYLE, scalar_field, tel, src, 2, true)
        AdaptiveOpticsSim._fill_telescope_field!(KA_CPU_STYLE, ka_field, tel, src, 2, true)
        mark_ka_cpu_kernel!(:fill_telescope_field_kernel!)
        @test ka_cpu_close(ka_field, scalar_field)

        scalar_ef = ElectricField(tel, src; zero_padding=2)
        ka_ef = ElectricField(tel, src; zero_padding=2)
        opd = reshape(collect(range(0.0, 1e-8; length=64)), 8, 8)
        phase = reshape(collect(range(0.0, 0.1; length=64)), 8, 8)
        amp = reshape(collect(range(0.9, 1.1; length=64)), 8, 8)
        AdaptiveOpticsSim._apply_phase!(SCALAR_CPU_STYLE, scalar_ef, opd, :opd)
        AdaptiveOpticsSim._apply_phase!(KA_CPU_STYLE, ka_ef, opd, :opd)
        mark_ka_cpu_kernel!(:apply_phase_opd_kernel!)
        @test ka_cpu_close(ka_ef.state.field, scalar_ef.state.field)
        AdaptiveOpticsSim._apply_phase!(SCALAR_CPU_STYLE, scalar_ef, phase, :phase)
        AdaptiveOpticsSim._apply_phase!(KA_CPU_STYLE, ka_ef, phase, :phase)
        mark_ka_cpu_kernel!(:apply_phase_rad_kernel!)
        @test ka_cpu_close(ka_ef.state.field, scalar_ef.state.field)
        AdaptiveOpticsSim._apply_amplitude!(SCALAR_CPU_STYLE, scalar_ef, amp)
        AdaptiveOpticsSim._apply_amplitude!(KA_CPU_STYLE, ka_ef, amp)
        mark_ka_cpu_kernel!(:apply_amplitude_kernel!)
        @test ka_cpu_close(ka_ef.state.field, scalar_ef.state.field)

        scalar_intensity = similar(scalar_ef.state.intensity)
        ka_intensity = similar(ka_ef.state.intensity)
        AdaptiveOpticsSim._intensity!(SCALAR_CPU_STYLE, scalar_intensity, scalar_ef.state.field)
        AdaptiveOpticsSim._intensity!(KA_CPU_STYLE, ka_intensity, ka_ef.state.field)
        mark_ka_cpu_kernel!(:intensity_kernel!)
        @test ka_cpu_close(ka_intensity, scalar_intensity)
        AdaptiveOpticsSim._accumulate_intensity!(SCALAR_CPU_STYLE, scalar_intensity, scalar_ef.state.field)
        AdaptiveOpticsSim._accumulate_intensity!(KA_CPU_STYLE, ka_intensity, ka_ef.state.field)
        mark_ka_cpu_kernel!(:accumulate_abs2_kernel!)
        @test ka_cpu_close(ka_intensity, scalar_intensity)

        scalar_centered = copy(scalar_field)
        ka_centered = copy(scalar_field)
        AdaptiveOpticsSim.apply_centering_phase!(SCALAR_CPU_STYLE, scalar_centered, 0.2)
        AdaptiveOpticsSim.apply_centering_phase!(KA_CPU_STYLE, ka_centered, 0.2)
        mark_ka_cpu_kernel!(:apply_centering_phase_kernel!)
        @test ka_cpu_close(ka_centered, scalar_centered)

        input = reshape(ComplexF64.(collect(1:16), collect(17:32)), 4, 4)
        scalar_copy = similar(input)
        ka_copy = similar(input)
        AdaptiveOpticsSim._complex_scale_copy!(SCALAR_CPU_STYLE, scalar_copy, input, 0.25)
        AdaptiveOpticsSim._complex_scale_copy!(KA_CPU_STYLE, ka_copy, input, 0.25)
        mark_ka_cpu_kernel!(:complex_scale_copy_kernel!)
        @test ka_copy == scalar_copy
        weights = fill(2.0 + 0.5im, 4, 4)
        AdaptiveOpticsSim._complex_hadamard!(SCALAR_CPU_STYLE, scalar_copy, weights)
        AdaptiveOpticsSim._complex_hadamard!(KA_CPU_STYLE, ka_copy, weights)
        mark_ka_cpu_kernel!(:complex_hadamard_kernel!)
        @test ka_copy == scalar_copy

        freqs = collect(range(0.0, 1.0; length=4))
        scalar_transfer = similar(input)
        ka_transfer = similar(input)
        AdaptiveOpticsSim.build_fresnel_transfer!(SCALAR_CPU_STYLE, scalar_transfer, freqs, -0.1)
        AdaptiveOpticsSim.build_fresnel_transfer!(KA_CPU_STYLE, ka_transfer, freqs, -0.1)
        mark_ka_cpu_kernel!(:fresnel_transfer_kernel!)
        @test ka_cpu_close(ka_transfer, scalar_transfer)

        pupil = trues(4, 4)
        basis = reshape(collect(range(0.0, 1.0; length=32)), 4, 4, 2)
        coeffs = [0.25, -0.5]
        scalar_opd = Matrix{Float64}(undef, 4, 4)
        ka_opd = similar(scalar_opd)
        AdaptiveOpticsSim.combine_basis!(SCALAR_CPU_STYLE, scalar_opd, basis, coeffs, pupil)
        AdaptiveOpticsSim.combine_basis!(KA_CPU_STYLE, ka_opd, basis, coeffs, pupil)
        mark_ka_cpu_kernel!(:combine_basis_kernel!)
        @test ka_cpu_close(ka_opd, scalar_opd)

        dm_scalar = DeformableMirror(tel; n_act=3, influence_width=0.3)
        dm_ka = DeformableMirror(tel; n_act=3, influence_width=0.3)
        fill!(dm_scalar.state.modes, 0.0)
        fill!(dm_ka.state.modes, 0.0)
        AdaptiveOpticsSim.build_influence_functions!(
            SCALAR_CPU_STYLE, dm_scalar, tel, dm_scalar.params.influence_model, AdaptiveOpticsSim.topology(dm_scalar))
        AdaptiveOpticsSim.build_influence_functions!(
            KA_CPU_STYLE, dm_ka, tel, dm_ka.params.influence_model, AdaptiveOpticsSim.topology(dm_ka))
        mark_ka_cpu_kernel!(:dm_mode_kernel!)
        @test ka_cpu_close(dm_ka.state.modes, dm_scalar.state.modes)

        dm_scalar.state.coefs .= collect(range(-0.1, 0.1; length=length(dm_scalar.state.coefs)))
        dm_ka.state.coefs .= dm_scalar.state.coefs
        AdaptiveOpticsSim.prepare_actuator_commands!(dm_scalar)
        AdaptiveOpticsSim.prepare_actuator_commands!(dm_ka)
        AdaptiveOpticsSim._apply_opd_separable!(SCALAR_CPU_STYLE, dm_scalar, tel)
        AdaptiveOpticsSim._apply_opd_separable!(KA_CPU_STYLE, dm_ka, tel)
        mark_ka_cpu_kernel!(:dm_separable_tmp_kernel!, :dm_separable_finalize_kernel!)
        @test ka_cpu_close(dm_ka.state.opd, dm_scalar.state.opd)
    end

    @testset "Detector kernels" begin
        frame = reshape(collect(1.0:25.0), 5, 5)
        scratch = similar(frame)
        scalar_frame = copy(frame)
        ka_frame = copy(frame)
        scalar_scratch = similar(frame)
        ka_scratch = similar(frame)
        gaussian = GaussianPixelResponse(response_width_px=0.5)
        AdaptiveOpticsSim.apply_response!(SCALAR_CPU_STYLE, gaussian, scalar_frame, scalar_scratch)
        AdaptiveOpticsSim.apply_response!(KA_CPU_STYLE, gaussian, ka_frame, ka_scratch)
        mark_ka_cpu_kernel!(:separable_response_rows_kernel!, :separable_response_cols_kernel!)
        @test ka_cpu_close(ka_frame, scalar_frame)

        sampled = SampledFrameResponse([0.0 0.1 0.0; 0.1 0.6 0.1; 0.0 0.1 0.0])
        scalar_frame .= frame
        ka_frame .= frame
        AdaptiveOpticsSim.apply_response!(SCALAR_CPU_STYLE, sampled, scalar_frame, scalar_scratch)
        AdaptiveOpticsSim.apply_response!(KA_CPU_STYLE, sampled, ka_frame, ka_scratch)
        mark_ka_cpu_kernel!(:sampled_response_kernel!)
        @test ka_cpu_close(ka_frame, scalar_frame)

        cube = reshape(collect(1.0:50.0), 2, 5, 5)
        scalar_cube = copy(cube)
        ka_cube = copy(cube)
        scalar_cube_scratch = similar(cube)
        ka_cube_scratch = similar(cube)
        AdaptiveOpticsSim._batched_apply_response!(SCALAR_CPU_STYLE, gaussian, scalar_cube, scalar_cube_scratch)
        AdaptiveOpticsSim._batched_apply_response!(KA_CPU_STYLE, gaussian, ka_cube, ka_cube_scratch)
        mark_ka_cpu_kernel!(:separable_response_stack_rows_kernel!, :separable_response_stack_cols_kernel!)
        @test ka_cpu_close(ka_cube, scalar_cube)

        scalar_cube .= cube
        ka_cube .= cube
        AdaptiveOpticsSim._batched_apply_response!(SCALAR_CPU_STYLE, sampled, scalar_cube, scalar_cube_scratch)
        AdaptiveOpticsSim._batched_apply_response!(KA_CPU_STYLE, sampled, ka_cube, ka_cube_scratch)
        mark_ka_cpu_kernel!(:sampled_response_stack_kernel!)
        @test ka_cpu_close(ka_cube, scalar_cube)

        noise = reshape(collect(range(-0.2, 0.2; length=25)), 5, 5)
        scalar_frame .= frame
        ka_frame .= frame
        AdaptiveOpticsSim.apply_column_noise!(SCALAR_CPU_STYLE, scalar_frame, noise, 0.5)
        AdaptiveOpticsSim.apply_column_noise!(KA_CPU_STYLE, ka_frame, noise, 0.5)
        mark_ka_cpu_kernel!(:add_column_noise_kernel!)
        @test ka_cpu_close(ka_frame, scalar_frame)

        output_model = StaticCMOSOutputPattern(2, [1.0, 1.1, 1.2], [0.0, 0.5, 1.0])
        scalar_frame .= frame
        ka_frame .= frame
        AdaptiveOpticsSim.apply_output_model!(SCALAR_CPU_STYLE, output_model, scalar_frame)
        AdaptiveOpticsSim.apply_output_model!(KA_CPU_STYLE, output_model, ka_frame)
        mark_ka_cpu_kernel!(:apply_cmos_output_pattern_kernel!)
        @test ka_cpu_close(ka_frame, scalar_frame)

        for model in (
            ReferencePixelCommonModeCorrection(1, 1),
            ReferenceRowCommonModeCorrection(1),
            ReferenceColumnCommonModeCorrection(1),
            ReferenceOutputCommonModeCorrection(3; edge_rows=1, edge_cols=1),
        )
            scalar_cube .= cube
            ka_cube .= cube
            fill!(scalar_cube_scratch, 0.0)
            fill!(ka_cube_scratch, 0.0)
            AdaptiveOpticsSim._batched_apply_readout_correction!(SCALAR_CPU_STYLE, model, scalar_cube, scalar_cube_scratch)
            AdaptiveOpticsSim._batched_apply_readout_correction!(KA_CPU_STYLE, model, ka_cube, ka_cube_scratch)
            @test ka_cpu_close(ka_cube, scalar_cube)
        end
        mark_ka_cpu_kernel!(
            :reference_pixel_bias_kernel!,
            :reference_row_bias_kernel!,
            :subtract_row_bias_kernel!,
            :reference_column_bias_kernel!,
            :subtract_column_bias_kernel!,
            :reference_output_bias_kernel!,
            :subtract_output_bias_kernel!,
        )
    end

    @testset "Pyramid kernels" begin
        tel = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
        wfs = PyramidWFS(tel; pupil_samples=4, modulation=2.0, modulation_points=3, mode=Diffractive())

        scalar_phasor = similar(wfs.state.phasor)
        ka_phasor = similar(wfs.state.phasor)
        AdaptiveOpticsSim._build_pyramid_phasor!(SCALAR_CPU_STYLE, scalar_phasor)
        AdaptiveOpticsSim._build_pyramid_phasor!(KA_CPU_STYLE, ka_phasor)
        mark_ka_cpu_kernel!(:pyramid_phasor_kernel!)
        @test ka_cpu_close(ka_phasor, scalar_phasor)

        scalar_mask = similar(wfs.state.pyramid_mask)
        ka_mask = similar(wfs.state.pyramid_mask)
        AdaptiveOpticsSim._build_pyramid_mask!(SCALAR_CPU_STYLE, scalar_mask, wfs, tel)
        AdaptiveOpticsSim._build_pyramid_mask!(KA_CPU_STYLE, ka_mask, wfs, tel)
        mark_ka_cpu_kernel!(:pyramid_mask_kernel!)
        @test ka_cpu_close(ka_mask, scalar_mask)

        scalar_phases = similar(wfs.state.modulation_phases)
        ka_phases = similar(wfs.state.modulation_phases)
        AdaptiveOpticsSim._build_modulation_phases!(SCALAR_CPU_STYLE, scalar_phases, wfs.params.modulation,
            (tel.params.resolution + 1) / 2, wfs.params.modulation_points, tel.params.resolution)
        AdaptiveOpticsSim._build_modulation_phases!(KA_CPU_STYLE, ka_phases, wfs.params.modulation,
            (tel.params.resolution + 1) / 2, wfs.params.modulation_points, tel.params.resolution)
        mark_ka_cpu_kernel!(:modulation_phases_kernel!)
        @test ka_cpu_close(ka_phases, scalar_phases)

        intensity = reshape(collect(1.0:(4 * 4 * 4 * 4)), 16, 16)
        scalar_slopes = similar(wfs.state.slopes)
        ka_slopes = similar(wfs.state.slopes)
        valid_mask = trues(4, 4)
        AdaptiveOpticsSim._pyramid_slopes!(SCALAR_CPU_STYLE, scalar_slopes, intensity, valid_mask, 2, 4, 16, 16,
            0, 0, 0, 8, 8, 0, 8, 8, (0, 0, 0, 0), (0, 0, 0, 0))
        AdaptiveOpticsSim._pyramid_slopes!(KA_CPU_STYLE, ka_slopes, intensity, valid_mask, 2, 4, 16, 16,
            0, 0, 0, 8, 8, 0, 8, 8, (0, 0, 0, 0), (0, 0, 0, 0))
        mark_ka_cpu_kernel!(:pyramid_slopes_kernel!)
        @test ka_slopes == scalar_slopes
    end

    @testset "BioEdge kernels" begin
        tel = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
        wfs = BioEdgeWFS(tel; pupil_samples=4, mode=Diffractive())

        scalar_edge_mask = similar(wfs.state.edge_mask)
        ka_edge_mask = similar(wfs.state.edge_mask)
        AdaptiveOpticsSim._update_edge_mask!(SCALAR_CPU_STYLE, scalar_edge_mask, tel.state.pupil, tel.params.resolution)
        AdaptiveOpticsSim._update_edge_mask!(KA_CPU_STYLE, ka_edge_mask, tel.state.pupil, tel.params.resolution)
        mark_ka_cpu_kernel!(:edge_mask_kernel!)
        @test ka_edge_mask == scalar_edge_mask

        scalar_phasor = similar(wfs.state.phasor)
        ka_phasor = similar(wfs.state.phasor)
        AdaptiveOpticsSim._build_bioedge_phasor!(SCALAR_CPU_STYLE, scalar_phasor)
        AdaptiveOpticsSim._build_bioedge_phasor!(KA_CPU_STYLE, ka_phasor)
        mark_ka_cpu_kernel!(:bioedge_phasor_kernel!)
        @test ka_cpu_close(ka_phasor, scalar_phasor)

        scalar_masks = similar(wfs.state.bioedge_masks)
        ka_masks = similar(wfs.state.bioedge_masks)
        AdaptiveOpticsSim._build_bioedge_masks!(SCALAR_CPU_STYLE, scalar_masks, Float64)
        AdaptiveOpticsSim._build_bioedge_masks!(KA_CPU_STYLE, ka_masks, Float64)
        mark_ka_cpu_kernel!(:bioedge_masks_kernel!)
        @test ka_cpu_close(ka_masks, scalar_masks)

        mask = falses(8, 8)
        mask[1:2:end, :] .= true
        scalar_binned = Matrix{Bool}(undef, 4, 4)
        ka_binned = similar(scalar_binned)
        AdaptiveOpticsSim._bin_edge_mask!(SCALAR_CPU_STYLE, scalar_binned, mask, 2, 4, 4)
        AdaptiveOpticsSim._bin_edge_mask!(KA_CPU_STYLE, ka_binned, mask, 2, 4, 4)
        mark_ka_cpu_kernel!(:bin_edge_mask_kernel!)
        @test ka_binned == scalar_binned
    end

    @testset "Elongation kernel" begin
        intensity = reshape(collect(range(0.0, 1.0; length=64)), 8, 8)
        scalar_intensity = copy(intensity)
        ka_intensity = copy(intensity)
        scalar_tmp = similar(intensity)
        ka_tmp = similar(intensity)
        scalar_kernel = Vector{Float64}(undef, 1)
        ka_kernel = similar(scalar_kernel)
        AdaptiveOpticsSim._apply_elongation!(SCALAR_CPU_STYLE, scalar_intensity, scalar_tmp,
            AdaptiveOpticsSim.apply_elongation!(scalar_intensity, 1.6, scalar_tmp, scalar_kernel), 1, 8, 8)
        ka_kernel = AdaptiveOpticsSim.apply_elongation!(ka_intensity, 1.6, ka_tmp, ka_kernel)
        scalar_baseline = copy(scalar_intensity)
        ka_probe = copy(intensity)
        ka_tmp2 = similar(ka_probe)
        AdaptiveOpticsSim._apply_elongation!(KA_CPU_STYLE, ka_probe, ka_tmp2, ka_kernel, 1, 8, 8)
        copyto!(ka_probe, ka_tmp2)
        mark_ka_cpu_kernel!(:elongation_apply_kernel!)
        @test ka_cpu_close(ka_probe, scalar_baseline)

        intensity_stack = reshape(collect(range(0.0, 1.0; length=128)), 2, 8, 8)
        scalar_stack = copy(intensity_stack)
        ka_stack = copy(intensity_stack)
        scalar_stack_tmp = similar(intensity_stack)
        ka_stack_tmp = similar(intensity_stack)
        AdaptiveOpticsSim._apply_elongation_stack!(SCALAR_CPU_STYLE, scalar_stack, scalar_stack_tmp, ka_kernel, 1, 2, 8, 8)
        AdaptiveOpticsSim._apply_elongation_stack!(KA_CPU_STYLE, ka_stack, ka_stack_tmp, ka_kernel, 1, 2, 8, 8)
        copyto!(scalar_stack, scalar_stack_tmp)
        copyto!(ka_stack, ka_stack_tmp)
        mark_ka_cpu_kernel!(:elongation_apply_stack_kernel!)
        @test ka_cpu_close(ka_stack, scalar_stack)
    end

    @testset "KA CPU kernel inventory" begin
        # Deferred kernels are still accounted for: most require larger
        # end-to-end wfs/tomography fixtures, while masked_sum2d_kernel! uses
        # atomics that are meaningful for GPU backends but not portable to the
        # KA CPU Array backend on Julia 1.12.
        deferred = Set([
            :accumulate_selected_block_kernel!,
            :accumulate_selected_block_transpose_kernel!,
            :apply_command_kernel!,
            :calibration_ramp_kernel!,
            :complex_to_real_scaled_stack_kernel!,
            :covariance_matrix_kernel!,
            :curvature_abs2_stack_kernel!,
            :curvature_branch_field_from_input_kernel!,
            :curvature_branch_field_stack_kernel!,
            :curvature_channel_pack_kernel!,
            :curvature_frame_pack_kernel!,
            :curvature_phasor_kernel!,
            :curvature_sample_branch_kernel!,
            :curvature_signal_from_channels_kernel!,
            :curvature_signal_from_frame_kernel!,
            :diagonal_matrix_kernel!,
            :fit_source_average_kernel!,
            :gather_bioedge_slopes_kernel!,
            :gather_pyramid_slopes_kernel!,
            :gather_stencil_data_kernel!,
            :gather_zernike_signal_kernel!,
            :guide_grid_kernel!,
            :guide_grid_stack_kernel!,
            :inject_column_negative_kernel!,
            :inject_column_positive_kernel!,
            :inject_row_negative_kernel!,
            :inject_row_positive_kernel!,
            :lift_gather_kernel!,
            :lift_scatter_update_kernel!,
            :masked_sum2d_kernel!,
            :multiply_kernel_fft_stack_kernel!,
            :real_to_complex_stack_kernel!,
            :scaled_shifted_coord_stack_kernel!,
            :selected_covariance_block_kernel!,
            :submatrix_extract_kernel!,
            :zernike_phasor_kernel!,
            :zernike_signal_kernel!,
        ])
        all_kernels = source_kernel_names()
        classified = union(KA_CPU_EXERCISED_KERNELS, deferred)
        @test isempty(setdiff(KA_CPU_EXERCISED_KERNELS, all_kernels))
        @test isempty(setdiff(deferred, all_kernels))
        @test isempty(setdiff(all_kernels, classified))
        @test length(KA_CPU_EXERCISED_KERNELS) >= 40
    end

    @testset "Python reference agreement" begin
        bundle = has_reference_bundle() ? load_reference_bundle() : nothing
        if bundle === nothing
            @test_skip "No reference bundle available"
        else
            for case in bundle.cases
                if case.id in ("shack_hartmann_geometric_ramp_xy", "shack_hartmann_geometric_ramp_y")
                    expected = load_reference_array(case)
                    scalar_actual = adapt_reference_actual(case, compute_reference_actual(case))
                    ka_actual = adapt_reference_actual(case, compute_reference_actual_ka_cpu(case))
                    @test isapprox(scalar_actual, expected; atol=case.atol, rtol=case.rtol)
                    @test isapprox(ka_actual, expected; atol=case.atol, rtol=case.rtol)
                    @test ka_cpu_close(ka_actual, scalar_actual; atol=case.atol, rtol=case.rtol)
                end
            end
        end
    end
end
