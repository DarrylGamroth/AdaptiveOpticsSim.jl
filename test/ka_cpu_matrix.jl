using KernelAbstractions

const KA_CPU_STYLE = AdaptiveOpticsSim.AcceleratorStyle(KernelAbstractions.CPU())
const SCALAR_CPU_STYLE = AdaptiveOpticsSim.ScalarCPUStyle()

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
        @test ka_fftshift == scalar_fftshift

        scalar_shift = similar(src)
        ka_shift = similar(src)
        AdaptiveOpticsSim._circshift2d!(SCALAR_CPU_STYLE, scalar_shift, src, (1, -1))
        AdaptiveOpticsSim._circshift2d!(KA_CPU_STYLE, ka_shift, src, (1, -1))
        @test ka_shift == scalar_shift

        scalar_bin = Matrix{Float64}(undef, 2, 2)
        ka_bin = similar(scalar_bin)
        AdaptiveOpticsSim._bin2d!(SCALAR_CPU_STYLE, scalar_bin, src, 2)
        AdaptiveOpticsSim._bin2d!(KA_CPU_STYLE, ka_bin, src, 2)
        @test ka_bin == scalar_bin

        scalar_freqs = Vector{Float64}(undef, 8)
        ka_freqs = similar(scalar_freqs)
        AdaptiveOpticsSim._fftfreq!(SCALAR_CPU_STYLE, scalar_freqs, 8, 0.25, 0.0)
        AdaptiveOpticsSim._fftfreq!(KA_CPU_STYLE, ka_freqs, 8, 0.25, 0.0)
        @test ka_freqs == scalar_freqs
    end

    @testset "Telescope and pupil kernels" begin
        params = TelescopeParams{Float64}(16, 8.0, 1e-3, 0.15, 0.0)
        scalar_pupil = Matrix{Bool}(undef, 16, 16)
        ka_pupil = similar(scalar_pupil)
        AdaptiveOpticsSim._generate_pupil!(SCALAR_CPU_STYLE, scalar_pupil, params)
        AdaptiveOpticsSim._generate_pupil!(KA_CPU_STYLE, ka_pupil, params)
        @test ka_pupil == scalar_pupil

        scalar_spiders = copy(scalar_pupil)
        ka_spiders = copy(ka_pupil)
        angles = [0.0, 45.0, 90.0]
        AdaptiveOpticsSim._apply_spiders!(SCALAR_CPU_STYLE, scalar_spiders, angles, 0.1, 0.0, 0.0, 8.5, 8.5, 8.0, 16)
        AdaptiveOpticsSim._apply_spiders!(KA_CPU_STYLE, ka_spiders, angles, 0.1, 0.0, 0.0, 8.5, 8.5, 8.0, 16)
        @test ka_spiders == scalar_spiders
    end

    @testset "Wavefront support kernels" begin
        pupil = trues(8, 8)
        pupil[1:2, :] .= false
        threshold = 0.4
        scalar_valid = Matrix{Bool}(undef, 2, 2)
        ka_valid = similar(scalar_valid)
        AdaptiveOpticsSim._set_valid_subapertures!(SCALAR_CPU_STYLE, scalar_valid, pupil, threshold, 4, 2)
        AdaptiveOpticsSim._set_valid_subapertures!(KA_CPU_STYLE, ka_valid, pupil, threshold, 4, 2)
        @test ka_valid == scalar_valid

        opd = reshape(collect(1.0:64.0), 8, 8)
        scalar_slopes = Vector{Float64}(undef, 8)
        ka_slopes = similar(scalar_slopes)
        AdaptiveOpticsSim._geometric_slopes!(SCALAR_CPU_STYLE, scalar_slopes, opd, scalar_valid, 4, 2, 4)
        AdaptiveOpticsSim._geometric_slopes!(KA_CPU_STYLE, ka_slopes, opd, ka_valid, 4, 2, 4)
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
        @test ka_edge_slopes == scalar_edge_slopes
    end

    @testset "Pyramid kernels" begin
        tel = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
        wfs = PyramidWFS(tel; n_subap=4, modulation=2.0, modulation_points=3, mode=Diffractive())

        scalar_phasor = similar(wfs.state.phasor)
        ka_phasor = similar(wfs.state.phasor)
        AdaptiveOpticsSim._build_pyramid_phasor!(SCALAR_CPU_STYLE, scalar_phasor)
        AdaptiveOpticsSim._build_pyramid_phasor!(KA_CPU_STYLE, ka_phasor)
        @test ka_cpu_close(ka_phasor, scalar_phasor)

        scalar_mask = similar(wfs.state.pyramid_mask)
        ka_mask = similar(wfs.state.pyramid_mask)
        AdaptiveOpticsSim._build_pyramid_mask!(SCALAR_CPU_STYLE, scalar_mask, wfs, tel)
        AdaptiveOpticsSim._build_pyramid_mask!(KA_CPU_STYLE, ka_mask, wfs, tel)
        @test ka_cpu_close(ka_mask, scalar_mask)

        scalar_phases = similar(wfs.state.modulation_phases)
        ka_phases = similar(wfs.state.modulation_phases)
        AdaptiveOpticsSim._build_modulation_phases!(SCALAR_CPU_STYLE, scalar_phases, wfs.params.modulation,
            (tel.params.resolution + 1) / 2, wfs.params.modulation_points, tel.params.resolution)
        AdaptiveOpticsSim._build_modulation_phases!(KA_CPU_STYLE, ka_phases, wfs.params.modulation,
            (tel.params.resolution + 1) / 2, wfs.params.modulation_points, tel.params.resolution)
        @test ka_cpu_close(ka_phases, scalar_phases)

        intensity = reshape(collect(1.0:(4 * 4 * 4 * 4)), 16, 16)
        scalar_slopes = similar(wfs.state.slopes)
        ka_slopes = similar(wfs.state.slopes)
        valid_mask = trues(4, 4)
        AdaptiveOpticsSim._pyramid_slopes!(SCALAR_CPU_STYLE, scalar_slopes, intensity, valid_mask, 2, 4, 16, 16,
            0, 0, 0, 8, 8, 0, 8, 8, (0, 0, 0, 0), (0, 0, 0, 0))
        AdaptiveOpticsSim._pyramid_slopes!(KA_CPU_STYLE, ka_slopes, intensity, valid_mask, 2, 4, 16, 16,
            0, 0, 0, 8, 8, 0, 8, 8, (0, 0, 0, 0), (0, 0, 0, 0))
        @test ka_slopes == scalar_slopes
    end

    @testset "BioEdge kernels" begin
        tel = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
        wfs = BioEdgeWFS(tel; n_subap=4, mode=Diffractive())

        scalar_edge_mask = similar(wfs.state.edge_mask)
        ka_edge_mask = similar(wfs.state.edge_mask)
        AdaptiveOpticsSim._update_edge_mask!(SCALAR_CPU_STYLE, scalar_edge_mask, tel.state.pupil, tel.params.resolution)
        AdaptiveOpticsSim._update_edge_mask!(KA_CPU_STYLE, ka_edge_mask, tel.state.pupil, tel.params.resolution)
        @test ka_edge_mask == scalar_edge_mask

        scalar_phasor = similar(wfs.state.phasor)
        ka_phasor = similar(wfs.state.phasor)
        AdaptiveOpticsSim._build_bioedge_phasor!(SCALAR_CPU_STYLE, scalar_phasor)
        AdaptiveOpticsSim._build_bioedge_phasor!(KA_CPU_STYLE, ka_phasor)
        @test ka_cpu_close(ka_phasor, scalar_phasor)

        scalar_masks = similar(wfs.state.bioedge_masks)
        ka_masks = similar(wfs.state.bioedge_masks)
        AdaptiveOpticsSim._build_bioedge_masks!(SCALAR_CPU_STYLE, scalar_masks, Float64)
        AdaptiveOpticsSim._build_bioedge_masks!(KA_CPU_STYLE, ka_masks, Float64)
        @test ka_cpu_close(ka_masks, scalar_masks)

        mask = falses(8, 8)
        mask[1:2:end, :] .= true
        scalar_binned = Matrix{Bool}(undef, 4, 4)
        ka_binned = similar(scalar_binned)
        AdaptiveOpticsSim._bin_edge_mask!(SCALAR_CPU_STYLE, scalar_binned, mask, 2, 4, 4)
        AdaptiveOpticsSim._bin_edge_mask!(KA_CPU_STYLE, ka_binned, mask, 2, 4, 4)
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
        @test ka_cpu_close(ka_probe, scalar_baseline)
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
