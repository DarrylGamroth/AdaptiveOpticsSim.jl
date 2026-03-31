@testset "Atmosphere propagation" begin
    tel = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    atm = KolmogorovAtmosphere(tel; r0=0.2, L0=25.0)
    advance!(atm, tel; rng=MersenneTwister(1))
    propagate!(atm, tel)
    @test size(atm.state.opd) == (32, 32)
    @test sum(abs.(tel.state.opd)) > 0

    delta = tel.params.diameter / tel.params.resolution
    ensure_psd!(atm, delta)
    psd_snapshot = copy(atm.state.psd)
    ensure_psd!(atm, delta)
    @test psd_snapshot == atm.state.psd

    atm = KolmogorovAtmosphere(tel; r0=0.1, L0=25.0)
    ensure_psd!(atm, delta)
    @test psd_snapshot != atm.state.psd
end

@testset "Infinite atmosphere boundary math" begin
    tel = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    @test AdaptiveOpticsSim.default_infinite_screen_resolution(tel.params.resolution) == 3 * tel.params.resolution
    @test AdaptiveOpticsSim.default_infinite_stencil_size(tel.params.resolution) >= 257

    stencil = AdaptiveOpticsSim.infinite_boundary_stencil(4, 0.25;
        stencil_size=9,
        orientation=:column,
        side=:positive,
        tail_stride=4,
        T=Float64,
    )
    @test size(stencil.boundary_coords) == (9, 2)
    @test size(stencil.stencil_positions, 2) == 2
    @test all(stencil.boundary_coords[:, 1] .== 0)

    row_stencil = AdaptiveOpticsSim.infinite_boundary_stencil(4, 0.25;
        stencil_size=9,
        orientation=:row,
        side=:negative,
        tail_stride=4,
        T=Float64,
    )
    @test all(row_stencil.boundary_coords[:, 2] .== 10)

    op = AdaptiveOpticsSim.boundary_injection_operator(stencil, 0.2, 25.0)
    @test size(op.predictor, 1) == size(stencil.boundary_coords, 1)
    @test size(op.predictor, 2) == size(stencil.stencil_coords, 1)
    @test size(op.residual_factor, 1) == size(stencil.boundary_coords, 1)
    @test isfinite(op.condition_ratio)
    @test op.condition_ratio > 0
    @test isapprox(op.residual_covariance, op.residual_factor * transpose(op.residual_factor); rtol=1e-10, atol=1e-10)

    rng = MersenneTwister(12)
    stencil_data = randn(rng, size(op.predictor, 2))
    nsamp = 4000
    samples = Matrix{Float64}(undef, size(op.predictor, 1), nsamp)
    expected_mean = op.predictor * stencil_data
    for i in 1:nsamp
        samples[:, i] = AdaptiveOpticsSim.sample_boundary_line(op, stencil_data, rng)
    end
    empirical_mean = vec(mean(samples; dims=2))
    centered = samples .- empirical_mean
    empirical_cov = centered * transpose(centered) / (nsamp - 1)
    @test isapprox(empirical_mean, expected_mean; rtol=0.08, atol=0.08)
    @test isapprox(empirical_cov, op.residual_covariance; rtol=0.12, atol=0.12)

    screen = InfinitePhaseScreen(tel; r0=0.2, L0=25.0, screen_resolution=33, stencil_size=35)
    @test size(screen.state.screen) == (35, 35)
    @test size(screen.state.extract_buffer) == (32, 32)
    @test size(screen.state.column_positive.operator.predictor, 1) == 35
    @test size(screen.state.row_negative.operator.predictor, 1) == 35

    atm = InfiniteMultiLayerAtmosphere(tel;
        r0=0.2,
        L0=25.0,
        fractional_cn2=[1.0],
        wind_speed=[0.0],
        wind_direction=[0.0],
        altitude=[0.0],
        screen_resolution=33,
        stencil_size=35,
    )
    @test length(atm.layers) == 1
    @test size(atm.state.opd) == (32, 32)
    advance!(atm, tel; rng=MersenneTwister(1))
    propagate!(atm, tel)
    @test size(atm.layers[1].screen.state.screen) == (35, 35)
    @test atm.layers[1].screen.state.initialized
    @test atm.layers[1].state.integer_shift_x == 0
    @test atm.layers[1].state.integer_shift_y == 0
    @test tel.state.opd == atm.state.opd .* tel.state.pupil

    @test_throws InvalidConfiguration AdaptiveOpticsSim.infinite_boundary_stencil(4, 0.25;
        stencil_size=8,
        orientation=:column,
        side=:positive,
    )
    @test_throws NumericalConditionError AdaptiveOpticsSim.boundary_injection_operator(stencil, 0.2, 25.0; conditioning_tol=1.0)
end

@testset "Infinite atmosphere stepping regressions" begin
    tel = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)

    normalized_correlation(a::AbstractMatrix, b::AbstractMatrix) =
        dot(vec(a), vec(b)) / sqrt(dot(vec(a), vec(a)) * dot(vec(b), vec(b)))

    stationary = InfiniteMultiLayerAtmosphere(tel;
        r0=0.2,
        L0=25.0,
        fractional_cn2=[1.0],
        wind_speed=[0.0],
        wind_direction=[0.0],
        altitude=[0.0],
        screen_resolution=33,
        stencil_size=35,
    )
    advance!(stationary, tel; rng=MersenneTwister(4))
    stationary_snapshot = copy(stationary.state.opd)
    advance!(stationary, tel; rng=MersenneTwister(5))
    @test stationary.state.opd == stationary_snapshot
    @test stationary.layers[1].state.integer_shift_x == 0
    @test stationary.layers[1].state.integer_shift_y == 0

    small_step = InfiniteMultiLayerAtmosphere(tel;
        r0=0.2,
        L0=25.0,
        fractional_cn2=[1.0],
        wind_speed=[50.0],
        wind_direction=[0.0],
        altitude=[0.0],
        screen_resolution=33,
        stencil_size=35,
    )
    rng = MersenneTwister(6)
    advance!(small_step, tel; rng=rng)
    opd_1 = copy(small_step.state.opd)
    advance!(small_step, tel; rng=rng)
    opd_2 = copy(small_step.state.opd)
    corr = normalized_correlation(opd_1, opd_2)
    @test corr > 0.95
    @test small_step.layers[1].state.integer_shift_x == 0

    row_step = InfiniteMultiLayerAtmosphere(tel;
        r0=0.2,
        L0=25.0,
        fractional_cn2=[1.0],
        wind_speed=[500.0],
        wind_direction=[90.0],
        altitude=[0.0],
        screen_resolution=33,
        stencil_size=35,
    )
    advance!(row_step, tel; rng=MersenneTwister(7))
    @test row_step.layers[1].state.integer_shift_y == 1
    @test row_step.layers[1].state.integer_shift_x == 0

    diagonal_speed = 250.0 * sqrt(2)
    diagonal_step = InfiniteMultiLayerAtmosphere(tel;
        r0=0.2,
        L0=25.0,
        fractional_cn2=[1.0],
        wind_speed=[diagonal_speed],
        wind_direction=[45.0],
        altitude=[0.0],
        screen_resolution=33,
        stencil_size=35,
    )
    diagonal_rng = MersenneTwister(71)
    advance!(diagonal_step, tel; rng=diagonal_rng)
    advance!(diagonal_step, tel; rng=diagonal_rng)
    @test diagonal_step.layers[1].state.integer_shift_x == 1
    @test diagonal_step.layers[1].state.integer_shift_y == 1
    @test abs(diagonal_step.layers[1].state.offset_x) < 1e-10
    @test abs(diagonal_step.layers[1].state.offset_y) < 1e-10

    function infinite_trace(; seed::Integer=73, steps::Int=10)
        tel_local = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
        src_local = Source(band=:I, magnitude=0.0, coordinates=(15.0, 30.0))
        atm = InfiniteMultiLayerAtmosphere(tel_local;
            r0=0.2,
            L0=25.0,
            fractional_cn2=[0.7, 0.3],
            wind_speed=[6.0, 3.0],
            wind_direction=[0.0, 120.0],
            altitude=[0.0, 5000.0],
            screen_resolution=33,
            stencil_size=35,
        )
        rng = MersenneTwister(seed)
        opd_trace = Matrix{Float64}[]
        for _ in 1:steps
            advance!(atm, tel_local; rng=rng)
            propagate!(atm, tel_local, src_local)
            push!(opd_trace, copy(Array(tel_local.state.opd)))
        end
        return (; opd_trace, screen=copy(Array(atm.layers[1].screen.state.screen)))
    end

    trace_a = infinite_trace()
    trace_b = infinite_trace()
    @test trace_a.opd_trace == trace_b.opd_trace
    @test trace_a.screen == trace_b.screen

    function trajectory_std_windows(; seed::Integer=79, steps::Int=24)
        tel_local = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
        atm = InfiniteMultiLayerAtmosphere(tel_local;
            r0=0.2,
            L0=25.0,
            fractional_cn2=[0.5, 0.5],
            wind_speed=[10.0, 6.0],
            wind_direction=[0.0, 90.0],
            altitude=[0.0, 5000.0],
            screen_resolution=33,
            stencil_size=35,
        )
        rng = MersenneTwister(seed)
        stds = Float64[]
        for _ in 1:steps
            advance!(atm, tel_local; rng=rng)
            push!(stds, std(vec(Array(atm.state.opd))))
        end
        midpoint = steps ÷ 2
        return mean(@view(stds[1:midpoint])), mean(@view(stds[midpoint + 1:end]))
    end

    early_std, late_std = trajectory_std_windows()
    @test isapprox(late_std, early_std; rtol=0.2)

    wind_speed_px = tel.params.diameter / tel.params.resolution / tel.params.sampling_time
    finite = MultiLayerAtmosphere(tel;
        r0=0.2,
        L0=25.0,
        fractional_cn2=[1.0],
        wind_speed=[wind_speed_px],
        wind_direction=[0.0],
        altitude=[0.0],
    )
    infinite = InfiniteMultiLayerAtmosphere(tel;
        r0=0.2,
        L0=25.0,
        fractional_cn2=[1.0],
        wind_speed=[wind_speed_px],
        wind_direction=[0.0],
        altitude=[0.0],
        screen_resolution=33,
        stencil_size=35,
    )
    finite_rng = MersenneTwister(8)
    infinite_rng = MersenneTwister(8)
    advance!(finite, tel; rng=finite_rng)
    advance!(infinite, tel; rng=infinite_rng)
    finite_snapshot = copy(finite.state.opd)
    infinite_snapshot = copy(infinite.state.opd)
    period = AdaptiveOpticsSim.moving_layer_screen_resolution(tel.params.resolution)
    for _ in 1:period
        advance!(finite, tel; rng=finite_rng)
        advance!(infinite, tel; rng=infinite_rng)
    end
    @test finite.state.opd == finite_snapshot
    @test !isapprox(infinite.state.opd, infinite_snapshot; rtol=1e-8, atol=1e-8)
    for _ in 1:(2 * period)
        advance!(finite, tel; rng=finite_rng)
        advance!(infinite, tel; rng=infinite_rng)
    end
    @test finite.state.opd == finite_snapshot
    @test normalized_correlation(infinite.state.opd, infinite_snapshot) < 0.9
end

@testset "Source-aware atmosphere extraction" begin
    tel = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    onaxis = Source(band=:I, magnitude=0.0, coordinates=(0.0, 0.0))
    offaxis = Source(band=:I, magnitude=0.0, coordinates=(100.0, 0.0))
    lgs = LGSSource(magnitude=0.0, coordinates=(100.0, 0.0), altitude=10_000.0)

    shift_x_ngs, shift_y_ngs, footprint_ngs = AdaptiveOpticsSim.layer_source_geometry(offaxis, 5000.0, tel, Float64)
    shift_x_lgs, shift_y_lgs, footprint_lgs = AdaptiveOpticsSim.layer_source_geometry(lgs, 5000.0, tel, Float64)
    @test shift_x_ngs ≈ shift_x_lgs
    @test shift_y_ngs ≈ shift_y_lgs
    @test shift_y_ngs ≈ 0.0
    @test footprint_ngs == 1.0
    @test footprint_lgs ≈ 0.5

    function fill_x_ramp!(screen::AbstractMatrix)
        @inbounds for j in axes(screen, 2), i in axes(screen, 1)
            screen[i, j] = j
        end
        return screen
    end

    function offaxis_shift_regression(constructor; seed::Integer=41, kwargs...)
        atm = constructor(tel;
            r0=0.2,
            L0=25.0,
            fractional_cn2=[1.0],
            wind_speed=[0.0],
            wind_direction=[0.0],
            altitude=[5000.0],
            kwargs...,
        )
        advance!(atm, tel; rng=MersenneTwister(seed))
        if atm isa MultiLayerAtmosphere
            fill_x_ramp!(atm.layers[1].generator.state.opd)
            atm.layers[1].state.offset_x = 0.0
            atm.layers[1].state.offset_y = 0.0
        else
            fill_x_ramp!(atm.layers[1].screen.state.screen)
            atm.layers[1].state.offset_x = 0.0
            atm.layers[1].state.offset_y = 0.0
        end
        propagate!(atm, tel, onaxis)
        onaxis_opd = copy(tel.state.opd)
        propagate!(atm, tel, offaxis)
        offaxis_opd = copy(tel.state.opd)
        propagate!(atm, tel, lgs)
        lgs_opd = copy(tel.state.opd)
        return (; onaxis_opd, offaxis_opd, lgs_opd)
    end

    finite_offaxis = offaxis_shift_regression(MultiLayerAtmosphere)
    infinite_offaxis = offaxis_shift_regression(InfiniteMultiLayerAtmosphere; screen_resolution=33, stencil_size=35)
    @test mean(finite_offaxis.offaxis_opd) > mean(finite_offaxis.onaxis_opd)
    @test mean(infinite_offaxis.offaxis_opd) > mean(infinite_offaxis.onaxis_opd)
    @test std(vec(finite_offaxis.lgs_opd)) < std(vec(finite_offaxis.offaxis_opd))
    @test std(vec(infinite_offaxis.lgs_opd)) < std(vec(infinite_offaxis.offaxis_opd))

    function onaxis_fast_path_regression(constructor; kwargs...)
        atm = constructor(tel;
            r0=0.2,
            L0=25.0,
            fractional_cn2=[0.7, 0.3],
            wind_speed=[6.0, 3.0],
            wind_direction=[0.0, 120.0],
            altitude=[0.0, 5000.0],
            kwargs...,
        )
        rng = MersenneTwister(52)
        advance!(atm, tel; rng=rng)
        propagate!(atm, tel)
        base_opd = copy(tel.state.opd)
        @test !atm.state.source_geometry.valid
        propagate!(atm, tel, onaxis)
        @test tel.state.opd == base_opd
        @test !atm.state.source_geometry.valid
        propagate!(atm, tel, offaxis)
        @test atm.state.source_geometry.valid
        @test atm.state.source_geometry.cached_x_arcsec ≈ coordinates_xy_arcsec(offaxis)[1]
        @test atm.state.source_geometry.cached_y_arcsec ≈ coordinates_xy_arcsec(offaxis)[2]
        return atm
    end

    onaxis_fast_path_regression(MultiLayerAtmosphere)
    onaxis_fast_path_regression(InfiniteMultiLayerAtmosphere; screen_resolution=33, stencil_size=35)

    function ensemble_std(constructor, fractions; nsamp::Int=12, kwargs...)
        acc = 0.0
        for s in 1:nsamp
            atm = constructor(tel;
                r0=0.2,
                L0=25.0,
                fractional_cn2=fractions,
                wind_speed=fill(0.0, length(fractions)),
                wind_direction=fill(0.0, length(fractions)),
                altitude=fill(0.0, length(fractions)),
                kwargs...,
            )
            advance!(atm, tel; rng=MersenneTwister(s))
            acc += std(vec(atm.state.opd))
        end
        return acc / nsamp
    end

    finite_single_std = ensemble_std(MultiLayerAtmosphere, [1.0])
    infinite_single_std = ensemble_std(InfiniteMultiLayerAtmosphere, [1.0]; screen_resolution=33, stencil_size=35)
    infinite_equal_std = ensemble_std(InfiniteMultiLayerAtmosphere, [0.5, 0.5]; screen_resolution=33, stencil_size=35)
    infinite_uneven_std = ensemble_std(InfiniteMultiLayerAtmosphere, [0.8, 0.2]; screen_resolution=33, stencil_size=35)
    @test isapprox(infinite_single_std, finite_single_std; rtol=0.2)
    @test isapprox(infinite_equal_std, infinite_single_std; rtol=0.2)
    @test isapprox(infinite_uneven_std, infinite_single_std; rtol=0.2)
end

@testset "Sub-harmonic phase screens" begin
    tel = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    atm = KolmogorovAtmosphere(tel; r0=0.2, L0=25.0)
    delta = tel.params.diameter / tel.params.resolution
    rng = MersenneTwister(11)
    ws = PhaseStatsWorkspace(32; T=Float64)
    phs_base = ft_sh_phase_screen(atm, 32, delta; rng=rng, ws=ws, subharmonics=false)
    rng = MersenneTwister(11)
    phs_sh = ft_sh_phase_screen(atm, 32, delta; rng=rng, ws=ws, subharmonics=true,
        mode=FidelitySubharmonics())
    @test sum(abs.(phs_sh .- phs_base)) > 0
    @test AdaptiveOpticsSim.resolve_subharmonic_levels(25.0, 8.0) == 4
    @test AdaptiveOpticsSim.resolve_subharmonic_levels(200.0, 8.0) > 4

    atm_large = KolmogorovAtmosphere(tel; r0=0.2, L0=200.0)
    rng = MersenneTwister(11)
    phs_default = ft_sh_phase_screen(atm_large, 32, delta; rng=rng, ws=ws,
        subharmonics=true, mode=FidelitySubharmonics())
    rng = MersenneTwister(11)
    phs_legacy = ft_sh_phase_screen(atm_large, 32, delta; rng=rng, ws=ws,
        subharmonics=true, mode=FastSubharmonics())
    @test sum(abs.(phs_default .- phs_legacy)) > 0

    metrics_25_none = subharmonic_metrics(atm, tel.params.diameter; subharmonics=false)
    metrics_25_legacy = subharmonic_metrics(atm, tel.params.diameter;
        subharmonics=true, mode=FastSubharmonics())
    metrics_25_default = subharmonic_metrics(atm, tel.params.diameter;
        subharmonics=true, mode=FidelitySubharmonics())
    @test metrics_25_none.tiptilt < metrics_25_legacy.tiptilt < metrics_25_default.tiptilt
    @test metrics_25_none.structure < metrics_25_legacy.structure < metrics_25_default.structure

    metrics_200_none = subharmonic_metrics(atm_large, tel.params.diameter; subharmonics=false)
    metrics_200_legacy = subharmonic_metrics(atm_large, tel.params.diameter;
        subharmonics=true, mode=FastSubharmonics())
    metrics_200_default = subharmonic_metrics(atm_large, tel.params.diameter;
        subharmonics=true, mode=FidelitySubharmonics())
@test metrics_200_none.tiptilt < metrics_200_legacy.tiptilt < metrics_200_default.tiptilt
@test metrics_200_none.structure < metrics_200_legacy.structure < metrics_200_default.structure
end

@testset "Fidelity profiles" begin
    @test default_fidelity_profile() isa ScientificProfile
    @test default_subharmonic_mode(ScientificProfile()) isa FidelitySubharmonics
    @test default_subharmonic_mode(FastProfile()) isa FastSubharmonics
    @test default_ncpa_basis(ScientificProfile()).method isa KLHHtPSD
    @test default_ncpa_basis(FastProfile()).method isa KLDMModes

    mixed = ProfileBundle(ScientificProfile(); lift=FastProfile(), tomography=FastProfile())
    @test atmosphere_profile(mixed) isa ScientificProfile
    @test calibration_profile(mixed) isa ScientificProfile
    @test detector_profile(mixed) isa ScientificProfile
    @test lift_profile(mixed) isa FastProfile
    @test tomography_profile(mixed) isa FastProfile
    @test default_subharmonic_mode(mixed) isa FidelitySubharmonics
    @test default_ncpa_basis(mixed).method isa KLHHtPSD

    fast_cal = ProfileBundle(ScientificProfile(); calibration=FastProfile())
    @test default_ncpa_basis(fast_cal).method isa KLDMModes

    tel = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    atm = KolmogorovAtmosphere(tel; r0=0.15, L0=200.0)
    scientific = ft_sh_phase_screen(atm, 16, 0.1; rng=MersenneTwister(3), profile=ScientificProfile())
    fast = ft_sh_phase_screen(atm, 16, 0.1; rng=MersenneTwister(3), profile=FastProfile())
    @test !all(scientific .== fast)

    dm = DeformableMirror(tel; n_act=4, influence_width=0.3)
    coeffs = zeros(4)
    ncpa_fast = NCPA(tel, dm, atm; profile=FastProfile(), coefficients=coeffs)
    ncpa_scientific = NCPA(tel, dm, atm; profile=ScientificProfile(), coefficients=coeffs)
    @test size(ncpa_fast.opd) == size(tel.state.opd)
    @test size(ncpa_scientific.opd) == size(tel.state.opd)
end

@testset "Multi-layer atmosphere" begin
    tel = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    delta = tel.params.diameter / tel.params.resolution
    one_pixel_speed = delta / tel.params.sampling_time
    quarter_pixel_speed = 0.25 * one_pixel_speed

    atm = MultiLayerAtmosphere(tel;
        r0=0.2,
        L0=25.0,
        fractional_cn2=[1.0],
        wind_speed=[one_pixel_speed],
        wind_direction=[0.0],
        altitude=[0.0],
    )
    advance!(atm, tel; rng=MersenneTwister(3))
    first = copy(atm.state.opd)
    advance!(atm, tel; rng=MersenneTwister(3))
    second = copy(atm.state.opd)
    propagate!(atm, tel)
    @test size(atm.state.opd) == (16, 16)
    @test sum(abs.(tel.state.opd)) > 0
    @test second[:, 2:end] ≈ first[:, 1:end-1]

    stationary = MultiLayerAtmosphere(tel;
        r0=0.2,
        L0=25.0,
        fractional_cn2=[1.0],
        wind_speed=[0.0],
        wind_direction=[0.0],
        altitude=[0.0],
    )
    advance!(stationary, tel; rng=MersenneTwister(4))
    static_first = copy(stationary.state.opd)
    advance!(stationary, tel; rng=MersenneTwister(4))
    @test stationary.state.opd ≈ static_first

    subpixel = MultiLayerAtmosphere(tel;
        r0=0.2,
        L0=25.0,
        fractional_cn2=[1.0],
        wind_speed=[quarter_pixel_speed],
        wind_direction=[0.0],
        altitude=[0.0],
    )
    advance!(subpixel, tel; rng=MersenneTwister(5))
    subpixel_first = copy(subpixel.state.opd)
    for _ in 1:4
        advance!(subpixel, tel; rng=MersenneTwister(5))
    end
    @test subpixel.state.opd[:, 2:end] ≈ subpixel_first[:, 1:end-1]

    function ensemble_std(fractions; nsamp::Int=16)
        acc = 0.0
        for s in 1:nsamp
            atm_local = MultiLayerAtmosphere(tel;
                r0=0.2,
                L0=25.0,
                fractional_cn2=fractions,
                wind_speed=fill(0.0, length(fractions)),
                wind_direction=fill(0.0, length(fractions)),
                altitude=fill(0.0, length(fractions)),
            )
            advance!(atm_local, tel; rng=MersenneTwister(s))
            acc += std(vec(atm_local.state.opd))
        end
        return acc / nsamp
    end

    single_std = ensemble_std([1.0])
    equal_std = ensemble_std([0.5, 0.5])
    uneven_std = ensemble_std([0.8, 0.2])
    @test isapprox(equal_std, single_std; rtol=0.15)
    @test isapprox(uneven_std, single_std; rtol=0.15)

    replay_a = moving_atmosphere_trace(
        seed=11,
        steps=5,
        resolution=16,
        fractional_cn2=[0.7, 0.3],
        wind_speed=[one_pixel_speed, quarter_pixel_speed],
        wind_direction=[0.0, 90.0],
        altitude=[0.0, 5000.0],
    )
    replay_b = moving_atmosphere_trace(
        seed=11,
        steps=5,
        resolution=16,
        fractional_cn2=[0.7, 0.3],
        wind_speed=[one_pixel_speed, quarter_pixel_speed],
        wind_direction=[0.0, 90.0],
        altitude=[0.0, 5000.0],
    )
    replay_c = moving_atmosphere_trace(
        seed=12,
        steps=5,
        resolution=16,
        fractional_cn2=[0.7, 0.3],
        wind_speed=[one_pixel_speed, quarter_pixel_speed],
        wind_direction=[0.0, 90.0],
        altitude=[0.0, 5000.0],
    )
    @test replay_a == replay_b
    @test any(!=(0.0), abs.(replay_a[1] .- replay_c[1]))
end

@testset "Moving-atmosphere regressions" begin
    slopes_a = moving_wfs_slope_trace(seed=21, steps=4)
    slopes_b = moving_wfs_slope_trace(seed=21, steps=4)
    slopes_c = moving_wfs_slope_trace(seed=22, steps=4)
    @test slopes_a == slopes_b
    @test any(norm(slopes_a[i + 1] - slopes_a[i]) > 0 for i in 1:length(slopes_a)-1)
    @test any(norm(slopes_a[i] - slopes_c[i]) > 0 for i in eachindex(slopes_a))

    loop_a = moving_closed_loop_trace(seed=31, steps=5)
    loop_b = moving_closed_loop_trace(seed=31, steps=5)
    loop_c = moving_closed_loop_trace(seed=32, steps=5)
    @test loop_a == loop_b
    @test any(diff(loop_a.slope_norms) .!= 0)
    @test all(isfinite, loop_a.command_norms)
    @test maximum(loop_a.command_norms) > 0
    @test all(>(0), loop_a.wfs_energy)
    @test all(>(0), loop_a.science_energy)
    @test any(abs.(loop_a.slope_norms .- loop_c.slope_norms) .> 0)
end

