# Analytic coordinate oracles shared by scalar and accelerator checks.
# Native optical arrays are [x,y]; unequal slopes distinguish a shared transpose.
function atmosphere_axis_checks(selector, backend, style;
    T=Float64, kernel_style=style)
    A = AdaptiveOpticsSim.Atmospheres
    screen_host = T[2i + 7j for i in 1:17, j in 1:17]
    screen = backend(screen_host)
    out = backend(zeros(T, 5, 5))
    for scale in (one(T), T(0.5)), (dx, dy) in
            ((T(0.25), T(-0.5)), (T(-0.75), T(0.125)))
        sx = T(9) - 2scale - dx
        sy = T(9) - 2scale - dy
        expected = T[2(sx + scale * (i - 1)) +
            7(sy + scale * (j - 1)) for i in 1:5, j in 1:5]
        A._extract_shifted_screen!(style, out, screen,
            sx, sy, scale, one(T), 5, 17)
        @test Array(out) ≈ expected
        A._accumulate_shifted_screen!(style, out, screen,
            sx, sy, scale, T(2), 5, 17)
        @test Array(out) ≈ 3expected
        starts_x = backend(T[sx])
        starts_y = backend(T[sy])
        for accumulate in (false, true)
            Backends.launch_kernel!(kernel_style, A.moving_layer_replay_kernel!,
                out, screen, starts_x, starts_y, 1, scale, one(T),
                accumulate, 5, 17; ndrange=size(out))
            @test Array(out) ≈ (accumulate ? 2expected : expected)
        end
    end

    # Whole-pixel periodic wrapping, with different x/y markers.
    marker = zeros(T, 5, 5)
    marker[2, 4] = one(T)
    marker[4, 2] = T(3)
    for (dx, dy) in ((1, 0), (0, 1), (-1, 0), (0, -1), (7, -6))
        A._extract_shifted_screen!(style, out, backend(marker),
            A.normalize_start_coordinate(T(1 - dx), 5),
            A.normalize_start_coordinate(T(1 - dy), 5),
            one(T), one(T), 5, 5)
        @test Array(out) == circshift(marker, (dx, dy))
    end

    source = Source(band=:I, magnitude=0, separation_arcsec=5,
        position_angle_deg=atan(4, 3) * 180 / pi, T=T)
    sx, sy, scale = A.layer_source_geometry(source, T(1000),
        (T(2), T(5)), T)
    @test sx ≈ T(3pi / (180 * 3600) * 1000 / 2)
    @test sy ≈ T(4pi / (180 * 3600) * 1000 / 5)
    @test scale == one(T)

    tel = Telescope(resolution=5, diameter=T(5),
        central_obstruction=zero(T), T=T, backend=selector)
    # 1.25 pixels crosses the infinite-screen integer/fractional boundary.
    for constructor in (MultiLayerAtmosphere, InfiniteMultiLayerAtmosphere),
            (angle, dx, dy) in ((0, 1.25, 0), (90, 0, 1.25),
                (180, -1.25, 0), (270, 0, -1.25))
        kwargs = constructor === InfiniteMultiLayerAtmosphere ?
            (; screen_resolution=9, stencil_size=17) : (;)
        atm = constructor(tel; r0=T(0.2),
            reference_wavelength_m=T(500e-9), L0=T(25),
            fractional_cn2=T[1], wind_speed=T[1],
            wind_direction_deg=T[angle], altitude=T[0], T=T, kwargs...)
        rng = MersenneTwister(42)
        advance_by!(atm, zero(T); rng)
        layer = first(atm.layers)
        phase = A._direction_batch_layer_screen(layer)
        m = size(phase, 1)
        initial = T[2i + 7j for i in 1:m, j in 1:m]
        copyto!(phase, backend(initial))
        advance_by!(atm, T(1.25); rng)
        A.render_layer!(out, layer, zero(T), zero(T), one(T))
        c = T(m + 1) / 2
        opd_scale = A._direction_batch_layer_opd_scale(layer)
        expected = T[opd_scale *
            (2(c + i - 3 - dx) + 7(c + j - 3 - dy))
            for i in 1:5, j in 1:5]
        @test Array(out) ≈ expected rtol=T(1e-5)
        if constructor === InfiniteMultiLayerAtmosphere
            evolved = Array(A._direction_batch_layer_screen(layer))
            model = if dx > 0
                layer.screen.state.column_positive
            elseif dx < 0
                layer.screen.state.column_negative
            elseif dy > 0
                layer.screen.state.row_positive
            else
                layer.screen.state.row_negative
            end
            coords = Array(model.stencil.stencil_coords)
            @test Array(layer.screen.state.stencil_buffer) ==
                T[initial[coords[k, 1], coords[k, 2]] for k in axes(coords, 1)]
            if dx > 0
                @test evolved[2:end, :] == initial[1:end-1, :]
            elseif dx < 0
                @test evolved[1:end-1, :] == initial[2:end, :]
            elseif dy > 0
                @test evolved[:, 2:end] == initial[:, 1:end-1]
            else
                @test evolved[:, 1:end-1] == initial[:, 2:end]
            end
        end
    end

    # The prepared direction batch has its own fused accelerator sampler.
    source_groups = ([Source(band=:I, magnitude=0, T=T),
        Source(band=:I, magnitude=0, separation_arcsec=10,
            position_angle_deg=30, T=T)],
        [
        LGSSource(magnitude=0, separation_arcsec=15,
            position_angle_deg=120, altitude=T(2000), T=T)])
    for constructor in (MultiLayerAtmosphere, InfiniteMultiLayerAtmosphere),
            sources in source_groups
        kwargs = constructor === InfiniteMultiLayerAtmosphere ?
            (; screen_resolution=9, stencil_size=17) : (;)
        atm = constructor(tel; r0=T(0.2),
            reference_wavelength_m=T(500e-9), L0=T(25),
            fractional_cn2=T[1], wind_speed=T[0],
            wind_direction_deg=T[0], altitude=T[1000], T=T, kwargs...)
        epoch = advance_by!(atm, zero(T); rng=MersenneTwister(43))
        layer = first(atm.layers)
        phase = A._direction_batch_layer_screen(layer)
        m = size(phase, 1)
        copyto!(phase, backend(T[2i + 7j for i in 1:m, j in 1:m]))
        output = backend(zeros(T, 5, 5, length(sources)))
        batch = prepare_atmosphere_direction_batch(atm, tel, sources, output)
        actual = Array(render_atmosphere_directions!(batch, atm, epoch))
        pupil = Array(pupil_mask(tel))
        for (direction, src) in enumerate(sources)
            ax, ay = AdaptiveOpticsSim.Optics.coordinates_xy_arcsec(src)
            height = AdaptiveOpticsSim.Optics.source_height_m(src)
            scale = isfinite(height) ? T((height - 1000) / height) : one(T)
            sx = T((m + 1) / 2 + ax * pi / (180 * 3600) * 1000)
            sy = T((m + 1) / 2 + ay * pi / (180 * 3600) * 1000)
            opd_scale = A._direction_batch_layer_opd_scale(layer)
            expected = T[pupil[i, j] * opd_scale *
                (2(sx + scale * (i - 3)) + 7(sy + scale * (j - 3)))
                for i in 1:5, j in 1:5]
            @test actual[:, :, direction] ≈ expected rtol=T(1e-5)
        end
    end

    # Gather uses the same (x,y) geometry as the conditional Gaussian operator.
    coords_host = [2 5; 7 3; 4 8]
    coords = backend(coords_host)
    gathered = backend(zeros(T, 3))
    Backends.launch_kernel!(kernel_style, A.gather_stencil_data_kernel!,
        gathered, screen, coords, 3; ndrange=3)
    @test Array(gathered) == T[39, 35, 64]
    boundary = backend(T.(101:117))
    injected = similar(screen)
    for (kernel, axis, sign) in
            ((A.inject_column_positive_kernel!, 1, 1),
             (A.inject_column_negative_kernel!, 1, -1),
             (A.inject_row_positive_kernel!, 2, 1),
             (A.inject_row_negative_kernel!, 2, -1))
        Backends.launch_kernel!(kernel_style, kernel, injected,
            screen, boundary, 17; ndrange=size(screen))
        expected = circshift(screen_host, axis == 1 ? (sign, 0) : (0, sign))
        selectdim(expected, axis, sign > 0 ? 1 : 17) .= T.(101:117)
        @test Array(injected) == expected
    end
    return nothing
end
