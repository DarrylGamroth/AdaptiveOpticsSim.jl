# Shared by the KA CPU and optional accelerator targets. The scalar implementation
# remains the independent adjacent-difference reference for exceptional inputs.
function test_geometric_slopes_numerics(array_backend, style)
    wfs = AdaptiveOpticsSim.WavefrontSensors
    scalar = AdaptiveOpticsSim.Backends.ScalarCPUStyle()
    @testset "geometric slopes numerical cases" begin
        for T in (Float32, Float64), sub in (1, 2, 3, 4, 5, 16, 32)
            n = 2sub
            mask = Bool[true true; false true]
            device_mask = array_backend(mask)
            output = array_backend(zeros(T, 8))
            scale_x, scale_y = T(0.25), T(4)
            plane = T[2x - 3y for x in 1:n, y in 1:n]
            spike = zeros(T, n, n)
            spike[clamp(sub ÷ 2 + 1, 1, sub), clamp(sub ÷ 2 + 1, 1, sub)] = T(16)
            checkerboard = T[isodd(x + y) ? -8 : 8 for x in 1:n, y in 1:n]
            cancellation = T[(isodd(x) ? 1024 : -1024) + y / 8 for x in 1:n, y in 1:n]
            offset_plane = T[2^20 + x / 4 - y / 2 for x in 1:n, y in 1:n]
            subnormal_plane = T[nextfloat(zero(T)) * (x - y) for x in 1:n, y in 1:n]
            for opd in (fill(T(3), n, n), plane, spike, checkerboard,
                    cancellation, offset_plane, subnormal_plane)
                wfs._geometric_slopes!(style, output, array_backend(opd), device_mask,
                    sub, 2, 4, scale_x, scale_y)
                actual = Array(output)
                # Dyadic fixtures have exact neighbor differences and exact sums
                # in both precisions. A BigFloat sum defines the mathematical
                # gradient independently of either production implementation.
                for j in 1:2, i in 1:2
                    idx = i + 2(j - 1)
                    if !mask[i, j] || sub == 1
                        @test actual[idx] == actual[idx + 4] == 0
                        continue
                    end
                    xs, ys = (i - 1) * sub + 1, (j - 1) * sub + 1
                    sx, sy = BigFloat(0), BigFloat(0)
                    for y in ys:(ys + sub - 1), x in xs:(xs + sub - 2)
                        sx += BigFloat(opd[x + 1, y]) - BigFloat(opd[x, y])
                    end
                    for y in ys:(ys + sub - 2), x in xs:(xs + sub - 1)
                        sy += BigFloat(opd[x, y + 1]) - BigFloat(opd[x, y])
                    end
                    @test actual[idx] == T(BigFloat(scale_x) * sx / (sub * (sub - 1)))
                    @test actual[idx + 4] == T(BigFloat(scale_y) * sy / (sub * (sub - 1)))
                end
            end
            # Invalid subapertures must suppress exceptional data; sub=1 must
            # return zero without reading the OPD, even for valid subapertures.
            for value in (T(NaN), T(Inf), T(-Inf)), location in ((1, 1), (min(3, sub), min(3, sub)))
                opd = zeros(T, n, n)
                opd[location...] = value
                opd[(sub + 1):end, 1:sub] .= value
                expected = zeros(T, 8)
                wfs._geometric_slopes!(scalar, expected, opd, mask, sub, 2, 4,
                    scale_x, scale_y)
                wfs._geometric_slopes!(style, output, array_backend(opd), device_mask,
                    sub, 2, 4, scale_x, scale_y)
                @test isequal(Array(output), expected)
            end
            if sub > 2
                # Finite differences can overflow a partial sum even when the
                # mathematical endpoint cancellation is zero.
                for pattern in (T[0, 1, 0, -1, 0], T[1, -1, 1, -1, 1])
                    opd = T[pattern[mod1(x, 5)] * (floatmax(T) / 2) for x in 1:n, y in 1:n]
                    expected = zeros(T, 8)
                    wfs._geometric_slopes!(scalar, expected, opd, mask, sub, 2, 4,
                        scale_x, scale_y)
                    wfs._geometric_slopes!(style, output, array_backend(opd), device_mask,
                        sub, 2, 4, scale_x, scale_y)
                    @test isequal(Array(output), expected)
                end
            end
            wfs._geometric_slopes!(style, output, array_backend(fill(T(NaN), n, n)),
                array_backend(falses(2, 2)), sub, 2, 4, scale_x, scale_y)
            @test all(iszero, Array(output))
        end
        @testset "non-dyadic cooperative accumulation" begin
            for T in (Float32, Float64), sub in (16, 32)
                # Equal periodic endpoints make the gradient small compared with
                # the sum of the interior differences that must cancel.
                for dc in (zero(T), T(1024))
                    opd = T[dc + sin(2π * (x - 1) / (sub - 1)) +
                        cos(2π * (y - 1) / (sub - 1)) + 1e-4x - 3e-4y
                        for x in 1:sub, y in 1:sub]
                    valid = trues(1, 1)
                    output = array_backend(zeros(T, 2))
                    serial_output = zeros(T, 2)
                    scales = (T(0.25), T(4))
                    wfs._geometric_slopes!(scalar, serial_output, opd, valid,
                        sub, 1, 1, scales...)
                    wfs._geometric_slopes!(style, output, array_backend(opd), array_backend(valid),
                        sub, 1, 1, scales...)
                    actual = Array(output)
                    count = sub * (sub - 1)
                    for axis in 1:2
                        exact = BigFloat(0)
                        magnitude = BigFloat(0)
                        for k in 0:(count - 1)
                            width = axis == 1 ? sub - 1 : sub
                            x, y = k % width + 1, k ÷ width + 1
                            a = axis == 1 ? opd[x + 1, y] : opd[x, y + 1]
                            b = opd[x, y]
                            exact += BigFloat(a) - BigFloat(b)
                            magnitude += abs(BigFloat(a - b))
                        end
                        exact *= BigFloat(scales[axis]) / count
                        # One subtraction, ceil(count/32) lane additions, five
                        # tree stages, and multiply/divide. Unit roundoff is ε/2.
                        steps = cld(count, 32) + 8
                        u = BigFloat(eps(T)) / 2
                        gamma = steps * u / (1 - steps * u)
                        budget = gamma * magnitude * abs(BigFloat(scales[axis])) / count
                        @test abs(BigFloat(actual[axis]) - exact) <= budget
                        # The serial reference has a longer accumulation chain.
                        serial_gamma = (count + 3) * u / (1 - (count + 3) * u)
                        serial_budget = serial_gamma * magnitude * abs(BigFloat(scales[axis])) / count
                        @test abs(BigFloat(serial_output[axis]) - exact) <= serial_budget
                    end
                end
            end
        end
        @testset "scale extremes and public validation" begin
            for T in (Float32, Float64)
                opd = T[x / 8 - y / 4 for x in 1:10, y in 1:10]
                valid = trues(2, 2)
                expected = zeros(T, 8)
                output = array_backend(zeros(T, 8))
                for scale in (floatmax(T), T(Inf), T(NaN), zero(T))
                    wfs._geometric_slopes!(scalar, expected, opd, valid, 5, 2, 4, scale, one(T))
                    wfs._geometric_slopes!(style, output, array_backend(opd), array_backend(valid),
                        5, 2, 4, scale, one(T))
                    @test isequal(Array(output), expected)
                end
                wfs.geometric_wavefront_slopes!(expected, opd, valid, (T(4), T(0.25)))
                @test expected[1] == T(1 / 32)
                @test expected[5] == T(-1)
                @test @allocated(wfs.geometric_wavefront_slopes!(expected, opd, valid,
                    (T(4), T(0.25)))) == 0
                @test_throws AdaptiveOpticsSim.InvalidConfiguration wfs.geometric_wavefront_slopes!(
                    expected, opd, valid, (zero(T), one(T)))
                @test_throws AdaptiveOpticsSim.DimensionMismatchError wfs.geometric_slopes!(
                    zeros(T, 7), opd, valid)
            end
        end
        @testset "output-grid crossover" begin
            for T in (Float32, Float64), n_sub in (32, 33)
                sub = 5
                n = sub * n_sub
                opd = T[x / 8 - y / 4 for x in 1:n, y in 1:n]
                valid = Bool[mod(i + 2j, 5) != 0 for i in 1:n_sub, j in 1:n_sub]
                expected = zeros(T, 2n_sub^2)
                output = array_backend(similar(expected))
                wfs._geometric_slopes!(scalar, expected, opd, valid,
                    sub, n_sub, n_sub^2, T(0.25), T(4))
                wfs._geometric_slopes!(style, output, array_backend(opd), array_backend(valid),
                    sub, n_sub, n_sub^2, T(0.25), T(4))
                @test Array(output) == expected
            end
        end
        @testset "mixed precision and edge masks" begin
            for T in (Float32, Float64), S in (Float32, Float64)
                sub = 5
                opd = S[sin(x + y / 3) for x in 1:10, y in 1:10]
                valid = trues(2, 2)
                expected = zeros(T, 8)
                output = array_backend(zeros(T, 8))
                wfs._geometric_slopes!(scalar, expected, opd, valid, sub, 2, 4,
                    T(0.25), T(4))
                wfs._geometric_slopes!(style, output, array_backend(opd), array_backend(valid),
                    sub, 2, 4, T(0.25), T(4))
                @test Array(output) ≈ expected rtol=8eps(T) atol=8eps(T)
                edges = Bool[isodd(x + y) for x in 1:10, y in 1:10]
                wfs._edge_geometric_slopes!(scalar, expected, opd, valid, edges, sub, 2, 4)
                wfs._edge_geometric_slopes!(style, output, array_backend(opd),
                    array_backend(valid), array_backend(edges), sub, 2, 4)
                @test Array(output) == expected
            end
        end
    end
end
