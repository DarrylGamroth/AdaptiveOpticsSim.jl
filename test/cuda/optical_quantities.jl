using Test
using CUDA
using AdaptiveOpticsSim
using AdaptiveOpticsSim.Optics
using AdaptiveOpticsSim.Backends: CUDABackend

CUDA.functional() || error("CUDA hardware is required for optical quantity tests")
CUDA.allowscalar(false)

@testset "CUDA nominal-axis optical quantities" begin
    for n in (7, 8), obstruction in (0f0, 0.4f0)
        cpu = Telescope(resolution=n, diameter=2, central_obstruction=obstruction, T=Float32)
        reflectivity = Float32[0.2 + 0.7 * (i + j) / (2n) for i in 1:n, j in 1:n]
        set_pupil_reflectivity!(cpu, reflectivity)
        host = PupilFunction(cpu)
        gpu = Telescope(resolution=n, diameter=2, central_obstruction=obstruction,
            T=Float32, backend=CUDABackend())
        set_pupil_reflectivity!(gpu, CuArray(reflectivity))
        pupil = PupilFunction(gpu)
        source = Source(band=:custom, wavelength=750f-9, photon_irradiance=3f0, T=Float32)
        prepared = prepare_direct_imaging(pupil, source; zero_padding=2)
        @test prepared.input.amplitude isa CuArray
        peak = diffraction_limited_peak(prepared)
        expected_energy = 3f0 * (2f0 / n)^2 * sum(abs2, host.amplitude)
        @test peak ≈ expected_energy * abs2(sum(host.amplitude)) /
            (length(prepared.field.values) * sum(abs2, host.amplitude)) rtol=2e-6
        for opd in (zeros(Float32, n, n), fill(250f-9, n, n),
            Float32[20f-9 * (i - (n + 1) / 2) for i in 1:n, j in 1:n])
            copyto!(pupil.opd, opd)
            # Independent host real/imaginary quadrature, no AOS CPU call.
            real_field = sum(host.amplitude .* cos.(2π .* Float64.(opd) ./ 750f-9))
            imag_field = sum(host.amplitude .* sin.(2π .* Float64.(opd) ./ 750f-9))
            expected = (real_field^2 + imag_field^2) / sum(host.amplitude)^2
            @test on_axis_strehl(pupil, 750f-9) ≈ expected rtol=3e-6
            @test diffraction_limited_peak(prepared) == peak
        end
        product = form_direct_image!(prepared)
        @test product.values isa CuArray
        @test sum(product.values) ≈ expected_energy rtol=5e-6
        @test product.metadata.normalization isa PhotonRateNormalization
        @test_throws InvalidConfiguration on_axis_strehl(pupil, -750f-9)
    end
end
