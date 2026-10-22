# Standalone: julia --project=<CUDA environment> test/cuda/modulation_batches.jl
using Test
using LinearAlgebra
using AdaptiveOpticsSim
using AdaptiveOpticsSim.Optics
using AdaptiveOpticsSim.Backends
using AdaptiveOpticsSim.WavefrontSensors
import CUDA
include(joinpath(@__DIR__, "..", "..", "benchmarks", "modulation_batch_override.jl"))

CUDA.allowscalar(false)
BLAS.set_num_threads(1)
Backends.set_fft_provider_threads!(1)

function batch_fixture(::Type{T}, selector, strategy, geometry, policy, batch) where {T}
    MODULATION_BATCH[] = batch
    telescope = Telescope(resolution=16, diameter=T(8), T=T, backend=selector)
    pupil = PupilFunction(telescope; T=T, backend=selector)
    amplitude = Array(pupil.amplitude)
    modes = Array{T}(undef, 16, 16, 3)
    for j in 1:16, i in 1:16
        x, y = T(2(i - 1) / 15 - 1), T(2(j - 1) / 15 - 1)
        amplitude[i, j] *= T(0.8) + T(0.1)*x + T(0.05)*y
        modes[i, j, 1] = x
        modes[i, j, 2] = y
        modes[i, j, 3] = x*y
    end
    copyto!(pupil.amplitude, amplitude)
    source = Source(band=:custom, wavelength=T(0.75e-6), photon_irradiance=T(10), T=T)
    sensor = PyramidWFS(telescope; pupil_samples=4, diffraction_padding=4,
        modulation_propagation_strategy=strategy, T=T, geometry..., policy...)
    front = PyramidOpticalFrontEnd(sensor, source)
    rate = pyramid_rate_map(front, pupil)
    prepared = prepare_wfs_optics(front, pupil, rate)
    return (; pupil, modes, source, sensor, rate, prepared)
end

function batch_rate(fixture, coefficients)
    # Coefficients are in nm. Every backend receives exactly the same host OPD.
    opd = dropdims(sum(fixture.modes .* reshape(coefficients, 1, 1, 3); dims=3); dims=3)
    opd .*= eltype(opd)(1e-9)
    copyto!(fixture.pupil.opd, opd)
    form_wfs_optical_products!(fixture.rate, fixture.pupil, fixture.prepared)
    return vec(Array(fixture.rate.values))
end

function batch_response(fixture)
    T = eltype(fixture.modes)
    response = Matrix{T}(undef, length(fixture.rate.values), 3)
    for mode in 1:3
        coefficients = zeros(T, 3)
        coefficients[mode] = T(2)
        plus = batch_rate(fixture, coefficients)
        coefficients[mode] = -T(2)
        response[:, mode] .= (plus .- batch_rate(fixture, coefficients)) ./ T(4)
    end
    return response
end

function batch_feedback(fixture, reference, reconstructor)
    T = eltype(fixture.modes)
    residual = T[12, -9, 6]
    trajectory = Matrix{T}(undef, 3, 5)
    for step in 1:5
        trajectory[:, step] .= residual
        residual .-= T(0.3) .* (reconstructor * (batch_rate(fixture, residual) .- reference))
    end
    return trajectory
end

function batch_replay!(captured, context)
    Backends._with_prepared_device_execution_context(context) do
        Backends._launch_prepared_device_graph!(captured, context)
        Backends._synchronize_prepared_device_execution_context_blocking!(context)
    end
    return nothing
end
batch_replay_bytes(captured, context) = @allocated batch_replay!(captured, context)

function batch_capture_checks(fixture, expected; rtol, atol)
    context = Backends._prepare_device_execution_context(fixture.rate.values)
    CUDA.synchronize()
    captured = Backends._with_prepared_device_execution_context(context) do
        WavefrontSensors.enqueue_wfs_optical_products!(fixture.rate, fixture.pupil, fixture.prepared)
        Backends._synchronize_prepared_device_execution_context_blocking!(context)
        Backends._capture_prepared_device_graph(context) do
            WavefrontSensors.enqueue_wfs_optical_products!(fixture.rate, fixture.pupil, fixture.prepared)
        end
    end
    output_pointer = pointer(fixture.rate.values)
    # Poison both the caller-visible output and FFT scratch. Replay must reset
    # accumulation and fully overwrite its scratch before consuming it.
    scratch = WavefrontSensors.pyramid_propagation_workspace(fixture.sensor).modulation_batch.field_stack
    for _ in 1:3
        fill!(fixture.rate.values, eltype(fixture.rate.values)(NaN))
        fill!(scratch, eltype(scratch)(NaN))
        CUDA.synchronize()
        batch_replay!(captured, context)
        @test vec(Array(fixture.rate.values)) ≈ expected rtol=rtol atol=atol
        @test pointer(fixture.rate.values) == output_pointer
    end
    for _ in 1:16
        batch_replay!(captured, context)
    end
    batch_replay_bytes(captured, context)
    @test batch_replay_bytes(captured, context) == 0
end

@testset "Pyramid modulation batch numerical and capture invariance" begin
    @test CUDA.functional()
    for T in (Float32, Float64)
        rtol, atol = T === Float32 ? (5f-5, 2f-6) : (3e-12, 2e-13)
        response_rtol = T === Float32 ? 8f-4 : 5e-11
        policies = (
            (modulation=T(0.7), modulation_points=6, modulation_phase_offset_rad=T(0.31)),
            (modulation=zero(T), user_modulation_path=(
                (T(1.3), T(0.7), T(0)), (T(-0.4), T(1.1), T(2)),
                (T(-0.9), T(-0.5), T(1)), (T(0.2), T(-1.2), T(3)),
                (T(0.6), T(0.3), T(0)), (T(-0.8), T(0.2), T(4)))),
        )
        cases = (
            (WavefrontSensors.PyramidPupilTiltStrategy(), (;)),
            (WavefrontSensors.PyramidPupilTiltStrategy(), (rooftop=T(0.25), phase_mask_rotation_rad=T(0.1))),
        )
        for (strategy, geometry) in cases, policy in policies
            @testset "$T strategy=$(typeof(strategy)) geometry=$geometry policy=$policy" begin
                cpu = batch_fixture(T, CPUBackend(), strategy, geometry, policy, 0)
                baseline = batch_fixture(T, CUDABackend(), strategy, geometry, policy, 1)
                coefficients = T[30, -20, 15]
                cpu_rate = batch_rate(cpu, coefficients)
                baseline_rate = batch_rate(baseline, coefficients)
                @test baseline_rate ≈ cpu_rate rtol=rtol atol=atol
                expected_flux = photon_irradiance(cpu.source) *
                    prod(cpu.pupil.metadata.sampling) * sum(abs2, Array(cpu.pupil.amplitude))
                @test sum(cpu_rate) ≈ expected_flux rtol=rtol
                reference = batch_rate(cpu, zeros(T, 3))
                cpu_response = batch_response(cpu)
                baseline_response = batch_response(baseline)
                @test baseline_response ≈ cpu_response rtol=response_rtol atol=atol
                # Freeze the same CPU calibration and gain for every batch.
                reconstructor = pinv(cpu_response)
                cpu_trajectory = batch_feedback(cpu, reference, reconstructor)
                baseline_trajectory = batch_feedback(baseline, reference, reconstructor)
                @test baseline_trajectory ≈ cpu_trajectory rtol=response_rtol atol=T(0.002)
                @test norm(cpu_trajectory[:, end]) < norm(cpu_trajectory[:, 1])
                for batch in (1, 2, 3, 6)
                    @testset "batch=$batch" begin
                        fixture = batch_fixture(T, CUDABackend(), strategy, geometry, policy, batch)
                        workspace = WavefrontSensors.pyramid_propagation_workspace(fixture.sensor).modulation_batch
                        @test workspace.batch_size == batch
                        actual = batch_rate(fixture, coefficients)
                        @test actual ≈ baseline_rate rtol=rtol atol=atol
                        @test actual ≈ cpu_rate rtol=rtol atol=atol
                        @test sum(actual) ≈ expected_flux rtol=rtol
                        @test batch_response(fixture) ≈ baseline_response rtol=response_rtol atol=atol
                        @test batch_feedback(fixture, reference, reconstructor) ≈ baseline_trajectory rtol=response_rtol atol=T(0.002)
                        expected = batch_rate(fixture, coefficients)
                        batch_capture_checks(fixture, expected; rtol, atol)
                    end
                end
            end
        end
    end
    MODULATION_BATCH[] = 0
end
