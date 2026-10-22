# Standalone focused checks; run in an environment containing CUDA and this checkout.
using Test
using AdaptiveOpticsSim
using AdaptiveOpticsSim.Optics
using FixedSizeArrays: FixedSizeVector
using AdaptiveOpticsSim.Backends
using AdaptiveOpticsSim.WavefrontSensors
import CUDA

CUDA.allowscalar(false)

const ModulationOptics = AdaptiveOpticsSim.Optics

function modulation_storage_fixture(::Type{T}, selector) where {T}
    telescope = Telescope(resolution=16, diameter=T(8), T=T, backend=selector)
    pupil = PupilFunction(telescope; T=T, backend=selector)
    amplitude = Array(pupil.amplitude)
    opd = similar(amplitude)
    for j in 1:16, i in 1:16
        x, y = T(2(i - 1) / 15 - 1), T(2(j - 1) / 15 - 1)
        amplitude[i, j] *= T(0.8) + T(0.1) * x + T(0.05) * y
        opd[i, j] = T(40e-9) * (sinpi(x) + T(0.4) * cospi(y) + T(0.2) * x * y)
    end
    copyto!(pupil.amplitude, amplitude)
    copyto!(pupil.opd, opd)
    source = Source(band=:custom, wavelength=T(0.75e-6),
        photon_irradiance=T(10), T=T)
    return telescope, pupil, source
end

function modulation_storage_sensor(telescope, ::Type{T}, strategy;
    geometry=(;), policy=(;)) where {T}
    return PyramidWFS(telescope; pupil_samples=4, diffraction_padding=4,
        modulation_propagation_strategy=strategy, T=T, geometry..., policy...)
end

function modulation_storage_diagnostic(sensor, pupil, source)
    out = similar(WavefrontSensors.pyramid_propagation_workspace(sensor).intensity)
    pyramid_modulation_frame!(out, sensor, pupil, source)
    return Array(out)
end

function modulation_storage_rate(sensor, pupil, source)
    front = PyramidOpticalFrontEnd(sensor, source)
    rate = pyramid_rate_map(front, pupil)
    prepared = prepare_wfs_optics(front, pupil, rate)
    form_wfs_optical_products!(rate, pupil, prepared)
    return Array(rate.values)
end

function modulation_storage_arrays(sensor)
    owner = WavefrontSensors.pyramid_operating_modulation(sensor)
    batch = WavefrontSensors.pyramid_propagation_workspace(sensor).modulation_batch
    masks = if hasproperty(batch, :shifted_masks)
        (batch.shifted_masks,)
    else
        (batch.axis_1_factors, batch.axis_2_factors)
    end
    return (owner.coordinates, owner.offsets, owner.amplitude_weights,
        batch.field_stack, batch.operating_weights, batch.axis_1_shifts_rad,
        batch.axis_2_shifts_rad, masks...)
end

@testset "Compact Pyramid modulation storage" begin
    @test CUDA.functional()
    for T in (Float32, Float64)
        cpu_tel, cpu_pupil, source = modulation_storage_fixture(T, CPUBackend())
        gpu_tel, gpu_pupil, _ = modulation_storage_fixture(T, CUDABackend())
        rtol = T === Float32 ? T(3e-5) : T(2e-12)
        atol = T === Float32 ? T(1e-6) : T(1e-13)
        policies = (
            (modulation=zero(T),),
            (modulation=T(0.7), modulation_points=4,
                modulation_phase_offset_rad=T(0.31)),
            (modulation=zero(T), user_modulation_path=(
                (T(1.3), T(0.7), T(0)), (T(-0.4), T(1.1), T(2)),
                (T(-0.9), T(-0.5), T(1)))),
        )
        geometries = ((;), (rooftop=T(0.25), phase_mask_rotation_rad=T(0.1)))
        for geometry in geometries, policy in policies
            @testset "$T geometry=$geometry policy=$policy" begin
                cpu = modulation_storage_sensor(cpu_tel, T, WavefrontSensors.PyramidShiftedMaskStrategy();
                    geometry, policy)
                gpu = modulation_storage_sensor(gpu_tel, T, WavefrontSensors.PyramidShiftedMaskStrategy();
                    geometry, policy)
                dense_cpu = modulation_storage_sensor(cpu_tel, T, WavefrontSensors.PyramidPupilTiltStrategy();
                    geometry, policy)
                dense_gpu = modulation_storage_sensor(gpu_tel, T, WavefrontSensors.PyramidPupilTiltStrategy();
                    geometry, policy)
                for (sensor, on_gpu) in ((cpu, false), (gpu, true))
                    owner = WavefrontSensors.pyramid_operating_modulation(sensor)
                    @test !hasproperty(owner, :phases)
                    @test size(owner.offsets) == (2, length(owner.amplitude_weights))
                    @test length(owner.coordinates) == 16
                    @test (owner.coordinates isa CUDA.CuArray) == on_gpu
                    @test (owner.offsets isa CUDA.CuArray) == on_gpu
                    @test owner.amplitude_weights isa FixedSizeVector{T}
                    @test Array(owner.coordinates) == collect(range(-T(pi), T(pi); length=16))
                    expected_offsets = hcat((collect(ModulationOptics.modulation_offset(
                        owner.policy, point, T)) for point in eachindex(owner.amplitude_weights))...)
                    @test Array(owner.offsets) == expected_offsets
                    @test sum(abs2, owner.amplitude_weights) ≈ one(T)
                    if haskey(policy, :user_modulation_path)
                        @test owner.amplitude_weights[1] == zero(T)
                        @test owner.amplitude_weights ≈ sqrt.(T[0, 2 / 3, 1 / 3])
                    end
                    batch = WavefrontSensors.pyramid_propagation_workspace(sensor).modulation_batch
                    @test (batch isa WavefrontSensors.PyramidSeparableShiftedMaskModulationWorkspace) == isempty(geometry)
                end
                cpu_batch = WavefrontSensors.pyramid_propagation_workspace(cpu).modulation_batch
                gpu_batch = WavefrontSensors.pyramid_propagation_workspace(gpu).modulation_batch
                for name in (isempty(geometry) ? (:axis_1_factors, :axis_2_factors) : (:shifted_masks,))
                    @test Array(getproperty(gpu_batch, name)) ≈
                        Array(getproperty(cpu_batch, name)) rtol=rtol atol=atol
                end
                cpu_diag = modulation_storage_diagnostic(cpu, cpu_pupil, source)
                gpu_diag = modulation_storage_diagnostic(gpu, gpu_pupil, source)
                @test cpu_diag ≈ modulation_storage_diagnostic(dense_cpu, cpu_pupil, source) rtol=rtol atol=atol
                @test gpu_diag ≈ modulation_storage_diagnostic(dense_gpu, gpu_pupil, source) rtol=rtol atol=atol
                @test gpu_diag ≈ cpu_diag rtol=rtol atol=atol
                @test modulation_storage_rate(gpu, gpu_pupil, source) ≈
                    modulation_storage_rate(cpu, cpu_pupil, source) rtol=rtol atol=atol
            end
        end

        for (tel, pupil) in ((cpu_tel, cpu_pupil), (gpu_tel, gpu_pupil)), geometry in geometries
            @testset "$T update $(backend(tel)) geometry=$geometry" begin
                policy = (modulation=T(0.7), modulation_points=4,
                    modulation_phase_offset_rad=T(0.31))
                sensor = modulation_storage_sensor(tel, T, WavefrontSensors.PyramidShiftedMaskStrategy(); geometry, policy)
                arrays = modulation_storage_arrays(sensor)
                addresses = map(pointer, arrays)
                weights = copy(WavefrontSensors.pyramid_operating_modulation(sensor).amplitude_weights)
                for (radius, enabled) in ((T(1.1), true), (T(0.5), false), (T(0.3), true))
                    @test ModulationOptics.update_cycle_averaged_circular_modulation!(
                        sensor, radius; enabled) === sensor
                    @test map(pointer, modulation_storage_arrays(sensor)) == addresses
                    @test all(a === b for (a, b) in zip(arrays, modulation_storage_arrays(sensor)))
                    @test WavefrontSensors.pyramid_operating_modulation(sensor).amplitude_weights == weights
                    fresh_policy = (modulation=enabled ? radius : zero(T), modulation_points=4,
                        modulation_phase_offset_rad=T(0.31))
                    # Zero-radius fresh construction selects NoModulation; the four
                    # identical updated points have the same normalized optical sum.
                    fresh = modulation_storage_sensor(tel, T, WavefrontSensors.PyramidShiftedMaskStrategy();
                        geometry, policy=fresh_policy)
                    @test modulation_storage_diagnostic(sensor, pupil, source) ≈
                        modulation_storage_diagnostic(fresh, pupil, source) rtol=rtol atol=atol
                    @test modulation_storage_rate(sensor, pupil, source) ≈
                        modulation_storage_rate(fresh, pupil, source) rtol=rtol atol=atol
                end
                for invalid in (T(-1), T(NaN), T(Inf), T === Float32 ? Float64(1) : Float32(1))
                    before = map(a -> copy(Array(a)), modulation_storage_arrays(sensor))
                    @test_throws InvalidConfiguration ModulationOptics.update_cycle_averaged_circular_modulation!(sensor, invalid)
                    @test isequal(map(a -> Array(a), modulation_storage_arrays(sensor)), before)
                    @test map(pointer, modulation_storage_arrays(sensor)) == addresses
                end
                @test_throws InvalidConfiguration ModulationOptics.update_cycle_averaged_circular_modulation!(
                    WavefrontSensors.pyramid_operating_modulation(sensor), T(0.4))
                for rejected_policy in (policies[1], policies[3])
                    rejected = modulation_storage_sensor(tel, T, WavefrontSensors.PyramidShiftedMaskStrategy();
                        geometry, policy=rejected_policy)
                    before = map(a -> copy(Array(a)), modulation_storage_arrays(rejected))
                    @test_throws InvalidConfiguration ModulationOptics.update_cycle_averaged_circular_modulation!(rejected, T(0.4))
                    @test isequal(map(a -> Array(a), modulation_storage_arrays(rejected)), before)
                end
            end
        end
    end
end


function modulation_storage_replay!(captured, context)
    Backends._with_prepared_device_execution_context(context) do
        Backends._launch_prepared_device_graph!(captured, context)
        Backends._synchronize_prepared_device_execution_context_blocking!(context)
    end
    return nothing
end

function modulation_storage_replay_bytes(captured, context)
    return @allocated modulation_storage_replay!(captured, context)
end

@testset "Captured Pyramid rate observes in-place modulation updates" begin
    T = Float32
    telescope, pupil, source = modulation_storage_fixture(T, CUDABackend())
    geometries = ((;), (rooftop=T(0.25), phase_mask_rotation_rad=T(0.1)))
    for geometry in geometries
        @testset "geometry=$geometry" begin
            policy = (modulation=T(0.7), modulation_points=4,
                modulation_phase_offset_rad=T(0.31))
            sensor = modulation_storage_sensor(telescope, T,
                WavefrontSensors.PyramidShiftedMaskStrategy(); geometry, policy)
            front = PyramidOpticalFrontEnd(sensor, source)
            rate = pyramid_rate_map(front, pupil)
            prepared = prepare_wfs_optics(front, pupil, rate)
            # Complete setup and compile every enqueue before beginning capture.
            form_wfs_optical_products!(rate, pupil, prepared)
            context = Backends._prepare_device_execution_context(rate.values)
            CUDA.synchronize()
            captured = Backends._with_prepared_device_execution_context(context) do
                WavefrontSensors.enqueue_wfs_optical_products!(rate, pupil, prepared)
                Backends._synchronize_prepared_device_execution_context_blocking!(context)
                Backends._capture_prepared_device_graph(context) do
                    WavefrontSensors.enqueue_wfs_optical_products!(rate, pupil, prepared)
                end
            end
            arrays = modulation_storage_arrays(sensor)
            addresses = map(pointer, arrays)
            output_address = pointer(rate.values)
            modulation_storage_replay!(captured, context)
            initial_rate = Array(rate.values)
            for (radius, enabled) in ((T(1.1), true), (T(0.5), false), (T(0.3), true))
                ModulationOptics.update_cycle_averaged_circular_modulation!(
                    sensor, radius; enabled)
                @test map(pointer, modulation_storage_arrays(sensor)) == addresses
                @test pointer(rate.values) == output_address
                # The graph was captured before the update; its buffer bindings
                # must now see the rebuilt shifted masks without recapture.
                modulation_storage_replay!(captured, context)
                actual = Array(rate.values)
                fresh_policy = (modulation=enabled ? radius : zero(T),
                    modulation_points=4, modulation_phase_offset_rad=T(0.31))
                fresh = modulation_storage_sensor(telescope, T,
                    WavefrontSensors.PyramidShiftedMaskStrategy();
                    geometry, policy=fresh_policy)
                expected = modulation_storage_rate(fresh, pupil, source)
                @test actual ≈ expected rtol=3f-5 atol=1f-6
                @test !isapprox(actual, initial_rate; rtol=1f-4, atol=1f-6)
            end
            for _ in 1:16
                modulation_storage_replay!(captured, context)
            end
            modulation_storage_replay_bytes(captured, context)
            @test modulation_storage_replay_bytes(captured, context) == 0
        end
    end
end
