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

function modulation_storage_sensor(telescope, ::Type{T};
    geometry=(;), policy=(;)) where {T}
    return PyramidWFS(telescope; pupil_samples=4, diffraction_padding=4,
        T=T, geometry..., policy...)
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

modulation_storage_batch_arrays(::WavefrontSensors.NoPyramidModulationBatchWorkspace) = ()
modulation_storage_fft_arrays(::WavefrontSensors.NoPyramidModulationFFTWorkspace) = ()
modulation_storage_fft_arrays(workspace::WavefrontSensors.PyramidTransposeFFTWorkspace) =
    (workspace.field_stack,)
modulation_storage_batch_arrays(batch::WavefrontSensors.PyramidModulationBatchWorkspace) =
    (batch.field_stack, batch.pupil_phase, batch.operating_weights,
        modulation_storage_fft_arrays(batch.fft_workspace)...)

modulation_storage_phase_arrays(owner::ModulationOptics.PreparedFocalPlaneModulation) =
    (owner.phases,)
modulation_storage_phase_arrays(owner::ModulationOptics.PreparedSeparableFocalPlaneModulation) =
    (owner.axis_1_factors, owner.axis_2_factors)

modulation_storage_dense_phases(owner::ModulationOptics.PreparedFocalPlaneModulation) =
    Array(owner.phases)
function modulation_storage_dense_phases(owner::ModulationOptics.PreparedSeparableFocalPlaneModulation)
    axis_1 = Array(owner.axis_1_factors)
    axis_2 = Array(owner.axis_2_factors)
    resolution, points = size(axis_1)
    return reshape(axis_1, resolution, 1, points) .*
        reshape(axis_2, 1, resolution, points)
end

function modulation_storage_check_phases(owner, ::Type{T};
    policy=owner.policy, rtol, atol) where {T}
    points = length(owner.amplitude_weights)
    coordinates = collect(range(-T(pi), T(pi); length=16))
    expected = Array{Complex{T}}(undef, 16, 16, points)
    for point in 1:points
        x, y = ModulationOptics.modulation_offset(policy, point, T)
        for j in 1:16, i in 1:16
            expected[i, j, point] = cis(x * coordinates[j] + y * coordinates[i])
        end
    end
    @test modulation_storage_dense_phases(owner) ≈ expected rtol=rtol atol=atol
end

function modulation_storage_arrays(sensor)
    owner = WavefrontSensors.pyramid_operating_modulation(sensor)
    batch = WavefrontSensors.pyramid_propagation_workspace(sensor).modulation_batch
    return (modulation_storage_phase_arrays(owner)..., owner.amplitude_weights,
        modulation_storage_batch_arrays(batch)...)
end

@testset "Physical Pyramid modulation phase storage" begin
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
        geometries = ((;), (rooftop=T(0.25), phase_mask_rotation_rad=T(0.1)),
            (phase_mask_rotation_rad=T(0.17),
                pupil_shift_x_pixels=(T(0.15),T(-0.1),T(0.2),T(-0.25)),
                pupil_shift_y_pixels=(T(-0.2),T(0.05),T(0.1),T(0.3))))
        for geometry in geometries, policy in policies
            @testset "$T geometry=$geometry policy=$policy" begin
                cpu = modulation_storage_sensor(cpu_tel, T; geometry, policy)
                gpu = modulation_storage_sensor(gpu_tel, T; geometry, policy)
                for (sensor, on_gpu) in ((cpu, false), (gpu, true))
                    owner = WavefrontSensors.pyramid_operating_modulation(sensor)
                    points = length(owner.amplitude_weights)
                    if on_gpu
                        @test owner isa ModulationOptics.PreparedSeparableFocalPlaneModulation
                        @test !hasproperty(owner, :phases)
                        @test size(owner.axis_1_factors) == (16, points)
                        @test size(owner.axis_2_factors) == (16, points)
                        @test owner.axis_1_factors isa CUDA.CuArray
                        @test owner.axis_2_factors isa CUDA.CuArray
                        @test length(owner.axis_1_factors) + length(owner.axis_2_factors) ==
                            2 * 16 * points
                    else
                        @test owner isa ModulationOptics.PreparedFocalPlaneModulation
                        @test size(owner.phases) == (16, 16, points)
                        @test owner.phases isa Array
                    end
                    @test owner.amplitude_weights isa FixedSizeVector{T}
                    modulation_storage_check_phases(owner, T; rtol, atol)
                    @test sum(abs2, owner.amplitude_weights) ≈ one(T)
                    if haskey(policy, :user_modulation_path)
                        @test owner.amplitude_weights[1] == zero(T)
                        @test owner.amplitude_weights ≈ sqrt.(T[0, 2 / 3, 1 / 3])
                    end
                end
                gpu_owner = WavefrontSensors.pyramid_operating_modulation(gpu)
                cpu_owner = WavefrontSensors.pyramid_operating_modulation(cpu)
                @test modulation_storage_dense_phases(gpu_owner) ≈
                    cpu_owner.phases rtol=rtol atol=atol
                gpu_batch = WavefrontSensors.pyramid_propagation_workspace(gpu).modulation_batch
                @test gpu_batch isa WavefrontSensors.PyramidModulationBatchWorkspace
                @test gpu_batch.field_stack isa CUDA.CuArray
                @test gpu_batch.pupil_phase isa CUDA.CuArray
                @test size(gpu_batch.pupil_phase) == (16, 16)
                @test gpu_batch.operating_weights isa CUDA.CuArray
                @test Array(gpu_batch.operating_weights) == collect(cpu_owner.amplitude_weights)
                @test modulation_storage_diagnostic(gpu, gpu_pupil, source) ≈
                    modulation_storage_diagnostic(cpu, cpu_pupil, source) rtol=rtol atol=atol
                @test modulation_storage_rate(gpu, gpu_pupil, source) ≈
                    modulation_storage_rate(cpu, cpu_pupil, source) rtol=rtol atol=atol
            end
        end

        for (tel, pupil) in ((cpu_tel, cpu_pupil), (gpu_tel, gpu_pupil)), geometry in geometries
            @testset "$T update $(backend(tel)) geometry=$geometry" begin
                policy = (modulation=T(0.7), modulation_points=4,
                    modulation_phase_offset_rad=T(0.31))
                sensor = modulation_storage_sensor(tel, T; geometry, policy)
                arrays = modulation_storage_arrays(sensor)
                addresses = map(pointer, arrays)
                weights = copy(WavefrontSensors.pyramid_operating_modulation(sensor).amplitude_weights)
                for (radius, enabled) in ((T(1.1), true), (T(0.5), false), (T(0.3), true))
                    @test ModulationOptics.update_cycle_averaged_circular_modulation!(
                        sensor, radius; enabled) === sensor
                    @test map(pointer, modulation_storage_arrays(sensor)) == addresses
                    @test all(a === b for (a, b) in zip(arrays, modulation_storage_arrays(sensor)))
                    @test WavefrontSensors.pyramid_operating_modulation(sensor).amplitude_weights == weights
                    updated_policy = CircularModulation(enabled ? radius : zero(T);
                        samples=4, phase_offset_rad=T(0.31), T=T)
                    modulation_storage_check_phases(
                        WavefrontSensors.pyramid_operating_modulation(sensor), T;
                        policy=updated_policy, rtol, atol)
                    fresh_policy = (modulation=enabled ? radius : zero(T), modulation_points=4,
                        modulation_phase_offset_rad=T(0.31))
                    # A centered cycle and a single unmodulated point have the
                    # same normalized optical sum.
                    fresh = modulation_storage_sensor(tel, T; geometry, policy=fresh_policy)
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
            end
        end
    end
end



function modulation_storage_field_rate(sensor, pupil, source, ::Type{T}) where {T}
    field = ElectricField(pupil, source; zero_padding=1, T=T)
    formation = prepare_pupil_field(pupil, source, field; center_even_grid=false)
    fill_electric_field!(field, pupil, formation)
    front = PyramidOpticalFrontEnd(sensor)
    rate = pyramid_rate_map(front, field)
    prepared = prepare_wfs_optics(front, field, rate)
    form_wfs_optical_products!(rate, field, prepared)
    return Array(rate.values)
end

@testset "Pyramid separable ElectricField and uncentered pupil parity" begin
    for T in (Float32, Float64), centered in (true, false)
        @testset "$T psf_centering=$centered" begin
            cpu_tel, cpu_pupil, source = modulation_storage_fixture(T, CPUBackend())
            gpu_tel, gpu_pupil, _ = modulation_storage_fixture(T, CUDABackend())
            geometry = (psf_centering=centered,)
            policy = (modulation=T(0.7), modulation_points=4,
                modulation_phase_offset_rad=T(0.31))
            cpu = modulation_storage_sensor(cpu_tel, T; geometry, policy)
            gpu = modulation_storage_sensor(gpu_tel, T; geometry, policy)
            batch = WavefrontSensors.pyramid_propagation_workspace(gpu).modulation_batch
            if centered
                @test batch isa WavefrontSensors.PyramidModulationBatchWorkspace
            else
                @test batch isa WavefrontSensors.NoPyramidModulationBatchWorkspace
            end
            @test WavefrontSensors.pyramid_operating_modulation(gpu) isa
                ModulationOptics.PreparedSeparableFocalPlaneModulation
            rtol = T === Float32 ? T(3e-5) : T(2e-12)
            atol = T === Float32 ? T(1e-6) : T(1e-13)
            @test modulation_storage_diagnostic(gpu, gpu_pupil, source) ≈
                modulation_storage_diagnostic(cpu, cpu_pupil, source) rtol=rtol atol=atol
            cpu_pupil_rate = modulation_storage_rate(cpu, cpu_pupil, source)
            gpu_pupil_rate = modulation_storage_rate(gpu, gpu_pupil, source)
            @test gpu_pupil_rate ≈ cpu_pupil_rate rtol=rtol atol=atol
            cpu_field_rate = modulation_storage_field_rate(cpu, cpu_pupil, source, T)
            gpu_field_rate = modulation_storage_field_rate(gpu, gpu_pupil, source, T)
            @test gpu_field_rate ≈ cpu_field_rate rtol=rtol atol=atol
            @test cpu_field_rate ≈ cpu_pupil_rate rtol=rtol atol=atol
            @test gpu_field_rate ≈ gpu_pupil_rate rtol=rtol atol=atol
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
    initial_opd = Array(pupil.opd)
    geometries = ((;), (rooftop=T(0.25), phase_mask_rotation_rad=T(0.1)),
            (phase_mask_rotation_rad=T(0.17),
                pupil_shift_x_pixels=(T(0.15),T(-0.1),T(0.2),T(-0.25)),
                pupil_shift_y_pixels=(T(-0.2),T(0.05),T(0.1),T(0.3))))
    for geometry in geometries
        @testset "geometry=$geometry" begin
            copyto!(pupil.opd, initial_opd)
            policy = (modulation=T(0.7), modulation_points=4,
                modulation_phase_offset_rad=T(0.31))
            sensor = modulation_storage_sensor(telescope, T; geometry, policy)
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
                # must now see the updated modulation factors without recapture.
                modulation_storage_replay!(captured, context)
                actual = Array(rate.values)
                fresh_policy = (modulation=enabled ? radius : zero(T),
                    modulation_points=4, modulation_phase_offset_rad=T(0.31))
                fresh = modulation_storage_sensor(telescope, T; geometry, policy=fresh_policy)
                expected = modulation_storage_rate(fresh, pupil, source)
                @test actual ≈ expected rtol=3f-5 atol=1f-6
                @test !isapprox(actual, initial_rate; rtol=1f-4, atol=1f-6)
            end
            # The same graph must recompute the phase cache from current OPD.
            # Poisoning the cache makes accidental reuse observable directly.
            previous_rate = Array(rate.values)
            changed_opd = initial_opd .+ [T(110e-9) *
                (T(i - 1) / T(15) - T(0.4) * T(j - 1) / T(15))
                for i in 1:16, j in 1:16]
            copyto!(pupil.opd, changed_opd)
            batch = WavefrontSensors.pyramid_propagation_workspace(sensor).modulation_batch
            fill!(batch.pupil_phase, Complex{T}(T(NaN), T(NaN)))
            CUDA.synchronize()
            modulation_storage_replay!(captured, context)
            refreshed_rate = Array(rate.values)
            @test all(isfinite, refreshed_rate)
            @test Array(batch.pupil_phase) ≈
                cispi.((T(2) / T(0.75e-6)) .* changed_opd) rtol=3f-5 atol=1f-6
            @test map(pointer, modulation_storage_arrays(sensor)) == addresses
            @test pointer(rate.values) == output_address
            cpu_tel, cpu_pupil, _ = modulation_storage_fixture(T, CPUBackend())
            copyto!(cpu_pupil.opd, changed_opd)
            cpu_policy = (modulation=T(0.3), modulation_points=4,
                modulation_phase_offset_rad=T(0.31))
            cpu_sensor = modulation_storage_sensor(cpu_tel, T; geometry, policy=cpu_policy)
            @test refreshed_rate ≈
                modulation_storage_rate(cpu_sensor, cpu_pupil, source) rtol=3f-5 atol=1f-6
            @test !isapprox(refreshed_rate, previous_rate; rtol=1f-4, atol=1f-6)
            for _ in 1:16
                modulation_storage_replay!(captured, context)
            end
            modulation_storage_replay_bytes(captured, context)
            @test modulation_storage_replay_bytes(captured, context) == 0
        end
    end
end
