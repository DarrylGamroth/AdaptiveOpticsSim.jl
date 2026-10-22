# Standalone focused checks; set AOS_TEST_SEPARABLE_MODULATION_CPU_ONLY=1 to
# validate the scalar and KernelAbstractions CPU paths without loading CUDA.
using Test
using AdaptiveOpticsSim
using AdaptiveOpticsSim.Backends
using FixedSizeArrays: FixedSizeVector
using KernelAbstractions

const SeparableModulationOptics = AdaptiveOpticsSim.Optics

function separable_modulation_update_bytes(prepared, radius)
    return @allocated SeparableModulationOptics.update_cycle_averaged_circular_modulation!(
        prepared, radius)
end

function check_separable_modulation(prepared, policy, ::Type{T}) where {T}
    O = SeparableModulationOptics
    first_factors = Array(prepared.axis_1_factors)
    second_factors = Array(prepared.axis_2_factors)
    resolution, points = size(first_factors)
    coordinates = collect(range(-T(pi), T(pi); length=resolution))
    dense = O.prepare_focal_plane_modulation(policy, resolution, zeros(T, 1), T)
    @test size(second_factors) == (resolution, points)
    @test O.modulation_point_count(prepared) == points
    @test prepared.amplitude_weights isa FixedSizeVector{T}
    @test prepared.amplitude_weights == dense.amplitude_weights
    @test sum(abs2, prepared.amplitude_weights) ≈ one(T)
    @test O._modulation_phase_storage(prepared) ===
        (prepared.axis_1_factors, prepared.axis_2_factors)
    @test O._modulation_phase_storage(dense) === dense.phases
    reconstructed = [O._modulation_phase((first_factors, second_factors), i, j, point)
        for i in 1:resolution, j in 1:resolution, point in 1:points]
    @test reconstructed ≈ dense.phases rtol=64eps(T) atol=16eps(T)
    @test O._modulation_phase(dense.phases, 1, 1, 1) == dense.phases[1, 1, 1]
    for point in 1:points
        offset_x, offset_y = O.modulation_offset(policy, point, T)
        @test first_factors[:, point] ≈ cis.(offset_y .* coordinates) rtol=16eps(T) atol=8eps(T)
        @test second_factors[:, point] ≈ cis.(offset_x .* coordinates) rtol=16eps(T) atol=8eps(T)
    end
end

function test_separable_modulation_storage(storage, ::Type{T}; host=false) where {T}
    O = SeparableModulationOptics
    policies = (
        O.NoModulation(),
        O.CircularModulation(zero(T); samples=5, T=T),
        O.CircularModulation(T(0.7); samples=7, phase_offset_rad=T(0.31), T=T),
        O.SampledModulation(((T(1.3), T(0.7)), (T(-0.4), T(1.1)),
            (T(-3.4), T(-2.8))); weights=T[0, 2, 1], T=T),
    )
    for resolution in (2, 16, 17), policy in policies
        prepared = O.prepare_separable_focal_plane_modulation(policy, resolution, storage, T)
        @test typeof(prepared.axis_1_factors) === typeof(similar(storage, Complex{T}, resolution, O.modulation_point_count(policy)))
        @test length(prepared.axis_1_factors) + length(prepared.axis_2_factors) ==
            2resolution * O.modulation_point_count(policy)
        check_separable_modulation(prepared, policy, T)
    end
    @test_throws InvalidConfiguration O.prepare_separable_focal_plane_modulation(
        O.NoModulation(), 0, storage, T)

    policy = O.CircularModulation(T(0.7); samples=7, phase_offset_rad=T(0.31), T=T)
    prepared = O.prepare_separable_focal_plane_modulation(policy, 17, storage, T)
    arrays = (prepared.axis_1_factors, prepared.axis_2_factors, prepared.amplitude_weights)
    addresses = map(pointer, arrays)
    weights = copy(prepared.amplitude_weights)
    for (radius, enabled) in ((T(1.1), true), (zero(T), true), (T(0.5), false), (T(0.3), true))
        @test O.update_cycle_averaged_circular_modulation!(prepared, radius; enabled) === prepared
        @test map(pointer, arrays) == addresses
        @test O._modulation_phase_storage(prepared) === arrays[1:2]
        @test prepared.amplitude_weights == weights
        current_policy = O.CircularModulation(enabled ? radius : zero(T);
            samples=7, phase_offset_rad=T(0.31), T=T)
        check_separable_modulation(prepared, current_policy, T)
    end
    for invalid in (T(-1), T(NaN), T(Inf), T === Float32 ? Float64(1) : Float32(1))
        before = map(Array, arrays)
        @test_throws InvalidConfiguration O.update_cycle_averaged_circular_modulation!(prepared, invalid)
        @test isequal(map(Array, arrays), before)
        @test map(pointer, arrays) == addresses
    end
    for malformed in (
        O.PreparedSeparableFocalPlaneModulation(policy, prepared.axis_1_factors,
            similar(storage, Complex{T}, 16, 7), prepared.amplitude_weights),
        O.PreparedSeparableFocalPlaneModulation(policy,
            similar(storage, Complex{T}, 17, 6), similar(storage, Complex{T}, 17, 6),
            prepared.amplitude_weights),
    )
        @test_throws InvalidConfiguration O.update_cycle_averaged_circular_modulation!(malformed, T(0.7))
    end
    if host
        separable_modulation_update_bytes(prepared, T(0.6))
        @test separable_modulation_update_bytes(prepared, T(0.6)) == 0
    end
end

@testset "Separable physical modulation factors" begin
    for T in (Float32, Float64)
        @testset "$T scalar host" begin
            test_separable_modulation_storage(zeros(T, 1), T; host=true)
        end
        @testset "$T accelerator kernel on CPU" begin
            O = SeparableModulationOptics
            style = Backends.AcceleratorStyle(KernelAbstractions.CPU())
            policy = O.SampledModulation(((T(1.3), T(0.7)), (T(-0.4), T(1.1)));
                weights=T[1, 2], T=T)
            prepared = O.prepare_separable_focal_plane_modulation(policy, 17, zeros(T, 1), T)
            O._prepare_separable_modulation_factors!(style, prepared.axis_1_factors,
                prepared.axis_2_factors, policy, T)
            check_separable_modulation(prepared, policy, T)
            circular = O.CircularModulation(T(0.7); samples=7, phase_offset_rad=T(0.31), T=T)
            prepared_circular = O.prepare_separable_focal_plane_modulation(circular, 17, zeros(T, 1), T)
            O._update_cycle_averaged_circular_modulation!(style, prepared_circular, T(1.1), T(0.31))
            updated = O.CircularModulation(T(1.1); samples=7, phase_offset_rad=T(0.31), T=T)
            check_separable_modulation(prepared_circular, updated, T)
        end
    end
end

if get(ENV, "AOS_TEST_SEPARABLE_MODULATION_CPU_ONLY", "0") != "1"
    @eval import CUDA
    CUDA.allowscalar(false)
    @testset "CUDA separable physical modulation factors" begin
        @test CUDA.functional()
        for T in (Float32, Float64)
            @testset "$T CUDA" begin
                test_separable_modulation_storage(CUDA.zeros(T, 1), T)
                O = SeparableModulationOptics
                policy = O.NoModulation()
                prepared = O.prepare_separable_focal_plane_modulation(policy, 16, CUDA.zeros(T, 1), T)
                target = Backends.compute_device(prepared.axis_1_factors)
                require_target = AdaptiveOpticsSim.WavefrontSensors._require_exact_focal_plane_modulation_target
                @test require_target(prepared, target, "test modulation") === prepared
                for mixed in (
                    O.PreparedSeparableFocalPlaneModulation(policy, Array(prepared.axis_1_factors),
                        prepared.axis_2_factors, prepared.amplitude_weights),
                    O.PreparedSeparableFocalPlaneModulation(policy, prepared.axis_1_factors,
                        Array(prepared.axis_2_factors), prepared.amplitude_weights),
                )
                    @test_throws AdaptiveOpticsSim.WavefrontSensors.WFSPreparationError require_target(mixed, target, "test modulation")
                end
            end
        end
    end
end
