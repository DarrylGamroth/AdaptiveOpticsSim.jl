using Test
import AdaptiveOpticsSim

const InfluenceOptics = AdaptiveOpticsSim.Optics
const InfluenceInvalid = AdaptiveOpticsSim.InvalidConfiguration
const InfluenceDimensions = AdaptiveOpticsSim.DimensionMismatchError

function influence_sampling_options()
    return (; grid_shape=(3, 3), sample_origin=(1.0, 1.0), sample_pitch=(1.0, 1.0),
        pupil_resolution=2, actuator_pitch_pupil_radius=1.0, command_indices=(2, 1),
        opd_scale_m=2e-6, reference_command_m=1e-6)
end

@testset "measured influence bilinear sampling" begin
    options = influence_sampling_options()
    for T in (Float32, Float64)
        samples = [T(100plane + 10x + y + x*y) for x in 1:3, y in 1:3, plane in 1:2]
        # x increases along the first array axis; source y decreases along
        # the second. The xy term also checks true bilinear interpolation.
        first_plane = Float32[121.25, 133.75, 118.75, 130.25]
        expected = 2.0f0 .* hcat(first_plane .+ 100.0f0, first_plane)
        modes = InfluenceOptics.resample_measured_influences(samples; options...)
        @test typeof(modes) === Matrix{Float32}
        @test modes == expected
        @test InfluenceOptics.resample_measured_influences(
            view(samples, :, :, :); options...) == modes
        backing = zeros(T, 6, 3, 2)
        backing[1:2:6, :, :] .= samples
        @test InfluenceOptics.resample_measured_influences(
            view(backing, 1:2:6, :, :); options...) == modes
        @test InfluenceOptics.resample_measured_influences(samples;
            merge(options, (; reference_command_m=-1e-6))...) == -modes
        @test InfluenceOptics.resample_measured_influences(samples;
            merge(options, (; command_indices=[1, 2]))...) == modes[:, [2, 1]]
        @test InfluenceOptics.resample_measured_influences(samples;
            merge(options, (; pupil_resolution=1))...) == reshape(Float32[452, 252], 1, 2)
        # Pupil samples on the source-grid endpoints remain admissible.
        @test size(InfluenceOptics.resample_measured_influences(samples;
            merge(options, (; actuator_pitch_pupil_radius=0.5))...)) == (4, 2)
    end
    precise = fill(1.0 + 2.0^-24, 3, 3, 2)
    precision_options = merge(options, (; opd_scale_m=3e-6))
    @test all(==(3.0f0), InfluenceOptics.resample_measured_influences(
        Float32.(precise); precision_options...))
    @test all(==(nextfloat(3.0f0)), InfluenceOptics.resample_measured_influences(
        precise; precision_options...))
end

@testset "measured influence input admission" begin
    options = influence_sampling_options()
    samples = ones(Float64, 3, 3, 2)
    for change in (
        (; grid_shape=(0, 3)), (; grid_shape=(true, 3)),
        (; sample_origin=(NaN, 1)), (; sample_pitch=(1, 0)),
        (; sample_pitch=(1, Inf)), (; pupil_resolution=0), (; pupil_resolution=true),
        (; pupil_resolution=typemax(Int)), (; actuator_pitch_pupil_radius=0),
        (; actuator_pitch_pupil_radius=NaN), (; actuator_pitch_pupil_radius=0.1),
        (; command_indices=(1, 1)), (; command_indices=(1, 3)),
        (; command_indices=(true, 2)), (; opd_scale_m=0), (; opd_scale_m=Inf),
        (; reference_command_m=0), (; reference_command_m=NaN),
        (; opd_scale_m=floatmax(Float64), reference_command_m=floatmin(Float64)),
    )
        @test_throws InfluenceInvalid InfluenceOptics.resample_measured_influences(
            samples; merge(options, change)...)
    end
    for change in ((; grid_shape=(3,)), (; sample_origin=(1,)),
        (; sample_pitch=(1,)), (; command_indices=(1,)))
        @test_throws InfluenceDimensions InfluenceOptics.resample_measured_influences(
            samples; merge(options, change)...)
    end
    for invalid in (ones(3, 3), ones(Int, 3, 3, 2))
        @test_throws InfluenceInvalid InfluenceOptics.resample_measured_influences(
            invalid; options...)
    end
    @test_throws InfluenceDimensions InfluenceOptics.resample_measured_influences(
        ones(1, 3, 2); options...)
    for value in (NaN, Inf, floatmax(Float64))
        @test_throws InfluenceInvalid InfluenceOptics.resample_measured_influences(
            fill(value, 3, 3, 2); options...)
    end
end
