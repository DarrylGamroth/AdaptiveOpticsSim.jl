function assert_prepared_four_pupil_lgs_conformance!(
    model::M,
    intensity::AbstractMatrix{T},
    scratch::AbstractMatrix{T},
    fft_buffer::AbstractMatrix{Complex{T}},
    fft_plan,
    ifft_buffer::AbstractMatrix{Complex{T}},
    ifft_plan,
    target::Backends.AbstractComputeDevice,
) where {M<:WavefrontSensors.AbstractPreparedFourPupilLGS,T<:AbstractFloat}
    @test @inferred(
        WavefrontSensors._require_exact_prepared_four_pupil_lgs_target(
            model, target)) === nothing
    result = @inferred WavefrontSensors.apply_prepared_four_pupil_lgs!(
        model, intensity, scratch, fft_buffer, fft_plan, ifft_buffer,
        ifft_plan)
    @test result === intensity
    @test all(isfinite, intensity)
    return intensity
end
