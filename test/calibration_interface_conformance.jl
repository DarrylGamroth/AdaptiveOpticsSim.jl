function assert_calibration_command_plan_conformance(
    plan::P,
    commands::AbstractMatrix{T},
    amplitude::T,
) where {P<:Calibration.AbstractCalibrationCommandPlan,T<:AbstractFloat}
    @test Calibration.calibration_command_count(plan) == size(commands, 2)
    coefficients = similar(commands, size(commands, 1))
    for column in axes(commands, 2)
        fill!(coefficients, T(NaN))
        result = @inferred Calibration.stage_calibration_command!(
            coefficients, plan, column, amplitude)
        @test result === coefficients
        @test coefficients == amplitude .* view(commands, :, column)
    end
    return nothing
end
