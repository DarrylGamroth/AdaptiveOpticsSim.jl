function captured_plant_time_error(f)
    try
        f()
    catch error
        return error
    end
    return nothing
end

function assert_plant_time_error(f, operation::Symbol, reason::Symbol)
    error = captured_plant_time_error(f)
    @test error isa PlantTimeError
    if error isa PlantTimeError
        @test error.operation === operation
        @test error.reason === reason
        @test !isempty(error.msg)
    end
    return error
end

function assert_plant_schedule_error(f, component::Symbol, reason::Symbol)
    error = captured_plant_time_error(f)
    @test error isa PlantScheduleError
    if error isa PlantScheduleError
        @test error.component === component
        @test error.reason === reason
        @test !isempty(error.msg)
    end
    return error
end

@inline function plant_time_hot_cycle(timestamp, duration, schedule)
    advanced = timestamp + duration
    elapsed = advanced - timestamp
    due = schedule_timestamp(schedule, 3, timestamp)
    return plant_nanoseconds(advanced) + plant_nanoseconds(elapsed) +
        plant_nanoseconds(due)
end

@testset "Canonical plant-time values" begin
    @test isbitstype(PlantTimestamp)
    @test isbitstype(PlantDuration)
    @test isbitstype(AdaptiveOpticsSim.PlantTimeOffset)
    @test sizeof(PlantTimestamp) == sizeof(Int64)
    @test sizeof(PlantDuration) == sizeof(Int64)
    @test sizeof(AdaptiveOpticsSim.PlantTimeOffset) == sizeof(Int64)

    timestamp = PlantTimestamp(1_250_000_000)
    duration = PlantDuration(UInt32(250_000_000))
    @test plant_nanoseconds(timestamp) == 1_250_000_000
    @test plant_nanoseconds(duration) == 250_000_000
    @test zero(PlantTimestamp) == PlantTimestamp(0)
    @test zero(timestamp) == PlantTimestamp(0)
    @test zero(PlantDuration) == PlantDuration(0)
    @test zero(duration) == PlantDuration(0)
    @test iszero(zero(timestamp))
    @test iszero(zero(duration))
    @test sprint(show, timestamp) == "PlantTimestamp(1250000000 ns)"
    @test sprint(show, duration) == "PlantDuration(250000000 ns)"
    offset = AdaptiveOpticsSim.PlantTimeOffset(-250_000_000)
    @test plant_nanoseconds(offset) == -250_000_000
    @test sprint(show, offset) == "PlantTimeOffset(-250000000 ns)"
    @test timestamp + offset == PlantTimestamp(1_000_000_000)
    @test offset + timestamp == PlantTimestamp(1_000_000_000)
    @test offset + AdaptiveOpticsSim.PlantTimeOffset(50_000_000) ==
        AdaptiveOpticsSim.PlantTimeOffset(-200_000_000)
    @test offset - AdaptiveOpticsSim.PlantTimeOffset(-50_000_000) ==
        AdaptiveOpticsSim.PlantTimeOffset(-200_000_000)
    @test zero(offset) == AdaptiveOpticsSim.PlantTimeOffset(0)
    @test_throws PlantTimeError PlantTimestamp(0) +
        AdaptiveOpticsSim.PlantTimeOffset(-1)
    @test_throws PlantTimeError PlantTimestamp(typemax(Int64)) +
        AdaptiveOpticsSim.PlantTimeOffset(1)
    @test_throws PlantTimeError AdaptiveOpticsSim.PlantTimeOffset(
        big(typemax(Int64)) + 1)
    @test_throws PlantTimeError AdaptiveOpticsSim.PlantTimeOffset(true)

    @test PlantTimestamp(1) < PlantTimestamp(2)
    @test PlantTimestamp(1) <= PlantTimestamp(1)
    @test PlantDuration(1) < PlantDuration(2)
    @test PlantDuration(1) <= PlantDuration(1)
    @test length(Set((PlantTimestamp(1), PlantTimestamp(1)))) == 1
    @test length(Set((PlantDuration(1), PlantDuration(1)))) == 1
    @test length(Set{Any}((PlantTimestamp(1), PlantDuration(1)))) == 2

    @test timestamp + duration == PlantTimestamp(1_500_000_000)
    @test duration + timestamp == PlantTimestamp(1_500_000_000)
    @test timestamp - duration == PlantTimestamp(1_000_000_000)
    @test timestamp - PlantTimestamp(250_000_000) ==
        PlantDuration(1_000_000_000)
    @test duration + PlantDuration(50_000_000) == PlantDuration(300_000_000)
    @test duration - PlantDuration(50_000_000) == PlantDuration(200_000_000)
    @test duration * 3 == PlantDuration(750_000_000)
    @test 3 * duration == PlantDuration(750_000_000)

    @test plant_time_seconds(timestamp) == 1.25
    @test plant_duration_seconds(duration) == 0.25
    @test plant_time_seconds(timestamp, Float32) === Float32(1.25)
    @test plant_duration_seconds(duration, Float32) === Float32(0.25)

    @test_throws MethodError timestamp + PlantTimestamp(1)
    @test_throws MethodError duration + 1
    @test_throws MethodError plant_time_seconds(duration)
    @test_throws MethodError plant_duration_seconds(timestamp)

    assert_plant_time_error(() -> PlantTimestamp(-1), :plant_timestamp,
        :negative)
    assert_plant_time_error(() -> PlantDuration(-1), :plant_duration,
        :negative)
    assert_plant_time_error(() -> PlantTimestamp(true), :plant_timestamp,
        :invalid_type)
    assert_plant_time_error(() -> PlantDuration(false), :plant_duration,
        :invalid_type)
    assert_plant_time_error(
        () -> PlantTimestamp(big(typemax(Int64)) + 1),
        :plant_timestamp,
        :overflow,
    )
    assert_plant_time_error(
        () -> PlantDuration(big(typemax(Int64)) + 1),
        :plant_duration,
        :overflow,
    )
    assert_plant_time_error(
        () -> PlantTimestamp(typemax(Int64)) + PlantDuration(1),
        :timestamp_add,
        :overflow,
    )
    assert_plant_time_error(
        () -> PlantDuration(typemax(Int64)) + PlantDuration(1),
        :duration_add,
        :overflow,
    )
    assert_plant_time_error(
        () -> PlantTimestamp(0) - PlantDuration(1),
        :timestamp_subtract,
        :negative_result,
    )
    assert_plant_time_error(
        () -> PlantTimestamp(0) - PlantTimestamp(1),
        :timestamp_difference,
        :negative_result,
    )
    assert_plant_time_error(
        () -> PlantDuration(0) - PlantDuration(1),
        :duration_subtract,
        :negative_result,
    )
    assert_plant_time_error(
        () -> PlantDuration(1) * -1,
        :duration_scale,
        :negative_factor,
    )
    assert_plant_time_error(
        () -> PlantDuration(1) * true,
        :duration_scale,
        :invalid_type,
    )
    assert_plant_time_error(
        () -> PlantDuration(typemax(Int64)) * 2,
        :duration_scale,
        :overflow,
    )
end

@testset "Periodic plant schedules" begin
    schedule = PeriodicSchedule(
        PlantDuration(500_000);
        phase=PlantDuration(125_000),
    )
    keyword_schedule = PeriodicSchedule(period_ns=500_000, phase_ns=125_000)
    @test isbitstype(PeriodicSchedule)
    @test !Base.ismutable(schedule)
    @test schedule_period(schedule) == PlantDuration(500_000)
    @test schedule_phase(schedule) == PlantDuration(125_000)
    @test schedule_period(keyword_schedule) == schedule_period(schedule)
    @test schedule_phase(keyword_schedule) == schedule_phase(schedule)
    @test schedule_timestamp(schedule, 1) == PlantTimestamp(125_000)
    @test schedule_timestamp(schedule, 2) == PlantTimestamp(625_000)
    @test schedule_timestamp(schedule, UInt8(3)) == PlantTimestamp(1_125_000)
    @test schedule_timestamp(schedule, 3, PlantTimestamp(2_000_000)) ==
        PlantTimestamp(3_125_000)

    assert_plant_schedule_error(
        () -> PeriodicSchedule(PlantDuration(0)),
        :periodic_schedule,
        :zero_period,
    )
    assert_plant_schedule_error(
        () -> schedule_timestamp(schedule, 0),
        :periodic_schedule,
        :invalid_sequence,
    )
    assert_plant_schedule_error(
        () -> schedule_timestamp(schedule, true),
        :periodic_schedule,
        :invalid_sequence,
    )
    assert_plant_time_error(
        () -> schedule_timestamp(
            PeriodicSchedule(PlantDuration(typemax(Int64))),
            2,
            PlantTimestamp(1),
        ),
        :timestamp_add,
        :overflow,
    )
end

@testset "Deterministic plant-event ordering" begin
    phases = (
        TriggerUpdatePhase,
        CommandApplicationPhase,
        AtmosphereEvolutionPhase,
        IntegrationBoundaryPhase,
        ExposureOpenPhase,
        OpticalSamplePhase,
        ReadoutCompletionPhase,
        AcquisitionReadyPhase,
    )
    @test UInt8.(phases) == Tuple(UInt8(1):UInt8(8))

    base = PlantEventKey(PlantTimestamp(10), TriggerUpdatePhase, 1, 1)
    later_time = PlantEventKey(PlantTimestamp(11), TriggerUpdatePhase, 1, 1)
    later_phase = PlantEventKey(
        PlantTimestamp(10),
        CommandApplicationPhase,
        1,
        1,
    )
    later_ordinal = PlantEventKey(PlantTimestamp(10), TriggerUpdatePhase, 2, 1)
    later_occurrence = PlantEventKey(
        PlantTimestamp(10),
        TriggerUpdatePhase,
        1,
        2,
    )
    @test base < later_time
    @test base < later_phase
    @test base < later_ordinal
    @test base < later_occurrence
    @test sort([later_occurrence, later_time, later_ordinal, later_phase, base]) ==
        [base, later_occurrence, later_ordinal, later_phase, later_time]

    assert_plant_schedule_error(
        () -> PlantEventKey(PlantTimestamp(0), TriggerUpdatePhase, 0, 1),
        :plant_event_key,
        :invalid_ordinal,
    )
    assert_plant_schedule_error(
        () -> PlantEventKey(
            PlantTimestamp(0),
            TriggerUpdatePhase,
            UInt32(0),
            UInt64(1),
        ),
        :plant_event_key,
        :invalid_ordinal,
    )
    assert_plant_schedule_error(
        () -> PlantEventKey(PlantTimestamp(0), TriggerUpdatePhase, true, 1),
        :plant_event_key,
        :invalid_ordinal,
    )
    assert_plant_schedule_error(
        () -> PlantEventKey(PlantTimestamp(0), TriggerUpdatePhase, 1, 0),
        :plant_event_key,
        :invalid_occurrence,
    )
    assert_plant_schedule_error(
        () -> PlantEventKey(
            PlantTimestamp(0),
            TriggerUpdatePhase,
            UInt32(1),
            UInt64(0),
        ),
        :plant_event_key,
        :invalid_occurrence,
    )
    assert_plant_schedule_error(
        () -> PlantEventKey(PlantTimestamp(0), TriggerUpdatePhase, 1, true),
        :plant_event_key,
        :invalid_occurrence,
    )

    @test !Base.isexported(AdaptiveOpticsSim, :PlantEventPhase)
    @test !Base.isexported(AdaptiveOpticsSim, :PlantEventKey)
end

@testset "Plant-time hot operations" begin
    timestamp = PlantTimestamp(1_000)
    duration = PlantDuration(250)
    schedule = PeriodicSchedule(PlantDuration(100); phase=PlantDuration(25))
    @test @inferred(timestamp + duration) == PlantTimestamp(1_250)
    @test @inferred(timestamp - PlantTimestamp(500)) == PlantDuration(500)
    @test @inferred(schedule_timestamp(schedule, 3, timestamp)) ==
        PlantTimestamp(1_225)
    @test @inferred(plant_time_hot_cycle(timestamp, duration, schedule)) == 2_725

    plant_time_hot_cycle(timestamp, duration, schedule)
    if coverage_instrumented()
        @test_skip "allocation gate disabled under coverage instrumentation"
    else
        @test @allocated(plant_time_hot_cycle(timestamp, duration, schedule)) == 0
    end
end
