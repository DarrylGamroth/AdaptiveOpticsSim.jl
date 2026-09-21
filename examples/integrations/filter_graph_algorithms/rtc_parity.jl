using Test

using FilterGraphAlgorithms

import FilterGraphAlgorithms: commit!, process!

const RtcRouteResult = Union{
    Nothing,
    ControllerToVdmFailure,
    VdmToPdmFailure,
    PdmCommandFailure,
    PdmFeedbackToVdmFailure,
    VdmFeedbackToControllerFailure,
}

# These numerical fixtures preserve the legacy AOS RTC oracle from
# d6a683008d9e40ca5bd90e5ccf748e5aa3384caa. They exercise the replacement FGA
# plans directly; no AOS RTC implementation participates at run time.

function prepare_rtc_projection_fixture()
    controller_to_vdm = Float32[
        1 0
        0 1
        1 -1
    ]
    active_to_full = Float32[
        1 0 0
        0 1 0
        0 0 1
        0.5 0.5 0
    ]
    vdm_to_pdm = Float32[
        1 0 0 0
        0 1 0 0
        0 0 1 0
        0 0 0 1
        1 0 0 1
    ]
    pdm_to_vdm = Float32[
        1 0 0 0 0
        0 1 0 0 0
        0 0 1 0 0
        0 0 0 1 0
    ]
    full_to_active = Float32[
        1 0 0 0
        0 1 0 0
        0 0 1 0
    ]
    vdm_to_controller = Float32[
        1 0 0
        0 0.5 0.5
    ]

    controller = prepare_algorithm(
        ControllerToVdmF32,
        ControllerToVdmF32Configuration(
            controller_command_extent=2,
            vdm_size=3,
        ),
    )
    replace_parameter!(controller, "controller-to-vdm", controller_to_vdm)
    vdm_to_pdm_prepared = prepare_algorithm(
        VdmToPdmF32,
        VdmToPdmF32Configuration(
            active_vdm_size=3,
            full_vdm_size=4,
            physical_actuator_extent=5,
        ),
    )
    replace_parameter!(vdm_to_pdm_prepared, "active-to-full", active_to_full)
    replace_parameter!(vdm_to_pdm_prepared, "vdm-to-pdm", vdm_to_pdm)
    pdm_command = prepare_algorithm(
        PdmCommandF32,
        PdmCommandF32Configuration(
            actuator_count=5,
            lower_limit=Float32[-1, -2, -3, 0, -4],
            upper_limit=Float32[1, 2, 2, 1, 2],
            slew_limit=nothing,
            slew_mode="independent",
            initial_current=zeros(Float32, 5),
            quantization_origin=nothing,
            quantization_step=nothing,
            dither_seed=nothing,
        ),
    )
    pdm_to_vdm_prepared = prepare_algorithm(
        PdmFeedbackToVdmF32,
        PdmFeedbackToVdmF32Configuration(
            active_vdm_size=3,
            full_vdm_size=4,
            physical_actuator_extent=5,
        ),
    )
    replace_parameter!(pdm_to_vdm_prepared, "full-to-active", full_to_active)
    replace_parameter!(pdm_to_vdm_prepared, "pdm-to-vdm", pdm_to_vdm)
    vdm_to_controller_prepared = prepare_algorithm(
        VdmFeedbackToControllerF32,
        VdmFeedbackToControllerF32Configuration(
            controller_command_extent=2,
            vdm_size=3,
        ),
    )
    replace_parameter!(
        vdm_to_controller_prepared,
        "vdm-to-controller",
        vdm_to_controller,
    )
    return (;
        controller,
        vdm_to_pdm_prepared,
        pdm_command,
        pdm_to_vdm_prepared,
        vdm_to_controller_prepared,
        controller_command=Float32[2, -1],
        vdm_command=zeros(Float32, 3),
        requested=zeros(Float32, 5),
        demanded=zeros(Float32, 5),
        pdm_feedback=zeros(Float32, 5),
        vdm_feedback=zeros(Float32, 3),
        controller_feedback=zeros(Float32, 2),
    )
end

function route_rtc_command!(fixture)
    result = process!(
        fixture.vdm_command,
        fixture.controller,
        fixture.controller_command,
    )
    result === nothing || return result
    result = process!(
        fixture.requested,
        fixture.vdm_to_pdm_prepared,
        fixture.vdm_command,
    )
    result === nothing || return result
    result = process!(
        fixture.demanded,
        fixture.pdm_feedback,
        fixture.pdm_command,
        fixture.requested,
    )
    result === nothing || return result
    result = process!(
        fixture.vdm_feedback,
        fixture.pdm_to_vdm_prepared,
        fixture.pdm_feedback,
    )
    result === nothing || return result
    result = process!(
        fixture.controller_feedback,
        fixture.vdm_to_controller_prepared,
        fixture.vdm_feedback,
    )
    result === nothing || return result
    return nothing
end

function prepare_closed_loop_fixture()
    plan = ClosedLoopCorrectionPlan(
        -0.3f0,
        0.99f0,
        1.0f0,
        0.0f0,
        2,
        zeros(Float32, 1, 2),
        zeros(Float32, 2, 1),
    )
    return (;
        plan,
        workspace=ClosedLoopCorrectionWorkspace(plan),
        correction=zeros(Float32, 2),
        controller_state=zeros(Float32, 2),
        residual_error=Float32[1, 2],
        constraint_feedback=Float32[0.1, -0.2],
    )
end

function apply_closed_loop_fixture!(fixture, feedback)
    result = process!(
        fixture.correction,
        fixture.controller_state,
        fixture.workspace,
        fixture.plan,
        fixture.residual_error,
        feedback,
    )
    result === nothing || return result
    commit!(fixture.workspace)
    return nothing
end

function prepare_fixed_delay_fixture()
    prepared = prepare_algorithm(
        FrameDelayF32,
        FrameDelayF32Configuration(
            extent=2,
            delay_frames=2,
            initial_value=0.0f0,
            input_schema=CONTROLLER_COMMAND_V1,
            output_schema=CONTROLLER_COMMAND_V1,
        ),
    )
    return (; prepared, input=zeros(Float32, 2), output=zeros(Float32, 2))
end

function delay_rtc_command!(fixture)
    return process!(fixture.output, fixture.prepared, fixture.input)
end

@testset "FGA RTC projection parity" begin
    fixture = prepare_rtc_projection_fixture()
    @test @inferred(RtcRouteResult, route_rtc_command!(fixture)) === nothing
    @test fixture.vdm_command == Float32[2, -1, 3]
    @test fixture.vdm_to_pdm_prepared.workspace.full_vdm ==
          Float32[2, -1, 3, 0.5]
    @test fixture.requested == Float32[2, -1, 3, 0.5, 2.5]
    @test fixture.demanded == Float32[1, -1, 2, 0.5, 2]
    @test fixture.pdm_feedback == Float32[1, 0, 1, 0, 0.5]
    @test fixture.vdm_feedback == Float32[1, 0, 1]
    @test fixture.controller_feedback == Float32[1, 0.5]

    route_rtc_command!(fixture)
    @test @allocated(route_rtc_command!(fixture)) == 0
end

@testset "FGA closed-loop correction parity" begin
    fixture = prepare_closed_loop_fixture()
    @test @inferred(Union{Nothing,ClosedLoopCorrectionFailure},
        apply_closed_loop_fixture!(fixture, nothing)) === nothing
    @test fixture.correction ≈ Float32[-0.3, -0.6]
    @test fixture.controller_state == Float32[0, 0]

    fixture.residual_error .= Float32[0.5, -1]
    @test @inferred(Union{Nothing,ClosedLoopCorrectionFailure},
        apply_closed_loop_fixture!(fixture, fixture.constraint_feedback)) ===
          nothing
    @test fixture.controller_state ≈ Float32[-0.4, -0.4]
    @test fixture.correction ≈ Float32[-0.546, -0.096]

    apply_closed_loop_fixture!(fixture, fixture.constraint_feedback)
    @test @allocated(apply_closed_loop_fixture!(fixture, fixture.constraint_feedback)) == 0
end

@testset "FGA fixed-frame delay parity" begin
    fixture = prepare_fixed_delay_fixture()
    expected = (
        (input=Float32[1, 10], output=Float32[0, 0]),
        (input=Float32[2, 20], output=Float32[0, 0]),
        (input=Float32[3, 30], output=Float32[1, 10]),
        (input=Float32[4, 40], output=Float32[2, 20]),
    )
    for sample in expected
        copyto!(fixture.input, sample.input)
        @test @inferred(Union{Nothing,FrameDelayFailure},
            delay_rtc_command!(fixture)) === nothing
        @test fixture.output == sample.output
    end

    copyto!(fixture.input, Float32[5, 50])
    delay_rtc_command!(fixture)
    @test @allocated(delay_rtc_command!(fixture)) == 0
end
