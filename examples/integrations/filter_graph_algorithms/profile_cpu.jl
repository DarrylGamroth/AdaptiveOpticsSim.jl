using Profile

include("s1_lockstep.jl")
include("s4_pyramid.jl")
using .AOSFGALockstep
using .AOSFGAPyramid

const PROFILE_STEPS = parse(Int,
    get(ENV, "AOS_FGA_PROFILE_STEPS", "500000"))

function run_steps!(prepared, step!, steps)
    for _ in 1:steps
        step!(prepared)
    end
    return prepared
end

function collect_allocation_profile(prepared, step!, steps)
    Profile.Allocs.clear()
    Profile.Allocs.start(sample_rate=1.0)
    try
        run_steps!(prepared, step!, steps)
    finally
        Profile.Allocs.stop()
    end
    return Profile.Allocs.fetch()
end

function profile_workload!(label, prepared, step!, command_values;
    steps::Int=PROFILE_STEPS)
    for _ in 1:32
        step!(prepared)
    end

    allocation_bytes = @allocated step!(prepared)
    allocation_bytes == 0 || error(
        "warmed $label step allocated $allocation_bytes Julia heap bytes",
    )
    all(isfinite, command_values(prepared)) || error(
        "$label command is not finite before profiling",
    )

    Profile.clear()
    Profile.init(n=10^7, delay=0.001)
    @profile run_steps!(prepared, step!, steps)

    println(label, "_CPU_PROFILE_BEGIN")
    Profile.print(format=:flat, sortedby=:count, mincount=10)
    println(label, "_CPU_PROFILE_END")

    allocation_profile = collect_allocation_profile(prepared, step!, 100)
    println(label, "_CPU_ALLOCATION_PROFILE_BEGIN")
    Profile.Allocs.print(stdout, allocation_profile)
    println(label, "_CPU_ALLOCATION_PROFILE_END")

    all(isfinite, command_values(prepared)) || error(
        "$label command is not finite after profiling",
    )
    println(label, "_profile_steps=", steps)
    println(label, "_final_sequence=", prepared.state.sequence)
    println(label, "_warmed_allocation_bytes=", allocation_bytes)
    return prepared
end

s1 = profile_workload!("S1", prepare_s1_lockstep(), step_lockstep!,
    prepared -> prepared.outputs.demanded)
abs(s1.outputs.demanded[1] + 3.0f-8) <= 1.0f-10 || error(
    "post-profile S1 demanded command did not retain the closed-loop oracle",
)

profile_workload!("S4", prepare_s4_pyramid(), step_s4_pyramid!,
    prepared -> prepared.adopted_command)
