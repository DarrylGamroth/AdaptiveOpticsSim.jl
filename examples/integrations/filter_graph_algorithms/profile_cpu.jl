using Profile

include("s1_lockstep.jl")
using .AOSFGALockstep

const PROFILE_STEPS = 500_000

function run_steps!(prepared, steps)
    for _ in 1:steps
        step_lockstep!(prepared)
    end
    return prepared
end

function collect_allocation_profile(prepared, steps)
    Profile.Allocs.clear()
    Profile.Allocs.start(sample_rate=1.0)
    try
        run_steps!(prepared, steps)
    finally
        Profile.Allocs.stop()
    end
    return Profile.Allocs.fetch()
end

prepared = prepare_s1_lockstep()
for _ in 1:32
    step_lockstep!(prepared)
end

allocation_bytes = @allocated step_lockstep!(prepared)
allocation_bytes == 0 || error(
    "warmed S1 lockstep allocated $allocation_bytes Julia heap bytes",
)
all(isfinite, prepared.outputs.demanded) || error(
    "pre-profile demanded command is not finite",
)

Profile.clear()
Profile.init(n=10^7, delay=0.001)
@profile run_steps!(prepared, PROFILE_STEPS)

println("S1_CPU_PROFILE_BEGIN")
Profile.print(format=:flat, sortedby=:count, mincount=10)
println("S1_CPU_PROFILE_END")

allocation_profile = collect_allocation_profile(prepared, 100)
println("S1_CPU_ALLOCATION_PROFILE_BEGIN")
Profile.Allocs.print(stdout, allocation_profile)
println("S1_CPU_ALLOCATION_PROFILE_END")

all(isfinite, prepared.outputs.demanded) || error(
    "post-profile demanded command is not finite",
)
abs(prepared.outputs.demanded[1] + 3.0f-8) <= 1.0f-10 || error(
    "post-profile demanded command did not retain the closed-loop oracle",
)

println("profile_steps=$PROFILE_STEPS")
println("final_sequence=$(prepared.state.sequence)")
println("final_demanded=$(prepared.outputs.demanded[1])")
println("warmed_allocation_bytes=$allocation_bytes")
