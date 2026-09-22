using Profile

include("profile_lift_aoc_common.jl")

function cpu_profile_region(f)
    Profile.clear()
    Profile.@profile f()
    println("cpu_profile_samples = ", length(Profile.fetch()))
    Profile.print(format=:flat, sortedby=:count)
    return nothing
end

profile_lift_adapter("cpu", Array, CPUBackend(), () -> nothing,
    f -> (f(); 0), () -> 0, cpu_profile_region)
