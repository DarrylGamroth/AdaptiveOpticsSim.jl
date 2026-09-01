include(joinpath(@__DIR__, "test_shared_memory.jl"))

if get(ENV, "AOS_PYRTC_PROCESS_TESTS", "0") == "1"
    include(joinpath(@__DIR__, "test_process_hil.jl"))
end
