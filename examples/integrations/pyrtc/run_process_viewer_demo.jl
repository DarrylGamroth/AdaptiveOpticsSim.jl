include(joinpath(@__DIR__, "pyrtc_process_hil.jl"))

system = isempty(ARGS) ? :pyramid : Symbol(ARGS[1])
duration = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 60.0
frame_rate = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : 10.0

PyRTCProcessHIL.viewer_main(
    system;
    duration,
    frame_rate,
)
