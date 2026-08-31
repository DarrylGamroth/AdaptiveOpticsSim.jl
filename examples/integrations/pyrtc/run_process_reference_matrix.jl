include(joinpath(@__DIR__, "pyrtc_process_hil.jl"))

for wavefront_sensor in (:shack_hartmann, :pyramid)
    PyRTCProcessHIL.main(wavefront_sensor, ARGS)
end
