include(joinpath(@__DIR__, "pyrtc_reference_hil.jl"))

for wavefront_sensor in (:shack_hartmann, :pyramid)
    PyRTCReferenceHIL.main(wavefront_sensor, ARGS)
end
