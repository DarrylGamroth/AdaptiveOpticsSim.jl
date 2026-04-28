#
# Closed-loop runtime execution and simulation I/O
#
# This file is intentionally kept as a small entry point. Implementation is
# split across the include files below by responsibility.
#

include("runtime/types.jl")
include("runtime/construction.jl")
include("runtime/outputs.jl")
include("runtime/execution.jl")
include("runtime/timing.jl")
