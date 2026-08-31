module AdaptiveOpticsSim

__precompile__(true)

using AbstractFFTs
import FFTW
using KernelAbstractions
using LinearAlgebra
using Logging
using Random
using SparseArrays
using SpecialFunctions
using Statistics

"""
AdaptiveOpticsSim.jl

Julia adaptive optics simulation toolkit (in development).
"""
const PROJECT_STATUS = :in_development

include("core/errors.jl")
include("core/profiles.jl")
include("core/random_services.jl")
include("backends/backends.jl")
include("core/utils.jl")

# Configuration serialization remains a root-owned cross-domain seam. Declare
# its generic before domain modules add methods, then install the shared
# fallbacks after the first physical owner has loaded.
function config_value end

include("optics/optics.jl")
include("core/types.jl")
include("core/config.jl")
include("core/kv56.jl")
include("core/workspace.jl")
include("core/telemetry.jl")

include("ensembles/ensembles.jl")
include("atmosphere/atmospheres.jl")
include("detectors/detectors.jl")
include("wfs/wavefront_sensors.jl")
include("calibration/calibration.jl")
include("control/control.jl")
include("tomography/tomography.jl")
include("algorithm_graphs/algorithm_graphs.jl")

include("exports.jl")

end # module AdaptiveOpticsSim
