module AdaptiveOpticsSim

__precompile__(true)

using AbstractFFTs
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
include("core/types.jl")
include("core/profiles.jl")
include("core/inverse_policies.jl")
include("core/backends.jl")
include("core/reductions.jl")
include("core/random_services.jl")
include("core/config.jl")
include("core/utils.jl")
include("core/kv56.jl")
include("core/workspace.jl")
include("core/parallel.jl")
include("core/telemetry.jl")

include("optics/aperture_masks.jl")
include("optics/telescope.jl")
include("optics/source.jl")
include("optics/spectrum.jl")
include("optics/electric_field.jl")
include("optics/propagation.jl")
include("optics/psf.jl")
include("optics/zernike.jl")
include("optics/misregistration.jl")
include("optics/controllable_optics.jl")
include("optics/deformable_mirror.jl")
include("detectors/detector.jl")
include("optics/asterism.jl")
include("optics/extended_source.jl")
include("optics/opd_map.jl")
include("optics/spatial_filter.jl")
include("atmosphere/source_geometry.jl")
include("optics/propagation_context.jl")
include("atmosphere/kolmogorov.jl")
include("atmosphere/infinite_screen_math.jl")
include("atmosphere/infinite_screen.jl")
include("atmosphere/multilayer.jl")
include("atmosphere/phase_stats.jl")
include("atmosphere/fast_atmosphere.jl")
include("optics/atmospheric_field_propagation.jl")
include("calibration/modal_basis.jl")
include("optics/ncpa.jl")
include("wfs/sensing_modes.jl")
include("wfs/grouped.jl")
include("wfs/calibration.jl")
include("wfs/elongation.jl")
include("wfs/subapertures.jl")
include("wfs/shack_hartmann.jl")
include("wfs/pyramid.jl")
include("wfs/bioedge.jl")
include("wfs/zernike.jl")
include("wfs/curvature.jl")
include("wfs/interface.jl")
include("wfs/lift.jl")
include("calibration/interaction_matrix.jl")
include("control/reconstructors.jl")
include("calibration/control_matrix.jl")
include("calibration/fitting_error.jl")
include("calibration/ao_calibration.jl")
include("simulation/assembly.jl")
include("calibration/gain_sensing_camera.jl")
include("calibration/misregistration_identification.jl")
include("control/controller.jl")
include("control/runtime_outputs.jl")
include("control/runtime.jl")
include("control/control_loop.jl")
include("tomography/parameters.jl")
include("tomography/geometry.jl")
include("tomography/fitting.jl")
include("tomography/reconstructors.jl")

include("exports.jl")

end # module AdaptiveOpticsSim
