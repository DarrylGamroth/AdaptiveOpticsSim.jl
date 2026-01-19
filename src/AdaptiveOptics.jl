module AdaptiveOptics

__precompile__(true)

using FFTW
using LinearAlgebra
using Logging
using Random
using Statistics

"""
AdaptiveOptics.jl

Julia adaptive optics simulation toolkit (in development).
"""
const PROJECT_STATUS = :in_development

include("Core/errors.jl")
include("Core/types.jl")
include("Core/utils.jl")
include("Core/workspace.jl")
include("Core/parallel.jl")

include("Optics/telescope.jl")
include("Optics/source.jl")
include("Optics/psf.jl")
include("Optics/zernike.jl")
include("Optics/misregistration.jl")
include("Optics/deformable_mirror.jl")
include("Optics/detector.jl")
include("Optics/asterism.jl")
include("Atmosphere/kolmogorov.jl")
include("Atmosphere/multilayer.jl")
include("WFS/shack_hartmann.jl")
include("Calibration/interaction_matrix.jl")
include("Calibration/reconstructor.jl")
include("Control/controller.jl")

export PROJECT_STATUS
export AdaptiveOpticsError, InvalidConfiguration, DimensionMismatchError
export Workspace, ensure_psf_buffers!
export bin2d, poisson_noise!, poisson_sample
export Telescope, TelescopeParams, TelescopeState, generate_pupil!, reset_opd!, apply_opd!
export set_pupil!, apply_spiders!
export Source, SourceParams, wavelength
export Asterism, coordinates_xy_arcsec, compute_psf!, psf_pixel_scale_arcsec
export compute_psf!, ensure_psf_state!
export ZernikeBasis, compute_zernike!, noll_to_nm
export KolmogorovAtmosphere, KolmogorovParams, KolmogorovState
export update_psd!, ensure_psd!, phase_screen_von_karman!
export MultiLayerAtmosphere, MultiLayerParams, MultiLayerState
export advance!, propagate!
export DeformableMirror, DeformableMirrorParams, DeformableMirrorState, build_influence_functions!, apply!
export Misregistration, apply_misregistration
export Detector, DetectorParams, DetectorState, capture!
export ShackHartmann, ShackHartmannParams, ShackHartmannState, update_valid_mask!, measure!
export InteractionMatrix, interaction_matrix
export ModalReconstructor, reconstruct!, reconstruct
export AbstractController, DiscreteIntegratorController, update!
export AbstractOpticalElement, AbstractSource, AbstractAtmosphere, AbstractWFS
export AbstractDetector, AbstractDeformableMirror, SensingMode, Diffractive, Geometric
export ParallelConfig, with_parallel_config

end # module AdaptiveOptics
