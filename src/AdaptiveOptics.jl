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
include("Optics/opd_map.jl")
include("Optics/ncpa.jl")
include("Optics/spatial_filter.jl")
include("Atmosphere/kolmogorov.jl")
include("Atmosphere/multilayer.jl")
include("Atmosphere/phase_stats.jl")
include("WFS/sensing_modes.jl")
include("WFS/shack_hartmann.jl")
include("WFS/pyramid.jl")
include("WFS/bioedge.jl")
include("WFS/lift.jl")
include("Calibration/interaction_matrix.jl")
include("Calibration/reconstructor.jl")
include("Calibration/calibration_vault.jl")
include("Calibration/modal_basis.jl")
include("Calibration/fitting_error.jl")
include("Calibration/ao_calibration.jl")
include("Calibration/fast_atmosphere.jl")
include("Calibration/initialization.jl")
include("Calibration/gain_sensing_camera.jl")
include("Calibration/misregistration_identification.jl")
include("Control/controller.jl")

export PROJECT_STATUS
export AdaptiveOpticsError, InvalidConfiguration, DimensionMismatchError
export Workspace, ensure_psf_buffers!
export bin2d, poisson_noise!, poisson_sample
export Telescope, TelescopeParams, TelescopeState, generate_pupil!, reset_opd!, apply_opd!
export set_pupil!, apply_spiders!
export Source, SourceParams, LGSSource, LGSSourceParams, wavelength
export lgs_elongation_factor
export Asterism, coordinates_xy_arcsec, compute_psf!, psf_pixel_scale_arcsec
export compute_psf!, ensure_psf_state!
export ZernikeBasis, compute_zernike!, noll_to_nm
export OPDMap
export NCPA, NCPABasis, KLBasis, ZernikeModalBasis, M2CBasis
export SpatialFilter, SpatialFilterShape, CircularFilter, SquareFilter, FoucaultFilter
export set_spatial_filter!, filter!
export KolmogorovAtmosphere, KolmogorovParams, KolmogorovState
export update_psd!, ensure_psd!, phase_screen_von_karman!
export MultiLayerAtmosphere, MultiLayerParams, MultiLayerState
export advance!, propagate!
export phase_variance, phase_covariance, phase_spectrum, covariance_matrix
export ft_phase_screen, ft_sh_phase_screen, PhaseStatsWorkspace
export DeformableMirror, DeformableMirrorParams, DeformableMirrorState, build_influence_functions!, apply!
export Misregistration, apply_misregistration
export Detector, DetectorParams, DetectorState, capture!
export ShackHartmann, ShackHartmannParams, ShackHartmannState, update_valid_mask!, measure!
export PyramidWFS, PyramidParams, PyramidState
export BioEdgeWFS, BioEdgeParams, BioEdgeState
export LiFT, lift_interaction_matrix
export InteractionMatrix, interaction_matrix
export CalibrationVault, with_truncation
export ModalBasis, dm_basis, kl_modal_basis, modal_basis, basis_from_m2c, basis_projector
export AOCalibration, ao_calibration
export fitting_error, fitting_error_dm
export fast_atmosphere
export AOSimulation, initialize_ao_pyramid, initialize_ao_shwfs
export GainSensingCamera, calibrate!, reset_calibration!, compute_optical_gains!
export MetaSensitivity, compute_meta_sensitivity_matrix, estimate_misregistration, SPRINT, estimate!
export ModalReconstructor, reconstruct!, reconstruct
export AbstractController, DiscreteIntegratorController, update!
export AbstractOpticalElement, AbstractSource, AbstractAtmosphere, AbstractWFS
export AbstractDetector, AbstractDeformableMirror, SensingMode, Diffractive, Geometric
export NoiseModel, NoiseNone, NoisePhoton, NoiseReadout, NoisePhotonReadout
export DMApplyMode, DMAdditive, DMReplace
export ParallelConfig, with_parallel_config

end # module AdaptiveOptics
