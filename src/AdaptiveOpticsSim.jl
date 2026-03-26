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

include("Core/errors.jl")
include("Core/types.jl")
include("Core/profiles.jl")
include("Core/inverse_policies.jl")
include("Core/backends.jl")
include("Core/config.jl")
include("Core/utils.jl")
include("Core/workspace.jl")
include("Core/parallel.jl")
include("Core/telemetry.jl")

include("Optics/telescope.jl")
include("Optics/source.jl")
include("Optics/psf.jl")
include("Optics/zernike.jl")
include("Optics/misregistration.jl")
include("Optics/deformable_mirror.jl")
include("Optics/detector.jl")
include("Optics/asterism.jl")
include("Optics/opd_map.jl")
include("Optics/spatial_filter.jl")
include("Atmosphere/kolmogorov.jl")
include("Atmosphere/multilayer.jl")
include("Atmosphere/phase_stats.jl")
include("Calibration/modal_basis.jl")
include("Optics/ncpa.jl")
include("WFS/sensing_modes.jl")
include("WFS/elongation.jl")
include("WFS/shack_hartmann.jl")
include("WFS/pyramid.jl")
include("WFS/bioedge.jl")
include("WFS/zernike.jl")
include("WFS/curvature.jl")
include("WFS/lift.jl")
include("Calibration/interaction_matrix.jl")
include("Calibration/reconstructor.jl")
include("Calibration/calibration_vault.jl")
include("Calibration/fitting_error.jl")
include("Calibration/ao_calibration.jl")
include("Calibration/fast_atmosphere.jl")
include("Calibration/initialization.jl")
include("Calibration/gain_sensing_camera.jl")
include("Calibration/misregistration_identification.jl")
include("Control/controller.jl")
include("Control/runtime.jl")
include("Tomography/parameters.jl")
include("Tomography/geometry.jl")
include("Tomography/fitting.jl")
include("Tomography/reconstructors.jl")

export PROJECT_STATUS
export AdaptiveOpticsSimError, InvalidConfiguration, DimensionMismatchError, UnsupportedAlgorithm
export FidelityProfile, ScientificProfile, FastProfile, ProfileBundle, default_fidelity_profile
export atmosphere_profile, calibration_profile, detector_profile, lift_profile, tomography_profile
export InversePolicy, ExactPseudoInverse, TSVDInverse, TikhonovInverse, InverseStats, inverse_operator
export default_modal_inverse_policy, default_calibration_inverse_policy, default_projector_inverse_policy
export BuildBackend, NativeBuildBackend, CPUBuildBackend, GPUArrayBuildBackend, default_build_backend
export Workspace, ensure_psf_buffers!
export Telemetry, TelemetryRow, record!
export ClosedLoopTrace, ClosedLoopTraceRow
export GSCClosedLoopTrace, GSCClosedLoopTraceRow
export GSCAtmosphereReplayTrace, GSCAtmosphereReplayTraceRow
export snapshot_config, write_config_toml, write_config_json
export write_telemetry_csv
export bin2d, poisson_noise!, poisson_sample
export Telescope, TelescopeParams, TelescopeState, generate_pupil!, reset_opd!, apply_opd!
export set_pupil!, set_pupil_reflectivity!, flux_map, apply_spiders!
export Source, SourceParams, LGSSource, LGSSourceParams, wavelength, optical_path, print_optical_path
export lgs_elongation_factor
export Asterism, coordinates_xy_arcsec, compute_psf!, psf_pixel_scale_arcsec
export ensure_psf_state!
export ZernikeBasis, compute_zernike!, noll_to_nm
export OPDMap
export NCPA, NCPABasis, KLBasis, ZernikeModalBasis, M2CBasis, default_ncpa_basis
export SpatialFilter, SpatialFilterShape, CircularFilter, SquareFilter, FoucaultFilter
export set_spatial_filter!, filter!
export KolmogorovAtmosphere, KolmogorovParams, KolmogorovState
export update_psd!, ensure_psd!, phase_screen_von_karman!
export MultiLayerAtmosphere, MultiLayerParams, MultiLayerState
export advance!, propagate!
export phase_variance, phase_covariance, phase_spectrum, covariance_matrix
export ft_phase_screen, ft_sh_phase_screen, PhaseStatsWorkspace
export SubharmonicMode, FastSubharmonics, FidelitySubharmonics, default_subharmonic_mode
export DeformableMirror, DeformableMirrorParams, DeformableMirrorState, build_influence_functions!, apply!
export Misregistration, apply_misregistration
export AbstractFrameDetector, AbstractCountingDetector
export Detector, DetectorParams, DetectorState, DetectorExportMetadata
export FrameResponseModel, NullFrameResponse, SeparableGaussianPixelResponse
export APDDetector, APDDetectorParams, APDDetectorState, CountingReadoutMetadata, CountingDetectorExportMetadata
export CountingDeadTimeModel, NoDeadTime, NonParalyzableDeadTime
export capture!, output_frame, channel_output, detector_export_metadata, readout_ready, reset_integration!
export supports_detector_mtf, supports_counting_noise, supports_dead_time, supports_channel_gain_map
export SensorType, FrameSensorType, CountingSensorType, CCDSensor, CMOSSensor, EMCCDSensor, APDSensor
export ShackHartmann, ShackHartmannParams, ShackHartmannState, update_valid_mask!, measure!
export PyramidWFS, PyramidParams, PyramidState
export pyramid_modulation_frame!
export BioEdgeWFS, BioEdgeParams, BioEdgeState
export ZernikeWFS, ZernikeWFSParams, ZernikeWFSState
export CurvatureReadoutModel, CurvatureFrameReadout, CurvatureCountingReadout, CurvatureBranchResponse
export CurvatureWFS, CurvatureWFSParams, CurvatureWFSState, ensure_curvature_calibration!
export apply_shift_wfs!, set_optical_gain!
export LiFT, lift_interaction_matrix, lift_interaction_matrix!
export LiFTSolveMode, LiFTSolveAuto, LiFTSolveQR, LiFTSolveNormalEquations
export LiFTDampingMode, LiFTDampingNone, LiFTLevenbergMarquardt, LiFTAdaptiveLevenbergMarquardt
export LiFTDiagnostics, diagnostics
export InteractionMatrix, interaction_matrix
export CalibrationVault, with_truncation
export ModalBasis, KLBasisMethod, KLDMModes, KLHHtPSD
export dm_basis, kl_modal_basis, modal_basis, basis_from_m2c, basis_projector
export AOCalibration, ao_calibration
export fitting_error, fitting_error_dm
export fast_atmosphere
export AOSimulation, initialize_ao_pyramid, initialize_ao_shwfs
export GainSensingCamera, GSCDetectorMetadata, detector_metadata, weak_mode_mask, attach_detector!, detach_detector!
export calibrate!, reset_calibration!, compute_optical_gains!
export MetaSensitivity, compute_meta_sensitivity_matrix, estimate_misregistration, SPRINT, estimate!
export AbstractReconstructorOperator, ModalReconstructor, MappedReconstructor, reconstruct!, reconstruct
export AbstractController, DiscreteIntegratorController, update!
export AbstractControlSimulation, AbstractExecutionPolicy
export SequentialExecution, ThreadedExecution, BackendStreamExecution
export VectorDelayLine, shift_delay!, prepare!, prepare_runtime_wfs!, init_execution_state
export supports_prepared_runtime, supports_detector_output, supports_stacked_sources, supports_grouped_execution
export ClosedLoopRuntime, SimulationInterface, CompositeSimulationInterface, SimulationReadout
export sense!, step!, set_command!, snapshot_outputs!
export simulation_readout, simulation_slopes, simulation_command, simulation_wfs_frame, simulation_science_frame
export simulation_wfs_metadata, simulation_science_metadata
export simulation_interface
export RuntimeTimingStats, RuntimePhaseTimingStats, runtime_timing, runtime_phase_timing
export with_reconstructor, with_reconstructors
export AbstractOpticalElement, AbstractSource, AbstractAtmosphere, AbstractWFS
export AbstractDetector, AbstractDeformableMirror, SensingMode, Diffractive, Geometric
export WFSNormalization, MeanValidFluxNormalization, IncidenceFluxNormalization
export NoiseModel, NoiseNone, NoisePhoton, NoiseReadout, NoisePhotonReadout
export DMApplyMode, DMAdditive, DMReplace
export ParallelConfig, with_parallel_config
export set_fft_provider_threads!
export gpu_backend_loaded, gpu_backend_array_type, gpu_backend_name, available_gpu_backends
export disable_scalar_backend!, backend_rand, backend_randn, backend_zeros, backend_fill
export GPUBackendTag, CUDABackendTag, MetalBackendTag, AMDGPUBackendTag
export GPUPrecisionPolicy, UnifiedGPUPrecision, SplitGPUPrecision
export gpu_runtime_type, gpu_build_type, default_gpu_precision_policy, high_accuracy_gpu_precision_policy
export TomographyAtmosphereParams, LGSAsterismParams, LGSWFSParams
export TomographyParams, TomographyDMParams, TomographyFitting
export AbstractTomographyMethod, ModelBasedTomography, InteractionMatrixTomography
export TomographyNoiseModel, RelativeSignalNoise, ScalarMeasurementNoise, DiagonalMeasurementNoise, PhotonReadoutSlopeNoise
export AbstractSlopeOrder, SimulationSlopes, InterleavedSlopes, InvertedSlopes
export TomographyOperators, TomographicReconstructor, TomographyCommandReconstructor
export build_reconstructor, assemble_reconstructor_and_fitting
export airmass, layer_altitude_m, wind_direction_rad, wind_velocity_components
export lgs_height_m, lgs_directions!, lgs_directions
export optimization_geometry!, optimization_geometry, direction_vectors!, direction_vectors
export n_valid_subapertures, valid_lenslet_support, dm_valid_support, support_diameter
export sparse_gradient_matrix, auto_correlation, cross_correlation
export influence_functions, fit_commands!, fit_commands, reconstruct_wavefront!, reconstruct_wavefront
export reconstruct_wavefront_map!, reconstruct_wavefront_map
export swap_xy_blocks, interleave_xy_columns, prepare_slope_order
export mask_actuators!, dm_commands!, dm_commands

end # module AdaptiveOpticsSim
