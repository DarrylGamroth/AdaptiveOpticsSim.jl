# Curated package API.
#
# Export only names used routinely after `using AdaptiveOpticsSim`. Declare
# stable advanced and developer-facing seams `public` so callers can use
# `AdaptiveOpticsSim.name` or import them explicitly without adding them to the
# ordinary user namespace.

export AdaptiveOpticsSimError, InvalidConfiguration, DimensionMismatchError
export UnsupportedAlgorithm, NumericalConditionError
export Backends, Optics, Detectors, Atmospheres, WavefrontSensors, Plant

export FidelityProfile, ScientificProfile, FastProfile, default_fidelity_profile
export runtime_rng, deterministic_reference_rng

export InversePolicy, ExactPseudoInverse, TSVDInverse, TikhonovInverse, default_modal_inverse_policy
export KLBasis, ZernikeModalBasis

export ShackHartmannDirectFrontEnd, ShackHartmannOpticalFrontEnd
export shack_hartmann_rate_map
export PyramidOpticalFrontEnd, pyramid_rate_map
export BioEdgeOpticalFrontEnd, bioedge_rate_map
export set_pyramid_calibration!, set_bioedge_calibration!
export ZernikeOpticalFrontEnd, zernike_rate_map
export CurvatureOpticalFrontEnd, curvature_rate_maps
export CurvaturePackedAcquisition
export set_zernike_calibration!, set_curvature_calibration!
export ShackHartmannWFS, PyramidWFS, BioEdgeWFS, ZernikeWFS, CurvatureWFS
export CurvatureReadoutModel, CurvatureFrameReadout
export CurvatureCountingReadout, CurvatureBranchResponse
export FluxThresholdValidSubapertures
export AbstractSlopeExtractionModel, CenterOfGravityExtraction
export SubapertureLayout, SubapertureCalibration
export subaperture_layout, subaperture_calibration, slope_extraction_model
export set_subaperture_calibration!
export valid_subaperture_indices
export pyramid_modulation_frame, pyramid_modulation_frame!
export shack_hartmann_detector_image, shack_hartmann_detector_image!
export n_valid_subapertures
export LiFT, PreparedLiFTForwardModel, LiFTObservation
export LiFTIdentityMapping, LiFTFrameMapping
export LiFTPhotonRate, LiFTExpectedCounts, LiFTNormalizedIntensity
export prepare_lift_forward_model, evaluate_lift_forward!
export predict_lift_observation!, lift_forward_output
export lift_observation_contract, diagnostics
export LiFTSolveAuto, LiFTSolveQR, LiFTSolveNormalEquations
export LiFTLevenbergMarquardt, LiFTAdaptiveLevenbergMarquardt

export InteractionMatrix, interaction_matrix
export ControlMatrix
export ModalBasis, KLDMModes, KLHHtPSD
export kl_modal_basis, modal_basis, basis_from_m2c
export AOCalibration, ao_calibration, control_matrix
export fitting_error
export GainSensingCamera, calibrate!, compute_optical_gains!

export NullReconstructor, ModalReconstructor, FactorizedReconstructor, MappedReconstructor
export ControlledReconstructor
export reconstruct!, reconstruct
export DiscreteIntegratorController, VectorDelayLine, shift_delay!

export AbstractExecutionPolicy
export SequentialExecution, ThreadedExecution, BackendStreamExecution
export DeterministicExecution, AcceleratedKernelsExecution, DaggerExecution
export SimulationEnsemble
export runtime_timing

public controller_output, reset_controller!, supports_controller_reset
public run_ensemble!, ensemble_members, execution_policy
public ensemble_ownership_roots, init_ensemble_scheduler, execute_ensemble!
public init_execution_state

export TomographyAtmosphereParams, LGSAsterismParams, LGSWFSParams
export TomographyParams, TomographyDMParams
export ModelBasedTomography, InteractionMatrixTomography
export SimulationSlopes, InterleavedSlopes, InvertedSlopes
export build_reconstructor, assemble_reconstructor_and_fitting
export zenith_angle_deg, wind_direction_deg
export reconstruct_wavefront_map
export dm_commands
