# Curated public API exported by `using AdaptiveOpticsSim`.
#
# Advanced and developer-facing names remain accessible as
# `AdaptiveOpticsSim.name`. Add names here only when they are part of ordinary
# user workflows, maintained examples, or stable extension seams.

export AdaptiveOpticsSimError, InvalidConfiguration, DimensionMismatchError
export UnsupportedAlgorithm, NumericalConditionError

export FidelityProfile, ScientificProfile, FastProfile, default_fidelity_profile
export runtime_rng, deterministic_reference_rng

export InversePolicy, ExactPseudoInverse, TSVDInverse, TikhonovInverse, default_modal_inverse_policy
export CPUBackend, CUDABackend, AMDGPUBackend, MetalBackend
export AbstractArrayBackend, backend

export CircularAperture, AnnularAperture, SpiderMask, RectangularROI, SubapertureGridMask
export build_mask!, apply_mask!

export AbstractAtmosphere
export Telescope, Source, LGSSource, Asterism
export wavelength, optical_path
export reset_opd!, apply_opd!, set_pupil!, set_pupil_reflectivity!
export pupil_mask, apply_spiders!
export compute_psf!, psf_pixel_scale_arcsec

export SpectralSample, SpectralBundle, SpectralSource, with_spectrum
export GaussianDiskSourceModel, PointCloudSourceModel, SampledImageSourceModel
export with_extended_source, extended_source_asterism

export ElectricField, FraunhoferPropagation, FresnelPropagation
export GeometricAtmosphericPropagation, LayeredFresnelAtmosphericPropagation
export AtmosphericFieldPropagation

export ZernikeBasis, compute_zernike!
export OPDMap, Misregistration, apply_misregistration
export NCPA, KLBasis, ZernikeModalBasis
export SpatialFilter, CircularFilter, SquareFilter, FoucaultFilter, filter!

export KolmogorovAtmosphere, MultiLayerAtmosphere, InfinitePhaseScreen, InfiniteMultiLayerAtmosphere
export advance!, propagate!

export ActuatorGridTopology, SampledActuatorTopology
export GaussianInfluenceWidth, GaussianMechanicalCoupling, DenseInfluenceMatrix, MeasuredInfluenceFunctions
export ClippedActuators, ActuatorHealthMap, CompositeDMActuatorModel
export DeformableMirror, apply!
export influence_model, influence_width
export mechanical_coupling, n_actuators
export FunctionModalBasis, MatrixModalBasis, ZernikeOpticBasis, CartesianTiltBasis
export CompositeControllableOptic, ModalControllableOptic, TipTiltMirror, FocusStage
export update_command!
export DMAdditive, DMReplace

export Detector, APDDetector, SPADArrayDetector
export NoiseModel, NoiseNone, NoisePhoton, NoiseReadout, NoisePhotonReadout
export SensorType, CCDSensor, CMOSSensor, EMCCDSensor, InGaAsSensor
export LinearEMMode, PhotonCountingEMMode, EMOutput, ConventionalOutput, emccd_snr
export QCMOSSensor, ORCAQuestSensor, ORCAQuest2Sensor, ORCAQuestIQSensor
export ORCAQuest, ORCAQuest2, ORCAQuestIQ
export QCMOSStandardScan, QCMOSUltraQuietScan, QCMOSPhotonNumberResolvingScan, QCMOSRawScan
export QCMOSDetector, ORCAQuestDetector, ORCAQuest2Detector, ORCAQuestIQDetector
export HgCdTeAvalancheArraySensor, APDSensor, SPADArraySensor
export FrameResponseModel, NullFrameResponse, GaussianPixelResponse, SampledFrameResponse
export RectangularPixelAperture, SeparablePixelMTF
export PixelResponseNonuniformity, DarkSignalNonuniformity, BadPixelMask, CompositeDetectorDefectModel
export RollingShutter, RollingExposure, GlobalResetExposure, CorrelatedDoubleSampling, FowlerSampling
export FunctionFrameSource, InPlaceFrameSource, FunctionExposureFrameSource, InPlaceExposureFrameSource
export FrameReadoutCorrectionModel, NullFrameReadoutCorrection
export ReferencePixelCommonModeCorrection, ReferenceRowCommonModeCorrection
export ReferenceColumnCommonModeCorrection, ReferenceOutputCommonModeCorrection
export CompositeFrameReadoutCorrection
export FrameReadoutProducts, NoFrameReadoutProducts, MultiReadFrameReadoutProducts, HgCdTeReadoutProducts
export SaturatingFrameNonlinearity, ExponentialPersistence
export AbstractDetectorThermalModel, NullDetectorThermalModel, FixedTemperature, FirstOrderThermalModel
export ArrheniusRateLaw, LinearTemperatureLaw, ExponentialTemperatureLaw
export CountingDeadTimeModel, NoDeadTime, NonParalyzableDeadTime, ParalyzableDeadTime
export DutyCycleGate, AfterpulsingModel, ChannelCrosstalkModel, CompositeCountingCorrelation
export capture!, output_frame, detector_export_metadata
export readout_ready, reset_integration!, thermal_model

export Diffractive, Geometric
export ShackHartmannWFS, PyramidWFS, BioEdgeWFS, ZernikeWFS, CurvatureWFS
export CurvatureReadoutModel, CurvatureCountingReadout, CurvatureBranchResponse
export FluxThresholdValidSubapertures
export CenterOfGravityExtraction, SubapertureLayout, SubapertureCalibration
export subaperture_layout, subaperture_calibration, slope_extraction_model
export valid_subaperture_indices
export MeanValidFluxNormalization, IncidenceFluxNormalization
export measure!, pyramid_modulation_frame!
export valid_subaperture_mask, camera_frame, wfs_detector_image
export shack_hartmann_detector_image, shack_hartmann_detector_image!
export n_valid_subapertures
export LiFT
export LiFTSolveAuto, LiFTSolveQR, LiFTSolveNormalEquations
export LiFTLevenbergMarquardt, LiFTAdaptiveLevenbergMarquardt

export InteractionMatrix, interaction_matrix
export ControlMatrix
export ModalBasis, KLDMModes, KLHHtPSD
export kl_modal_basis, modal_basis, basis_from_m2c
export AOCalibration, ao_calibration, control_matrix
export fitting_error
export GainSensingCamera, calibrate!, compute_optical_gains!

export NullReconstructor, ModalReconstructor, MappedReconstructor
export reconstruct!, reconstruct
export DiscreteIntegratorController

export AbstractControlSimulation, AbstractExecutionPolicy, SimulationReadout
export AOSimulation, SequentialExecution, ThreadedExecution, BackendStreamExecution
export ScientificRuntimeProfile, HILRuntimeProfile, RuntimeOutputRequirements, GroupedRuntimeOutputRequirements
export VectorDelayLine, shift_delay!
export prepare!, prepare_runtime_wfs!
export command_segments, command_segment_range
export ControlLoopBranch, SingleControlLoopConfig, GroupedControlLoopConfig, ControlLoopScenario
export build_control_loop_scenario, control_loop_name, control_loop_branch_labels
export sense!, step!, set_command!
export readout, command, slopes, wfs_frame, science_frame, grouped_wfs_stack
export runtime_timing

export TomographyAtmosphereParams, LGSAsterismParams, LGSWFSParams
export TomographyParams, TomographyDMParams
export ModelBasedTomography, InteractionMatrixTomography
export SimulationSlopes, InterleavedSlopes, InvertedSlopes
export build_reconstructor, assemble_reconstructor_and_fitting
export zenith_angle_deg, wind_direction_deg
export reconstruct_wavefront_map
export dm_commands
