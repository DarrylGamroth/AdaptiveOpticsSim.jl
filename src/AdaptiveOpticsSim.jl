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

export AdaptiveOpticsSimError, InvalidConfiguration, DimensionMismatchError, UnsupportedAlgorithm, NumericalConditionError
export FidelityProfile, ScientificProfile, FastProfile, ProfileBundle, default_fidelity_profile
export atmosphere_profile, calibration_profile, detector_profile, lift_profile, tomography_profile
export InversePolicy, ExactPseudoInverse, TSVDInverse, TikhonovInverse, InverseStats, inverse_operator
export CPUBackend, CUDABackend, AMDGPUBackend, MetalBackend
export backend, backend_type, same_backend, require_same_backend
export default_modal_inverse_policy, default_calibration_inverse_policy, default_projector_inverse_policy
export AbstractMaskPrimitive, MaskGrid
export CircularAperture, AnnularAperture, SpiderMask, RectangularROI, SubapertureGridMask
export default_mask_grid, pixel_mask_grid, build_mask!, apply_mask!
export AbstractTelescope, Telescope, TelescopeParams, TelescopeState, generate_pupil!, reset_opd!, apply_opd!
export set_pupil!, set_pupil_reflectivity!, flux_map, apply_spiders!, pupil_mask, opd_map
export Source, SourceParams, LGSSource, LGSSourceParams, wavelength, optical_path, print_optical_path
export SpectralSample, SpectralBundle, SpectralSource, with_spectrum
export spectral_bundle, spectral_reference_source, has_spectral_bundle, is_polychromatic, weighted_wavelength
export PointCloudSourceModel, GaussianDiskSourceModel, SampledImageSourceModel
export ExtendedSource, with_extended_source, has_extended_source_model, extended_source_model, extended_source_asterism
export lgs_elongation_factor
export ElectricField, ElectricFieldParams, ElectricFieldState
export fill_from_telescope!, fill_telescope_field!, apply_phase!, apply_amplitude!, intensity!
export AbstractPropagationModel, FraunhoferPropagation, FresnelPropagation, propagate_field!
export AbstractAtmosphericFieldModel, GeometricAtmosphericPropagation, LayeredFresnelAtmosphericPropagation
export AtmosphericFieldPropagationParams, AtmosphericFieldPropagationState, AtmosphericFieldSlice
export AtmosphericFieldPropagation, propagate_atmosphere_field!, atmospheric_intensity!
export Asterism, coordinates_xy_arcsec, compute_psf!, psf_pixel_scale_arcsec
export ensure_psf_state!
export ZernikeBasis, compute_zernike!, noll_to_nm
export OPDMap
export NCPA, NCPABasis, KLBasis, ZernikeModalBasis, M2CBasis, default_ncpa_basis
export SpatialFilter, SpatialFilterShape, CircularFilter, SquareFilter, FoucaultFilter
export set_spatial_filter!, filter!
export KolmogorovAtmosphere, KolmogorovParams, KolmogorovState
export InfinitePhaseScreenParams, InfinitePhaseScreenState, InfinitePhaseScreen
export InfiniteLayerParams, InfiniteLayerState, InfiniteAtmosphereLayer
export InfiniteMultiLayerParams, InfiniteMultiLayerState, InfiniteMultiLayerAtmosphere
export update_psd!, ensure_psd!, phase_screen_von_karman!
export MultiLayerAtmosphere, MultiLayerParams, MultiLayerState
export advance!, propagate!
export phase_variance, phase_covariance, phase_spectrum, covariance_matrix
export ft_phase_screen, ft_sh_phase_screen, PhaseStatsWorkspace
export SubharmonicMode, FastSubharmonics, FidelitySubharmonics, default_subharmonic_mode
export AbstractDMTopology, ActuatorGridTopology, SampledActuatorTopology
export AbstractDMInfluenceModel, GaussianInfluenceWidth, GaussianMechanicalCoupling, DenseInfluenceMatrix, MeasuredInfluenceFunctions
export AbstractDMActuatorModel, LinearStaticActuators, ClippedActuators, ActuatorHealthMap, CompositeDMActuatorModel
export DeformableMirror, DeformableMirrorParams, DeformableMirrorState, build_influence_functions!, apply!
export topology, topology_axis_count, topology_command_count, actuator_coordinates, valid_actuator_mask, active_actuator_indices, topology_metadata
export influence_model, actuator_model, influence_width, mechanical_coupling, influence_width_from_mechanical_coupling, n_actuators
export AbstractModalOpticBasis, FunctionModalBasis, MatrixModalBasis, ZernikeOpticBasis, CartesianTiltBasis, QuadraticFocusBasis
export CompositeControllableOptic, ModalControllableOptic, TipTiltMirror, FocusStage
export command_storage, n_control_dofs, controllable_surface_labels, supports_segmented_command, update_command!, surface_opd
export Misregistration, apply_misregistration, rotation_rad, rotation_deg, anamorphosis_angle_rad, anamorphosis_angle_deg
export AbstractFrameDetector, AbstractCountingDetector
export Detector, DetectorParams, DetectorState, DetectorExportMetadata
export AbstractDetectorResponse, AbstractFrameResponse, AbstractFrameMTF
export FrameResponseModel, NullFrameResponse, GaussianPixelResponse, SampledFrameResponse, RectangularPixelAperture, SeparablePixelMTF
export SeparableGaussianPixelResponse
export AbstractDetectorDefectModel, NullDetectorDefectModel, PixelResponseNonuniformity, DarkSignalNonuniformity, BadPixelMask
export CompositeDetectorDefectModel
export FrameWindow
export AbstractFrameTimingModel, GlobalShutter, RollingShutter
export FrameSamplingMode, SingleRead, AveragedNonDestructiveReads, CorrelatedDoubleSampling, FowlerSampling
export FrameReadoutCorrectionModel, NullFrameReadoutCorrection, ReferencePixelCommonModeCorrection
export ReferenceRowCommonModeCorrection, ReferenceColumnCommonModeCorrection, ReferenceOutputCommonModeCorrection
export CompositeFrameReadoutCorrection
export FrameReadoutProducts, NoFrameReadoutProducts, SampledFrameReadoutProducts, MultiReadFrameReadoutProducts, HgCdTeReadoutProducts
export AbstractFrameNonlinearityModel, NullFrameNonlinearity, SaturatingFrameNonlinearity
export AbstractPersistenceModel, NullPersistence, ExponentialPersistence
export AbstractDetectorThermalModel, NullDetectorThermalModel, FixedTemperature, FirstOrderThermalModel
export AbstractDetectorThermalState, NoThermalState, DetectorThermalState
export AbstractTemperatureLaw, NullTemperatureLaw, ArrheniusRateLaw, LinearTemperatureLaw, ExponentialTemperatureLaw
export APDDetector, APDDetectorParams, APDDetectorState, SPADArrayDetector, SPADArrayDetectorParams, SPADArrayDetectorState
export CountingReadoutMetadata, CountingDetectorExportMetadata
export CountingDeadTimeModel, NoDeadTime, NonParalyzableDeadTime, ParalyzableDeadTime
export AbstractCountingGateModel, NullCountingGate, DutyCycleGate
export AbstractCountingCorrelationModel, NullCountingCorrelation, AfterpulsingModel, ChannelCrosstalkModel
export CompositeCountingCorrelation
export capture!, output_frame, readout_products
export detector_reference_frame, detector_signal_frame, detector_combined_frame
export detector_reference_cube, detector_signal_cube, detector_read_cube, detector_read_times
export channel_output, detector_export_metadata, readout_ready, reset_integration!
export supports_detector_mtf, supports_clock_induced_charge, supports_column_readout_noise
export supports_avalanche_gain, supports_sensor_glow, supports_nondestructive_reads, supports_reference_read_subtraction
export supports_readout_correction, supports_read_cube
export supports_detector_defect_maps, supports_detector_persistence, supports_detector_nonlinearity, supports_shutter_timing
export supports_counting_noise, supports_dead_time, supports_channel_gain_map, supports_counting_gating
export supports_afterpulsing, supports_channel_crosstalk, supports_paralyzable_dead_time
export supports_detector_thermal_model, supports_dynamic_thermal_state
export supports_temperature_dependent_dark_current, supports_temperature_dependent_glow
export supports_temperature_dependent_persistence, supports_temperature_dependent_dark_counts
export response_family, response_application_domain, response_support
export is_shift_invariant, supports_frequency_domain_application, supports_separable_application, supports_subpixel_geometry
export default_response_model
export thermal_model, thermal_state, detector_temperature, advance_thermal!
export evaluate_temperature_law, effective_dark_current, effective_glow_rate, effective_cic_rate, effective_dark_count_rate
export SensorType, FrameSensorType, CountingSensorType, AvalancheFrameSensorType, HgCdTeAvalancheArraySensorType, SPADArraySensorType
export CCDSensor, CMOSSensor, EMCCDSensor, InGaAsSensor, HgCdTeAvalancheArraySensor, APDSensor, SPADArraySensor
export AbstractCMOSOutputModel, NullCMOSOutputModel, StaticCMOSOutputPattern
export AbstractEMGainModel, ExcessNoiseApproximation, StochasticMultiplicationRegister
export ShackHartmannWFS, ShackHartmannWFSParams, ShackHartmannWFSState, update_valid_mask!, measure!
export AbstractValidSubaperturePolicy, GeometryValidSubapertures, FluxThresholdValidSubapertures
export PyramidWFS, PyramidParams, PyramidState
export pyramid_modulation_frame!
export BioEdgeWFS, BioEdgeParams, BioEdgeState
export ZernikeWFS, ZernikeWFSParams, ZernikeWFSState
export CurvatureReadoutModel, CurvatureFrameReadout, CurvatureCountingReadout, CurvatureBranchResponse
export CurvatureWFS, CurvatureWFSParams, CurvatureWFSState, ensure_curvature_calibration!
export apply_shift_wfs!, set_optical_gain!
export valid_subaperture_mask, reference_signal, camera_frame, wfs_detector_image
export shack_hartmann_detector_image, shack_hartmann_detector_image!
export supports_valid_subaperture_mask, supports_reference_signal, supports_camera_frame
export LiFT, lift_interaction_matrix, lift_interaction_matrix!
export LiFTSolveMode, LiFTSolveAuto, LiFTSolveQR, LiFTSolveNormalEquations
export LiFTDampingMode, LiFTDampingNone, LiFTLevenbergMarquardt, LiFTAdaptiveLevenbergMarquardt
export LiFTDiagnostics, diagnostics
export InteractionMatrix, interaction_matrix
export ControlMatrix, with_truncation
export ModalBasis, KLBasisMethod, KLDMModes, KLHHtPSD
export dm_basis, kl_modal_basis, modal_basis, basis_from_m2c, basis_projector
export modal_to_command, sampled_basis, modal_projector
export AOCalibration, ao_calibration, control_matrix
export calibration_amplitude, inverse_policy, singular_values, condition_number, effective_rank, truncation_count
export AOSimulation
export forward_operator, inverse_operator_matrix
export fitting_error, fitting_error_dm
export GainSensingCamera
export calibrate!, reset_calibration!, compute_optical_gains!
export AbstractReconstructorOperator, NullReconstructor, ModalReconstructor, MappedReconstructor, reconstruct!, reconstruct
export AbstractController, DiscreteIntegratorController, update!, controller_output, reset_controller!, supports_controller_reset
export AbstractControlSimulation, AbstractExecutionPolicy
export SequentialExecution, ThreadedExecution, BackendStreamExecution
export AbstractRuntimeProfile, ScientificRuntimeProfile, HILRuntimeProfile
export RuntimeOutputRequirements, GroupedRuntimeOutputRequirements
export runtime_outputs, grouped_runtime_outputs
export RuntimeLatencyModel, default_runtime_profile, default_runtime_outputs, default_grouped_runtime_outputs
export VectorDelayLine, shift_delay!, prepare!, prepare_runtime_wfs!, init_execution_state
export supports_prepared_runtime, supports_detector_output, supports_stacked_sources, supports_grouped_execution
export ClosedLoopRuntime, SimulationInterface, CompositeSimulationInterface, SimulationReadout
export RuntimeCommandSegment, RuntimeCommandLayout, command_layout, branch_command_layout, branch_command_layouts, command_segments, command_segment_labels, command_segment_range
export AbstractControlLoopConfig, ControlLoopBranch, SingleControlLoopConfig, GroupedControlLoopConfig, ControlLoopScenario
export build_control_loop_scenario, control_loop_config, control_loop_boundary, control_loop_name, control_loop_branch_labels
export sense!, step!, set_command!, snapshot_outputs!
export readout, command, slopes, wfs_frame, science_frame, wfs_metadata, science_metadata
export grouped_wfs_stack, grouped_science_stack
export runtime_profile, runtime_latency
export RuntimeTimingStats, RuntimePhaseTimingStats, runtime_timing, runtime_phase_timing
export runtime_rng, deterministic_reference_rng
export with_reconstructor, with_reconstructors
export AbstractOpticalElement, AbstractSource, AbstractAtmosphere, AbstractWFS
export AbstractDetector, AbstractControllableOptic, AbstractDeformableMirror, SensingMode, Diffractive, Geometric
export AbstractSlopeExtractionModel, CenterOfGravityExtraction
export SubapertureLayout, SubapertureCalibration
export subaperture_layout, subaperture_calibration, slope_extraction_model
export valid_subaperture_indices, n_valid_subapertures
export WFSNormalization, MeanValidFluxNormalization, IncidenceFluxNormalization
export NoiseModel, NoiseNone, NoisePhoton, NoiseReadout, NoisePhotonReadout
export DMApplyMode, DMAdditive, DMReplace
export ParallelConfig, with_parallel_config
export AbstractArrayBackend, CPUBackend, CUDABackend, MetalBackend, AMDGPUBackend
export gpu_backend_loaded, gpu_backend_array_type, gpu_backend_name, available_gpu_backends
export array_backend_type, resolve_array_backend
export disable_scalar_backend!, backend_rand, backend_randn, backend_zeros, backend_fill
export GPUPrecisionPolicy, UnifiedGPUPrecision, SplitGPUPrecision
export gpu_runtime_type, gpu_build_type, default_gpu_precision_policy, high_accuracy_gpu_precision_policy
export TomographyAtmosphereParams, LGSAsterismParams, LGSWFSParams
export TomographyParams, TomographyDMParams, TomographyFitting
export AbstractTomographyMethod, ModelBasedTomography, InteractionMatrixTomography
export TomographyNoiseModel, RelativeSignalNoise, ScalarMeasurementNoise, DiagonalMeasurementNoise, PhotonReadoutSlopeNoise
export AbstractSlopeOrder, SimulationSlopes, InterleavedSlopes, InvertedSlopes
export TomographyOperators, TomographicReconstructor, TomographyCommandReconstructor
export build_reconstructor, assemble_reconstructor_and_fitting
export airmass, zenith_angle_rad, zenith_angle_deg, layer_altitude_m, wind_direction_rad, wind_direction_deg, wind_velocity_components
export lgs_height_m, lgs_directions!, lgs_directions
export optimization_geometry!, optimization_geometry, direction_vectors!, direction_vectors
export n_valid_subapertures, valid_lenslet_support, dm_valid_support, support_diameter
export sparse_gradient_matrix, auto_correlation, cross_correlation
export influence_functions, fit_commands!, fit_commands, reconstruct_wavefront!, reconstruct_wavefront
export reconstruct_wavefront_map!, reconstruct_wavefront_map
export swap_xy_blocks, interleave_xy_columns, prepare_slope_order
export mask_actuators!, dm_commands!, dm_commands

end # module AdaptiveOpticsSim
