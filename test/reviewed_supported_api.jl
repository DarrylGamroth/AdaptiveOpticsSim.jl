# Reviewed supported API ownership. Each binding belongs to exactly one rationale:
# exported_routine: normal user vocabulary; exported_integration: cross-package identity;
# public_advanced: qualified prepared/inspection API; public_extension: dispatch
# contracts implemented by extensions. Changes to this fixture require owner review.
# Julia's implicit module self-binding is not a declared package API member.
const REVIEWED_SUPPORTED_API = (
    AdaptiveOpticsSim = (
        exported_routine =
        (
            :AdaptiveOpticsSimError, :AlgorithmGraphs, :Atmospheres, :Backends, :Calibration,
            :Detectors, :DimensionMismatchError, :Ensembles, :FastProfile, :FidelityProfile,
            :InvalidConfiguration, :NumericalConditionError, :Optics, :ScientificProfile,
            :SplitMix64RNG, :Tomography, :UnsupportedAlgorithm, :WavefrontSensors,
            :default_fidelity_profile, :deterministic_reference_rng, :runtime_rng,
        ),
        exported_integration =
        (),
        public_advanced =
        (),
        public_extension =
        (),
    ),
    Backends = (
        exported_routine =
        (
            :AMDGPUBackend, :AbstractArrayBackend, :CPUBackend, :CUDABackend, :MetalBackend,
            :backend, :compute_device,
        ),
        exported_integration =
        (),
        public_advanced =
        (
            :AbstractComputeDevice, :AbstractComputeDeviceAvailability,
            :AcceleratorComputeDevice, :ComputeDeviceAvailable, :ComputeDeviceError,
            :ComputeDeviceUnavailable, :HostComputeDevice, :compute_device_backend,
            :compute_device_identifier,
        ),
        public_extension =
        (
            :allocate_device_array, :compute_device_availability,
            :compute_device_is_available, :compute_device_unavailable_reason,
        ),
    ),
    Optics = (
        exported_routine =
        (
            :AbstractCombinationPolicy, :AbstractFocalPlaneModulation,
            :AbstractOpticalNormalization, :AbstractOpticalPlaneKind, :AbstractOpticalProduct,
            :AbstractPlaneCoordinateDomain, :AbstractSourceRadiometry,
            :AbstractSpatialMeasure, :AchromaticSpectralCoordinate, :ActuatorGridTopology,
            :ActuatorHealthMap, :AngularCoordinates, :AnnularAperture, :Asterism,
            :BiOEdgeAmplitudeMask, :CartesianTiltBasis, :CellIntegratedMeasure,
            :CircularAperture, :CircularFilter, :CircularModulation, :ClippedActuators,
            :CoherentFieldCombination, :CompositeDMActuatorModel, :CurvatureDefocusPair,
            :DMAdditive, :DMReplace, :DeformableMirror, :DenseInfluenceMatrix, :DetectorPlane,
            :DimensionlessNormalization, :ElectricField, :FocalPlane, :FocusStage,
            :FoucaultFilter, :FraunhoferPropagation, :FresnelPropagation, :FunctionModalBasis,
            :GaussianDiskSourceModel, :GaussianInfluenceWidth, :GaussianMechanicalCoupling,
            :IncoherentIntensityAddition, :IntegratedSpectralChannel, :IntensityMap,
            :InterSampleCentered, :IntermediatePlane, :LGSSource, :MatrixModalBasis,
            :MeasuredInfluenceFunctions, :MetricCoordinates, :MicrolensArray,
            :MicrolensArrayParams, :Misregistration, :ModalControllableOptic,
            :MonochromaticChannel, :NCPA, :NoModulation, :NonCombinableProduct,
            :NormalizedPupilCoordinates, :NormalizedTestSource, :OPDMap,
            :OpticalPlaneMetadata, :OpticalProductBundle, :PhotonRateNormalization,
            :PhysicalPhotonIrradianceSource, :PlaneAxisOrientation, :PlaneCentering,
            :PointCloudSourceModel, :PointSampledMeasure, :PreparedIncoherentSum,
            :PupilFieldFormationPlan, :PupilFunction, :PupilPlane, :PyramidPhaseMask,
            :RectangularROI, :SampleCentered, :SampledActuatorTopology,
            :SampledImageSourceModel, :SampledModulation, :SodiumLayerProfile, :Source,
            :SpatialDensityMeasure, :SpatialFilter, :SpectralBundle, :SpectralSample,
            :SpectralSource, :SpiderMask, :SquareFilter, :SubapertureGridMask, :Telescope,
            :TipTiltMirror, :UnspecifiedCoherence, :UnspecifiedNormalization,
            :UnspecifiedSpatialMeasure, :UnspecifiedSpectralCoordinate, :ZernikeBasis,
            :ZernikeOpticBasis, :ZernikePhaseSpot, :accumulate_intensity!, :apply_mask!,
            :apply_misregistration, :apply_opd!, :apply_spiders!, :apply_surface!,
            :build_mask!, :compute_zernike!, :direct_imaging_components,
            :direct_imaging_output, :extended_source_asterism, :field_values,
            :fill_electric_field!, :fill_electric_field_async!, :filter!,
            :focal_plane_pixel_scale_arcsec, :form_direct_image!,
            :fraunhofer_intensity_from_field!, :influence_model, :influence_width,
            :intensity_values, :mechanical_coupling, :microlens_array, :n_actuators, :opd_map,
            :optical_path, :photon_irradiance, :plane_metadata, :prepare_direct_imaging,
            :prepare_direct_imaging_batch, :prepare_incoherent_sum,
            :prepare_microlens_propagation, :prepare_pupil_field, :prepare_spatial_filter,
            :propagate_field!, :propagation_output, :pupil_amplitude, :pupil_mask,
            :pupil_photon_rate_map, :pupil_reflectivity, :pupil_support, :reset_opd!,
            :set_command!, :set_pupil!, :set_pupil_reflectivity!, :sodium_layer_profile,
            :source_radiometric_value, :source_radiometry, :surface_opd, :update_surface!,
            :wavelength, :with_extended_source, :with_spectrum,
        ),
        exported_integration =
        (),
        public_advanced =
        (
            :DirectImagingPlan, :FraunhoferPropagationPlan, :FresnelPropagationPlan,
            :MicrolensPropagationPlan, :PreparedBundledDirectImaging, :PreparedDirectImaging,
            :PreparedDirectImagingBatch, :PreparedIncoherentDirectImaging,
            :PreparedMicrolensPropagation, :PreparedSpatialFilter, :SpatialFilterPlan,
            :StackedFraunhoferDirectImagingBatchCapability, :TelescopeDefinition,
            :UnsupportedDirectImagingBatchCapability, :direct_imaging_batch_count,
            :direct_imaging_batch_inputs, :direct_imaging_batch_products,
            :direct_imaging_batch_sources, :direct_imaging_plan, :microlens_propagation_plan,
            :prepare_telescope, :propagation_input_metadata, :propagation_output_metadata,
            :propagation_plan, :spatial_filter_output, :spatial_filter_plan,
            :validate_direct_imaging_batch,
        ),
        public_extension =
        (
            :AbstractDirectImagingBatchCapability, :AbstractPropagationModel,
            :AbstractPropagationPlan, :AbstractSource, :AbstractTelescope,
            :AbstractTelescopeDefinition, :direct_imaging_batch_capability,
            :validate_telescope_target,
        ),
    ),
    Atmospheres = (
        exported_routine =
        (
            :AbstractAtmosphere, :AtmosphereEpoch, :AtmosphereEpochError, :AtmosphereLayerID,
            :AtmosphereTimeError, :AtmosphericFieldPropagation,
            :GeometricAtmosphericPropagation, :InfiniteMultiLayerAtmosphere,
            :InfinitePhaseScreen, :KolmogorovAtmosphere,
            :LayeredFresnelAtmosphericPropagation, :MultiLayerAtmosphere, :advance!,
            :advance_by!, :advance_to!, :atmosphere_direction_output, :atmospheric_intensity!,
            :current_epoch, :direction_renderers, :epoch_sequence, :epoch_time,
            :prepare_atmosphere_direction_batch, :prepare_atmosphere_renderer,
            :prepare_atmosphere_renderers, :propagate!, :propagate_atmosphere_field!,
            :render_atmosphere!, :render_atmosphere_directions!,
        ),
        exported_integration =
        (),
        public_advanced =
        (
            :AtmosphereIdentity, :AtmosphereTimelineState,
            :InfiniteMultiLayerAtmosphereDefinition, :KolmogorovAtmosphereDefinition,
            :MultiLayerAtmosphereDefinition, :PreparedAtmosphereDirectionBatch,
            :atmosphere_direction_capacity, :atmosphere_direction_count,
            :atmosphere_direction_metadata, :atmosphere_identity, :atmosphere_timeline,
            :new_atmosphere_timeline, :prepare_timed_atmosphere,
            :validate_atmosphere_direction_batch,
        ),
        public_extension =
        (
            :AbstractTimedAtmosphere, :AbstractTimedAtmosphereDefinition,
            :atmosphere_numeric_type, :evolve_atmosphere!, :evolve_initial_atmosphere!,
            :initialize_atmosphere!, :validate_timed_atmosphere_target,
        ),
    ),
    Detectors = (
        exported_routine =
        (
            :AbstractDetectorThermalModel, :AbstractFrameResponse, :AbstractSensor,
            :ArrheniusRateLaw, :AveragedNonDestructiveReads, :BadPixelMask, :CCDSensor,
            :CMOSReadNoiseMap, :CMOSSensor, :CompositeCountingMeanResponse,
            :CompositeDetectorDefectModel, :CompositeFrameReadoutCorrection,
            :ConventionalOutput, :CorrelatedDoubleSampling, :CountingDeadTimeModel,
            :DarkSignalNonuniformity, :Detector, :DutyCycleGate, :EMCCDSensor, :EMOutput,
            :ExponentialPersistence, :ExponentialTemperatureLaw,
            :FirstOrderAfterpulseMeanResponse, :FirstOrderThermalModel, :FixedTemperature,
            :FowlerSampling, :FrameReadoutCorrectionModel, :FrameReadoutProducts,
            :FrameTransferAcquisition, :FunctionExposureFrameSource, :FunctionFrameSource,
            :GaussianPixelResponse, :GlobalResetExposure, :GlobalShutter,
            :HgCdTeAvalancheArraySensor, :HgCdTeSensor, :InGaAsSensor,
            :InPlaceExposureFrameSource, :InPlaceFrameSource, :InterpixelCapacitance,
            :LinearAPDChannelBank, :LinearAPDDetector, :LinearEMMode, :LinearTemperatureLaw,
            :MKIDArrayCharacteristics, :MKIDArrayDetector, :MKIDArraySensor,
            :MultiReadFrameReadoutProducts, :NearestNeighborCountRedistribution, :NoDeadTime,
            :NoFrameReadoutProducts, :NoiseModel, :NoiseNone, :NoisePhoton,
            :NoisePhotonReadout, :NoiseReadout, :NonParalyzableDeadTime,
            :NullDetectorThermalModel, :NullFrameReadoutCorrection, :NullFrameResponse,
            :ParalyzableDeadTime, :PhotonCountingEMMode, :PixelResponseNonuniformity,
            :RectangularPixelAperture, :ReferenceColumnCommonModeCorrection,
            :ReferenceOutputCommonModeCorrection, :ReferencePixelCommonModeCorrection,
            :ReferenceRowCommonModeCorrection, :RollingExposure, :RollingShutter,
            :SPADArrayDetector, :SPADArraySensor, :SampledFrameResponse,
            :SaturatingFrameNonlinearity, :SequentialAcquisition, :SingleElementLinearAPD,
            :SingleRead, :SkipperReadoutProducts, :SkipperSampling, :UpTheRampReadoutProducts,
            :UpTheRampSampling, :capture!, :channel_output, :detector_export_metadata,
            :detector_mtf, :detector_ramp_cube, :detector_ramp_intercept,
            :detector_ramp_read_offsets_s, :detector_ramp_slope, :emccd_snr, :output_frame,
            :prepare_detector_acquisition, :readout_ready, :reset_integration!,
            :thermal_model,
        ),
        exported_integration =
        (),
        public_advanced =
        (
            :ClippedGaussianAvalancheMultiplicationApproximation,
            :ClippedGaussianMultiplicationApproximation,
            :ConditionalGammaAvalancheMultiplication, :ConditionalGammaMultiplication,
            :DetectorAcquisitionPlan, :FrameWindow, :HgCdTeReadout, :MKIDArrayExportMetadata,
            :PreparedDetectorAcquisition, :StaticCMOSOutputPattern,
            :detector_acquisition_detector, :detector_acquisition_input,
            :detector_acquisition_plan, :detector_acquisition_products,
            :detector_acquisition_state, :detector_acquisition_workspace,
            :detector_ramp_acquisition,
        ),
        public_extension =
        (
            :AbstractEMGainModel, :AbstractHgCdTeAvalancheMultiplication,
        ),
    ),
    WavefrontSensors = (
        exported_routine =
        (
            :AbstractWFSMeasurementPath, :AcquiredObservationPath, :BiOEdgeOpticalFrontEnd,
            :BiOEdgeWFS, :CurvatureBranchResponse, :CurvatureChannelReadout,
            :CurvatureFrameReadout, :CurvatureOpticalFrontEnd, :CurvaturePackedAcquisition,
            :CurvatureReadoutModel, :CurvatureWFS, :Diffractive, :DirectMeasurementPath,
            :Geometric, :LiFTExpectedCounts, :LiFTForwardModel, :LiFTFrameMapping,
            :LiFTIdentityMapping, :LiFTNormalizedIntensity, :LiFTObservation, :LiFTPhotonRate,
            :PreparedLiFTForward, :PyramidOpticalFrontEnd, :PyramidWFS,
            :RelativeIlluminationValidSubapertures, :ShackHartmannOpticalFrontEnd,
            :ShackHartmannWFS, :SubapertureLayout, :WFSMeasurement, :WFSMeasurementMetadata,
            :WFSObservation, :WFSObservationMetadata, :WFSPreparationError,
            :ZernikeOpticalFrontEnd, :ZernikeWFS, :acquire_wfs_observation!,
            :bi_o_edge_rate_map, :curvature_rate_maps, :estimate_wfs_measurement!,
            :evaluate_lift_forward!, :form_wfs_optical_products!,
            :geometric_wavefront_slopes!, :lift_forward_output, :lift_observation_contract,
            :measurement_metadata, :measurement_storage, :measurement_units,
            :n_valid_subapertures, :observation_metadata, :observation_storage,
            :observation_units, :predict_lift_observation!, :prepare_lift_forward_model,
            :prepare_runtime_wfs!, :prepare_wfs_acquisition, :prepare_wfs_estimation,
            :prepare_wfs_optics, :pyramid_focal_mask, :pyramid_modulation_frame,
            :pyramid_modulation_frame!, :pyramid_rate_map, :set_valid_subapertures!,
            :shack_hartmann_optics, :shack_hartmann_rate_map, :subaperture_layout,
            :valid_subaperture_indices, :valid_subaperture_mask, :wfs_calibration_signature,
            :wfs_detector_image, :wfs_measurement_path, :zernike_rate_map,
        ),
        exported_integration =
        (),
        public_advanced =
        (
            :AbstractPyramidModulationPropagationStrategy, :BiOEdgeOpticsBundlePlan,
            :BiOEdgeOpticsPlan, :CurvatureOpticsPlan, :LiFTForwardPlan, :LiFTForwardWorkspace,
            :PreparedBiOEdgeOptics, :PreparedBiOEdgeOpticsBundle, :PreparedCurvatureOptics,
            :PreparedPyramidOptics, :PreparedPyramidOpticsBundle,
            :PreparedShackHartmannOptics, :PreparedShackHartmannOpticsBundle,
            :PreparedWFSCountingAcquisition, :PreparedWFSDetectorAcquisition,
            :PreparedWFSMultipleDetectorAcquisition, :PreparedZernikeOptics,
            :PyramidOpticsBundlePlan, :PyramidOpticsPlan, :PyramidPupilTiltStrategy,
            :PyramidShiftedMaskStrategy, :ShackHartmannOptics, :ShackHartmannOpticsBundlePlan,
            :ShackHartmannOpticsPlan, :WFSCountingAcquisitionPlan,
            :WFSDetectorAcquisitionPlan, :WFSMultipleDetectorAcquisitionPlan,
            :ZernikeOpticsPlan, :lift_forward_plan, :lift_forward_workspace,
            :wfs_acquisition_plan, :wfs_optical_products,
        ),
        public_extension =
        (
            :AbstractWFSAcquisitionPlan, :AbstractWFSEstimationPlan, :AbstractWFSOpticsPlan,
            :supports_detector_output, :supports_grouped_execution,
            :supports_prepared_runtime, :supports_stacked_sources, :validate_wfs_target,
        ),
    ),
    Calibration = (
        exported_routine =
        (
            :InteractionMatrix, :ModalBasis, :basis_from_m2c, :interaction_matrix,
            :modal_basis,
        ),
        exported_integration =
        (
            :KarhunenLoeveBasis,
        ),
        public_advanced =
        (
            :ModalOPDExpansionPlan, :combine_basis!,
        ),
        public_extension =
        (),
    ),
    Tomography = (
        exported_routine =
        (
            :InteractionMatrixTomography, :InterleavedSlopes, :InvertedSlopes,
            :LGSAsterismParams, :LGSWFSParams, :ModelBasedTomography, :SimulationSlopes,
            :TomographyAtmosphereParams, :TomographyDMParams, :TomographyParams,
            :assemble_reconstructor_and_fitting, :build_reconstructor, :wind_direction_deg,
            :zenith_angle_deg,
        ),
        exported_integration =
        (),
        public_advanced =
        (),
        public_extension =
        (),
    ),
    Ensembles = (
        exported_routine =
        (
            :AbstractExecutionPolicy, :AcceleratedKernelsExecution, :BackendStreamExecution,
            :DaggerExecution, :DeterministicExecution, :SequentialExecution,
            :SimulationEnsemble, :ThreadedExecution, :ensemble_members, :execution_policy,
            :run_ensemble!,
        ),
        exported_integration =
        (),
        public_advanced =
        (),
        public_extension =
        (
            :ensemble_ownership_roots, :execute_ensemble!, :init_ensemble_scheduler,
            :init_execution_state,
        ),
    ),
    AlgorithmGraphs = (
        exported_routine =
        (
            :AlgorithmGraphDefinition, :AlgorithmGraphError, :AlgorithmLink,
            :AlgorithmNodeDefinition, :CapturedGraphExecution, :CapturedModelTimestamp,
            :DelayedAlgorithmLink, :FixedStepModelTimeDriver, :GraphInputDefinition,
            :GraphOutputDefinition, :GraphStepTicket, :GroupedStreamGraphExecution,
            :ModelDuration, :ModelTimestamp, :PeriodicSchedule, :PreparedAlgorithmGraph,
            :PreparedBoundaryModelTimeDriver, :PreparedCapturedModelTimeDriver,
            :PreparedGraphHILBoundary, :StreamGraphExecution, :adopt_hil_command!,
            :advance_model_time!, :algorithm_graph, :algorithm_node,
            :builtin_graph_node_types, :capture_model_time_origin, :capture_model_timestamp,
            :captured_graph_node_count, :ccd_detector_acquisition_node,
            :cmos_detector_acquisition_node, :deformable_mirror_surface_node, :delayed_link,
            :emccd_detector_acquisition_node, :gaussian_deformable_mirror_surface_node,
            :graph_execution_policy, :graph_failed, :graph_input, :graph_name, :graph_output,
            :graph_step_pending, :graph_step_sequence,
            :grid_gaussian_deformable_mirror_surface_node, :hil_boundary_status,
            :hil_command_buffer, :hil_frame_buffer, :link, :load_algorithm_graph,
            :modal_opd_expansion_node, :model_duration_seconds, :model_nanoseconds,
            :model_time_exhausted, :model_time_provenance, :model_time_seconds,
            :model_time_sequence, :model_time_uncertainty, :model_timestamp,
            :multilayer_atmosphere_opd_node, :next_model_time_capture, :next_model_timestamp,
            :prepare_algorithm_graph, :prepare_boundary_model_time_driver,
            :prepare_captured_model_time_driver, :prepare_graph_hil_boundary,
            :pupil_opd_composition_node, :pyramid_rate_node, :reset_graph!,
            :reset_hil_boundary!, :reset_model_time!, :schedule_period, :schedule_phase,
            :schedule_timestamp, :shack_hartmann_rate_node, :sparse_parameter, :step_graph!,
            :step_graph_async!, :step_graph_at!, :step_hil_frame!, :step_hil_frame_at!,
            :wait_graph_step!,
        ),
        exported_integration =
        (),
        public_advanced =
        (
            :CCDDetectorAcquisitionNode, :CCDDetectorAcquisitionNodeConfig,
            :CMOSDetectorAcquisitionNode, :CMOSDetectorAcquisitionNodeConfig,
            :DeformableMirrorSurfaceNode, :DeformableMirrorSurfaceNodeConfig,
            :EMCCDDetectorAcquisitionNode, :EMCCDDetectorAcquisitionNodeConfig,
            :GaussianDeformableMirrorSurfaceNode, :GaussianDeformableMirrorSurfaceNodeConfig,
            :GridGaussianDeformableMirrorSurfaceNode,
            :GridGaussianDeformableMirrorSurfaceNodeConfig, :ModalOPDExpansionNode,
            :ModalOPDExpansionNodeConfig, :MultiLayerAtmosphereOPDNode,
            :MultiLayerAtmosphereOPDNodeConfig, :PupilOPDCompositionNode,
            :PupilOPDCompositionNodeConfig, :PyramidRateNode, :PyramidRateNodeConfig,
            :ShackHartmannRateNode, :ShackHartmannRateNodeConfig, :prepared_graph_node,
        ),
        public_extension =
        (
            :GraphNodeCaptureSafe, :GraphNodeCaptureUnsupported, :enqueue_graph_node!,
            :graph_node_capture_capability, :graph_node_ports, :graph_port_contract,
            :prepare_graph_node, :reset_graph_node!, :step_graph_node!,
        ),
    ),
)

@testset "Reviewed supported API allowlists" begin
    for (owner, contract) in pairs(REVIEWED_SUPPORTED_API)
        mod = owner === :AdaptiveOpticsSim ? AdaptiveOpticsSim :
            getfield(AdaptiveOpticsSim, owner)
        exported = vcat(collect(contract.exported_routine),
            collect(contract.exported_integration))
        qualified = vcat(collect(contract.public_advanced),
            collect(contract.public_extension))
        @test length(exported) == length(unique(exported))
        @test length(qualified) == length(unique(qualified))
        @test isempty(intersect(exported, qualified))
        @test Set(exported) == Set(filter(name ->
            name !== nameof(mod) && Base.isexported(mod, name),
            names(mod; all=true)))
        @test Set(qualified) == Set(filter(name ->
            Base.ispublic(mod, name) && !Base.isexported(mod, name),
            names(mod; all=true)))
    end
end
