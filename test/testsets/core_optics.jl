const TestAbstractFFTs = AdaptiveOpticsSim.AbstractFFTs

mutable struct MutableTelescopeTestDefinition <: AbstractTelescopeDefinition
end

struct TestExpandedSourceWrapper{S<:AdaptiveOpticsSim.Optics.AbstractSource} <:
    AdaptiveOpticsSim.Optics.AbstractSource
    source::S
end

AdaptiveOpticsSim.Optics.source_composition_style(::TestExpandedSourceWrapper) =
    AdaptiveOpticsSim.Optics.ExpandedSourceComposition()

# DenseArray is deliberately broad enough to reproduce the dispatch category
# used by accelerator array packages without depending on a GPU runtime.
struct TestBackendFFTArray{T,N} <: DenseArray{T,N}
    storage::Array{T,N}
end

Base.size(array::TestBackendFFTArray) = size(array.storage)
Base.getindex(array::TestBackendFFTArray, indices...) =
    getindex(array.storage, indices...)
Base.setindex!(array::TestBackendFFTArray, value, indices...) =
    setindex!(array.storage, value, indices...)
TestAbstractFFTs.plan_fft!(::TestBackendFFTArray) = :backend_fft
TestAbstractFFTs.plan_fft!(::TestBackendFFTArray, dims) = (:backend_fft, dims)
TestAbstractFFTs.plan_ifft!(::TestBackendFFTArray) = :backend_ifft
TestAbstractFFTs.plan_ifft!(::TestBackendFFTArray, dims) = (:backend_ifft, dims)

@testset "GPU backend registry" begin
    @test !gpu_backend_loaded(AdaptiveOpticsSim.Backends.CUDABackendTag)
    @test !gpu_backend_loaded(AdaptiveOpticsSim.Backends.MetalBackendTag)
    @test !gpu_backend_loaded(AdaptiveOpticsSim.Backends.AMDGPUBackendTag)
    @test gpu_backend_array_type(AdaptiveOpticsSim.Backends.CUDABackendTag) === nothing
    @test gpu_backend_array_type(AdaptiveOpticsSim.Backends.MetalBackendTag) === nothing
    @test gpu_backend_array_type(AdaptiveOpticsSim.Backends.AMDGPUBackendTag) === nothing
    @test gpu_backend_name(Matrix{Float64}) === nothing
    @test backend(zeros(2, 2)) isa CPUBackend
    @test backend(zeros(2)) isa CPUBackend
    @test available_gpu_backends() == ()
    @test Calibration.GPUArrayBuildBackend(AdaptiveOpticsSim.Backends.CUDABackendTag) isa
        Calibration.GPUArrayBuildBackend{AdaptiveOpticsSim.Backends.CUDABackendTag}
end

@testset "CPU FFT provider dispatch" begin
    for T in (ComplexF32, ComplexF64)
        original = T.(reshape(1:64, 8, 8))

        full = copy(original)
        full_fft = AdaptiveOpticsSim.Backends.plan_fft_backend!(full)
        full_ifft = AdaptiveOpticsSim.Backends.plan_ifft_backend!(full)
        @test full_fft isa TestAbstractFFTs.Plan
        @test full_ifft isa TestAbstractFFTs.Plan
        AdaptiveOpticsSim.Backends.execute_fft_plan!(full, full_fft)
        AdaptiveOpticsSim.Backends.execute_fft_plan!(full, full_ifft)
        @test full ≈ original

        first_dimension = copy(original)
        dim_fft = AdaptiveOpticsSim.Backends.plan_fft_backend!(first_dimension, 1)
        dim_ifft = AdaptiveOpticsSim.Backends.plan_ifft_backend!(first_dimension, 1)
        @test dim_fft isa TestAbstractFFTs.Plan
        @test dim_ifft isa TestAbstractFFTs.Plan
        AdaptiveOpticsSim.Backends.execute_fft_plan!(first_dimension, dim_fft)
        AdaptiveOpticsSim.Backends.execute_fft_plan!(first_dimension, dim_ifft)
        @test first_dimension ≈ original
    end

    backend_array = TestBackendFFTArray(zeros(ComplexF32, 4, 4))
    @test AdaptiveOpticsSim.Backends.plan_fft_backend!(backend_array) === :backend_fft
    @test AdaptiveOpticsSim.Backends.plan_ifft_backend!(backend_array) === :backend_ifft
    @test AdaptiveOpticsSim.Backends.plan_fft_backend!(backend_array, (1, 2)) ==
        (:backend_fft, (1, 2))
    @test AdaptiveOpticsSim.Backends.plan_ifft_backend!(backend_array, (1, 2)) ==
        (:backend_ifft, (1, 2))
end

function explicit_direct_image_cycle!(prepared)
    return form_direct_image!(prepared)
end

function explicit_spatial_filter_cycle!(prepared)
    return filter!(prepared)
end

@testset "API export curation" begin
    root_exported = filter(name -> Base.isexported(AdaptiveOpticsSim, name),
        names(AdaptiveOpticsSim))
    optics_exported = filter(name -> Base.isexported(Optics, name),
        names(Optics))
    detectors_exported = filter(name -> Base.isexported(Detectors, name),
        names(Detectors))
    atmospheres_exported = filter(name -> Base.isexported(Atmospheres, name),
        names(Atmospheres))
    plant_exported = filter(name -> Base.isexported(Plant, name),
        names(Plant))

    # The root package exposes the domain module, while the domain itself
    # distinguishes routine unqualified vocabulary from stable qualified API.
    @test length(root_exported) <= 500
    @test length(optics_exported) <= 160
    @test length(detectors_exported) <= 105
    @test length(atmospheres_exported) <= 35
    @test length(plant_exported) <= 100
    @test Base.isexported(AdaptiveOpticsSim, :Backends)
    @test Base.isexported(AdaptiveOpticsSim, :Optics)
    @test Base.isexported(AdaptiveOpticsSim, :Detectors)
    @test Base.isexported(AdaptiveOpticsSim, :Atmospheres)
    @test Base.isexported(AdaptiveOpticsSim, :Plant)
    @test Base.ispublic(AdaptiveOpticsSim, :Backends)
    @test Base.ispublic(AdaptiveOpticsSim, :Optics)
    @test Base.ispublic(AdaptiveOpticsSim, :Detectors)
    @test Base.ispublic(AdaptiveOpticsSim, :Atmospheres)
    @test Base.ispublic(AdaptiveOpticsSim, :Plant)
    @test AdaptiveOpticsSim.Backends === Backends
    @test AdaptiveOpticsSim.Optics === Optics
    @test AdaptiveOpticsSim.Detectors === Detectors
    @test AdaptiveOpticsSim.Atmospheres === Atmospheres
    @test AdaptiveOpticsSim.Plant === Plant

    # Every supported Plant-owned binding has one canonical owner. The root
    # module exposes only the domain module, never forwarding compatibility
    # bindings for exported or qualified-public Plant API.
    for name in names(Plant; imported=false)
        name === :Plant && continue
        @test !isdefined(AdaptiveOpticsSim, name)
    end

    for name in (
        :PlantTimestamp,
        :PlantDuration,
        :PeriodicSchedule,
        :OpticalPathID,
        :AcquisitionID,
        :ControllableOpticID,
        :CommandEndpointID,
        :PlantCommandSchema,
        :PlantCommand,
        :PlantDefinition,
        :PreparedPlant,
        :prepare_plant,
        :prepare_acquisition_selection,
        :execute_acquisition_selection!,
        :PlantEventLoopDefinition,
        :prepare_plant_event_loop,
    )
        @test !isdefined(AdaptiveOpticsSim, name)
        @test !Base.isexported(AdaptiveOpticsSim, name)
        @test Base.isexported(Plant, name)
        @test Base.ispublic(Plant, name)
    end

    for name in (
        :CommandValuePolicy,
        :CommandEndpointState,
        :CommandApplicationState,
        :PreparedPathExecutor,
        :PreparedTriggerTopology,
        :EventSchedulerWorkspace,
        :PlantEventKey,
        :RNGOwnerIdentity,
        :SingleIllumination,
        :FullOpticalProviderStyle,
        :validate_path_execution_binding,
    )
        @test !isdefined(AdaptiveOpticsSim, name)
        @test !Base.isexported(Plant, name)
        @test Base.ispublic(Plant, name)
    end

    for name in (
        :_PreparedTriggerSource,
        :_PLANT_COMMAND_CLAIM_TOKEN,
        :_PreparedPlantToken,
    )
        @test isdefined(Plant, name)
        @test !Base.isexported(Plant, name)
        @test !Base.ispublic(Plant, name)
    end

    @test parentmodule(Plant.PlantCommand) === Plant
    @test parentmodule(Plant.PlantTimestamp) === Plant
    @test parentmodule(Plant.PlantDefinitionError) === Plant
    @test Plant.advance_to! === AdaptiveOpticsSim.Atmospheres.advance_to!
    @test Plant.backend === AdaptiveOpticsSim.Backends.backend
    @test !Base.isexported(AdaptiveOpticsSim, :ShackHartmannWFS)
    @test Base.isexported(WavefrontSensors, :ShackHartmannWFS)
    for name in (
        :ZernikeBasis,
        :Misregistration,
        :SpatialFilter,
        :prepare_spatial_filter,
        :DeformableMirror,
        :influence_model,
        :ModalControllableOptic,
        :CircularModulation,
        :SampledModulation,
        :MicrolensArray,
        :microlens_array,
        :prepare_microlens_propagation,
        :PyramidPhaseMask,
        :BiOEdgeAmplitudeMask,
        :ZernikePhaseSpot,
        :CurvatureDefocusPair,
    )
        @test !Base.isexported(AdaptiveOpticsSim, name)
        @test !Base.ispublic(AdaptiveOpticsSim, name)
        @test Base.isexported(Optics, name)
        @test Base.ispublic(Optics, name)
        @test parentmodule(getfield(Optics, name)) === Optics
    end
    @test !Base.isexported(AdaptiveOpticsSim,
        :PreparedMicrolensPropagation)
    @test !Base.isexported(AdaptiveOpticsSim,
        :ShackHartmannDirectFrontEnd)
    @test Base.isexported(WavefrontSensors,
        :ShackHartmannDirectFrontEnd)
    @test !Base.isexported(AdaptiveOpticsSim,
        :ShackHartmannOpticalFrontEnd)
    @test Base.isexported(WavefrontSensors,
        :ShackHartmannOpticalFrontEnd)
    @test !Base.isexported(AdaptiveOpticsSim, :shack_hartmann_rate_map)
    @test Base.isexported(WavefrontSensors, :shack_hartmann_rate_map)
    for name in (
        :PyramidWFS,
        :BiOEdgeWFS,
        :PyramidOpticalFrontEnd,
        :BiOEdgeOpticalFrontEnd,
        :pyramid_rate_map,
        :bi_o_edge_rate_map,
        :set_pyramid_calibration!,
        :set_bi_o_edge_calibration!,
        :pyramid_modulation_frame,
        :pyramid_modulation_frame!,
    )
        @test !Base.isexported(AdaptiveOpticsSim, name)
        @test Base.isexported(WavefrontSensors, name)
    end
    for name in (
        :ZernikeOpticalFrontEnd,
        :CurvatureOpticalFrontEnd,
        :CurvaturePackedAcquisition,
    )
        @test !Base.isexported(AdaptiveOpticsSim, name)
        @test !Base.ispublic(AdaptiveOpticsSim, name)
        @test !isdefined(AdaptiveOpticsSim, name)
        @test Base.isexported(WavefrontSensors, name)
    end
    for name in (
        :PreparedLiFTForwardModel,
        :LiFTObservation,
        :LiFTIdentityMapping,
        :LiFTFrameMapping,
        :LiFTPhotonRate,
        :LiFTExpectedCounts,
        :LiFTNormalizedIntensity,
        :prepare_lift_forward_model,
        :evaluate_lift_forward!,
        :predict_lift_observation!,
        :lift_forward_output,
        :lift_observation_contract,
        :diagnostics,
    )
        @test !Base.isexported(AdaptiveOpticsSim, name)
        @test !Base.ispublic(AdaptiveOpticsSim, name)
        @test !isdefined(AdaptiveOpticsSim, name)
        @test Base.isexported(WavefrontSensors, name)
    end
    @test !Base.isexported(AdaptiveOpticsSim,
        :PreparedFocalPlaneModulation)
    @test !Base.isexported(AdaptiveOpticsSim,
        :set_subaperture_calibration!)
    @test Base.isexported(WavefrontSensors,
        :set_subaperture_calibration!)
    for name in (
        :Detector,
        :CMOSReadNoiseMap,
        :SkipperSampling,
        :GlobalShutter,
        :SingleRead,
        :AveragedNonDestructiveReads,
        :UpTheRampSampling,
        :FrameTransferAcquisition,
        :detector_ramp_cube,
        :detector_ramp_read_offsets_s,
        :InterpixelCapacitance,
        :LinearAPDDetector,
        :SPADArrayDetector,
        :MKIDArrayDetector,
        :SPADArraySensor,
        :MKIDArraySensor,
        :SingleElementLinearAPD,
        :LinearAPDChannelBank,
        :CountingDeadTimeModel,
        :DutyCycleGate,
        :NearestNeighborCountRedistribution,
        :channel_output,
    )
        @test !Base.isexported(AdaptiveOpticsSim, name)
        @test !Base.ispublic(AdaptiveOpticsSim, name)
        @test Base.isexported(Detectors, name)
        @test Base.ispublic(Detectors, name)
        @test parentmodule(getfield(Detectors, name)) === Detectors
    end
    @test Base.isexported(Control, :FactorizedReconstructor)
    @test Base.isexported(Control, :ControlledReconstructor)
    @test !Base.isexported(AdaptiveOpticsSim, :PyramidWFS)
    @test !Base.isexported(AdaptiveOpticsSim, :BiOEdgeWFS)
    @test Base.isexported(WavefrontSensors, :PyramidWFS)
    @test Base.isexported(WavefrontSensors, :BiOEdgeWFS)
    @test !Base.isexported(AdaptiveOpticsSim, :ZernikeWFS)
    @test !Base.isexported(AdaptiveOpticsSim, :CurvatureWFS)
    @test !isdefined(AdaptiveOpticsSim, :ZernikeWFS)
    @test !isdefined(AdaptiveOpticsSim, :CurvatureWFS)
    @test Base.isexported(WavefrontSensors, :ZernikeWFS)
    @test Base.isexported(WavefrontSensors, :CurvatureWFS)
    @test !Base.isexported(AdaptiveOpticsSim, :prepare_runtime_wfs!)
    @test !Base.ispublic(AdaptiveOpticsSim, :prepare_runtime_wfs!)
    @test Base.isexported(WavefrontSensors, :prepare_runtime_wfs!)
    @test parentmodule(WavefrontSensors.prepare_runtime_wfs!) ===
        WavefrontSensors
    @test !Base.isexported(AdaptiveOpticsSim, :subaperture_layout)
    @test Base.isexported(WavefrontSensors, :subaperture_layout)
    @test !Base.isexported(AdaptiveOpticsSim, :compute_device)
    @test !Base.ispublic(AdaptiveOpticsSim, :AbstractComputeDevice)
    @test Base.isexported(Backends, :compute_device)
    @test Base.ispublic(Backends, :AbstractComputeDevice)
    @test Base.ispublic(Backends, :HostComputeDevice)
    @test Base.ispublic(Backends, :AcceleratorComputeDevice)
    @test Base.ispublic(Backends, :compute_device_backend)
    @test Base.ispublic(Backends, :compute_device_identifier)
    for name in (
        :AbstractComputeDeviceAvailability,
        :ComputeDeviceAvailable,
        :ComputeDeviceUnavailable,
        :ComputeDeviceError,
        :compute_device_availability,
        :compute_device_is_available,
        :compute_device_unavailable_reason,
        :allocate_device_array,
    )
        @test !Base.isexported(AdaptiveOpticsSim, name)
        @test !Base.ispublic(AdaptiveOpticsSim, name)
        @test Base.ispublic(Backends, name)
        @test parentmodule(getfield(Backends, name)) === Backends
    end
    @test !Base.isexported(Backends, :allocate_array)
    @test !Base.ispublic(Backends, :allocate_array)
    for name in (
        :Telescope,
        :photon_irradiance,
        :OpticalPlaneMetadata,
        :MetricCoordinates,
        :AngularCoordinates,
        :AchromaticSpectralCoordinate,
        :MonochromaticChannel,
        :IntegratedSpectralChannel,
        :PupilFunction,
        :ElectricField,
        :IntensityMap,
        :PhotonRateNormalization,
        :DimensionlessNormalization,
        :PointSampledMeasure,
        :SpatialDensityMeasure,
        :CellIntegratedMeasure,
        :CoherentFieldCombination,
        :IncoherentIntensityAddition,
        :NonCombinableProduct,
        :PhysicalPhotonIrradianceSource,
        :NormalizedTestSource,
        :source_radiometry,
        :OpticalProductBundle,
        :prepare_incoherent_sum,
        :accumulate_intensity!,
        :prepare_pupil_field,
        :prepare_direct_imaging,
        :prepare_direct_imaging_batch,
        :form_direct_image!,
    )
        @test !Base.isexported(AdaptiveOpticsSim, name)
        @test !Base.ispublic(AdaptiveOpticsSim, name)
        @test Base.isexported(Optics, name)
        @test Base.ispublic(Optics, name)
    end
    for name in (
        :AbstractPropagationPlan,
        :FraunhoferPropagationPlan,
        :FraunhoferPropagationWorkspace,
        :FresnelPropagationPlan,
        :FresnelPropagationWorkspace,
        :DirectImagingPlan,
        :DirectImagingWorkspace,
        :PreparedDirectImaging,
        :SpatialFilterPlan,
        :SpatialFilterWorkspace,
        :PreparedSpatialFilter,
        :PreparedDirectImagingBatch,
        :DirectImagingBatchCompatibilitySignature,
        :DirectImagingBatchWorkspace,
        :direct_imaging_batch_capability,
        :direct_imaging_batch_signature,
        :validate_direct_imaging_batch,
    )
        @test !Base.isexported(AdaptiveOpticsSim, name)
        @test !Base.ispublic(AdaptiveOpticsSim, name)
        @test !Base.isexported(Optics, name)
        @test Base.ispublic(Optics, name)
    end
    for name in (
        :AbstractAtmosphere,
        :AtmosphereEpoch,
        :KolmogorovAtmosphere,
        :MultiLayerAtmosphere,
        :InfinitePhaseScreen,
        :AtmosphericFieldPropagation,
        :advance_by!,
        :advance_to!,
        :prepare_atmosphere_renderer,
        :prepare_atmosphere_direction_batch,
        :render_atmosphere_directions!,
        :atmosphere_direction_output,
        :render_atmosphere!,
    )
        @test !Base.isexported(AdaptiveOpticsSim, name)
        @test !Base.ispublic(AdaptiveOpticsSim, name)
        @test Base.isexported(Atmospheres, name)
        @test Base.ispublic(Atmospheres, name)
        @test parentmodule(getfield(Atmospheres, name)) === Atmospheres
    end
    for name in (
        :AbstractTimedAtmosphere,
        :PreparedAtmosphereDirectionBatch,
        :AbstractAtmosphereDirectionBatchCapability,
        :validate_atmosphere_direction_batch,
    )
        @test !Base.isexported(AdaptiveOpticsSim, name)
        @test !Base.ispublic(AdaptiveOpticsSim, name)
        @test !Base.isexported(Atmospheres, name)
        @test Base.ispublic(Atmospheres, name)
        @test parentmodule(getfield(Atmospheres, name)) === Atmospheres
    end
    @test !Base.isexported(AdaptiveOpticsSim, :TelescopeParams)
    @test !Base.isexported(AdaptiveOpticsSim, :TelescopeState)
    @test !Base.isexported(AdaptiveOpticsSim, :DetectorParams)
    @test !Base.isexported(AdaptiveOpticsSim, :DetectorState)
    @test !Base.isexported(AdaptiveOpticsSim, :supports_detector_mtf)
    @test !Base.isexported(AdaptiveOpticsSim, :supports_prepared_runtime)
    @test !Base.isexported(AdaptiveOpticsSim, :backend_rand)
    @test !Base.isexported(AdaptiveOpticsSim, :backend_zeros)
    @test !Base.isexported(AdaptiveOpticsSim, :BuildBackend)
    @test !Base.isexported(AdaptiveOpticsSim, :CPUBuildBackend)
    @test !Base.isexported(AdaptiveOpticsSim, :GPUArrayBuildBackend)
    @test !Base.isexported(AdaptiveOpticsSim, :default_build_backend)
    @test !Base.isexported(AdaptiveOpticsSim, :set_fft_provider_threads!)
    @test !Base.isexported(AdaptiveOpticsSim, :GPUBackendTag)
    @test !Base.isexported(AdaptiveOpticsSim, :AbstractRuntimeExecutionPlan)
    @test !Base.isexported(AdaptiveOpticsSim, :runtime_reconstructor_storage)
    for name in (
        :NullReconstructor,
        :ModalReconstructor,
        :FactorizedReconstructor,
        :MappedReconstructor,
        :ControlledReconstructor,
        :reconstruct!,
        :reconstruct,
        :DiscreteIntegratorController,
        :VectorDelayLine,
        :shift_delay!,
    )
        @test !Base.isexported(AdaptiveOpticsSim, name)
        @test !Base.ispublic(AdaptiveOpticsSim, name)
        @test Base.isexported(Control, name)
        @test Base.ispublic(Control, name)
        @test parentmodule(getfield(Control, name)) === Control
    end
    for name in (
        :controller_output,
        :reset_controller!,
        :supports_controller_reset,
    )
        @test !Base.isexported(AdaptiveOpticsSim, name)
        @test !Base.ispublic(AdaptiveOpticsSim, name)
        @test !Base.isexported(Control, name)
        @test Base.ispublic(Control, name)
        @test parentmodule(getfield(Control, name)) === Control
    end
    for name in (
        :TomographyAtmosphereParams,
        :LGSAsterismParams,
        :LGSWFSParams,
        :TomographyParams,
        :TomographyDMParams,
        :ModelBasedTomography,
        :InteractionMatrixTomography,
        :SimulationSlopes,
        :InterleavedSlopes,
        :InvertedSlopes,
        :build_reconstructor,
        :assemble_reconstructor_and_fitting,
        :zenith_angle_deg,
        :wind_direction_deg,
        :reconstruct_wavefront_map,
        :dm_commands,
    )
        @test !Base.isexported(AdaptiveOpticsSim, name)
        @test !Base.ispublic(AdaptiveOpticsSim, name)
        @test !isdefined(AdaptiveOpticsSim, name)
        @test Base.isexported(Tomography, name)
        @test Base.ispublic(Tomography, name)
        @test parentmodule(getfield(Tomography, name)) === Tomography
    end
    for name in (
        :AbstractExecutionPolicy,
        :SequentialExecution,
        :ThreadedExecution,
        :BackendStreamExecution,
        :DeterministicExecution,
        :AcceleratedKernelsExecution,
        :DaggerExecution,
        :SimulationEnsemble,
    )
        @test !Base.isexported(AdaptiveOpticsSim, name)
        @test !Base.ispublic(AdaptiveOpticsSim, name)
        @test !isdefined(AdaptiveOpticsSim, name)
        @test Base.isexported(Ensembles, name)
        @test Base.ispublic(Ensembles, name)
        @test parentmodule(getfield(Ensembles, name)) === Ensembles
    end
    for name in (
        :run_ensemble!,
        :ensemble_members,
        :execution_policy,
        :ensemble_ownership_roots,
        :init_ensemble_scheduler,
        :execute_ensemble!,
        :init_execution_state,
    )
        @test !Base.isexported(AdaptiveOpticsSim, name)
        @test !Base.ispublic(AdaptiveOpticsSim, name)
        @test !isdefined(AdaptiveOpticsSim, name)
        @test !Base.isexported(Ensembles, name)
        @test Base.ispublic(Ensembles, name)
        @test parentmodule(getfield(Ensembles, name)) === Ensembles
    end
    @test !Base.isexported(AdaptiveOpticsSim, :CUDABackendTag)
    @test !Base.isexported(AdaptiveOpticsSim, :MetalBackendTag)
    @test !Base.isexported(AdaptiveOpticsSim, :AMDGPUBackendTag)
    for removed_name in (
        :RuntimeCommandLayout,
        :RuntimeCommandSegment,
        :CompositeControllableOptic,
        :update_command!,
        :AbstractControlSimulation,
        :AOSimulation,
        :ClosedLoopRuntime,
        :ControlLoopScenario,
        :CPUHILExecutionPlan,
        :DeviceResidentExecutionPlan,
        :RuntimeOutputRequirements,
        :SharedOpticalRuntime,
        :wfs_source,
        :science_source,
        :synchronize_runtime!,
        :plane_device,
        :AbstractPlaneDevice,
        :HostPlaneDevice,
        :AcceleratorPlaneDevice,
    )
        @test !isdefined(AdaptiveOpticsSim, removed_name)
        @test !Base.isexported(AdaptiveOpticsSim, removed_name)
    end
    @test Calibration.CPUBuildBackend() isa Calibration.BuildBackend
end

@testset "Optical radiometry and combination contracts" begin
    for invalid_wavelength in (0.0, -1.0, Inf, NaN)
        @test_throws InvalidConfiguration MonochromaticChannel(
            invalid_wavelength)
    end
    @test_throws InvalidConfiguration IntegratedSpectralChannel(Symbol(""))

    metadata_storage = zeros(2, 2)
    for invalid_sampling in (
        (0.0, 1.0),
        (-1.0, 1.0),
        (Inf, 1.0),
        (NaN, 1.0),
    )
        @test_throws InvalidConfiguration OpticalPlaneMetadata(
            FocalPlane(), metadata_storage;
            coordinate_domain=AngularCoordinates(),
            sampling=invalid_sampling)
    end
    for invalid_origin in ((Inf, 0.0), (NaN, 0.0))
        @test_throws InvalidConfiguration OpticalPlaneMetadata(
            FocalPlane(), metadata_storage;
            coordinate_domain=AngularCoordinates(), sampling=(1.0, 1.0),
            origin=invalid_origin)
    end

    host_view = @view metadata_storage[:, :]
    host_wrappers = (
        host_view,
        reshape(host_view, 1, 4),
        transpose(host_view),
        PermutedDimsArray(host_view, (2, 1)),
        reinterpret(Float32, host_view),
    )
    host_device = compute_device(metadata_storage)
    @test host_device == AdaptiveOpticsSim.Backends.HostComputeDevice()
    @test AdaptiveOpticsSim.Backends.compute_device_backend(host_device) ==
        CPUBackend()
    @test isnothing(
        AdaptiveOpticsSim.Backends.compute_device_identifier(host_device))
    host_availability = @inferred(
        AdaptiveOpticsSim.Backends.compute_device_availability(host_device))
    @test typeof(host_availability) ===
        AdaptiveOpticsSim.Backends.ComputeDeviceAvailable
    @test AdaptiveOpticsSim.Backends.compute_device_is_available(
        host_availability)
    @test isnothing(
        AdaptiveOpticsSim.Backends.compute_device_unavailable_reason(
            host_availability))
    host_context =
        AdaptiveOpticsSim.Backends._prepare_device_execution_context(
            host_device)
    prepared_host_device =
        AdaptiveOpticsSim.Backends._prepared_device_execution_compute_device(
            host_context)
    @test prepared_host_device == host_device
    exact_host_storage = AdaptiveOpticsSim.Backends.allocate_device_array(
        host_device, Float32, 2, 3)
    @test exact_host_storage isa Matrix{Float32}
    @test size(exact_host_storage) == (2, 3)
    @test compute_device(exact_host_storage) == host_device
    cuda_device_0 = AdaptiveOpticsSim.Backends.AcceleratorComputeDevice(
        CUDABackend(), 0)
    cuda_device_1 = AdaptiveOpticsSim.Backends.AcceleratorComputeDevice(
        CUDABackend(), 1)
    amd_device_0 = AdaptiveOpticsSim.Backends.AcceleratorComputeDevice(
        AMDGPUBackend(), 0)
    @test cuda_device_0 != cuda_device_1
    @test cuda_device_0 != amd_device_0
    @test AdaptiveOpticsSim.Backends.compute_device_backend(cuda_device_0) ==
        CUDABackend()
    @test AdaptiveOpticsSim.Backends.compute_device_identifier(cuda_device_0) == 0
    unavailable_cuda = @inferred(
        AdaptiveOpticsSim.Backends.compute_device_availability(cuda_device_0))
    @test typeof(unavailable_cuda) ===
        AdaptiveOpticsSim.Backends.ComputeDeviceUnavailable
    @test !AdaptiveOpticsSim.Backends.compute_device_is_available(
        unavailable_cuda)
    @test AdaptiveOpticsSim.Backends.compute_device_unavailable_reason(
        unavailable_cuda) == :exact_device_selection_unavailable
    unavailable_error = try
        AdaptiveOpticsSim.Backends.allocate_device_array(
            cuda_device_0, Float32, 2, 3)
        nothing
    catch error
        error
    end
    @test typeof(unavailable_error) <:
        AdaptiveOpticsSim.Backends.ComputeDeviceError
    @test unavailable_error.operation == :select
    @test unavailable_error.reason == :exact_device_selection_unavailable
    @test unavailable_error.device == cuda_device_0
    unavailable_context_error = try
        AdaptiveOpticsSim.Backends._prepare_device_execution_context(
            cuda_device_0)
        nothing
    catch error
        error
    end
    @test typeof(unavailable_context_error) <:
        AdaptiveOpticsSim.Backends.ComputeDeviceError
    @test unavailable_context_error.operation == :prepare_context
    @test unavailable_context_error.reason ==
        :exact_device_selection_unavailable
    @test unavailable_context_error.device == cuda_device_0
    @test_throws InvalidConfiguration begin
        AdaptiveOpticsSim.Backends.AcceleratorComputeDevice(CPUBackend(), 0)
    end
    @test_throws InvalidConfiguration begin
        AdaptiveOpticsSim.Backends.AcceleratorComputeDevice(
            CUDABackend(), nothing)
    end
    @test_throws InvalidConfiguration begin
        AdaptiveOpticsSim.Backends.AcceleratorComputeDevice(
            CUDABackend(), "device-0")
    end
    for wrapper in host_wrappers
        @test compute_device(wrapper) == host_device
        @test backend(wrapper) isa CPUBackend
        @test AdaptiveOpticsSim.Backends.array_backend_selector(typeof(wrapper)) isa
            CPUBackend
        wrapper_metadata = OpticalPlaneMetadata(FocalPlane(), wrapper;
            coordinate_domain=AngularCoordinates(), sampling=(1.0, 1.0),
            normalization=PhotonRateNormalization(),
            spatial_measure=CellIntegratedMeasure(),
            coherence=IncoherentIntensityAddition())
        wrapper_map = IntensityMap(wrapper_metadata, wrapper)
        @test wrapper_map.values === wrapper
    end
    if coverage_instrumented()
        @test_skip "exact-device availability allocation gate disabled under coverage instrumentation"
    else
        AdaptiveOpticsSim.Backends.compute_device_availability(host_device)
        @test @allocated(
            AdaptiveOpticsSim.Backends.compute_device_availability(
                host_device)) == 0
    end

    tel = Telescope(resolution=8, diameter=8.0, central_obstruction=0.0,
        pupil_reflectivity=0.25)
    wavefront = PupilFunction(tel)
    @test wavefront.metadata.spectral isa AchromaticSpectralCoordinate

    physical_source = Source(band=:custom, wavelength=1.0e-6,
        photon_irradiance=3.0)
    physical_field = ElectricField(wavefront, physical_source;
        zero_padding=2)
    physical_formation = prepare_pupil_field(wavefront,
        physical_source, physical_field)
    fill_electric_field!(physical_field, wavefront, physical_formation)
    expected_rate = photon_irradiance(physical_source) *
        prod(wavefront.metadata.sampling) * sum(abs2, wavefront.amplitude)
    @test source_radiometry(physical_source) isa
        PhysicalPhotonIrradianceSource
    @test physical_field.metadata.normalization isa PhotonRateNormalization
    @test physical_field.metadata.spatial_measure isa CellIntegratedMeasure
    @test physical_field.metadata.coherence isa CoherentFieldCombination
    @test sum(abs2, physical_field.values) ≈ expected_rate

    physical_propagation = FraunhoferPropagation(physical_field)
    physical_intensity = IntensityMap(physical_field, physical_propagation)
    fraunhofer_intensity_from_field!(physical_intensity, physical_field,
        physical_propagation)
    @test physical_intensity.metadata.normalization isa
        PhotonRateNormalization
    @test physical_intensity.metadata.spatial_measure isa
        CellIntegratedMeasure
    @test physical_intensity.metadata.coherence isa
        IncoherentIntensityAddition
    @test sum(physical_intensity.values) ≈ expected_rate

    normalized_source = Source(band=:custom, wavelength=1.0e-6,
        normalized_power=2.5)
    normalized_field = ElectricField(wavefront, normalized_source;
        zero_padding=2)
    normalized_formation = prepare_pupil_field(wavefront,
        normalized_source, normalized_field)
    fill_electric_field!(normalized_field, wavefront, normalized_formation)
    @test source_radiometry(normalized_source) isa NormalizedTestSource
    @test_throws InvalidConfiguration photon_irradiance(normalized_source)
    @test normalized_field.metadata.normalization isa
        DimensionlessNormalization
    @test normalized_field.metadata.spatial_measure isa CellIntegratedMeasure
    @test sum(abs2, normalized_field.values) ≈ 2.5

    source_bundle = SpectralBundle(
        [wavelength(physical_source), 1.1 * wavelength(physical_source)],
        [0.5, 0.5])
    spectral_leaf = with_spectrum(physical_source, source_bundle)
    directional = Asterism([physical_source])
    extended = with_extended_source(physical_source,
        PointCloudSourceModel([(0.0, 0.0)], [1.0]))
    third_party_expansion = TestExpandedSourceWrapper(physical_source)
    @test_throws UnsupportedAlgorithm with_spectrum(spectral_leaf,
        source_bundle)
    @test_throws UnsupportedAlgorithm with_spectrum(directional,
        source_bundle)
    @test_throws UnsupportedAlgorithm with_spectrum(extended,
        source_bundle)
    @test_throws UnsupportedAlgorithm SpectralSource(physical_source,
        source_bundle)
    @test_throws UnsupportedAlgorithm SpectralSource(spectral_leaf,
        source_bundle)
    @test_throws UnsupportedAlgorithm Asterism([third_party_expansion])
    @test_throws UnsupportedAlgorithm Asterism(
        AdaptiveOpticsSim.Optics.AbstractSource[
            physical_source,
            third_party_expansion,
        ])

    normalized_propagation = FraunhoferPropagation(normalized_field)
    normalized_intensity = IntensityMap(normalized_field,
        normalized_propagation)
    fraunhofer_intensity_from_field!(normalized_intensity, normalized_field,
        normalized_propagation)
    @test normalized_intensity.metadata.normalization isa
        DimensionlessNormalization
    @test sum(normalized_intensity.values) ≈ 2.5

    normalized_direct = IntensityMap(normalized_field,
        normalized_propagation)
    @test_throws InvalidConfiguration prepare_direct_imaging(wavefront,
        normalized_source, normalized_field, normalized_direct)

    @test_throws InvalidConfiguration pupil_photon_rate_map(tel,
        normalized_source)
    @test_throws InvalidConfiguration prepare_direct_imaging(wavefront,
        normalized_source; zero_padding=2)

    float32_overflow = 2 * Float64(floatmax(Float32))
    for invalid_wavelength in (0.0, -1.0, Inf, -Inf, NaN,
        float32_overflow)
        @test_throws InvalidConfiguration Source(band=:custom,
            wavelength=invalid_wavelength, normalized_power=1.0,
            T=Float32)
        @test_throws InvalidConfiguration LGSSource(
            wavelength=invalid_wavelength, normalized_power=1.0,
            T=Float32)
    end
    for invalid_value in (-1.0, Inf, -Inf, NaN, float32_overflow)
        @test_throws InvalidConfiguration Source(band=:custom,
            wavelength=1.0e-6, photon_irradiance=invalid_value,
            T=Float32)
        @test_throws InvalidConfiguration Source(band=:custom,
            wavelength=1.0e-6, normalized_power=invalid_value,
            T=Float32)
        @test_throws InvalidConfiguration LGSSource(
            photon_irradiance=invalid_value, T=Float32)
    end
    @test source_radiometric_value(Source(band=:custom,
        wavelength=1.0e-6, normalized_power=0.0, T=Float32)) === 0.0f0

    for invalid_finite_value in (Inf, -Inf, NaN, float32_overflow)
        @test_throws InvalidConfiguration Source(band=:custom,
            wavelength=1.0e-6, magnitude=invalid_finite_value,
            normalized_power=1.0, T=Float32)
        @test_throws InvalidConfiguration Source(band=:custom,
            wavelength=1.0e-6,
            coordinates=(invalid_finite_value, 0.0),
            normalized_power=1.0, T=Float32)
        @test_throws InvalidConfiguration Source(band=:custom,
            wavelength=1.0e-6,
            coordinates=(1.0, invalid_finite_value),
            normalized_power=1.0, T=Float32)
        @test_throws InvalidConfiguration LGSSource(
            magnitude=invalid_finite_value, T=Float32)
        @test_throws InvalidConfiguration LGSSource(
            coordinates=(invalid_finite_value, 0.0), T=Float32)
        @test_throws InvalidConfiguration LGSSource(
            laser_coordinates=(invalid_finite_value, 0.0), T=Float32)
    end
    for invalid_altitude in (0.0, -1.0, Inf, -Inf, NaN,
        float32_overflow)
        @test_throws InvalidConfiguration LGSSource(
            altitude=invalid_altitude, T=Float32)
    end
    for invalid_nonnegative in (-1.0, Inf, -Inf, NaN,
        float32_overflow)
        @test_throws InvalidConfiguration LGSSource(
            elongation_factor=invalid_nonnegative, T=Float32)
        @test_throws InvalidConfiguration LGSSource(
            fwhm_spot_up=invalid_nonnegative, T=Float32)
    end
    @test LGSSource(elongation_factor=0.0,
        fwhm_spot_up=0.0).params.elongation_factor == 0.0

    valid_profile = SodiumLayerProfile(
        [80_000.0, 90_000.0, 100_000.0], [0.25, 0.5, 0.25])
    profile_source = LGSSource(sodium_layer_profile=valid_profile,
        T=Float32)
    @test profile_source.params.altitude ≈ 90_000.0f0
    @test profile_source.params.sodium_layer_profile == SodiumLayerProfile(
        Float32[80_000, 90_000, 100_000], Float32[0.25, 0.5, 0.25])
    @test sodium_layer_profile(profile_source) ===
        profile_source.params.sodium_layer_profile
    @test_throws DimensionMismatchError SodiumLayerProfile(
        [80_000.0, 90_000.0], [1.0])
    @test_throws InvalidConfiguration LGSSource(
        sodium_layer_profile=SodiumLayerProfile(Float64[], Float64[]),
        T=Float32)
    @test_throws InvalidConfiguration LGSSource(
        sodium_layer_profile=SodiumLayerProfile(
            [80_000.0, 90_000.0], [0.0, 0.0]), T=Float32)
    @test_throws InvalidConfiguration LGSSource(
        sodium_layer_profile=SodiumLayerProfile(
            [0.0, 90_000.0], [0.5, 0.5]), T=Float32)
    @test_throws InvalidConfiguration LGSSource(
        sodium_layer_profile=SodiumLayerProfile(
            [80_000.0, Inf], [0.5, 0.5]), T=Float32)
    @test_throws InvalidConfiguration LGSSource(
        sodium_layer_profile=SodiumLayerProfile(
            [80_000.0, float32_overflow], [0.5, 0.5]), T=Float32)
    @test_throws InvalidConfiguration LGSSource(
        sodium_layer_profile=SodiumLayerProfile(
            [80_000.0, 90_000.0], [-0.1, 1.1]), T=Float32)
    @test_throws InvalidConfiguration LGSSource(
        sodium_layer_profile=SodiumLayerProfile(
            [80_000.0, 90_000.0], [NaN, 1.0]), T=Float32)
    @test_throws InvalidConfiguration LGSSource(
        sodium_layer_profile=SodiumLayerProfile(
            [80_000.0, 90_000.0], [Inf, 1.0]), T=Float32)
    @test_throws InvalidConfiguration LGSSource(
        sodium_layer_profile=SodiumLayerProfile(
            [80_000.0, 90_000.0],
            [floatmax(Float32), floatmax(Float32)]), T=Float32)

    @test_throws InvalidConfiguration GeometricAtmosphericPropagation(
        chromatic_reference_wavelength=Inf, T=Float32)
    @test_throws InvalidConfiguration LayeredFresnelAtmosphericPropagation(
        chromatic_reference_wavelength=-1.0, T=Float32)
    for invalid_wavelength in (0.0, -1.0, Inf, -Inf, NaN,
        float32_overflow)
        @test_throws InvalidConfiguration SpectralBundle(
            [invalid_wavelength], [1.0]; T=Float32)
    end
    for invalid_weight in (-1.0, Inf, -Inf, NaN, float32_overflow)
        @test_throws InvalidConfiguration SpectralBundle(
            [1.0e-6], [invalid_weight]; T=Float32)
    end
    @test_throws InvalidConfiguration SpectralBundle(
        SpectralSample{Float32}[
            SpectralSample{Float32}(1.0f-6, floatmax(Float32)),
            SpectralSample{Float32}(1.1f-6, floatmax(Float32)),
        ])
    parameterized_samples = SpectralSample{Float64}[
        SpectralSample(1.0e-6, 0.25),
        SpectralSample(1.1e-6, 0.75),
    ]
    parameterized_bundle = SpectralBundle{Float64,
        typeof(parameterized_samples)}(parameterized_samples)
    @test parameterized_bundle.samples == parameterized_samples
    @test_throws InvalidConfiguration SpectralBundle{Float64,
        Vector{SpectralSample{Float64}}}([
            SpectralSample(NaN, -1.0),
        ])
    @test_throws InvalidConfiguration SpectralBundle{Float64,
        Vector{SpectralSample{Float64}}}([
            SpectralSample(1.0e-6, 2.0),
        ])

    density_field = ElectricField(wavefront, physical_source;
        normalization=DimensionlessNormalization(),
        spatial_measure=SpatialDensityMeasure(),
        coherence=CoherentFieldCombination())
    @test_throws InvalidConfiguration prepare_pupil_field(wavefront,
        physical_source, density_field; amplitude_scale=1.0)

    noncombinable_field = ElectricField(wavefront, physical_source;
        normalization=DimensionlessNormalization(),
        spatial_measure=CellIntegratedMeasure(),
        coherence=NonCombinableProduct())
    @test_throws InvalidConfiguration prepare_pupil_field(wavefront,
        physical_source, noncombinable_field; amplitude_scale=1.0)
    noncombinable_propagation = FraunhoferPropagation(noncombinable_field)
    @test_throws InvalidConfiguration IntensityMap(noncombinable_field,
        noncombinable_propagation)

    unspecified_values = zeros(ComplexF64, 8, 8)
    unspecified_metadata = OpticalPlaneMetadata(PupilPlane(),
        unspecified_values; coordinate_domain=MetricCoordinates(),
        sampling=(1.0, 1.0),
        normalization=PhotonRateNormalization(),
        spatial_measure=CellIntegratedMeasure(),
        coherence=CoherentFieldCombination())
    unspecified_field = ElectricField(unspecified_metadata,
        unspecified_values)
    @test_throws InvalidConfiguration FraunhoferPropagation(
        unspecified_field)
    @test_throws InvalidConfiguration FresnelPropagation(unspecified_field;
        distance_m=1.0)

    function intensity_map(values;
        sampling=(1.0, 1.0), wavelength_m=1.0e-6,
        normalization=PhotonRateNormalization(),
        spatial_measure=CellIntegratedMeasure(),
        coherence=IncoherentIntensityAddition())
        metadata = OpticalPlaneMetadata(FocalPlane(), values;
            coordinate_domain=AngularCoordinates(), sampling=sampling,
            spectral=MonochromaticChannel(wavelength_m),
            normalization=normalization, spatial_measure=spatial_measure,
            coherence=coherence)
        return IntensityMap(metadata, values)
    end

    output = intensity_map(zeros(4, 4))
    first_input = intensity_map(fill(1.0, 4, 4))
    second_input = intensity_map(fill(2.0, 4, 4))
    sum_plan = prepare_incoherent_sum(output, first_input, second_input)
    @test !applicable(PreparedIncoherentSum, output.metadata,
        (first_input.metadata, second_input.metadata))
    @test !applicable(
        PreparedIncoherentSum{typeof(output.metadata),
            typeof((first_input.metadata, second_input.metadata))},
        output.metadata, (first_input.metadata, second_input.metadata))
    @test @inferred(accumulate_intensity!(output,
        (first_input, second_input), sum_plan)) === output
    @test output.values == fill(3.0, 4, 4)
    @test @allocated(accumulate_intensity!(output,
        (first_input, second_input), sum_plan)) == 0

    integrated_metadata = OpticalPlaneMetadata(FocalPlane(), zeros(4, 4);
        coordinate_domain=AngularCoordinates(), sampling=(1.0, 1.0),
        spectral=IntegratedSpectralChannel(:science_passband),
        normalization=PhotonRateNormalization(),
        spatial_measure=CellIntegratedMeasure(),
        coherence=IncoherentIntensityAddition())
    integrated_output = IntensityMap(integrated_metadata, zeros(4, 4))
    integrated_first = IntensityMap(integrated_metadata, fill(1.0, 4, 4))
    integrated_second = IntensityMap(integrated_metadata, fill(2.0, 4, 4))
    integrated_plan = prepare_incoherent_sum(integrated_output,
        integrated_first, integrated_second)
    accumulate_intensity!(integrated_output,
        (integrated_first, integrated_second), integrated_plan)
    @test integrated_output.values == fill(3.0, 4, 4)

    achromatic_metadata = OpticalPlaneMetadata(FocalPlane(), zeros(4, 4);
        coordinate_domain=AngularCoordinates(), sampling=(1.0, 1.0),
        spectral=AchromaticSpectralCoordinate(),
        normalization=PhotonRateNormalization(),
        spatial_measure=CellIntegratedMeasure(),
        coherence=IncoherentIntensityAddition())
    achromatic_output = IntensityMap(achromatic_metadata, zeros(4, 4))
    achromatic_first = IntensityMap(achromatic_metadata, fill(1.0, 4, 4))
    achromatic_second = IntensityMap(achromatic_metadata, fill(2.0, 4, 4))
    achromatic_plan = prepare_incoherent_sum(achromatic_output,
        achromatic_first, achromatic_second)
    accumulate_intensity!(achromatic_output,
        (achromatic_first, achromatic_second), achromatic_plan)
    @test achromatic_output.values == fill(3.0, 4, 4)
    @test_throws InvalidConfiguration prepare_incoherent_sum(
        achromatic_output, first_input)

    other_integrated_metadata = OpticalPlaneMetadata(FocalPlane(),
        fill(1.0, 4, 4); coordinate_domain=AngularCoordinates(),
        sampling=(1.0, 1.0),
        spectral=IntegratedSpectralChannel(:other_passband),
        normalization=PhotonRateNormalization(),
        spatial_measure=CellIntegratedMeasure(),
        coherence=IncoherentIntensityAddition())
    other_integrated = IntensityMap(other_integrated_metadata,
        fill(1.0, 4, 4))
    @test_throws InvalidConfiguration prepare_incoherent_sum(
        integrated_output, other_integrated)

    unspecified_sum_values = zeros(4, 4)
    unspecified_sum_metadata = OpticalPlaneMetadata(FocalPlane(),
        unspecified_sum_values; coordinate_domain=AngularCoordinates(),
        sampling=(1.0, 1.0), normalization=PhotonRateNormalization(),
        spatial_measure=CellIntegratedMeasure(),
        coherence=IncoherentIntensityAddition())
    unspecified_sum = IntensityMap(unspecified_sum_metadata,
        unspecified_sum_values)
    @test_throws InvalidConfiguration prepare_incoherent_sum(
        unspecified_sum, unspecified_sum)
    identity_alias = IntensityMap{
        typeof(first_input.metadata),typeof(output.values),CPUBackend,
    }(first_input.metadata, output.values)
    view_values = @view output.values[:, :]
    view_alias = IntensityMap{
        typeof(first_input.metadata),typeof(view_values),CPUBackend,
    }(first_input.metadata, view_values)
    @test_throws InvalidConfiguration prepare_incoherent_sum(output,
        identity_alias)
    @test_throws InvalidConfiguration prepare_incoherent_sum(output,
        view_alias)
    fill!(output.values, 7.0)
    @test_throws InvalidConfiguration accumulate_intensity!(output,
        (view_alias, second_input), sum_plan)
    @test output.values == fill(7.0, 4, 4)
    unplanned_input = intensity_map(fill(1.0, 4, 4);
        sampling=(2.0, 2.0))
    @test_throws InvalidConfiguration accumulate_intensity!(output,
        (unplanned_input, second_input), sum_plan)
    @test output.values == fill(7.0, 4, 4)
    @test !applicable(accumulate_intensity!, output.values, physical_field)
    @test_throws DimensionMismatchError accumulate_intensity!(output,
        (first_input,), sum_plan)
    @test_throws InvalidConfiguration prepare_incoherent_sum(output)

    incompatible_normalization = intensity_map(fill(1.0, 4, 4);
        normalization=DimensionlessNormalization())
    incompatible_measure = intensity_map(fill(1.0, 4, 4);
        spatial_measure=SpatialDensityMeasure())
    incompatible_policy = intensity_map(fill(1.0, 4, 4);
        coherence=NonCombinableProduct())
    incompatible_sampling = intensity_map(fill(1.0, 4, 4);
        sampling=(2.0, 2.0))
    incompatible_spectral = intensity_map(fill(1.0, 4, 4);
        wavelength_m=1.1e-6)
    for incompatible in (
        incompatible_normalization,
        incompatible_measure,
        incompatible_policy,
        incompatible_sampling,
        incompatible_spectral,
    )
        @test_throws InvalidConfiguration prepare_incoherent_sum(output,
            incompatible)
    end

    bundle = OpticalProductBundle(first_input, incompatible_spectral)
    @test length(bundle) == 2
    @test bundle[1] === first_input
    @test bundle[2] === incompatible_spectral

    mixed_wavelengths = Asterism([
        Source(band=:custom, wavelength=1.0e-6, photon_irradiance=1.0),
        Source(band=:custom, wavelength=1.1e-6, photon_irradiance=1.0),
    ])
    @test_throws InvalidConfiguration prepare_direct_imaging(
        PupilFunction(tel), mixed_wavelengths)
end

@testset "Telescope and direct imaging" begin
    @test Base.ispublic(Optics, :AbstractTelescopeDefinition)
    @test Base.ispublic(Optics, :TelescopeDefinition)
    @test Base.ispublic(Optics, :prepare_telescope)
    @test !Base.isexported(Optics, :TelescopeDefinition)
    @test !Base.ispublic(AdaptiveOpticsSim, :TelescopeDefinition)

    definition = TelescopeDefinition(
        resolution=32,
        diameter=8.0,
        central_obstruction=0.2,
        fov_arcsec=10.0,
        pupil_reflectivity=0.25,
        revision=7,
        T=Float32,
    )
    @test definition isa AbstractTelescopeDefinition
    @test !ismutabletype(typeof(definition))
    @test fieldtypes(typeof(definition)) ==
        (Int, Float32, Float32, Float32, Float32, UInt)
    @test !hasproperty(definition, :backend)
    @test !hasproperty(definition, :device)
    @test all(!(getfield(definition, name) isa AbstractArray) for
        name in fieldnames(typeof(definition)))

    host_target = AdaptiveOpticsSim.Backends.HostComputeDevice()
    prepared_telescope = @inferred prepare_telescope(definition, host_target)
    standalone_telescope = Telescope(
        resolution=32,
        diameter=8.0,
        central_obstruction=0.2,
        fov_arcsec=10.0,
        pupil_reflectivity=0.25,
        T=Float32,
    )
    @test prepared_telescope isa Telescope
    @test backend(prepared_telescope) == CPUBackend()
    @test compute_device(pupil_mask(prepared_telescope)) == host_target
    @test compute_device(pupil_reflectivity(prepared_telescope)) == host_target
    @test aperture_revision(prepared_telescope) == UInt(7)
    @test pupil_mask(prepared_telescope) == pupil_mask(standalone_telescope)
    @test pupil_reflectivity(prepared_telescope) ==
        pupil_reflectivity(standalone_telescope)

    @test_throws UndefKeywordError TelescopeDefinition(
        resolution=8, diameter=8.0)
    for invalid_resolution in (0, -1)
        @test_throws InvalidConfiguration TelescopeDefinition(
            resolution=invalid_resolution, diameter=8.0, revision=0)
        @test_throws InvalidConfiguration Telescope(
            resolution=invalid_resolution, diameter=8.0)
    end
    for invalid_diameter in (0.0, -1.0, Inf, NaN)
        @test_throws InvalidConfiguration TelescopeDefinition(
            resolution=8, diameter=invalid_diameter, revision=0)
        @test_throws InvalidConfiguration Telescope(
            resolution=8, diameter=invalid_diameter)
    end
    for invalid_obstruction in (-0.1, 1.0, Inf, NaN)
        @test_throws InvalidConfiguration TelescopeDefinition(
            resolution=8, diameter=8.0,
            central_obstruction=invalid_obstruction, revision=0)
        @test_throws InvalidConfiguration Telescope(
            resolution=8, diameter=8.0,
            central_obstruction=invalid_obstruction)
    end
    for invalid_fov in (-0.1, Inf, NaN)
        @test_throws InvalidConfiguration TelescopeDefinition(
            resolution=8, diameter=8.0, fov_arcsec=invalid_fov,
            revision=0)
        @test_throws InvalidConfiguration Telescope(
            resolution=8, diameter=8.0, fov_arcsec=invalid_fov)
    end
    @test_throws InvalidConfiguration TelescopeDefinition(
        resolution=8, diameter=8.0, revision=-1)
    @test_throws InvalidConfiguration prepare_telescope(
        MutableTelescopeTestDefinition(), host_target)
    @test_throws InvalidConfiguration TelescopeDefinition(
        resolution=true, diameter=8.0, revision=0)
    @test_throws InvalidConfiguration TelescopeDefinition(
        resolution=8, diameter=8.0, revision=true)
    @test_throws InvalidConfiguration TelescopeDefinition(
        resolution=8, diameter=big(10.0)^1000, revision=0, T=Float32)
    @test_throws TypeError TelescopeDefinition(
        resolution=8, diameter=8.0, pupil_reflectivity=ones(8, 8),
        revision=0)

    unavailable_target = AdaptiveOpticsSim.Backends.AcceleratorComputeDevice(
        CUDABackend(), 0)
    unavailable_error = try
        prepare_telescope(definition, unavailable_target)
        nothing
    catch error
        error
    end
    @test unavailable_error isa AdaptiveOpticsSim.Backends.ComputeDeviceError
    @test unavailable_error.operation == :select
    @test unavailable_error.device == unavailable_target

    for invalid_reflectivity in (-0.1, 1.1, Inf, NaN)
        @test_throws InvalidConfiguration TelescopeDefinition(resolution=8,
            diameter=8.0, pupil_reflectivity=invalid_reflectivity,
            revision=0)
        @test_throws InvalidConfiguration Telescope(resolution=8,
            diameter=8.0, pupil_reflectivity=invalid_reflectivity)
    end
    invalid_reflectivity_map = ones(8, 8)
    invalid_reflectivity_map[1, 1] = NaN
    @test_throws InvalidConfiguration Telescope(resolution=8,
        diameter=8.0, pupil_reflectivity=invalid_reflectivity_map)

    tel = Telescope(resolution=32, diameter=8.0, central_obstruction=0.2)
    src = Source(band=:I, magnitude=0.0)
    wavefront = PupilFunction(tel)
    field = ElectricField(wavefront, src; zero_padding=2)
    formation = prepare_pupil_field(wavefront, src, field)
    fill_electric_field!(field, wavefront, formation)
    @test size(field.values) == (64, 64)
    @test field.metadata.coordinate_domain isa MetricCoordinates
    @test field.metadata.spectral == MonochromaticChannel(wavelength(src))
    fraunhofer = FraunhoferPropagation(field)
    @test Optics.propagation_output_metadata(fraunhofer).coordinate_domain isa
        AngularCoordinates
    centered_rate = similar(field.values, Float64)
    fraunhofer_intensity_from_field!(centered_rate, field, fraunhofer)
    direct = prepare_direct_imaging(wavefront, src;
        zero_padding=2)
    rate_map = form_direct_image!(direct)
    rate = intensity_values(rate_map)
    @test size(rate) == (64, 64)
    @test maximum(rate) > 0
    @test isfinite(sum(rate))
    @test centered_rate ≈ rate
    @test !hasproperty(tel, :state)

    tel_dim = Telescope(resolution=32, diameter=8.0, central_obstruction=0.2,
        pupil_reflectivity=0.25)
    dim_pupil = PupilFunction(tel_dim)
    dim_direct = prepare_direct_imaging(dim_pupil, src;
        zero_padding=2)
    rate_dim = intensity_values(form_direct_image!(dim_direct))
    @test sum(rate_dim) ≈ 0.25 * sum(rate)
    photon_rate = pupil_photon_rate_map(tel_dim, src)
    @test size(photon_rate) == size(pupil_mask(tel_dim))
    @test maximum(photon_rate) > 0
    @test optical_path(src, tel_dim) == "source(I) -> telescope"
    reflectivity_before_invalid_set = copy(pupil_reflectivity(tel_dim))
    @test_throws InvalidConfiguration set_pupil_reflectivity!(tel_dim, -0.1)
    @test pupil_reflectivity(tel_dim) == reflectivity_before_invalid_set

    tel_simple = Telescope(resolution=8, diameter=8.0, central_obstruction=0.0)
    src_simple = Source(band=:I, magnitude=0.0)
    wavefront_simple = PupilFunction(tel_simple)
    field_simple = ElectricField(wavefront_simple, src_simple;
        zero_padding=1)
    formation_simple = prepare_pupil_field(wavefront_simple,
        src_simple, field_simple)
    fill_electric_field!(field_simple, wavefront_simple, formation_simple)
    phase_map = fill(pi / 2, tel_simple.params.resolution, tel_simple.params.resolution)
    baseline = copy(field_simple.values)
    apply_phase!(field_simple, phase_map; units=:phase)
    @test field_simple.values ≈ baseline .* cis.(phase_map)

    fill_electric_field!(field_simple, wavefront_simple, formation_simple)
    quarter_wave_opd = fill(eltype(wavefront_simple.opd)(
        wavelength(src_simple) / 4), tel_simple.params.resolution,
        tel_simple.params.resolution)
    baseline = copy(field_simple.values)
    apply_phase!(field_simple, quarter_wave_opd; units=:opd)
    @test field_simple.values ≈ baseline .* cispi.(2 .* quarter_wave_opd ./ wavelength(src_simple))

    fill_electric_field!(field_simple, wavefront_simple, formation_simple)
    amplitude_map = fill(0.5, tel_simple.params.resolution, tel_simple.params.resolution)
    baseline = copy(field_simple.values)
    apply_amplitude!(field_simple, amplitude_map)
    @test field_simple.values ≈ baseline .* amplitude_map

    intensity_buffer = similar(field_simple.values, Float64)
    intensity!(intensity_buffer, field_simple)
    @test intensity_buffer ≈ abs2.(field_simple.values)

    propagated = propagation_output(field, fraunhofer)
    propagate_field!(propagated, field, fraunhofer)
    propagated_intensity = similar(propagated.values, Float64)
    @. propagated_intensity = abs2(propagated.values)
    @test propagated_intensity ≈ centered_rate atol=1e-10 rtol=1e-10
    @test Optics.propagation_plan(fraunhofer).params.output_sampling_rad ≈
        wavelength(src) /
        (field.metadata.dimensions[1] * field.metadata.sampling[1])
    @test propagated.metadata.kind isa FocalPlane

    fresnel = FresnelPropagation(field; distance_m=25.0)
    fresnel_out = propagation_output(field, fresnel)
    propagate_field!(fresnel_out, field, fresnel)
    @test size(fresnel_out.values) == size(field.values)
    @test isfinite(sum(abs, fresnel_out.values))
    reverse_fresnel = FresnelPropagation(fresnel_out; distance_m=-25.0,
        output_kind=PupilPlane())
    reverse_output = propagation_output(fresnel_out, reverse_fresnel)
    propagate_field!(reverse_output, fresnel_out, reverse_fresnel)
    @test reverse_output.values ≈ field.values atol=1e-9 rtol=1e-9

    zero_fresnel = FresnelPropagation(field; distance_m=0.0)
    zero_out = propagation_output(field, zero_fresnel)
    propagate_field!(zero_out, field, zero_fresnel)
    @test zero_out.values ≈ field.values atol=1e-10 rtol=1e-10

    atm_tel = Telescope(resolution=16, diameter=8.0, central_obstruction=0.0)
    atm_src = Source(band=:I, magnitude=0.0)
    atm = MultiLayerAtmosphere(atm_tel;
        r0=0.2,
        L0=25.0,
        fractional_cn2=[0.7, 0.3],
        wind_speed=[8.0, 4.0],
        wind_direction=[0.0, 90.0],
        altitude=[0.0, 5000.0],
    )
    advance_by!(atm, TEST_ATMOSPHERE_STEP; rng=MersenneTwister(1))
    atm_pupil = PupilFunction(atm_tel)
    geom_prop = AtmosphericFieldPropagation(atm, atm_pupil, atm_src;
        model=GeometricAtmosphericPropagation(T=Float64),
        zero_padding=1,
        T=Float64)
    @test AdaptiveOpticsSim.Atmospheres.atmospheric_field_execution_strategy(
        AdaptiveOpticsSim.Backends.execution_style(first(geom_prop.state.slices).field.values),
        geom_prop.params.model,
    ) isa AdaptiveOpticsSim.Atmospheres.GeometricFieldSynchronousStrategy
    geom_field = propagate_atmosphere_field!(geom_prop, atm,
        current_epoch(atm))
    tel_geom = Telescope(resolution=16, diameter=8.0, central_obstruction=0.0)
    collapsed_wavefront = PupilFunction(tel_geom)
    collapsed_renderer = prepare_atmosphere_renderer(atm, tel_geom, atm_src)
    render_atmosphere!(collapsed_wavefront, collapsed_renderer, atm,
        current_epoch(atm))
    collapsed = ElectricField(collapsed_wavefront, atm_src; zero_padding=1)
    collapsed_plan = prepare_pupil_field(collapsed_wavefront,
        atm_src, collapsed)
    fill_electric_field!(collapsed, collapsed_wavefront, collapsed_plan)
    @test geom_field.values ≈ collapsed.values atol=1e-8 rtol=1e-8

    fresnel_atm = MultiLayerAtmosphere(atm_tel;
        r0=0.2,
        L0=25.0,
        fractional_cn2=[1.0],
        wind_speed=[0.0],
        wind_direction=[0.0],
        altitude=[0.0],
    )
    advance_by!(fresnel_atm, TEST_ATMOSPHERE_STEP; rng=MersenneTwister(2))
    fresnel_prop = AtmosphericFieldPropagation(fresnel_atm, atm_pupil,
        atm_src;
        model=LayeredFresnelAtmosphericPropagation(T=Float64),
        zero_padding=1,
        T=Float64)
    @test AdaptiveOpticsSim.Atmospheres.atmospheric_field_execution_strategy(
        AdaptiveOpticsSim.Backends.execution_style(first(fresnel_prop.state.slices).field.values),
        fresnel_prop.params.model,
    ) isa AdaptiveOpticsSim.Atmospheres.LayeredFresnelFieldSynchronousStrategy
    fresnel_field = propagate_atmosphere_field!(fresnel_prop, fresnel_atm,
        current_epoch(fresnel_atm))
    geom_single = AtmosphericFieldPropagation(fresnel_atm, atm_pupil,
        atm_src;
        model=GeometricAtmosphericPropagation(T=Float64),
        zero_padding=1,
        T=Float64)
    geom_single_field = propagate_atmosphere_field!(geom_single, fresnel_atm,
        current_epoch(fresnel_atm))
    @test fresnel_field.values ≈ geom_single_field.values atol=1e-8 rtol=1e-8

    spectral = with_spectrum(atm_src, SpectralBundle([wavelength(atm_src), 1.1 * wavelength(atm_src)], [0.75, 0.25]))
    spectral_prop = AtmosphericFieldPropagation(atm, atm_pupil, spectral;
        model=GeometricAtmosphericPropagation(T=Float64),
        zero_padding=2,
        T=Float64)
    spectral_intensity = atmospheric_intensity!(spectral_prop, atm,
        current_epoch(atm))
    mono_prop = AtmosphericFieldPropagation(atm, atm_pupil, atm_src;
        model=GeometricAtmosphericPropagation(T=Float64),
        zero_padding=2,
        T=Float64)
    mono_intensity = atmospheric_intensity!(mono_prop, atm,
        current_epoch(atm))
    single_spectral = with_spectrum(atm_src, SpectralBundle([wavelength(atm_src)], [1.0]))
    single_spectral_prop = AtmosphericFieldPropagation(atm, atm_pupil,
        single_spectral;
        model=GeometricAtmosphericPropagation(T=Float64),
        zero_padding=2,
        T=Float64)
    @test atmospheric_intensity!(single_spectral_prop, atm,
        current_epoch(atm)) ≈ mono_intensity atol=1e-8 rtol=1e-8
    @test size(spectral_intensity) == size(mono_intensity)
    @test sum(spectral_intensity) > 0

    normalized_atm_src = Source(band=:custom,
        wavelength=wavelength(atm_src), normalized_power=2.0)
    normalized_spectral = with_spectrum(normalized_atm_src,
        SpectralBundle([wavelength(atm_src), 1.1 * wavelength(atm_src)],
            [0.75, 0.25]))
    normalized_spectral_prop = AtmosphericFieldPropagation(atm, atm_pupil,
        normalized_spectral;
        model=GeometricAtmosphericPropagation(T=Float64),
        zero_padding=1,
        T=Float64)
    @test all(slice.field.metadata.normalization isa
        DimensionlessNormalization for slice in
        normalized_spectral_prop.state.slices)
    normalized_spectral_intensity = atmospheric_intensity!(
        normalized_spectral_prop, atm, current_epoch(atm))
    @test sum(normalized_spectral_intensity) ≈ 2.0 rtol=1e-10
end

@testset "Direct-imaging exact target ownership" begin
    optics = AdaptiveOpticsSim.Optics
    backends = AdaptiveOpticsSim.Backends
    target = backends.HostComputeDevice()
    telescope = Telescope(
        resolution=8,
        diameter=8.0,
        central_obstruction=0.0,
    )
    pupil = PupilFunction(telescope)
    first_source = Source(
        band=:custom,
        wavelength=1.0e-6,
        photon_irradiance=1.0,
    )
    second_source = Source(
        band=:custom,
        wavelength=1.0e-6,
        photon_irradiance=2.0,
    )

    leaf = prepare_direct_imaging(pupil, first_source)
    @test @inferred(optics._require_exact_direct_imaging_target(
        leaf.plan, target)) === leaf.plan
    @test @inferred(optics._require_exact_direct_imaging_target(
        leaf, target)) === leaf

    preformed_field = ElectricField(pupil, first_source)
    preformed = prepare_direct_imaging(first_source, preformed_field)
    @test @inferred(optics._require_exact_direct_imaging_target(
        preformed, target)) === preformed

    asterism = Asterism([first_source, second_source])
    incoherent = prepare_direct_imaging(pupil, asterism)
    @test @inferred(optics._require_exact_direct_imaging_target(
        incoherent, target)) === incoherent

    spectral = with_spectrum(first_source, SpectralBundle(
        [1.0e-6, 1.1e-6],
        [0.5, 0.5],
    ))
    bundled = prepare_direct_imaging(pupil, spectral)
    @test @inferred(optics._require_exact_direct_imaging_target(
        bundled, target)) === bundled

    batch = prepare_direct_imaging_batch(pupil, asterism)
    @test @inferred(optics._require_exact_direct_imaging_target(
        batch, target)) === batch

    replacement_workspace = optics.DirectImagingWorkspace(
        optics.FraunhoferPropagationWorkspace(
            leaf.plan.propagation,
            leaf.field.values,
        ),
        copy(leaf.workspace.unshifted_intensity),
    )
    replacement_owner = optics.PreparedDirectImaging(
        optics._PREPARED_DIRECT_IMAGING_TOKEN,
        leaf.input,
        leaf.field,
        leaf.output,
        leaf.plan,
        replacement_workspace,
    )
    @test @inferred(optics._require_exact_direct_imaging_target(
        replacement_owner, target)) === replacement_owner

    wrong_target = backends.AcceleratorComputeDevice(CUDABackend(), 0)
    wrong_target_error = try
        optics._require_exact_direct_imaging_target(leaf, wrong_target)
        nothing
    catch error
        error
    end
    @test wrong_target_error isa backends.ComputeDeviceError
    @test wrong_target_error.operation == :validate_direct_imaging_target
    @test wrong_target_error.reason == :wrong_device
    @test wrong_target_error.device == wrong_target
    @test_throws InvalidConfiguration optics._require_exact_direct_imaging_target(
        nothing, target)
    @test !Base.ispublic(optics, :_require_exact_direct_imaging_target)
end

@testset "Explicit optical products and surfaces" begin
    tel = Telescope(resolution=16, diameter=8.0, central_obstruction=0.1)
    src = Source(band=:I, magnitude=0.0)
    aperture_before = copy(pupil_mask(tel))
    reflectivity_before = copy(pupil_reflectivity(tel))

    path_a = PupilFunction(tel)
    path_b = PupilFunction(tel)
    static_map = OPDMap(fill(2e-9, 16, 16))
    apply_surface!(path_a, static_map, DMAdditive())
    @test path_a.opd == static_map.opd
    @test all(iszero, path_b.opd)

    ncpa = NCPA(fill(3e-9, 16, 16))
    apply_surface!(path_b, ncpa, DMReplace())
    @test path_b.opd == ncpa.opd

    dm = DeformableMirror(tel; n_act=4, influence_width=0.3)
    set_command!(dm, fill(1e-8, n_actuators(dm)))
    update_surface!(dm)
    apply_surface!(path_a, dm, DMReplace())
    @test path_a.opd == surface_opd(dm)

    modal = TipTiltMirror(tel; scale=1e-8)
    set_command!(modal, [0.25, -0.5])
    update_surface!(modal)
    reset_opd!(path_b)
    apply_surface!(path_b, modal, DMAdditive())
    @test path_b.opd == surface_opd(modal)
    @test pupil_mask(tel) == aperture_before
    @test pupil_reflectivity(tel) == reflectivity_before

    reset_opd!(path_a)
    path_b_before_propagation = copy(path_b.opd)
    field = ElectricField(path_a, src; zero_padding=2)
    propagation = FraunhoferPropagation(field)
    output = IntensityMap(field, propagation)
    fill!(output.values, 0)
    prepared = prepare_direct_imaging(path_a, src, field, output)
    explicit_direct_image_cycle!(prepared)
    @test @allocated(explicit_direct_image_cycle!(prepared)) == 0
    @test path_b.opd == path_b_before_propagation
    @test !hasproperty(tel, :state)

    spatial_filter = SpatialFilter(tel; shape=SquareFilter(), diameter=5,
        zero_padding=2)
    spatial_field = ElectricField(path_a, src; zero_padding=2,
        normalization=DimensionlessNormalization(),
        spatial_measure=PointSampledMeasure(),
        coherence=CoherentFieldCombination())
    spatial_formation = prepare_pupil_field(path_a, src, spatial_field;
        center_even_grid=false, amplitude_scale=1)
    fill_electric_field!(spatial_field, path_a, spatial_formation)
    spatial_output = PupilFunction(tel)
    prepared_spatial_filter = prepare_spatial_filter(
        tel, spatial_filter, spatial_field,
        spatial_output)
    explicit_spatial_filter_cycle!(prepared_spatial_filter)
    @test @allocated(explicit_spatial_filter_cycle!(
        prepared_spatial_filter)) == 0
end

@testset "Optical-plane compatibility validation" begin
    tel = Telescope(resolution=8, diameter=8.0, central_obstruction=0.0)
    src = Source(band=:I, magnitude=0.0)
    wavefront = PupilFunction(tel)
    field = ElectricField(wavefront, src; zero_padding=1)
    values = similar(field.values)

    sampling_metadata = OpticalPlaneMetadata(PupilPlane(), values;
        coordinate_domain=field.metadata.coordinate_domain,
        sampling=(2.0, 2.0), spectral=field.metadata.spectral)
    sampling_field = ElectricField(sampling_metadata, values)
    @test_throws InvalidConfiguration prepare_pupil_field(wavefront,
        src, sampling_field)

    coordinate_metadata = OpticalPlaneMetadata(PupilPlane(), values;
        coordinate_domain=AngularCoordinates(),
        sampling=field.metadata.sampling, origin=field.metadata.origin,
        spectral=field.metadata.spectral)
    coordinate_field = ElectricField(coordinate_metadata, values)
    @test_throws InvalidConfiguration prepare_pupil_field(wavefront,
        src, coordinate_field)
    @test_throws InvalidConfiguration FraunhoferPropagation(coordinate_field)
    spatial_filter = SpatialFilter(tel; zero_padding=1)
    @test_throws InvalidConfiguration prepare_spatial_filter(tel,
        spatial_filter, coordinate_field, wavefront)

    kind_metadata = OpticalPlaneMetadata(FocalPlane(), values;
        coordinate_domain=field.metadata.coordinate_domain,
        sampling=field.metadata.sampling, origin=field.metadata.origin,
        spectral=field.metadata.spectral)
    kind_field = ElectricField(kind_metadata, values)
    @test_throws InvalidConfiguration prepare_pupil_field(wavefront,
        src, kind_field)

    origin_metadata = OpticalPlaneMetadata(PupilPlane(), values;
        coordinate_domain=field.metadata.coordinate_domain,
        sampling=field.metadata.sampling, origin=(0.0, 0.0),
        spectral=field.metadata.spectral)
    origin_field = ElectricField(origin_metadata, values)
    @test_throws InvalidConfiguration prepare_pupil_field(wavefront,
        src, origin_field)

    centering_metadata = OpticalPlaneMetadata(PupilPlane(), values;
        coordinate_domain=field.metadata.coordinate_domain,
        sampling=field.metadata.sampling, origin=field.metadata.origin,
        centering=(SampleCentered, SampleCentered),
        spectral=field.metadata.spectral)
    centering_field = ElectricField(centering_metadata, values)
    @test_throws InvalidConfiguration prepare_pupil_field(wavefront,
        src, centering_field)

    orientation_metadata = OpticalPlaneMetadata(PupilPlane(), values;
        coordinate_domain=field.metadata.coordinate_domain,
        sampling=field.metadata.sampling, origin=field.metadata.origin,
        orientation=PlaneAxisOrientation((:y, :x)),
        spectral=field.metadata.spectral)
    orientation_field = ElectricField(orientation_metadata, values)
    @test_throws InvalidConfiguration prepare_pupil_field(wavefront,
        src, orientation_field)

    spectral_metadata = OpticalPlaneMetadata(PupilPlane(), values;
        coordinate_domain=field.metadata.coordinate_domain,
        sampling=field.metadata.sampling, origin=field.metadata.origin,
        spectral=MonochromaticChannel(1.1 * wavelength(src)))
    spectral_field = ElectricField(spectral_metadata, values)
    @test_throws InvalidConfiguration prepare_pupil_field(wavefront,
        src, spectral_field)

    float32_values = Matrix{ComplexF32}(undef, size(values))
    numeric_metadata = OpticalPlaneMetadata(PupilPlane(), float32_values;
        coordinate_domain=field.metadata.coordinate_domain,
        sampling=field.metadata.sampling, origin=field.metadata.origin,
        spectral=field.metadata.spectral)
    numeric_field = ElectricField(numeric_metadata, float32_values)
    @test_throws InvalidConfiguration prepare_pupil_field(wavefront,
        src, numeric_field)

    wrong_size_values = Matrix{ComplexF64}(undef, 9, 9)
    dimension_metadata = OpticalPlaneMetadata(PupilPlane(),
        wrong_size_values; coordinate_domain=field.metadata.coordinate_domain,
        sampling=field.metadata.sampling,
        spectral=field.metadata.spectral)
    dimension_field = ElectricField(dimension_metadata, wrong_size_values)
    @test_throws DimensionMismatchError prepare_pupil_field(wavefront,
        src, dimension_field)

    propagation = FraunhoferPropagation(field)
    wrong_destination_values = Matrix{Float64}(undef, 9, 9)
    wrong_destination_metadata = OpticalPlaneMetadata(FocalPlane(),
        wrong_destination_values;
        coordinate_domain=Optics.propagation_output_metadata(
            propagation).coordinate_domain,
        sampling=Optics.propagation_output_metadata(propagation).sampling,
        spectral=Optics.propagation_output_metadata(propagation).spectral)
    wrong_destination = IntensityMap(wrong_destination_metadata,
        wrong_destination_values)
    @test_throws DimensionMismatchError prepare_direct_imaging(wavefront,
        src, field, wrong_destination)

    declared_device_metadata = OpticalPlaneMetadata(PupilPlane(), values;
        coordinate_domain=field.metadata.coordinate_domain,
        sampling=field.metadata.sampling, origin=field.metadata.origin,
        spectral=field.metadata.spectral,
        device=AdaptiveOpticsSim.Backends.AcceleratorComputeDevice(
            CUDABackend(), 1))
    @test_throws InvalidConfiguration ElectricField(
        declared_device_metadata, values)
end

@testset "Aperture masks" begin
    bool_mask = falses(32, 32)
    build_mask!(bool_mask, CircularAperture(radius=1.0))
    @test bool_mask[16, 16]
    @test !bool_mask[1, 1]

    annulus = falses(32, 32)
    build_mask!(annulus, AnnularAperture(inner_radius=0.25, outer_radius=1.0))
    @test !annulus[16, 16]
    @test count(annulus) < count(bool_mask)

    weighted = fill(-1.0, 16, 16)
    build_mask!(weighted, RectangularROI(5:8, 6:10); inside=2.0, outside=-1.0)
    @test all(weighted[5:8, 6:10] .== 2.0)
    @test weighted[4, 6] == -1.0

    spider = trues(32, 32)
    apply_mask!(spider, SpiderMask(thickness=0.08, angle_rad=pi / 2))
    @test count(spider) < length(spider)

    tel = Telescope(resolution=32, diameter=8.0, central_obstruction=0.25)
    expected = falses(32, 32)
    build_mask!(expected, AnnularAperture(inner_radius=0.25, outer_radius=1.0))
    @test pupil_mask(tel) == expected
    apply_spiders!(tel; thickness=0.4, angles=[0.0, 90.0])
    manual = copy(expected)
    apply_mask!(manual, SpiderMask(thickness=0.1, angle_rad=0.0))
    apply_mask!(manual, SpiderMask(thickness=0.1, angle_rad=pi / 2))
    @test pupil_mask(tel) == manual

    sf_tel = Telescope(resolution=8, diameter=8.0, central_obstruction=0.0)
    sf = SpatialFilter(sf_tel; shape=SquareFilter(), diameter=4, zero_padding=2)
    @test count(x -> !iszero(x), sf.mask) > 0
    foucault = SpatialFilter(sf_tel; shape=FoucaultFilter(), diameter=4, zero_padding=2)
    foucault_count = count(x -> !iszero(x), foucault.mask)
    @test 0 < foucault_count < length(foucault.mask)

    pupil = falses(8, 8)
    pupil[1:4, 1:4] .= true
    pupil[5:8, 5:8] .= true
    valid = falses(2, 2)
    build_mask!(valid, SubapertureGridMask(threshold=0.5), pupil)
    @test valid == Bool[true false; false true]
    valid2 = similar(valid)
    AdaptiveOpticsSim.WavefrontSensors.set_valid_subapertures!(
        valid2, pupil, 0.5)
    @test valid2 == valid
end

@testset "Zernike basis" begin
    tel = Telescope(resolution=32, diameter=8.0, central_obstruction=0.0)
    zb = ZernikeBasis(tel, 5)
    compute_zernike!(zb, tel)
    @test size(zb.modes) == (32, 32, 5)
    @test noll_to_nm(1) == (0, 0)
    @test noll_to_nm(2) == (1, -1)
    @test noll_to_nm(3) == (1, 1)
    @test noll_to_nm(4) == (2, -2)
end
