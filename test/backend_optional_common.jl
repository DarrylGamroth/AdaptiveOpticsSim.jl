backend_package_name(::Type{CUDABackendTag}) = "CUDA"
backend_package_name(::Type{AMDGPUBackendTag}) = "AMDGPU"

backend_label(::Type{CUDABackendTag}) = "CUDA"
backend_label(::Type{AMDGPUBackendTag}) = "AMDGPU"

backend_full_smoke_env(::Type{CUDABackendTag}) = "ADAPTIVEOPTICS_TEST_FULL_CUDA"
backend_full_smoke_env(::Type{AMDGPUBackendTag}) = "ADAPTIVEOPTICS_TEST_FULL_AMDGPU"

function build_optional_platform_branch(::Type{T}, BackendArray, label::Symbol; sensor::Symbol=:sh, seed::Integer=1) where {T<:AbstractFloat}
    tel = Telescope(resolution=16, diameter=T(8.0), sampling_time=T(1e-3),
        central_obstruction=T(0.0), T=T, backend=BackendArray)
    src = Source(band=:I, magnitude=0.0, T=T)
    atm = KolmogorovAtmosphere(tel; r0=T(0.2), L0=T(25.0), T=T, backend=BackendArray)
    dm = DeformableMirror(tel; n_act=4, influence_width=T(0.3), T=T, backend=BackendArray)
    wfs = sensor == :sh ?
        ShackHartmann(tel; n_subap=4, mode=Diffractive(), T=T, backend=BackendArray) :
        PyramidWFS(tel; n_subap=4, modulation=T(1.0), mode=Diffractive(), T=T, backend=BackendArray)
    det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0), binning=1, T=T, backend=BackendArray)
    sim = AdaptiveOpticsSim.AOSimulation(tel, atm, src, dm, wfs)
    imat = interaction_matrix(dm, wfs, tel, src; amplitude=T(0.05))
    recon = ModalReconstructor(imat; gain=T(0.5))
    return ClosedLoopBranchConfig(label, sim, recon; wfs_detector=det, rng=MersenneTwister(seed))
end

function import_backend_package!(::Type{CUDABackendTag})
    @eval import CUDA
    return nothing
end

function import_backend_package!(::Type{AMDGPUBackendTag})
    @eval import AMDGPU
    return nothing
end

backend_functional(::Type{CUDABackendTag}) = CUDA.functional()
backend_functional(::Type{AMDGPUBackendTag}) = AMDGPU.functional()

run_optional_backend_plan_checks(::Type{<:GPUBackendTag}, tel, backend) = nothing

function run_optional_backend_plan_checks(::Type{AMDGPUBackendTag}, tel, backend)
    T = Float32
    sh = ShackHartmann(tel; n_subap=4, mode=Diffractive(), T=T, backend=backend)
    pyr = PyramidWFS(tel; n_subap=4, modulation=T(1.0), mode=Diffractive(), T=T, backend=backend)
    bio = BioEdgeWFS(tel; n_subap=4, modulation=T(1.0), mode=Diffractive(), T=T, backend=backend)
    det = Detector(noise=NoiseReadout(T(1.0)), qe=1.0, sensor=HgCdTeAvalancheArraySensor(T=T), T=T, backend=backend)
    det_capture = Detector(noise=NoiseReadout(T(1.0)), qe=1.0, bits=12, full_well=T(100),
        sensor=CMOSSensor(T=T), T=T, backend=backend)
    src = Source(band=:I, magnitude=0.0, T=T)
    atm = MultiLayerAtmosphere(tel;
        r0=T(0.2),
        L0=T(25.0),
        fractional_cn2=T[1.0],
        wind_speed=T[0.0],
        wind_direction=T[0.0],
        altitude=T[0.0],
        T=T,
        backend=backend,
    )
    geom_prop = AtmosphericFieldPropagation(atm, tel, src;
        model=GeometricAtmosphericPropagation(T=T),
        zero_padding=1,
        T=T)
    fresnel_prop = AtmosphericFieldPropagation(atm, tel, src;
        model=LayeredFresnelAtmosphericPropagation(T=T),
        zero_padding=1,
        T=T)
    @test AdaptiveOpticsSim.grouped_accumulation_plan(AdaptiveOpticsSim.execution_style(pyr.state.intensity), pyr) isa AdaptiveOpticsSim.GroupedStaged2DPlan
    @test AdaptiveOpticsSim.grouped_accumulation_plan(AdaptiveOpticsSim.execution_style(bio.state.intensity), bio) isa AdaptiveOpticsSim.GroupedStaged2DPlan
    @test AdaptiveOpticsSim.sh_sensing_execution_plan(AdaptiveOpticsSim.execution_style(sh.state.slopes), sh) isa AdaptiveOpticsSim.ShackHartmannBatchedPlan
    @test AdaptiveOpticsSim.detector_execution_plan(typeof(AdaptiveOpticsSim.execution_style(det.state.frame)), typeof(det)) isa AdaptiveOpticsSim.DetectorHostMirrorPlan
    capture_psf = backend{T}(undef, 4, 4)
    fill!(capture_psf, T(10))
    captured = capture!(det_capture, capture_psf; rng=MersenneTwister(2))
    @test captured isa backend
    @test maximum(Array(captured)) <= Float64(exp2(T(12)) - one(T))
    @test AdaptiveOpticsSim.reduction_execution_plan(pyr.state.intensity) isa AdaptiveOpticsSim.HostMirrorReductionPlan
    randn_method = which(
        AdaptiveOpticsSim._randn_frame_noise!,
        (
            AdaptiveOpticsSim.DetectorHostMirrorPlan,
            typeof(det),
            typeof(MersenneTwister(1)),
            typeof(det.state.frame),
        ),
    )
    @test occursin("AdaptiveOpticsSimAMDGPUExt", String(randn_method.file))
    @test AdaptiveOpticsSim.atmospheric_field_execution_plan(
        AdaptiveOpticsSim.execution_style(first(geom_prop.state.slices).field.state.field),
        geom_prop.params.model,
    ) isa AdaptiveOpticsSim.GeometricFieldAsyncPlan
    @test AdaptiveOpticsSim.atmospheric_field_execution_plan(
        AdaptiveOpticsSim.execution_style(first(fresnel_prop.state.slices).field.state.field),
        fresnel_prop.params.model,
    ) isa AdaptiveOpticsSim.LayeredFresnelFieldAsyncPlan
    cpu_tel = Telescope(resolution=16, diameter=8.0f0, sampling_time=1.0f-3,
        central_obstruction=0.0f0, T=T, backend=Array)
    gpu_tel = Telescope(resolution=16, diameter=8.0f0, sampling_time=1.0f-3,
        central_obstruction=0.0f0, T=T, backend=backend)
    cpu_src = Source(band=:I, magnitude=0.0, T=T)
    gpu_src = Source(band=:I, magnitude=0.0, T=T)
    cpu_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
        sensor=CMOSSensor(T=T), response_model=NullFrameResponse(), T=T, backend=Array)
    gpu_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
        sensor=CMOSSensor(T=T), response_model=NullFrameResponse(), T=T, backend=backend)
    cpu_sh_stats = ShackHartmann(cpu_tel; n_subap=4, mode=Diffractive(), T=T, backend=Array,
        valid_subaperture_policy=FluxThresholdValidSubapertures(light_ratio=0.5f0))
    gpu_sh_stats = ShackHartmann(gpu_tel; n_subap=4, mode=Diffractive(), T=T, backend=backend,
        valid_subaperture_policy=FluxThresholdValidSubapertures(light_ratio=0.5f0))
    measure!(cpu_sh_stats, cpu_tel, cpu_src, cpu_det; rng=MersenneTwister(3))
    measure!(gpu_sh_stats, gpu_tel, gpu_src, gpu_det; rng=MersenneTwister(3))
    cpu_peak = AdaptiveOpticsSim.sh_safe_peak_value(cpu_sh_stats.state.spot_cube)
    cpu_cutoff = AdaptiveOpticsSim.centroid_threshold(cpu_sh_stats) * cpu_peak
    AdaptiveOpticsSim.sh_signal_from_spots!(cpu_sh_stats, cpu_cutoff)
    gpu_peak = AdaptiveOpticsSim.sh_safe_peak_value(gpu_sh_stats.state.spot_cube)
    gpu_cutoff = AdaptiveOpticsSim.centroid_threshold(gpu_sh_stats) * gpu_peak
    AdaptiveOpticsSim.sh_signal_from_spots_device_stats!(
        AdaptiveOpticsSim.execution_style(gpu_sh_stats.state.slopes),
        gpu_sh_stats,
        gpu_cutoff,
    )
    @test isapprox(Array(gpu_sh_stats.state.slopes), cpu_sh_stats.state.slopes; rtol=1f-5, atol=1f-4)

    correction_models = (
        ReferencePixelCommonModeCorrection(1, 1),
        ReferenceRowCommonModeCorrection(1),
        ReferenceColumnCommonModeCorrection(1),
        ReferenceOutputCommonModeCorrection(2; edge_rows=1, edge_cols=1),
        CompositeFrameReadoutCorrection((
            ReferenceRowCommonModeCorrection(1),
            ReferenceColumnCommonModeCorrection(1),
        )),
    )
    corr_input = reshape(T.(1:96), 2, 6, 8)
    cpu_quant_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
        bits=8, full_well=T(100), sensor=CMOSSensor(T=T),
        response_model=NullFrameResponse(), T=T, backend=Array)
    gpu_quant_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
        bits=8, full_well=T(100), sensor=CMOSSensor(T=T),
        response_model=NullFrameResponse(), T=T, backend=backend)
    quant_input = reshape(T.(1:48), 6, 8) .* T(3)
    cpu_quant_frame = capture!(cpu_quant_det, copy(quant_input); rng=MersenneTwister(13))
    gpu_quant_frame = capture!(gpu_quant_det, backend(copy(quant_input)); rng=MersenneTwister(13))
    @test gpu_quant_frame isa backend
    @test isapprox(Array(gpu_quant_frame), cpu_quant_frame; rtol=1f-5, atol=1f-4)

    for correction_model in correction_models
        cpu_frame_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
            sensor=HgCdTeAvalancheArraySensor(T=T),
            response_model=NullFrameResponse(),
            correction_model=correction_model,
            T=T,
            backend=Array)
        gpu_frame_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
            sensor=HgCdTeAvalancheArraySensor(T=T),
            response_model=NullFrameResponse(),
            correction_model=correction_model,
            T=T,
            backend=backend)
        frame_input = reshape(T.(1:48), 6, 8)
        cpu_frame = capture!(cpu_frame_det, frame_input; rng=MersenneTwister(12))
        gpu_frame = capture!(gpu_frame_det, backend(copy(frame_input)); rng=MersenneTwister(12))
        @test gpu_frame isa backend
        @test isapprox(Array(gpu_frame), cpu_frame; rtol=1f-5, atol=1f-4)

        cpu_corr_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
            sensor=HgCdTeAvalancheArraySensor(T=T),
            response_model=NullFrameResponse(),
            correction_model=correction_model,
            T=T,
            backend=Array)
        gpu_corr_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
            sensor=HgCdTeAvalancheArraySensor(T=T),
            response_model=NullFrameResponse(),
            correction_model=correction_model,
            T=T,
            backend=backend)
        cpu_corr_cube = copy(corr_input)
        gpu_corr_cube = backend(copy(corr_input))
        cpu_corr_scratch = similar(cpu_corr_cube)
        gpu_corr_scratch = similar(gpu_corr_cube)
        AdaptiveOpticsSim.capture_stack!(cpu_corr_det, cpu_corr_cube, cpu_corr_scratch; rng=MersenneTwister(11))
        AdaptiveOpticsSim.capture_stack!(gpu_corr_det, gpu_corr_cube, gpu_corr_scratch; rng=MersenneTwister(11))
        @test gpu_corr_cube isa backend
        @test isapprox(Array(gpu_corr_cube), cpu_corr_cube; rtol=1f-5, atol=1f-4)
    end

    generalized_correction = CompositeFrameReadoutCorrection((
        ReferenceRowCommonModeCorrection(1),
        ReferenceColumnCommonModeCorrection(1),
    ))
    cpu_gen_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
        bits=8, full_well=T(100), readout_window=FrameWindow(2:5, 3:7), output_precision=UInt16,
        sensor=HgCdTeAvalancheArraySensor(T=T),
        response_model=NullFrameResponse(),
        correction_model=generalized_correction,
        T=T,
        backend=Array)
    gpu_gen_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
        bits=8, full_well=T(100), readout_window=FrameWindow(2:5, 3:7), output_precision=UInt16,
        sensor=HgCdTeAvalancheArraySensor(T=T),
        response_model=NullFrameResponse(),
        correction_model=generalized_correction,
        T=T,
        backend=backend)
    cpu_gen_input = copy(corr_input)
    gpu_gen_input = backend(copy(corr_input))
    cpu_gen_out = Array{UInt16}(undef, 2, 4, 5)
    gpu_gen_out = backend{UInt16}(undef, 2, 4, 5)
    AdaptiveOpticsSim.capture_stack!(cpu_gen_det, cpu_gen_out, cpu_gen_input; rng=MersenneTwister(14))
    AdaptiveOpticsSim.capture_stack!(gpu_gen_det, gpu_gen_out, gpu_gen_input; rng=MersenneTwister(14))
    @test gpu_gen_out isa backend
    @test Array(gpu_gen_out) == cpu_gen_out

    cpu_windowed_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0), binning=1,
        gain=T(1.0),
        sensor=HgCdTeAvalancheArraySensor(T=T, sampling_mode=CorrelatedDoubleSampling()),
        response_model=NullFrameResponse(),
        readout_window=FrameWindow(2:3, 2:3),
        T=T,
        backend=Array)
    gpu_windowed_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0), binning=1,
        gain=T(1.0),
        sensor=HgCdTeAvalancheArraySensor(T=T, sampling_mode=CorrelatedDoubleSampling()),
        response_model=NullFrameResponse(),
        readout_window=FrameWindow(2:3, 2:3),
        T=T,
        backend=backend)
    windowed_input = fill(T(10), 4, 4)
    cpu_windowed_frame = capture!(cpu_windowed_det, windowed_input; rng=MersenneTwister(15))
    gpu_windowed_frame = capture!(gpu_windowed_det, backend(copy(windowed_input)); rng=MersenneTwister(15))
    @test isapprox(Array(gpu_windowed_frame), cpu_windowed_frame; rtol=1f-5, atol=1f-4)
    @test isapprox(Array(detector_reference_frame(gpu_windowed_det)), detector_reference_frame(cpu_windowed_det); rtol=1f-5, atol=1f-4)
    @test isapprox(Array(detector_signal_frame(gpu_windowed_det)), detector_signal_frame(cpu_windowed_det); rtol=1f-5, atol=1f-4)
    @test isapprox(Array(detector_combined_frame(gpu_windowed_det)), detector_combined_frame(cpu_windowed_det); rtol=1f-5, atol=1f-4)
    @test isapprox(Array(detector_reference_cube(gpu_windowed_det)), detector_reference_cube(cpu_windowed_det); rtol=1f-5, atol=1f-4)
    @test isapprox(Array(detector_signal_cube(gpu_windowed_det)), detector_signal_cube(cpu_windowed_det); rtol=1f-5, atol=1f-4)
    @test isapprox(Array(detector_read_cube(gpu_windowed_det)), detector_read_cube(cpu_windowed_det); rtol=1f-5, atol=1f-4)
    @test detector_read_times(gpu_windowed_det) == detector_read_times(cpu_windowed_det)
    return nothing
end

function run_optional_backend_plan_checks(::Type{CUDABackendTag}, tel, backend)
    T = Float32
    sh = ShackHartmann(tel; n_subap=4, mode=Diffractive(), T=T, backend=backend)
    pyr = PyramidWFS(tel; n_subap=4, modulation=T(1.0), mode=Diffractive(), T=T, backend=backend)
    bio = BioEdgeWFS(tel; n_subap=4, modulation=T(1.0), mode=Diffractive(), T=T, backend=backend)
    det = Detector(noise=NoiseReadout(T(1.0)), qe=1.0, sensor=HgCdTeAvalancheArraySensor(T=T), T=T, backend=backend)
    src = Source(band=:I, magnitude=0.0, T=T)
    atm = MultiLayerAtmosphere(tel;
        r0=T(0.2),
        L0=T(25.0),
        fractional_cn2=T[1.0],
        wind_speed=T[0.0],
        wind_direction=T[0.0],
        altitude=T[0.0],
        T=T,
        backend=backend,
    )
    geom_prop = AtmosphericFieldPropagation(atm, tel, src;
        model=GeometricAtmosphericPropagation(T=T),
        zero_padding=1,
        T=T)
    fresnel_prop = AtmosphericFieldPropagation(atm, tel, src;
        model=LayeredFresnelAtmosphericPropagation(T=T),
        zero_padding=1,
        T=T)
    @test AdaptiveOpticsSim.grouped_accumulation_plan(AdaptiveOpticsSim.execution_style(pyr.state.intensity), pyr) isa AdaptiveOpticsSim.GroupedStackReducePlan
    @test AdaptiveOpticsSim.grouped_accumulation_plan(AdaptiveOpticsSim.execution_style(bio.state.intensity), bio) isa AdaptiveOpticsSim.GroupedStackReducePlan
    @test AdaptiveOpticsSim.sh_sensing_execution_plan(AdaptiveOpticsSim.execution_style(sh.state.slopes), sh) isa AdaptiveOpticsSim.ShackHartmannBatchedPlan
    @test AdaptiveOpticsSim.detector_execution_plan(typeof(AdaptiveOpticsSim.execution_style(det.state.frame)), typeof(det)) isa AdaptiveOpticsSim.DetectorDirectPlan
    @test AdaptiveOpticsSim.reduction_execution_plan(pyr.state.intensity) isa AdaptiveOpticsSim.DirectReductionPlan
    @test AdaptiveOpticsSim.atmospheric_field_execution_plan(
        AdaptiveOpticsSim.execution_style(first(geom_prop.state.slices).field.state.field),
        geom_prop.params.model,
    ) isa AdaptiveOpticsSim.GeometricFieldAsyncPlan
    @test AdaptiveOpticsSim.atmospheric_field_execution_plan(
        AdaptiveOpticsSim.execution_style(first(fresnel_prop.state.slices).field.state.field),
        fresnel_prop.params.model,
    ) isa AdaptiveOpticsSim.LayeredFresnelFieldAsyncPlan
    cpu_tel = Telescope(resolution=16, diameter=8.0f0, sampling_time=1.0f-3,
        central_obstruction=0.0f0, T=T, backend=Array)
    gpu_tel = Telescope(resolution=16, diameter=8.0f0, sampling_time=1.0f-3,
        central_obstruction=0.0f0, T=T, backend=backend)
    cpu_src = Source(band=:I, magnitude=0.0, T=T)
    gpu_src = Source(band=:I, magnitude=0.0, T=T)
    cpu_sh = ShackHartmann(cpu_tel; n_subap=4, mode=Diffractive(), T=T, backend=Array)
    gpu_sh = ShackHartmann(gpu_tel; n_subap=4, mode=Diffractive(), T=T, backend=backend)
    cpu_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
        sensor=CMOSSensor(T=T), response_model=NullFrameResponse(), T=T, backend=Array)
    gpu_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
        sensor=CMOSSensor(T=T), response_model=NullFrameResponse(), T=T, backend=backend)
    measure!(cpu_sh, cpu_tel, cpu_src, cpu_det; rng=MersenneTwister(3))
    measure!(gpu_sh, gpu_tel, gpu_src, gpu_det; rng=MersenneTwister(3))
    cpu_export = Array(AdaptiveOpticsSim.sh_exported_spot_cube(cpu_sh))
    gpu_export = Array(AdaptiveOpticsSim.sh_exported_spot_cube(gpu_sh))
    cpu_frame = Array(AdaptiveOpticsSim.wfs_output_frame(cpu_sh, cpu_det))
    gpu_frame = Array(AdaptiveOpticsSim.wfs_output_frame(gpu_sh, gpu_det))
    @test size(gpu_export) == size(cpu_export)
    @test isapprox(gpu_export, cpu_export; rtol=1f-5, atol=1f-4)
    @test size(gpu_frame) == size(cpu_frame)
    @test isapprox(gpu_frame, cpu_frame; rtol=1f-5, atol=1f-4)
    cpu_sh_stats = ShackHartmann(cpu_tel; n_subap=4, mode=Diffractive(), T=T, backend=Array,
        valid_subaperture_policy=FluxThresholdValidSubapertures(light_ratio=0.5f0))
    gpu_sh_stats = ShackHartmann(gpu_tel; n_subap=4, mode=Diffractive(), T=T, backend=backend,
        valid_subaperture_policy=FluxThresholdValidSubapertures(light_ratio=0.5f0))
    measure!(cpu_sh_stats, cpu_tel, cpu_src, cpu_det; rng=MersenneTwister(3))
    measure!(gpu_sh_stats, gpu_tel, gpu_src, gpu_det; rng=MersenneTwister(3))
    cpu_peak = AdaptiveOpticsSim.sh_safe_peak_value(cpu_sh_stats.state.spot_cube)
    cpu_cutoff = AdaptiveOpticsSim.centroid_threshold(cpu_sh_stats) * cpu_peak
    AdaptiveOpticsSim.sh_signal_from_spots!(cpu_sh_stats, cpu_cutoff)
    gpu_peak = AdaptiveOpticsSim.sh_safe_peak_value(gpu_sh_stats.state.spot_cube)
    gpu_cutoff = AdaptiveOpticsSim.centroid_threshold(gpu_sh_stats) * gpu_peak
    AdaptiveOpticsSim.sh_signal_from_spots_device_stats!(
        AdaptiveOpticsSim.execution_style(gpu_sh_stats.state.slopes),
        gpu_sh_stats,
        gpu_cutoff,
    )
    @test isapprox(Array(gpu_sh_stats.state.slopes), cpu_sh_stats.state.slopes; rtol=1f-5, atol=1f-4)

    correction_models = (
        ReferencePixelCommonModeCorrection(1, 1),
        ReferenceRowCommonModeCorrection(1),
        ReferenceColumnCommonModeCorrection(1),
        ReferenceOutputCommonModeCorrection(2; edge_rows=1, edge_cols=1),
        CompositeFrameReadoutCorrection((
            ReferenceRowCommonModeCorrection(1),
            ReferenceColumnCommonModeCorrection(1),
        )),
    )
    corr_input = reshape(T.(1:96), 2, 6, 8)
    cpu_quant_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
        bits=8, full_well=T(100), sensor=CMOSSensor(T=T),
        response_model=NullFrameResponse(), T=T, backend=Array)
    gpu_quant_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
        bits=8, full_well=T(100), sensor=CMOSSensor(T=T),
        response_model=NullFrameResponse(), T=T, backend=backend)
    quant_input = reshape(T.(1:48), 6, 8) .* T(3)
    cpu_quant_frame = capture!(cpu_quant_det, copy(quant_input); rng=MersenneTwister(13))
    gpu_quant_frame = capture!(gpu_quant_det, backend(copy(quant_input)); rng=MersenneTwister(13))
    @test gpu_quant_frame isa backend
    @test isapprox(Array(gpu_quant_frame), cpu_quant_frame; rtol=1f-5, atol=1f-4)

    for correction_model in correction_models
        cpu_frame_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
            sensor=HgCdTeAvalancheArraySensor(T=T),
            response_model=NullFrameResponse(),
            correction_model=correction_model,
            T=T,
            backend=Array)
        gpu_frame_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
            sensor=HgCdTeAvalancheArraySensor(T=T),
            response_model=NullFrameResponse(),
            correction_model=correction_model,
            T=T,
            backend=backend)
        frame_input = reshape(T.(1:48), 6, 8)
        cpu_frame = capture!(cpu_frame_det, frame_input; rng=MersenneTwister(12))
        gpu_frame = capture!(gpu_frame_det, backend(copy(frame_input)); rng=MersenneTwister(12))
        @test gpu_frame isa backend
        @test isapprox(Array(gpu_frame), cpu_frame; rtol=1f-5, atol=1f-4)

        cpu_corr_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
            sensor=HgCdTeAvalancheArraySensor(T=T),
            response_model=NullFrameResponse(),
            correction_model=correction_model,
            T=T,
            backend=Array)
        gpu_corr_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
            sensor=HgCdTeAvalancheArraySensor(T=T),
            response_model=NullFrameResponse(),
            correction_model=correction_model,
            T=T,
            backend=backend)
        cpu_corr_cube = copy(corr_input)
        gpu_corr_cube = backend(copy(corr_input))
        cpu_corr_scratch = similar(cpu_corr_cube)
        gpu_corr_scratch = similar(gpu_corr_cube)
        AdaptiveOpticsSim.capture_stack!(cpu_corr_det, cpu_corr_cube, cpu_corr_scratch; rng=MersenneTwister(11))
        AdaptiveOpticsSim.capture_stack!(gpu_corr_det, gpu_corr_cube, gpu_corr_scratch; rng=MersenneTwister(11))
        @test gpu_corr_cube isa backend
        @test isapprox(Array(gpu_corr_cube), cpu_corr_cube; rtol=1f-5, atol=1f-4)
    end

    generalized_correction = CompositeFrameReadoutCorrection((
        ReferenceRowCommonModeCorrection(1),
        ReferenceColumnCommonModeCorrection(1),
    ))
    cpu_gen_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
        bits=8, full_well=T(100), readout_window=FrameWindow(2:5, 3:7), output_precision=UInt16,
        sensor=HgCdTeAvalancheArraySensor(T=T),
        response_model=NullFrameResponse(),
        correction_model=generalized_correction,
        T=T,
        backend=Array)
    gpu_gen_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
        bits=8, full_well=T(100), readout_window=FrameWindow(2:5, 3:7), output_precision=UInt16,
        sensor=HgCdTeAvalancheArraySensor(T=T),
        response_model=NullFrameResponse(),
        correction_model=generalized_correction,
        T=T,
        backend=backend)
    cpu_gen_input = copy(corr_input)
    gpu_gen_input = backend(copy(corr_input))
    cpu_gen_out = Array{UInt16}(undef, 2, 4, 5)
    gpu_gen_out = backend{UInt16}(undef, 2, 4, 5)
    AdaptiveOpticsSim.capture_stack!(cpu_gen_det, cpu_gen_out, cpu_gen_input; rng=MersenneTwister(14))
    AdaptiveOpticsSim.capture_stack!(gpu_gen_det, gpu_gen_out, gpu_gen_input; rng=MersenneTwister(14))
    @test gpu_gen_out isa backend
    @test Array(gpu_gen_out) == cpu_gen_out

    cpu_windowed_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0), binning=1,
        gain=T(1.0),
        sensor=HgCdTeAvalancheArraySensor(T=T, sampling_mode=CorrelatedDoubleSampling()),
        response_model=NullFrameResponse(),
        readout_window=FrameWindow(2:3, 2:3),
        T=T,
        backend=Array)
    gpu_windowed_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0), binning=1,
        gain=T(1.0),
        sensor=HgCdTeAvalancheArraySensor(T=T, sampling_mode=CorrelatedDoubleSampling()),
        response_model=NullFrameResponse(),
        readout_window=FrameWindow(2:3, 2:3),
        T=T,
        backend=backend)
    windowed_input = fill(T(10), 4, 4)
    cpu_windowed_frame = capture!(cpu_windowed_det, windowed_input; rng=MersenneTwister(15))
    gpu_windowed_frame = capture!(gpu_windowed_det, backend(copy(windowed_input)); rng=MersenneTwister(15))
    @test isapprox(Array(gpu_windowed_frame), cpu_windowed_frame; rtol=1f-5, atol=1f-4)
    @test isapprox(Array(detector_reference_frame(gpu_windowed_det)), detector_reference_frame(cpu_windowed_det); rtol=1f-5, atol=1f-4)
    @test isapprox(Array(detector_signal_frame(gpu_windowed_det)), detector_signal_frame(cpu_windowed_det); rtol=1f-5, atol=1f-4)
    @test isapprox(Array(detector_combined_frame(gpu_windowed_det)), detector_combined_frame(cpu_windowed_det); rtol=1f-5, atol=1f-4)
    @test isapprox(Array(detector_reference_cube(gpu_windowed_det)), detector_reference_cube(cpu_windowed_det); rtol=1f-5, atol=1f-4)
    @test isapprox(Array(detector_signal_cube(gpu_windowed_det)), detector_signal_cube(cpu_windowed_det); rtol=1f-5, atol=1f-4)
    @test isapprox(Array(detector_read_cube(gpu_windowed_det)), detector_read_cube(cpu_windowed_det); rtol=1f-5, atol=1f-4)
    @test detector_read_times(gpu_windowed_det) == detector_read_times(cpu_windowed_det)

    return nothing
end

function run_optional_backend_smoke(::Type{B}) where {B<:GPUBackendTag}
    pkg = backend_package_name(B)
    pkg_path = Base.find_package(pkg)
    if pkg_path === nothing
        @info "Skipping $(backend_label(B)) smoke: $(pkg).jl is not available in this environment"
        @test true
        return nothing
    end

    import_backend_package!(B)
    if !backend_functional(B)
        @info "Skipping $(backend_label(B)) smoke: backend runtime/device is not functional on this host"
        @test true
        return nothing
    end

    AdaptiveOpticsSim.disable_scalar_backend!(B)
    backend = AdaptiveOpticsSim.gpu_backend_array_type(B)
    @test backend !== nothing

    if get(ENV, backend_full_smoke_env(B), "0") == "1"
        include(joinpath(dirname(@__DIR__), "scripts", "gpu_smoke_contract.jl"))
        run_gpu_smoke_matrix(B)
        @test true
        return nothing
    end

    T = Float32
    rng = MersenneTwister(7)
    BackendArray = backend

    tel = Telescope(resolution=16, diameter=8.0f0, sampling_time=1.0f-3,
        central_obstruction=0.0f0, T=T, backend=BackendArray)
    src = Source(band=:I, magnitude=0.0, T=T)

    atm = MultiLayerAtmosphere(tel;
        r0=T(0.2),
        L0=T(25.0),
        fractional_cn2=T[0.7, 0.3],
        wind_speed=T[8.0, 4.0],
        wind_direction=T[0.0, 90.0],
        altitude=T[0.0, 5000.0],
        T=T,
        backend=BackendArray,
    )
    advance!(atm, tel; rng=rng)
    propagate!(atm, tel, src)
    @test atm.state.opd isa BackendArray
    @test tel.state.opd isa BackendArray

    inf_atm = InfiniteMultiLayerAtmosphere(tel;
        r0=T(0.2),
        L0=T(25.0),
        fractional_cn2=T[0.7, 0.3],
        wind_speed=T[8.0, 4.0],
        wind_direction=T[0.0, 90.0],
        altitude=T[0.0, 5000.0],
        screen_resolution=33,
        stencil_size=35,
        T=T,
        backend=BackendArray,
    )
    advance!(inf_atm, tel; rng=rng)
    propagate!(inf_atm, tel, src)
    @test inf_atm.state.opd isa BackendArray
    @test inf_atm.layers[1].screen.state.screen isa BackendArray

    prop = AtmosphericFieldPropagation(atm, tel, src;
        model=GeometricAtmosphericPropagation(T=T),
        zero_padding=2,
        T=T)
    field = propagate_atmosphere_field!(prop, atm, tel, src)
    @test field.state.field isa BackendArray
    intensity = atmospheric_intensity!(prop, atm, tel, src)
    @test intensity isa BackendArray

    bundle = SpectralBundle(T[0.9 * wavelength(src), 1.1 * wavelength(src)], T[0.4, 0.6]; T=T)
    poly = with_spectrum(src, bundle)
    sh = ShackHartmann(tel; n_subap=4, mode=Diffractive(), T=T, backend=BackendArray)
    slopes = measure!(sh, tel, poly)
    @test slopes isa BackendArray

    run_optional_backend_plan_checks(B, tel, BackendArray)

    curv = CurvatureWFS(tel; n_subap=4, T=T, backend=BackendArray)
    curv_slopes = measure!(curv, tel, src, atm)
    @test curv_slopes isa BackendArray

    single_platform = build_platform_scenario(
        SinglePlatformConfig(
            name=:optional_backend_single,
            branch_label=:main,
            products=RuntimeProductRequirements(slopes=true, wfs_pixels=true, science_pixels=false),
        ),
        build_optional_platform_branch(T, BackendArray, :main; sensor=:sh, seed=21),
    )
    prepare!(single_platform)
    step!(single_platform)
    @test simulation_interface(single_platform) isa SimulationInterface
    @test simulation_wfs_frame(single_platform) isa BackendArray
    @test platform_branch_labels(single_platform) == (:main,)

    grouped_platform = build_platform_scenario(
        GroupedPlatformConfig(
            (:hi, :lo);
            name=:optional_backend_grouped,
            products=GroupedRuntimeProductRequirements(wfs_frames=true, science_frames=false, wfs_stack=true, science_stack=false),
        ),
        build_optional_platform_branch(T, BackendArray, :hi; sensor=:sh, seed=31),
        build_optional_platform_branch(T, BackendArray, :lo; sensor=:sh, seed=32),
    )
    prepare!(grouped_platform)
    step!(grouped_platform)
    @test simulation_interface(grouped_platform) isa CompositeSimulationInterface
    @test simulation_grouped_wfs_stack(grouped_platform) isa BackendArray
    @test size(simulation_grouped_wfs_stack(grouped_platform), ndims(simulation_grouped_wfs_stack(grouped_platform))) == 2
    @test platform_branch_labels(grouped_platform) == (:hi, :lo)
    return nothing
end
