backend_package_name(::Type{AdaptiveOpticsSim.CUDABackendTag}) = "CUDA"
backend_package_name(::Type{AdaptiveOpticsSim.AMDGPUBackendTag}) = "AMDGPU"

backend_label(::Type{AdaptiveOpticsSim.CUDABackendTag}) = "CUDA"
backend_label(::Type{AdaptiveOpticsSim.AMDGPUBackendTag}) = "AMDGPU"

backend_full_smoke_env(::Type{AdaptiveOpticsSim.CUDABackendTag}) = "ADAPTIVEOPTICS_TEST_FULL_CUDA"
backend_full_smoke_env(::Type{AdaptiveOpticsSim.AMDGPUBackendTag}) = "ADAPTIVEOPTICS_TEST_FULL_AMDGPU"

backend_selector(::Type{AdaptiveOpticsSim.CUDABackendTag}) = AdaptiveOpticsSim.CUDABackend()
backend_selector(::Type{AdaptiveOpticsSim.AMDGPUBackendTag}) = AdaptiveOpticsSim.AMDGPUBackend()

struct OptionalStaticAtmosphere{A,B<:AbstractArrayBackend} <: AdaptiveOpticsSim.AbstractAtmosphere
    screen::A
end

AdaptiveOpticsSim.backend(::OptionalStaticAtmosphere{<:Any,B}) where {B} = B()
AdaptiveOpticsSim.advance!(atm::OptionalStaticAtmosphere, tel::Telescope, rng::AbstractRNG) = atm
AdaptiveOpticsSim.advance!(atm::OptionalStaticAtmosphere, tel::Telescope; rng::AbstractRNG=Random.default_rng()) = atm

function AdaptiveOpticsSim.propagate!(atm::OptionalStaticAtmosphere, tel::Telescope)
    copyto!(tel.state.opd, atm.screen)
    return tel
end

function OptionalStaticAtmosphere(tel::Telescope; T::Type{<:AbstractFloat}=Float32, backend::AbstractArrayBackend=backend(tel))
    selector = AdaptiveOpticsSim.require_same_backend(tel, AdaptiveOpticsSim._resolve_backend_selector(backend))
    array_backend = AdaptiveOpticsSim._resolve_array_backend(selector)
    host = zeros(T, tel.params.resolution, tel.params.resolution)
    host .*= Array(tel.state.pupil)
    screen = array_backend{T}(undef, size(host)...)
    copyto!(screen, host)
    return OptionalStaticAtmosphere{typeof(screen),typeof(selector)}(screen)
end

function run_optional_backend_selector_smoke(::Type{B}, BackendArray) where {B<:AdaptiveOpticsSim.GPUBackendTag}
    selector = backend_selector(B)
    T = Float32
    tel = Telescope(resolution=8, diameter=T(1), sampling_time=T(1e-3), central_obstruction=T(0), T=T, backend=selector)
    dm = DeformableMirror(tel; n_act=2, influence_width=T(0.3), T=T, backend=selector)
    dm_dense = DeformableMirror(tel; n_act=2, influence_model=DenseInfluenceMatrix(Array(dm.state.modes)), T=T, backend=selector)
    zernike_modal = ModalControllableOptic(tel, ZernikeOpticBasis([2, 3]); T=T, backend=selector)
    cartesian_modal = ModalControllableOptic(tel, CartesianTiltBasis(; scale=T(0.1)); T=T, backend=selector)
    wfs = ShackHartmann(tel; n_lenslets=2, mode=Diffractive(), T=T, backend=selector)
    det = Detector(noise=NoiseNone(), integration_time=T(1), qe=T(1), binning=1, T=T, backend=selector)
    @test tel.state.opd isa BackendArray
    @test dm.state.coefs isa BackendArray
    @test dm_dense.state.coefs isa BackendArray
    @test dm_dense.state.modes isa BackendArray
    @test Array(dm_dense.state.modes) ≈ Array(dm.state.modes) atol=0 rtol=0
    @test zernike_modal.state.coefs isa BackendArray
    @test zernike_modal.state.modes isa BackendArray
    @test cartesian_modal.state.coefs isa BackendArray
    @test cartesian_modal.state.modes isa BackendArray
    @test wfs.state.slopes isa BackendArray
    @test det.state.frame isa BackendArray
    return nothing
end

function build_optional_platform_branch(::Type{T}, backend, label::Symbol; sensor::Symbol=:sh, seed::Integer=1) where {T<:AbstractFloat}
    tel = Telescope(resolution=16, diameter=T(8.0), sampling_time=T(1e-3),
        central_obstruction=T(0.0), T=T, backend=backend)
    src = Source(band=:I, magnitude=0.0, T=T)
    atm = KolmogorovAtmosphere(tel; r0=T(0.2), L0=T(25.0), T=T, backend=backend)
    dm = DeformableMirror(tel; n_act=4, influence_width=T(0.3), T=T, backend=backend)
    wfs = sensor == :sh ?
        ShackHartmann(tel; n_lenslets=4, mode=Diffractive(), T=T, backend=backend) :
        PyramidWFS(tel; pupil_samples=4, modulation=T(1.0), mode=Diffractive(), T=T, backend=backend)
    det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0), binning=1, T=T, backend=backend)
    sim = AOSimulation(tel, src, atm, dm, wfs)
    imat = interaction_matrix(dm, wfs, tel, src; amplitude=T(0.05))
    recon = ModalReconstructor(imat; gain=T(0.5))
    return RuntimeBranch(label, sim, recon; wfs_detector=det, rng=MersenneTwister(seed))
end

@inline _optional_low_order_label(::Val{:tiptilt}) = :tiptilt
@inline _optional_low_order_label(::Val{:steering}) = :steering
@inline _optional_low_order_label(::Val{:focus}) = :focus

function _build_optional_low_order_wfs(tel::Telescope, backend, ::Type{T}, ::Val{:sh}) where {T<:AbstractFloat}
    return ShackHartmann(tel; n_lenslets=4, mode=Diffractive(), T=T, backend=backend)
end

function _build_optional_low_order_wfs(tel::Telescope, backend, ::Type{T}, ::Val{:pyr}) where {T<:AbstractFloat}
    return PyramidWFS(tel; pupil_samples=4, modulation=T(1.0), mode=Diffractive(), T=T, backend=backend)
end

function _build_optional_low_order_wfs(tel::Telescope, backend, ::Type{T}, ::Val{:bio}) where {T<:AbstractFloat}
    return BioEdgeWFS(tel; pupil_samples=4, modulation=T(1.0), mode=Diffractive(), T=T, backend=backend)
end

function _build_optional_composite_optic_case(backend, ::Type{T}, ::Val{:tiptilt}, wfs_case::Val{W}=Val(:sh)) where {T<:AbstractFloat,W}
    tel = Telescope(resolution=16, diameter=T(8.0), sampling_time=T(1e-3),
        central_obstruction=T(0.0), T=T, backend=backend)
    src = Source(band=:I, magnitude=0.0, T=T)
    atm = OptionalStaticAtmosphere(tel; T=T, backend=backend)
    tiptilt = TipTiltMirror(tel; scale=T(0.1), T=T, backend=backend, label=:tiptilt)
    dm = DeformableMirror(tel; n_act=4, influence_width=T(0.3), T=T, backend=backend)
    optic = CompositeControllableOptic(:tiptilt => tiptilt, :dm => dm)
    wfs = _build_optional_low_order_wfs(tel, backend, T, wfs_case)
    det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0), binning=1, T=T, backend=backend)
    sim = AOSimulation(tel, src, atm, optic, wfs)
    return build_runtime_scenario(
        SingleRuntimeConfig(name=Symbol(:optional_backend_composite_tiptilt_, W), branch_label=:main,
            products=RuntimeProductRequirements(slopes=true, wfs_pixels=true, science_pixels=false)),
        RuntimeBranch(:main, sim, NullReconstructor(); wfs_detector=det, rng=MersenneTwister(91)),
    )
end

function _build_optional_composite_optic_case(backend, ::Type{T}, ::Val{:steering}, ::Val{:sh}=Val(:sh)) where {T<:AbstractFloat}
    tel = Telescope(resolution=16, diameter=T(8.0), sampling_time=T(1e-3),
        central_obstruction=T(0.0), T=T, backend=backend)
    src = Source(band=:I, magnitude=0.0, T=T)
    atm = OptionalStaticAtmosphere(tel; T=T, backend=backend)
    steering = ModalControllableOptic(tel, CartesianTiltBasis(; scale=T(0.1));
        labels=:steering, T=T, backend=backend)
    dm = DeformableMirror(tel; n_act=4, influence_width=T(0.3), T=T, backend=backend)
    optic = CompositeControllableOptic(:steering => steering, :dm => dm)
    wfs = ShackHartmann(tel; n_lenslets=4, mode=Diffractive(), T=T, backend=backend)
    det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0), binning=1, T=T, backend=backend)
    sim = AOSimulation(tel, src, atm, optic, wfs)
    return build_runtime_scenario(
        SingleRuntimeConfig(name=:optional_backend_composite_steering, branch_label=:main,
            products=RuntimeProductRequirements(slopes=true, wfs_pixels=true, science_pixels=false)),
        RuntimeBranch(:main, sim, NullReconstructor(); wfs_detector=det, rng=MersenneTwister(91)),
    )
end

function _build_optional_composite_optic_case(backend, ::Type{T}, ::Val{:focus}, ::Val{:sh}=Val(:sh)) where {T<:AbstractFloat}
    tel = Telescope(resolution=16, diameter=T(8.0), sampling_time=T(1e-3),
        central_obstruction=T(0.0), T=T, backend=backend)
    src = Source(band=:I, magnitude=0.0, T=T)
    atm = OptionalStaticAtmosphere(tel; T=T, backend=backend)
    focus = FocusStage(tel; scale=T(0.1), T=T, backend=backend, label=:focus)
    dm = DeformableMirror(tel; n_act=4, influence_width=T(0.3), T=T, backend=backend)
    optic = CompositeControllableOptic(:focus => focus, :dm => dm)
    wfs = ShackHartmann(tel; n_lenslets=4, mode=Diffractive(), T=T, backend=backend)
    det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0), binning=1, T=T, backend=backend)
    sim = AOSimulation(tel, src, atm, optic, wfs)
    return build_runtime_scenario(
        SingleRuntimeConfig(name=:optional_backend_composite_focus, branch_label=:main,
            products=RuntimeProductRequirements(slopes=true, wfs_pixels=true, science_pixels=false)),
        RuntimeBranch(:main, sim, NullReconstructor(); wfs_detector=det, rng=MersenneTwister(91)),
    )
end

function _optional_low_order_commands(::Type{T}, ::Val{:focus}) where {T<:AbstractFloat}
    return fill(T(0.0125), 1), fill(T(0.025), 1), fill(T(0.02), 16)
end

function _optional_low_order_commands(::Type{T}, ::Union{Val{:tiptilt},Val{:steering}}) where {T<:AbstractFloat}
    return fill(T(0.0125), 2), fill(T(0.025), 2), fill(T(0.02), 16)
end

function _optional_low_order_tolerances(::Val{:focus}, ::Val{:sh}=Val(:sh))
    return (slopes_rtol=5f-3, slopes_atol=6f-3, frame_rtol=6f-3, frame_atol=1f6)
end

function _optional_low_order_tolerances(::Union{Val{:tiptilt},Val{:steering}}, ::Union{Val{:sh},Val{:pyr}})
    return (slopes_rtol=3f-3, slopes_atol=4f-3, frame_rtol=4f-3, frame_atol=1f6)
end

function _optional_low_order_tolerances(::Val{:tiptilt}, ::Val{:bio})
    return (slopes_rtol=1.5f-1, slopes_atol=4f-3, frame_rtol=6f-1, frame_atol=1f6)
end

function _run_optional_composite_optic_case(::Type{B}, case::Val{K}, wfs_case::Val{W}=Val(:sh)) where {B<:AdaptiveOpticsSim.GPUBackendTag,K,W}
    T = Float32
    cpu = _build_optional_composite_optic_case(CPUBackend(), T, case, wfs_case)
    gpu = _build_optional_composite_optic_case(backend_selector(B), T, case, wfs_case)
    prepare!(cpu)
    prepare!(gpu)
    initial_low_order, updated_low_order, dm_cmd = _optional_low_order_commands(T, case)
    label = _optional_low_order_label(case)
    tol = _optional_low_order_tolerances(case, wfs_case)
    set_command!(cpu, NamedTuple{(label, :dm)}((initial_low_order, dm_cmd)))
    set_command!(gpu, NamedTuple{(label, :dm)}((initial_low_order, dm_cmd)))
    sense!(cpu)
    sense!(gpu)
    AdaptiveOpticsSim.synchronize_backend!(AdaptiveOpticsSim.execution_style(command(gpu)))
    AdaptiveOpticsSim.synchronize_backend!(AdaptiveOpticsSim.execution_style(slopes(gpu)))
    AdaptiveOpticsSim.synchronize_backend!(AdaptiveOpticsSim.execution_style(wfs_frame(gpu)))
    @test command_segment_labels(command_layout(gpu)) == (label, :dm)
    @test isapprox(Array(command(gpu)), Array(command(cpu)); rtol=1f-6, atol=1f-6)
    @test isapprox(Array(slopes(gpu)), Array(slopes(cpu)); rtol=tol.slopes_rtol, atol=tol.slopes_atol)
    @test isapprox(Array(wfs_frame(gpu)), Array(wfs_frame(cpu)); rtol=tol.frame_rtol, atol=tol.frame_atol)

    update_command!(cpu, NamedTuple{(label,)}((updated_low_order,)))
    update_command!(gpu, NamedTuple{(label,)}((updated_low_order,)))
    sense!(cpu)
    sense!(gpu)
    AdaptiveOpticsSim.synchronize_backend!(AdaptiveOpticsSim.execution_style(command(gpu)))
    AdaptiveOpticsSim.synchronize_backend!(AdaptiveOpticsSim.execution_style(slopes(gpu)))
    AdaptiveOpticsSim.synchronize_backend!(AdaptiveOpticsSim.execution_style(wfs_frame(gpu)))
    @test isapprox(Array(command(gpu)), Array(command(cpu)); rtol=1f-6, atol=1f-6)
    @test isapprox(Array(slopes(gpu)), Array(slopes(cpu)); rtol=tol.slopes_rtol, atol=tol.slopes_atol)
    @test isapprox(Array(wfs_frame(gpu)), Array(wfs_frame(cpu)); rtol=tol.frame_rtol, atol=tol.frame_atol)
    return nothing
end

function run_optional_composite_optic_parity(::Type{B}, BackendArray) where {B<:AdaptiveOpticsSim.GPUBackendTag}
    _run_optional_composite_optic_case(B, Val(:tiptilt), Val(:sh))
    _run_optional_composite_optic_case(B, Val(:tiptilt), Val(:pyr))
    _run_optional_composite_optic_case(B, Val(:tiptilt), Val(:bio))
    _run_optional_composite_optic_case(B, Val(:steering), Val(:sh))
    _run_optional_composite_optic_case(B, Val(:focus), Val(:sh))
    return nothing
end

function import_backend_package!(::Type{AdaptiveOpticsSim.CUDABackendTag})
    @eval import CUDA
    return nothing
end

function import_backend_package!(::Type{AdaptiveOpticsSim.AMDGPUBackendTag})
    @eval import AMDGPU
    return nothing
end

function backend_functional(::Type{AdaptiveOpticsSim.CUDABackendTag})
    return Base.invokelatest(getproperty(getfield(Main, :CUDA), :functional))
end

function backend_functional(::Type{AdaptiveOpticsSim.AMDGPUBackendTag})
    return Base.invokelatest(getproperty(getfield(Main, :AMDGPU), :functional))
end

run_optional_backend_plan_checks(::Type{<:AdaptiveOpticsSim.GPUBackendTag}, tel, backend) = nothing

function run_optional_backend_plan_checks(::Type{AdaptiveOpticsSim.AMDGPUBackendTag}, tel, backend)
    T = Float32
    array_backend = AdaptiveOpticsSim._resolve_array_backend(backend)
    sh = ShackHartmann(tel; n_lenslets=4, mode=Diffractive(), T=T, backend=backend)
    pyr = PyramidWFS(tel; pupil_samples=4, modulation=T(1.0), mode=Diffractive(), T=T, backend=backend)
    bio = BioEdgeWFS(tel; pupil_samples=4, modulation=T(1.0), mode=Diffractive(), T=T, backend=backend)
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
    @test typeof(AdaptiveOpticsSim.sh_sensing_execution_plan(AdaptiveOpticsSim.execution_style(sh.state.slopes), sh)) === AdaptiveOpticsSim.ShackHartmannRocmSafePlan
    @test AdaptiveOpticsSim.detector_execution_plan(typeof(AdaptiveOpticsSim.execution_style(det.state.frame)), typeof(det)) isa AdaptiveOpticsSim.DetectorHostMirrorPlan
    capture_psf = array_backend{T}(undef, 4, 4)
    fill!(capture_psf, T(10))
    captured = capture!(det_capture, capture_psf; rng=MersenneTwister(2))
    @test captured isa array_backend
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
        central_obstruction=0.0f0, T=T, backend=CPUBackend())
    gpu_tel = Telescope(resolution=16, diameter=8.0f0, sampling_time=1.0f-3,
        central_obstruction=0.0f0, T=T, backend=backend)
    cpu_src = Source(band=:I, magnitude=0.0, T=T)
    gpu_src = Source(band=:I, magnitude=0.0, T=T)
    cpu_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
        sensor=CMOSSensor(T=T), response_model=NullFrameResponse(), T=T, backend=CPUBackend())
    gpu_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
        sensor=CMOSSensor(T=T), response_model=NullFrameResponse(), T=T, backend=backend)
    cpu_sh_stats = ShackHartmann(cpu_tel; n_lenslets=4, mode=Diffractive(), T=T, backend=CPUBackend(),
        valid_subaperture_policy=FluxThresholdValidSubapertures(light_ratio=0.5f0))
    gpu_sh_stats = ShackHartmann(gpu_tel; n_lenslets=4, mode=Diffractive(), T=T, backend=backend,
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
        response_model=NullFrameResponse(), T=T, backend=CPUBackend())
    gpu_quant_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
        bits=8, full_well=T(100), sensor=CMOSSensor(T=T),
        response_model=NullFrameResponse(), T=T, backend=backend)
    quant_input = reshape(T.(1:48), 6, 8) .* T(3)
    cpu_quant_frame = capture!(cpu_quant_det, copy(quant_input); rng=MersenneTwister(13))
    gpu_quant_frame = capture!(gpu_quant_det, array_backend(copy(quant_input)); rng=MersenneTwister(13))
    @test gpu_quant_frame isa array_backend
    @test isapprox(Array(gpu_quant_frame), cpu_quant_frame; rtol=1f-5, atol=1f-4)

    for correction_model in correction_models
        cpu_frame_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
            sensor=HgCdTeAvalancheArraySensor(T=T),
            response_model=NullFrameResponse(),
            correction_model=correction_model,
            T=T,
            backend=CPUBackend())
        gpu_frame_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
            sensor=HgCdTeAvalancheArraySensor(T=T),
            response_model=NullFrameResponse(),
            correction_model=correction_model,
            T=T,
            backend=backend)
        frame_input = reshape(T.(1:48), 6, 8)
        cpu_frame = capture!(cpu_frame_det, frame_input; rng=MersenneTwister(12))
        gpu_frame = capture!(gpu_frame_det, array_backend(copy(frame_input)); rng=MersenneTwister(12))
        @test gpu_frame isa array_backend
        @test isapprox(Array(gpu_frame), cpu_frame; rtol=1f-5, atol=1f-4)

        cpu_corr_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
            sensor=HgCdTeAvalancheArraySensor(T=T),
            response_model=NullFrameResponse(),
            correction_model=correction_model,
            T=T,
            backend=CPUBackend())
        gpu_corr_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
            sensor=HgCdTeAvalancheArraySensor(T=T),
            response_model=NullFrameResponse(),
            correction_model=correction_model,
            T=T,
            backend=backend)
        cpu_corr_cube = copy(corr_input)
        gpu_corr_cube = array_backend(copy(corr_input))
        cpu_corr_scratch = similar(cpu_corr_cube)
        gpu_corr_scratch = similar(gpu_corr_cube)
        AdaptiveOpticsSim.capture_stack!(cpu_corr_det, cpu_corr_cube, cpu_corr_scratch; rng=MersenneTwister(11))
        AdaptiveOpticsSim.capture_stack!(gpu_corr_det, gpu_corr_cube, gpu_corr_scratch; rng=MersenneTwister(11))
        @test gpu_corr_cube isa array_backend
        @test isapprox(Array(gpu_corr_cube), cpu_corr_cube; rtol=1f-5, atol=1f-4)
    end

    generalized_correction = CompositeFrameReadoutCorrection((
        ReferenceRowCommonModeCorrection(1),
        ReferenceColumnCommonModeCorrection(1),
    ))
    cpu_gen_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
        bits=8, full_well=T(100), readout_window=FrameWindow(2:5, 3:7), output_type=UInt16,
        sensor=HgCdTeAvalancheArraySensor(T=T),
        response_model=NullFrameResponse(),
        correction_model=generalized_correction,
        T=T,
        backend=CPUBackend())
    gpu_gen_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
        bits=8, full_well=T(100), readout_window=FrameWindow(2:5, 3:7), output_type=UInt16,
        sensor=HgCdTeAvalancheArraySensor(T=T),
        response_model=NullFrameResponse(),
        correction_model=generalized_correction,
        T=T,
        backend=backend)
    cpu_gen_input = copy(corr_input)
    gpu_gen_input = array_backend(copy(corr_input))
    cpu_gen_out = Array{UInt16}(undef, 2, 4, 5)
    gpu_gen_out = array_backend{UInt16}(undef, 2, 4, 5)
    AdaptiveOpticsSim.capture_stack!(cpu_gen_det, cpu_gen_out, cpu_gen_input; rng=MersenneTwister(14))
    AdaptiveOpticsSim.capture_stack!(gpu_gen_det, gpu_gen_out, gpu_gen_input; rng=MersenneTwister(14))
    @test gpu_gen_out isa array_backend
    @test Array(gpu_gen_out) == cpu_gen_out

    cpu_windowed_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0), binning=1,
        gain=T(1.0),
        sensor=HgCdTeAvalancheArraySensor(T=T, sampling_mode=CorrelatedDoubleSampling()),
        response_model=NullFrameResponse(),
        readout_window=FrameWindow(2:3, 2:3),
        T=T,
        backend=CPUBackend())
    gpu_windowed_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0), binning=1,
        gain=T(1.0),
        sensor=HgCdTeAvalancheArraySensor(T=T, sampling_mode=CorrelatedDoubleSampling()),
        response_model=NullFrameResponse(),
        readout_window=FrameWindow(2:3, 2:3),
        T=T,
        backend=backend)
    windowed_input = fill(T(10), 4, 4)
    cpu_windowed_frame = capture!(cpu_windowed_det, windowed_input; rng=MersenneTwister(15))
    gpu_windowed_frame = capture!(gpu_windowed_det, array_backend(copy(windowed_input)); rng=MersenneTwister(15))
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

function run_optional_backend_plan_checks(::Type{AdaptiveOpticsSim.CUDABackendTag}, tel, backend)
    T = Float32
    array_backend = AdaptiveOpticsSim._resolve_array_backend(backend)
    sh = ShackHartmann(tel; n_lenslets=4, mode=Diffractive(), T=T, backend=backend)
    pyr = PyramidWFS(tel; pupil_samples=4, modulation=T(1.0), mode=Diffractive(), T=T, backend=backend)
    bio = BioEdgeWFS(tel; pupil_samples=4, modulation=T(1.0), mode=Diffractive(), T=T, backend=backend)
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
        central_obstruction=0.0f0, T=T, backend=CPUBackend())
    gpu_tel = Telescope(resolution=16, diameter=8.0f0, sampling_time=1.0f-3,
        central_obstruction=0.0f0, T=T, backend=backend)
    cpu_src = Source(band=:I, magnitude=0.0, T=T)
    gpu_src = Source(band=:I, magnitude=0.0, T=T)
    cpu_sh = ShackHartmann(cpu_tel; n_lenslets=4, mode=Diffractive(), T=T, backend=CPUBackend())
    gpu_sh = ShackHartmann(gpu_tel; n_lenslets=4, mode=Diffractive(), T=T, backend=backend)
    cpu_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
        sensor=CMOSSensor(T=T), response_model=NullFrameResponse(), T=T, backend=CPUBackend())
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
    cpu_sh_stats = ShackHartmann(cpu_tel; n_lenslets=4, mode=Diffractive(), T=T, backend=CPUBackend(),
        valid_subaperture_policy=FluxThresholdValidSubapertures(light_ratio=0.5f0))
    gpu_sh_stats = ShackHartmann(gpu_tel; n_lenslets=4, mode=Diffractive(), T=T, backend=backend,
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
        response_model=NullFrameResponse(), T=T, backend=CPUBackend())
    gpu_quant_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
        bits=8, full_well=T(100), sensor=CMOSSensor(T=T),
        response_model=NullFrameResponse(), T=T, backend=backend)
    quant_input = reshape(T.(1:48), 6, 8) .* T(3)
    cpu_quant_frame = capture!(cpu_quant_det, copy(quant_input); rng=MersenneTwister(13))
    gpu_quant_frame = capture!(gpu_quant_det, array_backend(copy(quant_input)); rng=MersenneTwister(13))
    @test gpu_quant_frame isa array_backend
    @test isapprox(Array(gpu_quant_frame), cpu_quant_frame; rtol=1f-5, atol=1f-4)

    for correction_model in correction_models
        cpu_frame_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
            sensor=HgCdTeAvalancheArraySensor(T=T),
            response_model=NullFrameResponse(),
            correction_model=correction_model,
            T=T,
            backend=CPUBackend())
        gpu_frame_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
            sensor=HgCdTeAvalancheArraySensor(T=T),
            response_model=NullFrameResponse(),
            correction_model=correction_model,
            T=T,
            backend=backend)
        frame_input = reshape(T.(1:48), 6, 8)
        cpu_frame = capture!(cpu_frame_det, frame_input; rng=MersenneTwister(12))
        gpu_frame = capture!(gpu_frame_det, array_backend(copy(frame_input)); rng=MersenneTwister(12))
        @test gpu_frame isa array_backend
        @test isapprox(Array(gpu_frame), cpu_frame; rtol=1f-5, atol=1f-4)

        cpu_corr_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
            sensor=HgCdTeAvalancheArraySensor(T=T),
            response_model=NullFrameResponse(),
            correction_model=correction_model,
            T=T,
            backend=CPUBackend())
        gpu_corr_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
            sensor=HgCdTeAvalancheArraySensor(T=T),
            response_model=NullFrameResponse(),
            correction_model=correction_model,
            T=T,
            backend=backend)
        cpu_corr_cube = copy(corr_input)
        gpu_corr_cube = array_backend(copy(corr_input))
        cpu_corr_scratch = similar(cpu_corr_cube)
        gpu_corr_scratch = similar(gpu_corr_cube)
        AdaptiveOpticsSim.capture_stack!(cpu_corr_det, cpu_corr_cube, cpu_corr_scratch; rng=MersenneTwister(11))
        AdaptiveOpticsSim.capture_stack!(gpu_corr_det, gpu_corr_cube, gpu_corr_scratch; rng=MersenneTwister(11))
        @test gpu_corr_cube isa array_backend
        @test isapprox(Array(gpu_corr_cube), cpu_corr_cube; rtol=1f-5, atol=1f-4)
    end

    generalized_correction = CompositeFrameReadoutCorrection((
        ReferenceRowCommonModeCorrection(1),
        ReferenceColumnCommonModeCorrection(1),
    ))
    cpu_gen_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
        bits=8, full_well=T(100), readout_window=FrameWindow(2:5, 3:7), output_type=UInt16,
        sensor=HgCdTeAvalancheArraySensor(T=T),
        response_model=NullFrameResponse(),
        correction_model=generalized_correction,
        T=T,
        backend=CPUBackend())
    gpu_gen_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
        bits=8, full_well=T(100), readout_window=FrameWindow(2:5, 3:7), output_type=UInt16,
        sensor=HgCdTeAvalancheArraySensor(T=T),
        response_model=NullFrameResponse(),
        correction_model=generalized_correction,
        T=T,
        backend=backend)
    cpu_gen_input = copy(corr_input)
    gpu_gen_input = array_backend(copy(corr_input))
    cpu_gen_out = Array{UInt16}(undef, 2, 4, 5)
    gpu_gen_out = array_backend{UInt16}(undef, 2, 4, 5)
    AdaptiveOpticsSim.capture_stack!(cpu_gen_det, cpu_gen_out, cpu_gen_input; rng=MersenneTwister(14))
    AdaptiveOpticsSim.capture_stack!(gpu_gen_det, gpu_gen_out, gpu_gen_input; rng=MersenneTwister(14))
    @test gpu_gen_out isa array_backend
    @test Array(gpu_gen_out) == cpu_gen_out

    cpu_windowed_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0), binning=1,
        gain=T(1.0),
        sensor=HgCdTeAvalancheArraySensor(T=T, sampling_mode=CorrelatedDoubleSampling()),
        response_model=NullFrameResponse(),
        readout_window=FrameWindow(2:3, 2:3),
        T=T,
        backend=CPUBackend())
    gpu_windowed_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0), binning=1,
        gain=T(1.0),
        sensor=HgCdTeAvalancheArraySensor(T=T, sampling_mode=CorrelatedDoubleSampling()),
        response_model=NullFrameResponse(),
        readout_window=FrameWindow(2:3, 2:3),
        T=T,
        backend=backend)
    windowed_input = fill(T(10), 4, 4)
    cpu_windowed_frame = capture!(cpu_windowed_det, windowed_input; rng=MersenneTwister(15))
    gpu_windowed_frame = capture!(gpu_windowed_det, array_backend(copy(windowed_input)); rng=MersenneTwister(15))
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

function run_optional_backend_smoke(::Type{B}) where {B<:AdaptiveOpticsSim.GPUBackendTag}
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
    run_optional_backend_selector_smoke(B, backend)

    if get(ENV, backend_full_smoke_env(B), "0") == "1"
        include(joinpath(dirname(@__DIR__), "scripts", "gpu_smoke_contract.jl"))
        run_gpu_smoke_matrix(B)
        @test true
        return nothing
    end

    T = Float32
    rng = MersenneTwister(7)
    BackendArray = backend
    selector = backend_selector(B)

    tel = Telescope(resolution=16, diameter=8.0f0, sampling_time=1.0f-3,
        central_obstruction=0.0f0, T=T, backend=selector)
    src = Source(band=:I, magnitude=0.0, T=T)

    atm = MultiLayerAtmosphere(tel;
        r0=T(0.2),
        L0=T(25.0),
        fractional_cn2=T[0.7, 0.3],
        wind_speed=T[8.0, 4.0],
        wind_direction=T[0.0, 90.0],
        altitude=T[0.0, 5000.0],
        T=T,
        backend=selector,
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
        backend=selector,
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
    sh = ShackHartmann(tel; n_lenslets=4, mode=Diffractive(), T=T, backend=selector)
    slopes = measure!(sh, tel, poly)
    @test slopes isa BackendArray

    run_optional_backend_plan_checks(B, tel, selector)
    run_optional_composite_optic_parity(B, BackendArray)

    curv = CurvatureWFS(tel; pupil_samples=4, T=T, backend=selector)
    curv_slopes = measure!(curv, tel, src, atm)
    @test curv_slopes isa BackendArray

    single_platform = build_runtime_scenario(
        SingleRuntimeConfig(
            name=:optional_backend_single,
            branch_label=:main,
            products=RuntimeProductRequirements(slopes=true, wfs_pixels=true, science_pixels=false),
        ),
        build_optional_platform_branch(T, selector, :main; sensor=:sh, seed=21),
    )
    prepare!(single_platform)
    step!(single_platform)
    @test AdaptiveOpticsSim.simulation_interface(single_platform) isa SimulationInterface
    @test wfs_frame(single_platform) isa BackendArray
    @test platform_branch_labels(single_platform) == (:main,)

    grouped_platform = build_runtime_scenario(
        GroupedRuntimeConfig(
            (:hi, :lo);
            name=:optional_backend_grouped,
            products=GroupedRuntimeProductRequirements(wfs_frames=true, science_frames=false, wfs_stack=true, science_stack=false),
        ),
        build_optional_platform_branch(T, selector, :hi; sensor=:sh, seed=31),
        build_optional_platform_branch(T, selector, :lo; sensor=:sh, seed=32),
    )
    prepare!(grouped_platform)
    step!(grouped_platform)
    @test AdaptiveOpticsSim.simulation_interface(grouped_platform) isa CompositeSimulationInterface
    @test grouped_wfs_stack(grouped_platform) isa BackendArray
    @test size(grouped_wfs_stack(grouped_platform), ndims(grouped_wfs_stack(grouped_platform))) == 2
    @test platform_branch_labels(grouped_platform) == (:hi, :lo)
    return nothing
end
