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
    host .*= Array(pupil_mask(tel))
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
    wfs = ShackHartmannWFS(tel; n_lenslets=2, mode=Diffractive(), T=T, backend=selector)
    det = Detector(noise=NoiseNone(), integration_time=T(1), qe=T(1), binning=1, T=T, backend=selector)
    @test tel.state.opd isa BackendArray
    @test dm.state.coefs isa BackendArray
    @test dm.state.modes isa AdaptiveOpticsSim.GaussianInfluenceOperator
    @test typeof(backend(dm.state.modes)) === typeof(selector)
    @test similar(dm.state.modes, T, 2, 2) isa BackendArray
    @test AdaptiveOpticsSim.materialize_influence_matrix(dm) isa BackendArray
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

function run_optional_plane_product_checks(tel::Telescope,
    src::AdaptiveOpticsSim.AbstractSource,
    selector::AdaptiveOpticsSim.AbstractArrayBackend, BackendArray,
    ::Type{T}) where {T<:AbstractFloat}
    wavefront = PupilFunction(tel; T=T, backend=selector)
    apply_opd!(wavefront, opd_map(tel))
    field = ElectricField(wavefront, src; zero_padding=2, T=T)
    formation = prepare_pupil_field(tel, wavefront, src, field)
    fill_electric_field!(field, wavefront, formation)
    @test wavefront.amplitude isa BackendArray
    @test wavefront.opd isa BackendArray
    @test field.values isa BackendArray
    @test field.metadata.device == plane_device(field.values)

    prepared = prepare_direct_psf(tel, wavefront, src, field)
    compute_psf!(prepared.output, field, wavefront, prepared.plan,
        prepared.workspace)
    legacy_psf = copy(compute_psf!(tel, src; zero_padding=2))
    @test prepared.output.values isa BackendArray
    @test Array(prepared.output.values) ≈ Array(legacy_psf) atol=T(2e-5) rtol=T(2e-5)

    fraunhofer = FraunhoferPropagation(field)
    propagated = propagation_output(field, fraunhofer)
    propagate_field!(propagated, field, fraunhofer)
    @test propagated.values isa BackendArray
    @test propagated.metadata.kind isa FocalPlane

    spatial_filter = SpatialFilter(tel; shape=SquareFilter(), diameter=5,
        zero_padding=2, T=T, backend=selector)
    spatial_formation = prepare_pupil_field(tel, wavefront, src, field;
        center_even_grid=false, amplitude_scale=1)
    fill_electric_field!(field, wavefront, spatial_formation)
    filtered = PupilFunction(tel; T=T, backend=selector)
    spatial_plan = prepare_spatial_filter(tel, spatial_filter, field,
        filtered)
    spatial_workspace = SpatialFilterWorkspace(spatial_filter)
    filter!(filtered, field, spatial_filter, spatial_plan,
        spatial_workspace)
    @test filtered.opd isa BackendArray
    @test filtered.amplitude isa BackendArray
    @test all(isfinite, Array(filtered.opd))
    @test all(isfinite, Array(filtered.amplitude))

    telescope_opd = copy(Array(opd_map(tel)))
    dm = DeformableMirror(tel; n_act=4, influence_width=T(0.3), T=T,
        backend=selector)
    fill!(dm.state.coefs, T(1e-8))
    update_surface!(dm, tel)
    apply_surface!(wavefront, dm, DMReplace())
    @test Array(wavefront.opd) ≈ Array(surface_opd(dm)) atol=zero(T) rtol=zero(T)
    @test Array(opd_map(tel)) == telescope_opd
    return nothing
end

function build_optional_control_loop_branch(::Type{T}, backend, label::Symbol; sensor::Symbol=:sh, seed::Integer=1) where {T<:AbstractFloat}
    tel = Telescope(resolution=16, diameter=T(8.0), sampling_time=T(1e-3),
        central_obstruction=T(0.0), T=T, backend=backend)
    src = Source(band=:I, magnitude=0.0, T=T)
    atm = KolmogorovAtmosphere(tel; r0=T(0.2), L0=T(25.0), T=T, backend=backend)
    dm = DeformableMirror(tel; n_act=4, influence_width=T(0.3), T=T, backend=backend)
    wfs = sensor == :sh ?
        ShackHartmannWFS(tel; n_lenslets=4, mode=Diffractive(), T=T, backend=backend) :
        PyramidWFS(tel; pupil_samples=4, modulation=T(1.0), mode=Diffractive(), T=T, backend=backend)
    det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0), binning=1, T=T, backend=backend)
    sim = AOSimulation(tel, src, atm, dm, wfs)
    imat = interaction_matrix(dm, wfs, tel, src; amplitude=T(0.05))
    recon = ModalReconstructor(imat; gain=T(0.5))
    return ControlLoopBranch(label, sim, recon; wfs_detector=det, rng=MersenneTwister(seed))
end

function run_optional_scalable_reconstructor_checks(::Type{T}, selector,
    BackendArray) where {T<:AbstractFloat}
    rng = MersenneTwister(20260713)
    n_slopes = 12
    n_commands = 8
    D_host = randn(rng, T, n_slopes, n_commands)
    slopes_host = randn(rng, T, n_slopes)
    D = BackendArray(D_host)
    input = BackendArray(slopes_host)
    interaction = InteractionMatrix(D, T(0.1))
    dense = ModalReconstructor(interaction; gain=T(0.5))
    factorized = FactorizedReconstructor(interaction; gain=T(0.5))
    dense_out = BackendArray(zeros(T, n_commands))
    factorized_out = BackendArray(zeros(T, n_commands))
    reconstruct!(dense_out, dense, input)
    reconstruct!(factorized_out, factorized, input)
    AdaptiveOpticsSim.synchronize_backend!(
        AdaptiveOpticsSim.execution_style(factorized_out))
    @test Array(factorized_out) ≈ Array(dense_out) atol=T(2e-4) rtol=T(2e-4)
    @test AdaptiveOpticsSim.factorized_rank(factorized) == n_commands
    @test all(storage -> storage isa BackendArray,
        AdaptiveOpticsSim.runtime_reconstructor_storage(factorized))

    compact = FactorizedReconstructor(interaction; gain=T(0.5), max_rank=3)
    @test AdaptiveOpticsSim.factorized_rank(compact) == 3
    @test sum(length,
        AdaptiveOpticsSim.runtime_reconstructor_storage(compact)) <
        sum(length,
            AdaptiveOpticsSim.runtime_reconstructor_storage(factorized))

    controller = DiscreteIntegratorController(n_commands;
        gain=T(0.7), tau=T(0.01), T=T, backend=selector)
    controlled = ControlledReconstructor(factorized, controller; dt=T(1e-3))
    reconstruct!(factorized_out, controlled, input)
    AdaptiveOpticsSim.synchronize_backend!(
        AdaptiveOpticsSim.execution_style(factorized_out))
    @test all(isfinite, Array(factorized_out))
    @test all(storage -> storage isa BackendArray,
        AdaptiveOpticsSim.runtime_reconstructor_storage(controlled))
    @test AdaptiveOpticsSim.reset_controller!(controlled) === controlled
    @test_throws InvalidConfiguration ControlledReconstructor(
        factorized,
        DiscreteIntegratorController(n_commands; T=T, backend=CPUBackend());
        dt=T(1e-3),
    )
    return nothing
end

@inline _optional_low_order_label(::Val{:tiptilt}) = :tiptilt
@inline _optional_low_order_label(::Val{:steering}) = :steering
@inline _optional_low_order_label(::Val{:focus}) = :focus

function _build_optional_low_order_wfs(tel::Telescope, backend, ::Type{T}, ::Val{:sh}) where {T<:AbstractFloat}
    return ShackHartmannWFS(tel; n_lenslets=4, mode=Diffractive(), T=T, backend=backend)
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
    return build_control_loop_scenario(
        SingleControlLoopConfig(name=Symbol(:optional_backend_composite_tiptilt_, W), branch_label=:main,
            outputs=RuntimeOutputRequirements(slopes=true, wfs_pixels=true, science_pixels=false)),
        ControlLoopBranch(:main, sim, NullReconstructor(); wfs_detector=det, rng=MersenneTwister(91)),
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
    wfs = ShackHartmannWFS(tel; n_lenslets=4, mode=Diffractive(), T=T, backend=backend)
    det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0), binning=1, T=T, backend=backend)
    sim = AOSimulation(tel, src, atm, optic, wfs)
    return build_control_loop_scenario(
        SingleControlLoopConfig(name=:optional_backend_composite_steering, branch_label=:main,
            outputs=RuntimeOutputRequirements(slopes=true, wfs_pixels=true, science_pixels=false)),
        ControlLoopBranch(:main, sim, NullReconstructor(); wfs_detector=det, rng=MersenneTwister(91)),
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
    wfs = ShackHartmannWFS(tel; n_lenslets=4, mode=Diffractive(), T=T, backend=backend)
    det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0), binning=1, T=T, backend=backend)
    sim = AOSimulation(tel, src, atm, optic, wfs)
    return build_control_loop_scenario(
        SingleControlLoopConfig(name=:optional_backend_composite_focus, branch_label=:main,
            outputs=RuntimeOutputRequirements(slopes=true, wfs_pixels=true, science_pixels=false)),
        ControlLoopBranch(:main, sim, NullReconstructor(); wfs_detector=det, rng=MersenneTwister(91)),
    )
end

function _optional_low_order_commands(::Type{T}, ::Val{:focus}) where {T<:AbstractFloat}
    return fill(T(1.25e-7), 1), fill(T(2.5e-7), 1), fill(T(2e-8), 16)
end

function _optional_low_order_commands(::Type{T}, ::Union{Val{:tiptilt},Val{:steering}}) where {T<:AbstractFloat}
    # These coefficients produce nanometre-scale OPD. Centimetre-scale test
    # commands make diffraction output chaotic under otherwise sub-ulp
    # Float32 backend differences and do not represent an AO operating point.
    return fill(T(1.25e-7), 2), fill(T(2.5e-7), 2), fill(T(2e-8), 16)
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
    @test AdaptiveOpticsSim.command_segment_labels(AdaptiveOpticsSim.command_layout(gpu)) == (label, :dm)
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

function run_optional_counting_detector_parity(::Type{B}, BackendArray) where {B<:AdaptiveOpticsSim.GPUBackendTag}
    T = Float32
    selector = backend_selector(B)
    correlation = ChannelCrosstalkModel(T(0.4))
    sensor_kwargs = (
        qe=T(1),
        dark_count_rate=T(0),
        fill_factor=T(1),
        energy_resolution=T(12),
        timing_jitter_s=T(2e-6),
        wavelength_range_m=(T(0.8e-6), T(1.4e-6)),
        correlation_model=correlation,
        T=T,
    )
    cpu = MKIDArrayDetector(
        noise=NoiseNone(),
        sensor=MKIDArraySensor(; sensor_kwargs...),
        output_type=UInt16,
        T=T,
        backend=CPUBackend(),
    )
    gpu = MKIDArrayDetector(
        noise=NoiseNone(),
        sensor=MKIDArraySensor(; sensor_kwargs...),
        output_type=UInt16,
        T=T,
        backend=selector,
    )
    input = zeros(T, 5, 5)
    input[3, 3] = T(10)
    cpu_output = capture!(cpu, input, MersenneTwister(19))
    gpu_output = capture!(gpu, BackendArray(input), MersenneTwister(19))
    AdaptiveOpticsSim.synchronize_backend!(AdaptiveOpticsSim.execution_style(gpu_output))
    @test gpu_output isa BackendArray
    @test isapprox(Array(gpu_output), cpu_output; rtol=1f-6, atol=1f-6)
    @test isapprox(sum(Array(gpu_output)), sum(cpu_output); rtol=1f-6, atol=1f-6)

    inside_band = Source(band=:custom, wavelength=T(1.0e-6), T=T)
    outside_band = Source(band=:custom, wavelength=T(0.55e-6), T=T)
    gpu_inside = capture!(gpu, BackendArray(fill(T(2), 2, 2)), inside_band,
        MersenneTwister(20))
    AdaptiveOpticsSim.synchronize_backend!(AdaptiveOpticsSim.execution_style(gpu_inside))
    inside_host = Array(gpu_inside)
    gpu_outside = capture!(gpu, BackendArray(fill(T(2), 2, 2)), outside_band,
        MersenneTwister(20))
    AdaptiveOpticsSim.synchronize_backend!(AdaptiveOpticsSim.execution_style(gpu_outside))
    outside_host = Array(gpu_outside)
    @test inside_host == fill(T(2), 2, 2)
    @test outside_host == zeros(T, 2, 2)
    return nothing
end

function run_optional_avalanche_detector_parity(::Type{B}, BackendArray) where {B<:AdaptiveOpticsSim.GPUBackendTag}
    T = Float32
    selector = backend_selector(B)
    input = BackendArray(fill(one(T), 128, 128))

    pc_detector = Detector(noise=NoiseNone(), qe=one(T), gain=T(10),
        sensor=EMCCDSensor(operating_mode=PhotonCountingEMMode(
            threshold=T(5), detection_efficiency=T(0.5)), T=T),
        response_model=NullFrameResponse(), T=T, backend=selector)
    pc_output = capture!(pc_detector, input; rng=MersenneTwister(2026))
    AdaptiveOpticsSim.synchronize_backend!(
        AdaptiveOpticsSim.execution_style(pc_output))
    pc_host = Array(pc_output)
    @test all(x -> x == zero(T) || x == one(T), pc_host)
    @test T(0.47) <= sum(pc_host) / length(pc_host) <= T(0.53)

    stochastic_detector = Detector(noise=NoiseNone(), qe=one(T), gain=T(5),
        sensor=EMCCDSensor(excess_noise_factor=T(1.4),
            multiplication_model=AdaptiveOpticsSim.StochasticMultiplicationRegister(T(0.6)),
            T=T), response_model=NullFrameResponse(), T=T, backend=selector)
    stochastic_input = BackendArray(fill(T(50), 128, 128))
    stochastic_output = capture!(stochastic_detector, stochastic_input;
        rng=MersenneTwister(2027))
    AdaptiveOpticsSim.synchronize_backend!(
        AdaptiveOpticsSim.execution_style(stochastic_output))
    stochastic_host = Array(stochastic_output)
    @test all(x -> x >= zero(T), stochastic_host)
    @test isapprox(sum(stochastic_host) / length(stochastic_host), T(250);
        rtol=T(0.03))

    ramp_detector = Detector(integration_time=T(2), noise=NoiseNone(),
        qe=one(T), gain=one(T),
        sensor=HgCdTeAvalancheArraySensor(read_time=zero(T),
            sampling_mode=UpTheRampSampling(5), T=T),
        response_model=NullFrameResponse(), T=T, backend=selector)
    ramp_output = capture!(ramp_detector,
        BackendArray(fill(T(3), 32, 32)); rng=MersenneTwister(2028))
    AdaptiveOpticsSim.synchronize_backend!(
        AdaptiveOpticsSim.execution_style(ramp_output))
    @test Array(ramp_output) == fill(T(6), 32, 32)
    @test Array(detector_ramp_slope(ramp_detector)) == fill(T(3), 32, 32)
    @test maximum(abs, Array(detector_ramp_intercept(ramp_detector))) <= eps(T)
    @test size(detector_ramp_cube(ramp_detector)) == (32, 32, 5)
    @test detector_ramp_times(ramp_detector) == T[0, 0.5, 1, 1.5, 2]

    linear_apd = LinearAPDDetector(topology=APDChannelBank(4),
        integration_time=T(0.5), qe=T(0.5), avalanche_gain=T(4),
        conversion_gain=T(2), noise=NoiseNone(), T=T, backend=selector)
    linear_output = capture!(linear_apd, BackendArray(fill(T(10), 4));
        rng=MersenneTwister(2029))
    AdaptiveOpticsSim.synchronize_backend!(
        AdaptiveOpticsSim.execution_style(linear_output))
    @test linear_output isa BackendArray
    @test Array(linear_output) == fill(T(20), 4)
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
    sh = ShackHartmannWFS(tel; n_lenslets=4, mode=Diffractive(), T=T, backend=backend)
    sh_large = ShackHartmannWFS(tel; n_lenslets=8, mode=Diffractive(), T=T, backend=backend)
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
    @test typeof(AdaptiveOpticsSim.sh_sensing_execution_plan(
        AdaptiveOpticsSim.execution_style(sh.state.slopes), sh)) ===
        AdaptiveOpticsSim.ShackHartmannWFSRocmHostStatsPlan
    @test typeof(AdaptiveOpticsSim.sh_sensing_execution_plan(
        AdaptiveOpticsSim.execution_style(sh_large.state.slopes), sh_large)) ===
        AdaptiveOpticsSim.ShackHartmannWFSRocmHostStatsPlan
    @test AdaptiveOpticsSim.detector_execution_plan(typeof(AdaptiveOpticsSim.execution_style(det.state.frame)), typeof(det)) isa AdaptiveOpticsSim.DetectorHostMirrorPlan
    capture_psf = array_backend{T}(undef, 4, 4)
    fill!(capture_psf, T(10))
    captured = capture!(det_capture, capture_psf; rng=MersenneTwister(2))
    @test captured isa array_backend
    @test maximum(Array(captured)) <= Float64(exp2(T(12)) - one(T))
    cpu_poisson_det = Detector(noise=NoisePhoton(), integration_time=T(1.0), qe=T(1.0),
        sensor=CMOSSensor(T=T), response_model=NullFrameResponse(), T=T,
        backend=CPUBackend())
    gpu_poisson_det = Detector(noise=NoisePhoton(), integration_time=T(1.0), qe=T(1.0),
        sensor=CMOSSensor(T=T), response_model=NullFrameResponse(), T=T,
        backend=backend)
    poisson_input = fill(T(10), 4, 4)
    cpu_poisson = capture!(cpu_poisson_det, copy(poisson_input); rng=MersenneTwister(23))
    gpu_poisson = capture!(gpu_poisson_det, array_backend(copy(poisson_input));
        rng=MersenneTwister(23))
    @test Array(gpu_poisson) == cpu_poisson
    poisson_method = which(
        AdaptiveOpticsSim._poisson_noise_frame!,
        (
            AdaptiveOpticsSim.DetectorHostMirrorPlan,
            typeof(gpu_poisson_det),
            typeof(MersenneTwister(1)),
            typeof(gpu_poisson_det.state.frame),
        ),
    )
    @test occursin("AdaptiveOpticsSimAMDGPUExt", String(poisson_method.file))
    @test AdaptiveOpticsSim.reduction_execution_plan(pyr.state.intensity) isa AdaptiveOpticsSim.HostMirrorReductionPlan
    @test AdaptiveOpticsSim.backend_sum_value(capture_psf) == T(160)

    phase_freqs = T[-0.2, -0.1, 0.1, 0.2]
    cpu_phase_psd = zeros(T, 4, 4)
    gpu_phase_freqs = array_backend(phase_freqs)
    gpu_phase_psd = array_backend(zeros(T, 4, 4))
    phase_args = (T(0.02), T(4pi^2), T(0.01), T(-11 / 6), zero(T), 4)
    AdaptiveOpticsSim._fill_phase_psd!(AdaptiveOpticsSim.ScalarCPUStyle(),
        cpu_phase_psd, phase_freqs, phase_args...)
    phase_style = AdaptiveOpticsSim.execution_style(gpu_phase_psd)
    AdaptiveOpticsSim._fill_phase_psd!(phase_style, gpu_phase_psd,
        gpu_phase_freqs, phase_args...)
    @test isapprox(Array(gpu_phase_psd), cpu_phase_psd; rtol=1f-6, atol=1f-7)
    phase_method = which(
        AdaptiveOpticsSim._fill_phase_psd!,
        (typeof(phase_style), typeof(gpu_phase_psd), typeof(gpu_phase_freqs),
            T, T, T, T, T, Int),
    )
    @test occursin("AdaptiveOpticsSimAMDGPUExt", String(phase_method.file))

    host_opd = reshape(T.(1:256), 16, 16) .* T(1e-9)
    host_valid = trues(4, 4)
    cpu_geometric_slopes = zeros(T, 32)
    gpu_geometric_slopes = array_backend(zeros(T, 32))
    gpu_opd = array_backend(host_opd)
    gpu_valid = array_backend(host_valid)
    AdaptiveOpticsSim.geometric_slopes!(cpu_geometric_slopes, host_opd, host_valid)
    AdaptiveOpticsSim.geometric_slopes!(gpu_geometric_slopes, gpu_opd, gpu_valid)
    @test isapprox(Array(gpu_geometric_slopes), cpu_geometric_slopes;
        rtol=1f-6, atol=1f-8)
    geometric_method = which(
        AdaptiveOpticsSim._geometric_slopes!,
        (typeof(AdaptiveOpticsSim.execution_style(gpu_geometric_slopes)),
            typeof(gpu_geometric_slopes), typeof(gpu_opd), typeof(gpu_valid),
            Int, Int, Int),
    )
    @test occursin("AdaptiveOpticsSimAMDGPUExt", String(geometric_method.file))
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
        AdaptiveOpticsSim.execution_style(first(geom_prop.state.slices).field.values),
        geom_prop.params.model,
    ) isa AdaptiveOpticsSim.GeometricFieldAsyncPlan
    @test AdaptiveOpticsSim.atmospheric_field_execution_plan(
        AdaptiveOpticsSim.execution_style(first(fresnel_prop.state.slices).field.values),
        fresnel_prop.params.model,
    ) isa AdaptiveOpticsSim.LayeredFresnelFieldAsyncPlan
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
        bits=8, full_well=T(100), readout_window=AdaptiveOpticsSim.FrameWindow(2:5, 3:7), output_type=UInt16,
        sensor=HgCdTeAvalancheArraySensor(T=T),
        response_model=NullFrameResponse(),
        correction_model=generalized_correction,
        T=T,
        backend=CPUBackend())
    gpu_gen_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
        bits=8, full_well=T(100), readout_window=AdaptiveOpticsSim.FrameWindow(2:5, 3:7), output_type=UInt16,
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
        readout_window=AdaptiveOpticsSim.FrameWindow(2:3, 2:3),
        T=T,
        backend=CPUBackend())
    gpu_windowed_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0), binning=1,
        gain=T(1.0),
        sensor=HgCdTeAvalancheArraySensor(T=T, sampling_mode=CorrelatedDoubleSampling()),
        response_model=NullFrameResponse(),
        readout_window=AdaptiveOpticsSim.FrameWindow(2:3, 2:3),
        T=T,
        backend=backend)
    windowed_input = fill(T(10), 4, 4)
    cpu_windowed_frame = capture!(cpu_windowed_det, windowed_input; rng=MersenneTwister(15))
    gpu_windowed_frame = capture!(gpu_windowed_det, array_backend(copy(windowed_input)); rng=MersenneTwister(15))
    @test isapprox(Array(gpu_windowed_frame), cpu_windowed_frame; rtol=1f-5, atol=1f-4)
    @test isapprox(Array(AdaptiveOpticsSim.detector_reference_frame(gpu_windowed_det)), AdaptiveOpticsSim.detector_reference_frame(cpu_windowed_det); rtol=1f-5, atol=1f-4)
    @test isapprox(Array(AdaptiveOpticsSim.detector_signal_frame(gpu_windowed_det)), AdaptiveOpticsSim.detector_signal_frame(cpu_windowed_det); rtol=1f-5, atol=1f-4)
    @test isapprox(Array(AdaptiveOpticsSim.detector_combined_frame(gpu_windowed_det)), AdaptiveOpticsSim.detector_combined_frame(cpu_windowed_det); rtol=1f-5, atol=1f-4)
    @test isapprox(Array(AdaptiveOpticsSim.detector_reference_cube(gpu_windowed_det)), AdaptiveOpticsSim.detector_reference_cube(cpu_windowed_det); rtol=1f-5, atol=1f-4)
    @test isapprox(Array(AdaptiveOpticsSim.detector_signal_cube(gpu_windowed_det)), AdaptiveOpticsSim.detector_signal_cube(cpu_windowed_det); rtol=1f-5, atol=1f-4)
    @test isapprox(Array(AdaptiveOpticsSim.detector_read_cube(gpu_windowed_det)), AdaptiveOpticsSim.detector_read_cube(cpu_windowed_det); rtol=1f-5, atol=1f-4)
    @test AdaptiveOpticsSim.detector_read_times(gpu_windowed_det) == AdaptiveOpticsSim.detector_read_times(cpu_windowed_det)

    gsc_mask = array_backend(fill(one(T), 8, 8))
    gsc_basis = array_backend(reshape(T.(1:192), 8, 8, 3) .* T(1e-3))
    gsc_frame = array_backend(reshape(T.(1:64), 8, 8))
    gsc = GainSensingCamera(gsc_mask, gsc_basis; T=T)
    calibrate!(gsc, gsc_frame)
    optical_gains = compute_optical_gains!(gsc, gsc_frame)
    @test optical_gains isa array_backend
    @test all(isfinite, Array(optical_gains))

    lift_tel = Telescope(resolution=8, diameter=T(8), sampling_time=T(1e-3),
        central_obstruction=zero(T), T=T, backend=backend)
    lift_src = Source(band=:I, magnitude=zero(T), T=T)
    lift_det = Detector(noise=NoiseNone(), psf_sampling=1, T=T,
        backend=backend)
    lift_basis = array_backend(rand(MersenneTwister(29), T, 8, 8, 3))
    lift_diversity = array_backend(zeros(T, 8, 8))
    lift = LiFT(lift_tel, lift_src, lift_basis, lift_det;
        diversity_opd=lift_diversity, iterations=2, img_resolution=8,
        solve_mode=LiFTSolveAuto())
    lift_psf = compute_psf!(lift_tel, lift_src; zero_padding=1)
    lift_coeffs = reconstruct(lift, lift_psf, [1, 2])
    @test lift_coeffs isa array_backend
    @test all(isfinite, Array(lift_coeffs))
    return nothing
end

function run_optional_backend_plan_checks(::Type{AdaptiveOpticsSim.CUDABackendTag}, tel, backend)
    T = Float32
    array_backend = AdaptiveOpticsSim._resolve_array_backend(backend)
    sh = ShackHartmannWFS(tel; n_lenslets=4, mode=Diffractive(), T=T, backend=backend)
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
    @test AdaptiveOpticsSim.sh_sensing_execution_plan(AdaptiveOpticsSim.execution_style(sh.state.slopes), sh) isa AdaptiveOpticsSim.ShackHartmannWFSBatchedPlan
    @test AdaptiveOpticsSim.detector_execution_plan(typeof(AdaptiveOpticsSim.execution_style(det.state.frame)), typeof(det)) isa AdaptiveOpticsSim.DetectorDirectPlan
    @test AdaptiveOpticsSim.reduction_execution_plan(pyr.state.intensity) isa AdaptiveOpticsSim.DirectReductionPlan
    @test AdaptiveOpticsSim.atmospheric_field_execution_plan(
        AdaptiveOpticsSim.execution_style(first(geom_prop.state.slices).field.values),
        geom_prop.params.model,
    ) isa AdaptiveOpticsSim.GeometricFieldAsyncPlan
    @test AdaptiveOpticsSim.atmospheric_field_execution_plan(
        AdaptiveOpticsSim.execution_style(first(fresnel_prop.state.slices).field.values),
        fresnel_prop.params.model,
    ) isa AdaptiveOpticsSim.LayeredFresnelFieldAsyncPlan
    cpu_tel = Telescope(resolution=16, diameter=8.0f0, sampling_time=1.0f-3,
        central_obstruction=0.0f0, T=T, backend=CPUBackend())
    gpu_tel = Telescope(resolution=16, diameter=8.0f0, sampling_time=1.0f-3,
        central_obstruction=0.0f0, T=T, backend=backend)
    cpu_src = Source(band=:I, magnitude=0.0, T=T)
    gpu_src = Source(band=:I, magnitude=0.0, T=T)
    cpu_sh = ShackHartmannWFS(cpu_tel; n_lenslets=4, mode=Diffractive(), T=T, backend=CPUBackend())
    gpu_sh = ShackHartmannWFS(gpu_tel; n_lenslets=4, mode=Diffractive(), T=T, backend=backend)
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
    cpu_sh_stats = ShackHartmannWFS(cpu_tel; n_lenslets=4, mode=Diffractive(), T=T, backend=CPUBackend(),
        valid_subaperture_policy=FluxThresholdValidSubapertures(light_ratio=0.5f0))
    gpu_sh_stats = ShackHartmannWFS(gpu_tel; n_lenslets=4, mode=Diffractive(), T=T, backend=backend,
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
        bits=8, full_well=T(100), readout_window=AdaptiveOpticsSim.FrameWindow(2:5, 3:7), output_type=UInt16,
        sensor=HgCdTeAvalancheArraySensor(T=T),
        response_model=NullFrameResponse(),
        correction_model=generalized_correction,
        T=T,
        backend=CPUBackend())
    gpu_gen_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0),
        bits=8, full_well=T(100), readout_window=AdaptiveOpticsSim.FrameWindow(2:5, 3:7), output_type=UInt16,
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
        readout_window=AdaptiveOpticsSim.FrameWindow(2:3, 2:3),
        T=T,
        backend=CPUBackend())
    gpu_windowed_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0), binning=1,
        gain=T(1.0),
        sensor=HgCdTeAvalancheArraySensor(T=T, sampling_mode=CorrelatedDoubleSampling()),
        response_model=NullFrameResponse(),
        readout_window=AdaptiveOpticsSim.FrameWindow(2:3, 2:3),
        T=T,
        backend=backend)
    windowed_input = fill(T(10), 4, 4)
    cpu_windowed_frame = capture!(cpu_windowed_det, windowed_input; rng=MersenneTwister(15))
    gpu_windowed_frame = capture!(gpu_windowed_det, array_backend(copy(windowed_input)); rng=MersenneTwister(15))
    @test isapprox(Array(gpu_windowed_frame), cpu_windowed_frame; rtol=1f-5, atol=1f-4)
    @test isapprox(Array(AdaptiveOpticsSim.detector_reference_frame(gpu_windowed_det)), AdaptiveOpticsSim.detector_reference_frame(cpu_windowed_det); rtol=1f-5, atol=1f-4)
    @test isapprox(Array(AdaptiveOpticsSim.detector_signal_frame(gpu_windowed_det)), AdaptiveOpticsSim.detector_signal_frame(cpu_windowed_det); rtol=1f-5, atol=1f-4)
    @test isapprox(Array(AdaptiveOpticsSim.detector_combined_frame(gpu_windowed_det)), AdaptiveOpticsSim.detector_combined_frame(cpu_windowed_det); rtol=1f-5, atol=1f-4)
    @test isapprox(Array(AdaptiveOpticsSim.detector_reference_cube(gpu_windowed_det)), AdaptiveOpticsSim.detector_reference_cube(cpu_windowed_det); rtol=1f-5, atol=1f-4)
    @test isapprox(Array(AdaptiveOpticsSim.detector_signal_cube(gpu_windowed_det)), AdaptiveOpticsSim.detector_signal_cube(cpu_windowed_det); rtol=1f-5, atol=1f-4)
    @test isapprox(Array(AdaptiveOpticsSim.detector_read_cube(gpu_windowed_det)), AdaptiveOpticsSim.detector_read_cube(cpu_windowed_det); rtol=1f-5, atol=1f-4)
    @test AdaptiveOpticsSim.detector_read_times(gpu_windowed_det) == AdaptiveOpticsSim.detector_read_times(cpu_windowed_det)

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

    run_optional_counting_detector_parity(B, BackendArray)
    run_optional_avalanche_detector_parity(B, BackendArray)

    tel = Telescope(resolution=16, diameter=8.0f0, sampling_time=1.0f-3,
        central_obstruction=0.0f0, T=T, backend=selector)
    src = Source(band=:I, magnitude=0.0, T=T)
    run_optional_plane_product_checks(tel, src, selector, BackendArray, T)

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
    field = AdaptiveOpticsSim.propagate_atmosphere_field!(prop, atm, tel, src)
    @test field.values isa BackendArray
    intensity = AdaptiveOpticsSim.atmospheric_intensity!(prop, atm, tel, src)
    @test intensity isa BackendArray

    bundle = SpectralBundle(T[0.9 * wavelength(src), 1.1 * wavelength(src)], T[0.4, 0.6]; T=T)
    poly = with_spectrum(src, bundle)
    sh = ShackHartmannWFS(tel; n_lenslets=4, mode=Diffractive(), T=T, backend=selector)
    slopes = measure!(sh, tel, poly)
    @test slopes isa BackendArray

    science_src = Source(band=:K, magnitude=1.0, coordinates=(4.0, 90.0), T=T)
    split_dm = DeformableMirror(tel; n_act=4, influence_width=T(0.3), T=T,
        backend=selector)
    split_det = Detector(noise=NoiseNone(), integration_time=one(T), qe=one(T),
        binning=1, T=T, backend=selector)
    split_sim = AOSimulation(tel, src, atm, split_dm, sh;
        science_source=science_src)
    split_runtime = AdaptiveOpticsSim.ClosedLoopRuntime(
        split_sim,
        NullReconstructor();
        science_detector=split_det,
        outputs=RuntimeOutputRequirements(
            slopes=true,
            wfs_pixels=false,
            science_pixels=true,
        ),
        science_zero_padding=1,
        rng=MersenneTwister(17),
    )
    @test runtime_execution_plan(split_runtime) isa DeviceResidentExecutionPlan
    @test_throws InvalidConfiguration AdaptiveOpticsSim.ClosedLoopRuntime(
        split_sim,
        NullReconstructor();
        science_detector=split_det,
        execution_plan=CPUHILExecutionPlan(),
    )
    sense!(split_runtime)
    @test synchronize_runtime!(split_runtime) === split_runtime
    @test wfs_source(split_runtime) === src
    @test science_source(split_runtime) === science_src
    @test split_runtime.science_path isa AdaptiveOpticsSim.RepropagateScienceOpticalPath
    @test science_frame(split_runtime) isa BackendArray
    @test all(isfinite, Array(science_frame(split_runtime)))

    shared_detector_a = Detector(noise=NoiseNone(), integration_time=one(T),
        qe=one(T), binning=1, T=T, backend=selector)
    shared_detector_b = Detector(noise=NoiseNone(), integration_time=one(T),
        qe=one(T), binning=1, T=T, backend=selector)
    shared_arm = SharedOpticalArm(
        :shared_science,
        science_src;
        science_detectors=(shared_detector_a, shared_detector_b),
        science_zero_padding=1,
    )
    shared_runtime = SharedOpticalRuntime(split_runtime, shared_arm)
    sense!(shared_runtime)
    @test runtime_execution_plan(shared_runtime) isa DeviceResidentExecutionPlan
    @test all(frame -> frame isa BackendArray, science_frames(shared_arm))
    @test Array(science_frames(shared_arm)[1]) ≈
        Array(science_frames(shared_arm)[2]) atol=0 rtol=0
    @test synchronize_runtime!(shared_runtime) === shared_runtime

    run_optional_backend_plan_checks(B, tel, selector)
    run_optional_composite_optic_parity(B, BackendArray)
    run_optional_scalable_reconstructor_checks(T, selector, BackendArray)

    curv = CurvatureWFS(tel; pupil_samples=4, T=T, backend=selector)
    curv_slopes = measure!(curv, tel, src, atm)
    @test curv_slopes isa BackendArray

    single_control_loop = build_control_loop_scenario(
        SingleControlLoopConfig(
            name=:optional_backend_single,
            branch_label=:main,
            outputs=RuntimeOutputRequirements(slopes=true, wfs_pixels=true, science_pixels=false),
        ),
        build_optional_control_loop_branch(T, selector, :main; sensor=:sh, seed=21),
    )
    prepare!(single_control_loop)
    step!(single_control_loop)
    @test runtime_execution_plan(single_control_loop) isa DeviceResidentExecutionPlan
    @test synchronize_runtime!(single_control_loop) === single_control_loop
    single_runtime = AdaptiveOpticsSim.simulation_interface(single_control_loop).runtime
    @test all(storage -> typeof(AdaptiveOpticsSim.backend(storage)) === typeof(selector),
        AdaptiveOpticsSim.runtime_reconstructor_storage(single_runtime.reconstructor))
    @test AdaptiveOpticsSim.simulation_interface(single_control_loop) isa AdaptiveOpticsSim.SimulationInterface
    @test wfs_frame(single_control_loop) isa BackendArray
    @test control_loop_branch_labels(single_control_loop) == (:main,)

    grouped_control_loop = build_control_loop_scenario(
        GroupedControlLoopConfig(
            (:hi, :lo);
            name=:optional_backend_grouped,
            outputs=GroupedRuntimeOutputRequirements(wfs_frames=true, science_frames=false, wfs_stack=true, science_stack=false),
        ),
        build_optional_control_loop_branch(T, selector, :hi; sensor=:sh, seed=31),
        build_optional_control_loop_branch(T, selector, :lo; sensor=:sh, seed=32),
    )
    prepare!(grouped_control_loop)
    step!(grouped_control_loop)
    @test AdaptiveOpticsSim.simulation_interface(grouped_control_loop) isa AdaptiveOpticsSim.CompositeSimulationInterface
    @test grouped_wfs_stack(grouped_control_loop) isa BackendArray
    @test size(grouped_wfs_stack(grouped_control_loop), ndims(grouped_wfs_stack(grouped_control_loop))) == 2
    @test control_loop_branch_labels(grouped_control_loop) == (:hi, :lo)
    return nothing
end
