using AdaptiveOpticsSim
using Random
using Statistics
using LinearAlgebra
using TOML

const OUTDIR = joinpath(@__DIR__, "..", "benchmarks", "results", "platform")
const OUTFILE = joinpath(OUTDIR, "2026-04-13-multi-optic-hil.toml")
const MANIFEST = joinpath(OUTDIR, "manifest.toml")

struct StaticAtmosphere{A,B<:AbstractArrayBackend} <: AdaptiveOpticsSim.AbstractAtmosphere
    screen::A
end

AdaptiveOpticsSim.backend(::StaticAtmosphere{<:Any,B}) where {B} = B()
AdaptiveOpticsSim.advance!(atm::StaticAtmosphere, tel::Telescope, rng::AbstractRNG) = atm
AdaptiveOpticsSim.advance!(atm::StaticAtmosphere, tel::Telescope; rng::AbstractRNG=Random.default_rng()) = atm
function AdaptiveOpticsSim.propagate!(atm::StaticAtmosphere, tel::Telescope)
    copyto!(tel.state.opd, atm.screen)
    return tel
end

function deterministic_phase_screen(tel::Telescope, ::Type{T}) where {T<:AbstractFloat}
    host = zeros(T, tel.params.resolution, tel.params.resolution)
    host .*= Array(tel.state.pupil)
    return host
end

function StaticAtmosphere(tel::Telescope; T::Type{<:AbstractFloat}=Float32, backend::AbstractArrayBackend=CPUBackend())
    selector = AdaptiveOpticsSim.require_same_backend(tel, AdaptiveOpticsSim._resolve_backend_selector(backend))
    array_backend = AdaptiveOpticsSim._resolve_array_backend(selector)
    host = deterministic_phase_screen(tel, T)
    screen = array_backend{T}(undef, size(host)...)
    copyto!(screen, host)
    return StaticAtmosphere{typeof(screen),typeof(selector)}(screen)
end

function build_scenario(; seed::Integer=91)
    T = Float32
    tel = Telescope(resolution=16, diameter=T(8.0), sampling_time=T(1e-3),
        central_obstruction=T(0.0), T=T, backend=CPUBackend())
    src = Source(band=:I, magnitude=0.0, T=T)
    atm = KolmogorovAtmosphere(tel; r0=T(0.2), L0=T(25.0), T=T, backend=CPUBackend())
    tiptilt = TipTiltMirror(tel; scale=T(0.1), T=T, backend=CPUBackend(), label=:tiptilt)
    dm = DeformableMirror(tel; n_act=4, influence_width=T(0.3), T=T, backend=CPUBackend())
    optic = CompositeControllableOptic(:tiptilt => tiptilt, :dm => dm)
    wfs = ShackHartmann(tel; n_subap=4, mode=Diffractive(), T=T, backend=CPUBackend())
    det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0), binning=1, T=T, backend=CPUBackend())
    sim = AOSimulation(tel, src, atm, optic, wfs)
    scenario = build_runtime_scenario(
        SingleRuntimeConfig(name=:multi_optic_hil, branch_label=:main,
            products=RuntimeProductRequirements(slopes=true, wfs_pixels=true, science_pixels=false)),
        RuntimeBranch(:main, sim, NullReconstructor(); wfs_detector=det, rng=MersenneTwister(seed)),
    )
    prepare!(scenario)
    return scenario
end

@inline low_order_label(::Val{:tiptilt}) = :tiptilt
@inline low_order_label(::Val{:steering}) = :steering
@inline low_order_label(::Val{:focus}) = :focus

@inline wfs_family_label(::Val{:sh}) = :shack_hartmann
@inline wfs_family_label(::Val{:pyr}) = :pyramid
@inline wfs_family_label(::Val{:bio}) = :bioedge

function build_runtime_wfs(tel::Telescope, ::Val{:sh}; T::Type{<:AbstractFloat}=Float32,
    backend::AbstractArrayBackend=CPUBackend())
    return ShackHartmann(tel; n_subap=4, mode=Diffractive(), T=T, backend=backend)
end

function build_runtime_wfs(tel::Telescope, ::Val{:pyr}; T::Type{<:AbstractFloat}=Float32,
    backend::AbstractArrayBackend=CPUBackend())
    return PyramidWFS(tel; n_subap=4, modulation=T(1.0), mode=Diffractive(), T=T, backend=backend)
end

function build_runtime_wfs(tel::Telescope, ::Val{:bio}; T::Type{<:AbstractFloat}=Float32,
    backend::AbstractArrayBackend=CPUBackend())
    return BioEdgeWFS(tel; n_subap=4, modulation=T(1.0), mode=Diffractive(), T=T, backend=backend)
end

@inline build_runtime_detector_response(::Val{:null}, ::Type{T}, backend::AbstractArrayBackend) where {T<:AbstractFloat} = NullFrameResponse()
@inline build_runtime_detector_response(::Val{:gaussian}, ::Type{T}, backend::AbstractArrayBackend) where {T<:AbstractFloat} =
    default_response_model(CMOSSensor(T=T); T=T, backend=backend)

function build_runtime_detector(::Val{R}, ::Type{T}, backend::AbstractArrayBackend) where {R,T<:AbstractFloat}
    return Detector(
        noise=NoiseNone(),
        integration_time=T(1.0),
        qe=T(1.0),
        binning=1,
        sensor=CMOSSensor(T=T),
        response_model=build_runtime_detector_response(Val(R), T, backend),
        T=T,
        backend=backend,
    )
end

function build_behavior_optic(tel::Telescope, ::Val{:tiptilt}; T::Type{<:AbstractFloat}=Float32)
    return CompositeControllableOptic(
        :tiptilt => TipTiltMirror(tel; scale=T(0.1), T=T, backend=CPUBackend(), label=:tiptilt),
        :dm => DeformableMirror(tel; n_act=4, influence_width=T(0.3), T=T, backend=CPUBackend()),
    )
end

function build_behavior_optic(tel::Telescope, ::Val{:steering}; T::Type{<:AbstractFloat}=Float32)
    return CompositeControllableOptic(
        :steering => SteeringMirror(tel; scale=T(0.1), T=T, backend=CPUBackend(), label=:steering),
        :dm => DeformableMirror(tel; n_act=4, influence_width=T(0.3), T=T, backend=CPUBackend()),
    )
end

function build_behavior_optic(tel::Telescope, ::Val{:focus}; T::Type{<:AbstractFloat}=Float32)
    return CompositeControllableOptic(
        :focus => FocusStage(tel; scale=T(0.1), T=T, backend=CPUBackend(), label=:focus),
        :dm => DeformableMirror(tel; n_act=4, influence_width=T(0.3), T=T, backend=CPUBackend()),
    )
end

function build_behavior_scenario(::Val{K}; seed::Integer=91, wfs_family::Symbol=:sh,
    detector_profile::Symbol=:null) where {K}
    T = Float32
    tel = Telescope(resolution=16, diameter=T(8.0), sampling_time=T(1e-3),
        central_obstruction=T(0.0), T=T, backend=CPUBackend())
    src = Source(band=:I, magnitude=0.0, T=T)
    atm = StaticAtmosphere(tel; T=T, backend=CPUBackend())
    optic = build_behavior_optic(tel, Val(K); T=T)
    wfs = build_runtime_wfs(tel, Val(wfs_family); T=T, backend=CPUBackend())
    det = build_runtime_detector(Val(detector_profile), T, CPUBackend())
    sim = AOSimulation(tel, src, atm, optic, wfs)
    scenario = build_runtime_scenario(
        SingleRuntimeConfig(name=Symbol(:multi_optic_behavior_, K, :_, wfs_family, :_, detector_profile), branch_label=:main,
            products=RuntimeProductRequirements(slopes=true, wfs_pixels=true, science_pixels=false)),
        RuntimeBranch(:main, sim, NullReconstructor(); wfs_detector=det, rng=MersenneTwister(seed)),
    )
    prepare!(scenario)
    return scenario
end

function measure_behavior_response(::Val{K}, low_order_cmd::AbstractVector{<:AbstractFloat};
    wfs_family::Symbol=:sh, detector_profile::Symbol=:null) where {K}
    scenario = build_behavior_scenario(Val(K); wfs_family=wfs_family, detector_profile=detector_profile)
    label = low_order_label(Val(K))
    set_command!(scenario, NamedTuple{(label, :dm)}((collect(Float32.(low_order_cmd)), fill(Float32(0), 16))))
    sense!(scenario)
    slopes_vec = Array(slopes(scenario))
    frame = Array(wfs_frame(scenario))
    n = length(slopes_vec) ÷ 2
    return Dict(
        "command" => collect(Float32.(low_order_cmd)),
        "axis_1_mean" => mean(slopes_vec[1:n]),
        "axis_2_mean" => mean(slopes_vec[(n + 1):end]),
        "axis_1_norm" => norm(slopes_vec[1:n]),
        "axis_2_norm" => norm(slopes_vec[(n + 1):end]),
        "slopes" => slopes_vec,
        "wfs_frame" => frame,
        "frame_sum" => sum(frame),
    )
end

function response_delta_metrics(::Val{K}, low_order_cmd::AbstractVector{<:AbstractFloat};
    wfs_family::Symbol=:sh, detector_profile::Symbol=:null) where {K}
    active = measure_behavior_response(Val(K), low_order_cmd; wfs_family=wfs_family, detector_profile=detector_profile)
    baseline = measure_behavior_response(Val(K), zeros(Float32, length(low_order_cmd)); wfs_family=wfs_family, detector_profile=detector_profile)
    return Dict(
        "command" => collect(Float32.(low_order_cmd)),
        "slopes_delta_l2" => norm(active["slopes"] .- baseline["slopes"]),
        "frame_delta_l2" => norm(active["wfs_frame"] .- baseline["wfs_frame"]),
        "baseline_frame_sum" => baseline["frame_sum"],
        "active_frame_sum" => active["frame_sum"],
    )
end

function build_richer_runtime(::Val{S}; seed::Integer=91, wfs_family::Symbol=:sh,
    detector_profile::Symbol=:null) where {S}
    T = Float32
    tel = Telescope(resolution=16, diameter=T(8.0), sampling_time=T(1e-3),
        central_obstruction=T(0.0), T=T, backend=CPUBackend())
    src = Source(band=:I, magnitude=0.0, T=T)
    atm = StaticAtmosphere(tel; T=T, backend=CPUBackend())
    optic = if S === :tiptilt_focus_dm
        CompositeControllableOptic(
            :tiptilt => TipTiltMirror(tel; scale=T(0.1), T=T, backend=CPUBackend(), label=:tiptilt),
            :focus => FocusStage(tel; scale=T(0.1), T=T, backend=CPUBackend(), label=:focus),
            :dm => DeformableMirror(tel; n_act=4, influence_width=T(0.3), T=T, backend=CPUBackend()),
        )
    elseif S === :steering_focus_dm
        CompositeControllableOptic(
            :steering => SteeringMirror(tel; scale=T(0.1), T=T, backend=CPUBackend(), label=:steering),
            :focus => FocusStage(tel; scale=T(0.1), T=T, backend=CPUBackend(), label=:focus),
            :dm => DeformableMirror(tel; n_act=4, influence_width=T(0.3), T=T, backend=CPUBackend()),
        )
    else
        throw(InvalidConfiguration("unsupported richer runtime spec $(S)"))
    end
    wfs = build_runtime_wfs(tel, Val(wfs_family); T=T, backend=CPUBackend())
    det = build_runtime_detector(Val(detector_profile), T, CPUBackend())
    sim = AOSimulation(tel, src, atm, optic, wfs)
    scenario = build_runtime_scenario(
        SingleRuntimeConfig(name=Symbol(:richer_runtime_, S, :_, wfs_family, :_, detector_profile), branch_label=:main,
            products=RuntimeProductRequirements(slopes=true, wfs_pixels=true, science_pixels=false)),
        RuntimeBranch(:main, sim, NullReconstructor(); wfs_detector=det, rng=MersenneTwister(seed)),
    )
    prepare!(scenario)
    return scenario
end

runtime_snapshot(scenario) = Dict(
    "command" => copy(Array(command(scenario))),
    "slopes" => copy(Array(slopes(scenario))),
    "wfs_frame" => copy(Array(wfs_frame(scenario))),
)

function low_order_opd(::Val{K}, low_order_cmd::AbstractVector{<:AbstractFloat};
    dm_cmd::AbstractVector{<:AbstractFloat}=fill(Float32(0), 16), composite::Bool=false) where {K}
    T = Float32
    tel = Telescope(resolution=16, diameter=T(8.0), sampling_time=T(1e-3),
        central_obstruction=T(0.0), T=T, backend=CPUBackend())
    low_order = K === :tiptilt ? TipTiltMirror(tel; scale=T(0.1), T=T, backend=CPUBackend(), label=:tiptilt) :
        K === :steering ? SteeringMirror(tel; scale=T(0.1), T=T, backend=CPUBackend(), label=:steering) :
        FocusStage(tel; scale=T(0.1), T=T, backend=CPUBackend(), label=:focus)
    dm = DeformableMirror(tel; n_act=4, influence_width=T(0.3), T=T, backend=CPUBackend())
    fill!(tel.state.opd, zero(T))
    if composite
        optic = CompositeControllableOptic(low_order_label(Val(K)) => low_order, :dm => dm)
        set_command!(optic, NamedTuple{(low_order_label(Val(K)), :dm)}((collect(Float32.(low_order_cmd)), collect(Float32.(dm_cmd)))))
        apply!(optic, tel, DMAdditive())
    else
        set_command!(low_order, collect(Float32.(low_order_cmd)))
        set_command!(dm, collect(Float32.(dm_cmd)))
        apply!(low_order, tel, DMAdditive())
        apply!(dm, tel, DMAdditive())
    end
    return copy(Array(tel.state.opd))
end

function vector_stats(x)
    return Dict(
        "length" => length(x),
        "l2_norm" => norm(x),
        "mean" => mean(x),
        "maxabs" => maximum(abs, x),
        "sum" => sum(x),
    )
end

function frame_stats(frame)
    return Dict(
        "shape" => collect(size(frame)),
        "sum" => sum(frame),
        "mean" => mean(frame),
        "max" => maximum(frame),
        "min" => minimum(frame),
        "l2_norm" => norm(frame),
    )
end

function build_report()
    scenario_a = build_scenario(seed=91)
    scenario_b = build_scenario(seed=91)

    initial_tip = fill(Float32(0.0125), 2)
    initial_dm = fill(Float32(0.02), 16)
    updated_tip = fill(Float32(0.025), 2)

    initial_cmd = (; tiptilt=initial_tip, dm=initial_dm)
    update_cmd = (; tiptilt=updated_tip)

    set_command!(scenario_a, initial_cmd)
    set_command!(scenario_b, initial_cmd)
    sense!(scenario_a)
    sense!(scenario_b)

    initial_command = copy(Array(command(scenario_a)))
    initial_slopes_a = copy(Array(slopes(scenario_a)))
    initial_slopes_b = copy(Array(slopes(scenario_b)))
    initial_frame_a = copy(Array(wfs_frame(scenario_a)))
    initial_frame_b = copy(Array(wfs_frame(scenario_b)))

    update_command!(scenario_a, update_cmd)
    update_command!(scenario_b, update_cmd)
    sense!(scenario_a)
    sense!(scenario_b)

    updated_command = Array(command(scenario_a))
    updated_slopes_a = copy(Array(slopes(scenario_a)))
    updated_slopes_b = copy(Array(slopes(scenario_b)))
    updated_frame_a = copy(Array(wfs_frame(scenario_a)))
    updated_frame_b = copy(Array(wfs_frame(scenario_b)))

    layout = command_layout(scenario_a)
    segments = command_segments(layout)

    tip_plus = measure_behavior_response(Val(:tiptilt), Float32[0.0125, 0.0])
    tip_minus = measure_behavior_response(Val(:tiptilt), Float32[-0.0125, 0.0])
    tilt_plus = measure_behavior_response(Val(:tiptilt), Float32[0.0, 0.0125])
    tilt_minus = measure_behavior_response(Val(:tiptilt), Float32[0.0, -0.0125])
    tip_zero = measure_behavior_response(Val(:tiptilt), Float32[0.0, 0.0])
    steering_plus = measure_behavior_response(Val(:steering), Float32[0.0125, 0.0])
    steering_minus = measure_behavior_response(Val(:steering), Float32[-0.0125, 0.0])
    steering_zero = measure_behavior_response(Val(:steering), Float32[0.0, 0.0])
    focus_plus = measure_behavior_response(Val(:focus), Float32[0.0125])
    focus_minus = measure_behavior_response(Val(:focus), Float32[-0.0125])
    focus_zero = measure_behavior_response(Val(:focus), Float32[0.0])
    tip_seq = low_order_opd(Val(:tiptilt), Float32[0.01, -0.02]; dm_cmd=fill(Float32(0.03), 16), composite=false)
    tip_comp = low_order_opd(Val(:tiptilt), Float32[0.01, -0.02]; dm_cmd=fill(Float32(0.03), 16), composite=true)
    steer_seq = low_order_opd(Val(:steering), Float32[0.01, 0.02]; dm_cmd=fill(Float32(0.02), 16), composite=false)
    steer_comp = low_order_opd(Val(:steering), Float32[0.01, 0.02]; dm_cmd=fill(Float32(0.02), 16), composite=true)
    focus_seq = low_order_opd(Val(:focus), Float32[0.02]; dm_cmd=fill(Float32(0.01), 16), composite=false)
    focus_comp = low_order_opd(Val(:focus), Float32[0.02]; dm_cmd=fill(Float32(0.01), 16), composite=true)
    tip_a = low_order_opd(Val(:tiptilt), Float32[0.005, 0.0])
    tip_b = low_order_opd(Val(:tiptilt), Float32[0.0, 0.005])
    tip_ab = low_order_opd(Val(:tiptilt), Float32[0.005, 0.005])
    tip_2a = low_order_opd(Val(:tiptilt), Float32[0.01, 0.0])
    steer_a = low_order_opd(Val(:steering), Float32[0.005, 0.0])
    steer_2a = low_order_opd(Val(:steering), Float32[0.01, 0.0])
    focus_a = low_order_opd(Val(:focus), Float32[0.005])
    focus_2a = low_order_opd(Val(:focus), Float32[0.01])
    pyr_tip_plus = measure_behavior_response(Val(:tiptilt), Float32[0.0125, 0.0]; wfs_family=:pyr)
    pyr_tip_minus = measure_behavior_response(Val(:tiptilt), Float32[-0.0125, 0.0]; wfs_family=:pyr)
    pyr_tip_zero = measure_behavior_response(Val(:tiptilt), Float32[0.0, 0.0]; wfs_family=:pyr)
    bio_tip_plus = measure_behavior_response(Val(:tiptilt), Float32[0.0125, 0.0]; wfs_family=:bio)
    bio_tip_minus = measure_behavior_response(Val(:tiptilt), Float32[-0.0125, 0.0]; wfs_family=:bio)
    sh_sweep = map(cmd -> response_delta_metrics(Val(:tiptilt), Float32[cmd, 0.0]), Float32[0.0125, 0.025, 0.05])
    pyr_sweep = map(cmd -> response_delta_metrics(Val(:tiptilt), Float32[cmd, 0.0]; wfs_family=:pyr), Float32[0.0125, 0.025, 0.05])
    bio_sweep = map(cmd -> response_delta_metrics(Val(:tiptilt), Float32[cmd, 0.0]; wfs_family=:bio), Float32[0.0125, 0.025, 0.05])
    sh_det_plus = response_delta_metrics(Val(:tiptilt), Float32[0.0125, 0.0]; detector_profile=:gaussian)
    sh_det_minus = response_delta_metrics(Val(:tiptilt), Float32[-0.0125, 0.0]; detector_profile=:gaussian)
    pyr_det_plus = response_delta_metrics(Val(:tiptilt), Float32[0.0125, 0.0]; wfs_family=:pyr, detector_profile=:gaussian)
    pyr_det_minus = response_delta_metrics(Val(:tiptilt), Float32[-0.0125, 0.0]; wfs_family=:pyr, detector_profile=:gaussian)
    bio_det_plus = response_delta_metrics(Val(:tiptilt), Float32[0.0125, 0.0]; wfs_family=:bio, detector_profile=:gaussian)
    bio_det_minus = response_delta_metrics(Val(:tiptilt), Float32[-0.0125, 0.0]; wfs_family=:bio, detector_profile=:gaussian)
    richer = build_richer_runtime(Val(:tiptilt_focus_dm))
    richer_initial = (; tiptilt=Float32[0.01, -0.005], focus=Float32[0.015], dm=fill(Float32(0.02), 16))
    richer_tip_update = (; tiptilt=Float32[-0.02, 0.01])
    richer_focus_update = (; focus=Float32[-0.01])
    set_command!(richer, richer_initial)
    sense!(richer)
    richer_step_1 = runtime_snapshot(richer)
    update_command!(richer, richer_tip_update)
    sense!(richer)
    richer_step_2 = runtime_snapshot(richer)
    update_command!(richer, richer_focus_update)
    sense!(richer)
    richer_step_3 = runtime_snapshot(richer)

    return Dict(
        "artifact_id" => "MULTI-OPTIC-HIL-2026-04-13",
        "generated_on" => "2026-04-13",
        "scope" => Dict(
            "artifact_kind" => "multi_optic_hil_validation",
            "backend" => "cpu",
            "families" => ["runtime_orchestration", "composite_controllable_optic", "tiptilt_plus_dm", "shack_hartmann"],
        ),
        "surface" => Dict(
            "command_segment_labels" => collect(String.(command_segment_labels(layout))),
            "command_segment_lengths" => [length(command_segment_range(seg)) for seg in segments],
            "command_length" => length(updated_command),
            "wfs_family" => "shack_hartmann_diffractive",
            "science_pixels_enabled" => false,
        ),
        "determinism" => Dict(
            "initial_slopes_exact" => initial_slopes_a == initial_slopes_b,
            "initial_wfs_frame_exact" => initial_frame_a == initial_frame_b,
            "updated_slopes_exact" => updated_slopes_a == updated_slopes_b,
            "updated_wfs_frame_exact" => updated_frame_a == updated_frame_b,
        ),
        "initial_state" => Dict(
            "command" => Dict(
                "tiptilt_mean" => mean(initial_tip),
                "dm_mean" => mean(initial_dm),
                "packed_sum" => sum(initial_command),
            ),
            "slopes" => vector_stats(initial_slopes_a),
            "wfs_frame" => frame_stats(initial_frame_a),
        ),
        "updated_state" => Dict(
            "command" => Dict(
                "tiptilt_mean" => mean(updated_command[1:2]),
                "dm_mean" => mean(updated_command[3:end]),
                "packed_sum" => sum(updated_command),
                "dm_segment_unchanged" => updated_command[3:end] == initial_dm,
            ),
            "slopes" => vector_stats(updated_slopes_a),
            "wfs_frame" => frame_stats(updated_frame_a),
        ),
        "delta" => Dict(
            "slopes_delta_l2" => norm(updated_slopes_a .- initial_slopes_a),
            "wfs_frame_delta_l2" => norm(updated_frame_a .- initial_frame_a),
        ),
        "behavior" => Dict(
            "static_surface" => "zero_phase_static_atmosphere",
            "tip_plus" => Dict(
                "axis_1_mean" => tip_plus["axis_1_mean"],
                "axis_2_mean" => tip_plus["axis_2_mean"],
                "axis_1_norm" => tip_plus["axis_1_norm"],
                "axis_2_norm" => tip_plus["axis_2_norm"],
                "frame_sum" => tip_plus["frame_sum"],
            ),
            "tip_minus" => Dict(
                "axis_1_mean" => tip_minus["axis_1_mean"],
                "axis_2_mean" => tip_minus["axis_2_mean"],
                "axis_1_norm" => tip_minus["axis_1_norm"],
                "axis_2_norm" => tip_minus["axis_2_norm"],
                "frame_sum" => tip_minus["frame_sum"],
            ),
            "tilt_plus" => Dict(
                "axis_1_mean" => tilt_plus["axis_1_mean"],
                "axis_2_mean" => tilt_plus["axis_2_mean"],
                "axis_1_norm" => tilt_plus["axis_1_norm"],
                "axis_2_norm" => tilt_plus["axis_2_norm"],
            ),
            "tilt_minus" => Dict(
                "axis_1_mean" => tilt_minus["axis_1_mean"],
                "axis_2_mean" => tilt_minus["axis_2_mean"],
                "axis_1_norm" => tilt_minus["axis_1_norm"],
                "axis_2_norm" => tilt_minus["axis_2_norm"],
            ),
            "tip_antisymmetry_l2" => norm(tip_plus["slopes"] .+ tip_minus["slopes"]),
            "tilt_antisymmetry_l2" => norm(tilt_plus["slopes"] .+ tilt_minus["slopes"]),
            "tip_vs_tilt_l2" => norm(tip_plus["slopes"] .- tilt_plus["slopes"]),
            "tip_dominant_axis" => 2,
            "tilt_dominant_axis" => 1,
            "tip_axis_ratio" => tip_plus["axis_2_norm"] / max(tip_plus["axis_1_norm"], eps(Float32)),
            "tilt_axis_ratio" => tilt_plus["axis_1_norm"] / max(tilt_plus["axis_2_norm"], eps(Float32)),
            "zero_slopes_l2" => norm(tip_zero["slopes"]),
        ),
        "steering_behavior" => Dict(
            "steering_plus" => Dict(
                "axis_1_mean" => steering_plus["axis_1_mean"],
                "axis_2_mean" => steering_plus["axis_2_mean"],
                "axis_1_norm" => steering_plus["axis_1_norm"],
                "axis_2_norm" => steering_plus["axis_2_norm"],
                "frame_sum" => steering_plus["frame_sum"],
            ),
            "steering_minus" => Dict(
                "axis_1_mean" => steering_minus["axis_1_mean"],
                "axis_2_mean" => steering_minus["axis_2_mean"],
                "axis_1_norm" => steering_minus["axis_1_norm"],
                "axis_2_norm" => steering_minus["axis_2_norm"],
                "frame_sum" => steering_minus["frame_sum"],
            ),
            "antisymmetry_l2" => norm(steering_plus["slopes"] .+ steering_minus["slopes"]),
            "zero_slopes_l2" => norm(steering_zero["slopes"]),
            "matches_tiptilt_exactly" => steering_plus["slopes"] == tip_plus["slopes"],
        ),
        "focus_behavior" => Dict(
            "focus_plus" => Dict(
                "axis_1_mean" => focus_plus["axis_1_mean"],
                "axis_2_mean" => focus_plus["axis_2_mean"],
                "axis_1_norm" => focus_plus["axis_1_norm"],
                "axis_2_norm" => focus_plus["axis_2_norm"],
                "frame_sum" => focus_plus["frame_sum"],
            ),
            "focus_minus" => Dict(
                "axis_1_mean" => focus_minus["axis_1_mean"],
                "axis_2_mean" => focus_minus["axis_2_mean"],
                "axis_1_norm" => focus_minus["axis_1_norm"],
                "axis_2_norm" => focus_minus["axis_2_norm"],
                "frame_sum" => focus_minus["frame_sum"],
            ),
            "antisymmetry_l2" => norm(focus_plus["slopes"] .+ focus_minus["slopes"]),
            "zero_slopes_l2" => norm(focus_zero["slopes"]),
            "axis_balance_ratio" => focus_plus["axis_1_norm"] / max(focus_plus["axis_2_norm"], eps(Float32)),
        ),
        "pyramid_behavior" => Dict(
            "tip_plus_axis_ratio" => pyr_tip_plus["axis_2_norm"] / max(pyr_tip_plus["axis_1_norm"], eps(Float32)),
            "tip_antisymmetry_l2" => norm(pyr_tip_plus["slopes"] .+ pyr_tip_minus["slopes"]),
            "zero_slopes_l2" => norm(pyr_tip_zero["slopes"]),
            "tip_plus_frame_sum" => pyr_tip_plus["frame_sum"],
            "tip_minus_frame_sum" => pyr_tip_minus["frame_sum"],
        ),
        "bioedge_behavior" => Dict(
            "tip_plus_slopes_l2" => norm(bio_tip_plus["slopes"]),
            "tip_minus_slopes_l2" => norm(bio_tip_minus["slopes"]),
            "tip_frame_sum_delta" => abs(bio_tip_plus["frame_sum"] - bio_tip_minus["frame_sum"]),
        ),
        "sensitivity_envelope" => Dict(
            "shack_hartmann" => Dict(
                "commands" => [entry["command"][1] for entry in sh_sweep],
                "slopes_delta_l2" => [entry["slopes_delta_l2"] for entry in sh_sweep],
                "frame_delta_l2" => [entry["frame_delta_l2"] for entry in sh_sweep],
            ),
            "pyramid" => Dict(
                "commands" => [entry["command"][1] for entry in pyr_sweep],
                "slopes_delta_l2" => [entry["slopes_delta_l2"] for entry in pyr_sweep],
                "frame_delta_l2" => [entry["frame_delta_l2"] for entry in pyr_sweep],
            ),
            "bioedge" => Dict(
                "commands" => [entry["command"][1] for entry in bio_sweep],
                "frame_delta_l2" => [entry["frame_delta_l2"] for entry in bio_sweep],
            ),
        ),
        "detector_in_loop" => Dict(
            "shack_hartmann" => Dict(
                "plus_frame_delta_l2" => sh_det_plus["frame_delta_l2"],
                "minus_frame_delta_l2" => sh_det_minus["frame_delta_l2"],
                "frame_sum_delta" => abs(sh_det_plus["active_frame_sum"] - sh_det_minus["active_frame_sum"]),
            ),
            "pyramid" => Dict(
                "plus_frame_delta_l2" => pyr_det_plus["frame_delta_l2"],
                "minus_frame_delta_l2" => pyr_det_minus["frame_delta_l2"],
                "frame_sum_delta" => abs(pyr_det_plus["active_frame_sum"] - pyr_det_minus["active_frame_sum"]),
            ),
            "bioedge" => Dict(
                "plus_frame_delta_l2" => bio_det_plus["frame_delta_l2"],
                "minus_frame_delta_l2" => bio_det_minus["frame_delta_l2"],
                "frame_sum_delta" => abs(bio_det_plus["active_frame_sum"] - bio_det_minus["active_frame_sum"]),
            ),
        ),
        "stateful_sequence" => Dict(
            "surface" => "tiptilt_focus_dm_shack_hartmann",
            "step_1_command_sum" => sum(richer_step_1["command"]),
            "step_2_command_sum" => sum(richer_step_2["command"]),
            "step_3_command_sum" => sum(richer_step_3["command"]),
            "step_1_slopes_l2" => norm(richer_step_1["slopes"]),
            "step_2_slopes_l2" => norm(richer_step_2["slopes"]),
            "step_3_slopes_l2" => norm(richer_step_3["slopes"]),
            "step_1_frame_l2" => norm(richer_step_1["wfs_frame"]),
            "step_2_frame_l2" => norm(richer_step_2["wfs_frame"]),
            "step_3_frame_l2" => norm(richer_step_3["wfs_frame"]),
        ),
        "composition" => Dict(
            "tiptilt_dm_opd_l2" => norm(tip_seq .- tip_comp),
            "steering_dm_opd_l2" => norm(steer_seq .- steer_comp),
            "focus_dm_opd_l2" => norm(focus_seq .- focus_comp),
        ),
        "opd_linearity" => Dict(
            "tiptilt_additivity_l2" => norm(tip_ab .- (tip_a .+ tip_b)),
            "tiptilt_scaling_l2" => norm(tip_2a .- (2 .* tip_a)),
            "steering_scaling_l2" => norm(steer_2a .- (2 .* steer_a)),
            "focus_scaling_l2" => norm(focus_2a .- (2 .* focus_a)),
        ),
    )
end

function update_manifest!(artifact_path::AbstractString)
    mkpath(dirname(MANIFEST))
    manifest = isfile(MANIFEST) ? TOML.parsefile(MANIFEST) : Dict{String,Any}()
    artifacts = get!(manifest, "artifacts", Any[])
    kept = Any[item for item in artifacts if get(item, "id", "") != "MULTI-OPTIC-HIL-2026-04-13"]
    push!(kept, Dict(
        "purpose" => "multi-optic HIL validation artifact for low-order composite runtime self-check surfaces",
        "id" => "MULTI-OPTIC-HIL-2026-04-13",
        "path" => basename(artifact_path),
    ))
    manifest["artifacts"] = kept
    open(MANIFEST, "w") do io
        TOML.print(io, manifest)
    end
end

function main()
    mkpath(OUTDIR)
    report = build_report()
    open(OUTFILE, "w") do io
        TOML.print(io, report)
    end
    update_manifest!(OUTFILE)
    println(OUTFILE)
end

main()
