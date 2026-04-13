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

function build_behavior_scenario(; seed::Integer=91)
    T = Float32
    tel = Telescope(resolution=16, diameter=T(8.0), sampling_time=T(1e-3),
        central_obstruction=T(0.0), T=T, backend=CPUBackend())
    src = Source(band=:I, magnitude=0.0, T=T)
    atm = StaticAtmosphere(tel; T=T, backend=CPUBackend())
    tiptilt = TipTiltMirror(tel; scale=T(0.1), T=T, backend=CPUBackend(), label=:tiptilt)
    dm = DeformableMirror(tel; n_act=4, influence_width=T(0.3), T=T, backend=CPUBackend())
    optic = CompositeControllableOptic(:tiptilt => tiptilt, :dm => dm)
    wfs = ShackHartmann(tel; n_subap=4, mode=Diffractive(), T=T, backend=CPUBackend())
    det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0), binning=1, T=T, backend=CPUBackend())
    sim = AOSimulation(tel, src, atm, optic, wfs)
    scenario = build_runtime_scenario(
        SingleRuntimeConfig(name=:multi_optic_behavior, branch_label=:main,
            products=RuntimeProductRequirements(slopes=true, wfs_pixels=true, science_pixels=false)),
        RuntimeBranch(:main, sim, NullReconstructor(); wfs_detector=det, rng=MersenneTwister(seed)),
    )
    prepare!(scenario)
    return scenario
end

function measure_behavior_response(tiptilt_cmd::AbstractVector{<:AbstractFloat})
    scenario = build_behavior_scenario()
    set_command!(scenario, (; tiptilt=tiptilt_cmd, dm=fill(Float32(0), 16)))
    sense!(scenario)
    slopes_vec = Array(slopes(scenario))
    frame = Array(wfs_frame(scenario))
    n = length(slopes_vec) ÷ 2
    return Dict(
        "tiptilt_command" => collect(tiptilt_cmd),
        "axis_1_mean" => mean(slopes_vec[1:n]),
        "axis_2_mean" => mean(slopes_vec[(n + 1):end]),
        "axis_1_norm" => norm(slopes_vec[1:n]),
        "axis_2_norm" => norm(slopes_vec[(n + 1):end]),
        "slopes" => slopes_vec,
        "wfs_frame" => frame,
    )
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

    tip_plus = measure_behavior_response(Float32[0.0125, 0.0])
    tip_minus = measure_behavior_response(Float32[-0.0125, 0.0])
    tilt_plus = measure_behavior_response(Float32[0.0, 0.0125])
    tilt_minus = measure_behavior_response(Float32[0.0, -0.0125])

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
            ),
            "tip_minus" => Dict(
                "axis_1_mean" => tip_minus["axis_1_mean"],
                "axis_2_mean" => tip_minus["axis_2_mean"],
                "axis_1_norm" => tip_minus["axis_1_norm"],
                "axis_2_norm" => tip_minus["axis_2_norm"],
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
        ),
    )
end

function update_manifest!(artifact_path::AbstractString)
    mkpath(dirname(MANIFEST))
    manifest = isfile(MANIFEST) ? TOML.parsefile(MANIFEST) : Dict{String,Any}()
    artifacts = get!(manifest, "artifacts", Any[])
    kept = Any[item for item in artifacts if get(item, "id", "") != "MULTI-OPTIC-HIL-2026-04-13"]
    push!(kept, Dict(
        "purpose" => "multi-optic HIL validation artifact for tiptilt + dm runtime surface",
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
