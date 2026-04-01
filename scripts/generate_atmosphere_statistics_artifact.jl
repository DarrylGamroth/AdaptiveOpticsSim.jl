using AdaptiveOpticsSim
using LinearAlgebra
using Random
using Statistics
using TOML

const OUTDIR = joinpath(@__DIR__, "..", "benchmarks", "results", "atmosphere")
const OUTFILE = joinpath(OUTDIR, "2026-04-01-phase1-pvp02.toml")
const MANIFEST = joinpath(OUTDIR, "manifest.toml")

normalized_correlation(a::AbstractMatrix, b::AbstractMatrix) =
    dot(vec(a), vec(b)) / sqrt(dot(vec(a), vec(a)) * dot(vec(b), vec(b)))

function ensemble_std(tel::Telescope, constructor, fractions; nsamp::Int=16, kwargs...)
    acc = 0.0
    for s in 1:nsamp
        atm = constructor(tel;
            r0=0.2,
            L0=25.0,
            fractional_cn2=fractions,
            wind_speed=fill(0.0, length(fractions)),
            wind_direction=fill(0.0, length(fractions)),
            altitude=fill(0.0, length(fractions)),
            kwargs...,
        )
        advance!(atm, tel; rng=MersenneTwister(s))
        acc += std(vec(Array(atm.state.opd)))
    end
    return acc / nsamp
end

function trajectory_std_windows(; seed::Integer=79, steps::Int=24)
    tel = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    atm = InfiniteMultiLayerAtmosphere(tel;
        r0=0.2,
        L0=25.0,
        fractional_cn2=[0.5, 0.5],
        wind_speed=[10.0, 6.0],
        wind_direction=[0.0, 90.0],
        altitude=[0.0, 5000.0],
        screen_resolution=33,
        stencil_size=35,
    )
    rng = MersenneTwister(seed)
    stds = Float64[]
    for _ in 1:steps
        advance!(atm, tel; rng=rng)
        push!(stds, std(vec(Array(atm.state.opd))))
    end
    midpoint = steps ÷ 2
    return (; early=mean(@view(stds[1:midpoint])), late=mean(@view(stds[midpoint + 1:end])))
end

function periodicity_metrics()
    tel = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    delta = tel.params.diameter / tel.params.resolution
    wind_speed_px = delta / tel.params.sampling_time
    finite = MultiLayerAtmosphere(tel;
        r0=0.2,
        L0=25.0,
        fractional_cn2=[1.0],
        wind_speed=[wind_speed_px],
        wind_direction=[0.0],
        altitude=[0.0],
    )
    infinite = InfiniteMultiLayerAtmosphere(tel;
        r0=0.2,
        L0=25.0,
        fractional_cn2=[1.0],
        wind_speed=[wind_speed_px],
        wind_direction=[0.0],
        altitude=[0.0],
        screen_resolution=33,
        stencil_size=35,
    )
    finite_rng = MersenneTwister(8)
    infinite_rng = MersenneTwister(8)
    advance!(finite, tel; rng=finite_rng)
    advance!(infinite, tel; rng=infinite_rng)
    finite_snapshot = copy(Array(finite.state.opd))
    infinite_snapshot = copy(Array(infinite.state.opd))
    period = AdaptiveOpticsSim.moving_layer_screen_resolution(tel.params.resolution)
    for _ in 1:period
        advance!(finite, tel; rng=finite_rng)
        advance!(infinite, tel; rng=infinite_rng)
    end
    infinite_period_corr = normalized_correlation(Array(infinite.state.opd), infinite_snapshot)
    finite_period_exact = Array(finite.state.opd) == finite_snapshot
    for _ in 1:(2 * period)
        advance!(finite, tel; rng=finite_rng)
        advance!(infinite, tel; rng=infinite_rng)
    end
    infinite_long_corr = normalized_correlation(Array(infinite.state.opd), infinite_snapshot)
    return (; period, finite_period_exact, infinite_period_corr, infinite_long_corr)
end

function build_report()
    tel = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    finite_single = ensemble_std(tel, MultiLayerAtmosphere, [1.0])
    infinite_single = ensemble_std(tel, InfiniteMultiLayerAtmosphere, [1.0];
        screen_resolution=33, stencil_size=35)
    infinite_equal = ensemble_std(tel, InfiniteMultiLayerAtmosphere, [0.5, 0.5];
        screen_resolution=33, stencil_size=35)
    infinite_uneven = ensemble_std(tel, InfiniteMultiLayerAtmosphere, [0.8, 0.2];
        screen_resolution=33, stencil_size=35)
    stationarity = trajectory_std_windows()
    periodicity = periodicity_metrics()

    return Dict(
        "artifact_id" => "ATM-STAT-2026-04-01",
        "generated_on" => "2026-04-01",
        "scope" => Dict(
            "finite_model" => "MultiLayerAtmosphere",
            "infinite_model" => "InfiniteMultiLayerAtmosphere",
            "resolution" => 16,
            "diameter_m" => 8.0,
            "sampling_time_s" => 1.0e-3,
            "r0_m" => 0.2,
            "L0_m" => 25.0,
        ),
        "ensemble_std" => Dict(
            "finite_single_layer" => finite_single,
            "infinite_single_layer" => infinite_single,
            "infinite_equal_two_layer" => infinite_equal,
            "infinite_uneven_two_layer" => infinite_uneven,
            "single_layer_ratio" => infinite_single / finite_single,
            "equal_ratio" => infinite_equal / infinite_single,
            "uneven_ratio" => infinite_uneven / infinite_single,
        ),
        "stationarity" => Dict(
            "early_window_std" => stationarity.early,
            "late_window_std" => stationarity.late,
            "late_to_early_ratio" => stationarity.late / stationarity.early,
        ),
        "periodicity" => Dict(
            "screen_period_steps" => periodicity.period,
            "finite_period_exact" => periodicity.finite_period_exact,
            "infinite_period_corr" => periodicity.infinite_period_corr,
            "infinite_long_corr" => periodicity.infinite_long_corr,
        ),
        "interpretation" => Dict(
            "finite_infinite_single_layer_close" => abs(infinite_single / finite_single - 1) <= 0.2,
            "windowed_stationarity_close" => abs(stationarity.late / stationarity.early - 1) <= 0.2,
            "finite_periodicity_exact" => periodicity.finite_period_exact,
            "infinite_nonperiodic_after_wrap" => periodicity.infinite_period_corr < 0.9,
            "infinite_nonperiodic_long_run" => periodicity.infinite_long_corr < 0.9,
        ),
    )
end

function update_manifest()
    manifest = isfile(MANIFEST) ? TOML.parsefile(MANIFEST) : Dict{String,Any}()
    entries = get!(manifest, "artifacts", Any[])
    push!(entries, Dict(
        "id" => "ATM-STAT-2026-04-01",
        "path" => "2026-04-01-phase1-pvp02.toml",
        "purpose" => "finite/infinite atmosphere statistics artifact for PVP-02",
    ))
    unique_entries = Dict{String,Any}()
    for entry in entries
        unique_entries[entry["id"]] = entry
    end
    manifest["artifacts"] = collect(values(unique_entries))
    open(MANIFEST, "w") do io
        TOML.print(io, manifest)
    end
end

mkpath(OUTDIR)
open(OUTFILE, "w") do io
    TOML.print(io, build_report())
end
update_manifest()
println("wrote ", OUTFILE)
