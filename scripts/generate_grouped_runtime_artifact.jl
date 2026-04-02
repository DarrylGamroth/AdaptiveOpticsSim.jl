using AdaptiveOpticsSim
using TOML

include(joinpath(@__DIR__, "profile_multi_source_multi_wfs_runtime.jl"))

const OUTDIR = joinpath(@__DIR__, "..", "benchmarks", "results", "grouped")
const OUTFILE = joinpath(OUTDIR, "2026-04-01-gr.toml")
const MANIFEST = joinpath(OUTDIR, "manifest.toml")

function build_report()
    medium = run_profile(; backend_name="cpu", scale_name="medium", samples=12, warmup=3)
    representative = run_profile(; backend_name="cpu", scale_name="representative", samples=6, warmup=2)
    return Dict(
        "artifact_id" => "GROUPED-VAL-2026-04-01",
        "generated_on" => "2026-04-01",
        "scope" => Dict(
            "backend" => "cpu",
            "artifact_kind" => "grouped_runtime_validation",
            "families" => ["grouped_runtime", "multi_source_multi_wfs"],
        ),
        "cases" => Dict(
            "medium" => Dict(String(k) => v for (k, v) in pairs(medium)),
            "representative" => Dict(String(k) => v for (k, v) in pairs(representative)),
        ),
    )
end

function update_manifest!(artifact_path::AbstractString)
    mkpath(dirname(MANIFEST))
    manifest = isfile(MANIFEST) ? TOML.parsefile(MANIFEST) : Dict{String,Any}()
    artifacts = get!(manifest, "artifacts", Any[])
    kept = Any[item for item in artifacts if get(item, "id", "") != "GROUPED-VAL-2026-04-01"]
    push!(kept, Dict(
        "purpose" => "grouped runtime validation artifact for GR-4",
        "id" => "GROUPED-VAL-2026-04-01",
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
