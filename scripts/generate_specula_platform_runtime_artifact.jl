using TOML

include(joinpath(@__DIR__, "profile_multi_source_multi_wfs_runtime.jl"))

const OUTDIR = joinpath(@__DIR__, "..", "benchmarks", "results", "platform")
const OUTFILE = joinpath(OUTDIR, "2026-04-02-phase2-psp07.toml")
const MANIFEST = joinpath(OUTDIR, "manifest.toml")

function build_report()
    medium = run_profile(; backend_name="cpu", scale_name="medium", samples=12, warmup=3)
    representative = run_profile(; backend_name="cpu", scale_name="representative", samples=6, warmup=2)
    return Dict(
        "artifact_id" => "SPECULA-PLATFORM-2026-04-02",
        "generated_on" => "2026-04-02",
        "scope" => Dict(
            "artifact_kind" => "specula_informed_platform_runtime",
            "backend" => "cpu",
            "reference_mode" => "specula_informed_julia_native",
            "families" => ["grouped_runtime", "multi_source_multi_wfs", "composite_runtime"],
        ),
        "specula_inspiration" => Dict(
            "processing_object_tests" => [
                "test_mmse_reconstructor.py::TestMMSEReconstructor::test_mmse_multiple_wfs",
                "test_sprint.py",
                "test_sprint_pyr.py",
            ],
            "contract_statement" => "The scenario measures Julia-native grouped orchestration and heterogeneous multi-WFS runtime behavior on maintained exported surfaces inspired by SPECULA platform composition patterns.",
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
    kept = Any[item for item in artifacts if get(item, "id", "") != "SPECULA-PLATFORM-2026-04-02"]
    push!(kept, Dict(
        "purpose" => "SPECULA-informed platform runtime validation artifact for PSP-07",
        "id" => "SPECULA-PLATFORM-2026-04-02",
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
