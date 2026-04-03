using TOML
import AMDGPU

include(joinpath(@__DIR__, "profile_platform_runtime.jl"))

const OUTDIR = joinpath(@__DIR__, "..", "benchmarks", "results", "platform")
const OUTFILE = joinpath(OUTDIR, "2026-04-03-phase5-psp15.toml")
const MANIFEST = joinpath(OUTDIR, "manifest.toml")

function build_report()
    cpu_medium = run_profile(; backend_name="cpu", scale_name="medium", samples=12, warmup=3)
    cpu_representative = run_profile(; backend_name="cpu", scale_name="representative", samples=6, warmup=2)
    amdgpu_medium = run_profile(; backend_name="amdgpu", scale_name="medium", samples=12, warmup=3)
    return Dict(
        "artifact_id" => "PLATFORM-ORCH-2026-04-03",
        "generated_on" => "2026-04-03",
        "scope" => Dict(
            "artifact_kind" => "platform_orchestration_runtime",
            "backends" => ["cpu", "amdgpu", "cuda"],
            "families" => ["single_platform_runtime", "compatible_grouped_runtime", "mixed_grouped_runtime"],
            "notes" => "CUDA result is recorded from the maintained spiders host run because this local environment does not provide a CUDA device.",
        ),
        "cases" => Dict(
            "cpu_medium" => Dict(String(k) => v for (k, v) in pairs(cpu_medium)),
            "cpu_representative" => Dict(String(k) => v for (k, v) in pairs(cpu_representative)),
            "amdgpu_medium" => Dict(String(k) => v for (k, v) in pairs(amdgpu_medium)),
            "cuda_medium" => Dict{String,Any}(),
        ),
    )
end

function update_manifest!(artifact_path::AbstractString)
    mkpath(dirname(MANIFEST))
    manifest = isfile(MANIFEST) ? TOML.parsefile(MANIFEST) : Dict{String,Any}()
    artifacts = get!(manifest, "artifacts", Any[])
    kept = Any[item for item in artifacts if get(item, "id", "") != "PLATFORM-ORCH-2026-04-03"]
    push!(kept, Dict(
        "purpose" => "platform orchestration runtime artifact for PSP-15",
        "id" => "PLATFORM-ORCH-2026-04-03",
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
