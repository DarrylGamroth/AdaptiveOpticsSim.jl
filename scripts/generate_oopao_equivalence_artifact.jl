using LinearAlgebra
using TOML

using AdaptiveOpticsSim

include(joinpath(@__DIR__, "..", "test", "reference_harness.jl"))

const OUTDIR = joinpath(@__DIR__, "..", "benchmarks", "results", "equivalence")
const OUTFILE = joinpath(OUTDIR, "2026-04-09-oopao-production-equivalence.toml")
const MANIFEST = joinpath(OUTDIR, "manifest.toml")
const CASE_IDS = (
    "pyramid_diffractive_ramp",
    "bioedge_diffractive_ramp",
)

function case_by_id(bundle::ReferenceBundle, case_id::AbstractString)
    for case in bundle.cases
        case.id == case_id && return case
    end
    throw(InvalidConfiguration("reference case '$case_id' not found in bundle '$(bundle.root)'"))
end

function finite_metrics(actual, expected)
    diff = Float64.(actual .- expected)
    expected_f = Float64.(expected)
    finite_mask = .!(isnan.(expected_f) .| isnan.(diff))
    if !any(finite_mask)
        return Dict(
            "max_abs_error" => 0.0,
            "max_rel_error" => 0.0,
            "l2_rel_error" => 0.0,
            "actual_norm" => 0.0,
            "expected_norm" => 0.0,
        )
    end
    adiff = diff[finite_mask]
    aexp = expected_f[finite_mask]
    expected_norm = norm(aexp)
    denom = max.(abs.(aexp), eps(Float64))
    return Dict(
        "max_abs_error" => maximum(abs.(adiff)),
        "max_rel_error" => maximum(abs.(adiff) ./ denom),
        "l2_rel_error" => norm(adiff) / max(expected_norm, eps(Float64)),
        "actual_norm" => norm(Float64.(actual[finite_mask])),
        "expected_norm" => expected_norm,
    )
end

function case_report(case::ReferenceCase)
    expected = load_reference_array(case)
    actual = adapt_reference_actual(case, compute_reference_actual(case))
    verdict = validate_reference_case(case)
    metrics = finite_metrics(actual, expected)
    return Dict(
        "id" => case.id,
        "baseline" => String(case.baseline),
        "kind" => String(case.kind),
        "shape" => collect(size(expected)),
        "atol" => case.atol,
        "rtol" => case.rtol,
        "within_tolerance" => verdict.ok,
        "max_abs_error" => metrics["max_abs_error"],
        "max_rel_error" => metrics["max_rel_error"],
        "l2_rel_error" => metrics["l2_rel_error"],
        "actual_norm" => metrics["actual_norm"],
        "expected_norm" => metrics["expected_norm"],
    )
end

function build_report()
    bundle = load_reference_bundle()
    cases = Dict{String,Any}()
    for case_id in CASE_IDS
        case = case_by_id(bundle, case_id)
        cases[case_id] = case_report(case)
    end
    return Dict(
        "artifact_id" => "OOPAO-EQUIV-2026-04-09",
        "generated_on" => "2026-04-09",
        "scope" => Dict(
            "artifact_kind" => "oopao_external_equivalence",
            "baseline" => "oopao",
            "reference_bundle_root" => bundle.root,
            "surfaces" => [
                "pyramid_diffractive_ramp",
                "bioedge_diffractive_ramp",
            ],
            "interpretation" => "deterministic frozen OOPAO reference surfaces intended to complement the HEART Shack-Hartmann equivalence baseline",
        ),
        "cases" => cases,
        "all_cases_within_tolerance" => all(cases[id]["within_tolerance"] for id in keys(cases)),
    )
end

function update_manifest!(artifact_path::AbstractString)
    mkpath(dirname(MANIFEST))
    manifest = isfile(MANIFEST) ? TOML.parsefile(MANIFEST) : Dict{String,Any}()
    artifacts = get!(manifest, "artifacts", Any[])
    kept = Any[item for item in artifacts if get(item, "id", "") != "OOPAO-EQUIV-2026-04-09"]
    push!(kept, Dict(
        "purpose" => "production external-equivalence artifact for frozen OOPAO pyramid and bioedge surfaces",
        "id" => "OOPAO-EQUIV-2026-04-09",
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
