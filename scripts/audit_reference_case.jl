using TOML
using LinearAlgebra
using Statistics

using AdaptiveOpticsSim

include(joinpath(@__DIR__, "..", "test", "reference_harness.jl"))

function summarize_diff(actual, expected)
    diff = actual .- expected
    rel = similar(diff, Float64)
    @. rel = abs(diff) / max(abs(expected), eps(Float64))
    println("shape: ", size(actual))
    println("maxabs: ", maximum(abs, diff))
    println("maxrel: ", maximum(rel))
    println("rms: ", sqrt(sum(abs2, diff) / length(diff)))

    if ndims(actual) == 2
        println("\nper-column maxabs/maxrel:")
        for j in axes(actual, 2)
            println(
                "  col ", j, ": ",
                maximum(abs.(diff[:, j])),
                " / ",
                maximum(rel[:, j]),
            )
        end

        println("\nper-row maxabs/maxrel:")
        for i in axes(actual, 1)
            println(
                "  row ", i, ": ",
                maximum(abs.(diff[i, :])),
                " / ",
                maximum(rel[i, :]),
            )
        end
    end
end

function main()
    root = length(ARGS) >= 1 ? ARGS[1] : default_reference_root()
    case_id = length(ARGS) >= 2 ? ARGS[2] : "gsc_atmosphere_replay_trace_bounded"

    manifest = TOML.parsefile(reference_manifest_path(root))
    raw_cases = get(manifest, "cases", Dict{String,Any}())
    haskey(raw_cases, case_id) || error("reference case '$case_id' not found in $(reference_manifest_path(root))")

    case = parse_reference_case(case_id, raw_cases[case_id], root)
    actual = compute_reference_actual(case)
    expected = load_reference_array(case)

    println("case: ", case.id)
    summarize_diff(actual, expected)
end

main()
