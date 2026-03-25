using TOML
using Statistics

using AdaptiveOpticsSim

include(joinpath(@__DIR__, "..", "test", "reference_harness.jl"))

function _maxrel(actual, expected)
    diff = abs.(actual .- expected)
    scale = max.(abs.(expected), eps(Float64))
    return maximum(diff ./ scale)
end

function _case_ids(bundle::ReferenceBundle, requested::Vector{String})
    if isempty(requested)
        return [case.id for case in bundle.cases]
    end
    ids = Set(requested)
    found = String[]
    for case in bundle.cases
        if case.id in ids
            push!(found, case.id)
        end
    end
    missing = setdiff(ids, Set(found))
    isempty(missing) || error("reference case(s) not found: $(join(sort!(collect(missing)), ", "))")
    return found
end

function main()
    root = isempty(ARGS) ? default_reference_root() : ARGS[1]
    requested = length(ARGS) <= 1 ? String[] : String.(ARGS[2:end])

    has_reference_bundle(root) || error("reference bundle not found at $(reference_manifest_path(root))")
    bundle = load_reference_bundle(root)
    selected = Set(_case_ids(bundle, requested))

    println("OOPAO reference audit")
    println("  root: ", root)
    println("  selected_cases: ", length(selected))

    n_ok = 0
    worst_case = nothing
    worst_score = -Inf

    for case in bundle.cases
        case.id in selected || continue
        result = validate_reference_case(case)
        maxrel = if size(result.actual) == size(result.expected)
            _maxrel(result.actual, result.expected)
        else
            Inf
        end
        score = result.maxabs / max(case.atol, eps(Float64))
        if score > worst_score
            worst_score = score
            worst_case = (id=case.id, kind=case.kind, maxabs=result.maxabs, maxrel=maxrel, atol=case.atol, rtol=case.rtol, ok=result.ok)
        end
        status = result.ok ? "ok" : "FAIL"
        println("  ", status, " :: ", case.id,
            " kind=", case.kind,
            " maxabs=", result.maxabs,
            " maxrel=", maxrel,
            " atol=", case.atol,
            " rtol=", case.rtol)
        n_ok += result.ok ? 1 : 0
    end

    println("  passed: ", n_ok, "/", length(selected))
    if worst_case !== nothing
        println("  worst_case: ", worst_case.id,
            " kind=", worst_case.kind,
            " maxabs=", worst_case.maxabs,
            " maxrel=", worst_case.maxrel,
            " atol=", worst_case.atol,
            " rtol=", worst_case.rtol,
            " ok=", worst_case.ok)
    end

    n_ok == length(selected) || error("OOPAO reference audit failed")
end

main()
