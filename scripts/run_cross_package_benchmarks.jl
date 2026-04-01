using AdaptiveOpticsSim
using Dates
using TOML
using LinearAlgebra

const _repo_root = normpath(joinpath(@__DIR__, ".."))

include(joinpath(_repo_root, "test", "reference_harness.jl"))

repo_join(parts...) = normpath(joinpath(_repo_root, parts...))

function relpath_or_identity(path::AbstractString)
    normed = normpath(path)
    if isabspath(normed)
        return relpath(normed, _repo_root)
    end
    return normed
end

function default_cross_package_result_path()
    stamp = Dates.format(Dates.now(), dateformat"yyyy-mm-dd")
    return repo_join("benchmarks", "results", "cross_package", "$(stamp)-phase5-baseline.toml")
end

function resolve_contract_path(path::AbstractString)
    return isabspath(path) ? path : repo_join(path)
end

function parse_profile_value(raw::AbstractString)
    value = strip(raw)
    value == "nothing" && return nothing
    value == "true" && return true
    value == "false" && return false
    int_value = tryparse(Int, value)
    int_value !== nothing && return int_value
    float_value = tryparse(Float64, value)
    float_value !== nothing && return float_value
    return value
end

function parse_profile_output(output::AbstractString)
    lines = split(output, '\n')
    root = nothing
    data = Dict{String,Any}()
    for line in lines
        stripped = rstrip(line)
        isempty(strip(stripped)) && continue
        if root === nothing
            root = strip(stripped)
            continue
        end
        startswith(stripped, "  ") || continue
        payload = strip(stripped)
        parts = split(payload, ": ", limit=2)
        length(parts) == 2 || continue
        key, value = parts
        data[key] = parse_profile_value(value)
    end
    root === nothing && error("profile output did not contain a root block")
    return root, data
end

function normalize_profile_field(key::AbstractString, value)
    if value isa AbstractString && (endswith(key, "_dir") || endswith(key, "_path"))
        return relpath_or_identity(value)
    end
    return value
end

function max_relative_error(actual, expected)
    actual_vec = vec(actual)
    expected_vec = vec(expected)
    finite_mask = .!(isnan.(actual_vec) .| isnan.(expected_vec))
    any(finite_mask) || return 0.0
    a = actual_vec[finite_mask]
    e = expected_vec[finite_mask]
    denom = max.(abs.(e), eps(Float64))
    return maximum(abs.(a .- e) ./ denom)
end

function l2_relative_error(actual, expected)
    actual_vec = vec(actual)
    expected_vec = vec(expected)
    finite_mask = .!(isnan.(actual_vec) .| isnan.(expected_vec))
    any(finite_mask) || return 0.0
    a = actual_vec[finite_mask]
    e = expected_vec[finite_mask]
    denom = norm(e)
    denom == 0 && return norm(a - e)
    return norm(a - e) / denom
end

function run_reference_baseline(baseline::Symbol)
    if baseline === :oopao
        root = default_reference_root()
        has_reference_bundle(root) || return Dict{String,Any}("status" => "skipped", "reason" => "missing_reference_bundle", "root" => relpath_or_identity(root))
    elseif baseline === :specula
        root = default_specula_reference_root()
        has_specula_reference_bundle(root) || return Dict{String,Any}("status" => "skipped", "reason" => "missing_reference_bundle", "root" => relpath_or_identity(root))
    else
        return Dict{String,Any}("status" => "skipped", "reason" => "unsupported_baseline", "baseline" => String(baseline))
    end

    bundle = load_reference_bundle(root)
    cases = reference_cases(bundle, baseline)
    isempty(cases) && return Dict{String,Any}("status" => "skipped", "reason" => "no_cases", "root" => relpath_or_identity(root))

    failed = 0
    max_abs = 0.0
    max_rel = 0.0
    max_l2_rel = 0.0
    for case in cases
        result = validate_reference_case(case)
        failed += result.ok ? 0 : 1
        max_abs = max(max_abs, Float64(result.maxabs))
        if size(result.actual) == size(result.expected)
            max_rel = max(max_rel, max_relative_error(result.actual, result.expected))
            max_l2_rel = max(max_l2_rel, l2_relative_error(result.actual, result.expected))
        else
            max_rel = Inf
            max_l2_rel = Inf
        end
    end

    return Dict{String,Any}(
        "status" => failed == 0 ? "passed" : "failed",
        "root" => relpath_or_identity(root),
        "case_count" => length(cases),
        "failed_case_count" => failed,
        "max_abs_error" => max_abs,
        "max_rel_error" => max_rel,
        "l2_rel_error" => max_l2_rel,
    )
end

function run_reference_scenario(name::AbstractString, cfg::Dict{String,Any})
    baselines = Symbol.(cfg["baselines"])
    baseline_results = Dict{String,Any}()
    failed = 0
    completed = 0
    for baseline in baselines
        result = run_reference_baseline(baseline)
        baseline_results[String(baseline)] = result
        if get(result, "status", "") == "failed"
            failed += 1
        elseif get(result, "status", "") == "passed"
            completed += 1
        end
    end
    status = failed > 0 ? "failed" : (completed > 0 ? "passed" : "skipped")
    return Dict{String,Any}(
        "family" => cfg["family"],
        "class" => cfg["class"],
        "kind" => cfg["kind"],
        "status" => status,
        "description" => cfg["description"],
        "baselines" => baseline_results,
    )
end

function run_profile_implementation(name::AbstractString, cfg::Dict{String,Any})
    cwd = resolve_contract_path(cfg["cwd"])
    cmd = Cmd(Cmd(cfg["command"]); dir=cwd)
    optional = Bool(get(cfg, "optional", false))
    output = ""
    ok = false
    error_message = nothing
    try
        output = read(pipeline(cmd, stderr=stdout), String)
        ok = true
    catch err
        output = sprint(showerror, err)
        error_message = output
    end

    if !ok
        return Dict{String,Any}(
            "label" => cfg["label"],
            "cwd" => relpath_or_identity(cwd),
            "status" => optional ? "skipped" : "failed",
            "optional" => optional,
            "error" => error_message,
        )
    end

    root, data = parse_profile_output(output)
    expected_root = cfg["profile_root"]
    root == expected_root || return Dict{String,Any}(
        "label" => cfg["label"],
        "cwd" => relpath_or_identity(cwd),
        "status" => optional ? "skipped" : "failed",
        "optional" => optional,
        "error" => "unexpected_profile_root",
        "actual_profile_root" => root,
        "expected_profile_root" => expected_root,
    )

    result = Dict{String,Any}(
        "label" => cfg["label"],
        "cwd" => relpath_or_identity(cwd),
        "status" => "passed",
        "profile_root" => root,
        "build_time_ns" => get(data, "build_time_ns", nothing),
        "total_mean_ns" => get(data, cfg["mean_key"], nothing),
        "total_p95_ns" => get(data, cfg["p95_key"], nothing),
        "frame_rate_hz" => get(data, cfg["frame_rate_key"], nothing),
        "total_alloc_bytes" => get(data, cfg["alloc_key"], nothing),
    )

    for field in get(cfg, "comparison_fields", String[])
        result[String(field)] = get(data, String(field), nothing)
    end
    for (alias, source_key) in get(cfg, "field_aliases", Dict{String,Any}())
        result[String(alias)] = get(data, String(source_key), nothing)
    end
    for (key, value) in data
        haskey(result, key) || (result[key] = normalize_profile_field(key, value))
    end
    return result
end

function runtime_ratio(numerator, denominator)
    if numerator === nothing || denominator === nothing
        return nothing
    elseif denominator == 0
        return nothing
    end
    return Float64(numerator) / Float64(denominator)
end

function run_runtime_scenario(name::AbstractString, cfg::Dict{String,Any})
    if get(cfg, "status", "") == "contract_only"
        return Dict{String,Any}(
            "family" => cfg["family"],
            "class" => cfg["class"],
            "kind" => cfg["kind"],
            "status" => "deferred",
            "description" => cfg["description"],
            "skip_reason" => cfg["skip_reason"],
        )
    end

    impl_results = Dict{String,Any}()
    required_failed = false
    required_passed = false
    for (impl_name, impl_cfg_any) in sort!(collect(pairs(cfg)); by=first)
        impl_name in ("family", "class", "kind", "description", "comparison_fields",
            "normalized_fields", "known_differences", "runtime_metrics", "status",
            "skip_reason") && continue
        impl_cfg = Dict{String,Any}(impl_cfg_any)
        if haskey(cfg, "comparison_fields")
            impl_cfg["comparison_fields"] = cfg["comparison_fields"]
        end
        result = run_profile_implementation(String(impl_name), impl_cfg)
        impl_results[String(impl_name)] = result
        optional = Bool(get(impl_cfg, "optional", false))
        if !optional
            required_passed |= get(result, "status", "") == "passed"
            required_failed |= get(result, "status", "") == "failed"
        end
    end

    main = get(impl_results, "main", Dict{String,Any}())
    revolt_real = get(impl_results, "revolt_real", Dict{String,Any}())
    comparison = Dict{String,Any}()
    for field in get(cfg, "comparison_fields", String[])
        field_name = String(field)
        comparison["$(field_name)_match"] = get(main, field_name, nothing) == get(revolt_real, field_name, nothing)
    end
    comparison["mean_time_ratio_main_to_revolt_real"] = runtime_ratio(get(main, "total_mean_ns", nothing), get(revolt_real, "total_mean_ns", nothing))
    comparison["frame_rate_ratio_main_to_revolt_real"] = runtime_ratio(get(main, "frame_rate_hz", nothing), get(revolt_real, "frame_rate_hz", nothing))
    comparison["alloc_ratio_main_to_revolt_real"] = runtime_ratio(get(main, "total_alloc_bytes", nothing), get(revolt_real, "total_alloc_bytes", nothing))

    status = required_failed ? "failed" : (required_passed ? "passed" : "skipped")
    return Dict{String,Any}(
        "family" => cfg["family"],
        "class" => cfg["class"],
        "kind" => cfg["kind"],
        "status" => status,
        "description" => cfg["description"],
        "normalized_fields" => get(cfg, "normalized_fields", String[]),
        "known_differences" => get(cfg, "known_differences", String[]),
        "implementations" => impl_results,
        "comparison" => comparison,
    )
end

function load_contract(path::AbstractString)
    raw = TOML.parsefile(path)
    return Dict{String,Any}(raw)
end

function run_contract(path::AbstractString)
    contract = load_contract(path)
    scenarios_cfg = Dict{String,Any}(contract["scenarios"])
    scenario_results = Dict{String,Any}()
    for name in sort!(collect(keys(scenarios_cfg)))
        cfg = Dict{String,Any}(scenarios_cfg[name])
        kind = String(cfg["kind"])
        if kind == "reference_fidelity"
            scenario_results[name] = run_reference_scenario(name, cfg)
        elseif kind == "runtime_profile"
            scenario_results[name] = run_runtime_scenario(name, cfg)
        else
            scenario_results[name] = Dict{String,Any}("status" => "skipped", "reason" => "unsupported_kind", "kind" => kind)
        end
    end
    return contract, scenario_results
end

function write_cross_package_manifest!(result_path::AbstractString)
    manifest_path = repo_join("benchmarks", "results", "cross_package", "manifest.toml")
    mkpath(dirname(manifest_path))
    rel_result = relpath(result_path, _repo_root)
    entries = String[]
    if isfile(manifest_path)
        raw = TOML.parsefile(manifest_path)
        existing = get(raw, "results", String[])
        append!(entries, String.(existing))
    end
    rel_result in entries || push!(entries, rel_result)
    sort!(entries)
    manifest = Dict{String,Any}(
        "version" => 1,
        "status" => "active",
        "latest_result" => rel_result,
        "results" => entries,
    )
    open(manifest_path, "w") do io
        TOML.print(io, manifest; sorted=true)
    end
    return manifest_path
end

function sanitize_toml_value(value)
    if value === nothing
        return "none"
    elseif value isa Dict
        return Dict{String,Any}(String(k) => sanitize_toml_value(v) for (k, v) in value)
    elseif value isa AbstractVector
        return [sanitize_toml_value(v) for v in value]
    end
    return value
end

function write_results(contract_path::AbstractString, output_path::AbstractString)
    contract, scenario_results = run_contract(contract_path)
    summary = Dict{String,Any}("passed" => 0, "failed" => 0, "skipped" => 0, "deferred" => 0)
    for result in values(scenario_results)
        status = String(get(result, "status", "skipped"))
        if haskey(summary, status)
            summary[status] += 1
        else
            summary[status] = 1
        end
    end

    payload = Dict{String,Any}(
        "version" => 1,
        "status" => "active",
        "generated_at" => string(Dates.now()),
        "contract_path" => relpath(contract_path, _repo_root),
        "summary" => summary,
        "scenarios" => scenario_results,
        "contract_metadata" => Dict{String,Any}(contract["metadata"]),
    )
    payload = sanitize_toml_value(payload)

    mkpath(dirname(output_path))
    open(output_path, "w") do io
        TOML.print(io, payload; sorted=true)
    end
    write_cross_package_manifest!(output_path)
    return output_path
end

function main()
    contract_path = length(ARGS) >= 1 ? resolve_contract_path(ARGS[1]) : repo_join("benchmarks", "contracts", "cross_package.toml")
    output_path = length(ARGS) >= 2 ? resolve_contract_path(ARGS[2]) : default_cross_package_result_path()
    result_path = write_results(contract_path, output_path)
    println("cross_package_benchmark_result")
    println("  contract: ", relpath(contract_path, _repo_root))
    println("  result: ", relpath(result_path, _repo_root))
end

main()
