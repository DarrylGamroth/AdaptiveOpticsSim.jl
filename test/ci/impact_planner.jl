#!/usr/bin/env julia

include(joinpath(@__DIR__, "..", "test_selection.jl"))

@enum AcceleratorScope begin
    NoAccelerator
    DetectorAccelerator
    FullAccelerator
end

mutable struct ValidationPlan
    selectors::Set{String}
    run_cpu::Bool
    run_cuda::Bool
    run_amdgpu::Bool
    accelerator_scope::AcceleratorScope
    run_aarch64::Bool
    run_grouped_cpu::Bool
    run_scheduler::Bool
    manual_gates::Set{String}
end

ValidationPlan() = ValidationPlan(
    Set{String}(), false, false, false, NoAccelerator,
    false, false, false, Set{String}())

const CLOSURE_SELECTORS = Tuple(first.(TEST_CI_SHARD_SPECS))
const REGISTERED_SELECTORS = Set((test_suite_names()..., test_group_names()...))
const SELECTOR_ORDER = (test_suite_names()..., test_group_names()...)

const DETECTOR_SHARED_SELECTORS = (
    "detector-parameter-ownership",
    "prepared-detector-ownership",
    "detector-shared",
    "detector-lifecycle",
    "detector-response",
    "detector-shared-qualification",
    "detector-thermal",
    "detector-artifacts",
    "plant-detector-transitions",
    "plant-device-model-matrix",
)

const DETECTOR_FILE_SELECTORS = Dict(
    "ccd.jl" => ("detector-ccd", "detector-skipper"),
    "cmos.jl" => ("detector-cmos",),
    "emccd.jl" => ("detector-emccd",),
    "hgcdte.jl" => ("detector-hgcdte",),
    "hgcdte_avalanche_array.jl" => ("detector-hgcdte-avalanche",),
    "ingaas.jl" => ("detector-ingaas",),
    "linear_apd.jl" => ("detector-linear-apd",),
    "mkid_array.jl" => ("detector-mkid",),
    "spad_array.jl" => ("detector-spad",),
)

const WFS_SELECTORS = (
    "wfs-acquisition-ownership",
    "wfs-common",
    "wfs-shack-hartmann",
    "wfs-pyramid-bi-o-edge",
    "wfs-zernike-curvature",
    "wfs-lift",
)

const GPU_RELEVANT_SUITES = Set((
    "ka-cpu",
    "backend-smoke",
    "direct-imaging-batch",
    "atmosphere-direction-batch",
    "plant-device-batching",
    "plant-device-model-matrix",
    WFS_SELECTORS...,
    DETECTOR_TEST_SUITE_NAMES...,
))

const DETECTOR_GPU_SUITES = Set(DETECTOR_TEST_SUITE_NAMES)

function registered_test_path_selectors()
    result = Dict{String,Set{String}}()
    for spec in TEST_SUITE_SPECS, path in (spec.paths..., spec.fixtures...)
        repository_path = replace(joinpath("test", path), '\\' => '/')
        push!(get!(Set{String}, result, repository_path), spec.name)
    end
    return result
end

const REGISTERED_TEST_PATH_SELECTORS = registered_test_path_selectors()

function normalized_repository_path(path::AbstractString)
    normalized = replace(normpath(String(path)), '\\' => '/')
    startswith(normalized, "./") && (normalized = normalized[3:end])
    (isempty(normalized) || isabspath(normalized) || normalized == ".." ||
        startswith(normalized, "../")) && throw(ArgumentError(
            "expected a repository-relative path, got '$path'"))
    return normalized
end

function add_selectors!(plan::ValidationPlan, selectors)
    for selector in selectors
        selector in REGISTERED_SELECTORS || error(
            "unknown test selector '$selector' in validation policy")
        push!(plan.selectors, selector)
    end
    plan.run_cpu = true
    return plan
end

function widen_accelerator_scope!(plan::ValidationPlan, scope::AcceleratorScope)
    Int(scope) > Int(plan.accelerator_scope) && (plan.accelerator_scope = scope)
    return plan
end

function request_cuda!(plan::ValidationPlan, scope::AcceleratorScope)
    plan.run_cuda = true
    return widen_accelerator_scope!(plan, scope)
end

function request_amdgpu!(plan::ValidationPlan, scope::AcceleratorScope)
    plan.run_amdgpu = true
    return widen_accelerator_scope!(plan, scope)
end

function request_accelerators!(plan::ValidationPlan, scope::AcceleratorScope)
    request_cuda!(plan, scope)
    request_amdgpu!(plan, scope)
    return plan
end

function all_validation!(plan::ValidationPlan)
    add_selectors!(plan, CLOSURE_SELECTORS)
    request_accelerators!(plan, FullAccelerator)
    plan.run_aarch64 = true
    plan.run_grouped_cpu = true
    plan.run_scheduler = true
    union!(plan.manual_gates, (
        "macos-selector", "windows-selector", "appleaccelerate",
        "metal", "coverage"))
    return plan
end

function add_registered_test_impacts!(plan::ValidationPlan, path::String)
    selectors = get(REGISTERED_TEST_PATH_SELECTORS, path, nothing)
    selectors === nothing && return false
    add_selectors!(plan, selectors)
    if !isempty(intersect(selectors, GPU_RELEVANT_SUITES))
        scope = issubset(selectors, DETECTOR_GPU_SUITES) ?
            DetectorAccelerator : FullAccelerator
        request_accelerators!(plan, scope)
    end
    if !isempty(intersect(selectors, Set((
            "plant-topology-growth", "plant-event-composition",
            "plant-cpu-execution"))))
        plan.run_grouped_cpu = true
    end
    return true
end

function add_detector_impacts!(plan::ValidationPlan, path::String)
    family = get(DETECTOR_FILE_SELECTORS, basename(path), nothing)
    family === nothing ?
        add_selectors!(plan, ("detectors", "plant-detector-transitions",
            "plant-device-model-matrix")) :
        add_selectors!(plan, (DETECTOR_SHARED_SELECTORS..., family...))
    request_accelerators!(plan, DetectorAccelerator)
    return plan
end

function add_wfs_impacts!(plan::ValidationPlan, path::String)
    selectors = if occursin("shack_hartmann", path)
        ("wfs-acquisition-ownership", "wfs-common", "wfs-shack-hartmann")
    elseif occursin("pyramid", path) || occursin("bi_o_edge", path) ||
            occursin("bioedge", path) || occursin("focal_plane_modulation", path)
        ("wfs-acquisition-ownership", "wfs-common", "wfs-pyramid-bi-o-edge")
    elseif occursin("zernike", path) || occursin("curvature", path)
        ("wfs-acquisition-ownership", "wfs-common", "wfs-zernike-curvature")
    elseif occursin("lift", path)
        ("wfs-acquisition-ownership", "wfs-common", "wfs-lift")
    else
        WFS_SELECTORS
    end
    add_selectors!(plan, selectors)
    request_accelerators!(plan, FullAccelerator)
    return plan
end

function add_plant_impacts!(plan::ValidationPlan, path::String)
    stem = first(splitext(basename(path)))
    candidate = "plant-" * replace(stem, '_' => '-')
    if candidate in REGISTERED_SELECTORS
        add_selectors!(plan, (candidate, "prepared-execution-contracts"))
    else
        add_selectors!(plan, ("ci-plant-runtime", "ci-plant-optics"))
    end
    if any(fragment -> occursin(fragment, stem), (
            "preparation", "resource", "partition", "device", "target",
            "sampled", "direct_measurement", "path_input"))
        request_accelerators!(plan, FullAccelerator)
    end
    if any(fragment -> occursin(fragment, stem), (
            "cpu_execution", "device_batching", "mixed_serial",
            "topology_growth"))
        plan.run_grouped_cpu = true
    end
    return plan
end

function add_source_impacts!(plan::ValidationPlan, path::String)
    startswith(path, "src/") || return false
    if path in ("src/AdaptiveOpticsSim.jl", "src/exports.jl") ||
            startswith(path, "src/core/") || startswith(path, "src/backends/")
        all_validation!(plan)
    elseif startswith(path, "src/detectors/")
        add_detector_impacts!(plan, path)
    elseif startswith(path, "src/wfs/")
        add_wfs_impacts!(plan, path)
    elseif startswith(path, "src/plant/")
        add_plant_impacts!(plan, path)
    elseif startswith(path, "src/optics/") || startswith(path, "src/atmosphere/")
        add_selectors!(plan,
            ("ci-foundations", "ci-sensors-control", "ci-plant-optics"))
        request_accelerators!(plan, FullAccelerator)
    elseif startswith(path, "src/control/")
        add_selectors!(plan, ("control", "plant-controller-routing"))
    elseif startswith(path, "src/calibration/")
        add_selectors!(plan, ("calibration", "control-reconstruction"))
    elseif startswith(path, "src/tomography/")
        add_selectors!(plan, ("tomography", "control-reconstruction"))
    elseif startswith(path, "src/ensembles/")
        add_selectors!(plan, ("ci-foundations",))
        plan.run_grouped_cpu = true
        plan.run_scheduler = true
    elseif startswith(path, "src/algorithm_graphs/")
        add_selectors!(plan, ("algorithm-graphs", "namespace-authority"))
    else
        all_validation!(plan)
    end
    return true
end

function add_extension_impacts!(plan::ValidationPlan, path::String)
    startswith(path, "ext/") || return false
    add_selectors!(plan, ("quality", "interface-conformance"))
    if occursin("CUDA", path)
        request_cuda!(plan, FullAccelerator)
    elseif occursin("AMDGPU", path)
        request_amdgpu!(plan, FullAccelerator)
    elseif occursin("AppleAccelerate", path)
        push!(plan.manual_gates, "appleaccelerate")
    elseif occursin("Metal", path)
        push!(plan.manual_gates, "metal")
    elseif occursin("AcceleratedKernels", path) || occursin("Dagger", path)
        plan.run_grouped_cpu = true
        plan.run_scheduler = true
    end
    return true
end

function add_test_impacts!(plan::ValidationPlan, path::String)
    startswith(path, "test/") || return false
    add_registered_test_impacts!(plan, path) && return true
    if startswith(path, "test/ci/")
        add_selectors!(plan, ("quality",))
    elseif startswith(path, "test/cuda/") || occursin("cuda", lowercase(path))
        request_cuda!(plan, occursin("detector", path) ?
            DetectorAccelerator : FullAccelerator)
    elseif startswith(path, "test/amdgpu/") || occursin("amdgpu", lowercase(path))
        request_amdgpu!(plan, occursin("detector", path) ?
            DetectorAccelerator : FullAccelerator)
    elseif startswith(path, "test/appleaccelerate/")
        push!(plan.manual_gates, "appleaccelerate")
    elseif startswith(path, "test/metal/")
        push!(plan.manual_gates, "metal")
    elseif startswith(path, "test/schedulers/")
        plan.run_scheduler = true
    else
        all_validation!(plan)
    end
    return true
end

is_documentation_path(path::String) =
    startswith(path, "docs/") || path in ("README.md", "LICENSE", "LICENSE.md")

function add_path_impacts!(plan::ValidationPlan, raw_path::AbstractString)
    path = normalized_repository_path(raw_path)
    if is_documentation_path(path)
        add_selectors!(plan, path == "docs/glossary.md" ?
            ("quality", "namespace-authority") : ("quality",))
        return plan
    end
    if path == "Project.toml"
        return all_validation!(plan)
    elseif path in ("AGENTS.md", "codecov.yml") ||
            startswith(path, ".github/workflows/")
        add_selectors!(plan, ("quality",))
        path == "codecov.yml" && push!(plan.manual_gates, "coverage")
        return plan
    elseif startswith(path, "examples/")
        add_selectors!(plan, ("reference-tutorials",))
        return plan
    elseif startswith(path, "scripts/") || startswith(path, "benchmarks/")
        add_selectors!(plan, ("quality",))
        return plan
    end
    add_source_impacts!(plan, path) && return plan
    add_extension_impacts!(plan, path) && return plan
    add_test_impacts!(plan, path) && return plan
    return all_validation!(plan)
end

function plan_validation(paths)
    plan = ValidationPlan()
    isempty(paths) && return all_validation!(plan)
    for path in paths
        add_path_impacts!(plan, path)
    end
    return plan
end

ordered_selectors(plan::ValidationPlan) =
    [selector for selector in SELECTOR_ORDER if selector in plan.selectors]

selector_arguments(plan::ValidationPlan) = join(ordered_selectors(plan), ' ')
selector_input(plan::ValidationPlan) = join(ordered_selectors(plan), ',')

function local_commands(plan::ValidationPlan)
    commands = Pair{String,String}[]
    selectors = selector_arguments(plan)
    plan.run_cpu && push!(commands, "local CPU" =>
        "julia --project=. --startup-file=no test/ci/run_selected_tests.jl $selectors")
    plan.run_grouped_cpu && push!(commands, "local grouped CPU" =>
        "julia --threads=4 --project=. --startup-file=no test/ci/run_selected_tests.jl gate6")
    plan.run_scheduler && push!(commands, "local scheduler extensions" =>
        "julia --threads=4 --project=test/schedulers --startup-file=no test/schedulers/runtests.jl")
    accelerator_target = plan.accelerator_scope == DetectorAccelerator ?
        "runtests_amdgpu_detectors.jl" : "runtests_amdgpu.jl"
    plan.run_amdgpu && push!(commands, "local AMDGPU" =>
        "julia --project=test/amdgpu --startup-file=no test/$accelerator_target")
    cuda_target = plan.accelerator_scope == DetectorAccelerator ?
        "runtests_cuda_detectors.jl" : "runtests_cuda.jl"
    plan.run_cuda && push!(commands, "WSL CUDA" =>
        "ssh wsl 'cd \"\${AOS_WSL_ROOT:-\$HOME/workspaces/codex/AdaptiveOpticsSim.jl}\" && julia --project=test/cuda --startup-file=no test/$cuda_target'")
    plan.run_aarch64 && push!(commands, "Raspberry Pi aarch64" =>
        "ssh raspberrypi 'cd \"\${AOS_RASPBERRYPI_ROOT:-\$HOME/workspaces/codex/AdaptiveOpticsSim.jl}\" && julia --project=. --startup-file=no test/ci/run_selected_tests.jl $selectors'")
    return commands
end

function manual_commands(plan::ValidationPlan)
    commands = Pair{String,String}[]
    selectors = selector_input(plan)
    ref = "--ref \"\$(git branch --show-current)\""
    for gate in ("macos-selector", "windows-selector")
        gate in plan.manual_gates && push!(commands, gate =>
            "gh workflow run \"CPU Validation\" $ref -f gate=$gate -f selectors=$selectors")
    end
    for gate in ("appleaccelerate", "metal")
        gate in plan.manual_gates && push!(commands, gate =>
            "gh workflow run \"CPU Validation\" $ref -f gate=$gate")
    end
    "coverage" in plan.manual_gates && push!(commands, "coverage" =>
        "gh workflow run Coverage $ref")
    return commands
end

function changed_paths(base_sha::AbstractString, head_sha::AbstractString="HEAD")
    isempty(strip(base_sha)) && error("a nonempty base SHA is required")
    command = Cmd([
        "git", "-c", "core.quotePath=false", "diff", "--name-status",
        "-z", "--find-renames", "$(String(base_sha))...$(String(head_sha))",
    ])
    fields = split(read(command, String), '\0'; keepempty=false)
    paths = String[]
    index = 1
    while index <= length(fields)
        status = fields[index]
        index += 1
        if !isempty(status) && first(status) in ('R', 'C')
            index + 1 <= length(fields) || error(
                "malformed git name-status output for status '$status'")
            push!(paths, fields[index], fields[index + 1])
            index += 2
        else
            index <= length(fields) || error(
                "malformed git name-status output for status '$status'")
            push!(paths, fields[index])
            index += 1
        end
    end
    return unique!(paths)
end

function validate_local_markdown_links(repository_root::AbstractString)
    root = abspath(repository_root)
    markdown_files = String[]
    for path in (joinpath(root, "README.md"), joinpath(root, "AGENTS.md"))
        isfile(path) && push!(markdown_files, path)
    end
    docs_root = joinpath(root, "docs")
    if isdir(docs_root)
        for (directory, _, files) in walkdir(docs_root), file in files
            endswith(file, ".md") && push!(markdown_files, joinpath(directory, file))
        end
    end

    inline_link = r"!?\[[^\]]*\]\(([^)]+)\)"
    reference_link = r"^\s*\[[^\]]+\]:\s*(\S+)"m
    failures = String[]
    for markdown in markdown_files
        text = read(markdown, String)
        destinations = String[]
        append!(destinations,
            (match.captures[1] for match in eachmatch(inline_link, text)))
        append!(destinations,
            (match.captures[1] for match in eachmatch(reference_link, text)))
        for raw_destination in destinations
            destination = strip(raw_destination)
            startswith(destination, '<') && endswith(destination, '>') &&
                (destination = destination[2:end-1])
            destination = first(split(destination, r"\s+\""; limit=2))
            isempty(destination) && continue
            startswith(destination, '#') && continue
            occursin(r"^[A-Za-z][A-Za-z0-9+.-]*:", destination) && continue
            startswith(destination, '/') && continue
            local_path = first(split(destination, ('#', '?'); limit=2))
            local_path = replace(local_path, "%20" => " ")
            isempty(local_path) && continue
            resolved = normpath(joinpath(dirname(markdown), local_path))
            ispath(resolved) || push!(failures,
                "$(relpath(markdown, root)) -> $destination")
        end
    end
    isempty(failures) || error(
        "maintained Markdown contains missing local links:\n" * join(failures, '\n'))
    return nothing
end

function validate_impact_policy(repository_root::AbstractString=pwd())
    validate_test_suite_registry()
    for (path, selectors) in REGISTERED_TEST_PATH_SELECTORS
        plan = plan_validation((path,))
        issubset(selectors, plan.selectors) || error(
            "registered test path '$path' does not select its owning suites")
    end
    root = abspath(repository_root)
    tracked = split(read(`git -C $root ls-files -z`, String), '\0'; keepempty=false)
    for path in tracked
        if any(prefix -> startswith(path, prefix), ("src/", "ext/", "test/")) ||
                path == "Project.toml"
            plan = plan_validation((path,))
            (plan.run_cpu || plan.run_cuda || plan.run_amdgpu ||
                plan.run_aarch64 || plan.run_grouped_cpu ||
                plan.run_scheduler || !isempty(plan.manual_gates)) || error(
                    "behavior-bearing path '$path' selects no validation target")
        end
    end
    validate_local_markdown_links(root)
    return nothing
end

function parse_cli(arguments)
    base_sha = nothing
    head_sha = "HEAD"
    paths = String[]
    select_all = false
    validate = false
    index = 1
    while index <= length(arguments)
        argument = arguments[index]
        if argument == "--all"
            select_all = true
        elseif argument == "--validate-repository"
            validate = true
        elseif argument in ("--base-sha", "--head-sha", "--path")
            index += 1
            index <= length(arguments) || error("$argument requires a value")
            value = arguments[index]
            argument == "--base-sha" && (base_sha = value)
            argument == "--head-sha" && (head_sha = value)
            argument == "--path" && push!(paths, value)
        else
            error("unknown argument '$argument'")
        end
        index += 1
    end
    inputs = Int(select_all) + Int(base_sha !== nothing) + Int(!isempty(paths))
    inputs == 1 || error("use exactly one of --all, --base-sha, or --path")
    return (; base_sha, head_sha, paths, select_all, validate)
end

function print_plan(io::IO, paths, plan::ValidationPlan)
    println(io, "Changed paths: ", join(paths, ", "))
    println(io, "Selectors: ", join(ordered_selectors(plan), ", "))
    println(io, "Local validation:")
    for (label, command) in local_commands(plan)
        println(io, "  [$label] $command")
    end
    manual = manual_commands(plan)
    if !isempty(manual)
        println(io, "Optional manual GitHub gates (require explicit approval):")
        for (label, command) in manual
            println(io, "  [$label] $command")
        end
    end
    return nothing
end

function main(arguments=ARGS)
    options = parse_cli(arguments)
    options.validate && validate_impact_policy()
    paths = if options.select_all
        String[]
    elseif options.base_sha !== nothing
        changed_paths(options.base_sha, options.head_sha)
    else
        options.paths
    end
    plan = options.select_all ? all_validation!(ValidationPlan()) :
        plan_validation(paths)
    print_plan(stdout, options.select_all ? ["<all>"] : paths, plan)
    return nothing
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
