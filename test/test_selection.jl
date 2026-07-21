struct TestSuiteSpec{P<:Tuple}
    name::String
    paths::P
end

TestSuiteSpec(name::AbstractString, paths::AbstractString...) =
    TestSuiteSpec(String(name), Tuple(String.(paths)))

# Registry order is the full-suite execution order. Keep it stable so bare
# `Pkg.test()` remains the deterministic composition gate.
const TEST_SUITE_SPECS = (
    TestSuiteSpec("ka-cpu", "ka_cpu_matrix.jl"),
    TestSuiteSpec("tomography", "tomography.jl"),
    TestSuiteSpec("quality", "testsets/quality.jl"),
    TestSuiteSpec("core-optics", "testsets/core_optics.jl"),
    TestSuiteSpec("direct-science", "testsets/direct_science.jl"),
    TestSuiteSpec("atmosphere", "testsets/atmosphere.jl"),
    TestSuiteSpec("plant-topology", "testsets/plant_topology.jl"),
    TestSuiteSpec("plant-time", "testsets/plant_time.jl"),
    TestSuiteSpec("runtime", "testsets/control_and_runtime.jl"),
    TestSuiteSpec(
        "detectors-wfs",
        "testsets/detectors.jl",
        "testsets/wfs_stage_contracts.jl",
        "testsets/shack_hartmann_and_sources.jl",
        "testsets/pyramid_bioedge_and_lgs.jl",
        "testsets/zernike_and_curvature.jl",
    ),
    TestSuiteSpec("plant-preparation", "testsets/plant_preparation.jl"),
    TestSuiteSpec("plant-providers", "testsets/plant_providers.jl"),
    TestSuiteSpec("plant-rng", "testsets/plant_rng.jl"),
    TestSuiteSpec("plant-illumination", "testsets/plant_illumination.jl"),
    TestSuiteSpec(
        "calibration-analysis",
        "testsets/calibration_and_analysis.jl",
    ),
    TestSuiteSpec(
        "reference-tutorials",
        "testsets/reference_and_tutorials.jl",
    ),
    TestSuiteSpec("gate0", "testsets/gate0_characterization.jl"),
    TestSuiteSpec(
        "backend-smoke",
        "backend_optional_common.jl",
        "optional_amdgpu_backends.jl",
        "optional_cuda_backends.jl",
    ),
)

const TEST_GROUP_SPECS = (
    "core" => ("core-optics", "direct-science", "atmosphere"),
    "plant" => (
        "plant-topology",
        "plant-time",
        "plant-preparation",
        "plant-providers",
        "plant-rng",
        "plant-illumination",
    ),
    "control" => ("tomography", "runtime", "calibration-analysis"),
    "sensors" => ("detectors-wfs",),
    "references" => ("reference-tutorials", "gate0"),
    "backends" => ("ka-cpu", "backend-smoke"),
)

test_suite_names() = Tuple(spec.name for spec in TEST_SUITE_SPECS)
test_group_names() = Tuple(first(group) for group in TEST_GROUP_SPECS)

function _test_suite_spec(name::AbstractString)
    for spec in TEST_SUITE_SPECS
        spec.name == name && return spec
    end
    return nothing
end

function _test_group_members(name::AbstractString)
    for (group_name, members) in TEST_GROUP_SPECS
        group_name == name && return members
    end
    return nothing
end

function validate_test_suite_registry()
    suite_names = test_suite_names()
    length(unique(suite_names)) == length(suite_names) ||
        throw(ArgumentError("test suite names must be unique"))

    group_names = test_group_names()
    length(unique(group_names)) == length(group_names) ||
        throw(ArgumentError("test group names must be unique"))
    isempty(intersect(Set(suite_names), Set(group_names))) || throw(
        ArgumentError("test suite and group names must not overlap"),
    )

    registered_paths = String[]
    for spec in TEST_SUITE_SPECS
        isempty(spec.paths) && throw(ArgumentError(
            "test suite '$(spec.name)' must register at least one path"))
        append!(registered_paths, spec.paths)
    end
    length(unique(registered_paths)) == length(registered_paths) ||
        throw(ArgumentError("test paths must belong to exactly one suite"))

    known_suites = Set(suite_names)
    for (group_name, members) in TEST_GROUP_SPECS
        isempty(members) && throw(ArgumentError(
            "test group '$group_name' must contain at least one suite"))
        for member in members
            member in known_suites || throw(ArgumentError(
                "test group '$group_name' names unknown suite '$member'"))
        end
    end
    return nothing
end

validate_test_suite_registry()

function _test_selector_error(selector::AbstractString)
    known = join(("all", test_suite_names()..., test_group_names()...), ", ")
    return ArgumentError(
        "unknown test selector '$selector'; expected one of: $known",
    )
end

"""
    resolve_test_suites(arguments=ARGS)

Resolve suite and group selectors in canonical registry order. Empty arguments
select the complete suite. `--list` returns `nothing` and must be used alone.
Duplicate or overlapping selectors are idempotent.
"""
function resolve_test_suites(arguments=ARGS)
    selectors = String[String(argument) for argument in arguments]
    isempty(selectors) && return TEST_SUITE_SPECS

    if "--list" in selectors
        all(==("--list"), selectors) || throw(ArgumentError(
            "--list cannot be combined with test suite selectors"))
        return nothing
    end
    if "all" in selectors
        all(==("all"), selectors) || throw(ArgumentError(
            "all cannot be combined with other test suite selectors"))
        return TEST_SUITE_SPECS
    end

    selected_names = Set{String}()
    for selector in selectors
        suite = _test_suite_spec(selector)
        if suite !== nothing
            push!(selected_names, suite.name)
            continue
        end

        members = _test_group_members(selector)
        members === nothing && throw(_test_selector_error(selector))
        union!(selected_names, members)
    end
    return Tuple(spec for spec in TEST_SUITE_SPECS
        if spec.name in selected_names)
end

function print_test_suite_help(io::IO=stdout)
    println(io, "Test suites:")
    for spec in TEST_SUITE_SPECS
        println(io, "  ", spec.name)
    end
    println(io, "Test groups:")
    for (name, members) in TEST_GROUP_SPECS
        println(io, "  ", name, " = ", join(members, ", "))
    end
    println(io, "Special selectors: all, --list")
    return nothing
end

function registered_testset_paths()
    paths = String[]
    for spec in TEST_SUITE_SPECS, path in spec.paths
        startswith(path, "testsets/") || continue
        push!(paths, normpath(joinpath(@__DIR__, path)))
    end
    return sort!(paths)
end
