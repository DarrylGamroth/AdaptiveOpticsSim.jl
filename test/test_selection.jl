struct TestSuiteSpec{F<:Tuple,P<:Tuple}
    name::String
    fixtures::F
    paths::P
end

TestSuiteSpec(name::AbstractString, paths::AbstractString...;
    fixtures=()) = TestSuiteSpec(
        String(name),
        Tuple(String.(fixtures)),
        Tuple(String.(paths)),
    )

const DETECTOR_TEST_SUITE_NAMES = (
    "detector-parameter-ownership",
    "prepared-detector-ownership",
    "detector-shared",
    "detector-lifecycle",
    "detector-response",
    "detector-shared-qualification",
    "detector-ccd",
    "detector-skipper",
    "detector-emccd",
    "detector-cmos",
    "detector-ingaas",
    "detector-thermal",
    "detector-hgcdte",
    "detector-hgcdte-avalanche",
    "detector-linear-apd",
    "detector-spad",
    "detector-mkid",
    "detector-artifacts",
)

# Registry order is the full-suite execution order. Keep it stable so bare
# `Pkg.test()` remains the deterministic composition gate.
const TEST_SUITE_SPECS = (
    TestSuiteSpec("ka-cpu", "ka_cpu_matrix.jl";
        fixtures=("ka_cpu_style_fixture.jl",)),
    TestSuiteSpec("tomography", "tomography.jl"),
    TestSuiteSpec("quality", "testsets/quality.jl"),
    TestSuiteSpec(
        "prepared-execution-interfaces",
        "testsets/prepared_execution_interfaces.jl";
        fixtures=(
            "calibration_interface_conformance.jl",
            "wfs_four_pupil_interface_conformance.jl",
            "backend_execution_context_conformance.jl",
        ),
    ),
    TestSuiteSpec(
        "execution-strategies",
        "testsets/execution_strategies.jl",
    ),
    TestSuiteSpec("algorithm-graphs", "testsets/algorithm_graphs.jl"),
    TestSuiteSpec("prepared-optics-ownership",
        "testsets/prepared_optics_ownership.jl"),
    TestSuiteSpec("core-optics", "testsets/core_optics.jl"),
    TestSuiteSpec("direct-science", "testsets/direct_science.jl"),
    TestSuiteSpec(
        "direct-imaging-batch",
        "testsets/direct_imaging_batch.jl",
        fixtures=("wfs_stage_contract_fixtures.jl",),
    ),
    TestSuiteSpec("atmosphere", "testsets/atmosphere.jl"),
    TestSuiteSpec(
        "atmosphere-direction-batch",
        "testsets/atmosphere_direction_batch.jl",
        fixtures=("wfs_stage_contract_fixtures.jl",),
    ),
    TestSuiteSpec("control-primitives", "testsets/control_primitives.jl"),
    TestSuiteSpec(
        "control-reconstruction",
        "testsets/control_reconstruction.jl",
    ),
    TestSuiteSpec(
        "detector-parameter-ownership",
        "testsets/detector_parameter_ownership.jl";
        fixtures=("detector_test_fixtures.jl",),
    ),
    TestSuiteSpec(
        "prepared-detector-ownership",
        "testsets/prepared_detector_ownership.jl";
        fixtures=("detector_test_fixtures.jl",),
    ),
    TestSuiteSpec(
        "detector-shared",
        "testsets/detector_shared.jl";
        fixtures=("detector_test_fixtures.jl",),
    ),
    TestSuiteSpec(
        "detector-lifecycle",
        "testsets/detector_lifecycle.jl";
        fixtures=("detector_test_fixtures.jl",),
    ),
    TestSuiteSpec(
        "detector-response",
        "testsets/detector_response.jl";
        fixtures=("detector_test_fixtures.jl",),
    ),
    TestSuiteSpec(
        "detector-shared-qualification",
        "testsets/detector_shared_qualification.jl";
        fixtures=("detector_test_fixtures.jl",),
    ),
    TestSuiteSpec(
        "detector-ccd",
        "testsets/detector_ccd.jl";
        fixtures=("detector_test_fixtures.jl",),
    ),
    TestSuiteSpec(
        "detector-skipper",
        "testsets/detector_skipper.jl";
        fixtures=("detector_test_fixtures.jl",),
    ),
    TestSuiteSpec(
        "detector-emccd",
        "testsets/detector_emccd.jl";
        fixtures=("detector_test_fixtures.jl",),
    ),
    TestSuiteSpec(
        "detector-cmos",
        "testsets/detector_cmos.jl";
        fixtures=("detector_test_fixtures.jl",),
    ),
    TestSuiteSpec(
        "detector-ingaas",
        "testsets/detector_ingaas.jl";
        fixtures=("detector_test_fixtures.jl",),
    ),
    TestSuiteSpec(
        "detector-thermal",
        "testsets/detector_thermal.jl";
        fixtures=("detector_test_fixtures.jl",),
    ),
    TestSuiteSpec(
        "detector-hgcdte",
        "testsets/detector_hgcdte.jl";
        fixtures=("detector_test_fixtures.jl",),
    ),
    TestSuiteSpec(
        "detector-hgcdte-avalanche",
        "testsets/detector_hgcdte_avalanche.jl";
        fixtures=("detector_test_fixtures.jl",),
    ),
    TestSuiteSpec(
        "detector-linear-apd",
        "testsets/detector_linear_apd.jl";
        fixtures=("detector_test_fixtures.jl",),
    ),
    TestSuiteSpec(
        "detector-spad",
        "testsets/detector_spad.jl";
        fixtures=("detector_test_fixtures.jl",),
    ),
    TestSuiteSpec(
        "detector-mkid",
        "testsets/detector_mkid.jl";
        fixtures=("detector_test_fixtures.jl",),
    ),
    TestSuiteSpec(
        "detector-artifacts",
        "testsets/detector_artifacts.jl";
        fixtures=("detector_test_fixtures.jl",),
    ),
    TestSuiteSpec(
        "wfs-acquisition-ownership",
        "testsets/wfs_acquisition_ownership.jl",
    ),
    TestSuiteSpec(
        "wfs-common",
        "testsets/wfs_common_and_parity.jl",
        "testsets/wfs_stage_contracts.jl",
        fixtures=("wfs_stage_contract_fixtures.jl",),
    ),
    TestSuiteSpec(
        "wfs-shack-hartmann",
        "testsets/shack_hartmann_and_sources.jl",
        fixtures=("ka_cpu_style_fixture.jl",),
    ),
    TestSuiteSpec(
        "wfs-pyramid-bi-o-edge",
        "testsets/pyramid_bi_o_edge_and_lgs.jl",
        fixtures=("ka_cpu_style_fixture.jl",),
    ),
    TestSuiteSpec(
        "wfs-zernike-curvature",
        "testsets/zernike_and_curvature.jl",
        fixtures=("ka_cpu_style_fixture.jl",),
    ),
    TestSuiteSpec("wfs-lift", "testsets/wfs_lift.jl"),
    TestSuiteSpec(
        "calibration-workflows",
        "testsets/calibration_workflows.jl",
    ),
    TestSuiteSpec("optics-ncpa", "testsets/optics_ncpa.jl"),
    TestSuiteSpec("optical-analysis", "testsets/optical_analysis.jl"),
    TestSuiteSpec(
        "interface-conformance",
        "testsets/interface_conformance.jl",
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
        fixtures=("wfs_stage_contract_fixtures.jl",),
    ),
)

const TEST_GROUP_SPECS = (
    "detectors" => DETECTOR_TEST_SUITE_NAMES,
    "core" => (
        "core-optics",
        "direct-science",
        "direct-imaging-batch",
        "atmosphere",
        "atmosphere-direction-batch",
    ),
    "control" => ("control-primitives", "control-reconstruction"),
    "calibration" => ("calibration-workflows",),
    "sensors" => (
        DETECTOR_TEST_SUITE_NAMES...,
        "wfs-acquisition-ownership",
        "wfs-common",
        "wfs-shack-hartmann",
        "wfs-pyramid-bi-o-edge",
        "wfs-zernike-curvature",
        "wfs-lift",
    ),
    "references" => ("reference-tutorials", "gate0"),
    "backends" => ("ka-cpu", "backend-smoke"),
)

const TEST_CI_SHARD_SPECS = (
    "ci-foundations" => (
        "ka-cpu",
        "tomography",
        "quality",
        "prepared-execution-interfaces",
        "execution-strategies",
        "algorithm-graphs",
        "prepared-optics-ownership",
        "core-optics",
        "direct-science",
        "direct-imaging-batch",
        "atmosphere",
        "atmosphere-direction-batch",
        "optics-ncpa",
        "optical-analysis",
        "reference-tutorials",
        "gate0",
        "backend-smoke",
    ),
    "ci-sensors-control" => (
        "control-primitives",
        "control-reconstruction",
        DETECTOR_TEST_SUITE_NAMES...,
        "wfs-acquisition-ownership",
        "wfs-common",
        "wfs-shack-hartmann",
        "wfs-pyramid-bi-o-edge",
        "wfs-zernike-curvature",
        "wfs-lift",
        "calibration-workflows",
        "interface-conformance",
    ),
)

const TEST_GROUP_SPECS_WITH_CI = (
    TEST_GROUP_SPECS...,
    TEST_CI_SHARD_SPECS...,
)

test_suite_names(specs=TEST_SUITE_SPECS) =
    Tuple(spec.name for spec in specs)
test_group_names(specs=TEST_GROUP_SPECS_WITH_CI) =
    Tuple(first(group) for group in specs)

function _test_suite_spec(name::AbstractString)
    for spec in TEST_SUITE_SPECS
        spec.name == name && return spec
    end
    return nothing
end

function _test_group_members(name::AbstractString)
    for (group_name, members) in TEST_GROUP_SPECS_WITH_CI
        group_name == name && return members
    end
    return nothing
end

function _registered_testset_paths(specs)
    paths = String[]
    for spec in specs, path in spec.paths
        startswith(path, "testsets/") || continue
        push!(paths, normpath(joinpath(@__DIR__, path)))
    end
    return sort!(paths)
end

function _discovered_testset_paths()
    root = joinpath(@__DIR__, "testsets")
    return sort!(normpath.(filter(
        path -> endswith(path, ".jl"),
        readdir(root; join=true),
    )))
end

function validate_test_suite_registry(
    suite_specs,
    group_specs;
    ci_shard_specs=(),
    require_complete::Bool=false,
)
    suite_names = test_suite_names(suite_specs)
    length(unique(suite_names)) == length(suite_names) ||
        throw(ArgumentError("test suite names must be unique"))

    group_names = test_group_names(group_specs)
    length(unique(group_names)) == length(group_names) ||
        throw(ArgumentError("test group names must be unique"))
    isempty(intersect(Set(suite_names), Set(group_names))) || throw(
        ArgumentError("test suite and group names must not overlap"),
    )

    registered_paths = String[]
    registered_fixtures = String[]
    for spec in suite_specs
        isempty(spec.paths) && throw(ArgumentError(
            "test suite '$(spec.name)' must register at least one path"))
        append!(registered_paths, spec.paths)
        append!(registered_fixtures, spec.fixtures)
        for path in spec.paths
            isfile(joinpath(@__DIR__, path)) || throw(ArgumentError(
                "registered test path '$path' does not exist"))
        end
    end
    length(unique(registered_paths)) == length(registered_paths) ||
        throw(ArgumentError("test paths must belong to exactly one suite"))
    isempty(intersect(Set(registered_paths), Set(registered_fixtures))) ||
        throw(ArgumentError(
            "test fixtures must not also be registered as suite paths"))
    for fixture in unique(registered_fixtures)
        isfile(joinpath(@__DIR__, fixture)) || throw(ArgumentError(
            "registered test fixture '$fixture' does not exist"))
    end

    known_suites = Set(suite_names)
    for (group_name, members) in group_specs
        isempty(members) && throw(ArgumentError(
            "test group '$group_name' must contain at least one suite"))
        for member in members
            member in known_suites || throw(ArgumentError(
                "test group '$group_name' names unknown suite '$member'"))
        end
    end

    if !isempty(ci_shard_specs)
        shard_members = String[]
        for (shard_name, members) in ci_shard_specs
            shard_name in group_names || throw(ArgumentError(
                "CI shard '$shard_name' is not a registered test group"))
            append!(shard_members, members)
        end
        length(unique(shard_members)) == length(shard_members) ||
            throw(ArgumentError(
                "CI closure shards must not contain duplicate suites"))
        Set(shard_members) == known_suites || throw(ArgumentError(
            "CI closure shards must include every suite exactly once"))
    end

    if require_complete
        registered_testsets = _registered_testset_paths(suite_specs)
        discovered_testsets = _discovered_testset_paths()
        registered_testsets == discovered_testsets ||
            throw(ArgumentError(
                "testsets directory and suite registry must match exactly"))
    end
    return nothing
end

validate_test_suite_registry() = validate_test_suite_registry(
    TEST_SUITE_SPECS,
    TEST_GROUP_SPECS_WITH_CI;
    ci_shard_specs=TEST_CI_SHARD_SPECS,
    require_complete=true,
)

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
    for (name, members) in TEST_GROUP_SPECS_WITH_CI
        println(io, "  ", name, " = ", join(members, ", "))
    end
    println(io, "Special selectors: all, --list")
    return nothing
end

function registered_testset_paths()
    return _registered_testset_paths(TEST_SUITE_SPECS)
end

function registered_test_fixture_paths()
    paths = String[]
    for spec in TEST_SUITE_SPECS, fixture in spec.fixtures
        push!(paths, normpath(joinpath(@__DIR__, fixture)))
    end
    return sort!(unique!(paths))
end
