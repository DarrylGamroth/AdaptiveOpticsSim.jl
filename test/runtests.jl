include("test_selection.jl")

requested_suites = resolve_test_suites()
if requested_suites === nothing
    print_test_suite_help()
else
    include("runtests_head.jl")
    for suite in requested_suites
        @info "Running test suite" suite=suite.name
        for path in suite.paths
            include(path)
        end
    end
end
