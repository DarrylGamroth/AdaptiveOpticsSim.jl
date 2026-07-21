include("test_selection.jl")

requested_suites = resolve_test_suites()
if requested_suites === nothing
    print_test_suite_help()
else
    include("runtests_head.jl")
    loaded_fixtures = Set{String}()
    for suite in requested_suites
        @info "Running test suite" suite=suite.name
        for fixture in suite.fixtures
            fixture in loaded_fixtures && continue
            @info "Loading test fixture" fixture
            include(fixture)
            push!(loaded_fixtures, fixture)
        end
        for path in suite.paths
            include(path)
        end
    end
end
