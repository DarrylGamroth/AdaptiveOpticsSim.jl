# Benchmark the same compiler and declarations as the instrument applications.
include(joinpath(INSTRUMENTS_ROOT, "tools", "simulation.jl"))

function instrument_system(directory; target, execution, simulation_graph="production",
        settings=NamedTuple())
    configuration = Simulation.Configuration
    declaration = configuration.load_simulation(joinpath(directory, "simulation", "definition.yaml"))
    configured = configuration.settings(declaration, settings)
    document = configuration.graph_document(declaration, configured; simulation_graph)
    return Simulation.prepare_system((; declaration, configured, document); target, execution)
end

function instrument_source_record(root)
    return Dict("root" => root, "revision" => readchomp(`git -C $root rev-parse HEAD`),
        "status" => readchomp(`git -C $root status --porcelain`))
end
