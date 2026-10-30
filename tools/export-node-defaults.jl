# Export source-owned defaults without loading AOS or an accelerator package.
using JSON3, SHA

module NodeDefaults
include(joinpath(@__DIR__, "..", "src", "algorithm_graphs", "node_defaults.jl"))
end

function export_node_defaults(; check=false)
    root = normpath(joinpath(@__DIR__, ".."))
    source = "src/algorithm_graphs/node_defaults.jl"
    document = (;schema_version=1, generated="Generated from AOS YAML-reader defaults.",
        source=(;path=source, sha256=bytes2hex(sha256(read(joinpath(root, source))))),
        regeneration="julia --project=tools tools/export-node-defaults.jl",
        defaults=NodeDefaults._NODE_DEFAULTS)
    output = joinpath(root, "data", "generated", "node-defaults.json")
    contents = JSON3.write(document) * "\n"
    if check
        isfile(output) && read(output, String) == contents || error("node defaults metadata is stale; run $(document.regeneration)")
    else
        mkpath(dirname(output))
        write(output, contents)
    end
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    ARGS in (String[], ["--check"]) || error("usage: export-node-defaults.jl [--check]")
    export_node_defaults(;check=ARGS == ["--check"])
end
