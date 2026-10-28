module AdaptiveOpticsSimYAMLExt

import AdaptiveOpticsSim.AlgorithmGraphs: AlgorithmGraphError, _read_algorithm_graph
import YAML

function _strict_mapping(constructor, node)
    # YAML.jl skips its duplicate check for merged mappings. Check explicit
    # keys before flattening so aliases can still supply defaults and overrides.
    explicit_keys = Set{String}()
    for (key_node, _) in node.value
        key_node.tag == "tag:yaml.org,2002:merge" && continue
        key = String(YAML.construct_object(constructor, key_node))
        key in explicit_keys && throw(AlgorithmGraphError(
            "duplicate YAML mapping key '$key'",
        ))
        push!(explicit_keys, key)
    end
    return YAML.construct_mapping(
        Dict{String,Any},
        constructor,
        node;
        strict_unique_keys=true,
    )
end

function _read_algorithm_graph(path::AbstractString)
    constructors = Dict("tag:yaml.org,2002:map" => _strict_mapping)
    documents = open(path, "r") do io
        collect(YAML.load_all(io, constructors))
    end
    length(documents) == 1 || throw(AlgorithmGraphError(
        "expected exactly one YAML graph document, got $(length(documents))",
    ))
    return only(documents)
end

include("zemax_preparation.jl")

end
