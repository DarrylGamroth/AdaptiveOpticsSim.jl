struct ProperOpticalTrainNodeFactory{T<:Union{Float32,Float64}} end

"""
    proper_optical_train_node_factory(; T=Float64)

Return the cold factory for a fixed scalar Proper optical train. Loading
`Proper` activates the optional implementation. The built-in declarative graph
map admits this factory as `proper_optical_train_f64` and
`proper_optical_train_f32`; no executable source is read from configuration.
"""
proper_optical_train_node_factory(; T::Type{<:Union{Float32,Float64}}=Float64) =
    ProperOpticalTrainNodeFactory{T}()

function (factory::ProperOpticalTrainNodeFactory{T})(name::Symbol,
        config::NamedTuple, props::NamedTuple) where {T}
    return _proper_optical_train_node(T, name, config, props)
end

function _proper_optical_train_node(::Type, name::Symbol,
        config::NamedTuple, props::NamedTuple)
    throw(AlgorithmGraphError("Proper optical train requires Proper.jl; load it with `using Proper`"))
end

"""
    optical_sampling_m(owner)

Return the physical sample pitch in metres on both axes of a prepared scalar
Proper optical train's raw output intensity. Read this after stepping the
graph; before its first step the wavefront has entrance sampling. Obtain the
owner with `prepared_graph_node(graph, node_name)`.
"""
function optical_sampling_m end
