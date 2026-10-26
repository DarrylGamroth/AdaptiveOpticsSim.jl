import ..Optics

"""Declaration of an independent-coordinate physical response graph node."""
struct SecondOrderResponseNode{T<:AbstractFloat} end

"""Cold configuration of a same-quantity independent-coordinate response."""
struct SecondOrderResponseNodeConfig
    natural_frequency_hz::Float64
    damping_ratio::Float64
    sample_period_s::Float64
    count::Int
    schema::String
end

_valid_response_count(count::Integer) = 0 < count <= typemax(Int)
_valid_response_count(::Bool) = false

"""
    second_order_response_node(name; natural_frequency_hz, damping_ratio,
        sample_period_s, count, schema, T=Float32)

Declare a generic independent-coordinate second-order physical response with
same-schema vector ports `:input` and `:output`. Each step holds the current
input for `sample_period_s` and emits the end-of-interval position. The caller
must match this period to its graph cadence. Preparation/reset starts at rest;
frequency and damping are explicit model parameters, with no instrument defaults.
"""
function second_order_response_node(name::Symbol; natural_frequency_hz::Real,
        damping_ratio::Real, sample_period_s::Real, count::Integer,
        schema::AbstractString, T::Type=Float32)
    _valid_response_count(count) || throw(AlgorithmGraphError(
        "second-order response count must be a positive machine integer"))
    isempty(schema) && throw(AlgorithmGraphError("second-order response schema must not be empty"))
    plan = Optics.SecondOrderResponsePlan(;
        natural_frequency_hz, damping_ratio, sample_period_s, T)
    config = SecondOrderResponseNodeConfig(plan.natural_frequency_hz,
        plan.damping_ratio, plan.sample_period_s, Int(count), String(schema))
    return algorithm_node(name, SecondOrderResponseNode{T}, config; props=NamedTuple())
end

function graph_node_ports(::Type{SecondOrderResponseNode{T}}, config::SecondOrderResponseNodeConfig) where {T}
    return (graph_port_contract(:input, :input, :data, T, (config.count,), config.schema, :column_major),
        graph_port_contract(:output, :output, :data, T, (config.count,), config.schema, :column_major))
end

struct _SecondOrderResponseOwner{P}
    prepared::P
end

function prepare_graph_node(::Type{SecondOrderResponseNode{T}}, config::SecondOrderResponseNodeConfig,
        ::NamedTuple{()}, inputs::NamedTuple{(:input,)}, outputs::NamedTuple{(:output,)},
        ::NamedTuple{()}, target) where {T}
    plan = Optics.SecondOrderResponsePlan(; natural_frequency_hz=config.natural_frequency_hz,
        damping_ratio=config.damping_ratio, sample_period_s=config.sample_period_s, T)
    return _SecondOrderResponseOwner(Optics.prepare_second_order_response(inputs.input,
        plan; output=outputs.output))
end

@inline function step_graph_node!(owner::_SecondOrderResponseOwner)
    Optics.advance_response!(owner.prepared)
    return nothing
end

@inline function enqueue_graph_node!(owner::_SecondOrderResponseOwner)
    Optics._enqueue_response!(execution_style(owner.prepared.input), owner.prepared)
    return nothing
end

@inline _second_order_response_capture_capability(::AbstractComputeDevice) = GraphNodeCaptureUnsupported()
@inline _second_order_response_capture_capability(::AcceleratorComputeDevice) = GraphNodeCaptureSafe()
@inline graph_node_capture_capability(owner::_SecondOrderResponseOwner) =
    _second_order_response_capture_capability(compute_device(owner.prepared.input))

function reset_graph_node!(owner::_SecondOrderResponseOwner)
    Optics.reset_response!(owner.prepared)
    return nothing
end
