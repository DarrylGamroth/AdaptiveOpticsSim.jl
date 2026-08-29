struct DiscreteIntegratorNode{T<:AbstractFloat} end

"""Construction and graph-rebuild values for one discrete-integrator node."""
struct DiscreteIntegratorNodeConfig{T<:AbstractFloat}
    extent::Int
    sample_period_s::T
    input_schema::String
    output_schema::String

    function DiscreteIntegratorNodeConfig{T}(
        extent::Integer,
        sample_period_s::Real,
        input_schema::AbstractString,
        output_schema::AbstractString,
    ) where {T<:AbstractFloat}
        extent > 0 || throw(AlgorithmGraphError(
            "discrete-integrator extent must be greater than zero",
        ))
        period = T(sample_period_s)
        isfinite(period) && period > zero(T) || throw(AlgorithmGraphError(
            "discrete-integrator sample_period_s must be finite and greater than zero",
        ))
        isempty(input_schema) && throw(AlgorithmGraphError(
            "discrete-integrator input_schema must not be empty",
        ))
        isempty(output_schema) && throw(AlgorithmGraphError(
            "discrete-integrator output_schema must not be empty",
        ))
        return new{T}(
            Int(extent),
            period,
            String(input_schema),
            String(output_schema),
        )
    end
end

"""Initial scalar props for one discrete-integrator node."""
struct DiscreteIntegratorNodeProps{T<:AbstractFloat}
    gain::T
    tau_s::T

    function DiscreteIntegratorNodeProps{T}(
        gain::Real,
        tau_s::Real,
    ) where {T<:AbstractFloat}
        converted_gain = T(gain)
        converted_tau_s = T(tau_s)
        try
            DiscreteIntegratorPlan(converted_gain, converted_tau_s)
        catch error
            throw(AlgorithmGraphError(sprint(showerror, error)))
        end
        return new{T}(converted_gain, converted_tau_s)
    end
end

"""
    discrete_integrator_node(name; extent, sample_period_s, input_schema,
                             output_schema, gain=0.3, tau_s=0.02,
                             T=Float32)

Declare one fixed-step graph node around the native AOS discrete-integrator
implementation. The node owns its plan, persistent controller state,
replaceable workspace, and exact graph-buffer bindings after preparation.
"""
function discrete_integrator_node(
    name::Symbol;
    extent::Integer,
    sample_period_s::Real,
    input_schema::AbstractString,
    output_schema::AbstractString,
    gain::Real=0.3,
    tau_s::Real=0.02,
    T::Type{<:AbstractFloat}=Float32,
)
    config = DiscreteIntegratorNodeConfig{T}(
        extent,
        sample_period_s,
        input_schema,
        output_schema,
    )
    props = DiscreteIntegratorNodeProps{T}(gain, tau_s)
    return algorithm_node(name, DiscreteIntegratorNode{T}, config; props=props)
end

function graph_node_ports(
    ::Type{DiscreteIntegratorNode{T}},
    config::DiscreteIntegratorNodeConfig{T},
) where {T}
    format(name, direction, schema) = graph_port_contract(
        name,
        direction,
        :data,
        T,
        (config.extent,),
        schema,
        :column_major,
    )
    return (
        format(:input, :input, config.input_schema),
        format(:output, :output, config.output_schema),
    )
end

struct _DiscreteIntegratorStepPlan{
    T<:AbstractFloat,
    P<:DiscreteIntegratorPlan{T},
}
    controller::P
    sample_period_s::T
end

struct _DiscreteIntegratorOwner{Plan,State,Workspace,Input,Output}
    plan::Plan
    state::State
    workspace::Workspace
    input::Input
    output::Output
end

function prepare_graph_node(
    ::Type{DiscreteIntegratorNode{T}},
    config::DiscreteIntegratorNodeConfig{T},
    props::DiscreteIntegratorNodeProps{T},
    inputs::NamedTuple{(:input,)},
    outputs::NamedTuple{(:output,)},
    ::NamedTuple{()},
    target,
) where {T}
    plan = _DiscreteIntegratorStepPlan(
        DiscreteIntegratorPlan(props.gain, props.tau_s),
        config.sample_period_s,
    )
    state = DiscreteIntegratorState(outputs.output)
    workspace = DiscreteIntegratorWorkspace(outputs.output)
    return _DiscreteIntegratorOwner(
        plan,
        state,
        workspace,
        inputs.input,
        outputs.output,
    )
end

@inline function step_graph_node!(owner::_DiscreteIntegratorOwner)
    update!(
        owner.output,
        owner.state,
        owner.workspace,
        owner.plan.controller,
        owner.input,
        owner.plan.sample_period_s,
    )
    return nothing
end

@inline function reset_graph_node!(owner::_DiscreteIntegratorOwner)
    reset_controller!(owner.state, owner.workspace)
    return nothing
end

struct ModalOPDExpansionNode{T<:AbstractFloat} end

"""Construction and graph-rebuild values for one modal OPD expansion node."""
struct ModalOPDExpansionNodeConfig
    pupil_rows::Int
    pupil_columns::Int
    mode_count::Int
    coefficients_schema::String
    opd_schema::String
    basis_schema::String
    pupil_support_schema::String

    function ModalOPDExpansionNodeConfig(
        pupil_rows::Integer,
        pupil_columns::Integer,
        mode_count::Integer,
        coefficients_schema::AbstractString,
        opd_schema::AbstractString,
        basis_schema::AbstractString,
        pupil_support_schema::AbstractString,
    )
        pupil_rows > 0 || throw(AlgorithmGraphError(
            "modal OPD pupil_rows must be greater than zero",
        ))
        pupil_columns > 0 || throw(AlgorithmGraphError(
            "modal OPD pupil_columns must be greater than zero",
        ))
        mode_count > 0 || throw(AlgorithmGraphError(
            "modal OPD mode_count must be greater than zero",
        ))
        schemas = (
            coefficients_schema,
            opd_schema,
            basis_schema,
            pupil_support_schema,
        )
        all(schema -> !isempty(schema), schemas) || throw(AlgorithmGraphError(
            "modal OPD scientific schemas must not be empty",
        ))
        return new(
            Int(pupil_rows),
            Int(pupil_columns),
            Int(mode_count),
            String(coefficients_schema),
            String(opd_schema),
            String(basis_schema),
            String(pupil_support_schema),
        )
    end
end

"""
    modal_opd_expansion_node(name; pupil_rows, pupil_columns, mode_count,
                             coefficients_schema, opd_schema, basis_schema,
                             pupil_support_schema, T=Float32)

Declare one stateless graph node around the native AOS modal OPD expansion.
The graph supplies the modal basis and pupil support as explicit startup sparse
parameters.
"""
function modal_opd_expansion_node(
    name::Symbol;
    pupil_rows::Integer,
    pupil_columns::Integer,
    mode_count::Integer,
    coefficients_schema::AbstractString,
    opd_schema::AbstractString,
    basis_schema::AbstractString,
    pupil_support_schema::AbstractString,
    T::Type{<:AbstractFloat}=Float32,
)
    config = ModalOPDExpansionNodeConfig(
        pupil_rows,
        pupil_columns,
        mode_count,
        coefficients_schema,
        opd_schema,
        basis_schema,
        pupil_support_schema,
    )
    return algorithm_node(
        name,
        ModalOPDExpansionNode{T},
        config;
        props=NamedTuple(),
    )
end

function graph_node_ports(
    ::Type{ModalOPDExpansionNode{T}},
    config::ModalOPDExpansionNodeConfig,
) where {T}
    return (
        graph_port_contract(
            :coefficients,
            :input,
            :data,
            T,
            (config.mode_count,),
            config.coefficients_schema,
            :column_major,
        ),
        graph_port_contract(
            :opd,
            :output,
            :data,
            T,
            (config.pupil_rows, config.pupil_columns),
            config.opd_schema,
            :column_major,
        ),
        graph_port_contract(
            :basis,
            :input,
            :parameter,
            T,
            (config.pupil_rows, config.pupil_columns, config.mode_count),
            config.basis_schema,
            :column_major,
        ),
        graph_port_contract(
            :pupil_support,
            :input,
            :parameter,
            Bool,
            (config.pupil_rows, config.pupil_columns),
            config.pupil_support_schema,
            :column_major,
        ),
    )
end

struct _ModalOPDExpansionOwner{Plan,Coefficients,OPD}
    plan::Plan
    coefficients::Coefficients
    opd::OPD
end

function prepare_graph_node(
    ::Type{ModalOPDExpansionNode{T}},
    ::ModalOPDExpansionNodeConfig,
    ::NamedTuple{()},
    inputs::NamedTuple{(:coefficients,)},
    outputs::NamedTuple{(:opd,)},
    parameters::NamedTuple{(:basis,:pupil_support)},
    target,
) where {T}
    plan = ModalOPDExpansionPlan(parameters.basis, parameters.pupil_support)
    return _ModalOPDExpansionOwner(plan, inputs.coefficients, outputs.opd)
end

@inline function step_graph_node!(owner::_ModalOPDExpansionOwner)
    combine_basis!(owner.opd, owner.plan, owner.coefficients)
    return nothing
end

@inline reset_graph_node!(::_ModalOPDExpansionOwner) = nothing
