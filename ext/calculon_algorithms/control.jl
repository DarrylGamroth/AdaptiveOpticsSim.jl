const _DISCRETE_INTEGRATOR_GAIN_PROPERTY = UInt32(0)
const _DISCRETE_INTEGRATOR_TAU_PROPERTY = UInt32(1)

"""Fixed-step Calculon plan around the native AOS controller plan."""
struct _DiscreteIntegratorStepPlan{
    T<:AbstractFloat,
    P<:Control.DiscreteIntegratorPlan{T},
} <: CA.AbstractAlgorithmPlan
    controller_plan::P
    sample_period_s::T
end

function _discrete_integrator_step_plan(
    controller_plan::Control.DiscreteIntegratorPlan{T},
    sample_period_s::T,
) where {T<:AbstractFloat}
    isfinite(sample_period_s) && sample_period_s > zero(T) || throw(
        ArgumentError("sample_period_s must be finite and greater than zero"),
    )
    return _DiscreteIntegratorStepPlan(
        controller_plan,
        sample_period_s,
    )
end

function CA.property_schema(
    ::Type{<:_DiscreteIntegratorStepPlan{T}},
) where {T}
    return (
        CA.runtime_property(
            _DISCRETE_INTEGRATOR_GAIN_PROPERTY,
            :gain,
            T;
            default=T(0.3),
            description="Discrete-integrator input gain",
            unit="1",
            constraint=CA.PropertyRange(-floatmax(T), floatmax(T)),
        ),
        CA.runtime_property(
            _DISCRETE_INTEGRATOR_TAU_PROPERTY,
            :tau_s,
            T;
            default=T(0.02),
            description="First-order DM lag time constant",
            unit="s",
            constraint=CA.PropertyRange(nextfloat(zero(T)), floatmax(T)),
        ),
    )
end

CA.property_values(plan::_DiscreteIntegratorStepPlan) = (
    gain=plan.controller_plan.gain,
    tau_s=plan.controller_plan.tau_s,
)

function CA.replace_properties(
    plan::_DiscreteIntegratorStepPlan{T},
    values::NamedTuple{(:gain,:tau_s),Tuple{T,T}},
) where {T}
    controller_plan = Control.DiscreteIntegratorPlan(values.gain, values.tau_s)
    return _discrete_integrator_step_plan(controller_plan, plan.sample_period_s)
end

CA.@calculon_algorithm DiscreteIntegratorControllerF32 begin
    label = "discrete-integrator-controller-f32"

    @configuration DiscreteIntegratorControllerF32Configuration begin
        extent::Int => graph_rebuild
        gain::Float32 => construction
        tau_s::Float32 => construction
        sample_period_s::Float32 => graph_rebuild
        input_schema::String => graph_rebuild
        output_schema::String => graph_rebuild
    end

    properties = :runtime

    ports(configuration, plan, owners) = (
        input(
            :input,
            Float32,
            (configuration.extent,),
            configuration.input_schema,
        ),
        output(
            :output,
            Float32,
            (configuration.extent,),
            configuration.output_schema;
            metadata_from=:input,
        ),
    )

    prepare(configuration) = begin
        controller = Control.DiscreteIntegratorController(
            configuration.extent;
            gain=configuration.gain,
            tau=configuration.tau_s,
            T=Float32,
        )
        plan = _discrete_integrator_step_plan(
            Control.discrete_integrator_plan(controller),
            configuration.sample_period_s,
        )
        owners = _StateWorkspace(
            Control.discrete_integrator_state(controller),
            Control.discrete_integrator_workspace(controller),
        )
        (plan, owners)
    end

    process(plan, owners) = begin
        Control.update!(
            output,
            owners.state,
            owners.workspace,
            plan.controller_plan,
            input,
            plan.sample_period_s,
        )
    end

    reset(plan, owners) =
        Control.reset_controller!(owners.state, owners.workspace)
end

public DiscreteIntegratorControllerF32
public DiscreteIntegratorControllerF32Configuration
