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

struct ClosedLoopCorrectionNode{T<:AbstractFloat} end

"""Construction values for one same-coordinate closed-loop correction node."""
struct ClosedLoopCorrectionNodeConfig
    extent::Int
    residual_error_schema::String
    constraint_feedback_schema::String
    correction_schema::String
    controller_state_schema::String

    function ClosedLoopCorrectionNodeConfig(
        extent::Integer,
        residual_error_schema::AbstractString,
        constraint_feedback_schema::AbstractString,
        correction_schema::AbstractString,
        controller_state_schema::AbstractString,
    )
        extent > 0 || throw(AlgorithmGraphError(
            "closed-loop correction extent must be positive",
        ))
        schemas = (
            residual_error_schema,
            constraint_feedback_schema,
            correction_schema,
            controller_state_schema,
        )
        all(schema -> !isempty(schema), schemas) || throw(AlgorithmGraphError(
            "closed-loop correction schemas must not be empty",
        ))
        return new(
            Int(extent),
            String(residual_error_schema),
            String(constraint_feedback_schema),
            String(correction_schema),
            String(controller_state_schema),
        )
    end
end

"""Initial scalar props for one closed-loop correction node."""
struct ClosedLoopCorrectionNodeProps{T<:AbstractFloat}
    gain::T
    pole::T
    anti_windup_gain::T

    function ClosedLoopCorrectionNodeProps{T}(
        gain::Real,
        pole::Real,
        anti_windup_gain::Real,
    ) where {T<:AbstractFloat}
        converted_gain = T(gain)
        converted_pole = T(pole)
        converted_anti_windup_gain = T(anti_windup_gain)
        try
            ClosedLoopCorrectionPlan(
                1,
                converted_gain,
                converted_pole,
                converted_anti_windup_gain,
            )
        catch error
            error isa AdaptiveOpticsSimError || rethrow()
            throw(AlgorithmGraphError(sprint(showerror, error)))
        end
        return new{T}(
            converted_gain,
            converted_pole,
            converted_anti_windup_gain,
        )
    end
end

"""
    closed_loop_correction_node(name; ...)

Declare one atomic, complete-frame leaky correction operation with delayed
same-coordinate constraint feedback. The first step starts from zero state and
ignores its placeholder feedback. Later steps subtract the supplied feedback
(`preceding demanded correction - realized correction`) from the preceding
correction before applying the pole and current residual error. Hidden-mode
removal, coordinate transforms, command constraints, and DM response remain
separate operations.
"""
function closed_loop_correction_node(
    name::Symbol;
    extent::Integer,
    residual_error_schema::AbstractString,
    constraint_feedback_schema::AbstractString,
    correction_schema::AbstractString,
    controller_state_schema::AbstractString,
    gain::Real=1.0,
    pole::Real=0.0,
    anti_windup_gain::Real=0.0,
    T::Type{<:AbstractFloat}=Float32,
)
    config = ClosedLoopCorrectionNodeConfig(
        extent,
        residual_error_schema,
        constraint_feedback_schema,
        correction_schema,
        controller_state_schema,
    )
    props = ClosedLoopCorrectionNodeProps{T}(
        gain,
        pole,
        anti_windup_gain,
    )
    return algorithm_node(name, ClosedLoopCorrectionNode{T}, config; props)
end

function graph_node_ports(
    ::Type{ClosedLoopCorrectionNode{T}},
    config::ClosedLoopCorrectionNodeConfig,
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
        format(:residual_error, :input, config.residual_error_schema),
        format(
            :constraint_feedback,
            :input,
            config.constraint_feedback_schema,
        ),
        format(:correction, :output, config.correction_schema),
        format(:controller_state, :output, config.controller_state_schema),
    )
end

struct _ClosedLoopCorrectionOwner{
    Plan,
    State,
    Workspace,
    ResidualError,
    ConstraintFeedback,
    Correction,
    ControllerState,
}
    plan::Plan
    state::State
    workspace::Workspace
    residual_error::ResidualError
    constraint_feedback::ConstraintFeedback
    correction::Correction
    controller_state::ControllerState
end

function prepare_graph_node(
    ::Type{ClosedLoopCorrectionNode{T}},
    config::ClosedLoopCorrectionNodeConfig,
    props::ClosedLoopCorrectionNodeProps{T},
    inputs::NamedTuple{(:residual_error,:constraint_feedback)},
    outputs::NamedTuple{(:correction,:controller_state)},
    ::NamedTuple{()},
    target,
) where {T}
    plan = ClosedLoopCorrectionPlan(
        config.extent,
        props.gain,
        props.pole,
        props.anti_windup_gain,
    )
    state = ClosedLoopCorrectionState(plan, outputs.correction)
    workspace = ClosedLoopCorrectionWorkspace(plan, outputs.correction)
    return _ClosedLoopCorrectionOwner(
        plan,
        state,
        workspace,
        inputs.residual_error,
        inputs.constraint_feedback,
        outputs.correction,
        outputs.controller_state,
    )
end

@inline function step_graph_node!(owner::_ClosedLoopCorrectionOwner)
    apply_closed_loop_correction!(
        owner.correction,
        owner.state,
        owner.workspace,
        owner.plan,
        owner.residual_error,
        owner.constraint_feedback,
    )
    copyto!(owner.controller_state, owner.state.integrator_state)
    return nothing
end

@inline function reset_graph_node!(owner::_ClosedLoopCorrectionOwner)
    reset_closed_loop_correction!(owner.state, owner.workspace)
    fill!(owner.correction, zero(eltype(owner.correction)))
    fill!(owner.controller_state, zero(eltype(owner.controller_state)))
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

struct MultiLayerAtmosphereOPDNode{T<:AbstractFloat} end

"""Construction values for one finite multilayer atmosphere OPD node."""
struct MultiLayerAtmosphereOPDNodeConfig{TD,AD,T<:AbstractFloat}
    telescope::TD
    atmosphere::AD
    atmosphere_step::T
    rng_seed::UInt64
    atmosphere_opd_schema::String
end

function _multilayer_atmosphere_opd_config(
    ::Type{T};
    resolution::Integer,
    telescope_diameter_m::Real,
    central_obstruction_ratio::Real,
    pupil_reflectivity::Real,
    aperture_revision::Integer,
    r0::Real,
    reference_wavelength_m::Real,
    L0::Real,
    fractional_cn2,
    wind_speed,
    wind_direction_deg,
    altitude,
    layer_ids,
    atmosphere_step::Real,
    rng_seed::Integer,
    atmosphere_opd_schema::AbstractString,
) where {T<:AbstractFloat}
    atmosphere_step isa Bool && throw(AlgorithmGraphError(
        "atmosphere_step must be a real duration, not Bool",
    ))
    step = T(atmosphere_step)
    isfinite(step) && step > zero(T) || throw(AlgorithmGraphError(
        "atmosphere_step must be finite and positive",
    ))
    rng_seed isa Bool && throw(AlgorithmGraphError(
        "multilayer atmosphere rng_seed must be an integer, not Bool",
    ))
    rng_seed >= 0 || throw(AlgorithmGraphError(
        "multilayer atmosphere rng_seed must be nonnegative",
    ))
    rng_seed <= typemax(UInt64) || throw(AlgorithmGraphError(
        "multilayer atmosphere rng_seed exceeds UInt64",
    ))
    isempty(atmosphere_opd_schema) && throw(AlgorithmGraphError(
        "multilayer atmosphere OPD schema must not be empty",
    ))

    telescope, atmosphere = try
        telescope = TelescopeDefinition(
            resolution=resolution,
            diameter=telescope_diameter_m,
            central_obstruction=central_obstruction_ratio,
            pupil_reflectivity=pupil_reflectivity,
            revision=aperture_revision,
            T=T,
        )
        atmosphere = MultiLayerAtmosphereDefinition(
            r0=r0,
            reference_wavelength_m=reference_wavelength_m,
            L0=L0,
            fractional_cn2=fractional_cn2,
            wind_speed=wind_speed,
            wind_direction_deg=wind_direction_deg,
            altitude=altitude,
            layer_ids=layer_ids,
            T=T,
        )
        (telescope, atmosphere)
    catch error
        error isa AdaptiveOpticsSimError || rethrow()
        throw(AlgorithmGraphError(sprint(showerror, error)))
    end
    return MultiLayerAtmosphereOPDNodeConfig(
        telescope,
        atmosphere,
        step,
        UInt64(rng_seed),
        String(atmosphere_opd_schema),
    )
end

"""
    multilayer_atmosphere_opd_node(name; resolution,
        telescope_diameter_m, r0, reference_wavelength_m, fractional_cn2, wind_speed,
        wind_direction_deg, altitude, layer_ids, atmosphere_step,
        rng_seed, atmosphere_opd_schema, ...)

Declare one finite moving-screen multilayer-atmosphere operation. Preparation
creates the stochastic phase screens and one on-axis infinite-source direction
renderer on the graph's exact compute device. Each invocation advances the
single atmosphere writer by `atmosphere_step` seconds and publishes one
complete pupil-plane atmospheric OPD map in metres.

The node owns its evolving atmosphere and RNG state. It does not infer cadence
from a detector exposure or wall time, apply a deformable-mirror surface, or
retain historical atmosphere epochs.
"""
function multilayer_atmosphere_opd_node(
    name::Symbol;
    resolution::Integer,
    telescope_diameter_m::Real,
    r0::Real,
    reference_wavelength_m::Real,
    fractional_cn2,
    wind_speed,
    wind_direction_deg,
    altitude,
    layer_ids,
    atmosphere_step::Real,
    rng_seed::Integer,
    atmosphere_opd_schema::AbstractString,
    L0::Real=25,
    central_obstruction_ratio::Real=0,
    pupil_reflectivity::Real=1,
    aperture_revision::Integer=0,
    T::Type{<:AbstractFloat}=Float32,
)
    config = _multilayer_atmosphere_opd_config(
        T;
        resolution,
        telescope_diameter_m,
        central_obstruction_ratio,
        pupil_reflectivity,
        aperture_revision,
        r0,
        reference_wavelength_m,
        L0,
        fractional_cn2,
        wind_speed,
        wind_direction_deg,
        altitude,
        layer_ids,
        atmosphere_step,
        rng_seed,
        atmosphere_opd_schema,
    )
    return algorithm_node(
        name,
        MultiLayerAtmosphereOPDNode{T},
        config;
        props=NamedTuple(),
    )
end

function graph_node_ports(
    ::Type{MultiLayerAtmosphereOPDNode{T}},
    config::MultiLayerAtmosphereOPDNodeConfig,
) where {T}
    resolution = config.telescope.resolution
    return (
        graph_port_contract(
            :atmosphere_opd,
            :output,
            :data,
            T,
            (resolution, resolution),
            config.atmosphere_opd_schema,
            :column_major,
        ),
    )
end

mutable struct _MultiLayerAtmosphereOPDOwner{
    D,
    TEL,
    TARGET,
    A,
    R,
    P,
    RNG,
    REPLAY,
    T<:AbstractFloat,
}
    definition::D
    telescope::TEL
    target::TARGET
    atmosphere::A
    renderer::R
    pupil::P
    rng::RNG
    replay::REPLAY
    atmosphere_step::T
    rng_seed::UInt64
end

function prepare_graph_node(
    ::Type{MultiLayerAtmosphereOPDNode{T}},
    config::MultiLayerAtmosphereOPDNodeConfig,
    ::NamedTuple{()},
    ::NamedTuple{()},
    outputs::NamedTuple{(:atmosphere_opd,)},
    ::NamedTuple{()},
    target,
) where {T}
    telescope = prepare_telescope(config.telescope, target)
    atmosphere = prepare_timed_atmosphere(
        config.atmosphere,
        telescope,
        target,
    )
    renderer = prepare_atmosphere_renderer(atmosphere, telescope)
    pupil = PupilFunction(telescope, outputs.atmosphere_opd)
    rng = _prepare_graph_rng(target, config.rng_seed)
    advance_by!(atmosphere, zero(T), rng)
    replay = _prepare_moving_screen_replay(
        atmosphere,
        renderer,
        config.atmosphere_step,
    )
    return _MultiLayerAtmosphereOPDOwner(
        config.atmosphere,
        telescope,
        target,
        atmosphere,
        renderer,
        pupil,
        rng,
        replay,
        config.atmosphere_step,
        config.rng_seed,
    )
end

@inline function step_graph_node!(owner::_MultiLayerAtmosphereOPDOwner)
    epoch = advance_by!(
        owner.atmosphere,
        owner.atmosphere_step,
        owner.rng,
    )
    render_atmosphere!(
        owner.pupil,
        owner.renderer,
        owner.atmosphere,
        epoch,
    )
    return nothing
end

@inline function enqueue_graph_node!(owner::_MultiLayerAtmosphereOPDOwner)
    epoch = advance_by!(
        owner.atmosphere,
        owner.atmosphere_step,
        owner.rng,
    )
    _render_atmosphere_async!(
        owner.pupil,
        owner.renderer,
        owner.atmosphere,
        epoch,
    )
    return nothing
end

@inline function enqueue_captured_graph_node!(
    owner::_MultiLayerAtmosphereOPDOwner,
)
    _enqueue_moving_screen_replay!(
        owner.pupil.opd,
        owner.renderer,
        owner.atmosphere,
        owner.replay,
    )
    return nothing
end

@inline function _preflight_captured_graph_node!(
    owner::_MultiLayerAtmosphereOPDOwner,
)
    _preflight_atmosphere_replay_step(
        owner.atmosphere,
        owner.atmosphere_step,
    )
    return nothing
end

@inline function _complete_captured_graph_node!(
    owner::_MultiLayerAtmosphereOPDOwner,
)
    advance_by!(
        owner.atmosphere,
        owner.atmosphere_step,
        owner.rng,
    )
    return nothing
end

@inline graph_node_capture_capability(
    owner::_MultiLayerAtmosphereOPDOwner,
) = _multilayer_atmosphere_capture_capability(owner.target)

@inline _multilayer_atmosphere_capture_capability(
    ::AcceleratorComputeDevice,
) = GraphNodeCaptureSafe()

@inline _multilayer_atmosphere_capture_capability(
    ::AbstractComputeDevice,
) = GraphNodeCaptureUnsupported()

function reset_graph_node!(owner::_MultiLayerAtmosphereOPDOwner)
    _reset_graph_rng!(owner.rng, owner.rng_seed)
    _reset_multilayer_atmosphere!(owner.atmosphere, owner.rng)
    _reset_moving_screen_replay!(owner.replay)
    fill!(owner.pupil.opd, zero(eltype(owner.pupil.opd)))
    return nothing
end

struct DeformableMirrorSurfaceNode{T<:AbstractFloat} end

"""Construction values for one measured deformable-mirror surface node."""
struct DeformableMirrorSurfaceNodeConfig{TD,AM}
    telescope::TD
    actuator_count::Int
    actuator_model::AM
    pdm_command_schema::String
    surface_opd_schema::String
    actuator_coordinates_schema::String
    influence_functions_schema::String
end

function _deformable_mirror_surface_config(
    ::Type{T};
    resolution::Integer,
    telescope_diameter_m::Real,
    central_obstruction_ratio::Real,
    pupil_reflectivity::Real,
    aperture_revision::Integer,
    actuator_count::Integer,
    actuator_model,
    pdm_command_schema::AbstractString,
    surface_opd_schema::AbstractString,
    actuator_coordinates_schema::AbstractString,
    influence_functions_schema::AbstractString,
) where {T<:AbstractFloat}
    count = Int(actuator_count)
    count > 0 || throw(AlgorithmGraphError(
        "deformable-mirror actuator_count must be greater than zero",
    ))
    schemas = (
        pdm_command_schema,
        surface_opd_schema,
        actuator_coordinates_schema,
        influence_functions_schema,
    )
    all(schema -> !isempty(schema), schemas) || throw(AlgorithmGraphError(
        "deformable-mirror scientific schemas must not be empty",
    ))
    telescope = try
        TelescopeDefinition(
            resolution=Int(resolution),
            diameter=telescope_diameter_m,
            central_obstruction=central_obstruction_ratio,
            pupil_reflectivity=pupil_reflectivity,
            revision=Int(aperture_revision),
            T=T,
        )
    catch error
        error isa AdaptiveOpticsSimError || rethrow()
        throw(AlgorithmGraphError(sprint(showerror, error)))
    end
    return DeformableMirrorSurfaceNodeConfig(
        telescope,
        count,
        actuator_model,
        String(pdm_command_schema),
        String(surface_opd_schema),
        String(actuator_coordinates_schema),
        String(influence_functions_schema),
    )
end

"""
    deformable_mirror_surface_node(name; resolution, telescope_diameter_m,
                                   actuator_count, pdm_command_schema,
                                   surface_opd_schema,
                                   actuator_coordinates_schema,
                                   influence_functions_schema, ...)

Declare one complete-command physical deformable-mirror operation. The graph
supplies sampled actuator coordinates and calibrated influence functions as
startup sparse parameters. Each step applies the complete PDM command through
the native [`DeformableMirror`](@ref) implementation and publishes one complete
surface OPD. No actuator layout or influence width is inferred.

The default actuator response is the native linear static response. Direct
Julia construction may supply an explicit native `actuator_model`; the
maintained TOML factory deliberately keeps that choice out of configuration
until a typed file representation is defined.
"""
function deformable_mirror_surface_node(
    name::Symbol;
    resolution::Integer,
    telescope_diameter_m::Real,
    actuator_count::Integer,
    pdm_command_schema::AbstractString,
    surface_opd_schema::AbstractString,
    actuator_coordinates_schema::AbstractString,
    influence_functions_schema::AbstractString,
    central_obstruction_ratio::Real=0,
    pupil_reflectivity::Real=1,
    aperture_revision::Integer=0,
    actuator_model=nothing,
    T::Type{<:AbstractFloat}=Float32,
)
    config = _deformable_mirror_surface_config(
        T;
        resolution,
        telescope_diameter_m,
        central_obstruction_ratio,
        pupil_reflectivity,
        aperture_revision,
        actuator_count,
        actuator_model,
        pdm_command_schema,
        surface_opd_schema,
        actuator_coordinates_schema,
        influence_functions_schema,
    )
    return algorithm_node(
        name,
        DeformableMirrorSurfaceNode{T},
        config;
        props=NamedTuple(),
    )
end

function graph_node_ports(
    ::Type{DeformableMirrorSurfaceNode{T}},
    config::DeformableMirrorSurfaceNodeConfig,
) where {T}
    resolution = config.telescope.resolution
    return (
        graph_port_contract(
            :pdm_command,
            :input,
            :data,
            T,
            (config.actuator_count,),
            config.pdm_command_schema,
            :column_major,
        ),
        graph_port_contract(
            :surface_opd,
            :output,
            :data,
            T,
            (resolution, resolution),
            config.surface_opd_schema,
            :column_major,
        ),
        graph_port_contract(
            :actuator_coordinates,
            :input,
            :parameter,
            T,
            (2, config.actuator_count),
            config.actuator_coordinates_schema,
            :column_major,
        ),
        graph_port_contract(
            :influence_functions,
            :input,
            :parameter,
            T,
            (resolution * resolution, config.actuator_count),
            config.influence_functions_schema,
            :column_major,
        ),
    )
end

struct _DeformableMirrorSurfaceOwner{DM,Command,Output}
    deformable_mirror::DM
    pdm_command::Command
    output::Output
end

@kernel function _count_nonfinite_graph_values_kernel!(count, values, n::Int)
    index = @index(Global, Linear)
    if index <= n && !isfinite(@inbounds(values[index]))
        @inbounds KernelAbstractions.@atomic count[1] = UInt32(1)
    end
end

function _require_finite_graph_values(
    ::ScalarCPUStyle,
    values::AbstractArray,
    role::AbstractString,
)
    all(isfinite, values) || throw(AlgorithmGraphError("$role must be finite"))
    return nothing
end

function _require_finite_graph_values(
    style::AcceleratorStyle,
    values::AbstractArray,
    role::AbstractString,
)
    nonfinite_count = similar(values, UInt32, 1)
    fill!(nonfinite_count, UInt32(0))
    launch_kernel!(
        style,
        _count_nonfinite_graph_values_kernel!,
        nonfinite_count,
        values,
        length(values);
        ndrange=length(values),
    )
    host_count = Vector{UInt32}(undef, 1)
    copyto!(host_count, nonfinite_count)
    iszero(only(host_count)) || throw(AlgorithmGraphError(
        "$role must be finite",
    ))
    return nothing
end

@inline function _require_finite_graph_values(
    values::AbstractArray,
    role::AbstractString,
)
    return _require_finite_graph_values(execution_style(values), values, role)
end

function prepare_graph_node(
    ::Type{DeformableMirrorSurfaceNode{T}},
    config::DeformableMirrorSurfaceNodeConfig,
    ::NamedTuple{()},
    inputs::NamedTuple{(:pdm_command,)},
    outputs::NamedTuple{(:surface_opd,)},
    parameters::NamedTuple{(:actuator_coordinates,:influence_functions)},
    target,
) where {T}
    actuator_coordinates = Array(parameters.actuator_coordinates)
    all(isfinite, actuator_coordinates) || throw(AlgorithmGraphError(
        "deformable-mirror actuator coordinates must be finite",
    ))
    _require_finite_graph_values(
        parameters.influence_functions,
        "deformable-mirror influence functions",
    )
    topology = SampledActuatorTopology(actuator_coordinates; T=T)
    metadata = (schema=config.influence_functions_schema,)
    influence_model = MeasuredInfluenceFunctions{
        T,
        typeof(parameters.influence_functions),
        typeof(metadata),
    }(
        parameters.influence_functions,
        metadata,
    )
    telescope = prepare_telescope(config.telescope, target)
    deformable_mirror = DeformableMirror(
        telescope;
        topology,
        influence_model,
        actuator_model=config.actuator_model,
        T=T,
        backend=compute_device_backend(target),
    )
    return _DeformableMirrorSurfaceOwner(
        deformable_mirror,
        inputs.pdm_command,
        outputs.surface_opd,
    )
end

@inline function step_graph_node!(owner::_DeformableMirrorSurfaceOwner)
    set_command!(owner.deformable_mirror, owner.pdm_command)
    update_surface!(owner.deformable_mirror)
    copyto!(owner.output, surface_opd(owner.deformable_mirror))
    return nothing
end

@inline function enqueue_graph_node!(owner::_DeformableMirrorSurfaceOwner)
    copyto_backend_async!(
        command_storage(owner.deformable_mirror),
        owner.pdm_command,
    )
    _update_dm_surface_async!(owner.deformable_mirror)
    copyto_backend_async!(
        owner.output,
        surface_opd(owner.deformable_mirror),
    )
    return nothing
end

@inline function reset_graph_node!(owner::_DeformableMirrorSurfaceOwner)
    fill!(surface_opd(owner.deformable_mirror), zero(eltype(owner.output)))
    fill!(owner.output, zero(eltype(owner.output)))
    return nothing
end

struct GaussianDeformableMirrorSurfaceNode{T<:AbstractFloat} end

"""Construction values for one analytic Gaussian deformable-mirror node."""
struct GaussianDeformableMirrorSurfaceNodeConfig{TD,T,AM}
    telescope::TD
    actuator_count::Int
    influence_width::T
    actuator_model::AM
    pdm_command_schema::String
    surface_opd_schema::String
    actuator_coordinates_schema::String
end

function _gaussian_deformable_mirror_surface_config(
    ::Type{T};
    resolution::Integer,
    telescope_diameter_m::Real,
    central_obstruction_ratio::Real,
    pupil_reflectivity::Real,
    aperture_revision::Integer,
    actuator_count::Integer,
    influence_width::Real,
    actuator_model,
    pdm_command_schema::AbstractString,
    surface_opd_schema::AbstractString,
    actuator_coordinates_schema::AbstractString,
) where {T<:AbstractFloat}
    count = Int(actuator_count)
    count > 0 || throw(AlgorithmGraphError(
        "Gaussian deformable-mirror actuator_count must be greater than zero",
    ))
    width = T(influence_width)
    isfinite(width) && width > zero(T) || throw(AlgorithmGraphError(
        "Gaussian deformable-mirror influence_width must be finite and " *
        "greater than zero",
    ))
    schemas = (
        pdm_command_schema,
        surface_opd_schema,
        actuator_coordinates_schema,
    )
    all(schema -> !isempty(schema), schemas) || throw(AlgorithmGraphError(
        "Gaussian deformable-mirror scientific schemas must not be empty",
    ))
    telescope = try
        TelescopeDefinition(
            resolution=Int(resolution),
            diameter=telescope_diameter_m,
            central_obstruction=central_obstruction_ratio,
            pupil_reflectivity=pupil_reflectivity,
            revision=Int(aperture_revision),
            T=T,
        )
    catch error
        error isa AdaptiveOpticsSimError || rethrow()
        throw(AlgorithmGraphError(sprint(showerror, error)))
    end
    return GaussianDeformableMirrorSurfaceNodeConfig(
        telescope,
        count,
        width,
        actuator_model,
        String(pdm_command_schema),
        String(surface_opd_schema),
        String(actuator_coordinates_schema),
    )
end

"""
    gaussian_deformable_mirror_surface_node(name; resolution,
        telescope_diameter_m, actuator_count, influence_width,
        pdm_command_schema, surface_opd_schema,
        actuator_coordinates_schema, ...)

Declare one complete-command analytic Gaussian deformable-mirror operation.
The graph snapshots caller-supplied actuator coordinates in normalized pupil
coordinates and evaluates the native matrix-free Gaussian influence model.
`influence_width` uses the same normalized-pupil coordinate system.

Each PDM command element is the unit-peak surface-OPD coefficient, in metres,
of its actuator's Gaussian influence function. The total output superposes all
actuator contributions and is also surface OPD in metres. A normalized hardware
command requires an explicit, separately qualified command-calibration adapter.
Use [`deformable_mirror_surface_node`](@ref) when sampled measured influence
functions are available.
"""
function gaussian_deformable_mirror_surface_node(
    name::Symbol;
    resolution::Integer,
    telescope_diameter_m::Real,
    actuator_count::Integer,
    influence_width::Real,
    pdm_command_schema::AbstractString,
    surface_opd_schema::AbstractString,
    actuator_coordinates_schema::AbstractString,
    central_obstruction_ratio::Real=0,
    pupil_reflectivity::Real=1,
    aperture_revision::Integer=0,
    actuator_model=nothing,
    T::Type{<:AbstractFloat}=Float32,
)
    config = _gaussian_deformable_mirror_surface_config(
        T;
        resolution,
        telescope_diameter_m,
        central_obstruction_ratio,
        pupil_reflectivity,
        aperture_revision,
        actuator_count,
        influence_width,
        actuator_model,
        pdm_command_schema,
        surface_opd_schema,
        actuator_coordinates_schema,
    )
    return algorithm_node(
        name,
        GaussianDeformableMirrorSurfaceNode{T},
        config;
        props=NamedTuple(),
    )
end

function graph_node_ports(
    ::Type{GaussianDeformableMirrorSurfaceNode{T}},
    config::GaussianDeformableMirrorSurfaceNodeConfig,
) where {T}
    resolution = config.telescope.resolution
    return (
        graph_port_contract(
            :pdm_command,
            :input,
            :data,
            T,
            (config.actuator_count,),
            config.pdm_command_schema,
            :column_major,
        ),
        graph_port_contract(
            :surface_opd,
            :output,
            :data,
            T,
            (resolution, resolution),
            config.surface_opd_schema,
            :column_major,
        ),
        graph_port_contract(
            :actuator_coordinates,
            :input,
            :parameter,
            T,
            (2, config.actuator_count),
            config.actuator_coordinates_schema,
            :column_major,
        ),
    )
end

function prepare_graph_node(
    ::Type{GaussianDeformableMirrorSurfaceNode{T}},
    config::GaussianDeformableMirrorSurfaceNodeConfig,
    ::NamedTuple{()},
    inputs::NamedTuple{(:pdm_command,)},
    outputs::NamedTuple{(:surface_opd,)},
    parameters::NamedTuple{(:actuator_coordinates,)},
    target,
) where {T}
    actuator_coordinates = Array(parameters.actuator_coordinates)
    all(isfinite, actuator_coordinates) || throw(AlgorithmGraphError(
        "Gaussian deformable-mirror actuator coordinates must be finite",
    ))
    metadata = (
        schema=config.actuator_coordinates_schema,
        coordinate_system=NormalizedPupilCoordinates(),
    )
    topology = SampledActuatorTopology(
        actuator_coordinates;
        metadata,
        T=T,
    )
    telescope = prepare_telescope(config.telescope, target)
    deformable_mirror = DeformableMirror(
        telescope;
        topology,
        influence_model=GaussianInfluenceWidth(config.influence_width),
        actuator_model=config.actuator_model,
        T=T,
        backend=compute_device_backend(target),
    )
    return _DeformableMirrorSurfaceOwner(
        deformable_mirror,
        inputs.pdm_command,
        outputs.surface_opd,
    )
end

struct GridGaussianDeformableMirrorSurfaceNode{T<:AbstractFloat} end

"""Construction values for one separable regular-grid Gaussian DM node."""
struct GridGaussianDeformableMirrorSurfaceNodeConfig{TD,T<:AbstractFloat}
    telescope::TD
    actuator_count::Int
    actuator_axis_count::Int
    actuator_pitch::T
    influence_width::T
    pdm_command_schema::String
    surface_opd_schema::String
    actuator_grid_indices_schema::String
end

function _grid_gaussian_deformable_mirror_surface_config(
    ::Type{T};
    resolution::Integer,
    telescope_diameter_m::Real,
    central_obstruction_ratio::Real,
    pupil_reflectivity::Real,
    aperture_revision::Integer,
    actuator_count::Integer,
    actuator_axis_count::Integer,
    actuator_pitch::Real,
    influence_width::Real,
    pdm_command_schema::AbstractString,
    surface_opd_schema::AbstractString,
    actuator_grid_indices_schema::AbstractString,
) where {T<:AbstractFloat}
    actuator_count isa Bool && throw(AlgorithmGraphError(
        "grid Gaussian deformable-mirror actuator_count must be an " *
        "integer, not Bool",
    ))
    0 < actuator_count <= typemax(Int) || throw(AlgorithmGraphError(
        "grid Gaussian deformable-mirror actuator_count must be positive " *
        "and fit the host index range",
    ))
    count = Int(actuator_count)
    actuator_axis_count isa Bool && throw(AlgorithmGraphError(
        "grid Gaussian deformable-mirror actuator_axis_count must be an " *
        "integer, not Bool",
    ))
    1 < actuator_axis_count <= isqrt(typemax(Int32)) || throw(
        AlgorithmGraphError(
            "grid Gaussian deformable-mirror actuator_axis_count must be " *
            "greater than one and its square must fit Int32",
        ),
    )
    axis_count = Int(actuator_axis_count)
    count <= axis_count * axis_count || throw(AlgorithmGraphError(
        "grid Gaussian deformable-mirror actuator_count must not exceed " *
        "actuator_axis_count squared",
    ))
    actuator_pitch isa Bool && throw(AlgorithmGraphError(
        "grid Gaussian deformable-mirror actuator_pitch must be a real " *
        "value, not Bool",
    ))
    pitch = T(actuator_pitch)
    isfinite(pitch) && pitch > zero(T) || throw(AlgorithmGraphError(
        "grid Gaussian deformable-mirror actuator_pitch must be finite " *
        "and greater than zero",
    ))
    width = T(influence_width)
    isfinite(width) && width > zero(T) || throw(AlgorithmGraphError(
        "grid Gaussian deformable-mirror influence_width must be finite " *
        "and greater than zero",
    ))
    schemas = (
        pdm_command_schema,
        surface_opd_schema,
        actuator_grid_indices_schema,
    )
    all(schema -> !isempty(schema), schemas) || throw(AlgorithmGraphError(
        "grid Gaussian deformable-mirror scientific schemas must not be " *
        "empty",
    ))
    telescope = try
        TelescopeDefinition(
            resolution=Int(resolution),
            diameter=telescope_diameter_m,
            central_obstruction=central_obstruction_ratio,
            pupil_reflectivity=pupil_reflectivity,
            revision=Int(aperture_revision),
            T=T,
        )
    catch error
        error isa AdaptiveOpticsSimError || rethrow()
        throw(AlgorithmGraphError(sprint(showerror, error)))
    end
    return GridGaussianDeformableMirrorSurfaceNodeConfig(
        telescope,
        count,
        axis_count,
        pitch,
        width,
        String(pdm_command_schema),
        String(surface_opd_schema),
        String(actuator_grid_indices_schema),
    )
end

"""
    grid_gaussian_deformable_mirror_surface_node(name; resolution,
        telescope_diameter_m, actuator_count, actuator_axis_count,
        actuator_pitch, influence_width, pdm_command_schema,
        surface_opd_schema, actuator_grid_indices_schema, ...)

Declare an analytic Gaussian deformable-mirror operation whose active
actuators occupy a subset of a square regular grid. The startup
`actuator_grid_indices` parameter maps each complete PDM command element to
one column-major cell of that grid.

The prepared node scatters the active command into a zero-filled grid and
evaluates the Gaussian surface as `X * C * Y'`. This factorization is
mathematically equivalent to summing the declared Gaussian influence function
at every mapped actuator coordinate. It is not valid for an irregular grid or
for a registration that breaks separability.
"""
function grid_gaussian_deformable_mirror_surface_node(
    name::Symbol;
    resolution::Integer,
    telescope_diameter_m::Real,
    actuator_count::Integer,
    actuator_axis_count::Integer,
    actuator_pitch::Real,
    influence_width::Real,
    pdm_command_schema::AbstractString,
    surface_opd_schema::AbstractString,
    actuator_grid_indices_schema::AbstractString,
    central_obstruction_ratio::Real=0,
    pupil_reflectivity::Real=1,
    aperture_revision::Integer=0,
    T::Type{<:AbstractFloat}=Float32,
)
    config = _grid_gaussian_deformable_mirror_surface_config(
        T;
        resolution,
        telescope_diameter_m,
        central_obstruction_ratio,
        pupil_reflectivity,
        aperture_revision,
        actuator_count,
        actuator_axis_count,
        actuator_pitch,
        influence_width,
        pdm_command_schema,
        surface_opd_schema,
        actuator_grid_indices_schema,
    )
    return algorithm_node(
        name,
        GridGaussianDeformableMirrorSurfaceNode{T},
        config;
        props=NamedTuple(),
    )
end

function graph_node_ports(
    ::Type{GridGaussianDeformableMirrorSurfaceNode{T}},
    config::GridGaussianDeformableMirrorSurfaceNodeConfig,
) where {T}
    resolution = config.telescope.resolution
    return (
        graph_port_contract(
            :pdm_command,
            :input,
            :data,
            T,
            (config.actuator_count,),
            config.pdm_command_schema,
            :column_major,
        ),
        graph_port_contract(
            :surface_opd,
            :output,
            :data,
            T,
            (resolution, resolution),
            config.surface_opd_schema,
            :column_major,
        ),
        graph_port_contract(
            :actuator_grid_indices,
            :input,
            :parameter,
            Int32,
            (config.actuator_count,),
            config.actuator_grid_indices_schema,
            :column_major,
        ),
    )
end

struct _GridGaussianDeformableMirrorSurfaceOwner{
    DM,
    Command,
    Indices,
    Output,
}
    deformable_mirror::DM
    pdm_command::Command
    actuator_grid_indices::Indices
    output::Output
end

@kernel function _scatter_grid_dm_command_kernel!(
    grid_command,
    active_command,
    grid_indices,
    count::Int,
)
    index = @index(Global, Linear)
    if index <= count
        @inbounds grid_command[grid_indices[index]] = active_command[index]
    end
end

@inline function _scatter_grid_dm_command!(
    ::ScalarCPUStyle,
    grid_command,
    active_command,
    grid_indices,
)
    fill!(grid_command, zero(eltype(grid_command)))
    @inbounds for index in eachindex(active_command, grid_indices)
        grid_command[grid_indices[index]] = active_command[index]
    end
    return grid_command
end

@inline function _scatter_grid_dm_command!(
    style::AcceleratorStyle,
    grid_command,
    active_command,
    grid_indices,
)
    fill!(grid_command, zero(eltype(grid_command)))
    launch_kernel!(
        style,
        _scatter_grid_dm_command_kernel!,
        grid_command,
        active_command,
        grid_indices,
        length(active_command);
        ndrange=length(active_command),
    )
    return grid_command
end

@inline function _scatter_grid_dm_command_async!(
    ::ScalarCPUStyle,
    grid_command,
    active_command,
    grid_indices,
)
    return _scatter_grid_dm_command!(ScalarCPUStyle(), grid_command,
        active_command, grid_indices)
end

@inline function _scatter_grid_dm_command_async!(
    style::AcceleratorStyle,
    grid_command,
    active_command,
    grid_indices,
)
    fill!(grid_command, zero(eltype(grid_command)))
    launch_kernel_async!(
        style,
        _scatter_grid_dm_command_kernel!,
        grid_command,
        active_command,
        grid_indices,
        length(active_command);
        ndrange=length(active_command),
    )
    return grid_command
end

function prepare_graph_node(
    ::Type{GridGaussianDeformableMirrorSurfaceNode{T}},
    config::GridGaussianDeformableMirrorSurfaceNodeConfig,
    ::NamedTuple{()},
    inputs::NamedTuple{(:pdm_command,)},
    outputs::NamedTuple{(:surface_opd,)},
    parameters::NamedTuple{(:actuator_grid_indices,)},
    target,
) where {T}
    grid_indices_host = Int32.(Array(parameters.actuator_grid_indices))
    grid_extent = config.actuator_axis_count * config.actuator_axis_count
    all(index -> 1 <= index <= grid_extent, grid_indices_host) || throw(
        AlgorithmGraphError(
            "grid Gaussian deformable-mirror actuator_grid_indices must " *
            "address the declared actuator grid",
        ),
    )
    length(unique(grid_indices_host)) == config.actuator_count || throw(
        AlgorithmGraphError(
            "grid Gaussian deformable-mirror actuator_grid_indices must " *
            "be unique",
        ),
    )
    grid_indices = similar(
        inputs.pdm_command,
        Int32,
        config.actuator_count,
    )
    copyto!(grid_indices, grid_indices_host)
    topology = ActuatorGridTopology(
        config.actuator_axis_count;
        actuator_pitch=config.actuator_pitch,
        T=T,
    )
    telescope = prepare_telescope(config.telescope, target)
    deformable_mirror = DeformableMirror(
        telescope;
        topology,
        influence_model=GaussianInfluenceWidth(config.influence_width),
        T=T,
        backend=compute_device_backend(target),
    )
    return _GridGaussianDeformableMirrorSurfaceOwner(
        deformable_mirror,
        inputs.pdm_command,
        grid_indices,
        outputs.surface_opd,
    )
end

@inline function step_graph_node!(
    owner::_GridGaussianDeformableMirrorSurfaceOwner,
)
    grid_command = command_storage(owner.deformable_mirror)
    _scatter_grid_dm_command!(
        execution_style(grid_command),
        grid_command,
        owner.pdm_command,
        owner.actuator_grid_indices,
    )
    update_surface!(owner.deformable_mirror)
    copyto!(owner.output, surface_opd(owner.deformable_mirror))
    return nothing
end

@inline function enqueue_graph_node!(
    owner::_GridGaussianDeformableMirrorSurfaceOwner,
)
    grid_command = command_storage(owner.deformable_mirror)
    _scatter_grid_dm_command_async!(
        execution_style(grid_command),
        grid_command,
        owner.pdm_command,
        owner.actuator_grid_indices,
    )
    _update_dm_surface_async!(owner.deformable_mirror)
    copyto_backend_async!(
        owner.output,
        surface_opd(owner.deformable_mirror),
    )
    return nothing
end

@inline graph_node_capture_capability(
    ::_GridGaussianDeformableMirrorSurfaceOwner,
) = GraphNodeCaptureSafe()

@inline function reset_graph_node!(
    owner::_GridGaussianDeformableMirrorSurfaceOwner,
)
    fill!(
        command_storage(owner.deformable_mirror),
        zero(eltype(owner.pdm_command)),
    )
    fill!(
        surface_opd(owner.deformable_mirror),
        zero(eltype(owner.output)),
    )
    fill!(owner.output, zero(eltype(owner.output)))
    return nothing
end

struct PupilOPDCompositionNode{T<:AbstractFloat} end

"""Construction values for one additive pupil-OPD composition node."""
struct PupilOPDCompositionNodeConfig
    resolution::Int
    uncompensated_opd_schema::String
    surface_opd_schema::String
    pupil_opd_schema::String

    function PupilOPDCompositionNodeConfig(
        resolution::Integer,
        uncompensated_opd_schema::AbstractString,
        surface_opd_schema::AbstractString,
        pupil_opd_schema::AbstractString,
    )
        resolution isa Bool && throw(AlgorithmGraphError(
            "pupil-OPD composition resolution must be an integer, not Bool",
        ))
        0 < resolution <= typemax(Int) || throw(AlgorithmGraphError(
            "pupil-OPD composition resolution must be positive and fit the " *
            "host index range",
        ))
        schemas = (
            uncompensated_opd_schema,
            surface_opd_schema,
            pupil_opd_schema,
        )
        all(schema -> !isempty(schema), schemas) || throw(AlgorithmGraphError(
            "pupil-OPD composition scientific schemas must not be empty",
        ))
        return new(
            Int(resolution),
            String(uncompensated_opd_schema),
            String(surface_opd_schema),
            String(pupil_opd_schema),
        )
    end
end

"""
    pupil_opd_composition_node(name; resolution, uncompensated_opd_schema,
                               surface_opd_schema, pupil_opd_schema,
                               T=Float32)

Declare one same-grid additive optical-path operation. Each step adds one
complete deformable-mirror surface OPD to the uncompensated pupil OPD and
publishes one complete pupil-plane OPD map. The node does not infer surface
sign, path placement, conjugation, resampling, or command timing.
"""
function pupil_opd_composition_node(
    name::Symbol;
    resolution::Integer,
    uncompensated_opd_schema::AbstractString,
    surface_opd_schema::AbstractString,
    pupil_opd_schema::AbstractString,
    T::Type{<:AbstractFloat}=Float32,
)
    config = PupilOPDCompositionNodeConfig(
        resolution,
        uncompensated_opd_schema,
        surface_opd_schema,
        pupil_opd_schema,
    )
    return algorithm_node(
        name,
        PupilOPDCompositionNode{T},
        config;
        props=NamedTuple(),
    )
end

function graph_node_ports(
    ::Type{PupilOPDCompositionNode{T}},
    config::PupilOPDCompositionNodeConfig,
) where {T}
    format(name, direction, schema) = graph_port_contract(
        name,
        direction,
        :data,
        T,
        (config.resolution, config.resolution),
        schema,
        :column_major,
    )
    return (
        format(
            :uncompensated_opd,
            :input,
            config.uncompensated_opd_schema,
        ),
        format(:surface_opd, :input, config.surface_opd_schema),
        format(:pupil_opd, :output, config.pupil_opd_schema),
    )
end

struct _PupilOPDCompositionOwner{Uncompensated,Surface,Output}
    uncompensated_opd::Uncompensated
    surface_opd::Surface
    pupil_opd::Output
end

function prepare_graph_node(
    ::Type{PupilOPDCompositionNode{T}},
    ::PupilOPDCompositionNodeConfig,
    ::NamedTuple{()},
    inputs::NamedTuple{(:uncompensated_opd,:surface_opd)},
    outputs::NamedTuple{(:pupil_opd,)},
    ::NamedTuple{()},
    target,
) where {T}
    return _PupilOPDCompositionOwner(
        inputs.uncompensated_opd,
        inputs.surface_opd,
        outputs.pupil_opd,
    )
end

@inline function step_graph_node!(owner::_PupilOPDCompositionOwner)
    @. owner.pupil_opd = owner.uncompensated_opd + owner.surface_opd
    return nothing
end

@inline graph_node_capture_capability(
    ::_PupilOPDCompositionOwner,
) = GraphNodeCaptureSafe()

@inline function reset_graph_node!(owner::_PupilOPDCompositionOwner)
    fill!(owner.pupil_opd, zero(eltype(owner.pupil_opd)))
    return nothing
end

struct ShackHartmannRateNode{T<:AbstractFloat} end

"""Construction values for one diffractive Shack–Hartmann photon-rate node."""
struct ShackHartmannRateNodeConfig{T<:AbstractFloat,TD,S}
    telescope::TD
    source::S
    n_lenslets::Int
    n_pix_subap::Int
    diffraction_padding::Int
    pixel_scale_arcsec::T
    valid_subaperture_threshold::T
    threshold_convolution::T
    half_pixel_shift::Bool
    shannon_sampling::Bool
    opd_schema::String
    photon_rate_schema::String
end

function _shack_hartmann_rate_config(
    ::Type{T};
    resolution::Integer,
    telescope_diameter_m::Real,
    central_obstruction_ratio::Real,
    pupil_reflectivity::Real,
    aperture_revision::Integer,
    n_lenslets::Integer,
    n_pix_subap::Integer,
    diffraction_padding::Integer,
    pixel_scale_arcsec::Real,
    valid_subaperture_threshold::Real,
    threshold_convolution::Real,
    half_pixel_shift::Bool,
    shannon_sampling::Bool,
    source_band::Symbol,
    source_magnitude::Real,
    source_wavelength_m::Real,
    source_photon_irradiance_m2_s::Union{Nothing,Real},
    source_separation_arcsec::Real,
    source_position_angle_deg::Real,
    opd_schema::AbstractString,
    photon_rate_schema::AbstractString,
) where {T<:AbstractFloat}
    lenslets = Int(n_lenslets)
    pixels = Int(n_pix_subap)
    padding = Int(diffraction_padding)
    pupil_resolution = Int(resolution)
    lenslets > 0 || throw(AlgorithmGraphError(
        "Shack–Hartmann n_lenslets must be greater than zero",
    ))
    pupil_resolution > 0 || throw(AlgorithmGraphError(
        "Shack–Hartmann resolution must be greater than zero",
    ))
    pupil_resolution % lenslets == 0 || throw(AlgorithmGraphError(
        "Shack–Hartmann resolution must be divisible by n_lenslets",
    ))
    pixels > 0 && iseven(pixels) || throw(AlgorithmGraphError(
        "Shack–Hartmann n_pix_subap must be positive and even",
    ))
    padding > 0 || throw(AlgorithmGraphError(
        "Shack–Hartmann diffraction_padding must be greater than zero",
    ))
    isempty(opd_schema) && throw(AlgorithmGraphError(
        "Shack–Hartmann opd_schema must not be empty",
    ))
    isempty(photon_rate_schema) && throw(AlgorithmGraphError(
        "Shack–Hartmann photon_rate_schema must not be empty",
    ))

    pixel_scale = T(pixel_scale_arcsec)
    valid_threshold = T(valid_subaperture_threshold)
    convolution_threshold = T(threshold_convolution)
    isfinite(pixel_scale) && pixel_scale > zero(T) || throw(
        AlgorithmGraphError(
            "Shack–Hartmann pixel_scale_arcsec must be finite and positive",
        ))
    isfinite(valid_threshold) && zero(T) <= valid_threshold <= one(T) ||
        throw(AlgorithmGraphError(
            "Shack–Hartmann valid_subaperture_threshold must lie in [0, 1]",
        ))
    isfinite(convolution_threshold) &&
        zero(T) <= convolution_threshold <= one(T) || throw(
            AlgorithmGraphError(
                "Shack–Hartmann threshold_convolution must lie in [0, 1]",
            ))

    telescope, source = try
        definition = TelescopeDefinition(
            resolution=pupil_resolution,
            diameter=telescope_diameter_m,
            central_obstruction=central_obstruction_ratio,
            pupil_reflectivity=pupil_reflectivity,
            revision=aperture_revision,
            T=T,
        )
        guide_source = Source(
            band=source_band,
            magnitude=source_magnitude,
            coordinates=(source_separation_arcsec, source_position_angle_deg),
            wavelength=source_wavelength_m,
            photon_irradiance=source_photon_irradiance_m2_s,
            T=T,
        )
        photon_irradiance(guide_source)
        (definition, guide_source)
    catch error
        error isa AdaptiveOpticsSimError || rethrow()
        throw(AlgorithmGraphError(sprint(showerror, error)))
    end

    return ShackHartmannRateNodeConfig(
        telescope,
        source,
        lenslets,
        pixels,
        padding,
        pixel_scale,
        valid_threshold,
        convolution_threshold,
        half_pixel_shift,
        shannon_sampling,
        String(opd_schema),
        String(photon_rate_schema),
    )
end

"""
    shack_hartmann_rate_node(name; ...)

Declare one complete-frame diffractive Shack–Hartmann optics node. The node
consumes a pupil-plane OPD map and writes a cell-integrated detector-plane
photon-rate mosaic. Detector exposure/noise and slope estimation remain
separate graph operations.
"""
function shack_hartmann_rate_node(
    name::Symbol;
    resolution::Integer,
    telescope_diameter_m::Real,
    n_lenslets::Integer,
    n_pix_subap::Integer,
    pixel_scale_arcsec::Real,
    source_wavelength_m::Real,
    source_photon_irradiance_m2_s::Union{Nothing,Real}=nothing,
    opd_schema::AbstractString,
    photon_rate_schema::AbstractString,
    central_obstruction_ratio::Real=0,
    pupil_reflectivity::Real=1,
    aperture_revision::Integer=0,
    diffraction_padding::Integer=2,
    valid_subaperture_threshold::Real=0.1,
    threshold_convolution::Real=0.05,
    half_pixel_shift::Bool=false,
    shannon_sampling::Bool=true,
    source_band::Symbol=:custom,
    source_magnitude::Real=0,
    source_separation_arcsec::Real=0,
    source_position_angle_deg::Real=0,
    T::Type{<:AbstractFloat}=Float32,
)
    config = _shack_hartmann_rate_config(
        T;
        resolution,
        telescope_diameter_m,
        central_obstruction_ratio,
        pupil_reflectivity,
        aperture_revision,
        n_lenslets,
        n_pix_subap,
        diffraction_padding,
        pixel_scale_arcsec,
        valid_subaperture_threshold,
        threshold_convolution,
        half_pixel_shift,
        shannon_sampling,
        source_band,
        source_magnitude,
        source_wavelength_m,
        source_photon_irradiance_m2_s,
        source_separation_arcsec,
        source_position_angle_deg,
        opd_schema,
        photon_rate_schema,
    )
    return algorithm_node(
        name,
        ShackHartmannRateNode{T},
        config;
        props=NamedTuple(),
    )
end

function graph_node_ports(
    ::Type{ShackHartmannRateNode{T}},
    config::ShackHartmannRateNodeConfig{T},
) where {T}
    output_extent = config.n_lenslets * config.n_pix_subap
    return (
        graph_port_contract(
            :opd,
            :input,
            :data,
            T,
            (config.telescope.resolution, config.telescope.resolution),
            config.opd_schema,
            :column_major,
        ),
        graph_port_contract(
            :photon_rate,
            :output,
            :data,
            T,
            (output_extent, output_extent),
            config.photon_rate_schema,
            :column_major,
        ),
    )
end

struct _ShackHartmannRateOwner{P}
    prepared::P
end

function prepare_graph_node(
    ::Type{ShackHartmannRateNode{T}},
    config::ShackHartmannRateNodeConfig{T},
    ::NamedTuple{()},
    inputs::NamedTuple{(:opd,)},
    outputs::NamedTuple{(:photon_rate,)},
    ::NamedTuple{()},
    target,
) where {T}
    telescope = prepare_telescope(config.telescope, target)
    pupil = PupilFunction(telescope, inputs.opd)
    sensor = ShackHartmannWFS(
        telescope;
        n_lenslets=config.n_lenslets,
        threshold=config.valid_subaperture_threshold,
        threshold_convolution=config.threshold_convolution,
        half_pixel_shift=config.half_pixel_shift,
        diffraction_padding=config.diffraction_padding,
        pixel_scale_arcsec=config.pixel_scale_arcsec,
        n_pix_subap=config.n_pix_subap,
        shannon_sampling=config.shannon_sampling,
        mode=Diffractive(),
        T=T,
        backend=compute_device_backend(target),
    )
    optics = shack_hartmann_optics(sensor, config.source)
    photon_rate = shack_hartmann_rate_map(
        optics,
        pupil,
        outputs.photon_rate,
    )
    prepared = prepare_wfs_optics(optics, pupil, photon_rate)
    return _ShackHartmannRateOwner(prepared)
end

@inline function step_graph_node!(owner::_ShackHartmannRateOwner)
    prepared = owner.prepared
    form_wfs_optical_products!(
        prepared.output,
        prepared.input,
        prepared,
    )
    return nothing
end

@inline reset_graph_node!(::_ShackHartmannRateOwner) = nothing

struct PyramidRateNode{T<:AbstractFloat} end

"""Construction values for one diffractive Pyramid photon-rate node."""
struct PyramidRateNodeConfig{T<:AbstractFloat,TD,S,MS}
    telescope::TD
    source::S
    pupil_samples::Int
    threshold::T
    modulation::T
    modulation_points::Union{Nothing,Int}
    modulation_propagation_strategy::MS
    light_ratio::T
    diffraction_padding::Int
    psf_centering::Bool
    n_pix_separation::Union{Nothing,Int}
    n_pix_edge::Union{Nothing,Int}
    binning::Int
    opd_schema::String
    photon_rate_schema::String
end

function _pyramid_rate_config(
    ::Type{T};
    resolution::Integer,
    telescope_diameter_m::Real,
    central_obstruction_ratio::Real,
    pupil_reflectivity::Real,
    aperture_revision::Integer,
    pupil_samples::Integer,
    threshold::Real,
    modulation::Real,
    modulation_points::Union{Nothing,Integer},
    modulation_propagation_strategy::
        AbstractPyramidModulationPropagationStrategy,
    light_ratio::Real,
    diffraction_padding::Integer,
    psf_centering::Bool,
    n_pix_separation::Union{Nothing,Integer},
    n_pix_edge::Union{Nothing,Integer},
    binning::Integer,
    source_band::Symbol,
    source_magnitude::Real,
    source_wavelength_m::Real,
    source_photon_irradiance_m2_s::Union{Nothing,Real},
    source_separation_arcsec::Real,
    source_position_angle_deg::Real,
    opd_schema::AbstractString,
    photon_rate_schema::AbstractString,
) where {T<:AbstractFloat}
    pupil_resolution = Int(resolution)
    samples = Int(pupil_samples)
    padding = Int(diffraction_padding)
    bin = Int(binning)
    points = isnothing(modulation_points) ? nothing : Int(modulation_points)
    separation = isnothing(n_pix_separation) ? nothing :
        Int(n_pix_separation)
    edge = isnothing(n_pix_edge) ? nothing : Int(n_pix_edge)

    pupil_resolution > 0 || throw(AlgorithmGraphError(
        "Pyramid resolution must be greater than zero",
    ))
    samples > 0 || throw(AlgorithmGraphError(
        "Pyramid pupil_samples must be greater than zero",
    ))
    pupil_resolution % samples == 0 || throw(AlgorithmGraphError(
        "Pyramid resolution must be divisible by pupil_samples",
    ))
    padding > 0 || throw(AlgorithmGraphError(
        "Pyramid diffraction_padding must be greater than zero",
    ))
    bin > 0 && samples % bin == 0 || throw(AlgorithmGraphError(
        "Pyramid binning must be positive and divide pupil_samples",
    ))
    isnothing(points) || points > 0 || throw(AlgorithmGraphError(
        "Pyramid modulation_points must be greater than zero",
    ))
    if isnothing(separation)
        isnothing(edge) || throw(AlgorithmGraphError(
            "Pyramid n_pix_edge requires n_pix_separation",
        ))
    else
        separation >= 0 || throw(AlgorithmGraphError(
            "Pyramid n_pix_separation must be nonnegative",
        ))
        separation % (2 * bin) == 0 || throw(AlgorithmGraphError(
            "Pyramid n_pix_separation must remain even after binning",
        ))
        resolved_edge = isnothing(edge) ? div(separation, 2) : edge
        resolved_edge >= 0 || throw(AlgorithmGraphError(
            "Pyramid n_pix_edge must be nonnegative",
        ))
        resolved_edge % bin == 0 || throw(AlgorithmGraphError(
            "Pyramid n_pix_edge must remain integral after binning",
        ))
    end

    typed_threshold = T(threshold)
    typed_modulation = T(modulation)
    typed_light_ratio = T(light_ratio)
    isfinite(typed_threshold) &&
        zero(T) <= typed_threshold <= one(T) || throw(
            AlgorithmGraphError("Pyramid threshold must lie in [0, 1]"),
        )
    isfinite(typed_modulation) && typed_modulation >= zero(T) || throw(
        AlgorithmGraphError(
            "Pyramid modulation must be finite and nonnegative",
        ),
    )
    isfinite(typed_light_ratio) &&
        zero(T) <= typed_light_ratio <= one(T) || throw(
            AlgorithmGraphError("Pyramid light_ratio must lie in [0, 1]"),
        )
    isempty(opd_schema) && throw(AlgorithmGraphError(
        "Pyramid opd_schema must not be empty",
    ))
    isempty(photon_rate_schema) && throw(AlgorithmGraphError(
        "Pyramid photon_rate_schema must not be empty",
    ))

    telescope, source = try
        definition = TelescopeDefinition(
            resolution=pupil_resolution,
            diameter=telescope_diameter_m,
            central_obstruction=central_obstruction_ratio,
            pupil_reflectivity=pupil_reflectivity,
            revision=aperture_revision,
            T=T,
        )
        guide_source = Source(
            band=source_band,
            magnitude=source_magnitude,
            coordinates=(source_separation_arcsec, source_position_angle_deg),
            wavelength=source_wavelength_m,
            photon_irradiance=source_photon_irradiance_m2_s,
            T=T,
        )
        photon_irradiance(guide_source)
        (definition, guide_source)
    catch error
        error isa AdaptiveOpticsSimError || rethrow()
        throw(AlgorithmGraphError(sprint(showerror, error)))
    end

    return PyramidRateNodeConfig(
        telescope,
        source,
        samples,
        typed_threshold,
        typed_modulation,
        points,
        modulation_propagation_strategy,
        typed_light_ratio,
        padding,
        psf_centering,
        separation,
        edge,
        bin,
        String(opd_schema),
        String(photon_rate_schema),
    )
end

"""
    pyramid_rate_node(name; ...)

Declare one complete-frame diffractive Pyramid optics node. The node consumes
a pupil-plane OPD map and writes the four re-imaged pupils as one
cell-integrated detector-plane photon-rate frame. Detector acquisition and
external RTC processing remain separate graph operations.
"""
function pyramid_rate_node(
    name::Symbol;
    resolution::Integer,
    telescope_diameter_m::Real,
    pupil_samples::Integer,
    source_wavelength_m::Real,
    opd_schema::AbstractString,
    photon_rate_schema::AbstractString,
    central_obstruction_ratio::Real=0,
    pupil_reflectivity::Real=1,
    aperture_revision::Integer=0,
    threshold::Real=0.1,
    modulation::Real=2,
    modulation_points::Union{Nothing,Integer}=nothing,
    modulation_propagation_strategy::
        AbstractPyramidModulationPropagationStrategy=
        PyramidPupilTiltStrategy(),
    light_ratio::Real=0,
    diffraction_padding::Integer=2,
    psf_centering::Bool=true,
    n_pix_separation::Union{Nothing,Integer}=nothing,
    n_pix_edge::Union{Nothing,Integer}=nothing,
    binning::Integer=1,
    source_band::Symbol=:custom,
    source_magnitude::Real=0,
    source_photon_irradiance_m2_s::Union{Nothing,Real}=nothing,
    source_separation_arcsec::Real=0,
    source_position_angle_deg::Real=0,
    T::Type{<:AbstractFloat}=Float32,
)
    config = _pyramid_rate_config(
        T;
        resolution,
        telescope_diameter_m,
        central_obstruction_ratio,
        pupil_reflectivity,
        aperture_revision,
        pupil_samples,
        threshold,
        modulation,
        modulation_points,
        modulation_propagation_strategy,
        light_ratio,
        diffraction_padding,
        psf_centering,
        n_pix_separation,
        n_pix_edge,
        binning,
        source_band,
        source_magnitude,
        source_wavelength_m,
        source_photon_irradiance_m2_s,
        source_separation_arcsec,
        source_position_angle_deg,
        opd_schema,
        photon_rate_schema,
    )
    return algorithm_node(
        name,
        PyramidRateNode{T},
        config;
        props=NamedTuple(),
    )
end

@inline function _pyramid_rate_output_extent(
    config::PyramidRateNodeConfig,
)
    if isnothing(config.n_pix_separation)
        return div(
            config.pupil_samples * config.diffraction_padding,
            config.binning,
        )
    end
    edge = isnothing(config.n_pix_edge) ?
        div(config.n_pix_separation, 2) : config.n_pix_edge
    return div(
        2 * config.pupil_samples + config.n_pix_separation + 2 * edge,
        config.binning,
    )
end

function graph_node_ports(
    ::Type{PyramidRateNode{T}},
    config::PyramidRateNodeConfig{T},
) where {T}
    output_extent = _pyramid_rate_output_extent(config)
    return (
        graph_port_contract(
            :opd,
            :input,
            :data,
            T,
            (config.telescope.resolution, config.telescope.resolution),
            config.opd_schema,
            :column_major,
        ),
        graph_port_contract(
            :photon_rate,
            :output,
            :data,
            T,
            (output_extent, output_extent),
            config.photon_rate_schema,
            :column_major,
        ),
    )
end

struct _PyramidRateOwner{P}
    prepared::P
end

function prepare_graph_node(
    ::Type{PyramidRateNode{T}},
    config::PyramidRateNodeConfig{T},
    ::NamedTuple{()},
    inputs::NamedTuple{(:opd,)},
    outputs::NamedTuple{(:photon_rate,)},
    ::NamedTuple{()},
    target,
) where {T}
    telescope = prepare_telescope(config.telescope, target)
    pupil = PupilFunction(telescope, inputs.opd)
    sensor = PyramidWFS(
        telescope;
        pupil_samples=config.pupil_samples,
        threshold=config.threshold,
        modulation=config.modulation,
        modulation_points=config.modulation_points,
        modulation_propagation_strategy=
            config.modulation_propagation_strategy,
        light_ratio=config.light_ratio,
        diffraction_padding=config.diffraction_padding,
        psf_centering=config.psf_centering,
        n_pix_separation=config.n_pix_separation,
        n_pix_edge=config.n_pix_edge,
        binning=config.binning,
        mode=Diffractive(),
        T=T,
        backend=compute_device_backend(target),
    )
    front_end = PyramidOpticalFrontEnd(sensor, config.source)
    photon_rate = pyramid_rate_map(front_end, pupil, outputs.photon_rate)
    prepared = prepare_wfs_optics(front_end, pupil, photon_rate)
    return _PyramidRateOwner(prepared)
end

@inline function step_graph_node!(owner::_PyramidRateOwner)
    prepared = owner.prepared
    form_wfs_optical_products!(
        prepared.output,
        prepared.input,
        prepared,
    )
    return nothing
end

@inline reset_graph_node!(::_PyramidRateOwner) = nothing

struct ShackHartmannCentroidNode{T<:AbstractFloat} end

"""Construction values for one calibrated Shack–Hartmann centroid node."""
struct ShackHartmannCentroidNodeConfig{T<:AbstractFloat,TD}
    telescope::TD
    n_lenslets::Int
    n_pix_subap::Int
    centroid_cutoff_fraction::T
    centroid_response::T
    calibration_wavelength_m::T
    calibration_signature::UInt
    frame_schema::String
    slopes_schema::String
    valid_subapertures_schema::String
    reference_signal_schema::String
end

function _shack_hartmann_centroid_config(
    ::Type{T};
    resolution::Integer,
    telescope_diameter_m::Real,
    n_lenslets::Integer,
    n_pix_subap::Integer,
    centroid_cutoff_fraction::Real,
    centroid_response::Real,
    calibration_wavelength_m::Real,
    calibration_signature::Integer,
    frame_schema::AbstractString,
    slopes_schema::AbstractString,
    valid_subapertures_schema::AbstractString,
    reference_signal_schema::AbstractString,
) where {T<:AbstractFloat}
    pupil_resolution = Int(resolution)
    lenslets = Int(n_lenslets)
    pixels = Int(n_pix_subap)
    pupil_resolution > 0 || throw(AlgorithmGraphError(
        "Shack–Hartmann centroid resolution must be greater than zero",
    ))
    lenslets > 0 || throw(AlgorithmGraphError(
        "Shack–Hartmann centroid n_lenslets must be greater than zero",
    ))
    pupil_resolution % lenslets == 0 || throw(AlgorithmGraphError(
        "Shack–Hartmann centroid resolution must be divisible by n_lenslets",
    ))
    pixels > 0 && iseven(pixels) || throw(AlgorithmGraphError(
        "Shack–Hartmann centroid n_pix_subap must be positive and even",
    ))

    cutoff_fraction = T(centroid_cutoff_fraction)
    response = T(centroid_response)
    wavelength = T(calibration_wavelength_m)
    isfinite(cutoff_fraction) &&
        zero(T) <= cutoff_fraction <= one(T) ||
        throw(AlgorithmGraphError(
            "Shack–Hartmann centroid_cutoff_fraction must lie in [0, 1]",
        ))
    isfinite(response) && response != zero(T) ||
        throw(AlgorithmGraphError(
            "Shack–Hartmann centroid_response must be finite and nonzero",
        ))
    isfinite(wavelength) && wavelength > zero(T) ||
        throw(AlgorithmGraphError(
            "Shack–Hartmann calibration_wavelength_m must be finite and positive",
        ))
    calibration_signature >= 0 || throw(AlgorithmGraphError(
        "Shack–Hartmann calibration_signature must be nonnegative",
    ))
    calibration_signature <= typemax(UInt) || throw(AlgorithmGraphError(
        "Shack–Hartmann calibration_signature exceeds UInt",
    ))
    schemas = (
        frame_schema,
        slopes_schema,
        valid_subapertures_schema,
        reference_signal_schema,
    )
    all(schema -> !isempty(schema), schemas) || throw(AlgorithmGraphError(
        "Shack–Hartmann centroid scientific schemas must not be empty",
    ))

    telescope = try
        TelescopeDefinition(
            resolution=pupil_resolution,
            diameter=telescope_diameter_m,
            revision=0,
            T=T,
        )
    catch error
        error isa AdaptiveOpticsSimError || rethrow()
        throw(AlgorithmGraphError(sprint(showerror, error)))
    end
    return ShackHartmannCentroidNodeConfig(
        telescope,
        lenslets,
        pixels,
        cutoff_fraction,
        response,
        wavelength,
        UInt(calibration_signature),
        String(frame_schema),
        String(slopes_schema),
        String(valid_subapertures_schema),
        String(reference_signal_schema),
    )
end

"""
    shack_hartmann_centroid_node(name; ...)

Declare one calibrated, complete-frame Shack–Hartmann centroid operation. The
node consumes a lenslet-mosaic detector frame and writes the full canonical
`[axis 1; axis 2]` centroid-slope vector. Its valid-subaperture mask and
reference signal are explicit startup sparse parameters; compact slope
selection remains a separate graph operation.
"""
function shack_hartmann_centroid_node(
    name::Symbol;
    resolution::Integer,
    telescope_diameter_m::Real,
    n_lenslets::Integer,
    n_pix_subap::Integer,
    centroid_cutoff_fraction::Real,
    centroid_response::Real,
    calibration_wavelength_m::Real,
    calibration_signature::Integer,
    frame_schema::AbstractString,
    slopes_schema::AbstractString,
    valid_subapertures_schema::AbstractString,
    reference_signal_schema::AbstractString,
    T::Type{<:AbstractFloat}=Float32,
)
    config = _shack_hartmann_centroid_config(
        T;
        resolution,
        telescope_diameter_m,
        n_lenslets,
        n_pix_subap,
        centroid_cutoff_fraction,
        centroid_response,
        calibration_wavelength_m,
        calibration_signature,
        frame_schema,
        slopes_schema,
        valid_subapertures_schema,
        reference_signal_schema,
    )
    return algorithm_node(
        name,
        ShackHartmannCentroidNode{T},
        config;
        props=NamedTuple(),
    )
end

function graph_node_ports(
    ::Type{ShackHartmannCentroidNode{T}},
    config::ShackHartmannCentroidNodeConfig{T},
) where {T}
    lenslet_count = config.n_lenslets * config.n_lenslets
    frame_extent = config.n_lenslets * config.n_pix_subap
    return (
        graph_port_contract(
            :frame,
            :input,
            :data,
            T,
            (frame_extent, frame_extent),
            config.frame_schema,
            :column_major,
        ),
        graph_port_contract(
            :slopes,
            :output,
            :data,
            T,
            (2 * lenslet_count,),
            config.slopes_schema,
            :column_major,
        ),
        graph_port_contract(
            :valid_subapertures,
            :input,
            :parameter,
            Bool,
            (config.n_lenslets, config.n_lenslets),
            config.valid_subapertures_schema,
            :column_major,
        ),
        graph_port_contract(
            :reference_signal,
            :input,
            :parameter,
            T,
            (lenslet_count, 2),
            config.reference_signal_schema,
            :column_major,
        ),
    )
end

struct _ShackHartmannCentroidOwner{P,O,M}
    prepared::P
    observation::O
    measurement::M
end

function prepare_graph_node(
    ::Type{ShackHartmannCentroidNode{T}},
    config::ShackHartmannCentroidNodeConfig{T},
    ::NamedTuple{()},
    inputs::NamedTuple{(:frame,)},
    outputs::NamedTuple{(:slopes,)},
    parameters::NamedTuple{(:valid_subapertures,:reference_signal)},
    target,
) where {T}
    telescope = prepare_telescope(config.telescope, target)
    sensor = ShackHartmannWFS(
        telescope;
        n_lenslets=config.n_lenslets,
        threshold_cog=config.centroid_cutoff_fraction,
        diffraction_padding=1,
        n_pix_subap=config.n_pix_subap,
        mode=Diffractive(),
        T=T,
        backend=compute_device_backend(target),
    )
    set_valid_subapertures!(sensor, parameters.valid_subapertures)
    set_subaperture_calibration!(
        subaperture_calibration(sensor),
        parameters.reference_signal;
        centroid_response=config.centroid_response,
        output_units=:pixel,
        wavelength=config.calibration_wavelength_m,
        signature=config.calibration_signature,
    )
    observation = WFSObservation(
        inputs.frame;
        units=:electron_count,
        layout=:lenslet_mosaic,
    )
    measurement = WFSMeasurement(
        outputs.slopes;
        units=:pixel,
        kind=:centroid_slopes,
    )
    prepared = prepare_wfs_estimation(sensor, observation, measurement)
    return _ShackHartmannCentroidOwner(
        prepared,
        observation,
        measurement,
    )
end

@inline function step_graph_node!(owner::_ShackHartmannCentroidOwner)
    estimate_wfs_measurement!(
        owner.measurement,
        owner.observation,
        owner.prepared,
    )
    return nothing
end

@inline function reset_graph_node!(owner::_ShackHartmannCentroidOwner)
    fill!(owner.measurement.storage, zero(eltype(owner.measurement.storage)))
    return nothing
end

struct ShackHartmannSlopeSelectionNode{T<:AbstractFloat} end

"""Construction values for one ordered Shack–Hartmann slope selector."""
struct ShackHartmannSlopeSelectionNodeConfig
    n_lenslets::Int
    selected_lenslet_count::Int
    full_slopes_schema::String
    selected_slopes_schema::String
    lenslet_order_schema::String
end

function ShackHartmannSlopeSelectionNodeConfig(
    n_lenslets::Integer,
    selected_lenslet_count::Integer,
    full_slopes_schema::AbstractString,
    selected_slopes_schema::AbstractString,
    lenslet_order_schema::AbstractString,
)
    lenslets = Int(n_lenslets)
    selected_count = Int(selected_lenslet_count)
    lenslets > 0 || throw(AlgorithmGraphError(
        "Shack–Hartmann slope-selection n_lenslets must be positive",
    ))
    full_count = Base.checked_mul(lenslets, lenslets)
    0 < selected_count <= full_count || throw(AlgorithmGraphError(
        "Shack–Hartmann selected_lenslet_count must lie in 1:$full_count",
    ))
    schemas = (
        full_slopes_schema,
        selected_slopes_schema,
        lenslet_order_schema,
    )
    all(schema -> !isempty(schema), schemas) || throw(AlgorithmGraphError(
        "Shack–Hartmann slope-selection schemas must not be empty",
    ))
    return ShackHartmannSlopeSelectionNodeConfig(
        lenslets,
        selected_count,
        String(full_slopes_schema),
        String(selected_slopes_schema),
        String(lenslet_order_schema),
    )
end

"""
    shack_hartmann_slope_selection_node(name; ...)

Declare one complete-frame slope selector. The node converts the full AOS
`[axis 1; axis 2]` slope vector into pair-interleaved `[axis 1, axis 2]`
pairs using an explicit startup lenslet order. Instrument ROI-file parsing
remains an application concern.
"""
function shack_hartmann_slope_selection_node(
    name::Symbol;
    n_lenslets::Integer,
    selected_lenslet_count::Integer,
    full_slopes_schema::AbstractString,
    selected_slopes_schema::AbstractString,
    lenslet_order_schema::AbstractString,
    T::Type{<:AbstractFloat}=Float32,
)
    config = ShackHartmannSlopeSelectionNodeConfig(
        n_lenslets,
        selected_lenslet_count,
        full_slopes_schema,
        selected_slopes_schema,
        lenslet_order_schema,
    )
    return algorithm_node(
        name,
        ShackHartmannSlopeSelectionNode{T},
        config;
        props=NamedTuple(),
    )
end

function graph_node_ports(
    ::Type{ShackHartmannSlopeSelectionNode{T}},
    config::ShackHartmannSlopeSelectionNodeConfig,
) where {T}
    full_count = config.n_lenslets * config.n_lenslets
    selected_count = config.selected_lenslet_count
    return (
        graph_port_contract(
            :full_slopes,
            :input,
            :data,
            T,
            (2 * full_count,),
            config.full_slopes_schema,
            :column_major,
        ),
        graph_port_contract(
            :selected_slopes,
            :output,
            :data,
            T,
            (2 * selected_count,),
            config.selected_slopes_schema,
            :column_major,
        ),
        graph_port_contract(
            :lenslet_order,
            :input,
            :parameter,
            UInt32,
            (selected_count,),
            config.lenslet_order_schema,
            :column_major,
        ),
    )
end

struct _ShackHartmannSlopeSelectionOwner{P,I,O}
    plan::P
    input::I
    output::O
end

function prepare_graph_node(
    ::Type{ShackHartmannSlopeSelectionNode{T}},
    config::ShackHartmannSlopeSelectionNodeConfig,
    ::NamedTuple{()},
    inputs::NamedTuple{(:full_slopes,)},
    outputs::NamedTuple{(:selected_slopes,)},
    parameters::NamedTuple{(:lenslet_order,)},
    target,
) where {T}
    plan = try
        ShackHartmannSlopeSelectionPlan(
            config.n_lenslets,
            parameters.lenslet_order,
        )
    catch error
        error isa AdaptiveOpticsSimError || rethrow()
        throw(AlgorithmGraphError(sprint(showerror, error)))
    end
    selected_lenslet_count(plan) == config.selected_lenslet_count || throw(
        AlgorithmGraphError(
            "Shack–Hartmann lenslet order does not match selected_lenslet_count",
        ),
    )
    return _ShackHartmannSlopeSelectionOwner(
        plan,
        inputs.full_slopes,
        outputs.selected_slopes,
    )
end

@inline function step_graph_node!(owner::_ShackHartmannSlopeSelectionOwner)
    select_shack_hartmann_slopes!(
        owner.output,
        owner.plan,
        owner.input,
    )
    return nothing
end

@inline function reset_graph_node!(owner::_ShackHartmannSlopeSelectionOwner)
    fill!(owner.output, zero(eltype(owner.output)))
    return nothing
end

struct ControlMatrixReconstructionNode{T<:AbstractFloat} end

"""Construction values for one calibrated dense reconstruction node."""
struct ControlMatrixReconstructionNodeConfig
    slope_count::Int
    reconstructed_count::Int
    slopes_schema::String
    reconstructed_schema::String
    control_matrix_schema::String

    function ControlMatrixReconstructionNodeConfig(
        slope_count::Integer,
        reconstructed_count::Integer,
        slopes_schema::AbstractString,
        reconstructed_schema::AbstractString,
        control_matrix_schema::AbstractString,
    )
        slope_count > 0 || throw(AlgorithmGraphError(
            "control-matrix reconstruction slope_count must be positive",
        ))
        reconstructed_count > 0 || throw(AlgorithmGraphError(
            "control-matrix reconstruction reconstructed_count must be positive",
        ))
        schemas = (
            slopes_schema,
            reconstructed_schema,
            control_matrix_schema,
        )
        all(schema -> !isempty(schema), schemas) || throw(AlgorithmGraphError(
            "control-matrix reconstruction schemas must not be empty",
        ))
        return new(
            Int(slope_count),
            Int(reconstructed_count),
            String(slopes_schema),
            String(reconstructed_schema),
            String(control_matrix_schema),
        )
    end
end

"""
    control_matrix_reconstruction_node(name; ...)

Declare one complete-frame dense reconstruction operation. The node snapshots
an already calibrated control matrix from a startup sparse parameter, consumes
one exact ordered slope vector, and writes reconstructed controller
coordinates. Slope ordering, reference subtraction, controller integration,
and command-basis expansion remain separate operations.
"""
function control_matrix_reconstruction_node(
    name::Symbol;
    slope_count::Integer,
    reconstructed_count::Integer,
    slopes_schema::AbstractString,
    reconstructed_schema::AbstractString,
    control_matrix_schema::AbstractString,
    T::Type{<:AbstractFloat}=Float32,
)
    config = ControlMatrixReconstructionNodeConfig(
        slope_count,
        reconstructed_count,
        slopes_schema,
        reconstructed_schema,
        control_matrix_schema,
    )
    return algorithm_node(
        name,
        ControlMatrixReconstructionNode{T},
        config;
        props=NamedTuple(),
    )
end

function graph_node_ports(
    ::Type{ControlMatrixReconstructionNode{T}},
    config::ControlMatrixReconstructionNodeConfig,
) where {T}
    return (
        graph_port_contract(
            :slopes,
            :input,
            :data,
            T,
            (config.slope_count,),
            config.slopes_schema,
            :column_major,
        ),
        graph_port_contract(
            :reconstructed,
            :output,
            :data,
            T,
            (config.reconstructed_count,),
            config.reconstructed_schema,
            :column_major,
        ),
        graph_port_contract(
            :control_matrix,
            :input,
            :parameter,
            T,
            (config.reconstructed_count, config.slope_count),
            config.control_matrix_schema,
            :column_major,
        ),
    )
end

struct _ControlMatrixReconstructionOwner{Plan,Slopes,Reconstructed}
    plan::Plan
    slopes::Slopes
    reconstructed::Reconstructed
end

function prepare_graph_node(
    ::Type{ControlMatrixReconstructionNode{T}},
    ::ControlMatrixReconstructionNodeConfig,
    ::NamedTuple{()},
    inputs::NamedTuple{(:slopes,)},
    outputs::NamedTuple{(:reconstructed,)},
    parameters::NamedTuple{(:control_matrix,)},
    target,
) where {T}
    plan = try
        ControlMatrixPlan(parameters.control_matrix)
    catch error
        error isa AdaptiveOpticsSimError || rethrow()
        throw(AlgorithmGraphError(sprint(showerror, error)))
    end
    return _ControlMatrixReconstructionOwner(
        plan,
        inputs.slopes,
        outputs.reconstructed,
    )
end

@inline function step_graph_node!(owner::_ControlMatrixReconstructionOwner)
    reconstruct!(owner.reconstructed, owner.plan, owner.slopes)
    return nothing
end

@inline function reset_graph_node!(owner::_ControlMatrixReconstructionOwner)
    fill!(owner.reconstructed, zero(eltype(owner.reconstructed)))
    return nothing
end

struct CCDDetectorAcquisitionNode{T<:AbstractFloat} end

"""Construction values for one complete-frame single-read CCD acquisition."""
struct CCDDetectorAcquisitionNodeConfig{T<:AbstractFloat,N<:NoiseModel}
    rows::Int
    columns::Int
    binning::Int
    pixel_scale_arcsec::T
    wavelength_m::T
    exposure_duration_s::T
    quantum_efficiency::T
    noise::N
    rng_seed::UInt64
    photon_rate_schema::String
    frame_schema::String
end

@inline function _frame_detector_noise(
    ::Type{T},
    detector_name::AbstractString,
    photon_noise::Bool,
    readout_noise::Bool,
    readout_noise_e::Real,
) where {T<:AbstractFloat}
    sigma = T(readout_noise_e)
    isfinite(sigma) && sigma >= zero(T) || throw(AlgorithmGraphError(
        "$detector_name readout_noise_e must be finite and nonnegative",
    ))
    if readout_noise
        sigma > zero(T) || throw(AlgorithmGraphError(
            "$detector_name readout_noise_e must be positive when readout noise is enabled",
        ))
        return photon_noise ? NoisePhotonReadout{T}(sigma) :
            NoiseReadout{T}(sigma)
    end
    iszero(sigma) || throw(AlgorithmGraphError(
        "$detector_name readout_noise_e must be zero when readout noise is disabled",
    ))
    return photon_noise ? NoisePhoton() : NoiseNone()
end

function _ccd_detector_acquisition_config(
    ::Type{T};
    rows::Integer,
    columns::Integer,
    binning::Integer,
    pixel_scale_arcsec::Real,
    wavelength_m::Real,
    exposure_duration_s::Real,
    quantum_efficiency::Real,
    photon_noise::Bool,
    readout_noise::Bool,
    readout_noise_e::Real,
    rng_seed::Integer,
    photon_rate_schema::AbstractString,
    frame_schema::AbstractString,
) where {T<:AbstractFloat}
    n_rows = Int(rows)
    n_columns = Int(columns)
    bin = Int(binning)
    n_rows > 0 && n_columns > 0 || throw(AlgorithmGraphError(
        "CCD photon-rate dimensions must be positive",
    ))
    bin > 0 || throw(AlgorithmGraphError(
        "CCD binning must be positive",
    ))
    n_rows % bin == 0 && n_columns % bin == 0 || throw(
        AlgorithmGraphError(
            "CCD photon-rate dimensions must be divisible by binning",
        ))

    pixel_scale = T(pixel_scale_arcsec)
    wavelength = T(wavelength_m)
    exposure = T(exposure_duration_s)
    qe = T(quantum_efficiency)
    isfinite(pixel_scale) && pixel_scale > zero(T) || throw(
        AlgorithmGraphError(
            "CCD pixel_scale_arcsec must be finite and positive",
        ))
    isfinite(wavelength) && wavelength > zero(T) || throw(
        AlgorithmGraphError(
            "CCD wavelength_m must be finite and positive",
        ))
    isfinite(exposure) && exposure > zero(T) || throw(
        AlgorithmGraphError(
            "CCD exposure_duration_s must be finite and positive",
        ))
    isfinite(qe) && zero(T) <= qe <= one(T) || throw(
        AlgorithmGraphError(
            "CCD quantum_efficiency must lie in [0, 1]",
        ))
    rng_seed >= 0 || throw(AlgorithmGraphError(
        "CCD rng_seed must be nonnegative",
    ))
    rng_seed <= typemax(UInt64) || throw(AlgorithmGraphError(
        "CCD rng_seed exceeds UInt64",
    ))
    isempty(photon_rate_schema) && throw(AlgorithmGraphError(
        "CCD photon_rate_schema must not be empty",
    ))
    isempty(frame_schema) && throw(AlgorithmGraphError(
        "CCD frame_schema must not be empty",
    ))

    noise = _frame_detector_noise(
        T,
        "CCD",
        photon_noise,
        readout_noise,
        readout_noise_e,
    )
    return CCDDetectorAcquisitionNodeConfig(
        n_rows,
        n_columns,
        bin,
        pixel_scale,
        wavelength,
        exposure,
        qe,
        noise,
        UInt64(rng_seed),
        String(photon_rate_schema),
        String(frame_schema),
    )
end

"""
    ccd_detector_acquisition_node(name; ...)

Declare one complete-frame, single-read CCD acquisition. The node consumes a
cell-integrated detector-plane photon-rate mosaic and writes a detector frame.
Its RNG is explicit persistent state and is restored by `reset_graph!`.
Partial exposure timing, rolling shutters, frame transfer, and readout
readiness are outside the complete-frame graph contract.
"""
function ccd_detector_acquisition_node(
    name::Symbol;
    rows::Integer,
    columns::Integer,
    pixel_scale_arcsec::Real,
    wavelength_m::Real,
    exposure_duration_s::Real,
    quantum_efficiency::Real,
    rng_seed::Integer,
    photon_rate_schema::AbstractString,
    frame_schema::AbstractString,
    binning::Integer=1,
    photon_noise::Bool=true,
    readout_noise::Bool=false,
    readout_noise_e::Real=0,
    T::Type{<:AbstractFloat}=Float32,
)
    config = _ccd_detector_acquisition_config(
        T;
        rows,
        columns,
        binning,
        pixel_scale_arcsec,
        wavelength_m,
        exposure_duration_s,
        quantum_efficiency,
        photon_noise,
        readout_noise,
        readout_noise_e,
        rng_seed,
        photon_rate_schema,
        frame_schema,
    )
    return algorithm_node(
        name,
        CCDDetectorAcquisitionNode{T},
        config;
        props=NamedTuple(),
    )
end

function graph_node_ports(
    ::Type{CCDDetectorAcquisitionNode{T}},
    config::CCDDetectorAcquisitionNodeConfig{T},
) where {T}
    return (
        graph_port_contract(
            :photon_rate,
            :input,
            :data,
            T,
            (config.rows, config.columns),
            config.photon_rate_schema,
            :column_major,
        ),
        graph_port_contract(
            :frame,
            :output,
            :data,
            T,
            (div(config.rows, config.binning),
                div(config.columns, config.binning)),
            config.frame_schema,
            :column_major,
        ),
    )
end

struct _FrameDetectorAcquisitionOwner{P,R,O,A}
    prepared::P
    rng::R
    output::O
    latent_buffer::A
    rng_seed::UInt64
end

@inline _arcsec_to_rad(value::T) where {T<:AbstractFloat} =
    value * T(pi / (180 * 3600))

function prepare_graph_node(
    ::Type{CCDDetectorAcquisitionNode{T}},
    config::CCDDetectorAcquisitionNodeConfig{T},
    ::NamedTuple{()},
    inputs::NamedTuple{(:photon_rate,)},
    outputs::NamedTuple{(:frame,)},
    ::NamedTuple{()},
    target,
) where {T}
    sampling = _arcsec_to_rad(config.pixel_scale_arcsec)
    metadata = OpticalPlaneMetadata(
        DetectorPlane(),
        inputs.photon_rate;
        coordinate_domain=AngularCoordinates(),
        sampling=(sampling, sampling),
        spectral=MonochromaticChannel(config.wavelength_m),
        normalization=PhotonRateNormalization(),
        spatial_measure=CellIntegratedMeasure(),
        coherence=IncoherentIntensityAddition(),
    )
    photon_rate = IntensityMap(metadata, inputs.photon_rate)
    detector = Detector(
        exposure_duration=config.exposure_duration_s,
        qe=config.quantum_efficiency,
        psf_sampling=1,
        binning=config.binning,
        noise=config.noise,
        sensor=CCDSensor(T=T),
        response_model=NullFrameResponse(),
        T=T,
        backend=compute_device_backend(target),
    )
    prepared = prepare_detector_acquisition(detector, photon_rate)
    size(output_frame(detector)) == size(outputs.frame) || throw(
        AlgorithmGraphError(
            "prepared CCD output shape does not match its graph frame port",
        ))
    rng = _prepare_graph_rng(target, config.rng_seed)
    state = detector_acquisition_state(prepared)
    return _FrameDetectorAcquisitionOwner(
        prepared,
        rng,
        outputs.frame,
        state.latent_buffer,
        config.rng_seed,
    )
end

@inline function step_graph_node!(owner::_FrameDetectorAcquisitionOwner)
    frame = capture!(owner.prepared, owner.rng)
    copyto!(owner.output, frame)
    return nothing
end

function reset_graph_node!(owner::_FrameDetectorAcquisitionOwner)
    detector = detector_acquisition_detector(owner.prepared)
    reset_integration!(detector)
    fill!(owner.latent_buffer, zero(eltype(owner.latent_buffer)))
    fill!(output_frame(detector), zero(eltype(output_frame(detector))))
    fill!(owner.output, zero(eltype(owner.output)))
    _reset_graph_rng!(owner.rng, owner.rng_seed)
    return nothing
end

struct CMOSDetectorAcquisitionNode{T<:AbstractFloat} end

@inline function _graph_adc_output_type(bits::Union{Nothing,Int})
    bits === nothing && return nothing
    bits <= 8 && return UInt8
    bits <= 16 && return UInt16
    bits <= 32 && return UInt32
    return UInt64
end

"""Construction values for one complete-frame global-shutter CMOS acquisition."""
struct CMOSDetectorAcquisitionNodeConfig{
    T<:AbstractFloat,
    N<:NoiseModel,
}
    rows::Int
    columns::Int
    binning::Int
    pixel_scale_arcsec::T
    wavelength_m::T
    exposure_duration_s::T
    quantum_efficiency::T
    gain::T
    dark_current_e_per_pixel_s::T
    bits::Union{Nothing,Int}
    full_well_e::Union{Nothing,T}
    noise::N
    column_readout_noise_e::T
    row_readout_noise_e::T
    rng_seed::UInt64
    photon_rate_schema::String
    frame_schema::String
end

function _cmos_detector_acquisition_config(
    ::Type{T};
    rows::Integer,
    columns::Integer,
    binning::Integer,
    pixel_scale_arcsec::Real,
    wavelength_m::Real,
    exposure_duration_s::Real,
    quantum_efficiency::Real,
    gain::Real,
    dark_current_e_per_pixel_s::Real,
    bits::Union{Nothing,Integer},
    full_well_e::Union{Nothing,Real},
    photon_noise::Bool,
    readout_noise::Bool,
    readout_noise_e::Real,
    column_readout_noise_e::Real,
    row_readout_noise_e::Real,
    rng_seed::Integer,
    photon_rate_schema::AbstractString,
    frame_schema::AbstractString,
) where {T<:AbstractFloat}
    n_rows = Int(rows)
    n_columns = Int(columns)
    bin = Int(binning)
    output_bits = isnothing(bits) ? nothing : Int(bits)
    n_rows > 0 && n_columns > 0 || throw(AlgorithmGraphError(
        "CMOS photon-rate dimensions must be positive",
    ))
    bin > 0 || throw(AlgorithmGraphError(
        "CMOS binning must be positive",
    ))
    n_rows % bin == 0 && n_columns % bin == 0 || throw(
        AlgorithmGraphError(
            "CMOS photon-rate dimensions must be divisible by binning",
        ),
    )
    isnothing(output_bits) || 1 <= output_bits <= 64 || throw(
        AlgorithmGraphError("CMOS bits must lie in 1:64"),
    )
    isnothing(output_bits) || !isnothing(full_well_e) || throw(
        AlgorithmGraphError("CMOS bits requires full_well_e"),
    )

    pixel_scale = T(pixel_scale_arcsec)
    wavelength = T(wavelength_m)
    exposure = T(exposure_duration_s)
    qe = T(quantum_efficiency)
    linear_gain = T(gain)
    dark_current = T(dark_current_e_per_pixel_s)
    full_well = isnothing(full_well_e) ? nothing : T(full_well_e)
    column_noise = T(column_readout_noise_e)
    row_noise = T(row_readout_noise_e)

    isfinite(pixel_scale) && pixel_scale > zero(T) || throw(
        AlgorithmGraphError(
            "CMOS pixel_scale_arcsec must be finite and positive",
        ),
    )
    isfinite(wavelength) && wavelength > zero(T) || throw(
        AlgorithmGraphError(
            "CMOS wavelength_m must be finite and positive",
        ),
    )
    isfinite(exposure) && exposure > zero(T) || throw(
        AlgorithmGraphError(
            "CMOS exposure_duration_s must be finite and positive",
        ),
    )
    isfinite(qe) && zero(T) <= qe <= one(T) || throw(
        AlgorithmGraphError(
            "CMOS quantum_efficiency must lie in [0, 1]",
        ),
    )
    isfinite(linear_gain) && linear_gain > zero(T) || throw(
        AlgorithmGraphError("CMOS gain must be finite and positive"),
    )
    isfinite(dark_current) && dark_current >= zero(T) || throw(
        AlgorithmGraphError(
            "CMOS dark_current_e_per_pixel_s must be finite and nonnegative",
        ),
    )
    isnothing(full_well) || isfinite(full_well) && full_well > zero(T) ||
        throw(AlgorithmGraphError(
            "CMOS full_well_e must be finite and positive",
        ))
    isfinite(column_noise) && column_noise >= zero(T) || throw(
        AlgorithmGraphError(
            "CMOS column_readout_noise_e must be finite and nonnegative",
        ),
    )
    isfinite(row_noise) && row_noise >= zero(T) || throw(
        AlgorithmGraphError(
            "CMOS row_readout_noise_e must be finite and nonnegative",
        ),
    )
    rng_seed >= 0 || throw(AlgorithmGraphError(
        "CMOS rng_seed must be nonnegative",
    ))
    rng_seed <= typemax(UInt64) || throw(AlgorithmGraphError(
        "CMOS rng_seed exceeds UInt64",
    ))
    isempty(photon_rate_schema) && throw(AlgorithmGraphError(
        "CMOS photon_rate_schema must not be empty",
    ))
    isempty(frame_schema) && throw(AlgorithmGraphError(
        "CMOS frame_schema must not be empty",
    ))

    noise = _frame_detector_noise(
        T,
        "CMOS",
        photon_noise,
        readout_noise,
        readout_noise_e,
    )
    return CMOSDetectorAcquisitionNodeConfig(
        n_rows,
        n_columns,
        bin,
        pixel_scale,
        wavelength,
        exposure,
        qe,
        linear_gain,
        dark_current,
        output_bits,
        full_well,
        noise,
        column_noise,
        row_noise,
        UInt64(rng_seed),
        String(photon_rate_schema),
        String(frame_schema),
    )
end

"""
    cmos_detector_acquisition_node(name; ...)

Declare one complete-frame global-shutter CMOS acquisition. The node consumes
a cell-integrated detector-plane photon-rate mosaic and writes one completed
frame. Independent, shared-column, and shared-row read-noise components are
configured separately. `gain` is a linear post-readout gain, not decibels.
When `bits` is present, the internal detector readout uses the smallest
unsigned integer type that contains the configured ADC range; the graph port
remains `T` and therefore carries integer-valued samples in its floating-point
buffer.
Its RNG is explicit persistent state and is restored by `reset_graph!`.
"""
function cmos_detector_acquisition_node(
    name::Symbol;
    rows::Integer,
    columns::Integer,
    pixel_scale_arcsec::Real,
    wavelength_m::Real,
    exposure_duration_s::Real,
    quantum_efficiency::Real,
    rng_seed::Integer,
    photon_rate_schema::AbstractString,
    frame_schema::AbstractString,
    binning::Integer=1,
    gain::Real=1,
    dark_current_e_per_pixel_s::Real=0,
    bits::Union{Nothing,Integer}=nothing,
    full_well_e::Union{Nothing,Real}=nothing,
    photon_noise::Bool=true,
    readout_noise::Bool=false,
    readout_noise_e::Real=0,
    column_readout_noise_e::Real=0,
    row_readout_noise_e::Real=0,
    T::Type{<:AbstractFloat}=Float32,
)
    config = _cmos_detector_acquisition_config(
        T;
        rows,
        columns,
        binning,
        pixel_scale_arcsec,
        wavelength_m,
        exposure_duration_s,
        quantum_efficiency,
        gain,
        dark_current_e_per_pixel_s,
        bits,
        full_well_e,
        photon_noise,
        readout_noise,
        readout_noise_e,
        column_readout_noise_e,
        row_readout_noise_e,
        rng_seed,
        photon_rate_schema,
        frame_schema,
    )
    return algorithm_node(
        name,
        CMOSDetectorAcquisitionNode{T},
        config;
        props=NamedTuple(),
    )
end

function graph_node_ports(
    ::Type{CMOSDetectorAcquisitionNode{T}},
    config::CMOSDetectorAcquisitionNodeConfig{T},
) where {T}
    return (
        graph_port_contract(
            :photon_rate,
            :input,
            :data,
            T,
            (config.rows, config.columns),
            config.photon_rate_schema,
            :column_major,
        ),
        graph_port_contract(
            :frame,
            :output,
            :data,
            T,
            (div(config.rows, config.binning),
                div(config.columns, config.binning)),
            config.frame_schema,
            :column_major,
        ),
    )
end

function prepare_graph_node(
    ::Type{CMOSDetectorAcquisitionNode{T}},
    config::CMOSDetectorAcquisitionNodeConfig{T},
    ::NamedTuple{()},
    inputs::NamedTuple{(:photon_rate,)},
    outputs::NamedTuple{(:frame,)},
    ::NamedTuple{()},
    target,
) where {T}
    sampling = _arcsec_to_rad(config.pixel_scale_arcsec)
    metadata = OpticalPlaneMetadata(
        DetectorPlane(),
        inputs.photon_rate;
        coordinate_domain=AngularCoordinates(),
        sampling=(sampling, sampling),
        spectral=MonochromaticChannel(config.wavelength_m),
        normalization=PhotonRateNormalization(),
        spatial_measure=CellIntegratedMeasure(),
        coherence=IncoherentIntensityAddition(),
    )
    photon_rate = IntensityMap(metadata, inputs.photon_rate)
    backend = compute_device_backend(target)
    sensor = CMOSSensor(
        column_readout_sigma=config.column_readout_noise_e,
        row_readout_sigma=config.row_readout_noise_e,
        timing_model=GlobalShutter(),
        T=T,
        backend=backend,
    )
    detector = Detector(
        exposure_duration=config.exposure_duration_s,
        qe=config.quantum_efficiency,
        psf_sampling=1,
        binning=config.binning,
        gain=config.gain,
        dark_current=config.dark_current_e_per_pixel_s,
        bits=config.bits,
        full_well=config.full_well_e,
        output_type=_graph_adc_output_type(config.bits),
        noise=config.noise,
        sensor=sensor,
        response_model=NullFrameResponse(),
        T=T,
        backend=backend,
    )
    prepared = prepare_detector_acquisition(detector, photon_rate)
    size(output_frame(detector)) == size(outputs.frame) || throw(
        AlgorithmGraphError(
            "prepared CMOS output shape does not match its graph frame port",
        ),
    )
    rng = _prepare_graph_rng(target, config.rng_seed)
    state = detector_acquisition_state(prepared)
    return _FrameDetectorAcquisitionOwner(
        prepared,
        rng,
        outputs.frame,
        state.latent_buffer,
        config.rng_seed,
    )
end

struct EMCCDDetectorAcquisitionNode{T<:AbstractFloat} end

"""Construction values for one complete-frame linear-mode EMCCD acquisition."""
struct EMCCDDetectorAcquisitionNodeConfig{
    T<:AbstractFloat,
    N<:NoiseModel,
}
    rows::Int
    columns::Int
    binning::Int
    normalized_pupil_sampling::T
    wavelength_m::T
    exposure_duration_s::T
    quantum_efficiency::T
    gain::T
    dark_current_e_per_pixel_s::T
    bits::Union{Nothing,Int}
    full_well_e::Union{Nothing,T}
    noise::N
    excess_noise_factor::T
    clock_induced_charge_e_per_pixel_frame::T
    register_full_well_e::Union{Nothing,T}
    em_gain_range::Tuple{T,T}
    rng_seed::UInt64
    photon_rate_schema::String
    frame_schema::String
end

function _emccd_detector_acquisition_config(
    ::Type{T};
    rows::Integer,
    columns::Integer,
    binning::Integer,
    normalized_pupil_sampling::Real,
    wavelength_m::Real,
    exposure_duration_s::Real,
    quantum_efficiency::Real,
    gain::Real,
    dark_current_e_per_pixel_s::Real,
    bits::Union{Nothing,Integer},
    full_well_e::Union{Nothing,Real},
    photon_noise::Bool,
    readout_noise::Bool,
    readout_noise_e::Real,
    excess_noise_factor::Real,
    clock_induced_charge_e_per_pixel_frame::Real,
    register_full_well_e::Union{Nothing,Real},
    em_gain_min::Real,
    em_gain_max::Real,
    rng_seed::Integer,
    photon_rate_schema::AbstractString,
    frame_schema::AbstractString,
) where {T<:AbstractFloat}
    n_rows = Int(rows)
    n_columns = Int(columns)
    bin = Int(binning)
    output_bits = isnothing(bits) ? nothing : Int(bits)
    n_rows > 0 && n_columns > 0 || throw(AlgorithmGraphError(
        "EMCCD photon-rate dimensions must be positive",
    ))
    bin > 0 || throw(AlgorithmGraphError(
        "EMCCD binning must be positive",
    ))
    n_rows % bin == 0 && n_columns % bin == 0 || throw(
        AlgorithmGraphError(
            "EMCCD photon-rate dimensions must be divisible by binning",
        ),
    )
    isnothing(output_bits) || 1 <= output_bits <= 64 || throw(
        AlgorithmGraphError("EMCCD bits must lie in 1:64"),
    )
    isnothing(output_bits) || !isnothing(full_well_e) || throw(
        AlgorithmGraphError("EMCCD bits requires full_well_e"),
    )

    sampling = T(normalized_pupil_sampling)
    wavelength = T(wavelength_m)
    exposure = T(exposure_duration_s)
    qe = T(quantum_efficiency)
    em_gain = T(gain)
    dark_current = T(dark_current_e_per_pixel_s)
    full_well = isnothing(full_well_e) ? nothing : T(full_well_e)
    excess_noise = T(excess_noise_factor)
    cic = T(clock_induced_charge_e_per_pixel_frame)
    register_full_well = isnothing(register_full_well_e) ? nothing :
        T(register_full_well_e)
    gain_range = (T(em_gain_min), T(em_gain_max))

    isfinite(sampling) && sampling > zero(T) || throw(AlgorithmGraphError(
        "EMCCD normalized_pupil_sampling must be finite and positive",
    ))
    isfinite(wavelength) && wavelength > zero(T) || throw(
        AlgorithmGraphError(
            "EMCCD wavelength_m must be finite and positive",
        ),
    )
    isfinite(exposure) && exposure > zero(T) || throw(AlgorithmGraphError(
        "EMCCD exposure_duration_s must be finite and positive",
    ))
    isfinite(qe) && zero(T) <= qe <= one(T) || throw(AlgorithmGraphError(
        "EMCCD quantum_efficiency must lie in [0, 1]",
    ))
    isfinite(em_gain) && em_gain > zero(T) || throw(AlgorithmGraphError(
        "EMCCD gain must be finite and positive",
    ))
    isfinite(dark_current) && dark_current >= zero(T) || throw(
        AlgorithmGraphError(
            "EMCCD dark_current_e_per_pixel_s must be finite and nonnegative",
        ),
    )
    isnothing(full_well) || isfinite(full_well) && full_well > zero(T) ||
        throw(AlgorithmGraphError(
            "EMCCD full_well_e must be finite and positive",
        ))
    isfinite(excess_noise) && excess_noise >= one(T) || throw(
        AlgorithmGraphError(
            "EMCCD excess_noise_factor must be finite and at least one",
        ),
    )
    isfinite(cic) && cic >= zero(T) || throw(AlgorithmGraphError(
        "EMCCD clock-induced charge must be finite and nonnegative",
    ))
    isnothing(register_full_well) ||
        isfinite(register_full_well) && register_full_well > zero(T) ||
        throw(AlgorithmGraphError(
            "EMCCD register_full_well_e must be finite and positive",
        ))
    isfinite(gain_range[1]) && gain_range[1] > zero(T) || throw(
        AlgorithmGraphError(
            "EMCCD em_gain_min must be finite and positive",
        ),
    )
    !isnan(gain_range[2]) && gain_range[2] >= gain_range[1] || throw(
        AlgorithmGraphError(
            "EMCCD em_gain_max must be at least em_gain_min",
        ),
    )
    gain_range[1] <= em_gain <= gain_range[2] || throw(
        AlgorithmGraphError("EMCCD gain must lie in its EM gain range"),
    )
    rng_seed >= 0 || throw(AlgorithmGraphError(
        "EMCCD rng_seed must be nonnegative",
    ))
    rng_seed <= typemax(UInt64) || throw(AlgorithmGraphError(
        "EMCCD rng_seed exceeds UInt64",
    ))
    isempty(photon_rate_schema) && throw(AlgorithmGraphError(
        "EMCCD photon_rate_schema must not be empty",
    ))
    isempty(frame_schema) && throw(AlgorithmGraphError(
        "EMCCD frame_schema must not be empty",
    ))

    noise = _frame_detector_noise(
        T,
        "EMCCD",
        photon_noise,
        readout_noise,
        readout_noise_e,
    )
    return EMCCDDetectorAcquisitionNodeConfig(
        n_rows,
        n_columns,
        bin,
        sampling,
        wavelength,
        exposure,
        qe,
        em_gain,
        dark_current,
        output_bits,
        full_well,
        noise,
        excess_noise,
        cic,
        register_full_well,
        gain_range,
        UInt64(rng_seed),
        String(photon_rate_schema),
        String(frame_schema),
    )
end

"""
    emccd_detector_acquisition_node(name; ...)

Declare one complete-frame linear-mode EMCCD acquisition. The node consumes a
cell-integrated four-pupil photon-rate frame and writes a detector frame. Its
RNG is explicit persistent state and is restored by `reset_graph!`.

Frame-transfer timing and readout readiness are outside the complete-frame
graph contract; this node performs one complete sequential acquisition per
graph step.
"""
function emccd_detector_acquisition_node(
    name::Symbol;
    rows::Integer,
    columns::Integer,
    normalized_pupil_sampling::Real,
    wavelength_m::Real,
    exposure_duration_s::Real,
    quantum_efficiency::Real,
    rng_seed::Integer,
    photon_rate_schema::AbstractString,
    frame_schema::AbstractString,
    binning::Integer=1,
    gain::Real=1,
    dark_current_e_per_pixel_s::Real=0,
    bits::Union{Nothing,Integer}=nothing,
    full_well_e::Union{Nothing,Real}=nothing,
    photon_noise::Bool=true,
    readout_noise::Bool=false,
    readout_noise_e::Real=0,
    excess_noise_factor::Real=1,
    clock_induced_charge_e_per_pixel_frame::Real=0,
    register_full_well_e::Union{Nothing,Real}=nothing,
    em_gain_min::Real=1,
    em_gain_max::Real=Inf,
    T::Type{<:AbstractFloat}=Float32,
)
    config = _emccd_detector_acquisition_config(
        T;
        rows,
        columns,
        binning,
        normalized_pupil_sampling,
        wavelength_m,
        exposure_duration_s,
        quantum_efficiency,
        gain,
        dark_current_e_per_pixel_s,
        bits,
        full_well_e,
        photon_noise,
        readout_noise,
        readout_noise_e,
        excess_noise_factor,
        clock_induced_charge_e_per_pixel_frame,
        register_full_well_e,
        em_gain_min,
        em_gain_max,
        rng_seed,
        photon_rate_schema,
        frame_schema,
    )
    return algorithm_node(
        name,
        EMCCDDetectorAcquisitionNode{T},
        config;
        props=NamedTuple(),
    )
end

function graph_node_ports(
    ::Type{EMCCDDetectorAcquisitionNode{T}},
    config::EMCCDDetectorAcquisitionNodeConfig{T},
) where {T}
    return (
        graph_port_contract(
            :photon_rate,
            :input,
            :data,
            T,
            (config.rows, config.columns),
            config.photon_rate_schema,
            :column_major,
        ),
        graph_port_contract(
            :frame,
            :output,
            :data,
            T,
            (div(config.rows, config.binning),
                div(config.columns, config.binning)),
            config.frame_schema,
            :column_major,
        ),
    )
end

function prepare_graph_node(
    ::Type{EMCCDDetectorAcquisitionNode{T}},
    config::EMCCDDetectorAcquisitionNodeConfig{T},
    ::NamedTuple{()},
    inputs::NamedTuple{(:photon_rate,)},
    outputs::NamedTuple{(:frame,)},
    ::NamedTuple{()},
    target,
) where {T}
    metadata = OpticalPlaneMetadata(
        DetectorPlane(),
        inputs.photon_rate;
        coordinate_domain=NormalizedPupilCoordinates(),
        sampling=(config.normalized_pupil_sampling,
            config.normalized_pupil_sampling),
        spectral=MonochromaticChannel(config.wavelength_m),
        normalization=PhotonRateNormalization(),
        spatial_measure=CellIntegratedMeasure(),
        coherence=IncoherentIntensityAddition(),
    )
    photon_rate = IntensityMap(metadata, inputs.photon_rate)
    sensor = EMCCDSensor(
        excess_noise_factor=config.excess_noise_factor,
        clock_induced_charge_per_frame=
            config.clock_induced_charge_e_per_pixel_frame,
        register_full_well=config.register_full_well_e,
        em_gain_range=config.em_gain_range,
        T=T,
    )
    detector = Detector(
        exposure_duration=config.exposure_duration_s,
        qe=config.quantum_efficiency,
        psf_sampling=1,
        binning=config.binning,
        noise=config.noise,
        gain=config.gain,
        dark_current=config.dark_current_e_per_pixel_s,
        bits=config.bits,
        full_well=config.full_well_e,
        sensor=sensor,
        response_model=NullFrameResponse(),
        T=T,
        backend=compute_device_backend(target),
    )
    prepared = prepare_detector_acquisition(detector, photon_rate)
    size(output_frame(detector)) == size(outputs.frame) || throw(
        AlgorithmGraphError(
            "prepared EMCCD output shape does not match its graph frame port",
        ),
    )
    rng = _prepare_graph_rng(target, config.rng_seed)
    state = detector_acquisition_state(prepared)
    return _FrameDetectorAcquisitionOwner(
        prepared,
        rng,
        outputs.frame,
        state.latent_buffer,
        config.rng_seed,
    )
end
