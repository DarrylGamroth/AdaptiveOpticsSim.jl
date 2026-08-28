module CalculonControlTests

using AdaptiveOpticsSim
using CalculonAlgorithms
using Test

using AdaptiveOpticsSim.AlgorithmGraphs
using AdaptiveOpticsSim.Control

const CalculonExtension = Base.get_extension(
    AdaptiveOpticsSim,
    :AdaptiveOpticsSimCalculonAlgorithmsExt,
)

@testset "AdaptiveOpticsSim Calculon controller declaration" begin
    @test !isnothing(CalculonExtension)
    declaration = CalculonExtension.DiscreteIntegratorControllerF32
    configuration = CalculonExtension.DiscreteIntegratorControllerF32Configuration(
        extent=4,
        gain=0.5f0,
        tau_s=0.2f0,
        sample_period_s=0.1f0,
        input_schema=CONTROLLER_RESIDUAL_ERROR_V1,
        output_schema=CONTROLLER_COMMAND_V1,
    )

    direct = DiscreteIntegratorController(4; gain=0.5f0, tau=0.2f0, T=Float32)
    declared = prepare_algorithm(declaration, configuration)
    @test algorithm_label(declaration) == "discrete-integrator-controller-f32"
    @test algorithm_property_policy(declaration) isa RuntimeDeclaredProperties
    @test property_values(declared.plan) == (gain=0.5f0, tau_s=0.2f0)
    @test declared.workspace.state isa Control.DiscreteIntegratorState
    @test declared.workspace.workspace isa Control.DiscreteIntegratorWorkspace
    @test !Base.mightalias(
        declared.workspace.state.integrated_command,
        declared.workspace.workspace.next_integrated_command,
    )
    @test !Base.mightalias(
        declared.workspace.state.command,
        declared.workspace.workspace.next_command,
    )

    input = ones(Float32, 4)
    direct_output = similar(input)
    declared_output = similar(input)
    for _ in 1:2
        copyto!(direct_output, Control.update!(direct, input, 0.1f0))
        @test CalculonAlgorithms.process!(declared_output, declared, input) ===
            declared_output
        @test declared_output == direct_output
        @test declared.workspace.state.integrated_command ==
            Control.discrete_integrator_state(direct).integrated_command
        @test declared.workspace.state.command == Control.controller_output(direct)
    end
    @test @allocated(CalculonAlgorithms.process!(declared_output, declared, input)) == 0
    @test @inferred(CalculonAlgorithms.process!(declared_output, declared, input)) ===
        declared_output

    schema = property_schema(declared.plan)
    gain = only(spec for spec in schema if spec.info.name === :gain)
    tau_s = only(spec for spec in schema if spec.info.name === :tau_s)
    replacement = CalculonAlgorithms.prepare_properties(
        declared.plan,
        (
            PropertyAssignment(gain.info.id, 0.25f0),
            PropertyAssignment(tau_s.info.id, 0.4f0),
        ),
    )
    @test property_values(replacement) == (gain=0.25f0, tau_s=0.4f0)
    @test replacement.sample_period_s == declared.plan.sample_period_s

    @test CalculonAlgorithms.reset!(declared) === declared
    @test all(iszero, declared.workspace.state.integrated_command)
    @test all(iszero, declared.workspace.state.command)
    @test all(iszero, declared.workspace.workspace.next_integrated_command)
    @test all(iszero, declared.workspace.workspace.next_command)
end

@testset "AdaptiveOpticsSim Calculon controller in portable graph" begin
    declaration = CalculonExtension.DiscreteIntegratorControllerF32
    configuration = CalculonExtension.DiscreteIntegratorControllerF32Configuration(
        extent=3,
        gain=0.5f0,
        tau_s=0.2f0,
        sample_period_s=0.1f0,
        input_schema=CONTROLLER_RESIDUAL_ERROR_V1,
        output_schema=CONTROLLER_COMMAND_V1,
    )
    input = ones(Float32, 3)
    graph = prepare_algorithm_graph(algorithm_graph(
        (algorithm_node(:controller, declaration, configuration),);
        inputs=(graph_input(:residual, :controller => :input, input),),
        outputs=(graph_output(:command, :controller => :output),),
    ))

    @test step_graph!(graph) === graph
    @test graph_output(graph, Val(:command)) == fill(0.025f0, 3)
    step_graph!(graph)
    @test graph_output(graph, Val(:command)) == fill(0.0625f0, 3)
    @test @allocated(step_graph!(graph)) == 0
    @test @inferred(step_graph!(graph)) === graph

    prepared = AlgorithmGraphs.prepared_algorithm(graph, Val(:controller))
    @test reset_graph!(graph) === graph
    @test all(iszero, prepared.workspace.state.command)
    @test all(iszero, prepared.workspace.workspace.next_command)
end

@testset "AdaptiveOpticsSim Calculon controller rejects invalid configuration" begin
    declaration = CalculonExtension.DiscreteIntegratorControllerF32
    configuration(; kwargs...) =
        CalculonExtension.DiscreteIntegratorControllerF32Configuration(;
            extent=2,
            gain=0.5f0,
            tau_s=0.2f0,
            sample_period_s=0.1f0,
            input_schema=CONTROLLER_RESIDUAL_ERROR_V1,
            output_schema=CONTROLLER_COMMAND_V1,
            kwargs...,
        )

    @test_throws InvalidConfiguration prepare_algorithm(
        declaration,
        configuration(extent=0),
    )
    @test_throws InvalidConfiguration prepare_algorithm(
        declaration,
        configuration(tau_s=0.0f0),
    )
    @test_throws ArgumentError prepare_algorithm(
        declaration,
        configuration(sample_period_s=0.0f0),
    )
    @test_throws ArgumentError prepare_algorithm(
        declaration,
        configuration(input_schema="unversioned"),
    )
end

end # module CalculonControlTests
