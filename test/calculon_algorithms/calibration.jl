module CalculonCalibrationTests

using AdaptiveOpticsSim
using CalculonAlgorithms
using Test

using AdaptiveOpticsSim.AlgorithmGraphs
using AdaptiveOpticsSim.Calibration

const CalculonExtension = Base.get_extension(
    AdaptiveOpticsSim,
    :AdaptiveOpticsSimCalculonAlgorithmsExt,
)

const COEFFICIENTS_SCHEMA = "org.adaptiveopticssim.modal-opd-coefficients/1"
const OPD_SCHEMA = "org.adaptiveopticssim.pupil-opd/1"
const BASIS_SCHEMA = "org.adaptiveopticssim.modal-opd-basis/1"
const PUPIL_SUPPORT_SCHEMA = "org.adaptiveopticssim.pupil-support/1"

function configuration(; kwargs...)
    return CalculonExtension.ModalOPDExpansionF32Configuration(;
        pupil_rows=2,
        pupil_columns=3,
        mode_count=2,
        coefficients_schema=COEFFICIENTS_SCHEMA,
        opd_schema=OPD_SCHEMA,
        basis_schema=BASIS_SCHEMA,
        pupil_support_schema=PUPIL_SUPPORT_SCHEMA,
        kwargs...,
    )
end

function fixture()
    basis = reshape(Float32.(1:12), 2, 3, 2)
    pupil_support = Bool[true false true; true true false]
    coefficients = Float32[0.25, -0.5]
    return basis, pupil_support, coefficients
end

@testset "AdaptiveOpticsSim modal OPD expansion plan" begin
    basis, pupil_support, coefficients = fixture()
    plan = ModalOPDExpansionPlan(basis, pupil_support)
    basis .= 0
    pupil_support .= false

    output = zeros(Float32, 2, 3)
    @test combine_basis!(output, plan, coefficients) === output
    @test output == Float32[-3.25 0.0 -4.25; -3.5 -4.0 0.0]
    @test any(!iszero, plan.basis)
    @test any(plan.pupil_support)
    @test @allocated(combine_basis!(output, plan, coefficients)) == 0
    @test @inferred(combine_basis!(output, plan, coefficients)) === output

    basis_parent = reshape(Float32.(1:36), 3, 4, 3)
    pupil_parent = trues(3, 4)
    coefficient_parent = Float32[9, 0.25, -0.5, 9]
    output_parent = zeros(Float32, 4, 5)
    view_plan = ModalOPDExpansionPlan(
        @view(basis_parent[1:2, 1:3, 1:2]),
        @view(pupil_parent[1:2, 1:3]),
    )
    output_view = @view(output_parent[2:3, 2:4])
    coefficient_view = @view(coefficient_parent[2:3])
    @test combine_basis!(output_view, view_plan, coefficient_view) === output_view
    @test output_view == Calibration.combine_basis(
        view_plan.basis,
        collect(coefficient_view),
        view_plan.pupil_support,
    )

    @test_throws DimensionMismatchError combine_basis!(
        zeros(Float32, 3, 2),
        plan,
        coefficients,
    )
    @test_throws DimensionMismatchError combine_basis!(
        output,
        plan,
        ones(Float32, 3),
    )
end

@testset "AdaptiveOpticsSim modal OPD Calculon declaration" begin
    declaration = CalculonExtension.ModalOPDExpansionF32
    prepared = prepare_algorithm(declaration, configuration())
    basis, pupil_support, coefficients = fixture()
    @test algorithm_label(declaration) == "modal-opd-expansion-f32"
    @test algorithm_property_policy(declaration) isa
        CalculonAlgorithms.NoDeclaredProperties
    @test replace_parameter!(prepared, :basis, basis) === prepared
    @test replace_parameter!(prepared, :pupil_support, pupil_support) === prepared
    basis .= 0
    pupil_support .= false

    output = zeros(Float32, 2, 3)
    @test CalculonAlgorithms.process!(output, prepared, coefficients) === output
    @test output == Float32[-3.25 0.0 -4.25; -3.5 -4.0 0.0]
    @test any(!iszero, prepared.plan.basis)
    @test any(prepared.plan.pupil_support)
    @test @allocated(
        CalculonAlgorithms.process!(output, prepared, coefficients)
    ) == 0
    @test @inferred(
        CalculonAlgorithms.process!(output, prepared, coefficients)
    ) === output

    formats = algorithm_formats(prepared)
    @test map(format -> format.shape, formats) == (
        (2,),
        (2, 3),
        (2, 3, 2),
        (2, 3),
    )
    @test map(format -> format.element_type, formats) == (
        Float32,
        Float32,
        Float32,
        Bool,
    )
    @test_throws DimensionMismatch replace_parameter!(
        prepared,
        :basis,
        zeros(Float32, 2, 3, 3),
    )
end

@testset "AdaptiveOpticsSim modal OPD expansion in portable graph" begin
    declaration = CalculonExtension.ModalOPDExpansionF32
    basis, pupil_support, coefficients = fixture()
    graph = prepare_algorithm_graph(algorithm_graph(
        (algorithm_node(:expansion, declaration, configuration()),);
        inputs=(graph_input(:coefficients, :expansion => :coefficients, coefficients),),
        outputs=(graph_output(:opd, :expansion => :opd),),
        parameters=(
            sparse_parameter(:expansion => :basis, basis),
            sparse_parameter(:expansion => :pupil_support, pupil_support),
        ),
    ))

    @test step_graph!(graph) === graph
    @test graph_output(graph, Val(:opd)) ==
        Float32[-3.25 0.0 -4.25; -3.5 -4.0 0.0]
    @test @allocated(step_graph!(graph)) == 0
    @test @inferred(step_graph!(graph)) === graph
    @test reset_graph!(graph) === graph
end

@testset "AdaptiveOpticsSim modal OPD declaration rejects invalid configuration" begin
    declaration = CalculonExtension.ModalOPDExpansionF32
    @test_throws ArgumentError prepare_algorithm(
        declaration,
        configuration(pupil_rows=0),
    )
    @test_throws ArgumentError prepare_algorithm(
        declaration,
        configuration(mode_count=0),
    )
    @test_throws ArgumentError prepare_algorithm(
        declaration,
        configuration(opd_schema="unversioned"),
    )
end

end # module CalculonCalibrationTests
