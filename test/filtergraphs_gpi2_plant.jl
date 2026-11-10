# Run with AdaptiveOpticsSim, Proper, FilterGraphAlgorithms, and FilterGraphs.
using AdaptiveOpticsSim, Proper, FilterGraphAlgorithms, FilterGraphs, Test
const A = AdaptiveOpticsSim.AlgorithmGraphs
const O = AdaptiveOpticsSim.Optics

@testset "GPI2 modal calibration admission" begin
    node = A.modal_opd_expansion_node(:tiptilt; pupil_rows=2, pupil_columns=3, mode_count=2,
        coefficients_schema="test.tiptilt/1", opd_schema="test.opd/1",
        basis_schema="test.basis/1", pupil_support_schema="test.pupil-support/1")
    basis = cat(fill(1f0, 2, 3), fill(2f0, 2, 3); dims=3)
    support = fill(true, 2, 3)
    function definition(basis, support)
        return A.algorithm_graph((node,);
            inputs=(A.graph_input(:tiptilt, :tiptilt=>:coefficients, zeros(Float32, 2)),),
            outputs=(A.graph_output(:opd, :tiptilt=>:opd),),
            parameters=(A.sparse_parameter(:tiptilt=>:pupil_support, support),
                A.sparse_parameter(:tiptilt=>:basis, basis)))
    end
    recipe = O.gpi2_plant_filter_graph(definition(basis, support))
    parameters = only(recipe.document.nodes).config.parameters
    @test keys(parameters) == (:basis, :pupil_support)
    @test parameters.basis == basis && parameters.basis !== basis
    @test_throws A.AlgorithmGraphError O.gpi2_plant_filter_graph(definition(Float64.(basis), support))
    @test_throws A.AlgorithmGraphError O.gpi2_plant_filter_graph(definition(basis, fill(true, 3, 2)))
    graph = FilterGraphs.prepare_graph(recipe.document; algorithms=recipe.algorithms)
    opd = zeros(Float32, 2, 3)
    @test FilterGraphs.process!((; opd), graph, (; tiptilt=Float32[3, 4])) == (; opd=true)
    @test opd == fill(11f0, 2, 3)
end
