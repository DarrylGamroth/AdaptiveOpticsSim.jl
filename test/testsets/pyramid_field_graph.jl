# Complex pupil fields carry relative amplitude here; source radiometry is applied
# once by the graph adapter before the existing photon-rate ElectricField seam.
using Test
using YAML
using AdaptiveOpticsSim
using AdaptiveOpticsSim.AlgorithmGraphs
using AdaptiveOpticsSim.AlgorithmGraphs: prepared_graph_node, graph_node_capture_capability
using AdaptiveOpticsSim.Optics
using AdaptiveOpticsSim.WavefrontSensors
using AdaptiveOpticsSim.Backends

function pyramid_field_graph(input, quantity; kwargs...)
    options = (resolution=size(input, 1), telescope_diameter_m=2.0,
        pupil_samples=size(input, 1), source_wavelength_m=750e-9,
        source_photon_irradiance_m2_s=9.0,
        modulation=2, modulation_points=4, diffraction_padding=2,
        opd_schema="test.pyramid.opd/1", field_schema="test.pyramid.field/1",
        photon_rate_schema="test.pyramid.rate/1")
    node = pyramid_rate_node(:pyramid; options..., input_quantity=quantity, kwargs...)
    definition = algorithm_graph((node,);
        inputs=(graph_input(:input, :pyramid => quantity, input),),
        outputs=(graph_output(:rate, :pyramid => :photon_rate),))
    return prepare_algorithm_graph(definition)
end

function pyramid_field_graph_allocation(graph)
    step_graph!(graph)
    return @allocated step_graph!(graph)
end

@testset "Pyramid graph complex pupil input" begin
    for T in (Float32, Float64)
        n = 8
        wavelength_m = T(750e-9)
        opd = T[20e-9 * (i - j) + 3e-9 * (i - 4)^2 for i in 1:n, j in 1:n]
        reference = pyramid_field_graph(opd, :opd; T=T)
        reference_owner = prepared_graph_node(reference, Val(:pyramid))
        pupil = reference_owner.prepared.input
        field = @. pupil.amplitude * cispi(T(2) * opd / wavelength_m)
        field_before = copy(field)
        graph = pyramid_field_graph(field, :electric_field; T=T)
        owner = prepared_graph_node(graph, Val(:pyramid))
        output = graph_output(graph, Val(:rate))
        @test owner.field_input === field
        @test owner.prepared.input.values !== field
        @test owner.prepared.output.values === output
        @test owner.prepared.front_end.source === nothing
        @test size(owner.prepared.input.values) == (n, n)
        @test size(owner.prepared.workspace.field) == (2n, 2n)
        @test owner.prepared.input.metadata.sampling == (T(2) / n, T(2) / n)
        @test owner.prepared.input.metadata.spectral.wavelength_m == wavelength_m
        @test graph_node_capture_capability(owner) isa AlgorithmGraphs.GraphNodeCaptureUnsupported
        @test @inferred(step_graph!(graph)) === graph
        step_graph!(reference)
        @test output ≈ graph_output(reference, Val(:rate)) rtol=100eps(T)
        @test field == field_before
        scale = T(3) * T(2) / n
        @test owner.prepared.input.values ≈ scale .* field rtol=eps(T)
        @test sum(output) ≈ sum(abs2, field) * scale^2 rtol=100eps(T)
        saved = copy(output)
        step_graph!(graph)
        @test output ≈ saved rtol=eps(T)
        @test pyramid_field_graph_allocation(graph) == 0

        # Piston is retained in the field, and cancels only in measured intensity.
        field .*= cis(T(0.73))
        step_graph!(graph)
        @test output ≈ saved rtol=100eps(T)
        field .= field_before
        conjugate_graph = pyramid_field_graph(conj.(field), :electric_field; T=T)
        step_graph!(conjugate_graph)
        @test !isapprox(graph_output(conjugate_graph, Val(:rate)), saved; rtol=T(1e-3))

        # Nonuniform illumination, including outside the nominal aperture, must
        # survive: applying the telescope mask again would fail this energy test.
        for j in 1:n, i in 1:n
            field[i, j] = T(0.3 + i / n) * cis(T(0.4 * i - 0.2 * j))
        end
        step_graph!(graph)
        @test owner.prepared.input.values ≈ scale .* field rtol=eps(T)
        @test sum(output) ≈ sum(abs2, field) * scale^2 rtol=100eps(T)
        @test pyramid_field_graph_allocation(graph) == 0
        @test_throws WavefrontSensors.WFSPreparationError prepare_wfs_optics(
            reference_owner.prepared.front_end,
            owner.prepared.input, owner.prepared.output)
        metadata = owner.prepared.input.metadata
        wrong_wavelength = OpticalPlaneMetadata(PupilPlane(), owner.prepared.input.values;
            coordinate_domain=MetricCoordinates(), sampling=metadata.sampling,
            spectral=MonochromaticChannel(T(500e-9)),
            normalization=PhotonRateNormalization(), spatial_measure=CellIntegratedMeasure(),
            coherence=CoherentFieldCombination())
        @test_throws WavefrontSensors.WFSPreparationError prepare_wfs_optics(
            owner.prepared.front_end,
            ElectricField(wrong_wavelength, owner.prepared.input.values), owner.prepared.output)
        replacement_field = ElectricField(metadata, copy(owner.prepared.input.values))
        @test_throws WavefrontSensors.WFSPreparationError form_wfs_optical_products!(
            owner.prepared.output, replacement_field, owner.prepared)
    end
end

@testset "Pyramid field declaration and file checks" begin
    options = (resolution=8, telescope_diameter_m=2.0, pupil_samples=8,
        source_wavelength_m=750e-9, source_photon_irradiance_m2_s=9.0,
        photon_rate_schema="test.pyramid.rate/1")
    @test_throws AlgorithmGraphError pyramid_rate_node(:pyramid; options...,
        input_quantity=:electric_field)
    @test_throws AlgorithmGraphError pyramid_rate_node(:pyramid; options...,
        input_quantity=:electric_field, field_schema="test.field/1",
        source_separation_arcsec=0.1)
    @test_throws AlgorithmGraphError pyramid_rate_node(:pyramid; options...,
        input_quantity=:other, field_schema="test.field/1")
    @test_throws AlgorithmGraphError pyramid_rate_node(:pyramid; options...)
    @test_throws AlgorithmGraphError pyramid_field_graph(zeros(ComplexF32, 8, 7), :electric_field)

    text = """
    schema_version: 1
    name: pyramid_field
    nodes:
    - name: pyramid
      type: pyramid_rate_f32
      config:
        resolution: 8
        telescope_diameter_m: 2.0
        pupil_samples: 8
        input_quantity: electric_field
        field_schema: test.pyramid.field/1
        source_wavelength_m: 7.5e-7
        source_photon_irradiance_m2_s: 9.0
        photon_rate_schema: test.pyramid.rate/1
    inputs:
    - name: field
      destination: pyramid.electric_field
      binding: field
    outputs:
    - name: rate
      source: pyramid.photon_rate
    """
    mktemp() do path, io
        write(io, text)
        close(io)
        field = ones(ComplexF32, 8, 8)
        definition = load_algorithm_graph(path; bindings=(; field))
        @test definition.nodes[1].config.input_quantity === :electric_field
        @test definition.nodes[1].config.field_schema == "test.pyramid.field/1"
        @test definition.nodes[1].config.opd_schema == ""
        graph = prepare_algorithm_graph(definition)
        @test prepared_graph_node(graph, Val(:pyramid)).field_input === field
        step_graph!(graph)
        @test sum(graph_output(graph, Val(:rate))) ≈ 36f0 rtol=1f-5
    end
end
