using Test
using AdaptiveOpticsSim
using FixedSizeArrays: FixedSizeArray

function snapshot_test_graph(target; execution=AdaptiveOpticsSim.AlgorithmGraphs.StreamGraphExecution())
    graphs = AdaptiveOpticsSim.AlgorithmGraphs
    backends = AdaptiveOpticsSim.Backends
    first = backends.allocate_device_array(target, Float32, 8, 8)
    second = backends.allocate_device_array(target, Float32, 8, 8)
    fill!(first, 1.0f0)
    fill!(second, 2.0f0)
    schema = "test.snapshot.opd.f32/1"
    definition = graphs.algorithm_graph((graphs.pupil_opd_composition_node(:pupil;
        resolution=8, uncompensated_opd_schema=schema, surface_opd_schema=schema,
        pupil_opd_schema=schema),);
        name=:snapshot_test,
        inputs=(graphs.graph_input(:first, :pupil => :uncompensated_opd, first),
            graphs.graph_input(:second, :pupil => :surface_opd, second)),
        outputs=(graphs.graph_output(:opd, :pupil => :pupil_opd),))
    return graphs.prepare_algorithm_graph(definition; target, execution)
end

function snapshot_wait!(snapshot)
    graphs = AdaptiveOpticsSim.AlgorithmGraphs
    deadline = time_ns() + 5_000_000_000
    while !graphs.poll_graph_output_snapshot!(snapshot)
        time_ns() < deadline || error("snapshot completion timed out")
        yield()
    end
    return nothing
end

function snapshot_cycle!(snapshot)
    graphs = AdaptiveOpticsSim.AlgorithmGraphs
    graphs.try_snapshot_graph_output!(snapshot) || error("snapshot is occupied")
    snapshot_wait!(snapshot)
    graphs.release_graph_output_snapshot!(snapshot)
    return nothing
end

function check_output_snapshot(target; execution=AdaptiveOpticsSim.AlgorithmGraphs.StreamGraphExecution())
    graphs = AdaptiveOpticsSim.AlgorithmGraphs
    graph = snapshot_test_graph(target; execution)
    snapshot = graphs.prepare_graph_output_snapshot(graph, Val(:opd))
    values = graphs.snapshot_values(snapshot)
    @test size(values) == (8, 8)
    @test !graphs.poll_graph_output_snapshot!(snapshot)
    @test_throws graphs.AlgorithmGraphError graphs.snapshot_metadata(snapshot)
    @test_throws graphs.AlgorithmGraphError graphs.release_graph_output_snapshot!(snapshot)
    # The actual consumer task can be warmed before admitting the first frame.
    @test graphs.try_snapshot_graph_output!(snapshot)
    snapshot_wait!(snapshot)
    @test graphs.snapshot_metadata(snapshot).sequence == 0
    graphs.release_graph_output_snapshot!(snapshot)
    graphs.step_graph!(graph)
    @test graphs.try_snapshot_graph_output!(snapshot;
        generation=UInt64(7), pts=Int64(123), model_time_ns=Int64(456))
    @test !graphs.try_snapshot_graph_output!(snapshot; generation=UInt64(99))
    # Advance the plant before the observer starts its transfer: the retained
    # sample must still contain frame one rather than the overwritten output.
    fill!(graphs.graph_input(graph, Val(:first)), 10.0f0)
    graphs.step_graph!(graph)
    snapshot_wait!(snapshot)
    @test graphs.snapshot_values(snapshot) === values
    @test all(==(3.0f0), values)
    @test graphs.snapshot_metadata(snapshot) ==
        (sequence=UInt64(1), generation=UInt64(7), pts=Int64(123), model_time_ns=Int64(456))
    @test graphs.poll_graph_output_snapshot!(snapshot)
    @test !graphs.try_snapshot_graph_output!(snapshot)
    graphs.release_graph_output_snapshot!(snapshot)
    @test graphs.try_snapshot_graph_output!(snapshot)
    snapshot_wait!(snapshot)
    @test all(==(12.0f0), values)
    @test graphs.snapshot_metadata(snapshot).sequence == 2
    graphs.release_graph_output_snapshot!(snapshot)

    snapshot_cycle!(snapshot)
    snapshot_cycle!(snapshot)
    @test (@allocated snapshot_cycle!(snapshot)) == 0
    @test graphs.try_snapshot_graph_output!(snapshot)
    graphs.reset_graph_output_snapshot!(snapshot)
    @test !graphs.poll_graph_output_snapshot!(snapshot)
    graphs.reset_graph!(graph)
    @test graphs.try_snapshot_graph_output!(snapshot; generation=UInt64(8))
    snapshot_wait!(snapshot)
    @test graphs.snapshot_metadata(snapshot).generation == 8
    @test graphs.snapshot_metadata(snapshot).sequence == 0
    graphs.release_graph_output_snapshot!(snapshot)
    @test graphs.try_snapshot_graph_output!(snapshot)
    graphs.poll_graph_output_snapshot!(snapshot)
    close(snapshot)
    close(snapshot)
    @test !graphs.try_snapshot_graph_output!(snapshot)
    @test !graphs.poll_graph_output_snapshot!(snapshot)
    @test_throws graphs.AlgorithmGraphError graphs.reset_graph_output_snapshot!(snapshot)
end

@testset "Prepared graph output snapshot CPU" begin
    check_output_snapshot(AdaptiveOpticsSim.Backends.HostComputeDevice())
end

function snapshot_imaging_allocations(imaging)
    AdaptiveOpticsSim.Optics.form_direct_image!(imaging)
    return @allocated AdaptiveOpticsSim.Optics.form_direct_image!(imaging)
end

function check_snapshot_imaging(values)
    optics = AdaptiveOpticsSim.Optics
    backends = AdaptiveOpticsSim.Backends
    @test backends.backend(values) == backends.CPUBackend()
    @test backends.compute_device(values) == backends.HostComputeDevice()
    telescope = optics.Telescope(resolution=8, diameter=1.22,
        central_obstruction=0.204918, T=Float32)
    source = optics.Source(band=:custom, wavelength=8.0f-7,
        photon_irradiance=2.0f0, T=Float32)
    pupil = optics.PupilFunction(telescope, values)
    @test optics.opd_map(pupil) === values
    imaging = optics.prepare_direct_imaging(pupil, source; zero_padding=2)
    reference = optics.prepare_direct_imaging(
        optics.PupilFunction(telescope, Array(values)), source; zero_padding=2)
    optics.form_direct_image!(imaging)
    optics.form_direct_image!(reference)
    @test optics.intensity_values(optics.direct_imaging_output(imaging)) ≈
        optics.intensity_values(optics.direct_imaging_output(reference))
    @test snapshot_imaging_allocations(imaging) == 0
    @test snapshot_imaging_allocations(reference) == 0
end

@testset "Host fixed-size snapshot optics" begin
    graphs = AdaptiveOpticsSim.AlgorithmGraphs
    graph = snapshot_test_graph(AdaptiveOpticsSim.Backends.HostComputeDevice())
    fill!(graphs.graph_input(graph, Val(:first)), 1.0f-8)
    fill!(graphs.graph_input(graph, Val(:second)), 2.0f-8)
    graphs.step_graph!(graph)
    snapshot = graphs.prepare_graph_output_snapshot(graph, Val(:opd))
    @test graphs.try_snapshot_graph_output!(snapshot)
    snapshot_wait!(snapshot)
    check_snapshot_imaging(graphs.snapshot_values(snapshot))
    # Explicit Vector backing is also host storage; default backing is selected
    # exclusively through FixedSizeArrays' public alias, without naming it here.
    check_snapshot_imaging(FixedSizeArray{Float32,2,Vector{Float32}}(
        graphs.snapshot_values(snapshot)))
    graphs.release_graph_output_snapshot!(snapshot)
    close(snapshot)
end

if get(ENV, "ADAPTIVEOPTICS_TEST_SNAPSHOT_CUDA", "false") == "true"
    @eval import CUDA
    @testset "Prepared graph output snapshot CUDA" begin
        @test CUDA.functional()
        @test_throws AdaptiveOpticsSim.InvalidConfiguration AdaptiveOpticsSim.Backends.array_backend_selector(
                FixedSizeArray{Float32,2,CUDA.CuArray{Float32,1,CUDA.DeviceMemory}})
        check_output_snapshot(
            AdaptiveOpticsSim.Backends.AcceleratorComputeDevice(
                AdaptiveOpticsSim.Backends.CUDABackend(), 0);
            execution=AdaptiveOpticsSim.AlgorithmGraphs.CapturedGraphExecution())
    end
end
