const AOG = AdaptiveOpticsSim.AlgorithmGraphs

function warmed_graph_step_allocation_bytes(graph)
    step_graph!(graph)
    return @allocated step_graph!(graph)
end

struct GraphTestGainDeclaration end

struct GraphTestGainConfiguration
    extent::Int
end

struct GraphTestGainOwner{G,I,O}
    gain::G
    input::I
    output::O
end

function AOG.graph_node_ports(
    ::Type{GraphTestGainDeclaration},
    configuration::GraphTestGainConfiguration,
)
    configuration.extent > 0 || throw(ArgumentError("extent must be positive"))
    shape = (configuration.extent,)
    schema = "test.graph.signal.f32/1"
    return (
        AOG.graph_port_contract(
            :input,
            :input,
            :data,
            Float32,
            shape,
            schema,
            :column_major,
        ),
        AOG.graph_port_contract(
            :output,
            :output,
            :data,
            Float32,
            shape,
            schema,
            :column_major,
        ),
        AOG.graph_port_contract(
            :gain,
            :input,
            :parameter,
            Float32,
            shape,
            "test.graph.gain.f32/1",
            :column_major,
        ),
    )
end

function AOG.prepare_graph_node(
    ::Type{GraphTestGainDeclaration},
    ::GraphTestGainConfiguration,
    props,
    inputs::NamedTuple,
    outputs::NamedTuple,
    parameters::NamedTuple,
    target,
)
    return GraphTestGainOwner(parameters.gain, inputs.input, outputs.output)
end

function AOG.step_graph_node!(owner::GraphTestGainOwner)
    @inbounds for index in eachindex(owner.output, owner.input, owner.gain)
        owner.output[index] = owner.input[index] * owner.gain[index]
    end
    return nothing
end

AOG.reset_graph_node!(::GraphTestGainOwner) = nothing

struct GraphTestAddDeclaration end

struct GraphTestAddConfiguration
    extent::Int
    schema::String
    increment::Float32
end

struct GraphTestAddOwner{I,O}
    increment::Float32
    input::I
    output::O
end

function AOG.graph_node_ports(
    ::Type{GraphTestAddDeclaration},
    configuration::GraphTestAddConfiguration,
)
    format(name, direction) = AOG.graph_port_contract(
        name,
        direction,
        :data,
        Float32,
        (configuration.extent,),
        configuration.schema,
        :column_major,
    )
    return (format(:input, :input), format(:output, :output))
end

function AOG.prepare_graph_node(
    ::Type{GraphTestAddDeclaration},
    configuration::GraphTestAddConfiguration,
    props,
    inputs::NamedTuple,
    outputs::NamedTuple,
    parameters::NamedTuple,
    target,
)
    return GraphTestAddOwner(
        configuration.increment,
        inputs.input,
        outputs.output,
    )
end

function AOG.step_graph_node!(owner::GraphTestAddOwner)
    @inbounds for index in eachindex(owner.output, owner.input)
        owner.output[index] = owner.input[index] + owner.increment
    end
    return nothing
end

AOG.reset_graph_node!(::GraphTestAddOwner) = nothing

struct GraphTestSourceDeclaration end

struct GraphTestSourceConfiguration
    extent::Int
    value::Float32
end

struct GraphTestSourceOwner{O}
    value::Float32
    output::O
end

function AOG.graph_node_ports(
    ::Type{GraphTestSourceDeclaration},
    configuration::GraphTestSourceConfiguration,
)
    return (
        AOG.graph_port_contract(
            :output,
            :output,
            :data,
            Float32,
            (configuration.extent,),
            "test.graph.signal.f32/1",
            :column_major,
        ),
    )
end

function AOG.prepare_graph_node(
    ::Type{GraphTestSourceDeclaration},
    configuration::GraphTestSourceConfiguration,
    props,
    ::NamedTuple{()},
    outputs::NamedTuple,
    parameters::NamedTuple,
    target,
)
    return GraphTestSourceOwner(configuration.value, outputs.output)
end

function AOG.step_graph_node!(owner::GraphTestSourceOwner)
    fill!(owner.output, owner.value)
    return nothing
end

AOG.reset_graph_node!(::GraphTestSourceOwner) = nothing

struct GraphTestExactBindingDeclaration end

struct GraphTestExactBindingOwner{I,O}
    input::I
    output::O
end

function AOG.graph_node_ports(
    ::Type{GraphTestExactBindingDeclaration},
    extent::Int,
)
    format(name, direction) = AOG.graph_port_contract(
        name,
        direction,
        :data,
        Float32,
        (extent,),
        "test.graph.signal.f32/1",
        :column_major,
    )
    return (format(:input, :input), format(:output, :output))
end

function AOG.prepare_graph_node(
    ::Type{GraphTestExactBindingDeclaration},
    extent::Int,
    props,
    inputs::NamedTuple,
    outputs::NamedTuple,
    parameters::NamedTuple,
    target,
)
    return GraphTestExactBindingOwner(inputs.input, outputs.output)
end

function AOG.step_graph_node!(owner::GraphTestExactBindingOwner)
    copyto!(owner.output, owner.input)
    return nothing
end

AOG.reset_graph_node!(::GraphTestExactBindingOwner) = nothing

struct GraphTestSinkDeclaration end

struct GraphTestSinkConfiguration
    extent::Int
end

mutable struct GraphTestSinkOwner{I}
    sum::Float32
    input::I
end

function AOG.graph_node_ports(
    ::Type{GraphTestSinkDeclaration},
    configuration::GraphTestSinkConfiguration,
)
    return (
        AOG.graph_port_contract(
            :input,
            :input,
            :data,
            Float32,
            (configuration.extent,),
            "test.graph.signal.f32/1",
            :column_major,
        ),
    )
end

function AOG.prepare_graph_node(
    ::Type{GraphTestSinkDeclaration},
    configuration::GraphTestSinkConfiguration,
    props,
    inputs::NamedTuple,
    ::NamedTuple{()},
    parameters::NamedTuple,
    target,
)
    return GraphTestSinkOwner(0.0f0, inputs.input)
end

function AOG.step_graph_node!(owner::GraphTestSinkOwner)
    total = 0.0f0
    @inbounds for value in owner.input
        total += value
    end
    owner.sum = total
    return nothing
end

function AOG.reset_graph_node!(owner::GraphTestSinkOwner)
    owner.sum = 0.0f0
    return nothing
end

struct GraphTestFailureDeclaration end

struct GraphTestFailureOwner{I,O}
    input::I
    output::O
end

function AOG.graph_node_ports(
    ::Type{GraphTestFailureDeclaration},
    extent::Int,
)
    format(name, direction) = AOG.graph_port_contract(
        name,
        direction,
        :data,
        Float32,
        (extent,),
        "test.graph.signal.f32/1",
        :column_major,
    )
    return (format(:input, :input), format(:output, :output))
end

function AOG.prepare_graph_node(
    ::Type{GraphTestFailureDeclaration},
    extent::Int,
    props,
    inputs::NamedTuple,
    outputs::NamedTuple,
    parameters::NamedTuple,
    target,
)
    return GraphTestFailureOwner(inputs.input, outputs.output)
end

function AOG.step_graph_node!(owner::GraphTestFailureOwner)
    any(<(0.0f0), owner.input) && throw(DomainError(
        owner.input,
        "negative graph fixture input",
    ))
    copyto!(owner.output, owner.input)
    return nothing
end

AOG.reset_graph_node!(::GraphTestFailureOwner) = nothing

@testset "portable algorithm graph direct links and sparse parameters" begin
    input = Float32[1, 2, 3]
    definition = algorithm_graph(
        (
            algorithm_node(
                :gain,
                GraphTestGainDeclaration,
                GraphTestGainConfiguration(3),
            ),
            algorithm_node(
                :sink,
                GraphTestSinkDeclaration,
                GraphTestSinkConfiguration(3),
            ),
        );
        inputs=(graph_input(:residual, :gain => :input, input),),
        outputs=(graph_output(:command, :gain => :output),),
        links=(link(:gain => :output, :sink => :input),),
        parameters=(sparse_parameter(:gain => :gain, Float32[2, 3, 4]),),
    )
    target = AdaptiveOpticsSim.Backends.HostComputeDevice()
    graph = prepare_algorithm_graph(definition; target=target)

    @test AdaptiveOpticsSim.Backends.compute_device(graph) == target
    @test AdaptiveOpticsSim.Backends.compute_device(
        graph_output(graph, Val(:command)),
    ) ==
        target
    @test graph_input(graph, Val(:residual)) === input
    @test graph_input(graph, :residual) === input
    @test prepared_graph_node(graph, Val(:gain)) isa GraphTestGainOwner
    @test iszero(graph_step_sequence(graph))
    @test all(iszero, graph_output(graph, Val(:command)))
    @test step_graph!(graph) === graph
    @test graph_output(graph, Val(:command)) == Float32[2, 6, 12]
    @test graph_output(graph, :command) === graph_output(graph, Val(:command))
    @test prepared_graph_node(graph, Val(:sink)).sum == 20.0f0
    @test graph_step_sequence(graph) == UInt64(1)
    @test !graph_failed(graph)

    @test warmed_graph_step_allocation_bytes(graph) == 0
    @test @inferred(step_graph!(graph)) === graph
    @test all(isconcretetype, fieldtypes(typeof(graph)))
end

@testset "portable algorithm graph delayed feedback and reset" begin
    configuration = GraphTestAddConfiguration(
        2,
        "test.graph.signal.f32/1",
        1.0f0,
    )
    definition = algorithm_graph(
        (
            algorithm_node(:first, GraphTestAddDeclaration, configuration),
            algorithm_node(:second, GraphTestAddDeclaration, configuration),
        );
        outputs=(graph_output(:result, :second => :output),),
        links=(link(:first => :output, :second => :input),),
        delayed_links=(delayed_link(
            :second => :output,
            :first => :input,
            zeros(Float32, 2),
        ),),
    )
    graph = prepare_algorithm_graph(definition)

    step_graph!(graph)
    @test graph_output(graph, Val(:result)) == Float32[2, 2]
    step_graph!(graph)
    @test graph_output(graph, Val(:result)) == Float32[4, 4]
    @test graph_step_sequence(graph) == UInt64(2)
    @test @allocated(step_graph!(graph)) == 0

    @test reset_graph!(graph) === graph
    @test graph_step_sequence(graph) == UInt64(0)
    @test !graph_failed(graph)
    step_graph!(graph)
    @test graph_output(graph, Val(:result)) == Float32[2, 2]
end

@testset "portable algorithm graph admits pure sources and sinks" begin
    source = prepare_algorithm_graph(algorithm_graph(
        (algorithm_node(
            :source,
            GraphTestSourceDeclaration,
            GraphTestSourceConfiguration(2, 7.0f0),
        ),);
        outputs=(graph_output(:sample, :source => :output),),
    ))
    step_graph!(source)
    @test graph_output(source, Val(:sample)) == Float32[7, 7]

    input = Float32[2, 5]
    sink = prepare_algorithm_graph(algorithm_graph(
        (algorithm_node(
            :sink,
            GraphTestSinkDeclaration,
            GraphTestSinkConfiguration(2),
        ),);
        inputs=(graph_input(:sample, :sink => :input, input),),
    ))
    step_graph!(sink)
    @test isempty(sink.outputs)
    @test prepared_graph_node(sink, Val(:sink)).sum == 7.0f0
end

@testset "portable algorithm graph binds exact node execution owners" begin
    input = Float32[3, 5]
    graph = prepare_algorithm_graph(algorithm_graph(
        (algorithm_node(:bound, GraphTestExactBindingDeclaration, 2),);
        inputs=(graph_input(:sample, :bound => :input, input),),
        outputs=(graph_output(:copy, :bound => :output),),
    ))

    owner = prepared_graph_node(graph, Val(:bound))
    @test owner isa GraphTestExactBindingOwner
    @test owner.input === input
    @test owner.output === graph_output(graph, Val(:copy))
    @test step_graph!(graph) === graph
    @test graph_output(graph, Val(:copy)) == input
    @test @allocated(step_graph!(graph)) == 0
    @test @inferred(step_graph!(graph)) === graph
end

@testset "portable algorithm graph failure is fail-stop" begin
    input = Float32[1, -1]
    graph = prepare_algorithm_graph(algorithm_graph(
        (algorithm_node(:checked, GraphTestFailureDeclaration, 2),);
        inputs=(graph_input(:input, :checked => :input, input),),
        outputs=(graph_output(:output, :checked => :output),),
    ))

    @test_throws DomainError step_graph!(graph)
    @test graph_failed(graph)
    @test graph_step_sequence(graph) == UInt64(0)
    @test_throws AlgorithmGraphError step_graph!(graph)
    fill!(input, 3.0f0)
    reset_graph!(graph)
    step_graph!(graph)
    @test graph_output(graph, Val(:output)) == input
end

@testset "native controller and modal OPD graph nodes" begin
    residual = Float32[1, 2]
    basis = zeros(Float32, 2, 2, 2)
    fill!(@view(basis[:, :, 1]), 1.0f0)
    fill!(@view(basis[:, :, 2]), 2.0f0)
    pupil_support = Bool[true false; true true]
    coefficient_schema = "test.graph.modal-coefficients.f32/1"

    definition = algorithm_graph(
        (
            discrete_integrator_node(
                :controller;
                extent=2,
                sample_period_s=0.1f0,
                input_schema="test.graph.residual.f32/1",
                output_schema=coefficient_schema,
                gain=2.0f0,
                tau_s=0.2f0,
            ),
            modal_opd_expansion_node(
                :modal_opd;
                pupil_rows=2,
                pupil_columns=2,
                mode_count=2,
                coefficients_schema=coefficient_schema,
                opd_schema="test.graph.opd.f32/1",
                basis_schema="test.graph.modal-basis.f32/1",
                pupil_support_schema="test.graph.pupil-support.bool/1",
            ),
        );
        name=:native_modal_control,
        inputs=(graph_input(:residual, :controller => :input, residual),),
        outputs=(
            graph_output(:command, :controller => :output),
            graph_output(:opd, :modal_opd => :opd),
        ),
        links=(link(:controller => :output, :modal_opd => :coefficients),),
        parameters=(
            sparse_parameter(:modal_opd => :basis, basis),
            sparse_parameter(:modal_opd => :pupil_support, pupil_support),
        ),
    )
    graph = prepare_algorithm_graph(definition)

    @test graph_name(definition) === :native_modal_control
    @test graph_name(graph) === :native_modal_control
    @test step_graph!(graph) === graph
    @test graph_output(graph, Val(:command)) ≈ Float32[0.1, 0.2]
    @test graph_output(graph, Val(:opd)) ≈ Float32[0.5 0.0; 0.5 0.5]
    step_graph!(graph)
    @test @allocated(step_graph!(graph)) == 0
    @test @inferred(step_graph!(graph)) === graph
    reset_graph!(graph)
    step_graph!(graph)
    @test graph_output(graph, Val(:command)) ≈ Float32[0.1, 0.2]
end

@testset "diffractive Shack-Hartmann photon-rate graph node" begin
    pupil_opd = zeros(Float32, 16, 16)
    node = shack_hartmann_rate_node(
        :shwfs;
        resolution=16,
        telescope_diameter_m=8.0,
        n_lenslets=4,
        n_pix_subap=4,
        pixel_scale_arcsec=0.1,
        source_wavelength_m=0.75e-6,
        source_photon_irradiance_m2_s=1.0,
        opd_schema="test.graph.pupil-opd.f32/1",
        photon_rate_schema="test.graph.shwfs-photon-rate.f32/1",
    )
    definition = algorithm_graph(
        (node,);
        name=:shwfs_photon_rate,
        inputs=(graph_input(:pupil_opd, :shwfs => :opd, pupil_opd),),
        outputs=(
            graph_output(:shwfs_photon_rate, :shwfs => :photon_rate),
        ),
    )
    graph = prepare_algorithm_graph(definition)
    owner = prepared_graph_node(graph, Val(:shwfs))
    output = graph_output(graph, Val(:shwfs_photon_rate))

    @test owner.prepared.input.opd === pupil_opd
    @test owner.prepared.output.values === output
    @test size(output) == (16, 16)
    step_graph!(graph)
    @test all(isfinite, output)
    @test sum(output) > 0
    @test warmed_graph_step_allocation_bytes(graph) == 0
    @test @inferred(step_graph!(graph)) === graph

    @test_throws AlgorithmGraphError shack_hartmann_rate_node(
        :invalid_shwfs;
        resolution=16,
        telescope_diameter_m=8.0,
        n_lenslets=4,
        n_pix_subap=3,
        pixel_scale_arcsec=0.1,
        source_wavelength_m=0.75e-6,
        source_photon_irradiance_m2_s=1.0,
        opd_schema="test.graph.pupil-opd.f32/1",
        photon_rate_schema="test.graph.shwfs-photon-rate.f32/1",
    )
end

@testset "single-read CCD detector-acquisition graph node" begin
    photon_rate = fill(1_000.0f0, 4, 4)
    node = ccd_detector_acquisition_node(
        :detector;
        rows=4,
        columns=4,
        pixel_scale_arcsec=0.1,
        wavelength_m=0.75e-6,
        exposure_duration_s=0.01,
        quantum_efficiency=0.5,
        rng_seed=0x1234,
        photon_rate_schema="test.graph.shwfs-photon-rate.f32/1",
        frame_schema="test.graph.shwfs-frame.f32/1",
        photon_noise=true,
        readout_noise=true,
        readout_noise_e=2.7,
    )
    definition = algorithm_graph(
        (node,);
        name=:ccd_detector_acquisition,
        inputs=(
            graph_input(
                :photon_rate,
                :detector => :photon_rate,
                photon_rate,
            ),
        ),
        outputs=(graph_output(:frame, :detector => :frame),),
    )
    graph = prepare_algorithm_graph(definition)
    owner = prepared_graph_node(graph, Val(:detector))
    frame = graph_output(graph, Val(:frame))

    @test owner.prepared.input.values === photon_rate
    @test owner.output === frame
    @test size(frame) == (4, 4)
    step_graph!(graph)
    first_frame = copy(frame)
    @test all(isfinite, frame)
    @test warmed_graph_step_allocation_bytes(graph) == 0
    @test @inferred(step_graph!(graph)) === graph
    @test frame != first_frame
    reset_graph!(graph)
    @test all(iszero, frame)
    step_graph!(graph)
    @test frame == first_frame

    deterministic = ccd_detector_acquisition_node(
        :detector;
        rows=4,
        columns=4,
        binning=2,
        pixel_scale_arcsec=0.1,
        wavelength_m=0.75e-6,
        exposure_duration_s=0.01,
        quantum_efficiency=0.5,
        rng_seed=1,
        photon_rate_schema="test.graph.shwfs-photon-rate.f32/1",
        frame_schema="test.graph.shwfs-frame.f32/1",
        photon_noise=false,
    )
    deterministic_graph = prepare_algorithm_graph(algorithm_graph(
        (deterministic,);
        inputs=(
            graph_input(
                :photon_rate,
                :detector => :photon_rate,
                photon_rate,
            ),
        ),
        outputs=(graph_output(:frame, :detector => :frame),),
    ))
    step_graph!(deterministic_graph)
    @test graph_output(deterministic_graph, Val(:frame)) ==
        fill(20.0f0, 2, 2)

    @test_throws AlgorithmGraphError ccd_detector_acquisition_node(
        :invalid_detector;
        rows=4,
        columns=4,
        pixel_scale_arcsec=0.1,
        wavelength_m=0.75e-6,
        exposure_duration_s=0.01,
        quantum_efficiency=0.5,
        rng_seed=1,
        photon_rate_schema="test.graph.shwfs-photon-rate.f32/1",
        frame_schema="test.graph.shwfs-frame.f32/1",
        photon_noise=false,
        readout_noise=false,
        readout_noise_e=2.7,
    )
    invalid_rate = copy(photon_rate)
    invalid_rate[1, 1] = -1
    invalid_definition = algorithm_graph(
        (node,);
        inputs=(
            graph_input(
                :photon_rate,
                :detector => :photon_rate,
                invalid_rate,
            ),
        ),
    )
    @test_throws InvalidConfiguration prepare_algorithm_graph(
        invalid_definition,
    )
end

@testset "calibrated Shack-Hartmann centroid graph node" begin
    frame = zeros(Float32, 4, 4)
    frame[1, 2] = 1
    frame[4, 1] = 1
    frame[2, 4] = 1
    frame[3, 3] = 1
    valid_subapertures = Bool[true false; true true]
    reference_signal = zeros(Float32, 4, 2)
    reference_signal[2, 1] = 0.25f0
    reference_signal[1, 2] = 0.5f0
    node = shack_hartmann_centroid_node(
        :centroid;
        resolution=4,
        telescope_diameter_m=2.0,
        n_lenslets=2,
        n_pix_subap=2,
        centroid_cutoff_fraction=0.1,
        centroid_response=2.0,
        calibration_wavelength_m=0.75e-6,
        calibration_signature=7,
        frame_schema="test.graph.shwfs-frame.f32/1",
        slopes_schema="test.graph.shwfs-centroid-slopes.f32/1",
        valid_subapertures_schema=
            "test.graph.shwfs-valid-subapertures.bool/1",
        reference_signal_schema=
            "test.graph.shwfs-centroid-reference.f32/1",
    )
    definition = algorithm_graph(
        (node,);
        name=:shwfs_centroid,
        inputs=(graph_input(:frame, :centroid => :frame, frame),),
        outputs=(graph_output(:slopes, :centroid => :slopes),),
        parameters=(
            sparse_parameter(
                :centroid => :valid_subapertures,
                valid_subapertures,
            ),
            sparse_parameter(
                :centroid => :reference_signal,
                reference_signal,
            ),
        ),
    )
    graph = prepare_algorithm_graph(definition)
    owner = prepared_graph_node(graph, Val(:centroid))
    slopes_output = graph_output(graph, Val(:slopes))

    fill!(valid_subapertures, false)
    fill!(reference_signal, 100)
    step_graph!(graph)

    @test owner.observation.storage === frame
    @test owner.measurement.storage === slopes_output
    @test slopes_output == Float32[0, 0.375, 0, 0, 0.25, 0, 0, 0]
    @test n_valid_subapertures(owner.prepared.sensor.front_end.layout) == 3
    @test warmed_graph_step_allocation_bytes(graph) == 0
    @test @inferred(step_graph!(graph)) === graph
    reset_graph!(graph)
    @test all(iszero, slopes_output)

    sensor = owner.prepared.sensor
    layout_revision = subaperture_layout_revision(sensor.front_end.layout)
    @test_throws DimensionMismatchError set_valid_subapertures!(
        sensor,
        fill(true, 1, 1),
    )
    mask_parent = fill(false, 4, 4)
    mask_parent[2:3, 2:3] .= Bool[true false; false true]
    mask_view = @view mask_parent[2:3, 2:3]
    @test set_valid_subapertures!(sensor, mask_view) === sensor
    @test subaperture_layout_revision(sensor.front_end.layout) ==
        layout_revision + UInt(1)
    @test n_valid_subapertures(sensor.front_end.layout) == 2
    @test !sensor.calibration.calibrated

    @test_throws AlgorithmGraphError shack_hartmann_centroid_node(
        :invalid_centroid;
        resolution=4,
        telescope_diameter_m=2.0,
        n_lenslets=2,
        n_pix_subap=3,
        centroid_cutoff_fraction=0.1,
        centroid_response=1.0,
        calibration_wavelength_m=0.75e-6,
        calibration_signature=0,
        frame_schema="test.graph.shwfs-frame.f32/1",
        slopes_schema="test.graph.shwfs-centroid-slopes.f32/1",
        valid_subapertures_schema=
            "test.graph.shwfs-valid-subapertures.bool/1",
        reference_signal_schema=
            "test.graph.shwfs-centroid-reference.f32/1",
    )
end

@testset "ordered Shack-Hartmann slope-selection graph node" begin
    full_slopes = Float32[1, 2, 3, 4, 11, 12, 13, 14]
    lenslet_order = UInt32[4, 1, 3]
    node = shack_hartmann_slope_selection_node(
        :selection;
        n_lenslets=2,
        selected_lenslet_count=3,
        full_slopes_schema="test.graph.shwfs-full-slopes.f32/1",
        selected_slopes_schema="test.graph.shwfs-selected-slopes.f32/1",
        lenslet_order_schema="test.graph.shwfs-lenslet-order.u32/1",
    )
    definition = algorithm_graph(
        (node,);
        name=:shwfs_slope_selection,
        inputs=(
            graph_input(
                :full_slopes,
                :selection => :full_slopes,
                full_slopes,
            ),
        ),
        outputs=(
            graph_output(
                :selected_slopes,
                :selection => :selected_slopes,
            ),
        ),
        parameters=(
            sparse_parameter(
                :selection => :lenslet_order,
                lenslet_order,
            ),
        ),
    )
    graph = prepare_algorithm_graph(definition)
    owner = prepared_graph_node(graph, Val(:selection))
    selected_slopes = graph_output(graph, Val(:selected_slopes))

    fill!(lenslet_order, UInt32(1))
    step_graph!(graph)
    @test owner.input === full_slopes
    @test owner.output === selected_slopes
    @test selected_slopes == Float32[4, 14, 1, 11, 3, 13]
    @test warmed_graph_step_allocation_bytes(graph) == 0
    @test @inferred(step_graph!(graph)) === graph
    reset_graph!(graph)
    @test all(iszero, selected_slopes)

    duplicate_definition = algorithm_graph(
        (node,);
        inputs=(graph_input(
            :full_slopes,
            :selection => :full_slopes,
            full_slopes,
        ),),
        outputs=(graph_output(
            :selected_slopes,
            :selection => :selected_slopes,
        ),),
        parameters=(sparse_parameter(
            :selection => :lenslet_order,
            UInt32[1, 1, 2],
        ),),
    )
    @test_throws AlgorithmGraphError prepare_algorithm_graph(
        duplicate_definition,
    )
end

@testset "TOML graph files compile to native graph definitions" begin
    residual = Float32[1, 2]
    basis = zeros(Float32, 2, 2, 2)
    fill!(@view(basis[:, :, 1]), 1.0f0)
    fill!(@view(basis[:, :, 2]), 2.0f0)
    pupil_support = Bool[true false; true true]
    path = joinpath(
        dirname(@__DIR__),
        "graph_files",
        "native_modal_control.toml",
    )

    definition = load_algorithm_graph(
        path;
        bindings=(;
            residual,
            basis,
            pupil_support,
        ),
    )
    @test graph_name(definition) === :native_modal_control
    @test isconcretetype(typeof(definition))
    @test all(isconcretetype, fieldtypes(typeof(definition)))
    graph = prepare_algorithm_graph(definition)
    step_graph!(graph)
    @test graph_output(graph, Val(:command)) ≈ Float32[0.1, 0.2]
    @test graph_output(graph, Val(:opd)) ≈ Float32[0.5 0.0; 0.5 0.5]
    step_graph!(graph)
    @test @allocated(step_graph!(graph)) == 0
    @test @inferred(step_graph!(graph)) === graph

    @test_throws AlgorithmGraphError load_algorithm_graph(
        path;
        bindings=(; residual, basis),
    )
    @test keys(builtin_graph_node_types()) == (
        :ccd_detector_acquisition_f32,
        :discrete_integrator_f32,
        :modal_opd_expansion_f32,
        :shack_hartmann_centroid_f32,
        :shack_hartmann_rate_f32,
        :shack_hartmann_slope_selection_f32,
    )

    mktemp() do invalid_path, io
        write(io, "schema_version = 2\nname = \"invalid_version\"\nnodes = []\n")
        close(io)
        @test_throws AlgorithmGraphError load_algorithm_graph(invalid_path)
    end
    mktemp() do invalid_path, io
        write(
            io,
            "schema_version = 1\nname = \"unknown_field\"\n" *
            "unexpected = true\nnodes = []\n",
        )
        close(io)
        @test_throws AlgorithmGraphError load_algorithm_graph(invalid_path)
    end
    mktemp() do invalid_path, io
        write(
            io,
            "schema_version = 1\nname = \"unknown_type\"\n" *
            "[[nodes]]\nname = \"node\"\ntype = \"not_registered\"\n" *
            "[nodes.config]\nextent = 2\n",
        )
        close(io)
        @test_throws AlgorithmGraphError load_algorithm_graph(invalid_path)
    end

    mktemp() do custom_path, io
        write(
            io,
            "schema_version = 1\nname = \"custom_type\"\n" *
            "[[nodes]]\nname = \"source\"\ntype = \"test_source\"\n" *
            "[nodes.config]\nextent = 2\nvalue = 7.0\n" *
            "[[outputs]]\nname = \"sample\"\nsource = \"source.output\"\n",
        )
        close(io)
        test_source = function (name, config, props)
            isempty(props) || error("test-source props must be empty")
            return algorithm_node(
                name,
                GraphTestSourceDeclaration,
                GraphTestSourceConfiguration(
                    Int(config.extent),
                    Float32(config.value),
                ),
            )
        end
        node_types = merge(
            builtin_graph_node_types(),
            (; test_source),
        )
        custom = prepare_algorithm_graph(load_algorithm_graph(
            custom_path;
            node_types,
        ))
        step_graph!(custom)
        @test graph_output(custom, Val(:sample)) == Float32[7, 7]
    end
end


@testset "REVOLT Classic SHWFS TOML path is executable" begin
    pupil_opd = zeros(Float32, 240, 240)
    valid_subapertures = fill(true, 16, 16)
    reference_signal = zeros(Float32, 16 * 16, 2)
    # Executable wiring fixture only. Applications bind the ROI-derived order.
    lenslet_order = UInt32.(1:188)
    path = joinpath(
        dirname(dirname(@__DIR__)),
        "examples",
        "graphs",
        "revolt_classic_shwfs.toml",
    )
    definition = load_algorithm_graph(
        path;
        bindings=(;
            pupil_opd,
            valid_subapertures,
            reference_signal,
            lenslet_order,
        ),
    )
    graph = prepare_algorithm_graph(definition)
    centroid_owner = prepared_graph_node(graph, Val(:centroid))
    step_graph!(graph)
    photon_rate = graph_output(graph, Val(:shwfs_photon_rate))
    frame = graph_output(graph, Val(:shwfs_frame))
    full_slopes = graph_output(graph, Val(:shwfs_full_slopes))
    selected_slopes = graph_output(graph, Val(:shwfs_selected_slopes))

    @test graph_name(graph) === :revolt_classic_shwfs
    @test size(photon_rate) == (352, 352)
    @test size(frame) == (352, 352)
    @test size(full_slopes) == (512,)
    @test size(selected_slopes) == (376,)
    @test size(centroid_owner.prepared.workspace.centroid_host) == (22, 22)
    @test all(isfinite, photon_rate)
    @test all(isfinite, frame)
    @test all(isfinite, full_slopes)
    @test all(isfinite, selected_slopes)
    @test selected_slopes[1:6] == Float32[
        full_slopes[1],
        full_slopes[257],
        full_slopes[2],
        full_slopes[258],
        full_slopes[3],
        full_slopes[259],
    ]
    @test sum(photon_rate) > 0
    @test sum(frame) > 0
end

@testset "portable algorithm graph rejects invalid topology" begin
    configuration = GraphTestAddConfiguration(
        2,
        "test.graph.signal.f32/1",
        1.0f0,
    )
    node = algorithm_node(:add, GraphTestAddDeclaration, configuration)

    @test_throws AlgorithmGraphError algorithm_graph(())
    @test_throws AlgorithmGraphError algorithm_graph((node, node))
    @test_throws AlgorithmGraphError prepare_algorithm_graph(algorithm_graph((node,)))
    @test_throws AlgorithmGraphError prepare_algorithm_graph(algorithm_graph(
        (node,);
        inputs=(graph_input(:input, :add => :missing, zeros(Float32, 2)),),
    ))
    @test_throws AlgorithmGraphError prepare_algorithm_graph(algorithm_graph(
        (node,);
        inputs=(graph_input(:input, :add => :input, zeros(Float32, 3)),),
    ))
    packed_parent = zeros(Float32, 3)
    wrapped_input = @view packed_parent[1:2]
    @test_throws AlgorithmGraphError prepare_algorithm_graph(algorithm_graph(
        (node,);
        inputs=(graph_input(:input, :add => :input, wrapped_input),),
    ))

    backward = algorithm_graph(
        (
            algorithm_node(:first, GraphTestAddDeclaration, configuration),
            algorithm_node(:second, GraphTestAddDeclaration, configuration),
        );
        links=(link(:second => :output, :first => :input),),
        delayed_links=(delayed_link(
            :first => :output,
            :second => :input,
            zeros(Float32, 2),
        ),),
    )
    @test_throws AlgorithmGraphError prepare_algorithm_graph(backward)

    mismatched = GraphTestAddConfiguration(
        2,
        "test.graph.other.f32/1",
        1.0f0,
    )
    schema_mismatch = algorithm_graph(
        (
            algorithm_node(:first, GraphTestAddDeclaration, configuration),
            algorithm_node(:second, GraphTestAddDeclaration, mismatched),
        );
        inputs=(graph_input(:input, :first => :input, zeros(Float32, 2)),),
        links=(link(:first => :output, :second => :input),),
    )
    @test_throws AlgorithmGraphError prepare_algorithm_graph(schema_mismatch)

    duplicate_destination = algorithm_graph(
        (node,);
        inputs=(graph_input(:input, :add => :input, zeros(Float32, 2)),),
        delayed_links=(delayed_link(
            :add => :output,
            :add => :input,
            zeros(Float32, 2),
        ),),
    )
    @test_throws AlgorithmGraphError prepare_algorithm_graph(duplicate_destination)
end

@testset "fixed-step model-time driver" begin
    driver = FixedStepModelTimeDriver(
        PeriodicSchedule(PlantDuration(100); phase=PlantDuration(25));
        origin=PlantTimestamp(1_000),
    )

    @test model_time_sequence(driver) == UInt64(0)
    @test !model_time_exhausted(driver)
    @test next_model_timestamp(driver) == PlantTimestamp(1_025)
    @test advance_model_time!(driver) == PlantTimestamp(1_025)
    @test model_time_sequence(driver) == UInt64(1)
    @test next_model_timestamp(driver) == PlantTimestamp(1_125)
    @test @inferred(advance_model_time!(driver)) == PlantTimestamp(1_125)
    @test @allocated(advance_model_time!(driver)) == 0

    @test reset_model_time!(driver) === driver
    @test model_time_sequence(driver) == UInt64(0)
    @test next_model_timestamp(driver) == PlantTimestamp(1_025)
end

@testset "prepared-boundary model-time driver" begin
    driver = prepare_boundary_model_time_driver((
        PlantTimestamp(10),
        PlantTimestamp(20),
        PlantTimestamp(35),
    ))

    @test next_model_timestamp(driver) == PlantTimestamp(10)
    @test advance_model_time!(driver) == PlantTimestamp(10)
    @test @inferred(advance_model_time!(driver)) == PlantTimestamp(20)
    @test @allocated(advance_model_time!(driver)) == 0
    @test model_time_sequence(driver) == UInt64(3)
    @test model_time_exhausted(driver)
    @test_throws AlgorithmGraphError next_model_timestamp(driver)
    @test_throws AlgorithmGraphError advance_model_time!(driver)
    reset_model_time!(driver)
    @test next_model_timestamp(driver) == PlantTimestamp(10)

    @test_throws AlgorithmGraphError prepare_boundary_model_time_driver(())
    @test_throws AlgorithmGraphError prepare_boundary_model_time_driver((
        PlantTimestamp(10),
        PlantTimestamp(10),
    ))
    @test_throws AlgorithmGraphError prepare_boundary_model_time_driver((
        PlantTimestamp(20),
        PlantTimestamp(10),
    ))
    @test_throws AlgorithmGraphError prepare_boundary_model_time_driver((1, 2))
end

@testset "periodic prepared-boundary model time" begin
    schedule = PeriodicSchedule(
        PlantDuration(1_000);
        phase=PlantDuration(100),
    )
    driver = prepare_boundary_model_time_driver(
        schedule,
        (PlantDuration(0), PlantDuration(250), PlantDuration(750)),
        2;
        origin=PlantTimestamp(10),
    )
    observed = ntuple(_ -> advance_model_time!(driver), 6)
    @test observed == (
        PlantTimestamp(110),
        PlantTimestamp(360),
        PlantTimestamp(860),
        PlantTimestamp(1_110),
        PlantTimestamp(1_360),
        PlantTimestamp(1_860),
    )
    @test model_time_exhausted(driver)

    @test_throws AlgorithmGraphError prepare_boundary_model_time_driver(
        schedule,
        (),
        1,
    )
    @test_throws AlgorithmGraphError prepare_boundary_model_time_driver(
        schedule,
        (PlantDuration(0), PlantDuration(1_000)),
        1,
    )
    @test_throws AlgorithmGraphError prepare_boundary_model_time_driver(
        schedule,
        (PlantDuration(250), PlantDuration(100)),
        1,
    )
    @test_throws AlgorithmGraphError prepare_boundary_model_time_driver(
        schedule,
        (PlantDuration(0),),
        0,
    )
    @test_throws AlgorithmGraphError prepare_boundary_model_time_driver(
        schedule,
        (PlantDuration(0),),
        true,
    )
end

@testset "captured model-time replay" begin
    captures = (
        CapturedModelTimestamp(
            PlantTimestamp(15),
            PlantDuration(2),
            (source=UInt32(7), sequence=UInt64(41)),
        ),
        CapturedModelTimestamp(
            PlantTimestamp(28),
            PlantDuration(3),
            (source=UInt32(7), sequence=UInt64(42)),
        ),
    )
    driver = prepare_captured_model_time_driver(captures)

    @test next_model_time_capture(driver) === captures[1]
    @test model_timestamp(next_model_time_capture(driver)) ==
        PlantTimestamp(15)
    @test model_time_uncertainty(next_model_time_capture(driver)) ==
        PlantDuration(2)
    @test model_time_provenance(next_model_time_capture(driver)) ==
        (source=UInt32(7), sequence=UInt64(41))
    @test advance_model_time!(driver) == PlantTimestamp(15)
    @test @inferred(next_model_time_capture(driver)) === captures[2]
    @test @allocated(advance_model_time!(driver)) == 0
    @test model_time_exhausted(driver)

    reset_model_time!(driver)
    @test next_model_timestamp(driver) == PlantTimestamp(15)
    @test_throws AlgorithmGraphError prepare_captured_model_time_driver(())
    @test_throws AlgorithmGraphError prepare_captured_model_time_driver((
        captures[1],
        CapturedModelTimestamp(
            PlantTimestamp(15),
            PlantDuration(3),
            (source=UInt32(7), sequence=UInt64(42)),
        ),
    ))
    @test_throws AlgorithmGraphError prepare_captured_model_time_driver((
        captures[1],
        CapturedModelTimestamp(
            PlantTimestamp(28),
            PlantDuration(3),
            (stream=UInt32(7), sequence=UInt64(42)),
        ),
    ))
    @test_throws AlgorithmGraphError CapturedModelTimestamp(
        PlantTimestamp(1),
        PlantDuration(0),
        "borrowed provenance",
    )
end

@testset "model-time graph stepping commits together" begin
    graph = prepare_algorithm_graph(algorithm_graph(
        (algorithm_node(
            :source,
            GraphTestSourceDeclaration,
            GraphTestSourceConfiguration(2, 4.0f0),
        ),);
        outputs=(graph_output(:sample, :source => :output),),
    ))
    driver = FixedStepModelTimeDriver(PeriodicSchedule(PlantDuration(50)))

    @test step_graph_at!(graph, driver) == PlantTimestamp(0)
    @test graph_output(graph, Val(:sample)) == Float32[4, 4]
    @test graph_step_sequence(graph) == model_time_sequence(driver) == UInt64(1)
    @test @allocated(step_graph_at!(graph, driver)) == 0

    step_graph!(graph)
    @test_throws AlgorithmGraphError step_graph_at!(graph, driver)

    reset_graph!(graph)
    reset_model_time!(driver)
    @test step_graph_at!(graph, driver) == PlantTimestamp(0)

    reset_graph!(graph)
    captured_driver = prepare_captured_model_time_driver((
        CapturedModelTimestamp(
            PlantTimestamp(75),
            PlantDuration(4),
            (source=UInt32(1), sequence=UInt64(1)),
        ),
    ))
    @test step_graph_at!(graph, captured_driver) == PlantTimestamp(75)
    @test model_time_exhausted(captured_driver)

    failing_input = Float32[-1, 1]
    failing_graph = prepare_algorithm_graph(algorithm_graph(
        (algorithm_node(:checked, GraphTestFailureDeclaration, 2),);
        inputs=(graph_input(:input, :checked => :input, failing_input),),
        outputs=(graph_output(:output, :checked => :output),),
    ))
    failing_driver = FixedStepModelTimeDriver(PeriodicSchedule(PlantDuration(50)))
    @test_throws DomainError step_graph_at!(failing_graph, failing_driver)
    @test model_time_sequence(failing_driver) == UInt64(0)
end
