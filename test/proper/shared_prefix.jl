using AdaptiveOpticsSim, Proper, Test
using AdaptiveOpticsSim.AlgorithmGraphs

function shared_prefix_test_graph(kind; T=Float64, wavelength_um=1.55,
        entrance_radius_m=0.006, sampled_phase=0.3, suffix_opd=false, suffix_role=:branch_opd,
        pixel_integration=false,dm_values=nothing,suffix_values=nothing,opd_values=nothing)
    n = 16
    amplitude = [T(0.7 + 0.2cos(i + j)) for i in 1:n, j in 1:n]
    phase = fill(T(sampled_phase), n, n)
    mask = [T(kind === :fringed ? (i < 13) :
        kind === :unfringed ? (i < 10) : (i >= 10)) for i in 1:n, j in 1:n]
    sampled(role; kwargs...) = (;op="sampled_mask", name=string(role),
        amplitude=string(role), sampling_m=(0.001,0.001),
        center_pixel=(9.0,9.0), outside_amplitude=0.0, outside_phase_rad=0.0, kwargs...)
    dm = (;op="add_opd", name="measured_dm", input="dm_opd", size=(n,n),
        axes=(1,2), sampling_m=(0.001,0.001), center_pixel=(9.0,9.0),
        schema="test.dm-opd-m/1", outside_opd_m=0.0)
    prefix = ((op="circular_aperture",radius_m=entrance_radius_m),
        (op="lens",focal_length_m=0.1), (op="propagate",distance_m=0.03), dm,
        sampled(:apodizer;phase="phase"), (op="propagate",distance_m=0.02))
    suffix = suffix_opd ? (merge(dm,(name=string(suffix_role),input=string(suffix_role),)),) : ()
    steps = (prefix..., sampled(:lyot),
        (op="lens",focal_length_m=kind === :llowfs ? 0.11 : 0.07),
        suffix...,
        (op="propagate",distance_m=kind === :llowfs ? 0.08 : 0.15))
    config = (;resolution=n, pupil_resolution=n, diameter_m=0.016, wavelength_um,
        pupil_grid=(;size=(n,n),axes=(1,2),sampling_m=(0.001,0.001),center_pixel=(9.0,9.0)),
        steps, apodizer="apodizer", phase="phase", lyot="lyot")
    node = proper_optical_train_node_factory(;T)(:optics, config, (;),
        (;apodizer=amplitude, phase, lyot=mask))
    inputs = (graph_input(:pupil_opd,:optics=>:pupil_opd,opd_values === nothing ? zeros(T,n,n) : opd_values),
            graph_input(:pupil_amplitude,:optics=>:pupil_amplitude,ones(T,n,n)),
            graph_input(:dm_opd,:optics=>:dm_opd,dm_values === nothing ? zeros(T,n,n) : dm_values))
    suffix_opd && (inputs = (inputs...,
        graph_input(suffix_role,:optics=>suffix_role,
            suffix_values === nothing ? zeros(T,n,n) : suffix_values)))
    if !pixel_integration
        return prepare_algorithm_graph(algorithm_graph((node,); inputs,
            outputs=(graph_output(:intensity,:optics=>:output),)))
    end
    pitch = optical_sampling_m(node.config)
    pixels = pixel_integration_node(:pixels;input_shape=(n,n),output_shape=(10,12),
        input_sampling_m=(pitch,pitch),output_sampling_m=(2pitch,2pitch),
        input_center_pixel=(9.0,9.0),output_center_pixel=(5.5,6.5),
        input_schema="org.adaptiveopticssim.optical-intensity/1",
        output_schema="test.relative-pixel-flux/1",T)
    return prepare_algorithm_graph(algorithm_graph((node,pixels); inputs,
        links=(link(:optics=>:output,:pixels=>:input),),
        outputs=(graph_output(:intensity,:pixels=>:output),)))
end

function shared_prefix_test_owner(;shared_prefix=true,T=Float64,llowfs_keywords=(;),suffix_opd=false,
        pixel_integration=false,alias_dm=false)
    scc = shared_prefix_test_graph(:fringed;T,suffix_opd)
    unfringed = shared_prefix_test_graph(:unfringed;T,suffix_opd)
    dm_values = alias_dm ? graph_output(scc,:intensity) : nothing
    llowfs = shared_prefix_test_graph(:llowfs;T,suffix_opd,pixel_integration,dm_values,llowfs_keywords...)
    owner = prepare_external_pupil((;scc,llowfs);
        exposure_durations=(scc=ModelDuration(4),llowfs=ModelDuration(2)),
        observations=(scc=(output=:intensity,sampling_node=:optics,
                chopping=(;unfringed,initial_phase=:fringed)),
            llowfs=(output=:intensity,sampling_node=pixel_integration ? :pixels : :optics)))
    shared_prefix && return owner
    # Preserve the original independent propagation as the numerical oracle.
    # This test-only reconstruction supplies no production optimization switch.
    return PreparedExternalPupil(owner.graphs,owner.exposures,owner.pupil_opd,
        owner.pupil_amplitude,owner.observations,owner.input_names,nothing,
        owner.exposure_execution,owner.state)
end

shared_prefix_branches(owner) = (owner.graphs.scc,
    owner.observations.scc.chopping.unfringed,owner.graphs.llowfs)

function shared_prefix_stage_dm!(owner, dm)
    for graph in shared_prefix_branches(owner)
        copyto!(graph_input(graph,:dm_opd),dm)
    end
    return nothing
end

shared_prefix_exposure_state(exposure) = (exposure.start_ns,exposure.covered_ns,
    exposure.next_ns,copy(exposure.integrated),copy(exposure.image),exposure.sequence,exposure.incomplete)

shared_prefix_feed!(owner,args...;kwargs...) =
    step_external_pupil!((name,image,metadata)->nothing,owner,args...;kwargs...)

@testset "Shared prefix retains terminal rectangular pixel integration" begin
    for T in (Float32,Float64)
        shared = shared_prefix_test_owner(;T,pixel_integration=true)
        independent = shared_prefix_test_owner(;T,pixel_integration=true,shared_prefix=false)
        joint = AlgorithmGraphs.optical_execution_graph(shared,:llowfs)
        @test joint !== shared.graphs.llowfs
        @test joint === AlgorithmGraphs.optical_execution_graph(shared,:scc)
        @test length(joint.nodes) == 2
        @test size(graph_output(joint,:path_3)) == (10,12)
        products = (NamedTuple[],NamedTuple[])
        for (ordinal,owner) in enumerate((shared,independent))
            for sequence in 1:3
                shared_prefix_stage_dm!(owner,fill(T(sequence*1e-8),16,16))
                @test step_external_pupil!((name,image,metadata)->
                    push!(products[ordinal],(;name,image=copy(image),metadata)),owner,
                    fill(T(sequence*1e-9),16,16),fill(T(0.7),16,16),
                    ModelTimestamp(2sequence-2),ModelDuration(2);
                    domain=ntuple(_->UInt8(7),16),generation=UInt64(1),sequence=UInt64(sequence))
            end
        end
        for (actual,reference) in zip(products...)
            @test actual.metadata == reference.metadata
            @test actual.image ≈ reference.image rtol=30eps(T)
        end
        actual_raw = shared.graphs.llowfs.nodes.optics.outputs.output
        reference_raw = independent.graphs.llowfs.nodes.optics.outputs.output
        @test actual_raw ≈ reference_raw rtol=30eps(T)
        @test graph_output(shared.graphs.llowfs,:intensity) ≈
            graph_output(independent.graphs.llowfs,:intensity) rtol=30eps(T)
    end
end

@testset "Suffix OPD inputs remain branch-specific" begin
    for llowfs_role in (:branch_opd,:other_opd)
        shared = shared_prefix_test_owner(;suffix_opd=true,llowfs_keywords=(suffix_role=llowfs_role,))
        independent = shared_prefix_test_owner(;suffix_opd=true,shared_prefix=false,
            llowfs_keywords=(suffix_role=llowfs_role,))
        @test shared.propagation !== nothing
        for owner in (shared,independent)
            for (index,graph) in enumerate(shared_prefix_branches(owner))
                role = index == 3 ? llowfs_role : :branch_opd
                graph_input(graph,role) .= [index*1e-8sin(i+j) for i in 1:16,j in 1:16]
            end
            @test shared_prefix_feed!(owner,zeros(16,16),ones(16,16),ModelTimestamp(0),ModelDuration(2);
                domain=ntuple(_->UInt8(5),16),generation=UInt64(1),sequence=UInt64(1))
        end
        for name in (:scc,:llowfs)
            @test graph_output(getproperty(shared.graphs,name),:intensity) ≈
                graph_output(getproperty(independent.graphs,name),:intensity) rtol=3e-13
        end
    end
end

@testset "Shared scalar prefix preserves F/U and independent LL exposures" begin
    for T in (Float32,Float64)
        shared = shared_prefix_test_owner(;T)
        independent = shared_prefix_test_owner(;T,shared_prefix=false)
        execution = AlgorithmGraphs.optical_execution_graph(shared,:scc)
        @test execution === AlgorithmGraphs.optical_execution_graph(shared,:llowfs)
        @test execution === AlgorithmGraphs.optical_execution_graph(shared,:scc;unfringed=true)
        @test independent.propagation === nothing
        @test AlgorithmGraphs.optical_execution_graph(independent,:scc) === independent.graphs.scc
        @test AlgorithmGraphs.optical_execution_graph(independent,:scc;unfringed=true) ===
            independent.observations.scc.chopping.unfringed
        @test_throws ArgumentError AlgorithmGraphs.optical_execution_graph(shared,:llowfs;unfringed=true)
        @test length(only(values(execution.nodes)).owner.prefix.train.plan.steps) == 6
        @test graph_execution_policy(execution) isa StreamGraphExecution
        @test captured_graph_node_count(execution) == 0
        products = (NamedTuple[],NamedTuple[])
        domain = ntuple(_->UInt8(3),16)
        schedule = ((0,1,1,0),(1,2,2,0),(3,2,3,0),(5,1,4,0),
            (6,2,5,0),(10,2,7,1),(12,2,8,1))
        for (index,(time,span,sequence,lost)) in enumerate(schedule)
            opd = [T(index*1e-8sin(i+2j)) for i in 1:16,j in 1:16]
            amplitude = [T(0.6+0.2sin(index+i-j)) for i in 1:16,j in 1:16]
            dm = [T(index*2e-8cos(2i+j)) for i in 1:16,j in 1:16]
            for (ordinal,owner) in enumerate((shared,independent))
                shared_prefix_stage_dm!(owner,dm)
                emit(name,image,metadata) = push!(products[ordinal],(;name,image=copy(image),metadata))
                @test step_external_pupil!(emit,owner,opd,amplitude,
                    ModelTimestamp(time),ModelDuration(span);domain,generation=UInt64(1),
                    sequence=UInt64(sequence),lost_samples=UInt64(lost))
            end
            @test graph_output(shared.graphs.llowfs,:intensity) ≈
                graph_output(independent.graphs.llowfs,:intensity) rtol=30eps(T)
        end
        @test length(products[1]) == length(products[2]) > 4
        for (actual,reference) in zip(products...)
            @test actual.name == reference.name
            @test actual.metadata == reference.metadata
            @test actual.image ≈ reference.image rtol=30eps(T) atol=eps(T)^2
        end
        @test Set(product.metadata.chopped_image_kind for product in products[1]
            if product.name === :scc) == Set((:fringed,:unfringed))
        @test graph_step_sequence(execution) == length(schedule)
        @test all(graph -> graph_step_sequence(graph) == 0,shared_prefix_branches(shared))
        for owner in (shared,independent)
            reset_external_pupil!(owner)
            @test all(exposure -> exposure.sequence == 0 && exposure.covered_ns == 0,
                values(owner.exposures))
            shared_prefix_stage_dm!(owner,zeros(T,16,16))
        end
        @test graph_step_sequence(execution) == 0
        empty!.(products)
        for (ordinal,owner) in enumerate((shared,independent))
            for (sequence,time) in enumerate((0,2,4,6))
                @test step_external_pupil!((name,image,metadata)->
                    push!(products[ordinal],(;name,image=copy(image),metadata)),owner,
                    zeros(T,16,16),ones(T,16,16),ModelTimestamp(time),ModelDuration(2);
                    domain,generation=UInt64(2),sequence=UInt64(sequence))
            end
        end
        for (actual,reference) in zip(products...)
            @test actual.metadata == reference.metadata
            @test actual.image ≈ reference.image rtol=30eps(T)
        end
    end
end

@testset "Unequal prefix inputs retain independent propagation" begin
    shared = shared_prefix_test_owner()
    independent = shared_prefix_test_owner(;shared_prefix=false)
    joint = AlgorithmGraphs.optical_execution_graph(shared,:scc)
    domain = ntuple(_->UInt8(4),16)
    for sequence in 1:3
        for owner in (shared,independent)
            dm = graph_input(owner.graphs.llowfs,:dm_opd)
            dm[3,4] = sequence == 2 ? 1e-8 : 0
            @test shared_prefix_feed!(owner,zeros(16,16),ones(16,16),
                ModelTimestamp(sequence-1),ModelDuration(1);domain,
                generation=UInt64(1),sequence=UInt64(sequence))
        end
        expected = sequence == 2 ? shared.graphs.scc : joint
        @test AlgorithmGraphs.optical_execution_graph(shared,:scc) === expected
        @test graph_output(shared.graphs.llowfs,:intensity) ≈
            graph_output(independent.graphs.llowfs,:intensity) rtol=3e-13
    end
    @test graph_step_sequence(joint) == 2
    @test graph_step_sequence(shared.graphs.llowfs) == 1
    reset_external_pupil!(shared)
    @test AlgorithmGraphs.optical_execution_graph(shared,:scc) === joint
end

@testset "Direct shared graph use preserves fail-stop admission" begin
    owner = shared_prefix_test_owner()
    graph = AlgorithmGraphs.optical_execution_graph(owner,:scc)
    domain = ntuple(_->UInt8(4),16)
    feed(time,sequence) = shared_prefix_feed!(owner,zeros(16,16),ones(16,16),
        ModelTimestamp(time),ModelDuration(1);domain,generation=UInt64(1),sequence=UInt64(sequence))
    @test feed(0,1)
    exposure = shared_prefix_exposure_state(owner.exposures.scc)
    outputs = map(branch->copy(graph_output(branch,:intensity)),shared_prefix_branches(owner))
    graph_input(owner.graphs.llowfs,:dm_opd)[3,4] = 1e-8
    @test_throws AlgorithmGraphError step_graph!(graph)
    @test graph_failed(graph)
    @test isequal(shared_prefix_exposure_state(owner.exposures.scc),exposure)
    @test all(zip(outputs,shared_prefix_branches(owner))) do (output,branch)
        output == graph_output(branch,:intensity)
    end
    graph_input(owner.graphs.llowfs,:dm_opd)[3,4] = 0
    @test_throws AlgorithmGraphError step_graph!(graph)
    # Unequal prefix inputs would select independent execution, but cannot
    # bypass the failed coordinated graph.
    graph_input(owner.graphs.llowfs,:dm_opd)[3,4] = 1e-8
    @test_throws AlgorithmGraphError feed(1,2)
    graph_input(owner.graphs.llowfs,:dm_opd)[3,4] = 0
    reset_external_pupil!(owner)
    @test !graph_failed(graph)
    @test feed(0,1)
end

@testset "Independent failure and pending work block shared selection" begin
    owner = shared_prefix_test_owner()
    joint = AlgorithmGraphs.optical_execution_graph(owner,:scc)
    original = owner.graphs.llowfs
    domain = ntuple(_->UInt8(6),16)
    feed() = shared_prefix_feed!(owner,zeros(16,16),ones(16,16),ModelTimestamp(0),ModelDuration(1);
        domain,generation=UInt64(1),sequence=UInt64(1))
    # Establish generation before injecting failures: a newer generation
    # intentionally resets graph histories as part of the public contract.
    @test feed()
    wf = only(values(original.nodes)).owner.train.workspace.wavefront
    field = wf.field
    try
        wf.field = similar(field,0,0)
        @test_throws ArgumentError step_graph!(original)
    finally
        wf.field = field
    end
    @test graph_failed(original)
    @test_throws AlgorithmGraphError shared_prefix_feed!(owner,zeros(16,16),ones(16,16),
        ModelTimestamp(1),ModelDuration(1);domain,generation=UInt64(1),sequence=UInt64(2))
    @test graph_step_sequence(joint) == 1
    reset_external_pupil!(owner)
    @test feed()
    ticket = step_graph_async!(original)
    @test_throws AlgorithmGraphError shared_prefix_feed!(owner,zeros(16,16),ones(16,16),
        ModelTimestamp(1),ModelDuration(1);domain,generation=UInt64(1),sequence=UInt64(2))
    wait_graph_step!(ticket)
    reset_external_pupil!(owner)
    @test feed()
end

@testset "Publication completion fails before graph sequence commits" begin
    owner = shared_prefix_test_owner()
    graph = AlgorithmGraphs.optical_execution_graph(owner,:scc)
    ticket = step_graph_async!(graph)
    @test_throws ErrorException AlgorithmGraphs._wait_graph_step!(ticket,
        ()->error("injected product publication failure"))
    @test graph_failed(graph)
    @test graph_step_sequence(graph) == 0
    @test !graph_step_pending(graph)
    @test_throws AlgorithmGraphError step_graph!(graph)
    reset_external_pupil!(owner)
    # External reset clears its retained host pupil (the first branch on CPU).
    # A direct graph caller must stage every bound branch input coherently.
    for branch in shared_prefix_branches(owner)
        fill!(graph_input(branch,:pupil_opd),0)
        fill!(graph_input(branch,:pupil_amplitude),1)
    end
    step_graph!(graph)
    @test graph_step_sequence(graph) == 1
end

@testset "Incompatible scalar prefixes retain independent execution" begin
    aliased = shared_prefix_test_owner(;alias_dm=true)
    @test graph_input(aliased.graphs.llowfs,:dm_opd) === graph_output(aliased.graphs.scc,:intensity)
    @test aliased.propagation === nothing
    for settings in ((wavelength_um=1.6,),(entrance_radius_m=0.005,))
        owner = shared_prefix_test_owner(;llowfs_keywords=settings)
        @test owner.propagation === nothing
        @test AlgorithmGraphs.optical_execution_graph(owner,:scc) === owner.graphs.scc
    end
    # A changed common sampled phase is a different prefix. It may share only
    # the exactly identical earlier operations, never the differing phase mask.
    owner = shared_prefix_test_owner(;llowfs_keywords=(sampled_phase=0.4,))
    joint = AlgorithmGraphs.optical_execution_graph(owner,:scc)
    @test length(only(values(joint.nodes)).owner.prefix.train.plan.steps) == 4
end

@testset "Runtime-input aliases preserve independent pupil-staging order" begin
    pupil = [1e-7sin(0.4i+0.8j) for i in 1:16,j in 1:16]
    for suffix_opd in (false,true)
        b = shared_prefix_test_graph(:llowfs;suffix_opd)
        aliased = graph_input(b,:pupil_opd)
        a = suffix_opd ? shared_prefix_test_graph(:fringed;suffix_opd,suffix_values=aliased) :
            shared_prefix_test_graph(:fringed;dm_values=aliased)
        owner = prepare_external_pupil((;a,b);
            exposure_durations=(a=ModelDuration(2),b=ModelDuration(2)),
            observations=(a=(output=:intensity,sampling_node=:optics),
                b=(output=:intensity,sampling_node=:optics)))
        @test owner.propagation === nothing
        reference = shared_prefix_test_graph(:fringed;suffix_opd)
        copyto!(graph_input(reference,:pupil_opd),pupil)
        step_graph!(reference)
        @test shared_prefix_feed!(owner,pupil,ones(16,16),ModelTimestamp(0),ModelDuration(2);
            domain=ntuple(_->UInt8(8),16),generation=UInt64(1),sequence=UInt64(1))
        @test graph_output(a,:intensity) ≈ graph_output(reference,:intensity) rtol=3e-13
        # Staging b before stepping a would feed its current pupil into the
        # aliased runtime screen, instead of the previous zero-valued pupil.
        role = suffix_opd ? :branch_opd : :dm_opd
        copyto!(graph_input(reference,role),pupil)
        step_graph!(reference)
        @test !isapprox(graph_output(a,:intensity),graph_output(reference,:intensity);rtol=1e-5)
    end
end

@testset "Distinct native pupil aliases retain independent staging" begin
    backing = zeros(256)
    first_pupil, second_pupil = reshape(backing,16,16), reshape(backing,16,16)
    @test first_pupil !== second_pupil
    @test Base.mightalias(first_pupil,second_pupil)
    a = shared_prefix_test_graph(:fringed;opd_values=first_pupil)
    b = shared_prefix_test_graph(:llowfs;opd_values=second_pupil)
    owner = prepare_external_pupil((;a,b);
        exposure_durations=(a=ModelDuration(2),b=ModelDuration(2)),
        observations=(a=(output=:intensity,sampling_node=:optics),
            b=(output=:intensity,sampling_node=:optics)))
    @test owner.propagation === nothing
end
