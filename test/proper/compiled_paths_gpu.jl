# Shared CPU-oracle, streamed and captured GPU checks; include from a backend selector.
function compare_compiled_input_capture(model,::Type{T};mode="reference_ratio",grid=nothing,same_backend_gate=nothing,resample_fields=false) where {T}
    grid===nothing && (grid=(size=(8,8),axes=(1,2),
        sampling_m=(model.fixture.dx,model.fixture.dx),center_pixel=(5.0,5.0)))
    reference=ones(T,grid.size)
    reference[1,1]=0
    dm=zeros(T,grid.size);opd=similar(dm);fill!(opd,0)
    amplitude=mode=="reference_ratio" ? copy(reference) : zeros(T,grid.size)
    device=map(COMPILED_ARRAY,(dm,opd,amplitude))
    kwargs=(;mode,reference,grid)
    cpu=compiled_paths_test_graph(model;dm,opd,amplitude,kwargs...)
    streamed=compiled_paths_test_graph(model;dm=device[1],opd=device[2],amplitude=device[3],
        target=COMPILED_TARGET,kwargs...)
    captured=compiled_paths_test_graph(model;dm=device[1],opd=device[2],amplitude=device[3],
        target=COMPILED_TARGET,execution=CapturedGraphExecution(),kwargs...)
    owner=AlgorithmGraphs.prepared_graph_node(captured,:optics)
    cpu_owner=AlgorithmGraphs.prepared_graph_node(cpu,:optics)
    stream_owner=AlgorithmGraphs.prepared_graph_node(streamed,:optics)
    @test captured_graph_node_count(captured)==1
    @test graph_execution_policy(captured) isa CapturedGraphExecution
    @test AlgorithmGraphs.graph_node_capture_capability(owner) isa AlgorithmGraphs.GraphNodeCaptureSafe
    @test Backends.compute_device(owner.propagation.inputs.dm_opd)==COMPILED_TARGET
    @test Backends.compute_device(owner.mapping.dm.row_indices)==COMPILED_TARGET
    @test mode!="reference_ratio" || Backends.compute_device(owner.reference)==COMPILED_TARGET
    retained=(field=owner.propagation.workspace.field,dm=owner.propagation.inputs.dm_opd,
        maps=owner.mapping,reference=owner.workspace.reference,flags=owner.workspace.validation.device)
    # Direct prepared-field capture follows the lowering actually selected for
    # this backend. The complete YAML graph below also tests ordinary backends.
    if resample_fields && owner.propagation.workspace.propagation[5]!==nothing
        Backends._with_prepared_device_execution_context(captured.context) do
            compare_compiled_resample_fields(owner.propagation;
                synchronize=COMPILED_SYNCHRONIZE,
                capture=operation->Backends._capture_prepared_device_graph(operation,captured.context),
                replay=recorded->Backends._launch_prepared_device_graph!(recorded,captured.context))
        end
    end
    first_images=nothing
    previous=nothing
    rtol,atol=T===Float32 ? (8e-4,8e-6) : (8e-10,8e-12)
    for frame in 1:4
        frame==2 && (dm .= T[1e-7sin(.3i+.7j) for i in axes(dm,1),j in axes(dm,2)])
        frame==3 && (opd .= T[1e-7cos(.4i-.2j) for i in axes(opd,1),j in axes(opd,2)])
        if frame==4
            transmission=T[.6+.3sin(.2i+.5j)^2 for i in axes(amplitude,1),j in axes(amplitude,2)]
            amplitude .= mode=="reference_ratio" ? reference .* transmission : -log.(transmission)
        end
        for (destination,source) in zip(device,(dm,opd,amplitude));copyto!(destination,source);end
        COMPILED_SYNCHRONIZE()
        step_graph!(cpu);step_graph!(streamed);step_graph!(captured)
        # Preserve an independent-path oracle for the shared-prefix fixture.
        CP._enqueue_compiled_paths!(cpu_owner.propagation,false,false)
        Backends._with_prepared_device_execution_context(streamed.context) do
            CP._enqueue_compiled_paths!(stream_owner.propagation,false,false,false)
            Backends._synchronize_prepared_device_execution_context!(streamed.context)
        end
        actual=map(name->Array(graph_output(captured,name)),(:fringed,:unfringed,:reference))
        for (name,values) in zip((:fringed,:unfringed,:reference),actual)
            @test values ≈ graph_output(cpu,name) rtol=rtol atol=atol
            @test values ≈ getproperty(cpu_owner.propagation.products,name) rtol=rtol atol=atol
            @test values ≈ Array(getproperty(stream_owner.propagation.products,name)) rtol=rtol atol=atol
            @test values ≈ Array(graph_output(streamed,name)) rtol=rtol atol=atol
            @test all(isfinite,values) && sum(values)>0
            if same_backend_gate!==nothing
                expected=Array(getproperty(stream_owner.propagation.products,name))
                @test norm(values-expected)/norm(expected)<=same_backend_gate
                @test maximum(abs,values-expected)/maximum(abs,expected)<=same_backend_gate
            end
        end
        if frame==1
            first_images=map(copy,actual)
        elseif model.metadata["mode"]=="sky" || frame==2
            @test !isapprox(first(actual),previous;rtol=1e-4,atol=atol)
        end
        previous=copy(first(actual))
        @test owner.propagation.workspace.field===retained.field
        @test owner.propagation.inputs.dm_opd===retained.dm
        @test owner.mapping===retained.maps
        @test owner.workspace.reference===retained.reference
        @test owner.workspace.validation.device===retained.flags
    end
    fill!(dm,0);fill!(opd,0)
    amplitude .= mode=="reference_ratio" ? reference : zero(T)
    for (destination,source) in zip(device,(dm,opd,amplitude));copyto!(destination,source);end
    COMPILED_SYNCHRONIZE()
    reset_graph!(captured)
    for name in (:fringed,:unfringed,:reference)
        @test all(iszero,Array(graph_output(captured,name)))
    end
    @test all(iszero,Array(owner.workspace.validation.device))
    step_graph!(captured)
    for (name,expected) in zip((:fringed,:unfringed,:reference),first_images)
        @test Array(graph_output(captured,name)) ≈ expected rtol=rtol atol=atol
    end
    return nothing
end

@testset "captured GPU compiled optical input grids" begin
    mktempdir() do directory
        model=compiled_paths_yaml_fixture(directory)
        for T in (Float32,Float64), mode in ("reference_ratio","log_attenuation")
            compare_compiled_input_capture(model,T;mode)
            grid=(size=(7,7),axes=(2,1),sampling_m=(2model.fixture.dx,1.5model.fixture.dx),
                center_pixel=(4.25,3.4))
            compare_compiled_input_capture(model,T;mode,grid)
        end
    end
    mktempdir() do directory
        model=compiled_paths_yaml_fixture(directory;source_mode=true)
        compare_compiled_input_capture(model,Float32)
    end
    mktempdir() do directory
        model=compiled_paths_yaml_fixture(directory;fixture=compiled_paths_shared_fixture())
        for T in (Float32,Float64)
            compare_compiled_input_capture(model,T)
        end
    end
    mktempdir() do directory
        model=compiled_paths_yaml_fixture(directory;fixture=compiled_paths_grouped_fixture())
        for T in (Float32,Float64)
            compare_compiled_input_capture(model,T)
        end
    end
    for n in (7,8),distance in (-0.004,0.0,0.004)
        mktempdir() do directory
            model=compiled_paths_yaml_fixture(directory;fixture=compiled_paths_resample_fixture(;n,distance))
            for T in (Float32,Float64)
                compare_compiled_input_capture(model,T;same_backend_gate=1e-5,resample_fields=true)
            end
        end
    end
    for image_size in (5,6)
        mktempdir() do directory
            model=compiled_paths_yaml_fixture(directory;
                fixture=compiled_paths_filtered_fixture(;image_size))
            for T in (Float32,Float64)
                compare_compiled_input_capture(model,T)
            end
        end
    end
end

@testset "compiled capture rejects invalid live input before step publication" begin
    mktempdir() do directory
        model=compiled_paths_yaml_fixture(directory)
        reference=ones(Float32,8,8);reference[1,1]=0;reference[4,4]=floatmin(Float32)
        host=(dm=zeros(Float32,8,8),opd=zeros(Float32,8,8),amplitude=copy(reference))
        device=map(COMPILED_ARRAY,host)
        gpu=compiled_paths_test_graph(model;mode="reference_ratio",reference,
            dm=device.dm,opd=device.opd,amplitude=device.amplitude,
            target=COMPILED_TARGET,execution=CapturedGraphExecution())
        for (role,index,value,message) in (
                (:dm,CartesianIndex(2,3),NaN32,r"DM OPD"),
                (:opd,CartesianIndex(2,3),Inf32,r"entrance OPD"),
                (:amplitude,CartesianIndex(2,3),-.1f0,r"finite and nonnegative"),
                (:amplitude,CartesianIndex(1,1),.2f0,r"reference is opaque"),
                (:amplitude,CartesianIndex(2,3),0f0,r"new opaque entrance"),
                (:amplitude,CartesianIndex(4,4),floatmax(Float32),r"not representable"))
            reset_graph!(gpu)
            step_graph!(gpu)
            saved=getproperty(host,role)[index]
            getproperty(host,role)[index]=value
            copyto!(getproperty(device,role),getproperty(host,role));COMPILED_SYNCHRONIZE()
            ticket=step_graph_async!(gpu)
            @test_throws message wait_graph_step!(ticket)
            @test graph_step_sequence(gpu)==1
            @test graph_failed(gpu)
            @test_throws AlgorithmGraphError step_graph!(gpu)
            getproperty(host,role)[index]=saved
            copyto!(getproperty(device,role),getproperty(host,role));COMPILED_SYNCHRONIZE()
            reset_graph!(gpu);step_graph!(gpu)
            @test graph_step_sequence(gpu)==1
            @test all(isfinite,Array(graph_output(gpu,:fringed)))
        end
    end
end
