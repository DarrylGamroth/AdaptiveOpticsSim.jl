# Two roots exercise cache-slot recycling. Each root has three branches:
# the later branches reuse different closed/only checkpoints of the first.
function compiled_paths_branch_fixture(;n=8,sign=1,dependent=false,zero_only=false)
    fixture=compiled_paths_shared_fixture(;n)
    metadata,arrays=fixture.metadata,fixture.arrays
    arrays["only_mask"]=arrays["mask"].*ComplexF64[cis(.03i+.07j) for i in 1:n,j in 1:n]
    zero_only && fill!(arrays["only_mask"],0)
    for branch in metadata["propagations"]
        pop!(branch["prefix"])
        wavelength=metadata["wavelengths"][branch["anchor_index"]+1]
        mult(name)=Dict("kind"=>"mult","array"=>name)
        ft(waist,d)=Dict("kind"=>"ft","dx"=>fixture.dx,"dz"=>d,"lam"=>wavelength,
            "edge"=>"edge","from_waist"=>waist,"forward"=>d>=0,
            "gouy_real"=>0.13,"gouy_imag"=>d>=0 ? -0.9 : 0.9)
        tf=Dict("kind"=>"tf","dx"=>fixture.dx,"dz"=>sign*.004,"lam"=>wavelength,"edge"=>"edge")
        closed=[mult("mask"),ft(true,sign*.013),ft(false,-sign*.009),
            mult("mask"),deepcopy(tf),ft(true,sign*.013),
            mult("mask"),ft(false,sign*.009),deepcopy(tf),
            Dict("kind"=>"resample","kernel"=>"identity","roll"=>"edge","post"=>"mask")]
        only=deepcopy(closed)
        only[1]=mult("only_mask")
        index=branch["path_index"]
        if index>0
            position=index==1 ? 4 : 7
            closed[position]=mult("branch_$index")
            only[position]=mult("branch_$index")
        elseif dependent && branch["anchor_index"]==0
            only[1]=mult("mask")
            only[4]=mult("only_mask")
        end
        branch["closed"],branch["only"]=closed,only
    end
    return fixture
end

function force_compiled_branch_batch(owner;image_batch=false)
    backend=AdaptiveOpticsSim.Backends
    w=owner.workspace
    # Force the generic numerical implementation on CPU/AMD as well as the
    # automatically selected CUDA class, without changing the backend policy.
    propagation=CP._compiled_paths_prepare_propagation(backend.CUDABackend(),owner.plan.paths,w.wavefront,w.context)
    pairs=CP._compiled_paths_branch_pairs(owner.plan,propagation)
    batch=CP._compiled_paths_branch_batch(w.field,pairs)
    images=image_batch ? CP._compiled_paths_prepare_image_resample(owner.plan,
        w.closed_sources,w.only_sources,w.mul_coefficients,Val(5)) : w.image_resample
    workspace=CP.CompiledPathsWorkspace(map(fieldnames(typeof(w))) do name
        name===:propagation ? propagation : name===:branch_batch ? batch :
            name===:image_resample ? images : getfield(w,name)
    end...)
    return CP.PreparedCompiledPaths(owner.plan,workspace,owner.inputs,owner.products,owner.execution_context,owner.synthesis,
        CP._compiled_paths_prepare_execution(owner.plan,propagation))
end

function check_compiled_branch_batch(array,target;captured=false)
    backend=AdaptiveOpticsSim.Backends
    context=backend._prepare_device_execution_context(target)
    backend._with_prepared_device_execution_context(context) do
        @testset "Paired optical branches $array" begin
            for T in (Float32,Float64),n in (7,8),sign in (-1,1),dependent in (false,true),zero_only in (false,true)
                dependent && zero_only && continue
                fixture=compiled_paths_branch_fixture(;n,sign,dependent,zero_only)
                dm=array(zeros(T,n,n))
                prepare()=AdaptiveOpticsSim.Optics.prepare_compiled_optical_paths(fixture.metadata,
                    key->fixture.arrays[key];target,T,dm_opd=dm)
                reference=prepare()
                # Cover the combined lowering on two existing cases rather
                # than repeating the full matrix. Twelve sources leave a final
                # partial group of two in the forced five-plane image batch.
                image_batch=(T===Float32 && n==7 && sign==1 && !dependent && !zero_only) ||
                    (T===Float64 && n==8 && sign==-1 && dependent && !zero_only)
                owner=force_compiled_branch_batch(prepare();image_batch)
                w=owner.workspace;batch=w.branch_batch
                if image_batch
                    @test w.image_resample!==nothing
                    @test size(w.image_resample.field,3)==5
                    @test length(w.image_resample.outputs)%size(w.image_resample.field,3)==2
                    @test any(pair->pair!==nothing && any(s->s.closed_cache==0,pair.segments),batch.pairs)
                end
                @test any(!isnothing,batch.pairs)
                if dependent
                    @test batch.pairs[1]===nothing
                else
                    @test all(!isnothing,batch.pairs)
                end
                @test length(w.shared_prefixes)==length(reference.workspace.shared_prefixes)
                # The second wavelength root reuses the same checkpoint slots.
                @test length(w.shared_prefixes)<2sum(length(s.prefix.checkpoints)+length(s.closed.checkpoints)+length(s.only.checkpoints)
                    for s in owner.plan.shared_prefixes[1:3])
                retained=Tuple(w.shared_prefixes)
                operation()=CP._enqueue_compiled_paths!(owner)
                # Preparation uses the owner's stream. Bind managed storage
                # and compile kernels on the test stream before capture.
                operation()
                backend._synchronize_prepared_device_execution_context!(context)
                recorded=captured ? backend._capture_prepared_device_graph(operation,context) : nothing
                first_products=nothing
                for state in (1,2,3,1)
                    copyto!(dm,array(T[state*2e-8sin(.2i+.3j) for i in 1:n,j in 1:n]))
                    # Poison both independent planes: every used value must be
                    # produced by this invocation, including shared consumers.
                    fill!(batch.field,Complex{T}(NaN))
                    if image_batch
                        fill!(w.image_resample.field,Complex{T}(NaN))
                        fill!(w.image_resample.temporary,Complex{T}(NaN))
                    end
                    foreach(field->fill!(field,Complex{T}(NaN)),w.shared_prefixes)
                    foreach(field->fill!(field,Complex{T}(NaN)),reference.workspace.shared_prefixes)
                    CP._enqueue_compiled_paths!(reference,true,true,false)
                    recorded===nothing ? operation() : backend._launch_prepared_device_graph!(recorded,context)
                    backend._synchronize_prepared_device_execution_context!(context)
                    tolerance=T===Float32 ? 1e-5 : 1e-11
                    for name in (:fringed,:unfringed,:reference)
                        a=Array(getproperty(owner.products,name));b=Array(getproperty(reference.products,name))
                        if iszero(norm(b))
                            @test all(iszero,a)
                        else
                            @test norm(a-b)/norm(b)<=tolerance
                            @test maximum(abs,a-b)/maximum(abs,b)<=tolerance
                        end
                    end
                    for (left,right) in ((w.closed_sources,reference.workspace.closed_sources),
                            (w.only_sources,reference.workspace.only_sources),(w.shared_prefixes,reference.workspace.shared_prefixes))
                        for (a,b) in zip(left,right)
                            expected=Array(b)
                            actual=Array(a)
                            @test iszero(norm(expected)) ? all(iszero,actual) : norm(actual-expected)/norm(expected)<=tolerance
                        end
                    end
                    @test all(i->w.shared_prefixes[i]===retained[i],eachindex(retained))
                    if image_batch
                        @test all(isfinite,Array(w.image_resample.field))
                        @test all(isfinite,Array(w.image_resample.temporary))
                    end
                    current=map(name->Array(getproperty(owner.products,name)),(:fringed,:unfringed,:reference))
                    first_products===nothing && (first_products=current)
                    state==1 && @test current==first_products
                end
            end
        end
    end
end

function check_compiled_branch_rejections()
    backend=AdaptiveOpticsSim.Backends
    fixture=compiled_paths_branch_fixture()
    prepare(f)=AdaptiveOpticsSim.Optics.prepare_compiled_optical_paths(f.metadata,key->f.arrays[key])
    owner=force_compiled_branch_batch(prepare(fixture))
    path=first(owner.plan.paths);sharing=first(owner.plan.shared_prefixes)
    propagation=owner.workspace.propagation;offset=length(path.prefix)
    prefix(count,cache,points)=CP.PathPrefixPlan(cache,count,CP._compiled_paths_seal(points))
    sharing_with(left,right)=CP.PathSharingPlan(sharing.prefix,left,right)
    @testset "Unsafe branch pairs retain serial execution" begin
        @test CP._compiled_paths_branch_pair(path,sharing,propagation,offset)!==nothing
        left,right=sharing.closed,sharing.only
        shifted=map(p->CP.PathPrefixCheckpoint(p.index+1,p.cache_index),collect(right.checkpoints))
        @test CP._compiled_paths_branch_pair(path,sharing_with(left,prefix(right.count,right.cache_index,shifted)),propagation,offset)===nothing
        overlap=collect(right.checkpoints)
        overlap[1]=CP.PathPrefixCheckpoint(overlap[1].index,left.checkpoints[1].cache_index)
        @test CP._compiled_paths_branch_pair(path,sharing_with(left,prefix(right.count,right.cache_index,overlap)),propagation,offset)===nothing
        # The only input cannot read a value whose closed producer has not
        # executed yet, even when the remaining operation stencils match.
        lcp=[left.checkpoints[end]];rcp=[right.checkpoints[end]]
        dependent=sharing_with(prefix(3,1,lcp),prefix(3,lcp[1].cache_index,rcp))
        @test CP._compiled_paths_branch_pair(path,dependent,propagation,offset)===nothing
        @test CP._compiled_paths_branch_pair(path,sharing_with(left,prefix(3,right.cache_index,collect(right.checkpoints))),propagation,offset)===nothing
        operations=collect(path.closed)
        old=operations[4]
        operations[4]=CP.PathOperation(CP.PathLiveScreen,old.array,old.auxiliary,old.post,
            old.dx,old.dz,old.wavelength_m,old.gouy,old.from_waist,old.forward,old.source)
        changed=(;closed=operations,only=path.only)
        @test CP._compiled_paths_branch_pair(changed,sharing,propagation,offset)===nothing
        identity=compiled_paths_branch_fixture()
        for p in identity.metadata["propagations"],part in ("closed","only")
            p[part][end-1]["dz"]=0.0
        end
        identity_owner=force_compiled_branch_batch(prepare(identity))
        @test all(isnothing,identity_owner.workspace.branch_batch.pairs)
    end
end

# Later branches fork after one/three transform pairs. Their remaining prefix
# contains ten/six transforms and starts directly from a retained checkpoint.
function check_compiled_input_groups(array,target;captured=false)
    backend=AdaptiveOpticsSim.Backends
    context=backend._prepare_device_execution_context(target)
    backend._with_prepared_device_execution_context(context) do
        @testset "Prepared prefix input groups $array" begin
            for T in (Float32,Float64),n in (7,8),sign in (-1,1)
                fixture=compiled_paths_branch_fixture(;n,sign)
                for branch in fixture.metadata["propagations"]
                    prefix=deepcopy(branch["prefix"][1:3])
                    wavelength=fixture.metadata["wavelengths"][branch["anchor_index"]+1]
                    fork=branch["path_index"]==1 ? 2 : branch["path_index"]==2 ? 4 : 0
                    for pair in 1:6
                        mask=pair==fork ? "branch_$(branch["path_index"])" : "mask"
                        push!(prefix,Dict("kind"=>"mult","array"=>mask))
                        for (waist,distance) in ((true,.013),(false,-.009))
                            push!(prefix,Dict("kind"=>"ft","dx"=>fixture.dx,"dz"=>sign*distance,
                                "lam"=>wavelength,"edge"=>"edge","from_waist"=>waist,
                                "forward"=>sign*distance>=0,"gouy_real"=>.13,"gouy_imag"=>-.9))
                        end
                    end
                    branch["prefix"]=prefix
                end
                dm=array(zeros(T,n,n))
                prepare()=AdaptiveOpticsSim.Optics.prepare_compiled_optical_paths(fixture.metadata,
                    key->fixture.arrays[key];target,T,dm_opd=dm)
                reference=prepare()
                forced=force_compiled_branch_batch(prepare())
                w=forced.workspace
                workspace=CP.CompiledPathsWorkspace(map(fieldnames(typeof(w))) do name
                    name===:branch_batch ? nothing : getfield(w,name)
                end...)
                owner=CP.PreparedCompiledPaths(forced.plan,workspace,forced.inputs,forced.products,
                    forced.execution_context,forced.synthesis,forced.execution)
                retained=Tuple(w.shared_prefixes)
                operation()=CP._enqueue_compiled_paths!(owner)
                operation();backend._synchronize_prepared_device_execution_context!(context)
                recorded=captured ? backend._capture_prepared_device_graph(operation,context) : nothing
                first_products=nothing
                for state in (1,2,1)
                    copyto!(dm,array(T[state*2e-8sin(.2i+.3j) for i in 1:n,j in 1:n]))
                    for work in (w,reference.workspace)
                        fill!(work.field,Complex{T}(NaN))
                        foreach(field->fill!(field,Complex{T}(NaN)),work.shared_prefixes)
                    end
                    CP._enqueue_compiled_paths!(reference,true,true,false)
                    recorded===nothing ? operation() : backend._launch_prepared_device_graph!(recorded,context)
                    backend._synchronize_prepared_device_execution_context!(context)
                    tolerance=T===Float32 ? 1e-5 : 1e-11
                    for name in (:fringed,:unfringed,:reference)
                        actual=Array(getproperty(owner.products,name));expected=Array(getproperty(reference.products,name))
                        @test norm(actual-expected)/norm(expected)<=tolerance
                        @test maximum(abs,actual-expected)/maximum(abs,expected)<=tolerance
                    end
                    for (left,right) in ((w.closed_sources,reference.workspace.closed_sources),
                            (w.only_sources,reference.workspace.only_sources),(w.shared_prefixes,reference.workspace.shared_prefixes))
                        for (a,b) in zip(left,right)
                            expected=Array(b)
                            @test norm(Array(a)-expected)/norm(expected)<=tolerance
                        end
                    end
                    @test all(i->w.shared_prefixes[i]===retained[i],eachindex(retained))
                    products=map(name->Array(getproperty(owner.products,name)),(:fringed,:unfringed,:reference))
                    first_products===nothing && (first_products=products)
                    state==1 && @test products==first_products
                end
            end
        end
    end
end
