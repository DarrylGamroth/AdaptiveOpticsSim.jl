# Adapted from SPIDERSSim b479871 test/physical15.jl. The oracle is the
# centred Fresnel equation and does not call the copied numerical executor.
using AdaptiveOpticsSim,Proper,Test,LinearAlgebra,FITSIO,YAML
using AdaptiveOpticsSim.AlgorithmGraphs
const CP=Base.get_extension(AdaptiveOpticsSim,:AdaptiveOpticsSimProperExt)
using FFTW
using SHA

# Redistributable compiler-contract fixture. The oracle uses the defining
# centred Fresnel equations, independently of PROPER and the prepared runtime.
function compiled_paths_synthetic_fixture(; forward=true, from_waist=true,n=8)
    m = 5
    lambda, dx, dz = 1.6e-6, 4e-5, forward ? 0.013 : -0.013
    q = collect(0:n-1) .- n÷2
    entrance = ComplexF64[exp(-0.07*(x*x+y*y))*cis(0.2x+0.03y*y) for y in q, x in q]
    edge = ComplexF64[exp(-0.012*(x*x+y*y)^2) for y in q, x in q]
    mask = ComplexF64[0.65+0.02x+0.01y for y in q, x in q]
    phase = ComplexF64[cis(0.07x*y) for y in -2:2, x in -2:2]
    output_coords = collect(-2:2)
    output_freq = FFTW.fftshift(FFTW.fftfreq(m))
    kernel = exp.(2im*pi .* (output_coords * transpose(output_freq)))
    arrays = Dict{String,Matrix{ComplexF64}}(
        "entrance"=>entrance, "identity"=>Matrix{ComplexF64}(I,n,n),
        "edge"=>edge, "mask"=>mask, "resample"=>kernel,
        "phase"=>phase, "one"=>ones(ComplexF64,m,m),
        "scalar"=>ones(ComplexF64,1,1), "pixel"=>ones(ComplexF64,m+2,m+2))
    mult(name) = Dict("kind"=>"mult","array"=>name)
    screen(name) = Dict("kind"=>"live_screen", "input"=>(name=="DM" ? "dm_opd" : "pupil_opd"),
        "kernel"=>"identity", "size"=>n, "dx"=>dx, "lam"=>lambda)
    prefix = [screen("Subaru_Pupil"), Dict("kind"=>"normalize"), screen("DM"),
        Dict("kind"=>"tf","dx"=>dx,"dz"=>0.004,"lam"=>lambda,"edge"=>"edge"),
        Dict("kind"=>"ft","dx"=>dx,"dz"=>dz,"lam"=>lambda,"edge"=>"edge",
            "from_waist"=>from_waist,"forward"=>forward,
            "gouy_real"=>0.0,"gouy_imag"=>forward ? -1.0 : 1.0)]
    metadata = Dict(
        "schema_version"=>1, "array_orientation"=>"yx", "grid"=>n,
        "image_size"=>m, "crop_start"=>1, "n_env"=>1,
        "wavelengths"=>[lambda], "channels"=>[lambda],
        "propagations"=>[Dict("path_index"=>0,"anchor_index"=>0,
            "entrance"=>"entrance","prefix"=>prefix,
            "closed"=>[mult("mask")],"only"=>[])],
        "fast"=>Dict("resample_kernel"=>"resample","main_phase"=>"phase",
            "basis"=>"scalar","coefficient"=>"scalar","visibility"=>["scalar"],
            "pinhole_gain_real"=>[0.3],"pinhole_gain_imag"=>[-0.2],
            "phantom_pad"=>1,"pixel_pad"=>1,"pixel_mtf"=>"pixel",
            "apply_pixel_mtf"=>false),
        "fast_paths"=>[Dict("source_index"=>0,"shift_transfer"=>"one",
            "phase_multiplier"=>["one"])])
    return (; metadata, arrays, n, m, lambda, dx, dz, forward, from_waist)
end

function compiled_paths_grouped_fixture()
    fixture=compiled_paths_synthetic_fixture()
    prefix=only(fixture.metadata["propagations"])["prefix"]
    for _ in 1:3
        append!(prefix,[Dict("kind"=>"mult","array"=>"mask"),
            deepcopy(prefix[4]),deepcopy(prefix[5])])
    end
    push!(prefix,deepcopy(prefix[3]))
    return fixture
end

function compiled_paths_resample_fixture(;n=8,distance=0.004)
    fixture=compiled_paths_synthetic_fixture(;n)
    prefix=only(fixture.metadata["propagations"])["prefix"]
    transfer=deepcopy(prefix[4])
    transfer["dz"]=distance
    append!(prefix,[transfer,Dict("kind"=>"resample","kernel"=>"identity",
        "roll"=>"edge","post"=>"mask"),deepcopy(prefix[3])])
    return fixture
end

function compare_compiled_resample_fields(owner;synchronize=()->nothing,capture=nothing,replay=nothing)
    w=owner.workspace
    path=only(owner.plan.paths)
    op,tf,resample=path.prefix[5:7]
    wf,ctx=w.wavefront,w.context
    wf.wavelength_m=op.wavelength_m;wf.sampling_m=op.dx
    first=Proper.prepare_propagation(wf,op.dz,ctx;kind=op.from_waist ? :wts : :stw)
    wf.wavelength_m=tf.wavelength_m;wf.sampling_m=tf.dx
    last=Proper.prepare_propagation(wf,tf.dz,ctx;kind=:ptp)
    work=w.operations[7]
    input=similar(w.field)
    T=real(eltype(input));n=size(input,1)
    ordinary()=begin
        AdaptiveOpticsSim.Backends.copyto_backend_async!(w.field,input)
        CP._compiled_paths_proper_pair!(w,op,tf,first,last)
        CP._compiled_paths_resample!(w,resample,work)
    end
    optimized()=begin
        AdaptiveOpticsSim.Backends.copyto_backend_async!(w.field,input)
        CP._compiled_paths_proper_resample!(w,op,tf,resample,work,first,last)
    end
    fill!(input,one(eltype(input)))
    ordinary();optimized();synchronize()
    captures=capture===nothing ? nothing : map(capture,(ordinary,optimized))
    # Native graph arguments retain addresses; explicitly retain all backing
    # arrays while comparing changed-input replays of the complex field.
    GC.@preserve first last input w begin
        for state in 1:3
            copyto!(input,Complex{T}[sin(.3y+.11x+state)+im*cos(.2y-.17x-state) for y in 1:n,x in 1:n])
            captures===nothing ? ordinary() : replay(captures[1])
            synchronize();expected=Array(w.field)
            captures===nothing ? optimized() : replay(captures[2])
            synchronize();actual=Array(w.field)
            gate=T===Float32 ? 1e-5 : 2e-12
            @test norm(actual-expected)/norm(expected)<=gate
            @test maximum(abs,actual-expected)/maximum(abs,expected)<=gate
            last.direction==0 && @test actual==expected
        end
    end
    return nothing
end

function compiled_paths_synthetic_oracle(fixture, dm, entrance_opd, attenuation)
    f = fixture
    fftc(a) = FFTW.fftshift(FFTW.fft(FFTW.ifftshift(a)))
    ifftc(a) = FFTW.fftshift(FFTW.ifft(FFTW.ifftshift(a)))
    field = f.arrays["entrance"] .* exp.(-attenuation) .* cis.(2pi/f.lambda .* entrance_opd)
    field ./= sqrt(sum(abs2,field))
    field .*= cis.(2pi/f.lambda .* dm)
    freq = FFTW.fftshift(FFTW.fftfreq(f.n,inv(f.dx)))
    rho2 = freq.^2 .+ transpose(freq.^2)
    field = ifftc(fftc(field .* f.arrays["edge"]) .* cis.(-pi*f.lambda*0.004 .* rho2))
    field .*= f.arrays["edge"]
    dx_out = f.lambda*abs(f.dz)/(f.n*f.dx)
    function chirp(dx)
        q = (collect(0:f.n-1) .- f.n÷2) .* dx
        cis.(pi/(f.lambda*f.dz) .* (q.^2 .+ transpose(q.^2)))
    end
    f.from_waist && (field .*= chirp(f.dx))
    field = f.forward ? fftc(field)/f.n : ifftc(field)*f.n
    field .*= f.forward ? -im : im
    !f.from_waist && (field .*= chirp(dx_out))
    closed = (field .* f.arrays["mask"])[2:6,2:6] .* f.arrays["phase"]
    only = field[2:6,2:6] .* f.arrays["phase"]
    return (fringed=abs2.(closed .+ (0.3-0.2im).*only),
        unfringed=abs2.(closed), reference=abs2.(only))
end

# Multiple live branches and spectral anchors, plus an independent ghost.
# Branches fork after the original fixture's full Fresnel prefix.
function compiled_paths_shared_fixture(;n=8)
    fixture = compiled_paths_synthetic_fixture(;n)
    metadata, arrays = fixture.metadata, fixture.arrays
    metadata["n_env"] = 3
    metadata["wavelengths"] = [fixture.lambda, 1.1fixture.lambda]
    metadata["channels"] = copy(metadata["wavelengths"])
    original = only(metadata["propagations"])
    paths = typeof(original)[]
    for path in 0:2, anchor in 0:1
        branch = deepcopy(original)
        branch["path_index"], branch["anchor_index"] = path, anchor
        for op in branch["prefix"]
            haskey(op, "lam") && (op["lam"] = metadata["wavelengths"][anchor + 1])
        end
        name = "branch_$path"
        arrays[name] = fixture.arrays["mask"] .* cis.(0.03path .* (1:fixture.n))
        push!(branch["prefix"], Dict("kind"=>"mult", "array"=>name))
        # A live suffix needs its own original operation workspace after the
        # shared operations have been skipped; it also checks phase metadata.
        branch["closed"] = Dict{String,Any}[op for op in branch["closed"]]
        push!(branch["closed"], deepcopy(branch["prefix"][3]))
        push!(branch["closed"], deepcopy(branch["prefix"][4]))
        push!(paths, branch)
    end
    metadata["propagations"] = paths
    arrays["basis"] = ComplexF64[0.8 0.2; 0.3 0.7]
    arrays["coefficients"] = ComplexF64[1 0.2im -0.1 0.03; 0.9 0.1 0.1im 0.02]
    arrays["visibility"] = ones(ComplexF64, 4, 4)
    arrays["ghost"] = 0.01 .* fixture.arrays["phase"]
    fast = metadata["fast"]
    fast["basis"], fast["coefficient"] = "basis", "coefficients"
    fast["visibility"] = ["visibility", "visibility"]
    fast["pinhole_gain_real"] = [0.3, 0.2]
    fast["pinhole_gain_imag"] = [-0.2, 0.1]
    metadata["fast_paths"] = [Dict("source_index"=>source,
        "shift_transfer"=>"one", "phase_multiplier"=>["one", "phase"],
        "ghost_field_anchor"=>source == -1 ? ["ghost", "ghost"] : String[])
        for source in (0, 1, 2, -1)]
    return fixture
end

# Nonseparable complex filters, odd/even image sizes, an independent ghost,
# and defocus on only part of the fast batch expose ordering and slot errors.
function compiled_paths_filtered_fixture(; image_size=5, phantom_pad=2)
    fixture = compiled_paths_shared_fixture()
    metadata, arrays = fixture.metadata, fixture.arrays
    m = image_size
    metadata["image_size"] = m
    metadata["crop_start"] = (fixture.n - m) ÷ 2
    metadata["fast"]["phantom_pad"] = phantom_pad
    q = collect(0:m-1) .- m÷2
    freq = FFTW.fftshift(FFTW.fftfreq(m))
    arrays["resample"] = exp.(2im*pi .* (q * transpose(freq)))
    arrays["phase"] = ComplexF64[cis(0.07x*y + 0.03y) for y in q, x in q]
    arrays["ghost"] = 0.01 .* arrays["phase"]
    arrays["pixel"] = ones(ComplexF64, m+2, m+2)
    for (k, path) in enumerate(metadata["fast_paths"])
        shift = "shift_$k"
        arrays[shift] = ComplexF64[(1+0.02y-0.01x)*cis(0.03k*x-0.05y+0.01x*y) for y in q, x in q]
        path["shift_transfer"] = shift
        phases = ["phase_$(k)_$a" for a in 1:2]
        for (a, name) in enumerate(phases)
            arrays[name] = ComplexF64[(1+0.01a*x)*cis(0.02k*y-0.04a*x+0.03x*y) for y in q, x in q]
        end
        path["phase_multiplier"] = phases
        if k in (1, 3)
            p = m + 2phantom_pad
            f = FFTW.fftshift(FFTW.fftfreq(p))
            transfers = ["defocus_$(k)_$a" for a in 1:2]
            for (a, name) in enumerate(transfers)
                arrays[name] = ComplexF64[cis(0.8a*(fy^2+1.3fx^2)+0.2k*fy*fx) for fy in f, fx in f]
            end
            path["defocus_transfer"] = transfers
        end
    end
    return merge(fixture, (; m))
end

function compiled_paths_fast_oracle(owner, sources, free_ghost)
    plan = owner.plan
    fftc(a) = FFTW.fftshift(FFTW.fft(FFTW.ifftshift(a)))
    ifftc(a) = FFTW.fftshift(FFTW.ifft(FFTW.ifftshift(a)))
    filter(a, transfer) = ifftc(fftc(a) .* transfer)
    return [begin
        if path.source_index == 0
            input = free_ghost ? copy(path.ghost[a]) : zeros(eltype(path.ghost[a]), plan.image_size, plan.image_size)
        else
            input = copy(sources[(path.source_index-1)*plan.n_anchor+a])
            if !isempty(path.defocus_transfer)
                pad = plan.phantom_pad
                padded = zeros(eltype(input), plan.image_size+2pad, plan.image_size+2pad)
                padded[pad+1:pad+plan.image_size, pad+1:pad+plan.image_size] .= input
                input = filter(padded, path.defocus_transfer[a])[pad+1:pad+plan.image_size, pad+1:pad+plan.image_size]
            end
        end
        filter(input, path.shift_transfer) .* path.phase_multiplier[a]
    end for path in plan.fast_paths for a in 1:plan.n_anchor]
end

function compiled_paths_allocated(owner)
    AdaptiveOpticsSim.Optics.propagate_compiled_optical_paths!(owner)
    return @allocated AdaptiveOpticsSim.Optics.propagate_compiled_optical_paths!(owner)
end

function compiled_paths_yaml_fixture(directory;source_mode=false,
        fixture=compiled_paths_synthetic_fixture())
    source_mode || (fixture.arrays["entrance"][2,3]=0) # Frozen opaque sky entrance sample.
    metadata=deepcopy(fixture.metadata)
    metadata["mode"]=source_mode ? "source" : "sky"
    metadata["dx_output"]=0.001
    metadata["live"]=Dict("dm_size"=>fixture.n,"dm_pitch"=>fixture.dx,
        "entrance_size"=>fixture.n,"entrance_pitch"=>fixture.dx)
    if source_mode
        popfirst!(metadata["propagations"][1]["prefix"])
        metadata["live"]["entrance_size"]=0
        metadata["live"]["entrance_pitch"]=0.0
    end
    metadata["arrays"]=Dict{String,Any}()
    for (name,values) in fixture.arrays
        path=joinpath(directory,name)
        header=FITSHeader(["GENERATE","SOURCE"],["compiled_paths.jl test fixture","synthetic Fresnel equations"],["Generated fixture","Regenerate by running tests"])
        FITS(path,"w") do file
            write(file,permutedims(real.(values));header)
            write(file,permutedims(imag.(values));header)
        end
        metadata["arrays"][name]=Dict("path"=>name,"sha256"=>bytes2hex(sha256(read(path))),
            "shape"=>collect(size(values)),"type"=>"ComplexF64","axes"=>[2,1],"real_hdu"=>1,"imag_hdu"=>2)
    end
    path=joinpath(directory,"model.yaml")
    write(path,YAML.write(metadata))
    return (;fixture,metadata,path,sha256=bytes2hex(sha256(read(path))))
end
function compiled_paths_test_graph(model;mode="log_attenuation",reference=nothing,
        dm=zeros(8,8),opd=zeros(8,8),amplitude=zeros(8,8),grid=nothing,
        target=AdaptiveOpticsSim.Backends.compute_device(dm),execution=StreamGraphExecution(),
        output_sampling_calibration=nothing)
    if grid===nothing
        grid=(size=(8,8),axes=(1,2),sampling_m=(model.fixture.dx,model.fixture.dx),center_pixel=(5.0,5.0))
    end
    config=(model_path=model.path,model_sha256=model.sha256,dm_grid=grid,pupil_grid=grid,
        entrance_amplitude_mode=mode)
    if output_sampling_calibration!==nothing
        config=merge(config,(;output_sampling_calibration))
    end
    bindings=NamedTuple()
    if mode=="reference_ratio"
        config=merge(config,(entrance_reference_amplitude="pupil_reference",))
        bindings=(pupil_reference=reference,)
    end
    node=AlgorithmGraphs.compiled_optical_paths_node_factory(T=eltype(dm))(:optics,config,(;),bindings)
    definition=algorithm_graph((node,);name=:compiled_paths,
        inputs=(graph_input(:dm,:optics=>:dm_opd,dm),graph_input(:opd,:optics=>:pupil_opd,opd),
            graph_input(:amplitude,:optics=>:pupil_amplitude,amplitude)),
        outputs=(graph_output(:fringed,:optics=>:fringed),graph_output(:unfringed,:optics=>:unfringed),
            graph_output(:reference,:optics=>:reference)))
    return prepare_algorithm_graph(definition;target,execution)
end
