# Diagnostic sweep of all exact divisors fitting the existing field-stack budget.
include("profile_physical_pyramid_cuda.jl")

function fine_batch_sweep(output)
    CUDA.allowscalar(false)
    BLAS.set_num_threads(1)
    Backends.set_fft_provider_threads!(1)
    result=Dict{String,Any}("revision"=>readchomp(`git -C $(@__DIR__) rev-parse HEAD`),
        "scope"=>"Fixed fine physical sampling, captured completed rate node, 2 randomized rounds; benchmark overrides only", "cases"=>Any[])
    persist()=open(io->TOML.print(io,result;sorted=true),output,"w")
    for name in ("gpi_fine","copper_fine")
        c=only(filter(c->c.name==name,PHYSICAL_PROFILE_CASES))
        baseline,reference=measure_batch(c,WavefrontSensors.PyramidPupilTiltStrategy(),-1,nothing;warmups=8,samples=16)
        bound=max(1,div(256*1024^2,baseline["padded_side"]^2*sizeof(ComplexF32)))
        candidates=filter(b->c.points%b==0,1:bound)
        record=Dict{String,Any}("name"=>name,"baseline"=>baseline,"trials"=>Any[])
        push!(result["cases"],record);persist()
        for round in 1:2, batch in shuffle(MersenneTwister(872+round),collect(candidates))
            GC.gc(true);CUDA.reclaim()
            trial,_=measure_batch(c,WavefrontSensors.PyramidPupilTiltStrategy(),batch,reference;warmups=8,samples=16)
            trial["round"]=round
            push!(record["trials"],trial);persist()
        end
    end
    result["complete"]=true;persist()
end
if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS)==1 || error("OUTPUT.toml")
    fine_batch_sweep(abspath(only(ARGS)))
end
