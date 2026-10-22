# GPI_ROOT OUTPUT_DIRECTORY; numerical inputs/RTC remain owned by GPI's harness.
using TOML

function fine_fft_gpi_preflight(n, points)
    target = AdaptiveOpticsSim.Backends.compute_device(CUDA.zeros(Float32,1))
    system = GPI2Sim.prepare_calibration_system(;
        parameters=GPI2Sim.GPI2Parameters(pupil_resolution=n,modulation_points=points,noise=false),
        target, execution=AdaptiveOpticsSim.AlgorithmGraphs.CapturedGraphExecution())
    batch = AdaptiveOpticsSim.AlgorithmGraphs.prepared_graph_node(system.graph,:pyramid).prepared.workspace.modulation_batch
    return Dict("pupil_resolution"=>n,"modulation_points"=>points,
        "batch_size"=>batch.batch_size,"fft_workspace_type"=>string(typeof(batch.fft_workspace)),"stack_shape"=>collect(size(batch.field_stack)),
        "field_stack_bytes"=>sizeof(eltype(batch.field_stack))*length(batch.field_stack))
end

function run_fine_fft_gpi(output)
    records=Dict{String,Any}()
    for (name,n,points) in (("default",240,24),("fine",720,192))
        records[name]=Base.invokelatest(fine_fft_gpi_preflight,n,points)
        GC.gc(true);CUDA.reclaim()
    end
    open(io->TOML.print(io,records;sorted=true),output*".preflight.toml","w")
    loop_main([output])
end

function fine_fft_gpi_main(args)
    length(args)==2 || error("GPI_ROOT OUTPUT_DIRECTORY")
    root=realpath(args[1]);output=abspath(args[2])
    include(joinpath(root,"benchmarks","gpi_closed_loop_cuda.jl"))
    Base.invokelatest(run_fine_fft_gpi,output)
end

if abspath(PROGRAM_FILE) == @__FILE__
    fine_fft_gpi_main(ARGS)
end
