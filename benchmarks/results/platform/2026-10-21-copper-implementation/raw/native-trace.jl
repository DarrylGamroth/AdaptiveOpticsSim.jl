# Usage: julia --project=ISOLATED_ENV native-trace.jl AOS_ROOT COPPER_ROOT
include(joinpath(ARGS[1], "benchmarks", "profile_copper_cuda.jl"))
function trace_native(args)
    @assert realpath(dirname(dirname(pathof(AdaptiveOpticsSim)))) == realpath(args[1])
    @assert realpath(dirname(dirname(pathof(REVOLTCopperSim)))) == realpath(args[2])
    CUDA.allowscalar(false)
    BLAS.set_num_threads(1)
    BACKENDS.set_fft_provider_threads!(1)
    CUDA.versioninfo()
    target=BACKENDS.compute_device(CUDA.zeros(Float32,1))
    for (label,execution) in (("stream",StreamGraphExecution()),("captured",CapturedGraphExecution()))
        system=prepare_hil_system(;target,execution)
        for _ in 1:32
            step_graph!(system.graph)
        end
        println("Native complete-frame timeline: ",label)
        profile=CUDA.@profile trace=true begin
            for _ in 1:4
                step_graph!(system.graph)
            end
        end
        show(stdout,profile);println()
        system=nothing;GC.gc(true);CUDA.reclaim()
    end
end
trace_native(ARGS)
