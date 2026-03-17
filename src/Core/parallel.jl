using LinearAlgebra
using Logging

struct ParallelConfig
    threads::Int
    fft_threads::Int
    blas_threads::Int
end

function with_parallel_config(cfg::ParallelConfig, f::Function)
    if Threads.nthreads() != cfg.threads
        @warn "start Julia with JULIA_NUM_THREADS=$(cfg.threads)" current=Threads.nthreads()
    end
    set_fft_provider_threads!(cfg.fft_threads)
    BLAS.set_num_threads(cfg.blas_threads)
    return f()
end
