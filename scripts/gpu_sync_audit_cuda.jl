using AdaptiveOpticsSim

try
    using CUDA
catch err
    error("gpu_sync_audit_cuda.jl requires CUDA.jl: $(sprint(showerror, err))")
end

include(joinpath(@__DIR__, "gpu_sync_audit_contract.jl"))
run_gpu_sync_audit(CUDABackendTag)
