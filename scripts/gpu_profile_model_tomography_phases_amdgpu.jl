using AdaptiveOpticsSim
using LinearAlgebra
using SparseArrays

try
    using AMDGPU
catch err
    error("gpu_profile_model_tomography_phases_amdgpu.jl requires AMDGPU.jl: $(sprint(showerror, err))")
end

include(joinpath(@__DIR__, "gpu_profile_model_tomography_phases_contract.jl"))
warm_runs = parse(Int, get(ENV, "AOS_TOMO_PROFILE_WARMED_RUNS", "3"))
warm_runs >= 3 || error("AOS_TOMO_PROFILE_WARMED_RUNS must be at least 3")
for run_label in ("first-use", ("warmed-$i" for i in 1:warm_runs)...)
    measured = @timed run_gpu_model_tomography_phase_profile(
        AdaptiveOpticsSim.Backends.AMDGPUBackendTag; run_label)
    println("  run_host_alloc_bytes: ", measured.bytes,
        "; run_gc_ms: ", round(measured.gctime * 1e3; digits=3))
end
