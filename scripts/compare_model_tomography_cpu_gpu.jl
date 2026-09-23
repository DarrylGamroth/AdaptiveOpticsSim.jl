"""
    julia --project=test/cuda -e 'using CUDA; include("scripts/compare_model_tomography_cpu_gpu.jl"); compare_model_tomography_cpu_gpu(AdaptiveOpticsSim.Backends.CUDABackendTag)'

Compare the CPU and one accelerator's Float32 cold tomography products for
`AOS_TOMO_PROFILE_LENSLETS` (default 3). GPU-to-host copies occur only after
the profiled build, for this diagnostic comparison.
"""

include(joinpath(@__DIR__, "profile_model_tomography_cpu_phases.jl"))
include(joinpath(@__DIR__, "gpu_profile_model_tomography_phases_contract.jl"))

function compare_model_tomography_cpu_gpu(::Type{B}) where {B<:AOS.Backends.GPUBackendTag}
    cpu = build_model(profile_case(Float32))
    gpu = run_gpu_model_tomography_phase_profile(B;
        run_label="comparison", return_products=true)
    cpu_products = (
        cxx=cpu.operators.cxx,
        cox=cpu.operators.cox,
        recstat=cpu.operators.recstat,
        recon=cpu.reconstructor,
    )
    println("CPU/GPU Float32 tomography product comparison")
    println("  n_lenslets: ", cpu.wfs.n_lenslets)
    for name in propertynames(cpu_products)
        cpu_matrix = getproperty(cpu_products, name)
        gpu_matrix = getproperty(gpu, name)
        host_matrix = Array(gpu_matrix)
        size(host_matrix) == size(cpu_matrix) || error("$name shape mismatch")
        relative_error = norm(host_matrix - cpu_matrix) / norm(cpu_matrix)
        tolerance = name === :cxx || name === :cox ? 1e-5 : 1e-3
        println("  ", name, ": shape=", size(cpu_matrix),
            " cpu_type=", typeof(cpu_matrix),
            " gpu_type=", typeof(gpu_matrix),
            " relative_Frobenius_error=", relative_error,
            " tolerance=", tolerance)
        isfinite(relative_error) && relative_error <= tolerance ||
            error("$name CPU/GPU relative Frobenius error exceeds $tolerance")
    end
    return nothing
end
