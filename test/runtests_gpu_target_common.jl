using Test
using AdaptiveOpticsSim
using LinearAlgebra
using Random
include("backend_optional_common.jl")

backend_target_branch_mode(::Type{AdaptiveOpticsSim.CUDABackendTag}) = BackendStreamExecution()
backend_target_branch_mode(::Type{AdaptiveOpticsSim.AMDGPUBackendTag}) = SequentialExecution()

function require_backend_target!(::Type{B}) where {B<:AdaptiveOpticsSim.GPUBackendTag}
    pkg = backend_package_name(B)
    pkg_path = Base.find_package(pkg)
    pkg_path === nothing && error("$(backend_label(B)) target requires $(pkg).jl in the active environment")
    import_backend_package!(B)
    backend_functional(B) || error("$(backend_label(B)) target requires a functional backend/device on this host")
    AdaptiveOpticsSim.disable_scalar_backend!(B)
    return nothing
end

function run_gpu_backend_target(::Type{B}) where {B<:AdaptiveOpticsSim.GPUBackendTag}
    require_backend_target!(B)
    @testset "$(backend_label(B)) hardware target" begin
        run_optional_backend_smoke(B)
        run_equivalence = getfield(Main, :run_gpu_runtime_equivalence)
        Base.invokelatest(run_equivalence, B; branch_mode=backend_target_branch_mode(B))
        run_equivalence_high_accuracy = getfield(Main, :run_gpu_runtime_equivalence_high_accuracy)
        Base.invokelatest(run_equivalence_high_accuracy, B; branch_mode=backend_target_branch_mode(B))
    end
    return nothing
end
