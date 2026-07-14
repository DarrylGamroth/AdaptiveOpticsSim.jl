using Test
using AdaptiveOpticsSim
using LinearAlgebra
using Random
include("backend_optional_common.jl")
include(normpath(joinpath(@__DIR__, "..", "benchmarks", "support", "revolt_like_hil_common.jl")))
include(normpath(joinpath(@__DIR__, "..", "scripts", "gpu_builder_contract.jl")))

BLAS.set_num_threads(1)
AdaptiveOpticsSim.set_fft_provider_threads!(1)

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
        run_gpu_builder_smoke(B)
        run_revolt_like_hil_backend_smoke(B)
    end
    return nothing
end

function run_revolt_like_hil_backend_smoke(::Type{B}) where {B<:AdaptiveOpticsSim.GPUBackendTag}
    config_dir = normpath(joinpath(@__DIR__, "..", "benchmarks", "assets", "revolt_like"))
    backend_name = lowercase(backend_label(B))
    cpu_ctx = build_revolt_like_hil_context(;
        backend_name="cpu",
        config_dir,
        sensor=CMOSSensor(T=Float32),
        T=Float32,
        rng=runtime_rng(20260713),
    )
    ctx = build_revolt_like_hil_context(;
        backend_name,
        config_dir,
        sensor=CMOSSensor(T=Float32),
        T=Float32,
        rng=runtime_rng(20260713),
    )
    revolt_like_step!(cpu_ctx)
    revolt_like_step!(ctx)
    @test size(ctx.tiled_frame) == (352, 352)
    @test all(isfinite, Array(ctx.tiled_frame))
    @test sum(Array(ctx.tiled_frame)) > 0
    @test isapprox(Array(ctx.tiled_frame), cpu_ctx.tiled_frame;
        rtol=1f-5, atol=1f-3)
    return nothing
end
