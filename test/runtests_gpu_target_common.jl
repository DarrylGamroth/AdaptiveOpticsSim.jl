using Test
using AdaptiveOpticsSim
using AdaptiveOpticsSim.Optics
import AdaptiveOpticsSim.Optics: filter!
using AdaptiveOpticsSim.Backends
using AdaptiveOpticsSim.Detectors
using AdaptiveOpticsSim.Atmospheres
using AdaptiveOpticsSim.WavefrontSensors
using AdaptiveOpticsSim.Calibration
using AdaptiveOpticsSim.Control
using AdaptiveOpticsSim.Tomography
using AdaptiveOpticsSim.Ensembles
using AdaptiveOpticsSim: Plant
using AdaptiveOpticsSim.Plant
using LinearAlgebra
using Random
using Statistics

# Hardware targets exercise qualified-public and internal plant contracts in
# addition to the routine exported workflow.
for name in names(Backends; all=true)
    s = String(name)
    if Base.isidentifier(s) && !startswith(s, "#") &&
            !isdefined(@__MODULE__, name)
        @eval const $(name) = getfield(Backends, $(QuoteNode(name)))
    end
end

for name in names(Detectors; all=true)
    s = String(name)
    if Base.isidentifier(s) && !startswith(s, "#") &&
            !isdefined(@__MODULE__, name)
        @eval const $(name) = getfield(Detectors, $(QuoteNode(name)))
    end
end

for name in names(Atmospheres; all=true)
    s = String(name)
    if Base.isidentifier(s) && !startswith(s, "#") &&
            !isdefined(@__MODULE__, name)
        @eval const $(name) = getfield(Atmospheres, $(QuoteNode(name)))
    end
end

for name in names(WavefrontSensors; all=true)
    s = String(name)
    if Base.isidentifier(s) && !startswith(s, "#") &&
            !isdefined(@__MODULE__, name)
        @eval const $(name) =
            getfield(WavefrontSensors, $(QuoteNode(name)))
    end
end

for name in names(Plant; all=true)
    s = String(name)
    if Base.isidentifier(s) && !startswith(s, "#") && !isdefined(@__MODULE__, name)
        @eval const $(name) = getfield(Plant, $(QuoteNode(name)))
    end
end

include("plant_device_batching_fixtures.jl")
include("plant_device_model_matrix_fixtures.jl")
include("backend_optional_common.jl")
include(normpath(joinpath(@__DIR__, "..", "benchmarks", "support", "revolt_like_hil_common.jl")))
include(normpath(joinpath(@__DIR__, "..", "scripts", "gpu_builder_contract.jl")))

BLAS.set_num_threads(1)
Backends.set_fft_provider_threads!(1)

backend_target_branch_mode(::Type{Backends.CUDABackendTag}) =
    BackendStreamExecution()
backend_target_branch_mode(::Type{Backends.AMDGPUBackendTag}) =
    SequentialExecution()

function require_backend_target!(::Type{B}) where {B<:Backends.GPUBackendTag}
    pkg = backend_package_name(B)
    pkg_path = Base.find_package(pkg)
    pkg_path === nothing && error("$(backend_label(B)) target requires $(pkg).jl in the active environment")
    import_backend_package!(B)
    backend_functional(B) || error("$(backend_label(B)) target requires a functional backend/device on this host")
    Backends.disable_scalar_backend!(B)
    return nothing
end

function run_gpu_backend_target(::Type{B}) where {B<:Backends.GPUBackendTag}
    require_backend_target!(B)
    @testset "$(backend_label(B)) hardware target" begin
        run_optional_backend_smoke(B)
        run_gpu_builder_smoke(B)
        run_revolt_like_hil_backend_smoke(B)
    end
    return nothing
end

"""
    run_gpu_detector_target(::Type{<:Backends.GPUBackendTag})

Run the maintained detector qualification surface on a required hardware
backend. This target deliberately excludes unrelated optics, WFS, control, and
orchestration checks so a failure in another subsystem cannot hide detector
qualification evidence. The complete hardware target remains the release
integration gate.
"""
function run_gpu_detector_target(::Type{B}) where {B<:Backends.GPUBackendTag}
    require_backend_target!(B)
    BackendArray = Backends.gpu_backend_array_type(B)
    @test BackendArray !== nothing

    @testset "$(backend_label(B)) detector hardware target" begin
        run_optional_detector_device_model_matrix_checks(B, BackendArray)
        run_optional_cmos_family_checks(B, BackendArray)
        run_optional_shared_detector_ipc_checks(B, BackendArray)
        run_optional_detector_event_checks(B, BackendArray)
        run_optional_avalanche_detector_parity(B, BackendArray)
        run_optional_skipper_ccd_checks(B, BackendArray)
        run_optional_ingaas_checks(B, BackendArray)
        run_optional_spad_qualification_checks(B, BackendArray)
        run_optional_counting_detector_parity(B, BackendArray)
    end
    return nothing
end

function run_revolt_like_hil_backend_smoke(
    ::Type{B},
) where {B<:Backends.GPUBackendTag}
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
