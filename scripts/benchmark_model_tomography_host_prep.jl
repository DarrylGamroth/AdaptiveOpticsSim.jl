"""
    julia --project=test/cuda -e 'using CUDA; include("scripts/benchmark_model_tomography_host_prep.jl"); compare_host_prep(AdaptiveOpticsSim.Backends.CUDABackendTag)'

Opt-in Float32 cold-preparation comparison for lenslet counts 3 and 8. The
CPU route uses the production AOC-backed model builder, then transfers the
retained tomography product arrays through AOS's GPU build backend. The GPU
route uses the production GPU model builder. Both timed routes synchronize
before returning. Julia host allocation bytes do not measure device memory.

Only retained products (K, gamma, mask, Cxx, Cox, Cnz, recstat) are copied;
temporary GPU covariance/solve intermediates are not part of the CPU route.
"""

include(joinpath(@__DIR__, "profile_model_tomography_cpu_phases.jl"))

const HOST_MIN_AVAILABLE_GIB = 6
const HOST_PREP_WARM_RUNS = 5

function available_host_gib()
    info = read("/proc/meminfo", String)
    found = match(r"(?m)^MemAvailable:\s+(\d+)\s+kB", info)
    isnothing(found) && error("cannot determine MemAvailable from /proc/meminfo")
    return parse(Int, found.captures[1]) / 1024^2
end

function ensure_host_memory!()
    available = available_host_gib()
    available >= HOST_MIN_AVAILABLE_GIB ||
        error("stopping: only $(round(available; digits=2)) GiB host memory available")
    return available
end

function synchronize_model!(model)
    style = AOS.Backends.execution_style(model.reconstructor)
    AOS.Backends.synchronize_backend!(style)
    return model
end

copy_product_to_device(backend, A::AbstractMatrix) =
    AOS.Calibration.materialize_build(backend, A)

# CPU Cnz is Diagonal, whereas the native GPU builder retains a dense matrix.
# Stage that representation explicitly so scalar GPU indexing stays disabled.
copy_product_to_device(backend, A::Diagonal) =
    AOS.Calibration.materialize_build(backend, Matrix(A))

function materialize_model(::Type{B}, cpu_model) where {B<:AOS.Backends.GPUBackendTag}
    backend = AOS.Calibration.GPUArrayBuildBackend(B)
    cpu_ops = cpu_model.operators
    gamma = copy_product_to_device(backend, cpu_ops.gamma)
    mask = copy_product_to_device(backend, cpu_model.grid_mask)
    cxx = copy_product_to_device(backend, cpu_ops.cxx)
    cox = copy_product_to_device(backend, cpu_ops.cox)
    cnz = copy_product_to_device(backend, cpu_ops.cnz)
    recstat = copy_product_to_device(backend, cpu_ops.recstat)
    K = copy_product_to_device(backend, cpu_model.reconstructor)
    operators = AOS.Tomography.TomographyOperators(
        gamma, mask, cxx, cox, cnz, recstat, cpu_ops.wavefront_to_meter)
    model = AOS.Tomography.TomographicReconstructor(
        cpu_model.method, K, mask, cpu_model.atmosphere,
        cpu_model.asterism, cpu_model.wfs, cpu_model.tomography,
        cpu_model.dm, cpu_model.fitting, operators)
    return synchronize_model!(model)
end

function build_gpu_model(::Type{B}, case) where {B<:AOS.Backends.GPUBackendTag}
    backend = AOS.Calibration.GPUArrayBuildBackend(B)
    model = build_reconstructor(
        ModelBasedTomography(), case.atmosphere, case.asterism, case.wfs,
        case.tomography, case.dm;
        noise_model=case.noise_model, build_backend=backend)
    return synchronize_model!(model)
end

function measure_routes(::Type{B}, case) where {B<:AOS.Backends.GPUBackendTag}
    ensure_host_memory!()
    t0 = time_ns()
    cpu = @timed build_model(case)
    ensure_host_memory!()
    transfer = @timed materialize_model(B, cpu.value)
    cpu_to_gpu_ms = (time_ns() - t0) / 1e6
    ensure_host_memory!()
    gpu = @timed build_gpu_model(B, case)
    return (; cpu, transfer, gpu, cpu_to_gpu_ms)
end

function product_arrays(model)
    return (
        K=model.reconstructor,
        gamma=model.operators.gamma,
        mask=model.grid_mask,
        cxx=model.operators.cxx,
        cox=model.operators.cox,
        cnz=model.operators.cnz,
        recstat=model.operators.recstat,
    )
end

function verify_products(::Type{B}, measured) where {B<:AOS.Backends.GPUBackendTag}
    BackendArray = AOS.Backends.gpu_backend_array_type(B)
    isnothing(BackendArray) && error("GPU backend $(B) is unavailable")
    materialized = product_arrays(measured.transfer.value)
    native = product_arrays(measured.gpu.value)
    cpu = product_arrays(measured.cpu.value)
    for name in propertynames(materialized)
        transferred = getproperty(materialized, name)
        built = getproperty(native, name)
        host = getproperty(cpu, name)
        transferred isa BackendArray || error("$name transfer is not device-resident")
        built isa BackendArray || error("$name native build is not device-resident")
        size(transferred) == size(built) == size(host) ||
            error("$name product shape mismatch")
        println("    ", name, " shape=", size(host),
            " transferred_type=", typeof(transferred),
            " native_type=", typeof(built))
    end
    cpu_K = cpu.K
    transferred_K = Array(materialized.K)
    native_K = Array(native.K)
    transfer_error = norm(transferred_K - cpu_K) / norm(cpu_K)
    native_error = norm(native_K - cpu_K) / norm(cpu_K)
    println("    K_transfer_relative_Frobenius_error=", transfer_error)
    println("    K_native_relative_Frobenius_error=", native_error)
    isfinite(transfer_error) && transfer_error <= 1e-6 ||
        error("CPU K changed unexpectedly during H2D materialization")
    isfinite(native_error) && native_error <= 1e-3 ||
        error("native GPU K disagrees with CPU K beyond 1e-3")
    return nothing
end

function print_result(label, measured)
    cpu = measured.cpu
    transfer = measured.transfer
    gpu = measured.gpu
    println("  ", label,
        " cpu_ms=", round(cpu.time * 1e3; digits=3),
        " h2d_ms=", round(transfer.time * 1e3; digits=3),
        " cpu_plus_h2d_ms=", round(measured.cpu_to_gpu_ms; digits=3),
        " native_gpu_ms=", round(gpu.time * 1e3; digits=3))
    println("    host_bytes: cpu=", cpu.bytes,
        " h2d=", transfer.bytes,
        " cpu_plus_h2d=", cpu.bytes + transfer.bytes,
        " native_gpu=", gpu.bytes,
        " gc_ms: cpu=", round(cpu.gctime * 1e3; digits=3),
        " h2d=", round(transfer.gctime * 1e3; digits=3),
        " native_gpu=", round(gpu.gctime * 1e3; digits=3))
end

function compare_case(::Type{B}, n_lenslets::Int) where {B<:AOS.Backends.GPUBackendTag}
    ENV["AOS_TOMO_PROFILE_LENSLETS"] = string(n_lenslets)
    case = profile_case(Float32)
    println("case: Float32, n_lenslets=", n_lenslets,
        ", 2 LGS, 2×2 fit, 2 layers, 2 DMs")
    println("  available_host_gib_start=", round(ensure_host_memory!(); digits=2))
    first = measure_routes(B, case)
    print_result("first-use", first)
    for i in 1:HOST_PREP_WARM_RUNS
        warmed = measure_routes(B, case)
        print_result("warmed-$i", warmed)
        i == HOST_PREP_WARM_RUNS && verify_products(B, warmed)
    end
    println("  available_host_gib_end=", round(ensure_host_memory!(); digits=2))
    return nothing
end

function compare_host_prep(::Type{B}) where {B<:AOS.Backends.GPUBackendTag}
    AOS.Backends.disable_scalar_backend!(B)
    println("AOS model tomography host-prep comparison")
    println("  Julia=", VERSION,
        " threads=", Threads.nthreads(),
        " BLAS_threads=", BLAS.get_num_threads(),
        " backend=", AOS.Backends.gpu_backend_name(B))
    for n_lenslets in (3, 8)
        compare_case(B, n_lenslets)
    end
    module_name = AOS.Backends.gpu_backend_name(B) === :cuda ? :CUDA : :AMDGPU
    isdefined(Main, module_name) || error("load $module_name before this script")
    println("backend runtime and device information:")
    getfield(Main, module_name).versioninfo()
    return nothing
end
