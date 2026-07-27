using AdaptiveOpticsSim

const BACKEND_NAME = length(ARGS) >= 1 ? lowercase(ARGS[1]) : "cpu"
const EXECUTION_MODE = length(ARGS) >= 2 ? lowercase(ARGS[2]) : "batch"
const SAMPLE_COUNT = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 3
const WARMUP_COUNT = length(ARGS) >= 4 ? parse(Int, ARGS[4]) : 2
const MEASURED_COUNT = length(ARGS) >= 5 ? parse(Int, ARGS[5]) : 5

if BACKEND_NAME == "cuda"
    import CUDA
elseif BACKEND_NAME == "amdgpu"
    import AMDGPU
elseif BACKEND_NAME != "cpu"
    error("backend must be cpu, cuda, or amdgpu")
end

EXECUTION_MODE in ("batch", "independent") ||
    error("execution mode must be batch or independent")
SAMPLE_COUNT >= 1 || error("sample count must be positive")
WARMUP_COUNT >= 1 || error("warmup count must be positive")
MEASURED_COUNT >= 1 || error("measured count must be positive")

function profile_backend(name::AbstractString)
    if name == "cpu"
        return CPUBackend()
    elseif name == "cuda"
        CUDA.functional() || error("CUDA is not functional")
        AdaptiveOpticsSim.disable_scalar_backend!(
            AdaptiveOpticsSim.CUDABackendTag,
        )
        return CUDABackend()
    end
    AMDGPU.functional() || error("AMDGPU is not functional")
    AdaptiveOpticsSim.disable_scalar_backend!(
        AdaptiveOpticsSim.AMDGPUBackendTag,
    )
    return AMDGPUBackend()
end

function profile_sources(count::Int, ::Type{T}) where {T<:AbstractFloat}
    sources = Vector{Source}(undef, count)
    @inbounds for index in 1:count
        fraction = count == 1 ? zero(T) : T(index - 1) / T(count - 1)
        sources[index] = Source(
            band=:custom,
            wavelength=T(650e-9) + fraction * T(300e-9),
            photon_irradiance=one(T) + fraction,
            coordinates=(fraction * T(0.12), fraction * T(210)),
            T=T,
        )
    end
    return sources
end

function synchronize_products!(products)
    first_product = first(products)
    AdaptiveOpticsSim.synchronize_backend!(
        AdaptiveOpticsSim.execution_style(first_product.values),
    )
    return products
end

function run_profile()
    T = Float32
    backend = profile_backend(BACKEND_NAME)
    telescope = Telescope(
        resolution=64,
        diameter=T(8),
        central_obstruction=zero(T),
        T=T,
        backend=backend,
    )
    pupil = PupilFunction(telescope; T=T, backend=backend)
    sources = profile_sources(SAMPLE_COUNT, T)

    execute! = if EXECUTION_MODE == "batch"
        prepared = prepare_direct_imaging_batch(
            fill(pupil, length(sources)),
            sources;
            zero_padding=2,
        )
        () -> synchronize_products!(form_direct_image!(prepared))
    else
        prepared = map(sources) do source
            prepare_direct_imaging(pupil, source; zero_padding=2)
        end
        () -> begin
            @inbounds for index in eachindex(prepared)
                form_direct_image!(prepared[index])
            end
            AdaptiveOpticsSim.synchronize_backend!(
                AdaptiveOpticsSim.execution_style(
                    direct_imaging_output(first(prepared)).values,
                ),
            )
            prepared
        end
    end

    for _ in 1:WARMUP_COUNT
        execute!()
    end
    GC.gc()
    for _ in 1:MEASURED_COUNT
        execute!()
    end

    fft_execute_calls_per_iteration =
        EXECUTION_MODE == "batch" ? 1 : SAMPLE_COUNT
    println("direct_imaging_batch_submission_profile")
    println("  backend: ", BACKEND_NAME)
    println("  execution_mode: ", EXECUTION_MODE)
    println("  sample_count: ", SAMPLE_COUNT)
    println("  pupil_resolution: 64")
    println("  padded_resolution: 128")
    println("  numeric_type: Float32")
    println("  warmup_count: ", WARMUP_COUNT)
    println("  measured_count: ", MEASURED_COUNT)
    println(
        "  fft_execute_calls_per_iteration: ",
        fft_execute_calls_per_iteration,
    )
    println(
        "  expected_total_fft_execute_calls: ",
        (WARMUP_COUNT + MEASURED_COUNT) *
        fft_execute_calls_per_iteration,
    )
    return nothing
end

run_profile()
