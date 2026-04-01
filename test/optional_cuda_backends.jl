@testset "Optional CUDA smoke" begin
    run_optional_backend_smoke(CUDABackendTag)
end
