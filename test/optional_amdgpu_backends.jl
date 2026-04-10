@testset "Optional AMDGPU smoke" begin
    run_optional_backend_smoke(AdaptiveOpticsSim.AMDGPUBackendTag)
end
