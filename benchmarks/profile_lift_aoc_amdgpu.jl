using AMDGPU

include("profile_lift_aoc_common.jl")

AMDGPU.functional() || error("AMDGPU is not functional on this host")
AMDGPU.allowscalar(false)

function amdgpu_device_allocated_bytes(f)
    before = AMDGPU.alloc_stats.alloc_bytes
    f()
    return AMDGPU.alloc_stats.alloc_bytes - before
end

function amdgpu_profile_region(f)
    status = ccall((:roctxProfilerResume, "librocprofiler-sdk-roctx.so"),
        Cint, (UInt64,), 0)
    iszero(status) || error("rocprofiler rejected LiFT region resume")
    try
        f()
    finally
        AMDGPU.synchronize()
        status = ccall((:roctxProfilerPause,
            "librocprofiler-sdk-roctx.so"), Cint, (UInt64,), 0)
        iszero(status) || error("rocprofiler rejected LiFT region pause")
    end
    return nothing
end

profile_lift_adapter("amdgpu", AMDGPU.ROCArray, Backends.AMDGPUBackend(),
    AMDGPU.synchronize, amdgpu_device_allocated_bytes, AMDGPU.used_memory,
    amdgpu_profile_region)
