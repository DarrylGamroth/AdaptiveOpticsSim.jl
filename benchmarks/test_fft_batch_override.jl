# Keep the historical baseline independent of the current production budget.
include("modulation_batch_override.jl")
using Test
@testset "FFT benchmark baseline and production budget selection" begin
    for (points,side,baseline) in ((192,2880,4),(256,3072,2))
        plane_bytes=side^2*sizeof(ComplexF32)
        MODULATION_BATCH[]=0
        @test ModulationBatchWFS._pyramid_modulation_batch_size(points,plane_bytes,16,1536*1024^2)==baseline
        MODULATION_BATCH[]=-1
        @test ModulationBatchWFS._pyramid_modulation_batch_size(points,plane_bytes,16,1536*1024^2)==16
        MODULATION_BATCH[]=8
        @test ModulationBatchWFS._pyramid_modulation_batch_size(points,plane_bytes,16,1536*1024^2)==8
    end
end
