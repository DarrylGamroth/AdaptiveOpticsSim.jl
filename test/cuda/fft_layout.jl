# Force the internal paired family on small fixtures; policy selection is tested
# separately without these preparation overrides in modulation_batch_policy.jl.
using AdaptiveOpticsSim, Test
using AdaptiveOpticsSim.Backends
import CUDA
const PairWFS = AdaptiveOpticsSim.WavefrontSensors
for T in (ComplexF32,ComplexF64)
    @eval PairWFS function _prepare_pyramid_modulation_fft(
        ::$(Backends.AcceleratorStyle{CUDA.CUDABackend}),
        stack::$(CUDA.CuArray{T,3}), points::Int)
        return _prepare_pyramid_transpose_fft(stack)
    end
end
include("modulation_storage.jl")

@testset "Paired FFT workspace ownership and resize" begin
    stack=CUDA.zeros(ComplexF32,17,17,2)
    forward,backward,workspace=PairWFS._prepare_pyramid_transpose_fft(stack)
    batch=PairWFS.PyramidModulationBatchWorkspace(stack,CUDA.zeros(ComplexF32,5,5),
        CUDA.ones(Float32,2),forward,backward,2,workspace)
    binding=PairWFS._pyramid_modulation_batch_workspace_binding(batch)
    resized=PairWFS._resize_pyramid_modulation_batch(batch,CUDA.zeros(ComplexF32,35,35))
    @test typeof(resized) === typeof(batch)
    @test size(resized.fft_workspace.field_stack)==(35,35,2)
    @test resized.pupil_phase === batch.pupil_phase
    @test resized.operating_weights === batch.operating_weights
    @test resized.fft_workspace.field_stack !== workspace.field_stack
    @test PairWFS._pyramid_modulation_batch_workspace_binding(resized) !== binding
    target=compute_device(stack)
    @test PairWFS._require_exact_pyramid_modulation_fft_target(workspace,target) === nothing
    wrong=PairWFS.PyramidTransposeFFTWorkspace(zeros(ComplexF32,17,17,2),
        workspace.fft_plan,workspace.bfft_plan)
    @test_throws PairWFS.WFSPreparationError PairWFS._require_exact_pyramid_modulation_fft_target(wrong,target)
end
