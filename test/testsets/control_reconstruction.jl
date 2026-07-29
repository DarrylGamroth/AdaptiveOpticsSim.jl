@testset "Control matrices and reconstructors" begin
    D = rand(4, 3)
    control_matrix = ControlMatrix(D)
    @test size(control_matrix.M) == (3, 4)
    @test control_matrix.cond > 0
    @test control_matrix.effective_rank == 3
    truncated_control_matrix = with_truncation(control_matrix, 1)
    @test truncated_control_matrix.n_trunc == 1

    D_sing = [1.0 0.0; 0.0 1e-12]
    exact_control_matrix = ControlMatrix(D_sing; policy=ExactPseudoInverse())
    tsvd_control_matrix = ControlMatrix(D_sing; policy=TSVDInverse(rtol=1e-9))
    @test exact_control_matrix.effective_rank == 2
    @test tsvd_control_matrix.effective_rank == 1
    @test maximum(abs, exact_control_matrix.M .- tsvd_control_matrix.M) > 0
    @test ControlMatrix(D_sing).policy isa TSVDInverse

    imat = InteractionMatrix(D_sing, 0.1)
    recon_exact = ModalReconstructor(imat; policy=ExactPseudoInverse())
    recon_tsvd = ModalReconstructor(imat; policy=TSVDInverse(rtol=1e-9))
    recon_factorized = FactorizedReconstructor(imat;
        policy=ExactPseudoInverse())
    recon_rank_one = FactorizedReconstructor(imat;
        policy=ExactPseudoInverse(), max_rank=1)
    @test recon_exact.effective_rank == 2
    @test recon_tsvd.effective_rank == 1
    @test ModalReconstructor(imat).policy isa TSVDInverse
    @test AdaptiveOpticsSim.factorized_rank(recon_factorized) == 2
    @test AdaptiveOpticsSim.factorized_rank(recon_rank_one) == 1
    probe_slopes = [0.3, -0.2]
    dense_command = zeros(2)
    factorized_command = zeros(2)
    reconstruct!(dense_command, recon_exact, probe_slopes)
    reconstruct!(factorized_command, recon_factorized, probe_slopes)
    @test factorized_command ≈ dense_command atol=1e-12 rtol=1e-12
    reconstruct!(factorized_command, recon_factorized, probe_slopes)
    if coverage_instrumented()
        @test_skip "allocation assertions are disabled under coverage instrumentation"
    else
        @test @allocated(reconstruct!(factorized_command,
            recon_factorized, probe_slopes)) == 0
    end
    @test sum(length,
        AdaptiveOpticsSim.runtime_reconstructor_storage(recon_rank_one)) <
        sum(length,
            AdaptiveOpticsSim.runtime_reconstructor_storage(recon_factorized))
end
