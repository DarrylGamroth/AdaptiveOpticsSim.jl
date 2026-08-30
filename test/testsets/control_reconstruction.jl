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

    matrix_parent = reshape(Float32.(1:30), 5, 6)
    matrix_view = @view matrix_parent[2:4, 2:5]
    expected_matrix = copy(matrix_view)
    prepared_matrix = @inferred ControlMatrixPlan(matrix_view)
    fill!(matrix_view, 0.0f0)
    @test prepared_matrix.matrix == expected_matrix
    @test Control.runtime_reconstructor_storage(prepared_matrix) ===
        (prepared_matrix.matrix,)

    slopes_parent = Float32[-1, 0.5, -0.25, 0.75, 1.25, 2]
    slopes_view = @view slopes_parent[2:5]
    output_parent = fill(-1.0f0, 5)
    output_view = @view output_parent[2:4]
    expected_output = expected_matrix * collect(slopes_view)
    @test (@inferred reconstruct!(
        output_view,
        prepared_matrix,
        slopes_view,
    )) === output_view
    @test output_view ≈ expected_output
    @test reconstruct(prepared_matrix, slopes_view) ≈ expected_output
    reconstruct!(output_view, prepared_matrix, slopes_view)
    if coverage_instrumented()
        @test_skip "allocation assertions are disabled under coverage instrumentation"
    else
        @test @allocated(reconstruct!(
            output_view,
            prepared_matrix,
            slopes_view,
        )) == 0
    end

    @test_throws InvalidConfiguration ControlMatrixPlan(zeros(Float32, 0, 2))
    @test_throws InvalidConfiguration ControlMatrixPlan(zeros(Float32, 2, 0))
    @test_throws DimensionMismatchError reconstruct!(
        zeros(Float32, 2),
        prepared_matrix,
        slopes_view,
    )
    @test_throws DimensionMismatchError reconstruct!(
        output_view,
        prepared_matrix,
        zeros(Float32, 3),
    )
    square_plan = ControlMatrixPlan(Float32[1 0; 0 1])
    aliased = zeros(Float32, 2)
    @test_throws InvalidConfiguration reconstruct!(
        aliased,
        square_plan,
        aliased,
    )
    @test_throws InvalidConfiguration reconstruct!(
        @view(square_plan.matrix[:, 1]),
        square_plan,
        Float32[1, 2],
    )

    imat = InteractionMatrix(D_sing, 0.1)
    recon_exact = @inferred ModalReconstructor(
        imat;
        policy=ExactPseudoInverse(),
    )
    recon_tsvd = ModalReconstructor(imat; policy=TSVDInverse(rtol=1e-9))
    recon_factorized = @inferred FactorizedReconstructor(imat;
        policy=ExactPseudoInverse())
    recon_rank_one = FactorizedReconstructor(imat;
        policy=ExactPseudoInverse(), max_rank=1)
    @test recon_exact.effective_rank == 2
    @test recon_tsvd.effective_rank == 1
    @test ModalReconstructor(imat).policy isa TSVDInverse
    @test Control.factorized_rank(recon_factorized) == 2
    @test Control.factorized_rank(recon_rank_one) == 1
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
        Control.runtime_reconstructor_storage(recon_rank_one)) <
        sum(length,
            Control.runtime_reconstructor_storage(recon_factorized))
end
