include(joinpath(@__DIR__, "common.jl"))
function main(; n_iter::Int=4)
    result = run_closed_loop_example(
        (tel, pupil_samples) -> ZernikeWFS(
            tel;
            pupil_samples=pupil_samples,
            diffraction_padding=2,
        );
        n_iter=n_iter,
        seed=5,
        amplitude=1e-8,
    )
    @info "Closed-loop ZernikeWFS tutorial complete" final_residual=result.residual_after[end]
    return result
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
