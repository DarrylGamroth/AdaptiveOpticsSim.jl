include(joinpath(@__DIR__, "common.jl"))

function main(; n_iter::Int=4)
    result = run_closed_loop_example(
        (tel, pupil_samples) -> PyramidWFS(tel; pupil_samples=pupil_samples, mode=Diffractive(), modulation=1.0,
            modulation_points=4, diffraction_padding=2);
        n_iter=n_iter,
    )
    @info "Closed-loop Pyramid tutorial complete" final_residual=result.residual_after[end]
    return result
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
