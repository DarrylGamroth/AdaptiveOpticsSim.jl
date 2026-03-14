include(joinpath(@__DIR__, "common.jl"))

function main(; n_iter::Int=4)
    result = run_closed_loop_example(
        (tel, n_subap) -> ShackHartmann(tel; n_subap=n_subap, mode=Diffractive(), pixel_scale=0.1, n_pix_subap=6);
        n_iter=n_iter,
    )
    @info "Closed-loop Shack-Hartmann tutorial complete" final_residual=result.residual_after[end]
    return result
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
