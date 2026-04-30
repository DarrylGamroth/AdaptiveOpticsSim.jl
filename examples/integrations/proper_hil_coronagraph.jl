Base.find_package("Proper") === nothing &&
    error("Proper.jl is not available in the active environment. Install it with `using Pkg; Pkg.add(\"Proper\")` before running this example.")

using AdaptiveOpticsSim
using Proper

include(joinpath(@__DIR__, "..", "support", "proper_hil_coronagraph_common.jl"))
using .ProperHILCoronagraphCommon

function main()
    backend_name = length(ARGS) >= 1 ? ARGS[1] : "cpu"
    n_steps = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 5
    backend = resolve_hil_backend(backend_name)
    ctx = build_proper_hil_context(; backend=backend)

    println("proper_hil_coronagraph_example")
    println("  backend: ", backend_name)
    println("  wavelength_um: ", ctx.wavelength_um)
    println("  command_length: ", length(command(ctx.scenario)))
    println("  slopes_length: ", length(slopes(ctx.scenario)))

    for _ in 1:n_steps
        psf, sampling = hil_step!(ctx)
        println("  step=", ctx.step_index,
            " slopes_l2=", sqrt(sum(abs2, slopes(ctx.scenario))),
            " psf_shape=", size(psf),
            " science_sampling=", sampling)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
