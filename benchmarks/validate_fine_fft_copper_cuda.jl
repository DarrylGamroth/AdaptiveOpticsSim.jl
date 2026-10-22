# Fine-sampling COPPER integration check for FFT layout experiments.
# COPPER_ROOT OUTPUT_DIR [FROZEN_BASELINE_DIR]; source packages stay unchanged.
include("validate_modulation_batches_copper_cuda.jl")

function fine_fft_copper_main(args)
    length(args) in (2, 3) || error("COPPER_ROOT OUTPUT_DIR [FROZEN_BASELINE_DIR]")
    copper = realpath(args[1])
    directory = abspath(args[2])
    include(joinpath(copper, "benchmarks", "copper_convergence.jl"))
    include(joinpath(copper, "benchmarks", "copper_atmospheric.jl"))
    Base.invokelatest(main, [copper, directory, "1440:256"])
    baseline = length(args) == 3 ? abspath(args[3]) : directory
    Base.invokelatest() do
        batch_closed_loop!(directory, baseline, Main.CopperConvergence,
            Main.CopperAtmospheric; n=1440, points=256, frames=320)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    fine_fft_copper_main(ARGS)
end
