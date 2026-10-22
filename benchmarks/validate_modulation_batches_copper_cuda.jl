# Instrument integration check for an AOS-only implementation change.
# COPPER_ROOT OUTPUT_DIR [FROZEN_BASELINE_DIR]; source packages stay unchanged.
include("validate_modulation_storage_cuda.jl")

function batch_closed_loop!(directory,baseline,convergence,atmospheric;n=240,points=96,frames=320)
    record=TOML.parsefile(joinpath(baseline,"validation.toml"))["cases"]["n$(n)-m$(points)"]
    matrix=reshape(reinterpret(Float64,read(joinpath(baseline,record["interaction"]["file"]))),1432,277)
    reference=reshape(reinterpret(Float32,read(joinpath(baseline,record["rates"]["file"]))),64,64,13)[:,:,1]
    divisor=record["divisor"]
    decomposition=svd(matrix;full=false)
    keep=decomposition.S .> .05maximum(decomposition.S)
    reconstructor=decomposition.V[:,keep]*Diagonal(1 ./ decomposition.S[keep])*decomposition.U[:,keep]'
    target=AOS.Backends.compute_device(CUDA.zeros(Float32,1))
    system=convergence.prepare_case(n,points;target)
    screens=CUDA.CuArray.(atmospheric.initialized_screens(823))
    opd=CUDA.zeros(Float32,n,n)
    science_opd=CUDA.zeros(Float32,480,480)
    fixture=atmospheric.science_fixture(480)
    flat_signal=convergence.raw_slopes(reference,divisor)
    command=zeros(Float32,277)
    commands=zeros(Float32,277,frames);slopes=zeros(Float64,1432,frames)
    diagnostics=zeros(Float64,2,frames)
    for frame in 1:frames
        atmospheric.render!(opd,screens,frame,.15)
        rate=atmospheric.expose(system,opd,command)
        signal=convergence.raw_slopes(rate,divisor).-flat_signal
        commands[:,frame].=command;slopes[:,frame].=signal
        atmospheric.render!(science_opd,screens,frame,.15)
        diagnostics[:,frame].=atmospheric.science(fixture,Array(science_opd),command)
        command .-= Float32.(.15 .* (reconstructor*signal))
        @assert all(isfinite,command) && maximum(abs,command)<1e-3
    end
    tail=(frames-63):frames
    result=Dict("frames"=>frames,"seed"=>823,"r0_m"=>.15,"gain"=>.15,
        "frozen_calibration"=>abspath(baseline),"rank"=>count(keep),
        "commands"=>write_array(directory,"closed-commands",commands),
        "slopes"=>write_array(directory,"closed-slopes",slopes),
        "diagnostics"=>write_array(directory,"closed-diagnostics",diagnostics),
        "diagnostic_rows"=>["piston_removed_residual_opd_rms_m","strehl"],
        "tail_residual_rms_m"=>sqrt(mean(abs2,diagnostics[1,tail])),
        "tail_mean_strehl"=>mean(diagnostics[2,tail]),
        "maximum_command_m"=>maximum(abs,commands),"all_finite"=>all(isfinite,diagnostics),
        "scope"=>"Noiseless 277-actuator frozen baseline calibration; native five-layer atmosphere rendered independently at fixed N/M; host RTC and diagnostics excluded from plant timings")
    open(io->TOML.print(io,result;sorted=true),joinpath(directory,"closed-loop.toml"),"w")
end

function batch_copper_main(args)
    length(args) in (2,3) || error("COPPER_ROOT OUTPUT_DIR [FROZEN_BASELINE_DIR]")
    copper=realpath(args[1]);directory=abspath(args[2])
    include(joinpath(copper,"benchmarks","copper_convergence.jl"))
    include(joinpath(copper,"benchmarks","copper_atmospheric.jl"))
    Base.invokelatest(main,[copper,directory,"240:96"])
    baseline=length(args)==3 ? abspath(args[3]) : directory
    # Includes introduce methods in a newer world; invoke the assembled cold driver.
    Base.invokelatest(batch_closed_loop!,directory,baseline,Main.CopperConvergence,Main.CopperAtmospheric)
end
if abspath(PROGRAM_FILE) == @__FILE__
    batch_copper_main(ARGS)
end
