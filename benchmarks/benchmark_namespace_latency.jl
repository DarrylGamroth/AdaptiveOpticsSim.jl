using TOML

const REPOSITORY_ROOT = normpath(joinpath(@__DIR__, ".."))
const CONTRACT_PATH =
    joinpath(@__DIR__, "contracts", "namespace_latency.toml")

function representative_namespace_workflow()
    tel = AdaptiveOpticsSim.Telescope(
        resolution=16,
        diameter=8.0,
        central_obstruction=0.14,
        T=Float64,
        backend=AdaptiveOpticsSim.Backends.CPUBackend(),
    )
    src = AdaptiveOpticsSim.Source(
        band=:custom,
        wavelength=8.0e-7,
        photon_irradiance=2.0,
    )
    pupil = AdaptiveOpticsSim.PupilFunction(tel)
    prepared = AdaptiveOpticsSim.prepare_direct_imaging(
        pupil,
        src;
        zero_padding=2,
    )
    output = AdaptiveOpticsSim.form_direct_image!(prepared)
    return sum(AdaptiveOpticsSim.intensity_values(output))
end

function run_probe()
    GC.gc()
    load_start = time_ns()
    @eval using AdaptiveOpticsSim
    load_ns = time_ns() - load_start

    @eval using LinearAlgebra
    Core.eval(Main, :(LinearAlgebra.BLAS.set_num_threads(1)))
    Core.eval(Main, :(AdaptiveOpticsSim.Backends.set_fft_provider_threads!(1)))

    GC.gc()
    first_start = time_ns()
    first_result = Base.invokelatest(representative_namespace_workflow)
    first_call_ns = time_ns() - first_start

    GC.gc()
    warmed_start = time_ns()
    warmed_result = Base.invokelatest(representative_namespace_workflow)
    warmed_call_ns = time_ns() - warmed_start

    correctness = isfinite(first_result) && first_result > 0 &&
        first_result == warmed_result
    TOML.print(stdout, Dict(
        "load_ns" => load_ns,
        "first_call_ns" => first_call_ns,
        "warmed_call_ns" => warmed_call_ns,
        "first_result" => first_result,
        "warmed_result" => warmed_result,
        "correctness_passed" => correctness,
    ))
    return correctness ? 0 : 1
end

function probe_command()
    julia = Base.julia_cmd()
    return `$julia --startup-file=no --history-file=no --threads=1 --project=$REPOSITORY_ROOT $(@__FILE__) --probe`
end

function run_fresh_probe()
    output = read(setenv(probe_command(), "JULIA_NUM_THREADS" => "1"), String)
    return TOML.parse(output)
end

function integer_summary(values::Vector{Int})
    sorted = sort(values)
    return Dict(
        "minimum_ns" => first(sorted),
        "median_ns" => sorted[cld(length(sorted), 2)],
        "maximum_ns" => last(sorted),
    )
end

function git_output(arguments...)
    return strip(read(`git -C $REPOSITORY_ROOT $(arguments)`, String))
end

function run_characterization(output_path::AbstractString)
    contract = TOML.parsefile(CONTRACT_PATH)
    runs = Int(contract["fresh_process_runs"])
    observations = [run_fresh_probe() for _ in 1:runs]
    all(observation["correctness_passed"] for observation in observations) ||
        error("a namespace-latency correctness probe failed")

    load_samples = Int[observation["load_ns"] for observation in observations]
    first_samples =
        Int[observation["first_call_ns"] for observation in observations]
    warmed_samples =
        Int[observation["warmed_call_ns"] for observation in observations]

    head = git_output("rev-parse", "HEAD")
    artifact = Dict(
        "schema_version" => 1,
        "name" => contract["name"],
        "source_contract" => relpath(CONTRACT_PATH, REPOSITORY_ROOT),
        "characterized_source_revision" => contract["source_revision"],
        "repository_head" => head,
        "repository_dirty" => !isempty(git_output(
            "status", "--porcelain", "--untracked-files=normal")),
        "package_source_dirty" => !isempty(git_output(
            "status",
            "--porcelain",
            "--untracked-files=normal",
            "--",
            "src",
            "Project.toml",
        )),
        "julia_version" => string(VERSION),
        "fresh_process_runs" => runs,
        "correctness_passed" => true,
        "load_samples_ns" => load_samples,
        "first_call_samples_ns" => first_samples,
        "warmed_call_samples_ns" => warmed_samples,
        "load_summary" => integer_summary(load_samples),
        "first_call_summary" => integer_summary(first_samples),
        "warmed_call_summary" => integer_summary(warmed_samples),
        "results" => observations,
    )

    mkpath(dirname(output_path))
    open(output_path, "w") do io
        TOML.print(io, artifact; sorted=true)
    end
    TOML.print(stdout, Dict(
        "output_path" => relpath(output_path, REPOSITORY_ROOT),
        "load_median_ns" => artifact["load_summary"]["median_ns"],
        "first_call_median_ns" =>
            artifact["first_call_summary"]["median_ns"],
        "warmed_call_median_ns" =>
            artifact["warmed_call_summary"]["median_ns"],
    ); sorted=true)
    return nothing
end

function parse_output_path(arguments)
    index = findfirst(==("--output"), arguments)
    index === nothing && return nothing
    index < length(arguments) ||
        throw(ArgumentError("--output requires a path"))
    return abspath(arguments[index + 1])
end

function main(arguments=ARGS)
    if arguments == ["--probe"]
        exit(run_probe())
    end
    contract = TOML.parsefile(CONTRACT_PATH)
    output = parse_output_path(arguments)
    if output === nothing
        output = joinpath(
            REPOSITORY_ROOT,
            String(contract["artifact_path"]),
        )
    end
    run_characterization(output)
    return nothing
end

main()
