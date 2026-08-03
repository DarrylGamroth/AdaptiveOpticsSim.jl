using TOML
using Pkg

const FIXED_SIZE_ARRAYS_PROBE_MODE = ARGS == ["--probe"]
FIXED_SIZE_ARRAYS_PROBE_MODE || (@eval using FixedSizeArrays)

const REPOSITORY_ROOT = normpath(joinpath(@__DIR__, ".."))
const CONTRACT_PATH =
    joinpath(@__DIR__, "contracts", "fixed_size_arrays.toml")
const ISOLATED_PROJECT = joinpath(@__DIR__, "fixed_size_arrays")

function integer_summary(values::Vector{Int})
    sorted = sort(values)
    return Dict(
        "minimum_ns" => first(sorted),
        "median_ns" => sorted[cld(length(sorted), 2)],
        "maximum_ns" => last(sorted),
    )
end

function run_probe()
    contract = TOML.parsefile(CONTRACT_PATH)
    count = Int(contract["vector_length"])

    GC.gc()
    load_start = time_ns()
    @eval using FixedSizeArrays
    load_ns = Int(time_ns() - load_start)

    GC.gc()
    first_start = time_ns()
    values = Base.invokelatest() do
        result = FixedSizeArrays.FixedSizeVectorDefault{Float64}(
            undef, count)
        fill!(result, 1.0)
        result
    end
    first_construction_ns = Int(time_ns() - first_start)

    correctness = Base.invokelatest(values) do loaded_values
        length(loaded_values) == count && sum(loaded_values) == count
    end
    TOML.print(stdout, Dict(
        "load_ns" => load_ns,
        "first_construction_ns" => first_construction_ns,
        "correctness_passed" => correctness,
    ); sorted=true)
    return correctness ? 0 : 1
end

@noinline function indexed_sum(values)
    result = zero(eltype(values))
    @inbounds for index in eachindex(values)
        result += values[index]
    end
    return result
end

@noinline function iterated_sum(values)
    result = zero(eltype(values))
    for value in values
        result += value
    end
    return result
end

@noinline copy_values!(destination, source) = copyto!(destination, source)

@inline project_armed_storage_eltype_allowed(::Type{T}) where {T} =
    isconcretetype(T)

function warmed_allocation(operation)
    operation()
    operation()
    GC.gc()
    return @allocated operation()
end

function construction_allocation(constructor)
    constructor()
    GC.gc()
    return @allocated constructor()
end

function probe_command()
    julia = Base.julia_cmd()
    return `$julia --startup-file=no --history-file=no --threads=1 --project=$ISOLATED_PROJECT $(@__FILE__) --probe`
end

function run_fresh_probe()
    command = setenv(probe_command(), "JULIA_NUM_THREADS" => "1")
    return TOML.parse(read(command, String))
end

function resolved_dependencies()
    dependencies = Pkg.dependencies()
    DependencyValue = Union{Bool,String}
    selected = Dict{String,DependencyValue}[]
    for (uuid, dependency) in pairs(dependencies)
        dependency.name in ("FixedSizeArrays", "Collects") || continue
        push!(selected, Dict{String,DependencyValue}(
            "direct" => dependency.is_direct_dep,
            "name" => dependency.name,
            "uuid" => string(uuid),
            "version" => string(dependency.version),
        ))
    end
    sort!(selected; by=entry -> entry["name"])
    return selected
end

function git_output(arguments...)
    return strip(read(`git -C $REPOSITORY_ROOT $(arguments)`, String))
end

function run_characterization(output_path::AbstractString)
    contract = TOML.parsefile(CONTRACT_PATH)
    count = Int(contract["vector_length"])
    shape = Tuple(Int.(contract["matrix_shape"]))
    cpu_information = Sys.cpu_info()
    cpu_model = isempty(cpu_information) ?
        "unknown" : first(cpu_information).model
    observations = [run_fresh_probe() for _ in
        1:Int(contract["fresh_process_runs"])]
    all(observation["correctness_passed"] for observation in observations) ||
        error("a FixedSizeArrays fresh-process probe failed")

    fixed = FixedSizeVectorDefault{Float64}(undef, count)
    source = FixedSizeVectorDefault{Float64}(undef, count)
    fill!(fixed, 1.0)
    fill!(source, 2.0)
    indexed_sum(fixed)
    iterated_sum(fixed)
    copy_values!(fixed, source)

    indexed_allocated = warmed_allocation(() -> indexed_sum(fixed))
    iterated_allocated = warmed_allocation(() -> iterated_sum(fixed))
    copy_allocated = warmed_allocation(() -> copy_values!(fixed, source))

    short = FixedSizeVectorDefault{Float64}(undef, 8)
    long = FixedSizeVectorDefault{Float64}(undef, 16)
    matrix = FixedSizeMatrixDefault{Float64}(undef, shape)
    alternate_matrix = FixedSizeMatrixDefault{Float64}(
        undef, (first(shape) + 1, last(shape)))
    any_values = FixedSizeVectorDefault{Any}(undef, 2)
    any_values[1] = 1
    any_values[2] = 2.0

    resize_rejected = try
        resize!(short, 4)
        false
    catch error
        error isa MethodError
    end
    push_rejected = try
        push!(short, 1.0)
        false
    catch error
        error isa MethodError
    end

    memory_constructor = () -> begin
        values = Memory{Float64}(undef, count)
        fill!(values, 1.0)
        values
    end
    vector_constructor = () -> fill(1.0, count)
    fixed_constructor = () -> begin
        values = FixedSizeVectorDefault{Float64}(undef, count)
        fill!(values, 1.0)
        values
    end

    memory_values = memory_constructor()
    vector_values = vector_constructor()
    artifact = Dict(
        "schema_version" => 1,
        "name" => contract["name"],
        "source_contract" => relpath(CONTRACT_PATH, REPOSITORY_ROOT),
        "repository_head" => git_output("rev-parse", "HEAD"),
        "repository_dirty" => !isempty(git_output(
            "status", "--porcelain", "--untracked-files=normal")),
        "kernel" => string(Sys.KERNEL),
        "architecture" => string(Sys.ARCH),
        "cpu_target" => string(Sys.CPU_NAME),
        "cpu_model" => cpu_model,
        "logical_cpu_threads" => Sys.CPU_THREADS,
        "julia_version" => string(VERSION),
        "julia_threads" => Threads.nthreads(),
        "julia_cpu_target_env" =>
            get(ENV, "JULIA_CPU_TARGET", "default"),
        "fixed_size_arrays_version" => string(Base.pkgversion(
            FixedSizeArrays)),
        "resolved_dependencies" => resolved_dependencies(),
        "fresh_process_runs" => length(observations),
        "correctness_passed" => true,
        "load_samples_ns" =>
            Int[entry["load_ns"] for entry in observations],
        "first_construction_samples_ns" =>
            Int[entry["first_construction_ns"] for entry in observations],
        "load_summary" => integer_summary(
            Int[entry["load_ns"] for entry in observations]),
        "first_construction_summary" => integer_summary(
            Int[entry["first_construction_ns"] for entry in observations]),
        "warmed_allocated_bytes" => Dict(
            "indexing" => indexed_allocated,
            "iteration" => iterated_allocated,
            "copyto" => copy_allocated,
        ),
        "inference" => Dict(
            "indexed_sum" => Core.Compiler.return_type(
                indexed_sum, Tuple{typeof(fixed)}) === Float64,
            "iterated_sum" => Core.Compiler.return_type(
                iterated_sum, Tuple{typeof(fixed)}) === Float64,
            "copyto" => Core.Compiler.return_type(
                copy_values!, Tuple{typeof(fixed),typeof(source)}) ===
                    typeof(fixed),
        ),
        "shape_contract" => Dict(
            "runtime_vector_lengths_share_type" =>
                typeof(short) === typeof(long),
            "runtime_matrix_shapes_share_type" =>
                typeof(matrix) === typeof(alternate_matrix),
            "matrix_dimensions" => collect(size(matrix)),
            "alternate_matrix_dimensions" =>
                collect(size(alternate_matrix)),
            "matrix_rank_in_type" => ndims(typeof(matrix)),
            "resize_rejected" => resize_rejected,
            "push_rejected" => push_rejected,
        ),
        "element_type_contract" => Dict(
            "any_is_constructible" => eltype(any_values) === Any,
            "concrete_is_accepted" =>
                project_armed_storage_eltype_allowed(eltype(fixed)),
            "any_is_rejected" =>
                !project_armed_storage_eltype_allowed(eltype(any_values)),
        ),
        "storage_comparison" => Dict(
            "fixed_size_array" => Dict(
                "construction_allocated_bytes" =>
                    construction_allocation(fixed_constructor),
                "summarysize_bytes" => Base.summarysize(fixed),
            ),
            "memory" => Dict(
                "construction_allocated_bytes" =>
                    construction_allocation(memory_constructor),
                "summarysize_bytes" => Base.summarysize(memory_values),
            ),
            "vector" => Dict(
                "construction_allocated_bytes" =>
                    construction_allocation(vector_constructor),
                "summarysize_bytes" => Base.summarysize(vector_values),
            ),
            "dependency_backing_inspected" => false,
        ),
        "results" => observations,
    )

    artifact["correctness_passed"] =
        all(iszero, values(artifact["warmed_allocated_bytes"])) &&
        all(identity, values(artifact["inference"])) &&
        artifact["shape_contract"]["runtime_vector_lengths_share_type"] &&
        artifact["shape_contract"]["runtime_matrix_shapes_share_type"] &&
        artifact["shape_contract"]["resize_rejected"] &&
        artifact["shape_contract"]["push_rejected"] &&
        artifact["element_type_contract"]["any_is_constructible"] &&
        artifact["element_type_contract"]["concrete_is_accepted"] &&
        artifact["element_type_contract"]["any_is_rejected"] &&
        !artifact["storage_comparison"]["dependency_backing_inspected"]
    artifact["correctness_passed"] ||
        error("FixedSizeArrays characterization contract failed")

    mkpath(dirname(output_path))
    open(output_path, "w") do io
        TOML.print(io, artifact; sorted=true)
    end
    TOML.print(stdout, Dict(
        "output_path" => relpath(output_path, REPOSITORY_ROOT),
        "load_median_ns" => artifact["load_summary"]["median_ns"],
        "first_construction_median_ns" =>
            artifact["first_construction_summary"]["median_ns"],
        "correctness_passed" => artifact["correctness_passed"],
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
    arguments == ["--probe"] && exit(run_probe())
    contract = TOML.parsefile(CONTRACT_PATH)
    output = parse_output_path(arguments)
    output === nothing && (output = joinpath(
        REPOSITORY_ROOT, String(contract["artifact_path"])))
    run_characterization(output)
    return nothing
end

main()
