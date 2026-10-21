#=
Compare CUDA workgroups for the existing 352-square Shack--Hartmann plant.

    JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 julia --startup-file=no \
        --project=test/cuda benchmarks/benchmark_sh_pupil_workgroups_cuda.jl cuda

Set AOS_SH_RESULTS to a TOML output path. AOS_SH_SAMPLES (default 256) and
AOS_SH_ROUNDS (default 4) control each half of an ABBA/BAAB comparison.
Set AOS_SH_MATRIX=1 to compare the default launch, an unchanged control, and
(44, 8, 1, 1) across representative pupil and lenslet dimensions. In matrix
mode AOS_SH_RESULTS is a filename prefix; each case adds -N<resolution>-L<lenslets>.toml.
Preparation, compilation, input changes, oracles, event allocation and result
serialization stay outside timing. This is concurrency-one closed-loop service
measurement, not a fixed-rate deadline or rare-tail qualification.

Only this process specializes the pupil kernel's launcher. Production dispatch
and every other kernel retain their existing workgroup choices. Captured plant
replay contains no workgroup selection or benchmark instrumentation.
=#
import CUDA
using TOML
using LinearAlgebra
using Statistics
using FixedSizeArrays: FixedSizeVectorDefault
using Test: @test, @inferred

include("profile_aos_gpu_plant.jl")

const SHWFS = AdaptiveOpticsSim.WavefrontSensors
const KA = AdaptiveOpticsSim.KernelAbstractions
const _SH_WORKGROUP = Ref{Union{Nothing,NTuple{4,Int}}}(nothing)
const _SH_SAMPLES = parse(Int, get(ENV, "AOS_SH_SAMPLES", "256"))
const _SH_ROUNDS = parse(Int, get(ENV, "AOS_SH_ROUNDS", "4"))
const _SH_KERNEL_REPLAYS = 64
const _SH_DEFAULT_GROUPS = (nothing, nothing, (32, 2, 1, 1), (32, 4, 1, 1),
    (32, 8, 1, 1), (16, 8, 1, 1), (16, 16, 1, 1), (32, 16, 1, 1),
    (44, 8, 1, 1), (44, 16, 1, 1), (8, 8, 2, 2))
const _SH_MATRIX_GROUPS = (nothing, nothing, (44, 8, 1, 1))
const _SH_MATRIX_CASES = ((128, 8), (256, 8), (256, 16), (352, 16), (384, 16),
    (512, 16), (512, 32), (704, 32), (768, 32))

struct SHWorkgroupBlock
    candidate::Int
    round::Int
    position::Int
    measured_index::Int
    mode::Symbol
    start_ns::UInt64
    end_ns::UInt64
    allocated_bytes::Int
    collections::Int
    samples_ns::Vector{Float64}
end

# The exact singleton dispatch isolates the pupil launch; other four-dimensional
# launches, including spot sampling, continue to use AOS's ordinary helper.
function Backends.launch_kernel_async!(style::Backends.AcceleratorStyle,
    kernel::typeof(SHWFS.sh_explicit_pupil_stack_kernel!), args...; ndrange)
    group = _SH_WORKGROUP[]
    if group === nothing
        kernel(style.backend)(args...; ndrange)
    else
        kernel(style.backend)(args...; ndrange, workgroupsize=group)
    end
    return nothing
end

function sh_pupil_arguments(prepared)
    optics, input = prepared.optics, prepared.input
    workspace = SHWFS.microlens_propagation_workspace(optics.propagation)
    n = input.metadata.dimensions[1]
    n_sub = SHWFS.n_lenslets(optics)
    sub = div(n, n_sub)
    pad = size(workspace.fft_stack, 1)
    offset = div(pad - sub, 2)
    T = eltype(workspace.intensity_stack)
    area = T(prod(input.metadata.sampling))
    amplitude_scale = sqrt(T(SHWFS.photon_irradiance(prepared.plan.source)) * area)
    opd_to_cycles = T(2) / T(prepared.plan.wavelength_m)
    args = (workspace.fft_stack, optics.front_end.layout.valid_mask,
        input.amplitude, input.opd, workspace.phasor, amplitude_scale,
        opd_to_cycles, n_sub, sub, offset, offset, n, pad)
    return args, (pad, pad, n_sub, n_sub)
end

function sh_launch_pupil!(context, args, ndrange)
    Backends._with_prepared_device_execution_context(context) do
        Backends.launch_kernel_async!(Backends.execution_style(first(args)),
            SHWFS.sh_explicit_pupil_stack_kernel!, args...; ndrange)
    end
    return nothing
end

function sh_capture_kernel(context, args, ndrange)
    sh_launch_pupil!(context, args, ndrange)
    Backends._synchronize_prepared_device_execution_context!(context)
    return Backends._with_prepared_device_execution_context(context) do
        Backends._capture_prepared_device_graph(context) do
            for _ in 1:_SH_KERNEL_REPLAYS
                Backends.launch_kernel_async!(Backends.execution_style(first(args)),
                    SHWFS.sh_explicit_pupil_stack_kernel!, args...; ndrange)
            end
        end
    end
end

function sh_native_geometry(graph)
    driver = CUDA.CUDACore
    count = Ref{Csize_t}(0)
    driver.cuGraphGetNodes(graph, C_NULL, count)
    nodes = Vector{driver.CUgraphNode}(undef, count[])
    driver.cuGraphGetNodes(graph, nodes, count)
    @assert length(nodes) == _SH_KERNEL_REPLAYS
    params = Ref{driver.CUDA_KERNEL_NODE_PARAMS}()
    driver.cuGraphKernelNodeGetParams_v2(first(nodes), params)
    p = params[]
    registers = Ref{Cint}(0)
    driver.cuFuncGetAttribute(registers, driver.CU_FUNC_ATTRIBUTE_NUM_REGS, p.func)
    return (registers=Int(registers[]), block=[Int(p.blockDimX), Int(p.blockDimY), Int(p.blockDimZ)],
        grid=[Int(p.gridDimX), Int(p.gridDimY), Int(p.gridDimZ)])
end

function sh_kernel_times!(times, captured, context, start_event, stop_event)
    Backends._with_prepared_device_execution_context(context) do
        for index in eachindex(times)
            CUDA.record(start_event, context.stream)
            Backends._launch_prepared_device_graph!(captured, context)
            CUDA.record(stop_event, context.stream)
            CUDA.synchronize(stop_event)
            times[index] = CUDA.elapsed(start_event, stop_event) * 1e9 /
                _SH_KERNEL_REPLAYS
        end
    end
    return nothing
end

function sh_plant_times!(times, graph)
    for index in eachindex(times)
        started = time_ns()
        _complete_step!(graph)
        times[index] = Float64(time_ns() - started)
    end
    return nothing
end

function sh_summary(times)
    sorted = sort(times)
    return (count=length(times), mean_ns=mean(times),
        p50_ns=median(times), p90_ns=sorted[cld(9length(times), 10)],
        p99_ns=sorted[cld(99length(times), 100)], max_ns=maximum(times))
end

function sh_fixture(definition, target, group, reference)
    _SH_WORKGROUP[] = group
    graph = prepare_algorithm_graph(definition; target,
        execution=CapturedGraphExecution())
    _replay!(graph, 256)
    _assert_reference_frame!(reference, graph, "workgroup $group")
    _warm_and_measure_allocation!(graph, "workgroup $group"; require_zero=true)
    prepared = AlgorithmGraphs.prepared_graph_node(graph, Val(:shwfs)).prepared
    args, ndrange = sh_pupil_arguments(prepared)
    context = Backends._prepare_device_execution_context(first(args))
    kernel = sh_capture_kernel(context, args, ndrange)
    return (; graph, kernel, context, args, ndrange,
        start_event=CUDA.CuEvent(), stop_event=CUDA.CuEvent())
end

function sh_pair!(records, fixtures, candidate, round, mode)
    # Alternate ABBA/BAAB and rotate candidate order between rounds.
    order = isodd(round) ? (1, candidate, candidate, 1) :
        (candidate, 1, 1, candidate)
    for (position, index) in enumerate(order)
        fixture = fixtures[index]
        times = Vector{Float64}(undef, mode === :plant ? _SH_SAMPLES : 32)
        _replay!(fixture.graph, 32)
        memory_before = CUDA.used_memory()
        gc_before = Base.gc_num()
        block_started = time_ns()
        if mode === :plant
            sh_plant_times!(times, fixture.graph)
        else
            sh_kernel_times!(times, fixture.kernel, fixture.context,
                fixture.start_event, fixture.stop_event)
        end
        block_ended = time_ns()
        gc_delta = Base.GC_Diff(Base.gc_num(), gc_before)
        @assert CUDA.used_memory() == memory_before
        push!(records, SHWorkgroupBlock(candidate, round, position, index, mode,
            block_started, block_ended, gc_delta.allocd, gc_delta.pause, times))
    end
end

function benchmark_sh_workgroups(;
    pupil_resolution::Int=_PUPIL_RESOLUTION,
    lenslets_per_axis::Int=_LENSLETS_PER_AXIS,
    pixels_per_lenslet::Int=div(pupil_resolution, lenslets_per_axis),
    groups=_SH_DEFAULT_GROUPS,
    output::AbstractString=get(ENV, "AOS_SH_RESULTS", "/tmp/aos-sh-workgroups.toml"))
    _SH_SAMPLES >= 100 || error("AOS_SH_SAMPLES must be at least 100")
    _SH_ROUNDS > 0 || error("AOS_SH_ROUNDS must be positive")
    CUDA.functional() || error("CUDA is not functional")
    CUDA.allowscalar(false)
    BLAS.set_num_threads(1)
    Backends.set_fft_provider_threads!(1)
    command, coordinates = _pdm_command(), _actuator_coordinates()
    uncompensated = zeros(Float32, pupil_resolution, pupil_resolution)
    cpu = prepare_algorithm_graph(_plant_definition(command, uncompensated,
        coordinates; name=:sh_workgroup_cpu, pupil_resolution, lenslets_per_axis,
        pixels_per_lenslet))
    _complete_step!(cpu)
    reference = copy(graph_output(cpu, Val(:frame)))
    gpu_command, gpu_coordinates = CUDA.CuArray(command), CUDA.CuArray(coordinates)
    gpu_opd = CUDA.CuArray(uncompensated)
    target = Backends.compute_device(gpu_command)
    definition = _plant_definition(gpu_command, gpu_opd, gpu_coordinates;
        name=:sh_workgroup_cuda, pupil_resolution, lenslets_per_axis,
        pixels_per_lenslet)
    first_fixture = sh_fixture(definition, target, first(groups), reference)
    builder = Vector{typeof(first_fixture)}(undef, length(groups))
    builder[1] = first_fixture
    for index in 2:length(groups)
        println("prepare_workgroup=$(groups[index])")
        flush(stdout)
        builder[index] = sh_fixture(definition, target, groups[index], reference)
    end
    fixtures = FixedSizeVectorDefault(builder)
    # Compile the exact measurement loops before any recorded block.
    for fixture in fixtures
        sh_plant_times!(zeros(2), fixture.graph)
        sh_kernel_times!(zeros(2), fixture.kernel, fixture.context,
            fixture.start_event, fixture.stop_event)
    end
    # Compare the raw complex pupil stack against the default launch, including
    # padding and invalid lenslets. The kernel contains no reductions.
    default_stack = reinterpret(UInt8, vec(Array(first(first(fixtures).args))))
    for fixture in fixtures
        @assert reinterpret(UInt8, vec(Array(first(fixture.args)))) == default_stack
    end
    # Changing-command CPU oracles precede timing and check full products.
    for scale in (0.0f0, -0.7f0, 1.3f0)
        copyto!(command, scale .* _pdm_command())
        copyto!(gpu_command, command)
        # Captured launches cannot discover a default-stream input write.
        CUDA.synchronize()
        _complete_step!(cpu)
        for fixture in fixtures
            _complete_step!(fixture.graph)
            _assert_reference_frame!(graph_output(cpu, Val(:frame)), fixture.graph,
                "changing command")
            Backends._with_prepared_device_execution_context(fixture.context) do
                Backends._launch_prepared_device_graph!(fixture.kernel, fixture.context)
                Backends._synchronize_prepared_device_execution_context!(fixture.context)
            end
        end
        changed_stack = reinterpret(UInt8, vec(Array(first(first(fixtures).args))))
        for fixture in fixtures
            @assert reinterpret(UInt8, vec(Array(first(fixture.args)))) == changed_stack
        end
    end
    copyto!(command, _pdm_command())
    copyto!(gpu_command, command)
    CUDA.synchronize()
    _complete_step!(cpu)
    for fixture in fixtures
        _replay!(fixture.graph, 256)
    end
    GC.gc(true)
    records = SHWorkgroupBlock[]
    for round in 1:_SH_ROUNDS
        candidates = circshift(collect(2:length(groups)), round - 1)
        for candidate in candidates
            for mode in (:plant, :kernel)
                sh_pair!(records, fixtures, candidate, round, mode)
            end
            println("completed_round=$round candidate=$(groups[candidate])")
            flush(stdout)
        end
    end
    for fixture in fixtures
        _complete_step!(fixture.graph)
        _assert_reference_frame!(graph_output(cpu, Val(:frame)), fixture.graph,
            "post measurement")
        _warm_and_measure_allocation!(fixture.graph, "post measurement"; require_zero=true)
    end
    metadata = Dict{String,Union{Int,String}}(
        "source_revision" => String(readchomp(`git rev-parse HEAD`)),
        "julia" => string(VERSION), "cuda_jl" => string(Base.pkgversion(CUDA)),
        "kernel_abstractions" => string(Base.pkgversion(KA)),
        "device" => CUDA.name(CUDA.device()), "runtime" => string(CUDA.runtime_version()),
        "pupil_resolution" => pupil_resolution, "lenslets_per_axis" => lenslets_per_axis,
        "pixels_per_lenslet" => pixels_per_lenslet,
        "samples_per_plant_block" => _SH_SAMPLES, "rounds" => _SH_ROUNDS,
        "compute_threads" => Threads.nthreads(:default),
        "interactive_threads" => Threads.nthreads(:interactive),
        "gc_threads" => Threads.ngcthreads(), "blas_threads" => BLAS.get_num_threads(),
        "correctness_changed_commands" => "0.0, -0.7, 1.3 times nominal deterministic command",
        "measurement_input" => "Unchanged nominal command within each block; observational repeated steady-state plant workload",
        "measurement_order" => "ABBA on odd rounds, BAAB on even rounds; candidate order rotates between rounds; plant then kernel per candidate",
        "affinity" => get(ENV, "AOS_SH_AFFINITY", "unbound"),
        "kernel_replays_per_event_sample" => _SH_KERNEL_REPLAYS,
        "kernel_boundary" => "CUDA events around one graph of 64 repeated pupil kernels; graph scheduling overhead amortized per kernel",
        "plant_boundary" => "host time_ns around complete captured plant step including wait; concurrency one",
        "oracles" => "CPU complete frames before/after and three changed commands; raw pupil stack bitwise equal to default")
    open(output, "w") do io
        TOML.print(io, metadata)
        for (index, group) in enumerate(groups)
            geometry = sh_native_geometry(fixtures[index].kernel.graph)
            println(io, "\n[[geometry]]")
            TOML.print(io, Dict{String,Union{Int,String,Vector{Int}}}(
                "index" => index, "workgroup" => group === nothing ? "default" : string(group),
                "role" => index == 1 ? "baseline" : index == 2 ? "unchanged_control" : "candidate",
                "ndrange" => collect(fixtures[index].ndrange),
                "registers_per_thread" => geometry.registers,
                "block" => geometry.block, "grid" => geometry.grid))
        end
        for record in records
            summary = sh_summary(record.samples_ns)
            println(io, "\n[[blocks]]")
            fields = Dict{String,Union{Int,Float64,String,Vector{Float64}}}(
                "candidate" => record.candidate, "round" => record.round,
                "position" => record.position, "measured_index" => record.measured_index,
                "mode" => string(record.mode), "samples_ns" => record.samples_ns,
                "start_monotonic_ns" => Int(record.start_ns),
                "end_monotonic_ns" => Int(record.end_ns),
                "gc_allocated_bytes" => record.allocated_bytes,
                "gc_collections" => record.collections)
            for (key, value) in pairs(summary)
                fields[string(key)] = value
            end
            TOML.print(io, fields)
        end
    end
    println("results=$output")
    return nothing
end

function sh_matrix_output(prefix::AbstractString, pupil_resolution::Int,
    lenslets_per_axis::Int)
    stem, extension = splitext(prefix)
    base = extension == ".toml" ? stem : prefix
    return "$base-N$pupil_resolution-L$lenslets_per_axis.toml"
end

function benchmark_sh_workgroup_matrix(;
    output_prefix::AbstractString=get(ENV, "AOS_SH_RESULTS", "/tmp/aos-sh-workgroups.toml"))
    CUDA.synchronize()
    GC.gc(true)
    for (pupil_resolution, lenslets_per_axis) in _SH_MATRIX_CASES
        println("case_pupil_resolution=$pupil_resolution lenslets_per_axis=$lenslets_per_axis")
        flush(stdout)
        benchmark_sh_workgroups(; pupil_resolution, lenslets_per_axis,
            groups=_SH_MATRIX_GROUPS,
            output=sh_matrix_output(output_prefix, pupil_resolution, lenslets_per_axis))
        # The completed case has returned without retaining its graph owners.
        # Collection and pool reclamation precede preparation of the next case.
        CUDA.synchronize()
        GC.gc(true)
        CUDA.reclaim()
    end
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    if get(ENV, "AOS_SH_MATRIX", "0") == "1"
        benchmark_sh_workgroup_matrix()
    else
        benchmark_sh_workgroups()
    end
end
