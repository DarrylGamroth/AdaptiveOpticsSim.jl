module PyRTCProcessHIL

using AdaptiveOpticsSim.AlgorithmGraphs
using LinearAlgebra

include(joinpath(@__DIR__, "pyrtc_shared_memory.jl"))
using .PyRTCSharedMemory

include(joinpath(@__DIR__, "..", "..", "support", "hil_reference_systems.jl"))
using .HILReferenceSystems

const PYRTC_STREAM_NAMES = ("wfs", "wfc", "signal", "signal2D")
const WORKER_PREFIX = "AOS_PYRTC_WORKER "

struct ProcessReferenceCase
    wavefront_sensor::Symbol
    frame_shape::Tuple{Int,Int}
    signal_shape::Tuple{Int}
    signal_2d_shape::Tuple{Int,Int}
    poke::Float32
    gain::Float32
    iterations::Int
    convergence_limit::Float32
end

@inline process_reference_case(::Val{:shack_hartmann}) = ProcessReferenceCase(
    :shack_hartmann,
    (64, 64),
    (104,),
    (16, 8),
    2.0f-8,
    0.4f0,
    15,
    1.0f-3,
)

@inline process_reference_case(::Val{:pyramid}) = ProcessReferenceCase(
    :pyramid,
    (36, 36),
    (344,),
    (16, 32),
    2.0f-8,
    0.4f0,
    15,
    1.0f-3,
)

function process_reference_case(wavefront_sensor::Symbol)
    wavefront_sensor in supported_wavefront_sensors() || throw(ArgumentError(
        "unsupported pyRTC process WFS '$wavefront_sensor'; expected one " *
        "of $(join(supported_wavefront_sensors(), ", "))",
    ))
    return process_reference_case(Val(wavefront_sensor))
end

@inline reference_label(::Val{:shack_hartmann}) = "Shack-Hartmann"
@inline reference_label(::Val{:pyramid}) = "Pyramid"

struct ProcessStreams{W,C,S,D}
    wfs::W
    wfc::C
    signal::S
    signal_2d::D
end

mutable struct PyRTCWorker
    process::Base.Process
    stopped::Bool
end

function require_available_stream_names()
    Sys.islinux() || error("the pyRTC shared-memory integration requires Linux")
    occupied = String[]
    for name in PYRTC_STREAM_NAMES, suffix in ("", "_meta", "_gpu_handle")
        path = joinpath("/dev/shm", name * suffix)
        ispath(path) && push!(occupied, path)
    end
    isempty(occupied) || error(
        "refusing to reuse active pyRTC shared-memory streams: " *
        join(occupied, ", "),
    )
    return nothing
end

function create_process_streams(case::ProcessReferenceCase)
    require_available_stream_names()
    wfs = nothing
    wfc = nothing
    signal = nothing
    signal_2d = nothing
    try
        wfs = create_stream("wfs", Float32, case.frame_shape)
        wfc = create_stream("wfc", Float32, (actuator_count(),))
        signal = create_stream("signal", Float32, case.signal_shape)
        signal_2d = create_stream(
            "signal2D",
            Float32,
            case.signal_2d_shape,
        )
        return ProcessStreams(wfs, wfc, signal, signal_2d)
    catch
        !isnothing(signal_2d) && close_and_unlink_noexcept!(signal_2d)
        !isnothing(signal) && close_and_unlink_noexcept!(signal)
        !isnothing(wfc) && close_and_unlink_noexcept!(wfc)
        !isnothing(wfs) && close_and_unlink_noexcept!(wfs)
        rethrow()
    end
end

function close_and_unlink_noexcept!(stream::PyRTCStream)
    try
        close(stream)
    catch
    end
    try
        unlink!(stream)
    catch
    end
    return nothing
end

function close_process_streams!(streams::ProcessStreams)
    close_and_unlink_noexcept!(streams.signal_2d)
    close_and_unlink_noexcept!(streams.signal)
    close_and_unlink_noexcept!(streams.wfc)
    close_and_unlink_noexcept!(streams.wfs)
    return nothing
end

function pyrtc_python()
    executable = get(
        ENV,
        "PYRTC_PYTHON",
        get(
            ENV,
            "JULIA_PYTHONCALL_EXE",
            something(Sys.which("python3"), ""),
        ),
    )
    isfile(executable) || error(
        "set PYRTC_PYTHON to the Python interpreter containing pyRTC dependencies",
    )
    return abspath(executable)
end

function worker_command(
    pyrtc_root::AbstractString,
    wavefront_sensor::Symbol,
    temporary_directory::AbstractString,
)
    return Cmd(String[
        pyrtc_python(),
        joinpath(@__DIR__, "pyrtc_process_worker.py"),
        "--sensor",
        String(wavefront_sensor),
        "--pyrtc-root",
        abspath(pyrtc_root),
        "--temporary-directory",
        abspath(temporary_directory),
    ])
end

function read_worker_message(process::Base.Process)
    while !eof(process)
        line = readline(process)
        startswith(line, WORKER_PREFIX) || continue
        return chop(line; head=length(WORKER_PREFIX), tail=0)
    end
    error("pyRTC worker exited without a control response")
end

function await_worker_message(
    worker::PyRTCWorker;
    timeout::Real=10.0,
)
    task = @async read_worker_message(worker.process)
    status = timedwait(() -> istaskdone(task), timeout; pollint=0.001)
    if status === :timed_out
        process_running(worker.process) && kill(worker.process)
        try
            wait(worker.process)
        catch
        end
        worker.stopped = true
        try
            wait(task)
        catch
        end
        error("pyRTC worker did not respond within $timeout seconds")
    end
    return fetch(task)
end

function start_worker(
    pyrtc_root::AbstractString,
    case::ProcessReferenceCase,
    temporary_directory::AbstractString,
)
    process = open(
        worker_command(
            pyrtc_root,
            case.wavefront_sensor,
            temporary_directory,
        ),
        "r+",
    )
    worker = PyRTCWorker(process, false)
    try
        ready = split(await_worker_message(worker; timeout=30.0))
        length(ready) == 3 && first(ready) == "READY" || error(
            "unexpected pyRTC worker startup response: $(join(ready, " "))",
        )
        parse(Int, ready[2]) == only(case.signal_shape) || error(
            "pyRTC worker signal length $(ready[2]) does not match " *
            "$(only(case.signal_shape))",
        )
        parse(Int, ready[3]) == actuator_count() || error(
            "pyRTC worker command length $(ready[3]) does not match " *
            "$(actuator_count())",
        )
        return worker
    catch
        stop_worker_noexcept!(worker)
        rethrow()
    end
end

function send_worker_command!(
    worker::PyRTCWorker,
    command::AbstractString,
    expected_response::AbstractString;
    timeout::Real=10.0,
)
    worker.stopped && error("pyRTC worker is stopped")
    write(worker.process, command, '\n')
    flush(worker.process)
    response = await_worker_message(worker; timeout)
    response == expected_response || error(
        "pyRTC worker responded '$response'; expected '$expected_response'",
    )
    return nothing
end

function stop_worker!(worker::PyRTCWorker)
    worker.stopped && return nothing
    if process_running(worker.process)
        try
            send_worker_command!(worker, "STOP", "STOPPED"; timeout=5.0)
        catch
            process_running(worker.process) && kill(worker.process)
        end
    end
    wait(worker.process)
    worker.stopped = true
    success(worker.process) || error("pyRTC worker exited unsuccessfully")
    return nothing
end

function stop_worker_noexcept!(worker::PyRTCWorker)
    worker.stopped && return nothing
    if process_running(worker.process)
        try
            write(worker.process, "STOP\n")
            flush(worker.process)
        catch
        end
        timedwait(
            () -> !process_running(worker.process),
            1.0;
            pollint=0.01,
        )
        process_running(worker.process) && kill(worker.process)
    end
    try
        wait(worker.process)
    catch
    end
    worker.stopped = true
    return nothing
end

function process_frame!(
    worker::PyRTCWorker,
    streams::ProcessStreams,
    frame::AbstractMatrix{Float32},
    signal::Vector{Float32};
    command::AbstractString="PROCESS",
    response::AbstractString="PROCESSED",
)
    publish!(streams.wfs, frame)
    send_worker_command!(worker, command, response)
    read_next!(signal, streams.signal; timeout=5.0)
    return signal
end

function set_flat_reference!(
    worker::PyRTCWorker,
    streams::ProcessStreams,
    prepared,
    signal::Vector{Float32},
)
    boundary = prepared.boundary
    sequence = step_hil_frame!(boundary)
    process_frame!(
        worker,
        streams,
        hil_frame_buffer(boundary),
        signal,
    )
    send_worker_command!(worker, "SET_REF", "REF_SET")
    process_frame!(
        worker,
        streams,
        hil_frame_buffer(boundary),
        signal,
    )
    return sequence, copy(signal)
end

function calibrate_interaction_matrix!(
    worker::PyRTCWorker,
    streams::ProcessStreams,
    prepared,
    signal::Vector{Float32};
    poke::Float32,
)
    boundary = prepared.boundary
    sequence, flat_signal = set_flat_reference!(
        worker,
        streams,
        prepared,
        signal,
    )
    interaction_matrix = Matrix{Float32}(
        undef,
        length(signal),
        actuator_count(),
    )
    positive_signal = similar(signal)

    for command_index in axes(interaction_matrix, 2)
        fill!(hil_command_buffer(boundary), 0.0f0)
        hil_command_buffer(boundary)[command_index] = poke
        adopt_hil_command!(boundary, sequence)
        sequence = step_hil_frame!(boundary)
        copyto!(
            positive_signal,
            process_frame!(
                worker,
                streams,
                hil_frame_buffer(boundary),
                signal,
            ),
        )

        fill!(hil_command_buffer(boundary), 0.0f0)
        hil_command_buffer(boundary)[command_index] = -poke
        adopt_hil_command!(boundary, sequence)
        sequence = step_hil_frame!(boundary)
        negative_signal = process_frame!(
            worker,
            streams,
            hil_frame_buffer(boundary),
            signal,
        )
        @views @. interaction_matrix[:, command_index] =
            (positive_signal - negative_signal) / (2 * poke)
    end

    fill!(hil_command_buffer(boundary), 0.0f0)
    adopt_hil_command!(boundary, sequence)
    reset_hil_boundary!(boundary)
    return flat_signal, interaction_matrix
end

function configure_worker_loop!(
    worker::PyRTCWorker,
    streams::ProcessStreams,
    interaction_matrix::Matrix{Float32},
    gain::Float32,
    temporary_directory::AbstractString,
)
    matrix_path = joinpath(temporary_directory, "interaction_matrix.f32")
    open(matrix_path, "w") do io
        write(io, interaction_matrix)
    end
    send_worker_command!(
        worker,
        "CONFIGURE $matrix_path $(Float64(gain))",
        "CONFIGURED";
        timeout=30.0,
    )
    send_worker_command!(worker, "FLATTEN", "FLATTENED")
    flat_command = zeros(Float32, actuator_count())
    read_next!(flat_command, streams.wfc; timeout=5.0)
    all(iszero, flat_command) || error(
        "pyRTC worker flatten command is nonzero",
    )
    return nothing
end

function controlled_disturbance!(prepared)
    graph = prepared.graph
    boundary = prepared.boundary
    disturbance_command = zeros(Float32, actuator_count())
    disturbance_command[8] = 3.0f-8
    disturbance_command[12] = -2.0f-8
    disturbance_command[14] = 2.0f-8
    disturbance_command[18] = -3.0f-8

    sequence = step_hil_frame!(boundary)
    copyto!(hil_command_buffer(boundary), disturbance_command)
    adopt_hil_command!(boundary, sequence)
    sequence = step_hil_frame!(boundary)
    disturbance_opd = copy(graph_output(graph, Val(:dm_surface_opd)))
    fill!(hil_command_buffer(boundary), 0.0f0)
    adopt_hil_command!(boundary, sequence)
    reset_hil_boundary!(boundary)
    copyto!(prepared.uncompensated_opd, disturbance_opd)
    return disturbance_command
end

function close_loop!(
    worker::PyRTCWorker,
    streams::ProcessStreams,
    prepared,
    signal::Vector{Float32};
    iterations::Int,
)
    boundary = prepared.boundary
    residual_norms = Vector{Float32}(undef, iterations)
    command = zeros(Float32, actuator_count())
    sequence = step_hil_frame!(boundary)

    for iteration in 1:iterations
        process_frame!(
            worker,
            streams,
            hil_frame_buffer(boundary),
            signal;
            command="STEP",
            response="STEPPED",
        )
        residual_norms[iteration] = norm(signal)
        read_next!(command, streams.wfc; timeout=5.0)
        copyto!(hil_command_buffer(boundary), command)
        adopt_hil_command!(boundary, sequence)
        iteration < iterations && (sequence = step_hil_frame!(boundary))
    end
    return residual_norms, command
end

function run_validation(
    pyrtc_root::AbstractString;
    wavefront_sensor::Symbol=:shack_hartmann,
)
    isdir(joinpath(pyrtc_root, "pyRTC")) || error(
        "expected a pyRTC source checkout containing pyRTC/ at '$pyrtc_root'",
    )
    case = process_reference_case(wavefront_sensor)
    prepared = prepare_hil_reference_system(case.wavefront_sensor)
    streams = create_process_streams(case)
    worker = nothing
    return mktempdir() do temporary_directory
        try
            worker = start_worker(pyrtc_root, case, temporary_directory)
            signal = zeros(Float32, case.signal_shape)
            flat_signal, interaction_matrix = calibrate_interaction_matrix!(
                worker,
                streams,
                prepared,
                signal;
                poke=case.poke,
            )
            norm(flat_signal) <= 1.0f-5 || error(
                "pyRTC process flat reference left a nonzero residual: " *
                "$(norm(flat_signal))",
            )
            all(isfinite, interaction_matrix) || error(
                "pyRTC process interaction matrix contains a non-finite value",
            )
            singular_values = svdvals(interaction_matrix)
            rank_tolerance = maximum(singular_values) * 1.0f-4
            interaction_rank = count(>(rank_tolerance), singular_values)
            interaction_rank == actuator_count() || error(
                "pyRTC process interaction matrix has rank $interaction_rank; " *
                "expected $(actuator_count())",
            )

            configure_worker_loop!(
                worker,
                streams,
                interaction_matrix,
                case.gain,
                temporary_directory,
            )
            disturbance_command = controlled_disturbance!(prepared)
            residual_norms, command = close_loop!(
                worker,
                streams,
                prepared,
                signal;
                iterations=case.iterations,
            )
            all(isfinite, residual_norms) || error(
                "pyRTC process loop produced a non-finite residual norm",
            )
            first(residual_norms) > 0 || error(
                "the controlled disturbance produced a zero initial residual",
            )
            convergence_ratio = last(residual_norms) / first(residual_norms)
            command_error = norm(command + disturbance_command)
            command_scale = norm(disturbance_command)
            convergence_ratio < case.convergence_limit || error(
                "pyRTC process loop did not converge: final/initial residual = " *
                "$convergence_ratio",
            )
            command_error <= max(5.0f-11, command_scale * 1.0f-3) || error(
                "pyRTC process command did not recover the disturbance: " *
                "error = $command_error",
            )
            stop_worker!(worker)

            return (
                wavefront_sensor=case.wavefront_sensor,
                signal_length=length(signal),
                interaction_rank,
                interaction_condition=maximum(singular_values) /
                    minimum(singular_values),
                initial_residual=first(residual_norms),
                final_residual=last(residual_norms),
                convergence_ratio,
                command_error,
            )
        finally
            !isnothing(worker) && stop_worker_noexcept!(worker)
            close_process_streams!(streams)
        end
    end
end

function main(wavefront_sensor::Symbol, args)
    pyrtc_root = if !isempty(args)
        first(args)
    else
        get(ENV, "PYRTC_ROOT", "")
    end
    isempty(pyrtc_root) && error(
        "pass the pyRTC checkout path or set PYRTC_ROOT",
    )
    result = run_validation(pyrtc_root; wavefront_sensor)
    println(
        "AOS/native-SHM pyRTC ",
        reference_label(Val(wavefront_sensor)),
        " process loop passed",
    )
    println("  signal length: ", result.signal_length)
    println("  interaction rank: ", result.interaction_rank)
    println("  interaction condition: ", result.interaction_condition)
    println("  initial residual: ", result.initial_residual)
    println("  final residual: ", result.final_residual)
    println("  convergence ratio: ", result.convergence_ratio)
    println("  command error: ", result.command_error)
    return nothing
end

end # module PyRTCProcessHIL
