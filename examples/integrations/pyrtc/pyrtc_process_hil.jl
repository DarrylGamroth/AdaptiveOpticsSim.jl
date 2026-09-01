module PyRTCProcessHIL

using AdaptiveOpticsSim.AlgorithmGraphs
using LinearAlgebra

include(joinpath(@__DIR__, "pyrtc_shared_memory.jl"))
using .PyRTCSharedMemory

include(joinpath(@__DIR__, "..", "..", "support", "hil_reference_systems.jl"))
using .HILReferenceSystems

const PYRTC_STREAM_NAMES = ("wfs", "wfc", "signal", "signal2D")
const VIEWER_STREAM_NAMES = (
    "wfc2D",
    "aosUncompensatedOpd",
    "aosDmSurfaceOpd",
    "aosResidualOpd",
    "aosOpenLoopPsf",
    "aosClosedLoopPsf",
)
const VIEWER_ITEMS = (
    "wfs",
    "signal2D",
    "wfc2D",
    "aosUncompensatedOpd",
    "aosDmSurfaceOpd",
    "aosResidualOpd",
    "aosOpenLoopPsf",
    "aosClosedLoopPsf",
)
const WORKER_PREFIX = "AOS_PYRTC_WORKER "
const REFERENCE_ATMOSPHERE_STEP_S = 1.0e-3
const ATMOSPHERE_VALIDATION_FRAMES = 20
const ATMOSPHERE_VALIDATION_BURN_IN_FRAMES = 10

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

struct ViewerStreams{C,U,D,R,O,L}
    command::C
    uncompensated_opd::U
    dm_surface_opd::D
    residual_opd::R
    open_loop_psf::O
    closed_loop_psf::L
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

function create_viewer_streams(prepared, science_diagnostics)
    occupied = String[]
    for name in VIEWER_STREAM_NAMES, suffix in ("", "_meta", "_gpu_handle")
        path = joinpath("/dev/shm", name * suffix)
        ispath(path) && push!(occupied, path)
    end
    isempty(occupied) || error(
        "refusing to reuse active viewer shared-memory streams: " *
        join(occupied, ", "),
    )
    actuator_axis_count = isqrt(actuator_count())
    actuator_axis_count^2 == actuator_count() || error(
        "the viewer requires a square deformable-mirror command layout",
    )
    command = nothing
    uncompensated_opd = nothing
    dm_surface_opd = nothing
    residual_opd = nothing
    open_loop_psf_stream = nothing
    closed_loop_psf_stream = nothing
    try
        graph = prepared.graph
        uncompensated_opd_values = graph_output(
            graph,
            Val(:atmosphere_opd),
        )
        dm_surface_opd_values = graph_output(graph, Val(:dm_surface_opd))
        residual_opd_values = graph_output(graph, Val(:pupil_opd))
        command = create_stream(
            VIEWER_STREAM_NAMES[1],
            Float32,
            (actuator_axis_count, actuator_axis_count),
        )
        uncompensated_opd = create_stream(
            VIEWER_STREAM_NAMES[2],
            eltype(uncompensated_opd_values),
            size(uncompensated_opd_values),
        )
        dm_surface_opd = create_stream(
            VIEWER_STREAM_NAMES[3],
            eltype(dm_surface_opd_values),
            size(dm_surface_opd_values),
        )
        residual_opd = create_stream(
            VIEWER_STREAM_NAMES[4],
            eltype(residual_opd_values),
            size(residual_opd_values),
        )
        open_loop_psf_values = open_loop_psf(science_diagnostics)
        closed_loop_psf_values = closed_loop_psf(science_diagnostics)
        open_loop_psf_stream = create_stream(
            VIEWER_STREAM_NAMES[5],
            eltype(open_loop_psf_values),
            size(open_loop_psf_values),
        )
        closed_loop_psf_stream = create_stream(
            VIEWER_STREAM_NAMES[6],
            eltype(closed_loop_psf_values),
            size(closed_loop_psf_values),
        )
        return ViewerStreams(
            command,
            uncompensated_opd,
            dm_surface_opd,
            residual_opd,
            open_loop_psf_stream,
            closed_loop_psf_stream,
        )
    catch
        !isnothing(closed_loop_psf_stream) &&
            close_and_unlink_noexcept!(closed_loop_psf_stream)
        !isnothing(open_loop_psf_stream) &&
            close_and_unlink_noexcept!(open_loop_psf_stream)
        !isnothing(residual_opd) && close_and_unlink_noexcept!(residual_opd)
        !isnothing(dm_surface_opd) &&
            close_and_unlink_noexcept!(dm_surface_opd)
        !isnothing(uncompensated_opd) &&
            close_and_unlink_noexcept!(uncompensated_opd)
        !isnothing(command) && close_and_unlink_noexcept!(command)
        rethrow()
    end
end

function close_viewer_streams!(streams::ViewerStreams)
    close_and_unlink_noexcept!(streams.closed_loop_psf)
    close_and_unlink_noexcept!(streams.open_loop_psf)
    close_and_unlink_noexcept!(streams.residual_opd)
    close_and_unlink_noexcept!(streams.dm_surface_opd)
    close_and_unlink_noexcept!(streams.uncompensated_opd)
    close_and_unlink_noexcept!(streams.command)
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
    wavefront_sensor::Symbol,
    temporary_directory::AbstractString,
)
    return Cmd(String[
        pyrtc_python(),
        joinpath(@__DIR__, "pyrtc_process_worker.py"),
        "--sensor",
        String(wavefront_sensor),
        "--temporary-directory",
        abspath(temporary_directory),
    ])
end

function viewer_refresh_rate(frame_rate::Real)
    requested = clamp(Float64(frame_rate) / 2, 1.0, 60.0)
    return floor(Int, requested)
end

function viewer_command(frame_rate::Real)
    executable = joinpath(dirname(pyrtc_python()), "pyrtc-view")
    isfile(executable) || error(
        "the selected pyRTC environment does not provide pyrtc-view; " *
        "install examples/integrations/pyrtc/requirements.txt",
    )
    command = Cmd(String[
        executable,
        VIEWER_ITEMS...,
        "--geometry",
        "2x4",
        "--fps",
        string(viewer_refresh_rate(frame_rate)),
        "--pixel-scale",
        "3",
    ])
    if !haskey(ENV, "QT_QPA_PLATFORM") &&
            !haskey(ENV, "DISPLAY") &&
            haskey(ENV, "WAYLAND_DISPLAY")
        command = addenv(command, "QT_QPA_PLATFORM" => "wayland")
    end
    return command
end

function start_viewer(frame_rate::Real)
    viewer = run(viewer_command(frame_rate); wait=false)
    sleep(0.5)
    if !process_running(viewer)
        wait(viewer)
        error("pyrtc-view exited during startup")
    end
    return viewer
end

function stop_viewer_noexcept!(viewer::Base.Process)
    process_running(viewer) && kill(viewer)
    try
        wait(viewer)
    catch
    end
    return nothing
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
    case::ProcessReferenceCase,
    temporary_directory::AbstractString,
)
    process = open(
        worker_command(
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

function publish_viewer_outputs!(
    streams::ViewerStreams,
    prepared,
    science_diagnostics,
    applied_command::Matrix{Float32},
)
    graph = prepared.graph
    publish!(streams.command, applied_command)
    publish!(
        streams.uncompensated_opd,
        graph_output(graph, Val(:atmosphere_opd)),
    )
    publish!(
        streams.dm_surface_opd,
        graph_output(graph, Val(:dm_surface_opd)),
    )
    publish!(
        streams.residual_opd,
        graph_output(graph, Val(:pupil_opd)),
    )
    publish!(streams.open_loop_psf, open_loop_psf(science_diagnostics))
    publish!(streams.closed_loop_psf, closed_loop_psf(science_diagnostics))
    return nothing
end

@inline function root_mean_square(values::AbstractArray)
    return sqrt(sum(abs2, values) / length(values))
end

function pupil_opd_rms(
    opd::AbstractMatrix{<:AbstractFloat},
    support::AbstractMatrix{Bool},
)
    axes(opd) == axes(support) || throw(DimensionMismatch(
        "the pupil OPD and support must have identical axes",
    ))
    sample_count = 0
    mean_opd = 0.0
    sum_squared_difference = 0.0
    for index in eachindex(opd, support)
        support[index] || continue
        sample_count += 1
        value = Float64(opd[index])
        difference = value - mean_opd
        mean_opd += difference / sample_count
        sum_squared_difference += difference * (value - mean_opd)
    end
    sample_count > 0 || error("the HIL reference pupil support is empty")
    return sqrt(sum_squared_difference / sample_count)
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

function mean_from(values::Vector{<:Real}, first_index::Int)
    first_index in eachindex(values) || throw(BoundsError(values, first_index))
    total = 0.0
    sample_count = 0
    for index in first_index:lastindex(values)
        total += Float64(values[index])
        sample_count += 1
    end
    return total / sample_count
end

function close_atmospheric_loop!(
    worker::PyRTCWorker,
    streams::ProcessStreams,
    case::ProcessReferenceCase,
    signal::Vector{Float32};
    frames::Int=ATMOSPHERE_VALIDATION_FRAMES,
)
    frames > ATMOSPHERE_VALIDATION_BURN_IN_FRAMES || throw(ArgumentError(
        "atmosphere validation requires more than " *
        "$(ATMOSPHERE_VALIDATION_BURN_IN_FRAMES) frames",
    ))
    prepared = prepare_atmospheric_hil_reference_system(
        case.wavefront_sensor;
        atmosphere_step=REFERENCE_ATMOSPHERE_STEP_S,
        rng_seed=1,
    )
    diagnostics = prepare_hil_science_diagnostics()
    boundary = prepared.boundary
    command = zeros(Float32, actuator_count())
    open_loop_values = Vector{Float32}(undef, frames)
    closed_loop_values = Vector{Float32}(undef, frames)
    sequence = step_hil_frame!(boundary)

    for frame_index in 1:frames
        atmosphere_opd = graph_output(prepared.graph, Val(:atmosphere_opd))
        residual_opd = graph_output(prepared.graph, Val(:pupil_opd))
        update_hil_science_diagnostics!(
            diagnostics,
            atmosphere_opd,
            residual_opd,
        )
        open_loop_values[frame_index] =
            open_loop_on_axis_strehl(diagnostics)
        closed_loop_values[frame_index] =
            closed_loop_on_axis_strehl(diagnostics)
        maximum(open_loop_psf(diagnostics)) <= 1.001f0 || error(
            "open-loop PSF exceeds its exact diffraction-limited peak",
        )
        maximum(closed_loop_psf(diagnostics)) <= 1.001f0 || error(
            "closed-loop PSF exceeds its exact diffraction-limited peak",
        )

        process_frame!(
            worker,
            streams,
            hil_frame_buffer(boundary),
            signal;
            command="STEP",
            response="STEPPED",
        )
        read_next!(command, streams.wfc; timeout=5.0)
        copyto!(hil_command_buffer(boundary), command)
        adopt_hil_command!(boundary, sequence)
        frame_index < frames && (sequence = step_hil_frame!(boundary))
    end

    first_steady_frame = ATMOSPHERE_VALIDATION_BURN_IN_FRAMES + 1
    mean_open_loop_on_axis_strehl =
        mean_from(open_loop_values, first_steady_frame)
    mean_closed_loop_on_axis_strehl =
        mean_from(closed_loop_values, first_steady_frame)
    return (;
        mean_open_loop_on_axis_strehl,
        mean_closed_loop_on_axis_strehl,
        improvement=mean_closed_loop_on_axis_strehl /
            mean_open_loop_on_axis_strehl,
    )
end

function run_validation(; wavefront_sensor::Symbol=:shack_hartmann)
    case = process_reference_case(wavefront_sensor)
    prepared = prepare_hil_reference_system(case.wavefront_sensor)
    streams = create_process_streams(case)
    worker = nothing
    return mktempdir() do temporary_directory
        try
            worker = start_worker(case, temporary_directory)
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

            configure_worker_loop!(
                worker,
                streams,
                interaction_matrix,
                case.gain,
                temporary_directory,
            )
            atmosphere = close_atmospheric_loop!(
                worker,
                streams,
                case,
                signal,
            )
            atmosphere.mean_closed_loop_on_axis_strehl > 0.5 || error(
                "pyRTC atmospheric loop did not produce a usable corrected " *
                "on-axis PSF: mean Strehl = " *
                "$(atmosphere.mean_closed_loop_on_axis_strehl)",
            )
            atmosphere.improvement > 10 || error(
                "pyRTC atmospheric loop did not improve the mean on-axis " *
                "Strehl sufficiently: ratio = $(atmosphere.improvement)",
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
                atmosphere...,
            )
        finally
            !isnothing(worker) && stop_worker_noexcept!(worker)
            close_process_streams!(streams)
        end
    end
end

function run_viewer_demo(;
    wavefront_sensor::Symbol=:pyramid,
    duration::Real=60.0,
    frame_rate::Real=10.0,
)
    isfinite(duration) && duration > 0 || throw(ArgumentError(
        "duration must be finite and positive",
    ))
    isfinite(frame_rate) && frame_rate > 0 || throw(ArgumentError(
        "frame_rate must be finite and positive",
    ))
    case = process_reference_case(wavefront_sensor)
    calibration = prepare_hil_reference_system(case.wavefront_sensor)
    streams = create_process_streams(case)
    viewer_streams = nothing
    worker = nothing
    viewer = nothing
    return mktempdir() do temporary_directory
        try
            worker = start_worker(case, temporary_directory)
            signal = zeros(Float32, case.signal_shape)
            command = zeros(Float32, actuator_count())
            actuator_axis_count = isqrt(actuator_count())
            applied_command = reshape(
                zeros(Float32, actuator_count()),
                actuator_axis_count,
                actuator_axis_count,
            )
            println("Calibrating the pyRTC loop against the flat reference...")
            _, interaction_matrix = calibrate_interaction_matrix!(
                worker,
                streams,
                calibration,
                signal;
                poke=case.poke,
            )
            configure_worker_loop!(
                worker,
                streams,
                interaction_matrix,
                case.gain,
                temporary_directory,
            )

            prepared = prepare_atmospheric_hil_reference_system(
                case.wavefront_sensor;
                atmosphere_step=REFERENCE_ATMOSPHERE_STEP_S,
                rng_seed=1,
            )
            science_diagnostics = prepare_hil_science_diagnostics()
            viewer_streams = create_viewer_streams(
                prepared,
                science_diagnostics,
            )
            boundary = prepared.boundary
            sequence = step_hil_frame!(boundary)
            atmosphere_opd = graph_output(
                prepared.graph,
                Val(:atmosphere_opd),
            )
            residual_opd = graph_output(prepared.graph, Val(:pupil_opd))
            update_hil_science_diagnostics!(
                science_diagnostics,
                atmosphere_opd,
                residual_opd,
            )
            publish_viewer_outputs!(
                viewer_streams,
                prepared,
                science_diagnostics,
                applied_command,
            )
            viewer = start_viewer(frame_rate)
            start_time = time()
            next_frame_time = start_time
            next_status_time = start_time
            frame_period = inv(Float64(frame_rate))

            println(
                "Running the ",
                reference_label(Val(wavefront_sensor)),
                " AOS/pyRTC live view for ",
                Float64(duration),
                " seconds; close the viewer to stop early.",
            )
            while time() - start_time < duration && process_running(viewer)
                publish_viewer_outputs!(
                    viewer_streams,
                    prepared,
                    science_diagnostics,
                    applied_command,
                )
                process_frame!(
                    worker,
                    streams,
                    hil_frame_buffer(boundary),
                    signal;
                    command="STEP",
                    response="STEPPED",
                )
                read_next!(command, streams.wfc; timeout=5.0)

                copyto!(hil_command_buffer(boundary), command)
                adopt_hil_command!(boundary, sequence)

                elapsed = time() - start_time

                if time() >= next_status_time
                    println(
                        "  t=", round(elapsed; digits=1),
                        " s, open-loop OPD=", round(
                            1.0e9 * pupil_opd_rms(
                                atmosphere_opd,
                                HILReferenceSystems.pupil_support(
                                    science_diagnostics,
                                ),
                            );
                            digits=3,
                        ),
                        " nm RMS, residual OPD=", round(
                            1.0e9 * pupil_opd_rms(
                                residual_opd,
                                HILReferenceSystems.pupil_support(
                                    science_diagnostics,
                                ),
                            );
                            digits=3,
                        ),
                        " nm RMS, on-axis Strehl=", round(
                            open_loop_on_axis_strehl(science_diagnostics);
                            digits=4,
                        ),
                        " -> ", round(
                            closed_loop_on_axis_strehl(science_diagnostics);
                            digits=4,
                        ),
                        ", command=", round(
                            1.0e9 * root_mean_square(applied_command);
                            digits=3,
                        ),
                        " nm RMS",
                    )
                    next_status_time += 1.0
                end

                next_frame_time += frame_period
                delay = next_frame_time - time()
                delay > 0 && sleep(delay)
                (time() - start_time >= duration ||
                    !process_running(viewer)) && break
                sequence = step_hil_frame!(boundary)
                copyto!(vec(applied_command), command)
                atmosphere_opd = graph_output(
                    prepared.graph,
                    Val(:atmosphere_opd),
                )
                residual_opd = graph_output(
                    prepared.graph,
                    Val(:pupil_opd),
                )
                update_hil_science_diagnostics!(
                    science_diagnostics,
                    atmosphere_opd,
                    residual_opd,
                )
            end
            return nothing
        finally
            !isnothing(viewer) && stop_viewer_noexcept!(viewer)
            !isnothing(worker) && stop_worker_noexcept!(worker)
            !isnothing(viewer_streams) &&
                close_viewer_streams!(viewer_streams)
            close_process_streams!(streams)
        end
    end
end

function main(wavefront_sensor::Symbol)
    result = run_validation(; wavefront_sensor)
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
    println(
        "  atmospheric mean open-loop on-axis Strehl: ",
        result.mean_open_loop_on_axis_strehl,
    )
    println(
        "  atmospheric mean closed-loop on-axis Strehl: ",
        result.mean_closed_loop_on_axis_strehl,
    )
    println("  atmospheric on-axis Strehl improvement: ", result.improvement)
    return nothing
end

function viewer_main(
    wavefront_sensor::Symbol;
    duration::Real=60.0,
    frame_rate::Real=10.0,
)
    run_viewer_demo(; wavefront_sensor, duration, frame_rate)
    return nothing
end

end # module PyRTCProcessHIL
