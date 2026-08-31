module PyRTCReferenceHIL

using AdaptiveOpticsSim.AlgorithmGraphs
using LinearAlgebra
using PythonCall

include(joinpath(@__DIR__, "..", "..", "support", "hil_reference_systems.jl"))
using .HILReferenceSystems

const PYRTC_STREAM_NAMES = ("wfs", "wfc", "signal", "signal2D")
const PYRAMID_PUPIL_SAMPLES = 16
const PYRAMID_PUPIL_SEPARATION_PIXELS = 2
const PYRAMID_FRAME_EDGE_PIXELS = 1
const PYRAMID_FRAME_EXTENT = 2 * PYRAMID_PUPIL_SAMPLES +
    PYRAMID_PUPIL_SEPARATION_PIXELS + 2 * PYRAMID_FRAME_EDGE_PIXELS

struct PyRTCReferenceCase
    wavefront_sensor::Symbol
    frame_shape::Tuple{Int,Int}
    expected_signal_length::Int
    poke::Float32
    gain::Float32
    iterations::Int
    convergence_limit::Float32
end

@inline reference_case(::Val{:shack_hartmann}) = PyRTCReferenceCase(
    :shack_hartmann,
    (64, 64),
    104,
    2.0f-8,
    0.4f0,
    15,
    1.0f-3,
)

@inline reference_case(::Val{:pyramid}) = PyRTCReferenceCase(
    :pyramid,
    (PYRAMID_FRAME_EXTENT, PYRAMID_FRAME_EXTENT),
    344,
    2.0f-8,
    0.4f0,
    15,
    1.0f-3,
)

function reference_case(wavefront_sensor::Symbol)
    wavefront_sensor in supported_wavefront_sensors() || throw(ArgumentError(
        "unsupported pyRTC reference WFS '$wavefront_sensor'; expected one " *
        "of $(join(supported_wavefront_sensors(), ", "))",
    ))
    return reference_case(Val(wavefront_sensor))
end

@inline reference_label(::Val{:shack_hartmann}) = "Shack-Hartmann"
@inline reference_label(::Val{:pyramid}) = "Pyramid"

@inline function pyramid_pupil_centres()
    lower = PYRAMID_FRAME_EDGE_PIXELS + div(PYRAMID_PUPIL_SAMPLES, 2)
    upper = lower + PYRAMID_PUPIL_SAMPLES +
        PYRAMID_PUPIL_SEPARATION_PIXELS
    return (
        (lower, lower),
        (upper, lower),
        (lower, upper),
        (upper, upper),
    )
end

struct PyRTCReferenceComponents
    numpy::Py
    shared_memory::Py
    wfs::Py
    wavefront_corrector::Py
    slopes::Py
    loop::Py
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

function clear_owned_streams!(shared_memory::Py)
    for name in PYRTC_STREAM_NAMES, suffix in ("", "_meta", "_gpu_handle")
        stream_name = name * suffix
        ispath(joinpath("/dev/shm", stream_name)) || continue
        handle = shared_memory.SharedMemory(; name=stream_name)
        try
            handle.unlink()
        finally
            handle.close()
        end
    end
    return nothing
end

function close_component_handles!(components::PyRTCReferenceComponents)
    components.loop.stop()
    components.slopes.stop()
    components.wavefront_corrector.stop()
    image_handles = (
        components.loop.signalShm,
        components.loop.wfcShm,
        components.slopes.wfsShm,
        components.slopes.signal,
        components.slopes.signal2D,
        components.wavefront_corrector.correctionVector,
        components.wfs,
    )
    for handle in image_handles
        handle.close()
        handle.metadataShm.close()
    end
    return nothing
end

function python_config(pybuiltins::Py, entries::Pair...)
    config = pybuiltins.dict()
    for entry in entries
        config[entry.first] = entry.second
    end
    return config
end

function prepare_slopes_config(
    ::Val{:shack_hartmann},
    pybuiltins::Py,
    numpy::Py,
    temporary_directory::AbstractString,
)
    valid_subapertures = shack_hartmann_valid_subapertures()
    numpy_mask = numpy.array(
        Py(valid_subapertures);
        dtype=numpy.bool_,
        order="C",
        copy=true,
    )
    valid_subapertures_path =
        joinpath(temporary_directory, "valid_subapertures.npy")
    numpy.save(
        valid_subapertures_path,
        numpy.concatenate((numpy_mask, numpy_mask); axis=0),
    )
    return python_config(
        pybuiltins,
        "type" => "SHWFS",
        "signalType" => "slopes",
        "subApSpacing" => 8,
        "subApOffsetX" => 0,
        "subApOffsetY" => 0,
        "imageNoise" => 0.0,
        "contrast" => 0,
        "validSubApsFile" => valid_subapertures_path,
    )
end

function prepare_slopes_config(
    ::Val{:pyramid},
    pybuiltins::Py,
    ::Py,
    ::AbstractString,
)
    pupil_centres = pybuiltins.list()
    for (row, column) in pyramid_pupil_centres()
        pupil_centres.append("$row,$column")
    end
    return python_config(
        pybuiltins,
        "type" => "PYWFS",
        "signalType" => "slopes",
        "imageNoise" => 0.0,
        "centralObscurationRatio" => 0.0,
        "flatNorm" => true,
        "pupils" => pupil_centres,
        "pupilsRadius" => div(PYRAMID_PUPIL_SAMPLES, 2),
    )
end

function prepare_pyrtc_components(
    pyrtc_root::AbstractString,
    case::PyRTCReferenceCase,
)
    isdir(joinpath(pyrtc_root, "pyRTC")) || error(
        "expected a pyRTC source checkout containing pyRTC/ at '$pyrtc_root'",
    )
    sys = pyimport("sys")
    sys.path.insert(0, abspath(pyrtc_root))
    numpy = pyimport("numpy")
    pybuiltins = pyimport("builtins")
    shared_memory = pyimport("multiprocessing.shared_memory")
    ImageSHM = pyimport("pyRTC.Pipeline").ImageSHM
    SlopesProcess = pyimport("pyRTC.SlopesProcess").SlopesProcess
    WavefrontCorrector =
        pyimport("pyRTC.WavefrontCorrector").WavefrontCorrector
    Loop = pyimport("pyRTC.Loop").Loop

    return mktempdir() do temporary_directory
        slopes_config = prepare_slopes_config(
            Val(case.wavefront_sensor),
            pybuiltins,
            numpy,
            temporary_directory,
        )

        wfs = ImageSHM(
            "wfs",
            case.frame_shape,
            numpy.float32;
            consumer=false,
        )
        wavefront_corrector = WavefrontCorrector(python_config(
            pybuiltins,
            "name" => "AOSReferenceDM",
            "numActuators" => actuator_count(),
            "numModes" => actuator_count(),
        ))
        slopes = SlopesProcess(slopes_config)
        loop = Loop(python_config(
            pybuiltins,
            "gain" => Float64(case.gain),
            "numDroppedModes" => 0,
        ))
        return PyRTCReferenceComponents(
            numpy,
            shared_memory,
            wfs,
            wavefront_corrector,
            slopes,
            loop,
        )
    end
end

function validate_sensor_geometry!(
    ::Val{:shack_hartmann},
    components::PyRTCReferenceComponents,
)
    valid_shape = pyconvert(
        Tuple{Int,Int},
        components.slopes.validSubAps.shape,
    )
    valid_shape == (16, 8) || error(
        "pyRTC SHWFS valid-signal mask has shape $valid_shape; expected " *
        "(16, 8)",
    )
    return nothing
end

function validate_sensor_geometry!(
    ::Val{:pyramid},
    components::PyRTCReferenceComponents,
)
    expected_locations = Tuple{Int,Int}[
        (column, row) for (row, column) in pyramid_pupil_centres()
    ]
    locations = pyconvert(
        Vector{Tuple{Int,Int}},
        components.slopes.pupilLocs,
    )
    locations == expected_locations || error(
        "pyRTC Pyramid pupil locations are $locations; expected " *
        "$expected_locations",
    )
    pupil_radius = pyconvert(Int, components.slopes.pupilRadius)
    pupil_radius == div(PYRAMID_PUPIL_SAMPLES, 2) || error(
        "pyRTC Pyramid pupil radius is $pupil_radius; expected " *
        "$(div(PYRAMID_PUPIL_SAMPLES, 2))",
    )
    valid_shape = pyconvert(
        Tuple{Int,Int},
        components.slopes.validSubAps.shape,
    )
    valid_shape == (PYRAMID_PUPIL_SAMPLES, 2 * PYRAMID_PUPIL_SAMPLES) ||
        error(
            "pyRTC Pyramid valid-signal mask has shape $valid_shape; " *
            "expected ($(PYRAMID_PUPIL_SAMPLES), " *
            "$(2 * PYRAMID_PUPIL_SAMPLES))",
        )
    return nothing
end

function validate_component_geometry!(
    components::PyRTCReferenceComponents,
    case::PyRTCReferenceCase,
)
    image_shape = Tuple(pyconvert(Vector{Int}, components.slopes.imageShape))
    image_shape == case.frame_shape || error(
        "pyRTC WFS image shape is $image_shape; expected $(case.frame_shape)",
    )
    validate_sensor_geometry!(Val(case.wavefront_sensor), components)
    return nothing
end

function publish_frame_and_compute_signal!(
    components::PyRTCReferenceComponents,
    frame::AbstractMatrix{Float32},
)
    python_frame = components.numpy.array(
        Py(frame);
        dtype=components.numpy.float32,
        order="C",
        copy=true,
    )
    write_status = pyconvert(Int, components.wfs.write(python_frame))
    write_status == 1 || error("pyRTC rejected the AOS WFS frame")
    components.slopes.computeSignal()
    return pyconvert(
        Vector{Float32},
        components.slopes.signal.read_noblock(),
    )
end

function set_flat_reference!(
    components::PyRTCReferenceComponents,
    boundary,
)
    sequence = step_hil_frame!(boundary)
    publish_frame_and_compute_signal!(components, hil_frame_buffer(boundary))
    components.slopes.setRefSlopes(
        components.slopes.signal2D.read_noblock(),
    )
    flat_signal = publish_frame_and_compute_signal!(
        components,
        hil_frame_buffer(boundary),
    )
    return sequence, flat_signal
end

function calibrate_interaction_matrix!(
    components::PyRTCReferenceComponents,
    prepared;
    poke::Float32,
)
    boundary = prepared.boundary
    sequence, flat_signal = set_flat_reference!(components, boundary)
    interaction_matrix = Matrix{Float32}(
        undef,
        length(flat_signal),
        actuator_count(),
    )
    positive_signal = similar(flat_signal)

    for command_index in axes(interaction_matrix, 2)
        fill!(hil_command_buffer(boundary), 0.0f0)
        hil_command_buffer(boundary)[command_index] = poke
        adopt_hil_command!(boundary, sequence)
        sequence = step_hil_frame!(boundary)
        copyto!(
            positive_signal,
            publish_frame_and_compute_signal!(
                components,
                hil_frame_buffer(boundary),
            ),
        )

        fill!(hil_command_buffer(boundary), 0.0f0)
        hil_command_buffer(boundary)[command_index] = -poke
        adopt_hil_command!(boundary, sequence)
        sequence = step_hil_frame!(boundary)
        negative_signal = publish_frame_and_compute_signal!(
            components,
            hil_frame_buffer(boundary),
        )
        @views @. interaction_matrix[:, command_index] =
            (positive_signal - negative_signal) / (2 * poke)
    end

    fill!(hil_command_buffer(boundary), 0.0f0)
    adopt_hil_command!(boundary, sequence)
    reset_hil_boundary!(boundary)
    return flat_signal, interaction_matrix
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

function configure_pyrtc_loop!(
    components::PyRTCReferenceComponents,
    interaction_matrix::AbstractMatrix{Float32},
    gain::Float32,
)
    components.loop.IM = components.numpy.array(
        Py(interaction_matrix);
        dtype=components.numpy.float32,
        order="C",
        copy=true,
    )
    components.loop.computeCM()
    components.loop.setGain(gain)
    components.loop.flatten()
    return nothing
end

function close_loop!(
    components::PyRTCReferenceComponents,
    prepared;
    iterations::Int,
)
    iterations > 1 || throw(ArgumentError("iterations must be greater than one"))
    boundary = prepared.boundary
    residual_norms = Vector{Float32}(undef, iterations)
    command = zeros(Float32, actuator_count())
    sequence = step_hil_frame!(boundary)

    for iteration in 1:iterations
        signal = publish_frame_and_compute_signal!(
            components,
            hil_frame_buffer(boundary),
        )
        residual_norms[iteration] = norm(signal)
        components.loop.standardIntegrator()
        copyto!(
            command,
            pyconvert(
                Vector{Float32},
                components.wavefront_corrector.read(),
            ),
        )
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
    case = reference_case(wavefront_sensor)
    prepared = prepare_hil_reference_system(case.wavefront_sensor)
    require_available_stream_names()
    components = nothing
    shared_memory = pyimport("multiprocessing.shared_memory")
    try
        components = prepare_pyrtc_components(pyrtc_root, case)
        validate_component_geometry!(components, case)
        flat_signal, interaction_matrix = calibrate_interaction_matrix!(
            components,
            prepared;
            poke=case.poke,
        )
        length(flat_signal) == case.expected_signal_length || error(
            "expected $(case.expected_signal_length) valid pyRTC " *
            "$(reference_label(Val(case.wavefront_sensor))) signals, got " *
            "$(length(flat_signal))",
        )
        norm(flat_signal) <= 1.0f-5 || error(
            "pyRTC flat reference left a nonzero residual: " *
            "$(norm(flat_signal))",
        )
        all(isfinite, interaction_matrix) || error(
            "pyRTC interaction matrix contains a non-finite value",
        )
        singular_values = svdvals(interaction_matrix)
        rank_tolerance = maximum(singular_values) * 1.0f-4
        interaction_rank = count(>(rank_tolerance), singular_values)
        interaction_rank == actuator_count() || error(
            "pyRTC interaction matrix has rank $interaction_rank; expected " *
            "$(actuator_count())",
        )

        configure_pyrtc_loop!(components, interaction_matrix, case.gain)
        disturbance_command = controlled_disturbance!(prepared)
        residual_norms, command = close_loop!(
            components,
            prepared;
            iterations=case.iterations,
        )
        all(isfinite, residual_norms) || error(
            "pyRTC loop produced a non-finite residual norm",
        )
        first(residual_norms) > 0 || error(
            "the controlled disturbance produced a zero initial residual",
        )
        convergence_ratio = last(residual_norms) / first(residual_norms)
        command_error = norm(command + disturbance_command)
        command_scale = norm(disturbance_command)
        convergence_ratio < case.convergence_limit || error(
            "pyRTC loop did not converge: final/initial residual = " *
            "$convergence_ratio",
        )
        command_error <= max(5.0f-11, command_scale * 1.0f-3) || error(
            "pyRTC command did not recover the controlled disturbance: " *
            "error = $command_error",
        )

        return (
            wavefront_sensor=case.wavefront_sensor,
            signal_length=length(flat_signal),
            interaction_rank,
            interaction_condition=maximum(singular_values) /
                minimum(singular_values),
            initial_residual=first(residual_norms),
            final_residual=last(residual_norms),
            convergence_ratio,
            command_error,
        )
    finally
        try
            if !isnothing(components)
                close_component_handles!(components)
                shared_memory = components.shared_memory
            end
        finally
            clear_owned_streams!(shared_memory)
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
        "AOS/pyRTC ",
        reference_label(Val(wavefront_sensor)),
        " reference loop passed",
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

end # module PyRTCReferenceHIL
