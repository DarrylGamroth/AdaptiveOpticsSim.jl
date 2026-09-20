"""
    profile_aos_gpu_plant.jl BACKEND [POLICY]

Profile the AOS physical Shack--Hartmann plant only:

    PDM command -> DM surface -> pupil OPD -> photon-rate mosaic -> CCD frame

`BACKEND` is `cuda` or `amdgpu`; `POLICY` is `stream`, `captured`, or `both`
(the default).  The detector is deliberately noiseless, so a host AOS graph is
an independent deterministic oracle.  No FilterGraphAlgorithms, reconstructor,
controller, or command adoption is constructed here: S1 remains a CPU bridge.

The profile region contains `AOS_PROFILE_REPLAYS` complete, synchronized graph
steps (default 1000).  Every policy receives exactly 32 cold-excluded warmups.
The checked-stream path reports its warmed Julia host-launch allocation; native
captured replay requires exactly zero Julia bytes, including completion.
The named ranges let an external tool select just the repeated plant work:

    nsys profile --trace=cuda,nvtx --capture-range=nvtx \
      --stop-on-range-end=true -- ./julia-command cuda captured

    rocprofv3 --runtime-trace --marker-trace --selected-regions \
      --output-directory rocprof-aos -- ./julia-command amdgpu captured

Use the maintained backend environments, for example:

    julia --project=test/cuda --startup-file=no \
      benchmarks/profile_aos_gpu_plant.jl cuda captured
    AOS_PROFILE_REPLAYS=2 julia --project=test/amdgpu --startup-file=no \
      benchmarks/profile_aos_gpu_plant.jl amdgpu both

External collection and PMU replay perturb timing.  The reported allocation
and device-memory gates run outside the marked measurement region.
"""

using AdaptiveOpticsSim
using AdaptiveOpticsSim.AlgorithmGraphs
using AdaptiveOpticsSim.Backends
using Profile
using Test: @inferred

const _PROFILE_BACKEND = isempty(ARGS) ? "" : lowercase(ARGS[1])
const _PROFILE_POLICY = length(ARGS) < 2 ? "both" : lowercase(ARGS[2])
const _PROFILE_REPLAYS = parse(Int, get(ENV, "AOS_PROFILE_REPLAYS", "1000"))
const _PROFILE_ALLOCATION_DIAGNOSTICS = get(
    ENV,
    "AOS_PROFILE_ALLOCATION_DIAGNOSTICS",
    "0",
) == "1"
const _PROFILE_WARMUPS = 32
const _PUPIL_RESOLUTION = 352
const _LENSLETS_PER_AXIS = 16
const _PIXELS_PER_LENSLET = 22
const _ACTUATORS_PER_AXIS = 5
const _ACTUATOR_COUNT = _ACTUATORS_PER_AXIS^2

_PROFILE_REPLAYS > 0 || error("AOS_PROFILE_REPLAYS must be positive")
_PROFILE_POLICY in ("stream", "captured", "both") || error(
    "POLICY must be stream, captured, or both",
)

if _PROFILE_BACKEND == "cuda"
    import CUDA
elseif _PROFILE_BACKEND == "amdgpu"
    import AMDGPU
else
    error("BACKEND must be cuda or amdgpu")
end

"""Return retained accelerator storage and its exact AOS compute-device target."""
function _backend_storage(backend::String)
    if backend == "cuda"
        CUDA.functional() || error("CUDA.jl has no functional CUDA device")
        Backends.disable_scalar_backend!(Backends.CUDABackendTag)
        return CUDA.CuArray, CUDABackend(), Backends.CUDABackendTag
    end
    AMDGPU.functional() || error("AMDGPU.jl has no functional ROCm device")
    Backends.disable_scalar_backend!(Backends.AMDGPUBackendTag)
    return AMDGPU.ROCArray, AMDGPUBackend(), Backends.AMDGPUBackendTag
end

function _actuator_coordinates()
    coordinates = Matrix{Float32}(undef, 2, _ACTUATOR_COUNT)
    axis = range(-0.8f0, 0.8f0; length=_ACTUATORS_PER_AXIS)
    index = 1
    @inbounds for y in axis, x in axis
        coordinates[1, index] = x
        coordinates[2, index] = y
        index += 1
    end
    return coordinates
end

function _pdm_command()
    command = Vector{Float32}(undef, _ACTUATOR_COUNT)
    @inbounds for index in eachindex(command)
        command[index] = 2.5f-8 * sinpi(Float32(index) / Float32(_ACTUATOR_COUNT))
    end
    return command
end

function _plant_definition(command, uncompensated_opd, coordinates; name::Symbol)
    command_schema = "org.adaptiveopticssim.profile.pdm-command.f32/1"
    surface_schema = "org.adaptiveopticssim.profile.dm-surface-opd.f32/1"
    pupil_schema = "org.adaptiveopticssim.profile.pupil-opd.f32/1"
    rate_schema = "org.adaptiveopticssim.profile.shwfs-photon-rate.f32/1"
    frame_schema = "org.adaptiveopticssim.profile.shwfs-electron-frame.f32/1"
    return algorithm_graph(
        (
            gaussian_deformable_mirror_surface_node(
                :dm;
                resolution=_PUPIL_RESOLUTION,
                telescope_diameter_m=8.0f0,
                actuator_count=_ACTUATOR_COUNT,
                influence_width=0.18f0,
                pdm_command_schema=command_schema,
                surface_opd_schema=surface_schema,
                actuator_coordinates_schema=
                    "org.adaptiveopticssim.profile.dm-coordinates.f32/1",
            ),
            pupil_opd_composition_node(
                :pupil;
                resolution=_PUPIL_RESOLUTION,
                uncompensated_opd_schema=pupil_schema,
                surface_opd_schema=surface_schema,
                pupil_opd_schema=pupil_schema,
            ),
            shack_hartmann_rate_node(
                :shwfs;
                resolution=_PUPIL_RESOLUTION,
                telescope_diameter_m=8.0f0,
                n_lenslets=_LENSLETS_PER_AXIS,
                n_pix_subap=_PIXELS_PER_LENSLET,
                pixel_scale_arcsec=0.2f0,
                source_wavelength_m=750.0f-9,
                source_photon_irradiance_m2_s=2.0f6,
                opd_schema=pupil_schema,
                photon_rate_schema=rate_schema,
            ),
            ccd_detector_acquisition_node(
                :detector;
                rows=_PUPIL_RESOLUTION,
                columns=_PUPIL_RESOLUTION,
                pixel_scale_arcsec=0.2f0,
                wavelength_m=750.0f-9,
                exposure_duration_s=1.0f-3,
                quantum_efficiency=0.8f0,
                photon_noise=false,
                readout_noise=false,
                rng_seed=0xA05,
                photon_rate_schema=rate_schema,
                frame_schema=frame_schema,
            ),
        );
        name,
        inputs=(
            graph_input(:pdm_command, :dm => :pdm_command, command),
            graph_input(
                :uncompensated_opd,
                :pupil => :uncompensated_opd,
                uncompensated_opd,
            ),
        ),
        outputs=(
            graph_output(:surface_opd, :dm => :surface_opd),
            graph_output(:pupil_opd, :pupil => :pupil_opd),
            graph_output(:photon_rate, :shwfs => :photon_rate),
            graph_output(:frame, :detector => :frame),
        ),
        links=(
            link(:dm => :surface_opd, :pupil => :surface_opd),
            link(:pupil => :pupil_opd, :shwfs => :opd),
            link(:shwfs => :photon_rate, :detector => :photon_rate),
        ),
        parameters=(sparse_parameter(:dm => :actuator_coordinates, coordinates),),
    )
end

@inline function _complete_step!(graph)
    wait_graph_step!(step_graph_async!(graph))
    return nothing
end

@inline function _replay!(graph, replays::Int)
    @inbounds for _ in 1:replays
        _complete_step!(graph)
    end
    return nothing
end

function _print_warmed_allocation_profile!(graph, label::AbstractString)
    Profile.Allocs.clear()
    Profile.Allocs.start(sample_rate=1.0)
    try
        _complete_step!(graph)
    finally
        Profile.Allocs.stop()
    end
    records = Profile.Allocs.fetch()
    println("warmed_allocation_profile_begin=$label")
    Profile.Allocs.print(stdout, records)
    println("warmed_allocation_profile_end=$label")
    return nothing
end

function _warm_and_measure_allocation!(
    graph,
    label::AbstractString;
    require_zero::Bool,
)
    _replay!(graph, _PROFILE_WARMUPS)
    _PROFILE_ALLOCATION_DIAGNOSTICS &&
        _print_warmed_allocation_profile!(graph, label)
    @inferred _complete_step!(graph)
    allocated = @allocated _complete_step!(graph)
    require_zero && allocated != 0 && error(
        "$label allocated $allocated warmed Julia bytes including completion",
    )
    return allocated
end

function _assert_reference_frame!(reference, graph, label::AbstractString)
    frame = Array(graph_output(graph, Val(:frame)))
    all(isfinite, frame) || error("$label produced a non-finite detector frame")
    sum(frame) > 0 || error("$label produced an empty detector frame")
    isapprox(frame, reference; rtol=5.0f-4, atol=1.0f-5) || error(
        "$label differs from the CPU AOS reference frame",
    )
    return nothing
end

_device_used_bytes(::Type{Backends.CUDABackendTag}) = Int(CUDA.used_memory())
_device_used_bytes(::Type{Backends.AMDGPUBackendTag}) = Int(AMDGPU.used_memory())

# These ranges are intentionally external-tool ABI calls rather than a
# dependency on a profiling package.  Nsight Systems records NVTX; rocprofv3
# records ROCTx.  They surround only the repeated selected plant region.
@inline function _profile_region!(
    ::Type{Backends.CUDABackendTag},
    label::String,
    f::F,
) where {F<:Function}
    ccall((:nvtxRangePushA, "libnvToolsExt"), Cint, (Cstring,), label)
    try
        return f()
    finally
        ccall((:nvtxRangePop, "libnvToolsExt"), Cint, ())
    end
end

@inline function _profile_region!(
    ::Type{Backends.AMDGPUBackendTag},
    label::String,
    f::F,
) where {F<:Function}
    ccall((:roctxRangePushA, "libroctx64"), Cint, (Cstring,), label)
    try
        return f()
    finally
        ccall((:roctxRangePop, "libroctx64"), Cint, ())
    end
end

function _profile_policy!(policy::Symbol, definition, target, backend_tag, reference)
    execution = policy === :stream ? StreamGraphExecution() : CapturedGraphExecution()
    println("stage=$(policy).prepare")
    graph = prepare_algorithm_graph(definition; target, execution)
    expected_captured_nodes = policy === :captured ? 4 : 0
    captured_graph_node_count(graph) == expected_captured_nodes || error(
        "$(policy) graph captured $(captured_graph_node_count(graph)) nodes; expected $expected_captured_nodes",
    )

    println("stage=$(policy).warm_infer_allocate")
    host_allocation_bytes = _warm_and_measure_allocation!(
        graph,
        string(policy);
        require_zero=(policy === :captured),
    )
    println("stage=$(policy).pre_profile_oracle")
    _assert_reference_frame!(reference, graph, "$(policy) pre-profile")

    println("stage=$(policy).device_memory_before")
    device_used_before = _device_used_bytes(backend_tag)
    region_name = "aos_plant_" *
        (policy === :stream ? "checked_stream" : "captured_graph")
    println("stage=$(policy).selected_region")
    start_ns = time_ns()
    _profile_region!(
        backend_tag,
        region_name,
        () -> _replay!(graph, _PROFILE_REPLAYS),
    )
    elapsed_ns = time_ns() - start_ns
    println("stage=$(policy).device_memory_after")
    device_used_after = _device_used_bytes(backend_tag)
    device_used_after == device_used_before || error(
        "$(policy) changed used device memory from $device_used_before to $device_used_after bytes",
    )
    println("stage=$(policy).post_profile_oracle")
    _assert_reference_frame!(reference, graph, "$(policy) post-profile")
    println("policy=$(policy) captured_nodes=$(captured_graph_node_count(graph)) " *
        "replays=$_PROFILE_REPLAYS total_ns=$elapsed_ns " *
        "mean_ns=$(elapsed_ns / _PROFILE_REPLAYS) " *
        "warmed_host_allocation_bytes=$host_allocation_bytes " *
        "device_used_bytes=$device_used_after")
    return nothing
end

function main()
    BackendArray, _, backend_tag = _backend_storage(_PROFILE_BACKEND)
    command = _pdm_command()
    coordinates = _actuator_coordinates()
    uncompensated_opd = zeros(Float32, _PUPIL_RESOLUTION, _PUPIL_RESOLUTION)

    println("stage=cpu_reference.prepare")
    cpu_graph = prepare_algorithm_graph(
        _plant_definition(command, uncompensated_opd, coordinates; name=:aos_plant_cpu_reference),
    )
    println("stage=cpu_reference.warm_infer_allocate")
    _replay!(cpu_graph, _PROFILE_WARMUPS)
    @inferred _complete_step!(cpu_graph)
    cpu_allocated = @allocated _complete_step!(cpu_graph)
    cpu_allocated == 0 || error("CPU reference allocated $cpu_allocated warmed Julia bytes")
    println("stage=cpu_reference.before")
    reference_before = copy(graph_output(cpu_graph, Val(:frame)))
    all(isfinite, reference_before) && sum(reference_before) > 0 || error(
        "CPU reference did not produce a finite nonempty detector frame",
    )

    gpu_command = BackendArray(command)
    gpu_uncompensated = BackendArray(uncompensated_opd)
    gpu_coordinates = BackendArray(coordinates)
    target = compute_device(gpu_command)
    definition = _plant_definition(
        gpu_command,
        gpu_uncompensated,
        gpu_coordinates;
        name=:aos_gpu_physical_plant,
    )

    println("backend=$_PROFILE_BACKEND frame_shape=$(_PUPIL_RESOLUTION)x$(_PUPIL_RESOLUTION) " *
        "warmups=$_PROFILE_WARMUPS replays=$_PROFILE_REPLAYS")
    _PROFILE_POLICY in ("stream", "both") &&
        _profile_policy!(:stream, definition, target, backend_tag, reference_before)
    _PROFILE_POLICY in ("captured", "both") &&
        _profile_policy!(:captured, definition, target, backend_tag, reference_before)

    println("stage=cpu_reference.after")
    reset_graph!(cpu_graph)
    _complete_step!(cpu_graph)
    reference_after = graph_output(cpu_graph, Val(:frame))
    reference_after == reference_before || error("CPU reference changed across the profile run")
    println("cpu_reference_before_after=passed")
    return nothing
end

main()
