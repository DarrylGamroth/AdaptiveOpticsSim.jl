module SubaruAO188Simulation

using AdaptiveOpticsSim
using LinearAlgebra
using Random
using Statistics

import AdaptiveOpticsSim: step!, runtime_timing, convert_noise, validate_noise, materialize_build,
    execution_style, synchronize_backend!, bin2d!, apply_command!, prepare_sampling!,
    ensure_sh_calibration!, wfs_output_frame, prepare!, prepare_runtime_wfs!, supports_prepared_runtime,
    supports_detector_output, supports_grouped_execution, simulation_interface, wfs_output_metadata,
    init_execution_state

export AO188ActuatorSupportModel, CircularActuatorSupport
export SubaruHighOrderWFSModel, OperationalShackHartmannModel, AO188CurvatureModel
export AO188ReplayMode, DirectReplayMode, PreparedReplayMode
export AO188LatencyModel, AO188DetectorConfig, AO188WFSDetectorConfig, AO188APDDetectorConfig
export AO188SimulationParams, AO188CurvatureSimulationParams
export AO188Simulation, subaru_ao188_simulation, subaru_ao188_curvature_simulation
export subaru_ao188_phase_timing, prepare_replay!

abstract type AO188ActuatorSupportModel end
struct CircularActuatorSupport <: AO188ActuatorSupportModel end

abstract type SubaruHighOrderWFSModel end
struct OperationalShackHartmannModel <: SubaruHighOrderWFSModel end
struct AO188CurvatureModel{T<:AbstractFloat,R<:CurvatureReadoutModel,B<:CurvatureBranchResponse{T}} <: SubaruHighOrderWFSModel
    defocus_rms_nm::T
    readout_model::R
    branch_response::B
    crop_samples_per_subap::Int
    readout_pixels_per_subap::Int
end
AO188CurvatureModel(; defocus_rms_nm::Real=500.0,
    readout_model::CurvatureReadoutModel=CurvatureCountingReadout(),
    branch_response::CurvatureBranchResponse=CurvatureBranchResponse(),
    crop_samples_per_subap::Integer=8,
    readout_pixels_per_subap::Integer=1,
    T::Type{<:AbstractFloat}=Float32) =
    AO188CurvatureModel{T, typeof(readout_model), CurvatureBranchResponse{T}}(
        T(defocus_rms_nm), readout_model, CurvatureBranchResponse(T=T,
            plus_throughput=branch_response.plus_throughput,
            minus_throughput=branch_response.minus_throughput,
            plus_background=branch_response.plus_background,
            minus_background=branch_response.minus_background),
        Int(crop_samples_per_subap),
        Int(readout_pixels_per_subap))

function ao188_curvature_readout_crop_resolution(resolution::Int, n_subap::Int, crop_samples_per_subap::Int)
    crop_samples_per_subap >= 1 ||
        throw(InvalidConfiguration("AO188 curvature crop_samples_per_subap must be >= 1"))
    max_divisible = resolution - mod(resolution, n_subap)
    max_divisible > 0 || throw(InvalidConfiguration("AO188 curvature resolution must be divisible by n_subap"))
    requested = n_subap * crop_samples_per_subap
    return min(requested, max_divisible)
end

function AO188CurvatureSimulationParams(; kwargs...)
    nt = (; kwargs...)
    T0 = get(nt, :T, Float32)
    sampling = get(nt, :sampling_time, 1e-3)
    high_detector = get(nt, :high_detector, AO188APDDetectorConfig(T=T0, integration_time=sampling))
    rest = Base.structdiff(nt, (; high_order_sensor_model=nothing, source_band=nothing, high_detector=nothing))
    return AO188SimulationParams(; source_band=:I, high_order_sensor_model=AO188CurvatureModel(T=T0),
        high_detector=high_detector, rest...)
end

function _build_high_order_wfs(::OperationalShackHartmannModel, tel::Telescope, params; backend=CPUBackend())
    T = eltype(tel.state.opd)
    return ShackHartmann(tel; n_subap=params.n_subap, mode=Diffractive(), T=T, backend=backend)
end

function _build_high_order_wfs(model::AO188CurvatureModel, tel::Telescope, params; backend=CPUBackend())
    T = eltype(tel.state.opd)
    readout_crop_resolution = ao188_curvature_readout_crop_resolution(
        tel.params.resolution, params.n_subap, model.crop_samples_per_subap)
    return CurvatureWFS(tel; n_subap=params.n_subap, defocus_rms_nm=model.defocus_rms_nm,
        readout_model=model.readout_model, branch_response=model.branch_response,
        readout_crop_resolution=readout_crop_resolution,
        readout_pixels_per_subap=model.readout_pixels_per_subap, T=T, backend=backend)
end

abstract type AO188ReplayMode end
struct DirectReplayMode <: AO188ReplayMode end
struct PreparedReplayMode <: AO188ReplayMode end

struct AO188LatencyModel
    high_measurement_delay_frames::Int
    low_measurement_delay_frames::Int
    reconstruction_delay_frames::Int
    dm_delay_frames::Int
end

function AO188LatencyModel(;
    high_measurement_delay_frames::Integer=1,
    low_measurement_delay_frames::Integer=1,
    reconstruction_delay_frames::Integer=0,
    dm_delay_frames::Integer=1,
)
    delays = (
        Int(high_measurement_delay_frames),
        Int(low_measurement_delay_frames),
        Int(reconstruction_delay_frames),
        Int(dm_delay_frames),
    )
    all(>=(0), delays) || throw(InvalidConfiguration("AO188 latency delays must be >= 0"))
    return AO188LatencyModel(delays...)
end

abstract type AO188DetectorConfig end

struct AO188WFSDetectorConfig{T<:AbstractFloat,N<:NoiseModel,S<:SensorType,R<:Union{Nothing,FrameResponseModel},C<:FrameReadoutCorrectionModel,TM<:Union{Nothing,AbstractDetectorThermalModel}} <: AO188DetectorConfig
    integration_time::T
    qe::T
    psf_sampling::Int
    binning::Int
    gain::T
    dark_current::T
    noise::N
    sensor::S
    response_model::R
    correction_model::C
    thermal_model::TM
end

function AO188WFSDetectorConfig(;
    T::Type{<:AbstractFloat}=Float32,
    integration_time::Real=1e-3,
    qe::Real=0.9,
    psf_sampling::Int=1,
    binning::Int=1,
    gain::Real=1.0,
    dark_current::Real=0.0,
    noise::NoiseModel=NoisePhotonReadout(0.3),
    sensor::SensorType=CCDSensor(),
    response_model::Union{Nothing,FrameResponseModel}=nothing,
    correction_model::FrameReadoutCorrectionModel=NullFrameReadoutCorrection(),
    thermal_model::Union{Nothing,AbstractDetectorThermalModel}=nothing,
)
    return AO188WFSDetectorConfig{T,typeof(convert_noise(noise, T)),typeof(sensor),typeof(response_model),typeof(correction_model),typeof(thermal_model)}(
        T(integration_time),
        T(qe),
        psf_sampling,
        binning,
        T(gain),
        T(dark_current),
        convert_noise(validate_noise(noise), T),
        sensor,
        response_model,
        correction_model,
        thermal_model,
    )
end

struct AO188APDDetectorConfig{T<:AbstractFloat,N<:NoiseModel,D<:CountingDeadTimeModel,GM,TM<:Union{Nothing,AbstractDetectorThermalModel}} <: AO188DetectorConfig
    integration_time::T
    qe::T
    gain::T
    dark_count_rate::T
    noise::N
    dead_time_model::D
    output_precision::Union{Nothing,DataType}
    channel_gain_map::GM
    thermal_model::TM
end

function AO188APDDetectorConfig(;
    T::Type{<:AbstractFloat}=Float32,
    integration_time::Real=1e-3,
    qe::Real=0.9,
    gain::Real=1.0,
    dark_count_rate::Real=0.0,
    noise::NoiseModel=NoisePhoton(),
    dead_time_model::CountingDeadTimeModel=NoDeadTime(),
    output_precision::Union{Nothing,DataType}=nothing,
    channel_gain_map=nothing,
    thermal_model::Union{Nothing,AbstractDetectorThermalModel}=nothing,
)
    return AO188APDDetectorConfig{T,typeof(convert_noise(noise, T)),typeof(dead_time_model),typeof(channel_gain_map),typeof(thermal_model)}(
        T(integration_time),
        T(qe),
        T(gain),
        T(dark_count_rate),
        convert_noise(noise, T),
        dead_time_model,
        output_precision,
        channel_gain_map,
        thermal_model,
    )
end

function detector_from_config(cfg::AO188WFSDetectorConfig{T}; backend=CPUBackend()) where {T<:AbstractFloat}
    return Detector(
        cfg.noise;
        integration_time=cfg.integration_time,
        qe=cfg.qe,
        psf_sampling=cfg.psf_sampling,
        binning=cfg.binning,
        gain=cfg.gain,
        dark_current=cfg.dark_current,
        sensor=cfg.sensor,
        response_model=cfg.response_model,
        correction_model=cfg.correction_model,
        thermal_model=cfg.thermal_model,
        T=T,
        backend=backend,
    )
end

function detector_from_config(cfg::AO188APDDetectorConfig{T}; backend=CPUBackend()) where {T<:AbstractFloat}
    return APDDetector(
        integration_time=cfg.integration_time,
        qe=cfg.qe,
        gain=cfg.gain,
        dark_count_rate=cfg.dark_count_rate,
        noise=cfg.noise,
        dead_time_model=cfg.dead_time_model,
        output_precision=cfg.output_precision,
        channel_gain_map=cfg.channel_gain_map,
        thermal_model=something(cfg.thermal_model, NullDetectorThermalModel()),
        T=T,
        backend=backend,
    )
end

function default_ao188_low_order_resolution(resolution::Integer, n_low_order_subap::Integer)
    resolution > 0 || throw(InvalidConfiguration("resolution must be > 0"))
    n_low_order_subap > 0 || throw(InvalidConfiguration("n_low_order_subap must be > 0"))
    target = max(n_low_order_subap, resolution ÷ 4)
    snapped = target - mod(target, n_low_order_subap)
    return max(n_low_order_subap, snapped)
end

struct AO188SimulationParams{
    T<:AbstractFloat,
    F<:FidelityProfile,
    S<:AO188ActuatorSupportModel,
    W<:SubaruHighOrderWFSModel,
    B<:AbstractExecutionPolicy,
    R<:AO188ReplayMode,
    H,
    L,
}
    diameter::T
    sampling_time::T
    central_obstruction::T
    resolution::Int
    n_act::Int
    n_active_actuators::Int
    n_control_modes::Int
    control_grid_side::Int
    n_subap::Int
    n_low_order_subap::Int
    low_order_resolution::Int
    n_low_order_modes::Int
    influence_width::T
    r0::T
    L0::T
    source_magnitude::T
    interaction_amplitude::T
    control_gain::T
    low_order_gain::T
    control_sign::T
    profile::F
    source_band::Symbol
    support_model::S
    high_order_sensor_model::W
    branch_execution::B
    replay_mode::R
    latency::AO188LatencyModel
    high_detector::H
    low_detector::L
end

function AO188SimulationParams(;
    T::Type{<:AbstractFloat}=Float32,
    diameter::Real=8.2,
    sampling_time::Real=1e-3,
    central_obstruction::Real=0.30,
    resolution::Int=112,
    n_act::Int=64,
    n_active_actuators::Int=3228,
    n_control_modes::Int=188,
    control_grid_side::Int=16,
    n_subap::Int=14,
    n_low_order_subap::Int=2,
    low_order_resolution::Union{Int,Nothing}=nothing,
    n_low_order_modes::Int=4,
    influence_width::Real=0.3,
    r0::Real=0.16,
    L0::Real=25.0,
    source_magnitude::Real=8.0,
    interaction_amplitude::Real=0.05,
    control_gain::Real=0.5,
    low_order_gain::Real=0.15,
    control_sign::Real=-1.0,
    profile::FidelityProfile=FastProfile(),
    source_band::Symbol=:I,
    support_model::AO188ActuatorSupportModel=CircularActuatorSupport(),
    high_order_sensor_model::SubaruHighOrderWFSModel=OperationalShackHartmannModel(),
    branch_execution::AbstractExecutionPolicy=SequentialExecution(),
    replay_mode::AO188ReplayMode=DirectReplayMode(),
    latency::AO188LatencyModel=AO188LatencyModel(),
    high_detector::Union{Nothing,AO188DetectorConfig}=AO188WFSDetectorConfig(
        T=T,
        integration_time=sampling_time,
        qe=0.9,
        psf_sampling=1,
        binning=1,
        gain=1.0,
        dark_current=0.0,
        noise=NoisePhotonReadout(0.3),
        sensor=CCDSensor(),
    ),
    low_detector::AO188WFSDetectorConfig=AO188WFSDetectorConfig(
        T=T,
        integration_time=sampling_time,
        qe=0.95,
        psf_sampling=1,
        binning=1,
        gain=1.0,
        dark_current=0.0,
        noise=NoisePhotonReadout(0.1),
        sensor=CCDSensor(),
    ),
)
    resolved_low_order_resolution = isnothing(low_order_resolution) ?
        default_ao188_low_order_resolution(resolution, n_low_order_subap) :
        low_order_resolution

    return AO188SimulationParams{T,typeof(profile),typeof(support_model),typeof(high_order_sensor_model),typeof(branch_execution),typeof(replay_mode),typeof(high_detector),typeof(low_detector)}(
        T(diameter),
        T(sampling_time),
        T(central_obstruction),
        resolution,
        n_act,
        n_active_actuators,
        n_control_modes,
        control_grid_side,
        n_subap,
        n_low_order_subap,
        resolved_low_order_resolution,
        n_low_order_modes,
        T(influence_width),
        T(r0),
        T(L0),
        T(source_magnitude),
        T(interaction_amplitude),
        T(control_gain),
        T(low_order_gain),
        T(control_sign),
        profile,
        source_band,
        support_model,
        high_order_sensor_model,
        branch_execution,
        replay_mode,
        latency,
        high_detector,
        low_detector,
    )
end

struct AO1883kPhaseTimingStats{T<:AbstractFloat}
    samples::Int
    high_sense_mean_ns::T
    low_sense_mean_ns::T
    high_reconstruct_mean_ns::T
    low_reconstruct_mean_ns::T
    delay_mean_ns::T
    apply_mean_ns::T
    total_mean_ns::T
    total_p95_ns::T
end

mutable struct AO188Simulation{
    P,
    TEL,
    LOWTEL,
    ATM,
    SRC,
    DM,
    HOWFS,
    LOWFS,
    HDET,
    LDET,
    HREC,
    LREC,
    CMD,
    MASK,
    IDX,
    OVL,
    M2CH,
    M2CL,
    DLH,
    DLL,
    DLR,
    DLD,
    BESTATE,
    RNG,
} <: AbstractControlSimulation
    params::P
    tel::TEL
    low_tel::LOWTEL
    atm::ATM
    src::SRC
    dm::DM
    high_wfs::HOWFS
    low_wfs::LOWFS
    high_detector::HDET
    low_detector::LDET
    high_reconstructor::HREC
    low_reconstructor::LREC
    high_command::CMD
    low_command::CMD
    combined_command::CMD
    command::CMD
    active_mask::MASK
    active_indices::IDX
    actuator_overlap::OVL
    high_M2C::M2CH
    low_M2C::M2CL
    high_measurement_delay::DLH
    low_measurement_delay::DLL
    reconstruction_delay::DLR
    dm_delay::DLD
    branch_execution_state::BESTATE
    rng::RNG
    replay_prepared::Bool
end

function _dm_grid_coordinates(n_act::Int, misregistration::Misregistration)
    xs = range(-1.0, 1.0; length=n_act)
    ys = range(-1.0, 1.0; length=n_act)
    x_coords = Vector{Float64}(undef, n_act * n_act)
    y_coords = Vector{Float64}(undef, n_act * n_act)
    idx = 1
    for x0 in xs, y0 in ys
        x_m, y_m = apply_misregistration(misregistration, x0, y0)
        x_coords[idx] = x_m
        y_coords[idx] = y_m
        idx += 1
    end
    return x_coords, y_coords
end

function _actuator_overlap_weights(tel::Telescope, dm::DeformableMirror)
    n = tel.params.resolution
    pupil = Array(tel.state.pupil)
    x_coords, y_coords = _dm_grid_coordinates(dm.params.n_act, dm.params.misregistration)
    cx = (n + 1) / 2
    scale = n / 2
    px = [((i - cx) / scale) for i in 1:n]
    py = [((j - cx) / scale) for j in 1:n]
    sigma2 = Float64(dm.params.influence_width)^2
    overlaps = Vector{Float64}(undef, length(x_coords))
    @inbounds for k in eachindex(x_coords)
        x0 = x_coords[k]
        y0 = y_coords[k]
        acc = 0.0
        for j in 1:n, i in 1:n
            if pupil[i, j]
                dx = px[i] - x0
                dy = py[j] - y0
                acc += exp(-(dx * dx + dy * dy) / (2 * sigma2))
            end
        end
        overlaps[k] = acc
    end
    return overlaps, x_coords, y_coords
end

function _active_actuator_support(::CircularActuatorSupport, overlaps::AbstractVector{<:Real},
    x_coords::AbstractVector{<:Real}, y_coords::AbstractVector{<:Real}, n_act::Int, n_active::Int)
    total = n_act * n_act
    0 < n_active <= total || throw(InvalidConfiguration("n_active_actuators must be in 1:$total"))
    radius2 = similar(overlaps, Float64)
    @inbounds for i in eachindex(radius2)
        radius2[i] = x_coords[i]^2 + y_coords[i]^2
    end
    order = sortperm(eachindex(overlaps); by=i -> (radius2[i], -Float64(overlaps[i])))
    active_indices = sort(order[1:n_active])
    active_mask = falses(n_act, n_act)
    active_mask[active_indices] .= true
    return active_mask, active_indices
end

function _grouped_command_basis(active_mask::AbstractMatrix{Bool}, overlaps::AbstractVector{<:Real},
    n_modes::Int, grid_side::Int, T::Type{<:AbstractFloat})
    n_act = size(active_mask, 1)
    size(active_mask, 1) == size(active_mask, 2) ||
        throw(InvalidConfiguration("AO1883k active support must be square"))
    grid_side > 0 || throw(InvalidConfiguration("control_grid_side must be positive"))
    groups = [Int[] for _ in 1:grid_side * grid_side]
    scores = zeros(Float64, grid_side * grid_side)
    @inbounds for linear_idx in eachindex(active_mask)
        active_mask[linear_idx] || continue
        row, col = Tuple(CartesianIndices(active_mask)[linear_idx])
        grid_row = min(grid_side, fld((row - 1) * grid_side, n_act) + 1)
        grid_col = min(grid_side, fld((col - 1) * grid_side, n_act) + 1)
        group_idx = grid_row + (grid_col - 1) * grid_side
        push!(groups[group_idx], linear_idx)
        scores[group_idx] += Float64(overlaps[linear_idx])
    end
    occupied = findall(!isempty, groups)
    length(occupied) >= n_modes || throw(InvalidConfiguration(
        "control_grid_side=$grid_side yields only $(length(occupied)) occupied groups; need $n_modes control modes"))
    selected = occupied[partialsortperm(scores[occupied], 1:n_modes; rev=true)]
    M2C = zeros(T, n_act * n_act, n_modes)
    @inbounds for (mode_idx, group_idx) in enumerate(selected)
        members = groups[group_idx]
        weight = inv(sqrt(T(length(members))))
        for actuator_idx in members
            M2C[actuator_idx, mode_idx] = weight
        end
    end
    return M2C
end

function _low_order_command_basis(dm::DeformableMirror, tel::Telescope, active_mask::AbstractMatrix{Bool},
    n_modes::Int, T::Type{<:AbstractFloat})
    M2C_native, _ = kl_modal_basis(KLDMModes(), dm, tel; n_modes=n_modes, remove_piston=true)
    M2C = Matrix{T}(Array(M2C_native[:, 1:n_modes]))
    inactive = .!vec(active_mask)
    @views M2C[inactive, :] .= zero(T)
    @inbounds for j in axes(M2C, 2)
        col = @view M2C[:, j]
        norm_col = sqrt(sum(abs2, col))
        if !iszero(norm_col)
            col ./= norm_col
        end
    end
    return M2C
end

function _full_command_reconstructor(M2C_host::AbstractMatrix{T}, imat::InteractionMatrix{T};
    gain::Real, policy::InversePolicy, inverse_build_backend::AdaptiveOpticsSim.BuildBackend,
    materialize_backend::AdaptiveOpticsSim.BuildBackend, ref::AbstractMatrix{T}) where {T<:AbstractFloat}
    return MappedReconstructor(
        M2C_host,
        imat;
        gain=gain,
        policy=policy,
        inverse_build_backend=inverse_build_backend,
        materialize_backend=materialize_backend,
        ref=ref,
    )
end

function _auto_build_backend(backend)
    backend === Array && return AdaptiveOpticsSim.NativeBuildBackend()
    for B in AdaptiveOpticsSim.available_gpu_backends()
        if backend === AdaptiveOpticsSim.gpu_backend_array_type(B)
            return AdaptiveOpticsSim.GPUArrayBuildBackend(B)
        end
    end
    return AdaptiveOpticsSim.NativeBuildBackend()
end

function _ao188_calibration_objects(params::AO188SimulationParams{T}) where {T<:AbstractFloat}
    tel = Telescope(
        resolution=params.resolution,
        diameter=params.diameter,
        sampling_time=params.sampling_time,
        central_obstruction=params.central_obstruction,
        T=T,
        backend=CPUBackend(),
    )
    low_tel = Telescope(
        resolution=params.low_order_resolution,
        diameter=params.diameter,
        sampling_time=params.sampling_time,
        central_obstruction=params.central_obstruction,
        T=T,
        backend=CPUBackend(),
    )
    src = Source(band=params.source_band, magnitude=params.source_magnitude, T=T)
    dm = DeformableMirror(tel; n_act=params.n_act, influence_width=params.influence_width, T=T, backend=CPUBackend())
    low_dm = DeformableMirror(low_tel; n_act=params.n_act, influence_width=params.influence_width, T=T, backend=CPUBackend())
    high_wfs = _build_high_order_wfs(params.high_order_sensor_model, tel, params; backend=CPUBackend())
    low_wfs = ShackHartmann(low_tel; n_subap=params.n_low_order_subap, mode=Diffractive(), T=T, backend=CPUBackend())
    return tel, low_tel, src, dm, low_dm, high_wfs, low_wfs
end

_high_order_wfs_frame(simulation::AO188Simulation) =
    isnothing(simulation.high_detector) ? simulation.high_wfs.state.camera_frame :
    wfs_output_frame(simulation.high_wfs, simulation.high_detector)
_high_order_wfs_metadata(simulation::AO188Simulation) =
    isnothing(simulation.high_detector) ? wfs_output_metadata(simulation.high_wfs) :
    detector_export_metadata(simulation.high_detector)

function _measure_high!(surrogate::AO188Simulation, rng::AbstractRNG)
    if isnothing(surrogate.high_detector)
        return measure!(surrogate.high_wfs, surrogate.tel, surrogate.src)
    end
    return measure!(surrogate.high_wfs, surrogate.tel, surrogate.src, surrogate.high_detector; rng=rng)
end

function _downsample_low_order_opd!(surrogate::AO188Simulation)
    high_res = surrogate.tel.params.resolution
    low_res = surrogate.low_tel.params.resolution
    factor = div(high_res, low_res)
    bin2d!(surrogate.low_tel.state.opd, surrogate.tel.state.opd, factor)
    surrogate.low_tel.state.opd .*= inv(eltype(surrogate.low_tel.state.opd)(factor * factor))
    surrogate.low_tel.state.opd .*= surrogate.low_tel.state.pupil
    return surrogate.low_tel
end

function _measure_low!(surrogate::AO188Simulation, rng::AbstractRNG)
    _downsample_low_order_opd!(surrogate)
    return measure!(surrogate.low_wfs, surrogate.low_tel, surrogate.src, surrogate.low_detector; rng=rng)
end

function _frame_rngs!(rng::AbstractRNG)
    return MersenneTwister(rand(rng, UInt64)), MersenneTwister(rand(rng, UInt64))
end

function measure_branches_backend!(::BackendStreamExecution, simulation::AO188Simulation, state)
    return _measure_branches!(SequentialExecution(), simulation)
end

function _measure_branches!(::SequentialExecution, simulation::AO188Simulation)
    _measure_high!(simulation, simulation.rng)
    _measure_low!(simulation, simulation.rng)
    return simulation
end

function _measure_branches!(::ThreadedExecution, simulation::AO188Simulation)
    Threads.nthreads() < 2 && return _measure_branches!(SequentialExecution(), simulation)
    high_rng, low_rng = _frame_rngs!(simulation.rng)
    high_task = Threads.@spawn _measure_high!(simulation, high_rng)
    low_task = Threads.@spawn _measure_low!(simulation, low_rng)
    fetch(high_task)
    fetch(low_task)
    return simulation
end

function _measure_branches!(::BackendStreamExecution, simulation::AO188Simulation)
    return measure_branches_backend!(BackendStreamExecution(), simulation, simulation.branch_execution_state)
end

function prepare_replay!(simulation::AO188Simulation)
    prepare_runtime_wfs!(simulation.high_wfs, simulation.tel, simulation.src)
    prepare_sampling!(simulation.low_wfs, simulation.low_tel, simulation.src)
    ensure_sh_calibration!(simulation.low_wfs, simulation.low_tel, simulation.src)
    simulation.replay_prepared = true
    return simulation
end

function subaru_ao188_simulation(; params::AO188SimulationParams=AO188SimulationParams(),
    backend=CPUBackend(), build_backend::Union{Nothing,AdaptiveOpticsSim.BuildBackend}=nothing, rng=MersenneTwister(0))
    resolved_materialize_backend = isnothing(build_backend) ? _auto_build_backend(backend) : build_backend
    resolved_calibration_backend = isnothing(build_backend) && backend !== Array ? AdaptiveOpticsSim.CPUBuildBackend() : resolved_materialize_backend
    T = typeof(params.diameter)
    if params.resolution % params.low_order_resolution != 0
        throw(InvalidConfiguration("low_order_resolution must evenly divide the main telescope resolution"))
    end
    if params.low_order_resolution % params.n_low_order_subap != 0
        throw(InvalidConfiguration("low_order_resolution must be divisible by n_low_order_subap"))
    end
    tel = Telescope(
        resolution=params.resolution,
        diameter=params.diameter,
        sampling_time=params.sampling_time,
        central_obstruction=params.central_obstruction,
        T=T,
        backend=backend,
    )
    low_tel = Telescope(
        resolution=params.low_order_resolution,
        diameter=params.diameter,
        sampling_time=params.sampling_time,
        central_obstruction=params.central_obstruction,
        T=T,
        backend=backend,
    )
    src = Source(band=params.source_band, magnitude=params.source_magnitude, T=T)
    atm = KolmogorovAtmosphere(tel; r0=params.r0, L0=params.L0, T=T, backend=backend)
    dm = DeformableMirror(tel; n_act=params.n_act, influence_width=params.influence_width, T=T, backend=backend)
    low_dm = DeformableMirror(low_tel; n_act=params.n_act, influence_width=params.influence_width, T=T, backend=backend)
    high_wfs = _build_high_order_wfs(params.high_order_sensor_model, tel, params; backend=backend)
    low_wfs = ShackHartmann(low_tel; n_subap=params.n_low_order_subap, mode=Diffractive(), T=T, backend=backend)
    high_detector = isnothing(params.high_detector) ? nothing : detector_from_config(params.high_detector; backend=backend)
    low_detector = detector_from_config(params.low_detector; backend=backend)

    calibration_tel = tel
    calibration_low_tel = low_tel
    calibration_src = src
    calibration_dm = dm
    calibration_low_dm = low_dm
    calibration_high_wfs = high_wfs
    calibration_low_wfs = low_wfs
    if resolved_calibration_backend isa AdaptiveOpticsSim.CPUBuildBackend && backend !== Array
        calibration_tel, calibration_low_tel, calibration_src,
        calibration_dm, calibration_low_dm, calibration_high_wfs, calibration_low_wfs =
            _ao188_calibration_objects(params)
    end

    overlap_weights, x_coords, y_coords = _actuator_overlap_weights(calibration_tel, calibration_dm)
    overlap_t = T.(overlap_weights)
    active_mask, active_indices = _active_actuator_support(
        params.support_model, overlap_weights, x_coords, y_coords, params.n_act, params.n_active_actuators)
    high_M2C_host = _grouped_command_basis(active_mask, overlap_weights, params.n_control_modes, params.control_grid_side, T)
    low_M2C_host = _low_order_command_basis(calibration_dm, calibration_tel, active_mask, params.n_low_order_modes, T)

    high_M2C = materialize_build(resolved_materialize_backend, dm.state.modes, high_M2C_host)
    low_M2C = materialize_build(resolved_materialize_backend, dm.state.modes, low_M2C_host)
    policy = default_modal_inverse_policy(T)

    high_imat = interaction_matrix(calibration_dm, calibration_high_wfs, calibration_tel, high_M2C_host,
        calibration_src; amplitude=params.interaction_amplitude)
    low_imat = interaction_matrix(calibration_low_dm, calibration_low_wfs, calibration_low_tel, low_M2C_host,
        calibration_src; amplitude=params.interaction_amplitude)
    high_recon = _full_command_reconstructor(high_M2C_host, high_imat;
        gain=params.control_gain, policy=policy, inverse_build_backend=resolved_calibration_backend,
        materialize_backend=resolved_materialize_backend, ref=dm.state.modes)
    low_recon = _full_command_reconstructor(low_M2C_host, low_imat;
        gain=params.low_order_gain, policy=policy, inverse_build_backend=resolved_calibration_backend,
        materialize_backend=resolved_materialize_backend, ref=dm.state.modes)

    high_command = similar(dm.state.coefs)
    low_command = similar(dm.state.coefs)
    combined_command = similar(dm.state.coefs)
    command = similar(dm.state.coefs)
    fill!(high_command, zero(T))
    fill!(low_command, zero(T))
    fill!(combined_command, zero(T))
    fill!(command, zero(T))
    branch_execution_state = init_execution_state(params.branch_execution, dm.state.coefs)

    surrogate = AO188Simulation(
        params,
        tel,
        low_tel,
        atm,
        src,
        dm,
        high_wfs,
        low_wfs,
        high_detector,
        low_detector,
        high_recon,
        low_recon,
        high_command,
        low_command,
        combined_command,
        command,
        active_mask,
        active_indices,
        overlap_t,
        high_M2C,
        low_M2C,
        VectorDelayLine(high_wfs.state.slopes, params.latency.high_measurement_delay_frames),
        VectorDelayLine(low_wfs.state.slopes, params.latency.low_measurement_delay_frames),
        VectorDelayLine(command, params.latency.reconstruction_delay_frames),
        VectorDelayLine(command, params.latency.dm_delay_frames),
        branch_execution_state,
        rng,
        false,
    )
    if params.replay_mode isa PreparedReplayMode
        prepare!(surrogate)
    end
    return surrogate
end

function subaru_ao188_curvature_simulation(; params::AO188SimulationParams=AO188CurvatureSimulationParams(),
    kwargs...)
    return subaru_ao188_simulation(; params=params, kwargs...)
end

function step!(surrogate::AO188Simulation)
    if surrogate.params.replay_mode isa PreparedReplayMode && !surrogate.replay_prepared
        prepare!(surrogate)
    end
    advance!(surrogate.atm, surrogate.tel, surrogate.rng)
    propagate!(surrogate.atm, surrogate.tel)
    apply!(surrogate.dm, surrogate.tel, DMAdditive())

    _measure_branches!(surrogate.params.branch_execution, surrogate)

    high_slopes = shift_delay!(surrogate.high_measurement_delay, surrogate.high_wfs.state.slopes)
    low_slopes = shift_delay!(surrogate.low_measurement_delay, surrogate.low_wfs.state.slopes)

    reconstruct!(surrogate.high_command, surrogate.high_reconstructor, high_slopes)
    reconstruct!(surrogate.low_command, surrogate.low_reconstructor, low_slopes)
    surrogate.combined_command .= surrogate.high_command .+ surrogate.low_command
    delayed_recon = shift_delay!(surrogate.reconstruction_delay, surrogate.combined_command)
    copyto!(surrogate.command, delayed_recon)
    delayed_dm = shift_delay!(surrogate.dm_delay, surrogate.command)
    apply_command!(surrogate.dm.state.coefs, delayed_dm, surrogate.params.control_sign)
    return surrogate.command
end

runtime_timing(surrogate::AO188Simulation; kwargs...) = runtime_timing(() -> step!(surrogate); kwargs...)

prepare!(simulation::AO188Simulation) = prepare_replay!(simulation)
supports_prepared_runtime(::Type{<:AO188Simulation}) = true
supports_detector_output(::Type{<:AO188Simulation}) = true
supports_grouped_execution(::Type{<:AO188Simulation}) = true

function simulation_interface(simulation::AO188Simulation)
    return SimulationReadout(
        simulation.command,
        (simulation.high_wfs.state.slopes, simulation.low_wfs.state.slopes),
        (
            _high_order_wfs_frame(simulation),
            wfs_output_frame(simulation.low_wfs, simulation.low_detector),
        ),
        nothing,
        (
            _high_order_wfs_metadata(simulation),
            detector_export_metadata(simulation.low_detector),
        ),
        nothing,
        nothing,
        nothing,
    )
end

function subaru_ao188_phase_timing(surrogate::AO188Simulation; warmup::Int=10, samples::Int=1000, gc_before::Bool=true)
    warmup >= 0 || throw(InvalidConfiguration("warmup must be >= 0"))
    samples > 0 || throw(InvalidConfiguration("samples must be positive"))
    for _ in 1:warmup
        step!(surrogate)
    end
    gc_before && GC.gc()
    high_sense_times = Vector{Int}(undef, samples)
    low_sense_times = Vector{Int}(undef, samples)
    high_reconstruct_times = Vector{Int}(undef, samples)
    low_reconstruct_times = Vector{Int}(undef, samples)
    delay_times = Vector{Int}(undef, samples)
    apply_times = Vector{Int}(undef, samples)
    total_times = Vector{Int}(undef, samples)
    @inbounds for i in 1:samples
        t0 = time_ns()
        advance!(surrogate.atm, surrogate.tel, surrogate.rng)
        propagate!(surrogate.atm, surrogate.tel)
        apply!(surrogate.dm, surrogate.tel, DMAdditive())
        if isnothing(surrogate.high_detector)
            measure!(surrogate.high_wfs, surrogate.tel, surrogate.src)
        else
            measure!(surrogate.high_wfs, surrogate.tel, surrogate.src, surrogate.high_detector; rng=surrogate.rng)
        end
        synchronize_backend!(execution_style(surrogate.high_wfs.state.slopes))
        t1 = time_ns()
        _downsample_low_order_opd!(surrogate)
        measure!(surrogate.low_wfs, surrogate.low_tel, surrogate.src, surrogate.low_detector; rng=surrogate.rng)
        synchronize_backend!(execution_style(surrogate.low_wfs.state.slopes))
        t2 = time_ns()

        high_slopes = shift_delay!(surrogate.high_measurement_delay, surrogate.high_wfs.state.slopes)
        reconstruct!(surrogate.high_command, surrogate.high_reconstructor, high_slopes)
        synchronize_backend!(execution_style(surrogate.high_command))
        t3 = time_ns()
        low_slopes = shift_delay!(surrogate.low_measurement_delay, surrogate.low_wfs.state.slopes)
        reconstruct!(surrogate.low_command, surrogate.low_reconstructor, low_slopes)
        synchronize_backend!(execution_style(surrogate.low_command))
        surrogate.combined_command .= surrogate.high_command .+ surrogate.low_command
        synchronize_backend!(execution_style(surrogate.combined_command))
        t4 = time_ns()

        delayed_recon = shift_delay!(surrogate.reconstruction_delay, surrogate.combined_command)
        copyto!(surrogate.command, delayed_recon)
        delayed_dm = shift_delay!(surrogate.dm_delay, surrogate.command)
        t5 = time_ns()

        apply_command!(surrogate.dm.state.coefs, delayed_dm, surrogate.params.control_sign)
        synchronize_backend!(execution_style(surrogate.dm.state.coefs))
        t6 = time_ns()
        high_sense_times[i] = t1 - t0
        low_sense_times[i] = t2 - t1
        high_reconstruct_times[i] = t3 - t2
        low_reconstruct_times[i] = t4 - t3
        delay_times[i] = t5 - t4
        apply_times[i] = t6 - t5
        total_times[i] = t6 - t0
    end
    sorted_total = sort(copy(total_times))
    p95_idx = clamp(round(Int, 0.95 * samples), 1, samples)
    return AO1883kPhaseTimingStats{Float64}(
        samples,
        mean(high_sense_times),
        mean(low_sense_times),
        mean(high_reconstruct_times),
        mean(low_reconstruct_times),
        mean(delay_times),
        mean(apply_times),
        mean(total_times),
        sorted_total[p95_idx],
    )
end

if isdefined(Main, :CUDA)
    const _AO188CUDA = Main.CUDA

    struct AO188CUDAStreamState
        high::_AO188CUDA.CuStream
        low::_AO188CUDA.CuStream
    end

    function init_execution_state(::BackendStreamExecution, ::_AO188CUDA.CuArray)
        return AO188CUDAStreamState(
            _AO188CUDA.CuStream(; flags=_AO188CUDA.STREAM_NON_BLOCKING),
            _AO188CUDA.CuStream(; flags=_AO188CUDA.STREAM_NON_BLOCKING),
        )
    end

    function measure_branches_backend!(::BackendStreamExecution, surrogate::AO188Simulation, state::AO188CUDAStreamState)
        high_rng, low_rng = _frame_rngs!(surrogate.rng)
        high_task = Threads.@spawn _AO188CUDA.stream!(state.high) do
            _measure_high!(surrogate, high_rng)
            _AO188CUDA.synchronize(state.high)
        end
        low_task = Threads.@spawn _AO188CUDA.stream!(state.low) do
            _measure_low!(surrogate, low_rng)
            _AO188CUDA.synchronize(state.low)
        end
        fetch(high_task)
        fetch(low_task)
        return surrogate
    end
end

end
