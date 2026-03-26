using Random

abstract type NoiseModel end
abstract type SensorType end
abstract type FrameSensorType <: SensorType end
abstract type CountingSensorType <: SensorType end
abstract type AbstractFrameDetector <: AbstractDetector end
abstract type AbstractCountingDetector <: AbstractDetector end
abstract type CountingDeadTimeModel end
abstract type FrameResponseModel end
abstract type BackgroundModel end
struct NoiseNone <: NoiseModel end
struct NoisePhoton <: NoiseModel end
struct NoiseReadout{T<:AbstractFloat} <: NoiseModel
    sigma::T
end
struct NoisePhotonReadout{T<:AbstractFloat} <: NoiseModel
    sigma::T
end
struct CCDSensor{T<:AbstractFloat} <: FrameSensorType
    clock_induced_charge_rate::T
end
struct CMOSSensor{T<:AbstractFloat} <: FrameSensorType
    column_readout_sigma::T
end
abstract type AvalancheFrameSensorType <: FrameSensorType end
struct EMCCDSensor{T<:AbstractFloat} <: AvalancheFrameSensorType
    excess_noise_factor::T
end
struct InGaAsSensor{T<:AbstractFloat} <: FrameSensorType
    glow_rate::T
end
struct SAPHIRASensor{T<:AbstractFloat} <: AvalancheFrameSensorType
    avalanche_gain::T
    excess_noise_factor::T
    glow_rate::T
end
struct APDSensor <: CountingSensorType end
struct NoDeadTime <: CountingDeadTimeModel end
struct NonParalyzableDeadTime{T<:AbstractFloat} <: CountingDeadTimeModel
    dead_time::T
end
struct NullFrameResponse <: FrameResponseModel end
struct SeparableGaussianPixelResponse{T<:AbstractFloat,V<:AbstractVector{T}} <: FrameResponseModel
    response_width_px::T
    kernel::V
end
struct NoBackground <: BackgroundModel end
struct ScalarBackground{T<:AbstractFloat} <: BackgroundModel
    level::T
end
struct BackgroundFrame{T<:AbstractFloat,A<:AbstractMatrix{T}} <: BackgroundModel
    map::A
end

struct DetectorExportMetadata{T<:AbstractFloat}
    integration_time::T
    qe::T
    psf_sampling::Int
    binning::Int
    gain::T
    dark_current::T
    bits::Union{Nothing,Int}
    full_well::Union{Nothing,T}
    sensor::Symbol
    noise::Symbol
    readout_sigma::Union{Nothing,T}
    output_precision::Union{Nothing,DataType}
    frame_size::Tuple{Int,Int}
    output_size::Tuple{Int,Int}
    frame_response::Symbol
    response_width_px::Union{Nothing,T}
end

struct CountingReadoutMetadata
    layout::Symbol
    output_size::Tuple{Int,Int}
    n_channels::Int
end

struct CountingDetectorExportMetadata{T<:AbstractFloat}
    integration_time::T
    qe::T
    gain::T
    dark_count_rate::T
    dead_time_model::Symbol
    dead_time::Union{Nothing,T}
    sensor::Symbol
    noise::Symbol
    output_precision::Union{Nothing,DataType}
    readout::CountingReadoutMetadata
end

NoiseReadout(sigma::Real) = NoiseReadout{Float64}(float(sigma))
NoisePhotonReadout(sigma::Real) = NoisePhotonReadout{Float64}(float(sigma))
NonParalyzableDeadTime(dead_time::Real) = NonParalyzableDeadTime{Float64}(float(dead_time))
function EMCCDSensor(; excess_noise_factor::Real=1.0, T::Type{<:AbstractFloat}=Float64)
    excess_noise_factor >= 1 || throw(InvalidConfiguration("EMCCDSensor excess_noise_factor must be >= 1"))
    return EMCCDSensor{T}(T(excess_noise_factor))
end
function CCDSensor(; clock_induced_charge_rate::Real=0.0, T::Type{<:AbstractFloat}=Float64)
    clock_induced_charge_rate >= 0 || throw(InvalidConfiguration("CCDSensor clock_induced_charge_rate must be >= 0"))
    return CCDSensor{T}(T(clock_induced_charge_rate))
end
function CMOSSensor(; column_readout_sigma::Real=0.0, T::Type{<:AbstractFloat}=Float64)
    column_readout_sigma >= 0 || throw(InvalidConfiguration("CMOSSensor column_readout_sigma must be >= 0"))
    return CMOSSensor{T}(T(column_readout_sigma))
end
function InGaAsSensor(; glow_rate::Real=0.0, T::Type{<:AbstractFloat}=Float64)
    glow_rate >= 0 || throw(InvalidConfiguration("InGaAsSensor glow_rate must be >= 0"))
    return InGaAsSensor{T}(T(glow_rate))
end
function SAPHIRASensor(; avalanche_gain::Real=1.0, excess_noise_factor::Real=1.0,
    glow_rate::Real=0.0, T::Type{<:AbstractFloat}=Float64)
    avalanche_gain >= 1 || throw(InvalidConfiguration("SAPHIRASensor avalanche_gain must be >= 1"))
    excess_noise_factor >= 1 || throw(InvalidConfiguration("SAPHIRASensor excess_noise_factor must be >= 1"))
    glow_rate >= 0 || throw(InvalidConfiguration("SAPHIRASensor glow_rate must be >= 0"))
    return SAPHIRASensor{T}(T(avalanche_gain), T(excess_noise_factor), T(glow_rate))
end

function SeparableGaussianPixelResponse(; response_width_px::Real=0.5, truncate_at::Real=3.0,
    T::Type{<:AbstractFloat}=Float64, backend=Array)
    response_width_px > 0 || throw(InvalidConfiguration("SeparableGaussianPixelResponse response_width_px must be > 0"))
    truncate_at > 0 || throw(InvalidConfiguration("SeparableGaussianPixelResponse truncate_at must be > 0"))
    radius = max(1, ceil(Int, truncate_at * response_width_px))
    host_kernel = Vector{T}(undef, 2 * radius + 1)
    inv_sigma2 = inv(T(response_width_px)^2)
    for (idx, offset) in enumerate(-radius:radius)
        host_kernel[idx] = exp(-T(0.5) * T(offset * offset) * inv_sigma2)
    end
    host_kernel ./= sum(host_kernel)
    kernel = backend{T}(undef, length(host_kernel))
    copyto!(kernel, host_kernel)
    return SeparableGaussianPixelResponse{T,typeof(kernel)}(T(response_width_px), kernel)
end

@kernel function separable_response_rows_kernel!(out, img, kernel, radius::Int, n::Int, m::Int, klen::Int)
    i, j = @index(Global, NTuple)
    if i <= n && j <= m
        acc = zero(eltype(out))
        @inbounds for kk in 1:klen
            jj = clamp(j + kk - radius - 1, 1, m)
            acc += kernel[kk] * img[i, jj]
        end
        @inbounds out[i, j] = acc
    end
end

@kernel function separable_response_cols_kernel!(out, img, kernel, radius::Int, n::Int, m::Int, klen::Int)
    i, j = @index(Global, NTuple)
    if i <= n && j <= m
        acc = zero(eltype(out))
        @inbounds for kk in 1:klen
            ii = clamp(i + kk - radius - 1, 1, n)
            acc += kernel[kk] * img[ii, j]
        end
        @inbounds out[i, j] = acc
    end
end

@kernel function add_column_noise_kernel!(frame, noise, sigma, n::Int, m::Int)
    i, j = @index(Global, NTuple)
    if i <= n && j <= m
        @inbounds frame[i, j] += sigma * noise[1, j]
    end
end

struct DetectorParams{T<:AbstractFloat,S<:SensorType,R<:FrameResponseModel}
    integration_time::T
    qe::T
    psf_sampling::Int
    binning::Int
    gain::T
    dark_current::T
    bits::Union{Nothing,Int}
    full_well::Union{Nothing,T}
    sensor::S
    response_model::R
    output_precision::Union{Nothing,DataType}
end

mutable struct DetectorState{T<:AbstractFloat,A<:AbstractMatrix{T},O}
    frame::A
    response_buffer::A
    bin_buffer::A
    noise_buffer::A
    accum_buffer::A
    output_buffer::O
    integrated_time::T
    readout_ready::Bool
end

struct Detector{N<:NoiseModel,P<:DetectorParams,S<:DetectorState,BF<:BackgroundModel,BM<:BackgroundModel} <: AbstractFrameDetector
    noise::N
    params::P
    state::S
    background_flux::BF
    background_map::BM
end

struct APDDetectorParams{T<:AbstractFloat,S<:APDSensor,D<:CountingDeadTimeModel}
    integration_time::T
    qe::T
    gain::T
    dark_count_rate::T
    dead_time_model::D
    sensor::S
    output_precision::Union{Nothing,DataType}
    layout::Symbol
end

mutable struct APDDetectorState{T<:AbstractFloat,A<:AbstractMatrix{T},O}
    channels::A
    noise_buffer::A
    output_buffer::O
end

struct APDDetector{N<:NoiseModel,P<:APDDetectorParams,S<:APDDetectorState,GM} <: AbstractCountingDetector
    noise::N
    params::P
    state::S
    channel_gain_map::GM
end

readout_ready(det::Detector) = det.state.readout_ready
output_frame(det::Detector) = det.state.output_buffer === nothing ? det.state.frame : det.state.output_buffer
readout_ready(det::APDDetector) = true
channel_output(det::APDDetector) = det.state.output_buffer === nothing ? det.state.channels : det.state.output_buffer
output_frame(det::APDDetector) = channel_output(det)
reset_integration!(det::APDDetector) = det

detector_noise_symbol(::NoiseNone) = :none
detector_noise_symbol(::NoisePhoton) = :photon
detector_noise_symbol(::NoiseReadout) = :readout
detector_noise_symbol(::NoisePhotonReadout) = :photon_readout

detector_sensor_symbol(::CCDSensor) = :ccd
detector_sensor_symbol(::CMOSSensor) = :cmos
detector_sensor_symbol(::EMCCDSensor) = :emccd
detector_sensor_symbol(::InGaAsSensor) = :ingaas
detector_sensor_symbol(::SAPHIRASensor) = :saphira
detector_sensor_symbol(::APDSensor) = :apd
counting_dead_time_symbol(::NoDeadTime) = :none
counting_dead_time_symbol(::NonParalyzableDeadTime) = :nonparalyzable
frame_response_symbol(::NullFrameResponse) = :none
frame_response_symbol(::SeparableGaussianPixelResponse) = :separable_gaussian

supports_detector_mtf(::AbstractFrameDetector) = false
supports_detector_mtf(det::Detector) = supports_detector_mtf(det.params.response_model)
supports_detector_mtf(::NullFrameResponse) = false
supports_detector_mtf(::SeparableGaussianPixelResponse) = true
supports_clock_induced_charge(::FrameSensorType) = false
supports_clock_induced_charge(::CCDSensor) = true
supports_column_readout_noise(::FrameSensorType) = false
supports_column_readout_noise(::CMOSSensor) = true
supports_avalanche_gain(::FrameSensorType) = false
supports_avalanche_gain(::AvalancheFrameSensorType) = true
supports_sensor_glow(::FrameSensorType) = false
supports_sensor_glow(::InGaAsSensor) = true
supports_sensor_glow(::SAPHIRASensor) = true

supports_counting_noise(::AbstractCountingDetector) = false
supports_counting_noise(::APDDetector) = true
supports_dead_time(::AbstractCountingDetector) = false
supports_dead_time(det::APDDetector) = supports_dead_time(det.params.dead_time_model)
supports_dead_time(::NoDeadTime) = false
supports_dead_time(::NonParalyzableDeadTime) = true
supports_channel_gain_map(::AbstractCountingDetector) = false
supports_channel_gain_map(det::APDDetector) = !isnothing(det.channel_gain_map)

sensor_is_frame_based(::SensorType) = false
sensor_is_frame_based(::FrameSensorType) = true
sensor_is_counting(::SensorType) = false
sensor_is_counting(::CountingSensorType) = true

detector_readout_sigma(::NoiseNone, ::Type{T}) where {T<:AbstractFloat} = nothing
detector_readout_sigma(::NoisePhoton, ::Type{T}) where {T<:AbstractFloat} = nothing
detector_readout_sigma(noise::NoiseReadout, ::Type{T}) where {T<:AbstractFloat} = T(noise.sigma)
detector_readout_sigma(noise::NoisePhotonReadout, ::Type{T}) where {T<:AbstractFloat} = T(noise.sigma)

function detector_export_metadata(det::Detector; T::Type{<:AbstractFloat}=eltype(det.state.frame))
    output = output_frame(det)
    full_well = det.params.full_well === nothing ? nothing : T(det.params.full_well)
    return DetectorExportMetadata{T}(
        T(det.params.integration_time),
        T(det.params.qe),
        det.params.psf_sampling,
        det.params.binning,
        T(det.params.gain),
        T(det.params.dark_current),
        det.params.bits,
        full_well,
        detector_sensor_symbol(det.params.sensor),
        detector_noise_symbol(det.noise),
        detector_readout_sigma(det.noise, T),
        det.params.output_precision,
        size(det.state.frame),
        size(output),
        frame_response_symbol(det.params.response_model),
        frame_response_width(det.params.response_model, T),
    )
end

function detector_export_metadata(det::APDDetector; T::Type{<:AbstractFloat}=eltype(det.state.channels))
    output = channel_output(det)
    return CountingDetectorExportMetadata{T}(
        T(det.params.integration_time),
        T(det.params.qe),
        T(det.params.gain),
        T(det.params.dark_count_rate),
        counting_dead_time_symbol(det.params.dead_time_model),
        counting_dead_time_value(det.params.dead_time_model, T),
        detector_sensor_symbol(det.params.sensor),
        detector_noise_symbol(det.noise),
        det.params.output_precision,
        CountingReadoutMetadata(det.params.layout, size(output), length(output)),
    )
end

function reset_integration!(det::Detector)
    fill!(det.state.accum_buffer, zero(eltype(det.state.accum_buffer)))
    det.state.integrated_time = zero(det.state.integrated_time)
    det.state.readout_ready = true
    return det
end

normalize_noise(noise::NoiseModel) = noise

function normalize_noise(noises::Tuple{Vararg{NoiseModel}})
    if isempty(noises)
        return NoiseNone()
    end
    has_photon = false
    has_readout = false
    sigma = 0.0
    sigma_set = false

    for noise in noises
        if noise isa NoiseNone
            continue
        elseif noise isa NoisePhoton
            has_photon = true
        elseif noise isa NoiseReadout
            has_readout = true
            sigma, sigma_set = merge_readout_sigma(sigma, sigma_set, noise.sigma)
        elseif noise isa NoisePhotonReadout
            has_photon = true
            has_readout = true
            sigma, sigma_set = merge_readout_sigma(sigma, sigma_set, noise.sigma)
        end
    end

    if has_photon && has_readout
        return NoisePhotonReadout(sigma)
    elseif has_photon
        return NoisePhoton()
    elseif has_readout
        return NoiseReadout(sigma)
    end
    return NoiseNone()
end

function merge_readout_sigma(current::Float64, set::Bool, new_sigma::Real)
    σ = float(new_sigma)
    if set && σ != current
        throw(InvalidConfiguration("conflicting readout noise values in noise tuple"))
    end
    return σ, true
end

convert_noise(noise::NoiseNone, ::Type{T}) where {T<:AbstractFloat} = NoiseNone()
convert_noise(noise::NoisePhoton, ::Type{T}) where {T<:AbstractFloat} = NoisePhoton()
convert_noise(noise::NoiseReadout, ::Type{T}) where {T<:AbstractFloat} = NoiseReadout{T}(T(noise.sigma))
convert_noise(noise::NoisePhotonReadout, ::Type{T}) where {T<:AbstractFloat} = NoisePhotonReadout{T}(T(noise.sigma))

validate_noise(noise::NoiseNone) = noise
validate_noise(noise::NoisePhoton) = noise
function validate_noise(noise::NoiseReadout)
    if noise.sigma <= 0
        throw(InvalidConfiguration("readout noise must be > 0 for NoiseReadout"))
    end
    return noise
end
function validate_noise(noise::NoisePhotonReadout)
    if noise.sigma <= 0
        throw(InvalidConfiguration("readout noise must be > 0 for NoisePhotonReadout"))
    end
    return noise
end

background_model(::Nothing; T::Type{<:AbstractFloat}, backend) = NoBackground()
background_model(level::Real; T::Type{<:AbstractFloat}, backend) = ScalarBackground{T}(T(level))
function background_model(map::AbstractMatrix; T::Type{<:AbstractFloat}, backend)
    background = backend{T}(undef, size(map)...)
    copyto!(background, T.(map))
    return BackgroundFrame{T, typeof(background)}(background)
end

counting_dead_time_value(::NoDeadTime, ::Type{T}) where {T<:AbstractFloat} = nothing
counting_dead_time_value(model::NonParalyzableDeadTime, ::Type{T}) where {T<:AbstractFloat} = T(model.dead_time)
frame_response_width(::NullFrameResponse, ::Type{T}) where {T<:AbstractFloat} = nothing
frame_response_width(model::SeparableGaussianPixelResponse, ::Type{T}) where {T<:AbstractFloat} = T(model.response_width_px)

convert_dead_time_model(::NoDeadTime, ::Type{T}) where {T<:AbstractFloat} = NoDeadTime()
convert_dead_time_model(model::NonParalyzableDeadTime, ::Type{T}) where {T<:AbstractFloat} =
    NonParalyzableDeadTime{T}(T(model.dead_time))
convert_frame_response_model(::NullFrameResponse, ::Type{T}, backend) where {T<:AbstractFloat} = NullFrameResponse()
function convert_frame_response_model(model::SeparableGaussianPixelResponse, ::Type{T}, backend) where {T<:AbstractFloat}
    kernel = backend{T}(undef, length(model.kernel))
    copyto!(kernel, T.(Array(model.kernel)))
    return SeparableGaussianPixelResponse{T,typeof(kernel)}(T(model.response_width_px), kernel)
end

validate_dead_time_model(model::NoDeadTime) = model
function validate_dead_time_model(model::NonParalyzableDeadTime)
    model.dead_time >= 0 || throw(InvalidConfiguration("NonParalyzableDeadTime dead_time must be >= 0"))
    return model
end
validate_frame_response_model(model::NullFrameResponse) = model
function validate_frame_response_model(model::SeparableGaussianPixelResponse)
    model.response_width_px > 0 || throw(InvalidConfiguration("SeparableGaussianPixelResponse response_width_px must be > 0"))
    length(model.kernel) > 0 || throw(InvalidConfiguration("SeparableGaussianPixelResponse kernel must not be empty"))
    isodd(length(model.kernel)) || throw(InvalidConfiguration("SeparableGaussianPixelResponse kernel length must be odd"))
    return model
end

validate_frame_detector_sensor(::FrameSensorType) = nothing
function validate_frame_detector_sensor(sensor::CountingSensorType)
    throw(InvalidConfiguration("Detector currently models frame-based sensors only; use a sensor-side readout model for $(detector_sensor_symbol(sensor))"))
end

function resolve_output_precision(bits::Union{Nothing,Int}, output_precision::Union{Nothing,DataType})
    output_precision !== nothing && return output_precision
    bits === 8 && return UInt8
    bits === 16 && return UInt16
    bits === 32 && return UInt32
    bits === 64 && return UInt64
    return nothing
end

function _build_detector(noise::NoiseModel; integration_time::Real, qe::Real,
    psf_sampling::Int, binning::Int, gain::Real, dark_current::Real,
    bits::Union{Nothing,Int}, full_well::Union{Nothing,Real}, sensor::SensorType,
    response_model::FrameResponseModel, output_precision::Union{Nothing,DataType}, background_flux, background_map,
    T::Type{<:AbstractFloat}, backend)
    validate_frame_detector_sensor(sensor)
    full_well_t = full_well === nothing ? nothing : T(full_well)
    flux_model = background_model(background_flux; T=T, backend=backend)
    map_model = background_model(background_map; T=T, backend=backend)
    response = validate_frame_response_model(convert_frame_response_model(response_model, T, backend))
    output_precision_t = resolve_output_precision(bits, output_precision)
    params = DetectorParams{T, typeof(sensor), typeof(response)}(
        T(integration_time),
        T(qe),
        psf_sampling,
        binning,
        T(gain),
        T(dark_current),
        bits,
        full_well_t,
        sensor,
        response,
        output_precision_t,
    )
    frame = backend{T}(undef, 1, 1)
    response_buffer = backend{T}(undef, 1, 1)
    bin_buffer = backend{T}(undef, 1, 1)
    noise_buffer = backend{T}(undef, 1, 1)
    accum_buffer = backend{T}(undef, 1, 1)
    output_buffer = output_precision_t === nothing ? nothing : backend{output_precision_t}(undef, 1, 1)
    fill!(frame, zero(T))
    fill!(response_buffer, zero(T))
    fill!(bin_buffer, zero(T))
    fill!(noise_buffer, zero(T))
    fill!(accum_buffer, zero(T))
    output_buffer === nothing || fill!(output_buffer, zero(eltype(output_buffer)))
    state = DetectorState{T, typeof(frame), typeof(output_buffer)}(
        frame,
        response_buffer,
        bin_buffer,
        noise_buffer,
        accum_buffer,
        output_buffer,
        zero(T),
        true,
    )
    return Detector{typeof(noise), typeof(params), typeof(state), typeof(flux_model), typeof(map_model)}(
        noise,
        params,
        state,
        flux_model,
        map_model,
    )
end

validate_apd_noise(noise::NoiseNone) = noise
validate_apd_noise(noise::NoisePhoton) = noise
validate_apd_noise(noise::NoiseReadout) =
    throw(InvalidConfiguration("APDDetector does not support additive readout noise; use NoiseNone or NoisePhoton"))
validate_apd_noise(noise::NoisePhotonReadout) =
    throw(InvalidConfiguration("APDDetector does not support frame-style readout noise; use NoisePhoton"))

function _build_apd_detector(noise::NoiseModel; integration_time::Real, qe::Real, gain::Real,
    dark_count_rate::Real, dead_time_model::CountingDeadTimeModel, sensor::APDSensor, output_precision::Union{Nothing,DataType},
    layout::Symbol, channel_gain_map,
    T::Type{<:AbstractFloat}, backend)
    gain >= 0 || throw(InvalidConfiguration("APDDetector gain must be >= 0"))
    dark_count_rate >= 0 || throw(InvalidConfiguration("APDDetector dark_count_rate must be >= 0"))
    converted = convert_noise(noise, T)
    validated = validate_apd_noise(converted)
    dead_time = validate_dead_time_model(convert_dead_time_model(dead_time_model, T))
    gain_map = channel_gain_map === nothing ? nothing : begin
        g = backend{T}(undef, size(channel_gain_map)...)
        copyto!(g, T.(channel_gain_map))
        g
    end
    params = APDDetectorParams{T,typeof(sensor),typeof(dead_time)}(
        T(integration_time),
        T(qe),
        T(gain),
        T(dark_count_rate),
        dead_time,
        sensor,
        output_precision,
        layout,
    )
    channels = backend{T}(undef, 1, 1)
    noise_buffer = backend{T}(undef, 1, 1)
    output_buffer = output_precision === nothing ? nothing : backend{output_precision}(undef, 1, 1)
    fill!(channels, zero(T))
    fill!(noise_buffer, zero(T))
    output_buffer === nothing || fill!(output_buffer, zero(eltype(output_buffer)))
    state = APDDetectorState{T,typeof(channels),typeof(output_buffer)}(channels, noise_buffer, output_buffer)
    return APDDetector{typeof(validated),typeof(params),typeof(state),typeof(gain_map)}(
        validated, params, state, gain_map)
end

function Detector(; integration_time::Real=1.0, qe::Real=1.0,
    psf_sampling::Int=1, binning::Int=1, noise=NoisePhoton(),
    gain::Real=1.0, dark_current::Real=0.0, bits::Union{Nothing,Int}=nothing,
    full_well::Union{Nothing,Real}=nothing, sensor::SensorType=CCDSensor(),
    response_model::FrameResponseModel=NullFrameResponse(),
    output_precision::Union{Nothing,DataType}=nothing, background_flux=nothing, background_map=nothing,
    T::Type{<:AbstractFloat}=Float64, backend=Array)
    normalized = normalize_noise(noise)
    return Detector(normalized; integration_time=integration_time, qe=qe,
        psf_sampling=psf_sampling, binning=binning, gain=gain,
        dark_current=dark_current, bits=bits, full_well=full_well,
        sensor=sensor, response_model=response_model, output_precision=output_precision,
        background_flux=background_flux, background_map=background_map,
        T=T, backend=backend)
end

function Detector(noise::NoiseModel; integration_time::Real=1.0, qe::Real=1.0,
    psf_sampling::Int=1, binning::Int=1, gain::Real=1.0, dark_current::Real=0.0,
    bits::Union{Nothing,Int}=nothing, full_well::Union{Nothing,Real}=nothing,
    sensor::SensorType=CCDSensor(), response_model::FrameResponseModel=NullFrameResponse(),
    output_precision::Union{Nothing,DataType}=nothing,
    background_flux=nothing, background_map=nothing,
    T::Type{<:AbstractFloat}=Float64, backend=Array)
    converted = convert_noise(noise, T)
    validated = validate_noise(converted)
    return _build_detector(validated; integration_time=integration_time, qe=qe,
        psf_sampling=psf_sampling, binning=binning, gain=gain,
        dark_current=dark_current, bits=bits, full_well=full_well,
        sensor=sensor, response_model=response_model, output_precision=output_precision,
        background_flux=background_flux, background_map=background_map,
        T=T, backend=backend)
end

function APDDetector(; integration_time::Real=1.0, qe::Real=1.0, noise::NoiseModel=NoisePhoton(),
    gain::Real=1.0, dark_count_rate::Real=0.0, output_precision::Union{Nothing,DataType}=nothing,
    layout::Symbol=:channels, channel_gain_map=nothing, dead_time_model::CountingDeadTimeModel=NoDeadTime(),
    T::Type{<:AbstractFloat}=Float64, backend=Array)
    return _build_apd_detector(noise; integration_time=integration_time, qe=qe, gain=gain,
        dark_count_rate=dark_count_rate, dead_time_model=dead_time_model, sensor=APDSensor(), output_precision=output_precision,
        layout=layout, channel_gain_map=channel_gain_map, T=T, backend=backend)
end

function ensure_channel_buffers!(det::APDDetector, dims::Tuple{Int,Int})
    if size(det.state.channels) != dims
        det.state.channels = similar(det.state.channels, dims...)
    end
    if size(det.state.noise_buffer) != dims
        det.state.noise_buffer = similar(det.state.noise_buffer, dims...)
    end
    if det.state.output_buffer !== nothing && size(det.state.output_buffer) != dims
        det.state.output_buffer = similar(det.state.output_buffer, dims...)
        fill!(det.state.output_buffer, zero(eltype(det.state.output_buffer)))
    end
    return det
end

function apply_apd_gain_map!(det::APDDetector)
    isnothing(det.channel_gain_map) && return det.state.channels
    size(det.channel_gain_map) == size(det.state.channels) ||
        throw(DimensionMismatchError("APDDetector channel_gain_map size must match counting-channel size"))
    det.state.channels .*= det.channel_gain_map
    return det.state.channels
end

function apply_apd_dark_counts!(det::APDDetector, exposure_time::Real)
    dark = det.params.dark_count_rate * exposure_time
    dark <= 0 && return det.state.channels
    det.state.channels .+= dark
    return det.state.channels
end

apply_apd_dead_time!(det::APDDetector) = apply_apd_dead_time!(det.params.dead_time_model, det)
apply_apd_dead_time!(::NoDeadTime, det::APDDetector) = det.state.channels
function apply_apd_dead_time!(model::NonParalyzableDeadTime, det::APDDetector)
    exposure_time = det.params.integration_time
    exposure_time > zero(exposure_time) || return det.state.channels
    scale = model.dead_time / exposure_time
    scale <= zero(scale) && return det.state.channels
    @. det.state.channels = det.state.channels / (1 + det.state.channels * scale)
    return det.state.channels
end

apply_apd_counting_noise!(det::APDDetector{NoiseNone}, rng::AbstractRNG) = det.state.channels
function apply_apd_counting_noise!(det::APDDetector{NoisePhoton}, rng::AbstractRNG)
    poisson_noise!(rng, det.state.channels)
    return det.state.channels
end

function apply_apd_gain!(det::APDDetector)
    det.state.channels .*= det.params.gain
    return det.state.channels
end

function write_channel_output!(det::APDDetector)
    output = det.state.output_buffer
    output === nothing && return det.state.channels
    output_eltype = eltype(output)
    if output_eltype <: Integer
        lo = typemin(output_eltype)
        hi = typemax(output_eltype)
        @. output = output_eltype(clamp(round(det.state.channels), lo, hi))
    else
        copyto!(output, det.state.channels)
    end
    return output
end

function capture!(det::APDDetector, channels::AbstractMatrix{T}; rng::AbstractRNG=Random.default_rng()) where {T<:AbstractFloat}
    ensure_channel_buffers!(det, size(channels))
    copyto!(det.state.channels, channels)
    det.state.channels .*= det.params.qe * det.params.integration_time
    apply_apd_gain_map!(det)
    apply_apd_dark_counts!(det, det.params.integration_time)
    apply_apd_counting_noise!(det, rng)
    apply_apd_dead_time!(det)
    apply_apd_gain!(det)
    return write_channel_output!(det)
end

function capture!(det::APDDetector, channels::AbstractMatrix{T}, rng::AbstractRNG) where {T<:AbstractFloat}
    return capture!(det, channels; rng=rng)
end

apply_frame_response!(::NullFrameResponse, det::Detector) = det.state.frame

function apply_frame_response!(model::SeparableGaussianPixelResponse, det::Detector)
    return apply_frame_response!(execution_style(det.state.frame), model, det)
end

function apply_frame_response!(::ScalarCPUStyle, model::SeparableGaussianPixelResponse, det::Detector)
    frame = det.state.frame
    scratch = det.state.response_buffer
    kernel = model.kernel
    n, m = size(frame)
    radius = fld(length(kernel), 2)
    @inbounds for i in 1:n, j in 1:m
        acc = zero(eltype(frame))
        for kk in eachindex(kernel)
            jj = clamp(j + kk - radius - 1, 1, m)
            acc += kernel[kk] * frame[i, jj]
        end
        scratch[i, j] = acc
    end
    @inbounds for i in 1:n, j in 1:m
        acc = zero(eltype(frame))
        for kk in eachindex(kernel)
            ii = clamp(i + kk - radius - 1, 1, n)
            acc += kernel[kk] * scratch[ii, j]
        end
        frame[i, j] = acc
    end
    return frame
end

function apply_frame_response!(style::AcceleratorStyle, model::SeparableGaussianPixelResponse, det::Detector)
    frame = det.state.frame
    scratch = det.state.response_buffer
    kernel = model.kernel
    n, m = size(frame)
    radius = fld(length(kernel), 2)
    launch_kernel!(style, separable_response_rows_kernel!, scratch, frame, kernel, radius, n, m, length(kernel); ndrange=(n, m))
    launch_kernel!(style, separable_response_cols_kernel!, frame, scratch, kernel, radius, n, m, length(kernel); ndrange=(n, m))
    return frame
end

function fill_frame!(det::Detector, psf::AbstractMatrix{T}, exposure_time::Real) where {T}
    n_in, m_in = size(psf)
    sampling = det.params.psf_sampling
    binning = det.params.binning
    if sampling < 1 || binning < 1
        throw(InvalidConfiguration("psf_sampling and binning must be >= 1"))
    end
    if n_in % sampling != 0 || m_in % sampling != 0
        throw(DimensionMismatchError("psf_sampling must evenly divide input dimensions"))
    end
    n_mid = div(n_in, sampling)
    m_mid = div(m_in, sampling)
    if n_mid % binning != 0 || m_mid % binning != 0
        throw(DimensionMismatchError("binning must evenly divide sampled dimensions"))
    end
    n_out = div(n_mid, binning)
    m_out = div(m_mid, binning)
    ensure_buffers!(det, n_mid, m_mid, n_out, m_out)

    if sampling > 1
        if binning > 1
            bin2d!(det.state.frame, psf, sampling * binning)
        else
            bin2d!(det.state.bin_buffer, psf, sampling)
            det.state.frame .= det.state.bin_buffer
        end
    else
        if binning > 1
            bin2d!(det.state.frame, psf, binning)
        else
            det.state.frame .= psf
        end
    end
    @. det.state.frame *= det.params.qe * exposure_time
    apply_frame_response!(det.params.response_model, det)
    return det.state.frame
end

function fill_frame!(det::Detector, psf::AbstractMatrix{T}) where {T}
    return fill_frame!(det, psf, det.params.integration_time)
end

function capture_signal!(det::Detector{NoiseNone}, psf::AbstractMatrix{T}, rng::AbstractRNG, exposure_time::Real) where {T}
    fill_frame!(det, psf, exposure_time)
    apply_background_flux!(det.background_flux, det, rng, exposure_time)
    return nothing
end

function capture_signal!(det::Detector{NoisePhoton}, psf::AbstractMatrix{T}, rng::AbstractRNG, exposure_time::Real) where {T}
    fill_frame!(det, psf, exposure_time)
    poisson_noise!(rng, det.state.frame)
    apply_background_flux!(det.background_flux, det, rng, exposure_time)
    return nothing
end

function capture_signal!(det::Detector{<:NoiseReadout}, psf::AbstractMatrix{T}, rng::AbstractRNG, exposure_time::Real) where {T}
    fill_frame!(det, psf, exposure_time)
    apply_background_flux!(det.background_flux, det, rng, exposure_time)
    return nothing
end

function capture_signal!(det::Detector{<:NoisePhotonReadout}, psf::AbstractMatrix{T}, rng::AbstractRNG, exposure_time::Real) where {T}
    fill_frame!(det, psf, exposure_time)
    poisson_noise!(rng, det.state.frame)
    apply_background_flux!(det.background_flux, det, rng, exposure_time)
    return nothing
end

apply_background_flux!(::NoBackground, det::Detector, rng::AbstractRNG, exposure_time::Real) = det.state.frame

function apply_background_flux!(background::ScalarBackground, det::Detector, rng::AbstractRNG, exposure_time::Real)
    fill!(det.state.noise_buffer, background.level * exposure_time)
    poisson_noise!(rng, det.state.noise_buffer)
    det.state.frame .+= det.state.noise_buffer
    return det.state.frame
end

function apply_background_flux!(background::BackgroundFrame, det::Detector, rng::AbstractRNG, exposure_time::Real)
    if size(background.map) != size(det.state.frame)
        throw(DimensionMismatchError("background_flux size must match detector frame size"))
    end
    copyto!(det.state.noise_buffer, background.map)
    det.state.noise_buffer .*= exposure_time
    poisson_noise!(rng, det.state.noise_buffer)
    det.state.frame .+= det.state.noise_buffer
    return det.state.frame
end

function apply_dark_current!(det::Detector, rng::AbstractRNG, exposure_time::Real)
    dark_signal = det.params.dark_current * exposure_time
    if dark_signal <= 0
        return det.state.frame
    end
    fill!(det.state.noise_buffer, dark_signal)
    poisson_noise!(rng, det.state.noise_buffer)
    det.state.frame .+= det.state.noise_buffer
    return det.state.frame
end

function sensor_saturation_limit(det::Detector)
    return sensor_saturation_limit(det.params.sensor, det)
end

sensor_saturation_limit(::FrameSensorType, det::Detector) = det.params.full_well
function sensor_saturation_limit(sensor::SAPHIRASensor, det::Detector)
    full_well = det.params.full_well
    full_well === nothing && return nothing
    return full_well / sensor.avalanche_gain
end

function apply_saturation!(det::Detector)
    full_well = sensor_saturation_limit(det)
    full_well === nothing && return det.state.frame
    clamp!(det.state.frame, zero(eltype(det.state.frame)), full_well)
    return det.state.frame
end

apply_sensor_statistics!(sensor::FrameSensorType, det::Detector, rng::AbstractRNG) = det.state.frame

function apply_sensor_statistics!(sensor::CCDSensor, det::Detector, rng::AbstractRNG)
    rate = sensor.clock_induced_charge_rate * det.params.integration_time
    rate <= zero(rate) && return det.state.frame
    fill!(det.state.noise_buffer, rate)
    poisson_noise!(rng, det.state.noise_buffer)
    det.state.frame .+= det.state.noise_buffer
    return det.state.frame
end

function apply_sensor_statistics!(sensor::CMOSSensor, det::Detector, rng::AbstractRNG)
    sigma = sensor.column_readout_sigma
    sigma <= zero(sigma) && return det.state.frame
    randn_backend!(rng, det.state.noise_buffer)
    return apply_cmos_column_noise!(execution_style(det.state.frame), det.state.frame, det.state.noise_buffer, sigma)
end

function apply_sensor_statistics!(sensor::InGaAsSensor, det::Detector, rng::AbstractRNG)
    rate = sensor.glow_rate * det.params.integration_time
    rate <= zero(rate) && return det.state.frame
    fill!(det.state.noise_buffer, rate)
    poisson_noise!(rng, det.state.noise_buffer)
    det.state.frame .+= det.state.noise_buffer
    return det.state.frame
end

function apply_avalanche_excess_noise!(factor, det::Detector, rng::AbstractRNG)
    factor <= one(factor) && return det.state.frame
    randn_backend!(rng, det.state.noise_buffer)
    scale2 = factor * factor - one(factor)
    zero_t = zero(eltype(det.state.frame))
    @. det.state.frame += sqrt(max(scale2 * det.state.frame, zero_t)) * det.state.noise_buffer
    return det.state.frame
end

apply_sensor_statistics!(sensor::EMCCDSensor, det::Detector, rng::AbstractRNG) =
    apply_avalanche_excess_noise!(sensor.excess_noise_factor, det, rng)

function apply_sensor_statistics!(sensor::SAPHIRASensor, det::Detector, rng::AbstractRNG)
    rate = sensor.glow_rate * det.params.integration_time
    if rate > zero(rate)
        fill!(det.state.noise_buffer, rate)
        poisson_noise!(rng, det.state.noise_buffer)
        det.state.frame .+= det.state.noise_buffer
    end
    return apply_avalanche_excess_noise!(sensor.excess_noise_factor, det, rng)
end

function apply_cmos_column_noise!(::ScalarCPUStyle, frame, noise, sigma)
    n, m = size(frame)
    @inbounds for j in 1:m
        offset = sigma * noise[1, j]
        for i in 1:n
            frame[i, j] += offset
        end
    end
    return frame
end

function apply_cmos_column_noise!(style::AcceleratorStyle, frame, noise, sigma)
    n, m = size(frame)
    launch_kernel!(style, add_column_noise_kernel!, frame, noise, sigma, n, m; ndrange=(n, m))
    return frame
end

apply_pre_readout_gain!(::CCDSensor, det::Detector) = det.state.frame
apply_pre_readout_gain!(::CMOSSensor, det::Detector) = det.state.frame
apply_pre_readout_gain!(::InGaAsSensor, det::Detector) = det.state.frame
function apply_pre_readout_gain!(::EMCCDSensor, det::Detector)
    det.state.frame .*= det.params.gain
    return det.state.frame
end
function apply_pre_readout_gain!(sensor::SAPHIRASensor, det::Detector)
    det.state.frame .*= sensor.avalanche_gain
    return det.state.frame
end

apply_readout_noise!(det::Detector{NoiseNone}, rng::AbstractRNG) = det.state.frame
apply_readout_noise!(det::Detector{NoisePhoton}, rng::AbstractRNG) = det.state.frame
function apply_readout_noise!(det::Detector{<:NoiseReadout}, rng::AbstractRNG)
    randn_backend!(rng, det.state.noise_buffer)
    det.state.frame .+= det.noise.sigma .* det.state.noise_buffer
    return det.state.frame
end
function apply_readout_noise!(det::Detector{<:NoisePhotonReadout}, rng::AbstractRNG)
    randn_backend!(rng, det.state.noise_buffer)
    det.state.frame .+= det.noise.sigma .* det.state.noise_buffer
    return det.state.frame
end

apply_post_readout_gain!(::EMCCDSensor, det::Detector) = det.state.frame
function apply_post_readout_gain!(::CCDSensor, det::Detector)
    det.state.frame .*= det.params.gain
    return det.state.frame
end
function apply_post_readout_gain!(::CMOSSensor, det::Detector)
    det.state.frame .*= det.params.gain
    return det.state.frame
end
function apply_post_readout_gain!(::InGaAsSensor, det::Detector)
    det.state.frame .*= det.params.gain
    return det.state.frame
end
function apply_post_readout_gain!(::SAPHIRASensor, det::Detector)
    det.state.frame .*= det.params.gain
    return det.state.frame
end

function apply_quantization!(det::Detector)
    bits = det.params.bits
    bits === nothing && return det.state.frame
    levels = exp2(eltype(det.state.frame)(bits))
    full_well = det.params.full_well
    if full_well === nothing
        peak = maximum(det.state.frame)
        if peak > 0
            det.state.frame .*= levels / peak
        end
    else
        det.state.frame .*= (levels - one(levels)) / full_well
        clamp!(det.state.frame, zero(eltype(det.state.frame)), levels - one(levels))
    end
    return det.state.frame
end

subtract_background_map!(::NoBackground, det::Detector) = det.state.frame

function subtract_background_map!(background::ScalarBackground, det::Detector)
    det.state.frame .-= background.level
    return det.state.frame
end

function subtract_background_map!(background::BackgroundFrame, det::Detector)
    if size(background.map) != size(det.state.frame)
        throw(DimensionMismatchError("background_map size must match detector frame size"))
    end
    det.state.frame .-= background.map
    return det.state.frame
end

function finalize_capture!(det::Detector, rng::AbstractRNG, exposure_time::Real)
    apply_dark_current!(det, rng, exposure_time)
    apply_saturation!(det)
    apply_sensor_statistics!(det.params.sensor, det, rng)
    apply_pre_readout_gain!(det.params.sensor, det)
    apply_readout_noise!(det, rng)
    apply_post_readout_gain!(det.params.sensor, det)
    apply_quantization!(det)
    subtract_background_map!(det.background_map, det)
    return det.state.frame
end

function write_output_frame!(det::Detector)
    output = det.state.output_buffer
    output === nothing && return det.state.frame
    output_eltype = eltype(output)
    if output_eltype <: Integer
        lo = typemin(output_eltype)
        hi = typemax(output_eltype)
        @. output = output_eltype(clamp(round(det.state.frame), lo, hi))
    else
        copyto!(output, det.state.frame)
    end
    return output
end

function capture!(det::Detector, psf::AbstractMatrix{T}, rng::AbstractRNG) where {T}
    capture_signal!(det, psf, rng, det.params.integration_time)
    finalize_capture!(det, rng, det.params.integration_time)
    return write_output_frame!(det)
end

function _require_batched_detector_compat(det::Detector, cube::AbstractArray, scratch::AbstractArray)
    size(cube) == size(scratch) ||
        throw(DimensionMismatchError("batched detector scratch must match cube size"))
    ndims(cube) == 3 || throw(DimensionMismatchError("batched detector input must be 3D"))
    det.params.psf_sampling == 1 ||
        throw(InvalidConfiguration("batched detector capture currently requires psf_sampling == 1"))
    det.params.binning == 1 ||
        throw(InvalidConfiguration("batched detector capture currently requires binning == 1"))
    det.params.output_precision === nothing ||
        throw(InvalidConfiguration("batched detector capture currently requires output_precision === nothing"))
    supports_detector_mtf(det) &&
        throw(InvalidConfiguration("batched detector capture currently requires NullFrameResponse()"))
    _require_batched_sensor_compat(det.params.sensor)
    return nothing
end

_require_batched_sensor_compat(::FrameSensorType) = nothing
function _require_batched_sensor_compat(sensor::CMOSSensor)
    sensor.column_readout_sigma <= zero(sensor.column_readout_sigma) ||
        throw(InvalidConfiguration("batched detector capture does not yet support CMOSSensor column_readout_sigma"))
    return nothing
end

_batched_background_map!(::NoBackground, cube::AbstractArray, scratch::AbstractArray) = cube
_batched_background_map!(background::ScalarBackground, cube::AbstractArray, scratch::AbstractArray) = (cube .-= background.level; cube)

function _batched_background_map!(background::BackgroundFrame, cube::AbstractArray, scratch::AbstractArray)
    size(background.map) == size(cube)[1:2] ||
        throw(DimensionMismatchError("background_map size must match detector frame size"))
    cube .-= reshape(background.map, size(background.map)..., 1)
    return cube
end

_batched_background_flux!(::NoBackground, det::Detector, cube::AbstractArray, scratch::AbstractArray, rng::AbstractRNG, exposure_time::Real) = cube

function _batched_background_flux!(background::ScalarBackground, det::Detector, cube::AbstractArray, scratch::AbstractArray,
    rng::AbstractRNG, exposure_time::Real)
    fill!(scratch, background.level * exposure_time)
    poisson_noise!(rng, scratch)
    cube .+= scratch
    return cube
end

function _batched_background_flux!(background::BackgroundFrame, det::Detector, cube::AbstractArray, scratch::AbstractArray,
    rng::AbstractRNG, exposure_time::Real)
    size(background.map) == size(cube)[1:2] ||
        throw(DimensionMismatchError("background_flux size must match detector frame size"))
    scratch .= reshape(background.map, size(background.map)..., 1)
    scratch .*= exposure_time
    poisson_noise!(rng, scratch)
    cube .+= scratch
    return cube
end

function _batched_dark_current!(det::Detector, cube::AbstractArray, scratch::AbstractArray, rng::AbstractRNG, exposure_time::Real)
    dark_signal = det.params.dark_current * exposure_time
    if dark_signal <= 0
        return cube
    end
    fill!(scratch, dark_signal)
    poisson_noise!(rng, scratch)
    cube .+= scratch
    return cube
end

_batched_pre_readout_gain!(::CCDSensor, det::Detector, cube::AbstractArray) = cube
_batched_pre_readout_gain!(::CMOSSensor, det::Detector, cube::AbstractArray) = cube
_batched_pre_readout_gain!(::InGaAsSensor, det::Detector, cube::AbstractArray) = cube
_batched_pre_readout_gain!(::EMCCDSensor, det::Detector, cube::AbstractArray) = (cube .*= det.params.gain; cube)
_batched_pre_readout_gain!(sensor::SAPHIRASensor, det::Detector, cube::AbstractArray) = (cube .*= sensor.avalanche_gain; cube)

_batched_sensor_statistics!(sensor::FrameSensorType, det::Detector, cube::AbstractArray, scratch::AbstractArray, rng::AbstractRNG) = cube
function _batched_sensor_statistics!(sensor::CCDSensor, det::Detector, cube::AbstractArray, scratch::AbstractArray, rng::AbstractRNG)
    rate = sensor.clock_induced_charge_rate * det.params.integration_time
    rate <= zero(rate) && return cube
    fill!(scratch, rate)
    poisson_noise!(rng, scratch)
    cube .+= scratch
    return cube
end
function _batched_sensor_statistics!(sensor::InGaAsSensor, det::Detector, cube::AbstractArray, scratch::AbstractArray, rng::AbstractRNG)
    rate = sensor.glow_rate * det.params.integration_time
    rate <= zero(rate) && return cube
    fill!(scratch, rate)
    poisson_noise!(rng, scratch)
    cube .+= scratch
    return cube
end
function _batched_avalanche_excess_noise!(factor, cube::AbstractArray, scratch::AbstractArray, rng::AbstractRNG)
    factor <= one(factor) && return cube
    randn_backend!(rng, scratch)
    scale2 = factor * factor - one(factor)
    zero_t = zero(eltype(cube))
    @. cube += sqrt(max(scale2 * cube, zero_t)) * scratch
    return cube
end
_batched_sensor_statistics!(sensor::EMCCDSensor, det::Detector, cube::AbstractArray, scratch::AbstractArray, rng::AbstractRNG) =
    _batched_avalanche_excess_noise!(sensor.excess_noise_factor, cube, scratch, rng)
function _batched_sensor_statistics!(sensor::SAPHIRASensor, det::Detector, cube::AbstractArray, scratch::AbstractArray, rng::AbstractRNG)
    rate = sensor.glow_rate * det.params.integration_time
    if rate > zero(rate)
        fill!(scratch, rate)
        poisson_noise!(rng, scratch)
        cube .+= scratch
    end
    return _batched_avalanche_excess_noise!(sensor.excess_noise_factor, cube, scratch, rng)
end

_batched_post_readout_gain!(::EMCCDSensor, det::Detector, cube::AbstractArray) = cube
_batched_post_readout_gain!(::CCDSensor, det::Detector, cube::AbstractArray) = (cube .*= det.params.gain; cube)
_batched_post_readout_gain!(::CMOSSensor, det::Detector, cube::AbstractArray) = (cube .*= det.params.gain; cube)
_batched_post_readout_gain!(::InGaAsSensor, det::Detector, cube::AbstractArray) = (cube .*= det.params.gain; cube)
_batched_post_readout_gain!(::SAPHIRASensor, det::Detector, cube::AbstractArray) = (cube .*= det.params.gain; cube)

_batched_readout_noise!(det::Detector{NoiseNone}, cube::AbstractArray, scratch::AbstractArray, rng::AbstractRNG) = cube
_batched_readout_noise!(det::Detector{NoisePhoton}, cube::AbstractArray, scratch::AbstractArray, rng::AbstractRNG) = cube

function _batched_readout_noise!(det::Detector{<:NoiseReadout}, cube::AbstractArray, scratch::AbstractArray, rng::AbstractRNG)
    randn_backend!(rng, scratch)
    cube .+= det.noise.sigma .* scratch
    return cube
end

function _batched_readout_noise!(det::Detector{<:NoisePhotonReadout}, cube::AbstractArray, scratch::AbstractArray, rng::AbstractRNG)
    randn_backend!(rng, scratch)
    cube .+= det.noise.sigma .* scratch
    return cube
end

function _batched_quantization!(det::Detector, cube::AbstractArray)
    bits = det.params.bits
    bits === nothing && return cube
    levels = exp2(eltype(cube)(bits))
    full_well = det.params.full_well
    if full_well === nothing
        peak = maximum(cube)
        if peak > 0
            cube .*= levels / peak
        end
    else
        cube .*= (levels - one(levels)) / full_well
        clamp!(cube, zero(eltype(cube)), levels - one(levels))
    end
    return cube
end

function capture_stack!(det::Detector, cube::AbstractArray{T,3}, scratch::AbstractArray{T,3};
    rng::AbstractRNG=Random.default_rng()) where {T<:AbstractFloat}
    _require_batched_detector_compat(det, cube, scratch)
    exposure_time = det.params.integration_time
    cube .*= det.params.qe * exposure_time
    if det isa Detector{NoisePhoton} || det isa Detector{<:NoisePhotonReadout}
        poisson_noise!(rng, cube)
    end
    _batched_background_flux!(det.background_flux, det, cube, scratch, rng, exposure_time)
    _batched_dark_current!(det, cube, scratch, rng, exposure_time)
    apply_saturation!(det, cube)
    _batched_sensor_statistics!(det.params.sensor, det, cube, scratch, rng)
    _batched_pre_readout_gain!(det.params.sensor, det, cube)
    _batched_readout_noise!(det, cube, scratch, rng)
    _batched_post_readout_gain!(det.params.sensor, det, cube)
    _batched_quantization!(det, cube)
    _batched_background_map!(det.background_map, cube, scratch)
    return cube
end

function apply_saturation!(det::Detector, cube::AbstractArray)
    full_well = det.params.full_well
    full_well === nothing && return cube
    clamp!(cube, zero(eltype(cube)), full_well)
    return cube
end

function capture!(det::Detector, psf::AbstractMatrix{T}; rng::AbstractRNG=Random.default_rng(),
    sample_time::Union{Nothing,Real}=nothing) where {T}
    if sample_time === nothing
        return capture!(det, psf, rng)
    end
    dt = eltype(det.state.frame)(sample_time)
    if dt <= 0
        throw(InvalidConfiguration("sample_time must be > 0"))
    end
    if det.params.integration_time < dt
        throw(InvalidConfiguration("sample_time must be <= detector integration_time"))
    end
    capture_signal!(det, psf, rng, dt)
    if size(det.state.accum_buffer) != size(det.state.frame)
        det.state.accum_buffer = similar(det.state.frame, size(det.state.frame)...)
        fill!(det.state.accum_buffer, zero(eltype(det.state.accum_buffer)))
    end
    det.state.accum_buffer .+= det.state.frame
    det.state.integrated_time += dt
    det.state.readout_ready = false
    if det.state.integrated_time + eps(dt) >= det.params.integration_time
        det.state.frame .= det.state.accum_buffer
        finalize_capture!(det, rng, det.state.integrated_time)
        fill!(det.state.accum_buffer, zero(eltype(det.state.accum_buffer)))
        det.state.integrated_time = zero(det.state.integrated_time)
        det.state.readout_ready = true
    end
    return write_output_frame!(det)
end

function ensure_buffers!(det::Detector, n_mid::Int, m_mid::Int, n_out::Int, m_out::Int)
    if size(det.state.frame) != (n_out, m_out)
        det.state.frame = similar(det.state.frame, n_out, m_out)
    end
    if size(det.state.response_buffer) != (n_out, m_out)
        det.state.response_buffer = similar(det.state.response_buffer, n_out, m_out)
    end
    if size(det.state.bin_buffer) != (n_mid, m_mid)
        det.state.bin_buffer = similar(det.state.bin_buffer, n_mid, m_mid)
    end
    if size(det.state.noise_buffer) != (n_out, m_out)
        det.state.noise_buffer = similar(det.state.noise_buffer, n_out, m_out)
    end
    if size(det.state.accum_buffer) != (n_out, m_out)
        det.state.accum_buffer = similar(det.state.accum_buffer, n_out, m_out)
        fill!(det.state.accum_buffer, zero(eltype(det.state.accum_buffer)))
        det.state.integrated_time = zero(det.state.integrated_time)
        det.state.readout_ready = true
    end
    if det.state.output_buffer !== nothing && size(det.state.output_buffer) != (n_out, m_out)
        det.state.output_buffer = similar(det.state.output_buffer, n_out, m_out)
        fill!(det.state.output_buffer, zero(eltype(det.state.output_buffer)))
    end
    return det
end
