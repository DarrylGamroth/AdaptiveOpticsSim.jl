using Random

abstract type NoiseModel end
abstract type SensorType end
abstract type BackgroundModel end
struct NoiseNone <: NoiseModel end
struct NoisePhoton <: NoiseModel end
struct NoiseReadout{T<:AbstractFloat} <: NoiseModel
    sigma::T
end
struct NoisePhotonReadout{T<:AbstractFloat} <: NoiseModel
    sigma::T
end
struct CCDSensor <: SensorType end
struct CMOSSensor <: SensorType end
struct EMCCDSensor <: SensorType end
struct NoBackground <: BackgroundModel end
struct ScalarBackground{T<:AbstractFloat} <: BackgroundModel
    level::T
end
struct BackgroundFrame{T<:AbstractFloat,A<:AbstractMatrix{T}} <: BackgroundModel
    map::A
end

NoiseReadout(sigma::Real) = NoiseReadout{Float64}(float(sigma))
NoisePhotonReadout(sigma::Real) = NoisePhotonReadout{Float64}(float(sigma))

struct DetectorParams{T<:AbstractFloat,S<:SensorType}
    integration_time::T
    qe::T
    psf_sampling::Int
    binning::Int
    gain::T
    dark_current::T
    bits::Union{Nothing,Int}
    full_well::Union{Nothing,T}
    sensor::S
    output_precision::Union{Nothing,DataType}
end

mutable struct DetectorState{T<:AbstractFloat,A<:AbstractMatrix{T},O}
    frame::A
    bin_buffer::A
    noise_buffer::A
    accum_buffer::A
    output_buffer::O
    integrated_time::T
    readout_ready::Bool
end

struct Detector{N<:NoiseModel,P<:DetectorParams,S<:DetectorState,BF<:BackgroundModel,BM<:BackgroundModel} <: AbstractDetector
    noise::N
    params::P
    state::S
    background_flux::BF
    background_map::BM
end

readout_ready(det::Detector) = det.state.readout_ready
output_frame(det::Detector) = det.state.output_buffer === nothing ? det.state.frame : det.state.output_buffer

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
    output_precision::Union{Nothing,DataType}, background_flux, background_map,
    T::Type{<:AbstractFloat}, backend)
    full_well_t = full_well === nothing ? nothing : T(full_well)
    flux_model = background_model(background_flux; T=T, backend=backend)
    map_model = background_model(background_map; T=T, backend=backend)
    output_precision_t = resolve_output_precision(bits, output_precision)
    params = DetectorParams{T, typeof(sensor)}(
        T(integration_time),
        T(qe),
        psf_sampling,
        binning,
        T(gain),
        T(dark_current),
        bits,
        full_well_t,
        sensor,
        output_precision_t,
    )
    frame = backend{T}(undef, 1, 1)
    bin_buffer = backend{T}(undef, 1, 1)
    noise_buffer = backend{T}(undef, 1, 1)
    accum_buffer = backend{T}(undef, 1, 1)
    output_buffer = output_precision_t === nothing ? nothing : backend{output_precision_t}(undef, 1, 1)
    fill!(frame, zero(T))
    fill!(bin_buffer, zero(T))
    fill!(noise_buffer, zero(T))
    fill!(accum_buffer, zero(T))
    output_buffer === nothing || fill!(output_buffer, zero(eltype(output_buffer)))
    state = DetectorState{T, typeof(frame), typeof(output_buffer)}(
        frame,
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

function Detector(; integration_time::Real=1.0, qe::Real=1.0,
    psf_sampling::Int=1, binning::Int=1, noise=NoisePhoton(),
    gain::Real=1.0, dark_current::Real=0.0, bits::Union{Nothing,Int}=nothing,
    full_well::Union{Nothing,Real}=nothing, sensor::SensorType=CCDSensor(),
    output_precision::Union{Nothing,DataType}=nothing, background_flux=nothing, background_map=nothing,
    T::Type{<:AbstractFloat}=Float64, backend=Array)
    normalized = normalize_noise(noise)
    return Detector(normalized; integration_time=integration_time, qe=qe,
        psf_sampling=psf_sampling, binning=binning, gain=gain,
        dark_current=dark_current, bits=bits, full_well=full_well,
        sensor=sensor, output_precision=output_precision,
        background_flux=background_flux, background_map=background_map,
        T=T, backend=backend)
end

function Detector(noise::NoiseModel; integration_time::Real=1.0, qe::Real=1.0,
    psf_sampling::Int=1, binning::Int=1, gain::Real=1.0, dark_current::Real=0.0,
    bits::Union{Nothing,Int}=nothing, full_well::Union{Nothing,Real}=nothing,
    sensor::SensorType=CCDSensor(), output_precision::Union{Nothing,DataType}=nothing,
    background_flux=nothing, background_map=nothing,
    T::Type{<:AbstractFloat}=Float64, backend=Array)
    converted = convert_noise(noise, T)
    validated = validate_noise(converted)
    return _build_detector(validated; integration_time=integration_time, qe=qe,
        psf_sampling=psf_sampling, binning=binning, gain=gain,
        dark_current=dark_current, bits=bits, full_well=full_well,
        sensor=sensor, output_precision=output_precision,
        background_flux=background_flux, background_map=background_map,
        T=T, backend=backend)
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
        bin2d!(det.state.bin_buffer, psf, sampling)
        if binning > 1
            bin2d!(det.state.frame, det.state.bin_buffer, binning)
        else
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

function apply_saturation!(det::Detector)
    full_well = det.params.full_well
    full_well === nothing && return det.state.frame
    clamp!(det.state.frame, zero(eltype(det.state.frame)), full_well)
    return det.state.frame
end

apply_pre_readout_gain!(::CCDSensor, det::Detector) = det.state.frame
apply_pre_readout_gain!(::CMOSSensor, det::Detector) = det.state.frame
function apply_pre_readout_gain!(::EMCCDSensor, det::Detector)
    det.state.frame .*= det.params.gain
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
