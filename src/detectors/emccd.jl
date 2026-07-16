abstract type AbstractEMGainModel end
abstract type AbstractEMCCDOperatingMode end
abstract type AbstractEMCCDOutputPath end
abstract type AbstractEMCCDAcquisitionMode end

struct ExcessNoiseApproximation <: AbstractEMGainModel end
struct LinearEMMode <: AbstractEMCCDOperatingMode end
struct EMOutput <: AbstractEMCCDOutputPath end
struct ConventionalOutput <: AbstractEMCCDOutputPath end
struct SequentialAcquisition <: AbstractEMCCDAcquisitionMode end

"""
    FrameTransferAcquisition(; transfer_time=0, T=Float64)

EMCCD acquisition timing in which the image area is rapidly shifted into a
storage area and storage readout overlaps the next integration. This affects
latency and steady-state frame period only; it does not alter detector optical
response, charge multiplication, or noise.
"""
struct FrameTransferAcquisition{T<:AbstractFloat} <:
    AbstractEMCCDAcquisitionMode
    transfer_time::T
end

FrameTransferAcquisition(transfer_time::Real) =
    validate_emccd_acquisition_mode(
        FrameTransferAcquisition{Float64}(float(transfer_time)))
FrameTransferAcquisition(; transfer_time::Real=0.0,
    T::Type{<:AbstractFloat}=Float64) =
    validate_emccd_acquisition_mode(FrameTransferAcquisition{T}(T(transfer_time)))

struct PhotonCountingEMMode{T<:AbstractFloat} <: AbstractEMCCDOperatingMode
    threshold::T
    detection_efficiency::T
end

function PhotonCountingEMMode(; threshold::Real, detection_efficiency::Real=1.0,
    T::Type{<:AbstractFloat}=Float64)
    threshold >= 0 || throw(InvalidConfiguration("PhotonCountingEMMode threshold must be >= 0"))
    zero(T) <= T(detection_efficiency) <= one(T) ||
        throw(InvalidConfiguration("PhotonCountingEMMode detection_efficiency must be in [0, 1]"))
    return PhotonCountingEMMode{T}(T(threshold), T(detection_efficiency))
end

struct StochasticMultiplicationRegister{T<:AbstractFloat} <: AbstractEMGainModel
    register_noise_factor::T
end

StochasticMultiplicationRegister(register_noise_factor::Real) = StochasticMultiplicationRegister{Float64}(float(register_noise_factor))

struct EMCCDSensor{T<:AbstractFloat,M<:AbstractEMGainModel,
    O<:AbstractEMCCDOperatingMode,P<:AbstractEMCCDOutputPath,
    A<:AbstractEMCCDAcquisitionMode} <: AvalancheFrameSensorType
    excess_noise_factor::T
    clock_induced_charge_per_frame::T
    multiplication_model::M
    register_full_well::Union{Nothing,T}
    operating_mode::O
    output_path::P
    acquisition_mode::A
    em_gain_range::Tuple{T,T}
    readout_rate_hz::Union{Nothing,T}
end

function EMCCDSensor(; excess_noise_factor::Real=1.0,
    clock_induced_charge_per_frame::Real=0.0,
    multiplication_model::AbstractEMGainModel=ExcessNoiseApproximation(),
    register_full_well::Union{Nothing,Real}=nothing,
    operating_mode::AbstractEMCCDOperatingMode=LinearEMMode(),
    output_path::AbstractEMCCDOutputPath=EMOutput(),
    acquisition_mode::AbstractEMCCDAcquisitionMode=SequentialAcquisition(),
    em_gain_range::Tuple{<:Real,<:Real}=(1.0, Inf),
    readout_rate_hz::Union{Nothing,Real}=nothing, T::Type{<:AbstractFloat}=Float64)
    excess_noise_factor >= 1 || throw(InvalidConfiguration("EMCCDSensor excess_noise_factor must be >= 1"))
    clock_induced_charge_per_frame >= 0 || throw(InvalidConfiguration(
        "EMCCDSensor clock_induced_charge_per_frame must be >= 0"))
    reg_full_well = register_full_well === nothing ? nothing : T(register_full_well)
    reg_full_well === nothing || reg_full_well > zero(reg_full_well) ||
        throw(InvalidConfiguration("EMCCDSensor register_full_well must be > 0"))
    em_range = (T(em_gain_range[1]), T(em_gain_range[2]))
    em_range[1] > zero(T) || throw(InvalidConfiguration("EMCCDSensor em_gain_range lower bound must be > 0"))
    em_range[2] >= em_range[1] ||
        throw(InvalidConfiguration("EMCCDSensor em_gain_range upper bound must be >= lower bound"))
    readout_rate = readout_rate_hz === nothing ? nothing : T(readout_rate_hz)
    readout_rate === nothing || readout_rate > zero(readout_rate) ||
        throw(InvalidConfiguration("EMCCDSensor readout_rate_hz must be > 0"))
    converted_model = convert_em_gain_model(multiplication_model, T)
    validated_model = validate_em_gain_model(converted_model)
    converted_mode = convert_emccd_operating_mode(operating_mode, T)
    validated_mode = validate_emccd_operating_mode(converted_mode)
    converted_acquisition = convert_emccd_acquisition_mode(acquisition_mode, T)
    validated_acquisition = validate_emccd_acquisition_mode(converted_acquisition)
    return EMCCDSensor{T,typeof(validated_model),typeof(validated_mode),
        typeof(output_path),typeof(validated_acquisition)}(
        T(excess_noise_factor), T(clock_induced_charge_per_frame), validated_model,
        reg_full_well, validated_mode,
        output_path, validated_acquisition, em_range, readout_rate)
end

detector_sensor_symbol(::EMCCDSensor) = :emccd
supports_clock_induced_charge(::EMCCDSensor) = true
supports_photon_number_resolving(sensor::EMCCDSensor) = supports_emccd_photon_counting(sensor.operating_mode)
supports_emccd_photon_counting(::AbstractEMCCDOperatingMode) = false
supports_emccd_photon_counting(::PhotonCountingEMMode) = true
configured_cic_rate(sensor::EMCCDSensor, ::Type{T}) where {T<:AbstractFloat} =
    T(sensor.clock_induced_charge_per_frame)
is_excess_noise_model(::AbstractEMGainModel) = false
is_excess_noise_model(::ExcessNoiseApproximation) = true

acquisition_mode_symbol(sensor::EMCCDSensor) =
    emccd_acquisition_mode_symbol(sensor.acquisition_mode)
emccd_acquisition_mode_symbol(::SequentialAcquisition) = :sequential
emccd_acquisition_mode_symbol(::FrameTransferAcquisition) = :frame_transfer

frame_transfer_time(sensor::EMCCDSensor, ::Type{T}) where {T<:AbstractFloat} =
    emccd_frame_transfer_time(sensor.acquisition_mode, T)
emccd_frame_transfer_time(::SequentialAcquisition,
    ::Type{T}) where {T<:AbstractFloat} = nothing
emccd_frame_transfer_time(mode::FrameTransferAcquisition,
    ::Type{T}) where {T<:AbstractFloat} = T(mode.transfer_time)

_emccd_readout_pixels(frame_size::Tuple{Int,Int}, ::Nothing) =
    frame_size[1] * frame_size[2]
_emccd_readout_pixels(frame_size::Tuple{Int,Int}, window::FrameWindow) =
    length(window.rows) * length(window.cols)

function sampling_read_time(sensor::EMCCDSensor, frame_size::Tuple{Int,Int},
    window::Union{Nothing,FrameWindow}, ::Type{T}) where {T<:AbstractFloat}
    sensor.readout_rate_hz === nothing && return nothing
    return T(_emccd_readout_pixels(frame_size, window)) /
        T(sensor.readout_rate_hz)
end

sampling_read_time(::EMCCDSensor, ::Type{T}) where {T<:AbstractFloat} = nothing

function sampling_wallclock_time(sensor::EMCCDSensor, integration_time,
    frame_size::Tuple{Int,Int}, window::Union{Nothing,FrameWindow},
    ::Type{T}) where {T<:AbstractFloat}
    readout_time = sampling_read_time(sensor, frame_size, window, T)
    readout_time === nothing && return nothing
    return emccd_first_frame_latency(sensor.acquisition_mode,
        T(integration_time), readout_time)
end

sampling_wallclock_time(::EMCCDSensor, integration_time,
    ::Type{T}) where {T<:AbstractFloat} = nothing

emccd_first_frame_latency(::SequentialAcquisition, integration_time,
    readout_time) = integration_time + readout_time
emccd_first_frame_latency(mode::FrameTransferAcquisition, integration_time,
    readout_time) = integration_time + mode.transfer_time + readout_time

function steady_state_frame_period(sensor::EMCCDSensor, integration_time,
    frame_size::Tuple{Int,Int}, window::Union{Nothing,FrameWindow},
    ::Type{T}) where {T<:AbstractFloat}
    readout_time = sampling_read_time(sensor, frame_size, window, T)
    readout_time === nothing && return nothing
    return emccd_steady_state_frame_period(sensor.acquisition_mode,
        T(integration_time), readout_time)
end

emccd_steady_state_frame_period(::SequentialAcquisition, integration_time,
    readout_time) = integration_time + readout_time
emccd_steady_state_frame_period(mode::FrameTransferAcquisition,
    integration_time, readout_time) =
    max(integration_time, readout_time) + mode.transfer_time

convert_em_gain_model(::ExcessNoiseApproximation, ::Type{T}) where {T<:AbstractFloat} = ExcessNoiseApproximation()
convert_em_gain_model(model::StochasticMultiplicationRegister, ::Type{T}) where {T<:AbstractFloat} =
    StochasticMultiplicationRegister{T}(T(model.register_noise_factor))
convert_emccd_operating_mode(::LinearEMMode, ::Type{T}) where {T<:AbstractFloat} = LinearEMMode()
convert_emccd_operating_mode(mode::PhotonCountingEMMode, ::Type{T}) where {T<:AbstractFloat} =
    PhotonCountingEMMode(threshold=T(mode.threshold), detection_efficiency=T(mode.detection_efficiency), T=T)
validate_em_gain_model(::ExcessNoiseApproximation) = ExcessNoiseApproximation()
validate_emccd_operating_mode(::LinearEMMode) = LinearEMMode()
validate_emccd_operating_mode(mode::PhotonCountingEMMode) = mode
convert_emccd_acquisition_mode(::SequentialAcquisition, ::Type{T}) where
    {T<:AbstractFloat} = SequentialAcquisition()
convert_emccd_acquisition_mode(mode::FrameTransferAcquisition,
    ::Type{T}) where {T<:AbstractFloat} =
    FrameTransferAcquisition{T}(T(mode.transfer_time))
validate_emccd_acquisition_mode(::SequentialAcquisition) = SequentialAcquisition()

function validate_emccd_acquisition_mode(mode::FrameTransferAcquisition)
    mode.transfer_time >= zero(mode.transfer_time) || throw(InvalidConfiguration(
        "FrameTransferAcquisition transfer_time must be >= 0"))
    return mode
end

function validate_em_gain_model(model::StochasticMultiplicationRegister)
    model.register_noise_factor >= zero(model.register_noise_factor) ||
        throw(InvalidConfiguration("StochasticMultiplicationRegister register_noise_factor must be >= 0"))
    return model
end

function apply_sensor_statistics!(sensor::EMCCDSensor, det::Detector,
    rng::AbstractRNG, exposure_time::Real)
    rate = effective_cic_rate(det)
    add_poisson_rate!(det.state.frame, det, rng, rate)
    return det.state.frame
end

function apply_pre_readout_gain!(sensor::EMCCDSensor, det::Detector, rng::AbstractRNG)
    return apply_pre_readout_gain!(sensor.output_path, sensor, det, rng)
end

function apply_pre_readout_gain!(::EMOutput, sensor::EMCCDSensor, det::Detector, rng::AbstractRNG)
    return apply_em_register_gain!(sensor.multiplication_model, sensor, det, rng)
end

apply_pre_readout_gain!(::ConventionalOutput, sensor::EMCCDSensor, det::Detector, rng::AbstractRNG) = det.state.frame

function apply_em_register_gain!(::ExcessNoiseApproximation, sensor::EMCCDSensor, det::Detector, rng::AbstractRNG)
    apply_avalanche_excess_noise!(sensor.excess_noise_factor, det, rng)
    det.state.frame .*= det.params.gain
    return det.state.frame
end

function apply_em_register_gain!(model::StochasticMultiplicationRegister, sensor::EMCCDSensor, det::Detector, rng::AbstractRNG)
    return _apply_stochastic_em_register!(execution_style(det.state.frame), model,
        sensor, det.state.frame, det.state.noise_buffer, det.params.gain, rng)
end

function apply_post_readout_gain!(sensor::EMCCDSensor, det::Detector)
    return det.state.frame
end

@inline function em_register_noise_factor(model::StochasticMultiplicationRegister,
    sensor::EMCCDSensor)
    sensor.excess_noise_factor <= one(sensor.excess_noise_factor) &&
        return model.register_noise_factor
    return max(model.register_noise_factor,
        sqrt(sensor.excess_noise_factor^2 - one(sensor.excess_noise_factor)))
end

@inline function _gamma_unit_scale(rng::AbstractRNG, shape::T) where {T<:AbstractFloat}
    if shape < one(T)
        return _gamma_unit_scale(rng, shape + one(T)) *
            rand(rng, T)^inv(shape)
    end

    d = shape - one(T) / T(3)
    c = inv(sqrt(T(9) * d))
    while true
        x = randn(rng, T)
        base = one(T) + c * x
        base <= zero(T) && continue
        v = base^3
        u = rand(rng, T)
        if u < one(T) - T(0.0331) * x^4 ||
            log(u) < T(0.5) * x^2 + d * (one(T) - v + log(v))
            return d * v
        end
    end
end

function _apply_stochastic_em_register!(::ScalarCPUStyle,
    model::StochasticMultiplicationRegister, sensor::EMCCDSensor,
    frame::AbstractArray{T}, scratch::AbstractArray, gain,
    rng::AbstractRNG) where {T<:AbstractFloat}
    factor = T(em_register_noise_factor(model, sensor))
    gain_t = T(gain)
    if factor <= zero(T)
        frame .*= gain_t
        return frame
    end
    factor2 = factor * factor
    @inbounds for i in eachindex(frame)
        charge = max(frame[i], zero(T))
        frame[i] = charge <= zero(T) ? zero(T) :
            gain_t * factor2 * _gamma_unit_scale(rng, charge / factor2)
    end
    return frame
end

function _apply_stochastic_em_register!(::AcceleratorStyle,
    model::StochasticMultiplicationRegister, sensor::EMCCDSensor,
    frame::AbstractArray{T}, scratch::AbstractArray, gain,
    rng::AbstractRNG) where {T<:AbstractFloat}
    randn_backend!(rng, scratch)
    factor = T(em_register_noise_factor(model, sensor))
    gain_t = T(gain)
    zero_t = zero(T)
    @. frame = gain_t * max(frame + factor * sqrt(max(frame, zero_t)) * scratch,
        zero_t)
    return frame
end

apply_detection_output!(sensor::EMCCDSensor, det::Detector,
    rng::AbstractRNG) = apply_emccd_detection_output!(sensor.operating_mode,
    det, rng)

apply_emccd_detection_output!(::LinearEMMode, det::Detector,
    rng::AbstractRNG) = det.state.frame

function apply_emccd_detection_output!(mode::PhotonCountingEMMode,
    det::Detector, rng::AbstractRNG)
    threshold = eltype(det.state.frame)(mode.threshold)
    efficiency = eltype(det.state.frame)(mode.detection_efficiency)
    zero_t = zero(eltype(det.state.frame))
    one_t = one(eltype(det.state.frame))
    if efficiency <= zero_t
        fill!(det.state.frame, zero_t)
    elseif efficiency >= one_t
        @. det.state.frame = ifelse(det.state.frame >= threshold, one_t, zero_t)
    else
        rand_uniform_backend!(rng, det.state.noise_buffer)
        @. det.state.frame = ifelse(det.state.frame >= threshold,
            ifelse(det.state.noise_buffer <= efficiency, one_t, zero_t), zero_t)
    end
    return det.state.frame
end

function _batched_sensor_statistics!(sensor::EMCCDSensor, det::Detector,
    cube::AbstractArray, scratch::AbstractArray, rng::AbstractRNG,
    exposure_time::Real)
    rate = effective_cic_rate(det)
    if rate > zero(rate)
        fill!(scratch, rate)
        poisson_noise_frame!(det, rng, scratch)
        cube .+= scratch
    end
    return cube
end

function _batched_pre_readout_gain!(sensor::EMCCDSensor, det::Detector, cube::AbstractArray,
    scratch::AbstractArray, rng::AbstractRNG)
    return _batched_pre_readout_gain!(sensor.output_path, sensor, det, cube, scratch, rng)
end

function _batched_pre_readout_gain!(::EMOutput, sensor::EMCCDSensor, det::Detector,
    cube::AbstractArray, scratch::AbstractArray, rng::AbstractRNG)
    if is_excess_noise_model(sensor.multiplication_model)
        _batched_avalanche_excess_noise!(sensor.excess_noise_factor, cube, scratch, rng)
        cube .*= det.params.gain
        return cube
    end
    return _apply_stochastic_em_register!(execution_style(cube),
        sensor.multiplication_model, sensor, cube, scratch, det.params.gain, rng)
end

_batched_pre_readout_gain!(::ConventionalOutput, sensor::EMCCDSensor, det::Detector,
    cube::AbstractArray, scratch::AbstractArray, rng::AbstractRNG) = cube

function _batched_post_readout_gain!(sensor::EMCCDSensor, det::Detector, cube::AbstractArray)
    return cube
end

_batched_detection_output!(sensor::EMCCDSensor, det::Detector,
    cube::AbstractArray, scratch::AbstractArray, rng::AbstractRNG) =
    _batched_emccd_detection_output!(sensor.operating_mode, cube, scratch, rng)

_batched_emccd_detection_output!(::LinearEMMode, cube::AbstractArray,
    scratch::AbstractArray, rng::AbstractRNG) = cube

function _batched_emccd_detection_output!(mode::PhotonCountingEMMode,
    cube::AbstractArray, scratch::AbstractArray, rng::AbstractRNG)
    threshold = eltype(cube)(mode.threshold)
    efficiency = eltype(cube)(mode.detection_efficiency)
    zero_t = zero(eltype(cube))
    one_t = one(eltype(cube))
    if efficiency <= zero_t
        fill!(cube, zero_t)
    elseif efficiency >= one_t
        @. cube = ifelse(cube >= threshold, one_t, zero_t)
    else
        rand_uniform_backend!(rng, scratch)
        @. cube = ifelse(cube >= threshold,
            ifelse(scratch <= efficiency, one_t, zero_t), zero_t)
    end
    return cube
end

function sensor_saturation_limit(sensor::EMCCDSensor, det::Detector)
    return sensor_saturation_limit(sensor.output_path, sensor, det)
end

function sensor_saturation_limit(::EMOutput, sensor::EMCCDSensor, det::Detector)
    sensor.register_full_well === nothing && return det.params.full_well
    return sensor.register_full_well / det.params.gain
end

sensor_saturation_limit(::ConventionalOutput, sensor::EMCCDSensor, det::Detector) = det.params.full_well

function emccd_snr(signal_electrons::Real; dark_electrons::Real=0, cic_electrons::Real=0,
    readout_noise::Real=0, gain::Real=1, excess_noise_factor::Real=sqrt(2),
    operating_mode::AbstractEMCCDOperatingMode=LinearEMMode(),
    output_path::AbstractEMCCDOutputPath=EMOutput(), T::Type{<:AbstractFloat}=Float64)
    signal = T(signal_electrons)
    dark = T(dark_electrons)
    cic = T(cic_electrons)
    read_noise = T(readout_noise)
    em_gain = T(gain)
    enf = T(excess_noise_factor)
    signal >= zero(T) || throw(InvalidConfiguration("emccd_snr signal_electrons must be >= 0"))
    dark >= zero(T) || throw(InvalidConfiguration("emccd_snr dark_electrons must be >= 0"))
    cic >= zero(T) || throw(InvalidConfiguration("emccd_snr cic_electrons must be >= 0"))
    read_noise >= zero(T) || throw(InvalidConfiguration("emccd_snr readout_noise must be >= 0"))
    em_gain > zero(T) || throw(InvalidConfiguration("emccd_snr gain must be > 0"))
    enf >= one(T) || throw(InvalidConfiguration("emccd_snr excess_noise_factor must be >= 1"))
    mode = convert_emccd_operating_mode(operating_mode, T)
    return _emccd_snr(mode, output_path, signal, dark, cic, read_noise, em_gain, enf)
end

function _emccd_snr(::LinearEMMode, ::EMOutput, signal, dark, cic, readout_noise, gain, excess_noise_factor)
    variance = excess_noise_factor^2 * (signal + dark + cic) + (readout_noise / gain)^2
    return _snr_ratio(signal, variance)
end

function _emccd_snr(::LinearEMMode, ::ConventionalOutput, signal, dark, cic, readout_noise, gain, excess_noise_factor)
    variance = signal + dark + cic + readout_noise^2
    return _snr_ratio(signal, variance)
end

function _emccd_snr(mode::PhotonCountingEMMode, ::AbstractEMCCDOutputPath, signal, dark, cic, readout_noise, gain, excess_noise_factor)
    detected_signal = mode.detection_efficiency * signal
    detected_background = mode.detection_efficiency * (dark + cic)
    return _snr_ratio(detected_signal, detected_signal + detected_background)
end

function _snr_ratio(signal, variance)
    variance > zero(variance) && return signal / sqrt(variance)
    return signal > zero(signal) ? oftype(signal / one(signal), Inf) : zero(signal)
end
