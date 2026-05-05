abstract type AbstractEMGainModel end
abstract type AbstractEMCCDOperatingMode end
abstract type AbstractEMCCDOutputPath end

struct ExcessNoiseApproximation <: AbstractEMGainModel end
struct LinearEMMode <: AbstractEMCCDOperatingMode end
struct EMOutput <: AbstractEMCCDOutputPath end
struct ConventionalOutput <: AbstractEMCCDOutputPath end

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

struct EMCCDSensor{T<:AbstractFloat,M<:AbstractEMGainModel,O<:AbstractEMCCDOperatingMode,P<:AbstractEMCCDOutputPath} <: AvalancheFrameSensorType
    excess_noise_factor::T
    cic_rate::T
    multiplication_model::M
    register_full_well::Union{Nothing,T}
    operating_mode::O
    output_path::P
    em_gain_range::Tuple{T,T}
    readout_rate_hz::Union{Nothing,T}
end

function EMCCDSensor(; excess_noise_factor::Real=1.0, cic_rate::Real=0.0,
    multiplication_model::AbstractEMGainModel=ExcessNoiseApproximation(),
    register_full_well::Union{Nothing,Real}=nothing,
    operating_mode::AbstractEMCCDOperatingMode=LinearEMMode(),
    output_path::AbstractEMCCDOutputPath=EMOutput(),
    em_gain_range::Tuple{<:Real,<:Real}=(1.0, Inf),
    readout_rate_hz::Union{Nothing,Real}=nothing, T::Type{<:AbstractFloat}=Float64)
    excess_noise_factor >= 1 || throw(InvalidConfiguration("EMCCDSensor excess_noise_factor must be >= 1"))
    cic_rate >= 0 || throw(InvalidConfiguration("EMCCDSensor cic_rate must be >= 0"))
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
    return EMCCDSensor{T,typeof(validated_model),typeof(validated_mode),typeof(output_path)}(
        T(excess_noise_factor), T(cic_rate), validated_model, reg_full_well, validated_mode,
        output_path, em_range, readout_rate)
end

detector_sensor_symbol(::EMCCDSensor) = :emccd
supports_clock_induced_charge(::EMCCDSensor) = true
supports_photon_number_resolving(sensor::EMCCDSensor) = supports_emccd_photon_counting(sensor.operating_mode)
supports_emccd_photon_counting(::AbstractEMCCDOperatingMode) = false
supports_emccd_photon_counting(::PhotonCountingEMMode) = true
configured_cic_rate(sensor::EMCCDSensor, ::Type{T}) where {T<:AbstractFloat} = T(sensor.cic_rate)
is_excess_noise_model(::AbstractEMGainModel) = false
is_excess_noise_model(::ExcessNoiseApproximation) = true

convert_em_gain_model(::ExcessNoiseApproximation, ::Type{T}) where {T<:AbstractFloat} = ExcessNoiseApproximation()
convert_em_gain_model(model::StochasticMultiplicationRegister, ::Type{T}) where {T<:AbstractFloat} =
    StochasticMultiplicationRegister{T}(T(model.register_noise_factor))
convert_emccd_operating_mode(::LinearEMMode, ::Type{T}) where {T<:AbstractFloat} = LinearEMMode()
convert_emccd_operating_mode(mode::PhotonCountingEMMode, ::Type{T}) where {T<:AbstractFloat} =
    PhotonCountingEMMode(threshold=T(mode.threshold), detection_efficiency=T(mode.detection_efficiency), T=T)
validate_em_gain_model(::ExcessNoiseApproximation) = ExcessNoiseApproximation()
validate_emccd_operating_mode(::LinearEMMode) = LinearEMMode()
validate_emccd_operating_mode(mode::PhotonCountingEMMode) = mode

function validate_em_gain_model(model::StochasticMultiplicationRegister)
    model.register_noise_factor >= zero(model.register_noise_factor) ||
        throw(InvalidConfiguration("StochasticMultiplicationRegister register_noise_factor must be >= 0"))
    return model
end

function apply_sensor_statistics!(sensor::EMCCDSensor, det::Detector, rng::AbstractRNG)
    rate = effective_cic_rate(det) * det.params.integration_time
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
    det.state.frame .*= det.params.gain
    apply_avalanche_excess_noise!(sensor.excess_noise_factor, det, rng)
    return det.state.frame
end

function apply_em_register_gain!(model::StochasticMultiplicationRegister, sensor::EMCCDSensor, det::Detector, rng::AbstractRNG)
    det.state.frame .*= det.params.gain
    randn_frame_noise!(det, rng, det.state.noise_buffer)
    factor = sensor.excess_noise_factor <= one(sensor.excess_noise_factor) ?
        model.register_noise_factor :
        max(model.register_noise_factor, sqrt(sensor.excess_noise_factor^2 - one(sensor.excess_noise_factor)))
    zero_t = zero(eltype(det.state.frame))
    @. det.state.frame = max(det.state.frame + factor * sqrt(max(det.state.frame, zero_t)) * det.state.noise_buffer, zero_t)
    return det.state.frame
end

function apply_post_readout_gain!(sensor::EMCCDSensor, det::Detector)
    return apply_emccd_operating_mode_output!(sensor.operating_mode, sensor, det)
end

apply_emccd_operating_mode_output!(::LinearEMMode, sensor::EMCCDSensor, det::Detector) = det.state.frame

function apply_emccd_operating_mode_output!(mode::PhotonCountingEMMode, sensor::EMCCDSensor, det::Detector)
    threshold = eltype(det.state.frame)(mode.threshold)
    detected = eltype(det.state.frame)(mode.detection_efficiency)
    @. det.state.frame = ifelse(det.state.frame >= threshold, detected, zero(det.state.frame))
    return det.state.frame
end

function _batched_sensor_statistics!(sensor::EMCCDSensor, det::Detector, cube::AbstractArray, scratch::AbstractArray, rng::AbstractRNG)
    rate = effective_cic_rate(det) * det.params.integration_time
    if rate > zero(rate)
        fill!(scratch, rate)
        poisson_noise_frame!(det, rng, scratch)
        cube .+= scratch
    end
    return cube
end

function _batched_pre_readout_gain!(sensor::EMCCDSensor, det::Detector, cube::AbstractArray, rng::AbstractRNG)
    return _batched_pre_readout_gain!(sensor.output_path, sensor, det, cube, rng)
end

function _batched_pre_readout_gain!(::EMOutput, sensor::EMCCDSensor, det::Detector, cube::AbstractArray, rng::AbstractRNG)
    cube .*= det.params.gain
    scratch = similar(cube)
    if is_excess_noise_model(sensor.multiplication_model)
        _batched_avalanche_excess_noise!(sensor.excess_noise_factor, cube, scratch, rng)
        return cube
    end
    randn_backend!(rng, scratch)
    factor = sensor.excess_noise_factor <= one(sensor.excess_noise_factor) ?
        sensor.multiplication_model.register_noise_factor :
        max(sensor.multiplication_model.register_noise_factor, sqrt(sensor.excess_noise_factor^2 - one(sensor.excess_noise_factor)))
    zero_t = zero(eltype(cube))
    @. cube = max(cube + factor * sqrt(max(cube, zero_t)) * scratch, zero_t)
    return cube
end

_batched_pre_readout_gain!(::ConventionalOutput, sensor::EMCCDSensor, det::Detector, cube::AbstractArray, rng::AbstractRNG) = cube

function _batched_post_readout_gain!(sensor::EMCCDSensor, det::Detector, cube::AbstractArray)
    return _batched_emccd_operating_mode_output!(sensor.operating_mode, sensor, det, cube)
end

_batched_emccd_operating_mode_output!(::LinearEMMode, sensor::EMCCDSensor, det::Detector, cube::AbstractArray) = cube

function _batched_emccd_operating_mode_output!(mode::PhotonCountingEMMode, sensor::EMCCDSensor, det::Detector, cube::AbstractArray)
    threshold = eltype(cube)(mode.threshold)
    detected = eltype(cube)(mode.detection_efficiency)
    @. cube = ifelse(cube >= threshold, detected, zero(cube))
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
