abstract type AbstractEMGainModel end

struct ExcessNoiseApproximation <: AbstractEMGainModel end

struct StochasticMultiplicationRegister{T<:AbstractFloat} <: AbstractEMGainModel
    register_noise_factor::T
end

StochasticMultiplicationRegister(register_noise_factor::Real) = StochasticMultiplicationRegister{Float64}(float(register_noise_factor))

struct EMCCDSensor{T<:AbstractFloat,M<:AbstractEMGainModel} <: AvalancheFrameSensorType
    excess_noise_factor::T
    cic_rate::T
    multiplication_model::M
    register_full_well::Union{Nothing,T}
end

function EMCCDSensor(; excess_noise_factor::Real=1.0, cic_rate::Real=0.0,
    multiplication_model::AbstractEMGainModel=ExcessNoiseApproximation(),
    register_full_well::Union{Nothing,Real}=nothing, T::Type{<:AbstractFloat}=Float64)
    excess_noise_factor >= 1 || throw(InvalidConfiguration("EMCCDSensor excess_noise_factor must be >= 1"))
    cic_rate >= 0 || throw(InvalidConfiguration("EMCCDSensor cic_rate must be >= 0"))
    reg_full_well = register_full_well === nothing ? nothing : T(register_full_well)
    reg_full_well === nothing || reg_full_well > zero(reg_full_well) ||
        throw(InvalidConfiguration("EMCCDSensor register_full_well must be > 0"))
    converted_model = convert_em_gain_model(multiplication_model, T)
    validated_model = validate_em_gain_model(converted_model)
    return EMCCDSensor{T,typeof(validated_model)}(T(excess_noise_factor), T(cic_rate), validated_model, reg_full_well)
end

detector_sensor_symbol(::EMCCDSensor) = :emccd
supports_clock_induced_charge(::EMCCDSensor) = true
configured_cic_rate(sensor::EMCCDSensor, ::Type{T}) where {T<:AbstractFloat} = T(sensor.cic_rate)
is_excess_noise_model(::AbstractEMGainModel) = false
is_excess_noise_model(::ExcessNoiseApproximation) = true

convert_em_gain_model(::ExcessNoiseApproximation, ::Type{T}) where {T<:AbstractFloat} = ExcessNoiseApproximation()
convert_em_gain_model(model::StochasticMultiplicationRegister, ::Type{T}) where {T<:AbstractFloat} =
    StochasticMultiplicationRegister{T}(T(model.register_noise_factor))
validate_em_gain_model(::ExcessNoiseApproximation) = ExcessNoiseApproximation()

function validate_em_gain_model(model::StochasticMultiplicationRegister)
    model.register_noise_factor >= zero(model.register_noise_factor) ||
        throw(InvalidConfiguration("StochasticMultiplicationRegister register_noise_factor must be >= 0"))
    return model
end

function apply_sensor_statistics!(sensor::EMCCDSensor, det::Detector, rng::AbstractRNG)
    rate = effective_cic_rate(det) * det.params.integration_time
    if rate > zero(rate)
        fill!(det.state.noise_buffer, rate)
        poisson_noise!(rng, det.state.noise_buffer)
        det.state.frame .+= det.state.noise_buffer
    end
    return det.state.frame
end

function apply_pre_readout_gain!(sensor::EMCCDSensor, det::Detector, rng::AbstractRNG)
    return apply_pre_readout_gain!(sensor.multiplication_model, sensor, det, rng)
end

function apply_pre_readout_gain!(::ExcessNoiseApproximation, sensor::EMCCDSensor, det::Detector, rng::AbstractRNG)
    det.state.frame .*= det.params.gain
    apply_avalanche_excess_noise!(sensor.excess_noise_factor, det, rng)
    return det.state.frame
end

function apply_pre_readout_gain!(model::StochasticMultiplicationRegister, sensor::EMCCDSensor, det::Detector, rng::AbstractRNG)
    det.state.frame .*= det.params.gain
    randn_backend!(rng, det.state.noise_buffer)
    factor = sensor.excess_noise_factor <= one(sensor.excess_noise_factor) ?
        model.register_noise_factor :
        max(model.register_noise_factor, sqrt(sensor.excess_noise_factor^2 - one(sensor.excess_noise_factor)))
    zero_t = zero(eltype(det.state.frame))
    @. det.state.frame = max(det.state.frame + factor * sqrt(max(det.state.frame, zero_t)) * det.state.noise_buffer, zero_t)
    return det.state.frame
end

function _batched_sensor_statistics!(sensor::EMCCDSensor, det::Detector, cube::AbstractArray, scratch::AbstractArray, rng::AbstractRNG)
    rate = effective_cic_rate(det) * det.params.integration_time
    if rate > zero(rate)
        fill!(scratch, rate)
        poisson_noise!(rng, scratch)
        cube .+= scratch
    end
    return cube
end

function _batched_pre_readout_gain!(sensor::EMCCDSensor, det::Detector, cube::AbstractArray, rng::AbstractRNG)
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

function sensor_saturation_limit(sensor::EMCCDSensor, det::Detector)
    sensor.register_full_well === nothing && return det.params.full_well
    return sensor.register_full_well / det.params.gain
end
