"""
    MKIDArraySensor(; qe=0.7, dark_count_rate=0, fill_factor=1,
        energy_resolution=10, timing_jitter_s=1e-6,
        wavelength_range_m=nothing, ...)

Physical parameters for an accumulated-image microwave kinetic inductance
detector (MKID). `energy_resolution` is the dimensionless resolving power
`E/ΔE`; `timing_jitter_s` is in seconds; and `wavelength_range_m`, when set, is
an inclusive `(minimum, maximum)` passband in meters.

The passband is applied by source-aware `capture!` calls. Matrix-only capture
assumes that its input has already been spectrally filtered. Energy resolution
and timing jitter are exported as metadata; this model does not generate
per-photon energy/timestamp events.
"""
struct MKIDArraySensor{
    T<:AbstractFloat,
    D<:CountingDeadTimeModel,
    C<:AbstractCountingCorrelationModel,
} <: MKIDArraySensorType
    qe::T
    dark_count_rate::T
    fill_factor::T
    energy_resolution::T
    timing_jitter_s::T
    wavelength_min_m::Union{Nothing,T}
    wavelength_max_m::Union{Nothing,T}
    dead_time_model::D
    correlation_model::C
end

function MKIDArraySensor(; qe::Real=0.7, dark_count_rate::Real=0.0, fill_factor::Real=1.0,
    energy_resolution::Real=10.0, timing_jitter_s::Real=1e-6,
    wavelength_range_m::Union{Nothing,Tuple{<:Real,<:Real}}=nothing,
    dead_time_model::CountingDeadTimeModel=NoDeadTime(),
    correlation_model::AbstractCountingCorrelationModel=NullCountingCorrelation(),
    T::Type{<:AbstractFloat}=Float64)
    typed_qe = T(qe)
    typed_dark_count_rate = T(dark_count_rate)
    typed_fill_factor = T(fill_factor)
    typed_energy_resolution = T(energy_resolution)
    typed_timing_jitter_s = T(timing_jitter_s)
    isfinite(typed_qe) && zero(T) <= typed_qe <= one(T) ||
        throw(InvalidConfiguration("MKIDArraySensor qe must lie in [0, 1]"))
    isfinite(typed_dark_count_rate) && typed_dark_count_rate >= zero(T) ||
        throw(InvalidConfiguration("MKIDArraySensor dark_count_rate must be >= 0"))
    isfinite(typed_fill_factor) && zero(T) < typed_fill_factor <= one(T) ||
        throw(InvalidConfiguration("MKIDArraySensor fill_factor must lie in (0, 1]"))
    isfinite(typed_energy_resolution) && typed_energy_resolution > zero(T) ||
        throw(InvalidConfiguration("MKIDArraySensor energy_resolution must be > 0"))
    isfinite(typed_timing_jitter_s) && typed_timing_jitter_s >= zero(T) ||
        throw(InvalidConfiguration("MKIDArraySensor timing_jitter_s must be >= 0"))
    wavelength_min_m, wavelength_max_m = if wavelength_range_m === nothing
        nothing, nothing
    else
        lo, hi = T.(wavelength_range_m)
        isfinite(lo) && lo > zero(T) ||
            throw(InvalidConfiguration("MKIDArraySensor wavelength_range_m lower bound must be finite and > 0"))
        isfinite(hi) && hi > lo ||
            throw(InvalidConfiguration("MKIDArraySensor wavelength_range_m upper bound must be finite and > lower bound"))
        lo, hi
    end
    dead_time = validate_dead_time_model(convert_dead_time_model(dead_time_model, T))
    correlation = validate_correlation_model(convert_correlation_model(correlation_model, T))
    return MKIDArraySensor{T,typeof(dead_time),typeof(correlation)}(
        typed_qe,
        typed_dark_count_rate,
        typed_fill_factor,
        typed_energy_resolution,
        typed_timing_jitter_s,
        wavelength_min_m,
        wavelength_max_m,
        dead_time,
        correlation,
    )
end

struct MKIDArrayDetectorParams{
    T<:AbstractFloat,
    S<:MKIDArraySensorType,
    G<:AbstractCountingGateModel,
    TM<:AbstractDetectorThermalModel,
}
    integration_time::T
    gate_model::G
    thermal_model::TM
    sensor::S
    output_type::Union{Nothing,DataType}
    layout::Symbol
end

mutable struct MKIDArrayDetectorState{
    T<:AbstractFloat,
    A<:AbstractMatrix{T},
    H<:AbstractMatrix{T},
    O,
    OH,
    TS<:AbstractDetectorThermalState,
}
    counts::A
    noise_buffer::A
    host_buffer::H
    output_buffer::O
    output_buffer_host::OH
    thermal_state::TS
end

struct MKIDArrayDetector{
    N<:NoiseModel,
    P<:MKIDArrayDetectorParams,
    S<:MKIDArrayDetectorState,
    B<:AbstractArrayBackend,
} <: AbstractCountingDetector
    noise::N
    params::P
    state::S
end

@inline backend(::MKIDArrayDetector{<:Any,<:Any,<:Any,B}) where {B} = B()

counting_sensor(det::MKIDArrayDetector) = det.params.sensor
counting_gate_model(det::MKIDArrayDetector) = det.params.gate_model
counting_dead_time_model(det::MKIDArrayDetector) = det.params.sensor.dead_time_model
counting_correlation_model(det::MKIDArrayDetector) = det.params.sensor.correlation_model
counting_integration_time(det::MKIDArrayDetector) = det.params.integration_time
counting_layout(det::MKIDArrayDetector) = det.params.layout
counting_output_type(det::MKIDArrayDetector) = det.params.output_type
counting_array(det::MKIDArrayDetector) = det.state.counts
counting_noise_buffer(det::MKIDArrayDetector) = det.state.noise_buffer
counting_host_buffer(det::MKIDArrayDetector) = det.state.host_buffer
counting_output_buffer(det::MKIDArrayDetector) = det.state.output_buffer
counting_output_host_buffer(det::MKIDArrayDetector) = det.state.output_buffer_host
set_counting_array!(det::MKIDArrayDetector, values) = (det.state.counts = values; det)
set_counting_noise_buffer!(det::MKIDArrayDetector, values) = (det.state.noise_buffer = values; det)
set_counting_host_buffer!(det::MKIDArrayDetector, values) = (det.state.host_buffer = values; det)
set_counting_output_buffer!(det::MKIDArrayDetector, values) = (det.state.output_buffer = values; det)
set_counting_output_host_buffer!(det::MKIDArrayDetector, values) =
    (det.state.output_buffer_host = values; det)
counting_qe(det::MKIDArrayDetector, ::Type{T}=eltype(counting_array(det))) where {T<:AbstractFloat} = T(det.params.sensor.qe)
counting_fill_factor(det::MKIDArrayDetector, ::Type{T}=eltype(counting_array(det))) where {T<:AbstractFloat} = T(det.params.sensor.fill_factor)
counting_reported_fill_factor(det::MKIDArrayDetector, ::Type{T}=eltype(counting_array(det))) where {T<:AbstractFloat} = T(det.params.sensor.fill_factor)
counting_dark_count_rate(det::MKIDArrayDetector, ::Type{T}=eltype(counting_array(det))) where {T<:AbstractFloat} = T(det.params.sensor.dark_count_rate)
counting_energy_resolution(det::MKIDArrayDetector, ::Type{T}=eltype(counting_array(det))) where {T<:AbstractFloat} = T(det.params.sensor.energy_resolution)
counting_timing_jitter_s(det::MKIDArrayDetector, ::Type{T}=eltype(counting_array(det))) where {T<:AbstractFloat} = T(det.params.sensor.timing_jitter_s)
counting_wavelength_min_m(det::MKIDArrayDetector, ::Type{T}=eltype(counting_array(det))) where {T<:AbstractFloat} =
    det.params.sensor.wavelength_min_m === nothing ? nothing : T(det.params.sensor.wavelength_min_m)
counting_wavelength_max_m(det::MKIDArrayDetector, ::Type{T}=eltype(counting_array(det))) where {T<:AbstractFloat} =
    det.params.sensor.wavelength_max_m === nothing ? nothing : T(det.params.sensor.wavelength_max_m)

@inline function mkid_wavelength_in_range(det::MKIDArrayDetector, wavelength_m, ::Type{T}) where {T<:AbstractFloat}
    lo = counting_wavelength_min_m(det, T)
    lo === nothing && return true
    hi = counting_wavelength_max_m(det, T)
    λ = T(wavelength_m)
    return lo <= λ <= hi
end

function counting_source_throughput(det::MKIDArrayDetector, src::AbstractSource,
    ::Type{T}=eltype(counting_array(det))) where {T<:AbstractFloat}
    return mkid_wavelength_in_range(det, wavelength(src), T) ? one(T) : zero(T)
end

function counting_source_throughput(det::MKIDArrayDetector, src::SpectralSource,
    ::Type{T}=eltype(counting_array(det))) where {T<:AbstractFloat}
    counting_wavelength_min_m(det, T) === nothing && return one(T)
    throughput = zero(T)
    @inbounds for sample in spectral_bundle(src)
        mkid_wavelength_in_range(det, sample.wavelength, T) || continue
        throughput += T(sample.weight)
    end
    return throughput
end

detector_sensor_symbol(::MKIDArraySensor) = :mkid_array
supports_photon_number_resolving(::MKIDArraySensor) = true
supports_energy_resolving(::MKIDArraySensor) = true

dark_count_law(::MKIDArrayDetector) = NullTemperatureLaw()

function _build_mkid_array_detector(noise::NoiseModel; integration_time::Real,
    gate_model::AbstractCountingGateModel,
    thermal_model::AbstractDetectorThermalModel,
    sensor::MKIDArraySensorType,
    output_type::Union{Nothing,DataType},
    layout::Symbol,
    T::Type{<:AbstractFloat},
    backend)
    integration_time > 0 || throw(InvalidConfiguration("MKIDArrayDetector integration_time must be > 0"))
    converted = convert_noise(noise, T)
    validated = validate_counting_noise(converted)
    gate = validate_gate_model(convert_gate_model(gate_model, T))
    thermal = validate_thermal_model(convert_thermal_model(thermal_model, T))
    typed_sensor = convert_mkid_sensor(sensor, T)
    params = MKIDArrayDetectorParams{T,typeof(typed_sensor),typeof(gate),typeof(thermal)}(
        T(integration_time),
        gate,
        thermal,
        typed_sensor,
        output_type,
        layout,
    )
    counts = backend{T}(undef, 1, 1)
    noise_buffer = backend{T}(undef, 1, 1)
    host_buffer = Matrix{T}(undef, 1, 1)
    output_buffer = output_type === nothing ? nothing : backend{output_type}(undef, 1, 1)
    output_buffer_host = output_type === nothing ? nothing : Matrix{output_type}(undef, 1, 1)
    fill!(counts, zero(T))
    fill!(noise_buffer, zero(T))
    fill!(host_buffer, zero(T))
    output_buffer === nothing || fill!(output_buffer, zero(eltype(output_buffer)))
    output_buffer_host === nothing || fill!(output_buffer_host,
        zero(eltype(output_buffer_host)))
    thermal_state = thermal_state_from_model(thermal, T)
    state = MKIDArrayDetectorState{T,typeof(counts),typeof(host_buffer),
        typeof(output_buffer),typeof(output_buffer_host),typeof(thermal_state)}(
        counts, noise_buffer, host_buffer, output_buffer, output_buffer_host,
        thermal_state)
    selector = _resolve_backend_selector(backend)
    return MKIDArrayDetector{typeof(validated),typeof(params),typeof(state),typeof(selector)}(
        validated, params, state)
end

"""
    MKIDArrayDetector(; sensor=MKIDArraySensor(), ...)

Accumulated-count image detector backed by `MKIDArraySensor`. The detector uses
the shared counting pipeline and preallocated buffers. It exports a frame, not a
per-photon event stream.
"""
function MKIDArrayDetector(; integration_time::Real=1.0, noise::NoiseModel=NoisePhoton(),
    sensor::MKIDArraySensor=MKIDArraySensor(),
    output_type::Union{Nothing,DataType}=nothing,
    layout::Symbol=:pixel_counts,
    gate_model::AbstractCountingGateModel=NullCountingGate(),
    thermal_model::AbstractDetectorThermalModel=NullDetectorThermalModel(),
    T::Type{<:AbstractFloat}=Float64,
    backend::AbstractArrayBackend=CPUBackend())
    backend = _resolve_array_backend(backend)
    return _build_mkid_array_detector(noise; integration_time=integration_time,
        gate_model=gate_model, thermal_model=thermal_model, sensor=sensor,
        output_type=output_type, layout=layout, T=T, backend=backend)
end

convert_mkid_sensor(sensor::MKIDArraySensor{T}, ::Type{T}) where {T<:AbstractFloat} = sensor

function convert_mkid_sensor(sensor::MKIDArraySensorType, ::Type{T}) where {T<:AbstractFloat}
    wavelength_range_m = sensor.wavelength_min_m === nothing ? nothing :
        (sensor.wavelength_min_m, sensor.wavelength_max_m)
    return MKIDArraySensor(
        qe=sensor.qe,
        dark_count_rate=sensor.dark_count_rate,
        fill_factor=sensor.fill_factor,
        energy_resolution=sensor.energy_resolution,
        timing_jitter_s=sensor.timing_jitter_s,
        wavelength_range_m=wavelength_range_m,
        dead_time_model=sensor.dead_time_model,
        correlation_model=sensor.correlation_model,
        T=T,
    )
end
