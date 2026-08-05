"""
    MKIDArrayCharacteristics(; energy_resolving_power=nothing,
        photon_arrival_time_resolution_s=nothing,
        wavelength_passband_m=nothing, ...)

Optional physical characteristics of an MKID array. Resolving power is the
dimensionless ratio `E/ΔE`; photon-arrival-time resolution is in seconds; and
the wavelength passband, when set, is an inclusive `(minimum, maximum)` tuple
in meters. These characteristics do not add energy estimates or photon-arrival
timestamps to the accumulated-count image produced by `MKIDArrayDetector`.
"""
struct MKIDArrayCharacteristics{T<:AbstractFloat}
    energy_resolving_power::Union{Nothing,T}
    photon_arrival_time_resolution_s::Union{Nothing,T}
    wavelength_passband_m::Union{Nothing,Tuple{T,T}}
end

function MKIDArrayCharacteristics(;
    energy_resolving_power::Union{Nothing,Real}=nothing,
    photon_arrival_time_resolution_s::Union{Nothing,Real}=nothing,
    wavelength_passband_m::Union{Nothing,Tuple{<:Real,<:Real}}=nothing,
    T::Type{<:AbstractFloat}=Float64)
    resolving_power = if energy_resolving_power === nothing
        nothing
    else
        value = T(energy_resolving_power)
        isfinite(value) && value > zero(T) || throw(InvalidConfiguration(
            "MKIDArrayCharacteristics energy_resolving_power must be finite and > 0"))
        value
    end
    arrival_time_resolution = if photon_arrival_time_resolution_s === nothing
        nothing
    else
        value = T(photon_arrival_time_resolution_s)
        isfinite(value) && value > zero(T) || throw(InvalidConfiguration(
            "MKIDArrayCharacteristics photon_arrival_time_resolution_s must be finite and > 0"))
        value
    end
    passband = if wavelength_passband_m === nothing
        nothing
    else
        lo, hi = T.(wavelength_passband_m)
        isfinite(lo) && lo > zero(T) || throw(InvalidConfiguration(
            "MKIDArrayCharacteristics wavelength_passband_m lower bound must be finite and > 0"))
        isfinite(hi) && hi > lo || throw(InvalidConfiguration(
            "MKIDArrayCharacteristics wavelength_passband_m upper bound must be finite and greater than the lower bound"))
        (lo, hi)
    end
    return MKIDArrayCharacteristics{T}(
        resolving_power, arrival_time_resolution, passband)
end

function convert_mkid_characteristics(
    characteristics::MKIDArrayCharacteristics{T}, ::Type{T}) where {
    T<:AbstractFloat}
    return characteristics
end

function convert_mkid_characteristics(
    characteristics::MKIDArrayCharacteristics, ::Type{T}) where {
    T<:AbstractFloat}
    return MKIDArrayCharacteristics(
        energy_resolving_power=characteristics.energy_resolving_power,
        photon_arrival_time_resolution_s=
            characteristics.photon_arrival_time_resolution_s,
        wavelength_passband_m=characteristics.wavelength_passband_m,
        T=T,
    )
end

"""
    MKIDArraySensor(; qe=0.7, dark_count_rate=0, fill_factor=1,
        characteristics=MKIDArrayCharacteristics(), ...)

Physical parameters for an accumulated-count microwave kinetic inductance
detector (MKID). Source-aware capture applies a configured characteristics
passband. Matrix-only capture assumes its input was spectrally prefiltered.
"""
struct MKIDArraySensor{
    T<:AbstractFloat,
    D<:CountingDeadTimeModel,
} <: AbstractMKIDArraySensor
    qe::T
    dark_count_rate::T
    fill_factor::T
    characteristics::MKIDArrayCharacteristics{T}
    dead_time_model::D
end

function MKIDArraySensor(; qe::Real=0.7, dark_count_rate::Real=0.0,
    fill_factor::Real=1.0,
    characteristics::MKIDArrayCharacteristics=MKIDArrayCharacteristics(),
    dead_time_model::CountingDeadTimeModel=NoDeadTime(),
    T::Type{<:AbstractFloat}=Float64)
    typed_qe = T(qe)
    typed_dark_count_rate = T(dark_count_rate)
    typed_fill_factor = T(fill_factor)
    isfinite(typed_qe) && zero(T) <= typed_qe <= one(T) ||
        throw(InvalidConfiguration(
            "MKIDArraySensor qe must be finite and lie in [0, 1]"))
    isfinite(typed_dark_count_rate) && typed_dark_count_rate >= zero(T) ||
        throw(InvalidConfiguration(
            "MKIDArraySensor dark_count_rate must be finite and >= 0"))
    isfinite(typed_fill_factor) && zero(T) < typed_fill_factor <= one(T) ||
        throw(InvalidConfiguration(
            "MKIDArraySensor fill_factor must be finite and lie in (0, 1]"))
    typed_characteristics = convert_mkid_characteristics(characteristics, T)
    dead_time = validate_dead_time_model(
        convert_dead_time_model(dead_time_model, T))
    return MKIDArraySensor{T,typeof(dead_time)}(
        typed_qe,
        typed_dark_count_rate,
        typed_fill_factor,
        typed_characteristics,
        dead_time,
    )
end

struct MKIDArrayExportMetadata{T<:AbstractFloat}
    counting::CountingDetectorExportMetadata{T}
    characteristics::MKIDArrayCharacteristics{T}
    observable::Symbol
    provides_energy_estimates::Bool
    provides_photon_arrival_timestamps::Bool
end

struct MKIDArrayDetectorParams{
    T<:AbstractFloat,
    S<:AbstractMKIDArraySensor,
    G<:AbstractCountingGateModel,
    TM<:AbstractDetectorThermalModel,
}
    exposure_duration::T
    gate_model::G
    thermal_model::TM
    sensor::S
    output_type::Union{Nothing,DataType}
end

struct MKIDArrayDetector{
    N<:NoiseModel,
    P<:MKIDArrayDetectorParams,
    S<:CountingDetectorState,
    W<:CountingDetectorWorkspace,
    R<:CountingDetectorProducts,
    B<:AbstractArrayBackend,
} <: AbstractCountingDetector
    noise::N
    params::P
    state::S
    workspace::W
    products::R
end

@inline backend(
    ::MKIDArrayDetector{<:Any,<:Any,<:Any,<:Any,<:Any,B}) where {B} = B()

counting_sensor(det::MKIDArrayDetector) = det.params.sensor
counting_gate_model(det::MKIDArrayDetector) = det.params.gate_model
counting_dead_time_model(det::MKIDArrayDetector) =
    det.params.sensor.dead_time_model
counting_mean_response_model(::MKIDArrayDetector) = NullCountingMeanResponse()
counting_exposure_duration(det::MKIDArrayDetector) = det.params.exposure_duration
counting_layout(::MKIDArrayDetector) = :pixel_counts
counting_output_type(det::MKIDArrayDetector) = det.params.output_type
counting_array(det::MKIDArrayDetector) = det.products.counts
counting_noise_buffer(det::MKIDArrayDetector) = det.workspace.noise_buffer
counting_host_buffer(det::MKIDArrayDetector) = det.workspace.host_buffer
counting_output_buffer(det::MKIDArrayDetector) = det.products.output_buffer
counting_output_host_buffer(det::MKIDArrayDetector) =
    det.workspace.output_buffer_host
set_counting_array!(det::MKIDArrayDetector, values) =
    (det.products.counts = values; det)
set_counting_noise_buffer!(det::MKIDArrayDetector, values) =
    (det.workspace.noise_buffer = values; det)
set_counting_host_buffer!(det::MKIDArrayDetector, values) =
    (det.workspace.host_buffer = values; det)
set_counting_output_buffer!(det::MKIDArrayDetector, values) =
    (det.products.output_buffer = values; det)
set_counting_output_host_buffer!(det::MKIDArrayDetector, values) =
    (det.workspace.output_buffer_host = values; det)
counting_detection_efficiency(det::MKIDArrayDetector,
    ::Type{T}=eltype(counting_array(det))) where {T<:AbstractFloat} =
    T(det.params.sensor.qe)
counting_fill_factor(det::MKIDArrayDetector,
    ::Type{T}=eltype(counting_array(det))) where {T<:AbstractFloat} =
    T(det.params.sensor.fill_factor)
counting_reported_fill_factor(det::MKIDArrayDetector,
    ::Type{T}=eltype(counting_array(det))) where {T<:AbstractFloat} =
    T(det.params.sensor.fill_factor)
counting_dark_count_rate(det::MKIDArrayDetector,
    ::Type{T}=eltype(counting_array(det))) where {T<:AbstractFloat} =
    T(det.params.sensor.dark_count_rate)

@inline function mkid_wavelength_in_passband(det::MKIDArrayDetector,
    wavelength_m, ::Type{T}) where {T<:AbstractFloat}
    passband = det.params.sensor.characteristics.wavelength_passband_m
    passband === nothing && return true
    lo, hi = passband
    λ = T(wavelength_m)
    return lo <= λ <= hi
end

function counting_source_throughput(det::MKIDArrayDetector,
    src::AbstractSource,
    ::Type{T}=eltype(counting_array(det))) where {T<:AbstractFloat}
    return mkid_wavelength_in_passband(det, wavelength(src), T) ?
        one(T) : zero(T)
end

function counting_source_throughput(det::MKIDArrayDetector,
    src::SpectralSource,
    ::Type{T}=eltype(counting_array(det))) where {T<:AbstractFloat}
    det.params.sensor.characteristics.wavelength_passband_m === nothing &&
        return one(T)
    throughput = zero(T)
    @inbounds for sample in spectral_bundle(src)
        mkid_wavelength_in_passband(det, sample.wavelength, T) || continue
        throughput += T(sample.weight)
    end
    return throughput
end

detector_sensor_symbol(::MKIDArraySensor) = :mkid_array
supports_photon_counting(::MKIDArraySensor) = true

function detector_export_metadata(det::MKIDArrayDetector;
    T::Type{<:AbstractFloat}=eltype(counting_array(det)))
    characteristics = convert_mkid_characteristics(
        det.params.sensor.characteristics, T)
    return MKIDArrayExportMetadata{T}(
        counting_detector_export_metadata(det, T),
        characteristics,
        :accumulated_count_image,
        false,
        false,
    )
end

dark_count_law(::MKIDArrayDetector) = NullTemperatureLaw()

function _build_mkid_array_detector(noise::NoiseModel;
    exposure_duration::Real,
    gate_model::AbstractCountingGateModel,
    thermal_model::AbstractDetectorThermalModel,
    sensor::AbstractMKIDArraySensor,
    output_type::Union{Nothing,DataType},
    T::Type{<:AbstractFloat},
    backend)
    typed_exposure_duration = T(exposure_duration)
    isfinite(typed_exposure_duration) && typed_exposure_duration > zero(T) ||
        throw(InvalidConfiguration(
            "MKIDArrayDetector exposure_duration must be finite and > 0"))
    converted = convert_noise(noise, T)
    validated = validate_counting_noise(converted)
    gate = validate_gate_model(convert_gate_model(gate_model, T))
    thermal = validate_thermal_model(convert_thermal_model(thermal_model, T))
    typed_sensor = convert_mkid_sensor(sensor, T)
    params = MKIDArrayDetectorParams{
        T,typeof(typed_sensor),typeof(gate),typeof(thermal)}(
        typed_exposure_duration,
        gate,
        thermal,
        typed_sensor,
        output_type,
    )
    counts = backend{T}(undef, 1, 1)
    noise_buffer = backend{T}(undef, 1, 1)
    host_buffer = Matrix{T}(undef, 1, 1)
    output_buffer = output_type === nothing ? nothing :
        backend{output_type}(undef, 1, 1)
    output_buffer_host = output_type === nothing ? nothing :
        Matrix{output_type}(undef, 1, 1)
    fill!(counts, zero(T))
    fill!(noise_buffer, zero(T))
    fill!(host_buffer, zero(T))
    output_buffer === nothing ||
        fill!(output_buffer, zero(eltype(output_buffer)))
    output_buffer_host === nothing || fill!(output_buffer_host,
        zero(eltype(output_buffer_host)))
    thermal_state = thermal_state_from_model(thermal, T)
    state = CountingDetectorState{typeof(thermal_state)}(thermal_state)
    workspace = CountingDetectorWorkspace{T,typeof(noise_buffer),
        typeof(host_buffer),typeof(output_buffer_host)}(
        noise_buffer, host_buffer, output_buffer_host)
    products = CountingDetectorProducts{T,typeof(counts),
        typeof(output_buffer)}(counts, output_buffer)
    selector = _resolve_backend_selector(backend)
    return MKIDArrayDetector{
        typeof(validated),typeof(params),typeof(state),typeof(workspace),
        typeof(products),typeof(selector)}(
        validated, params, state, workspace, products)
end

"""
    MKIDArrayDetector(; sensor=MKIDArraySensor(), ...)

Accumulated-count image detector backed by `MKIDArraySensor`. The detector uses
the shared counting pipeline and preallocated buffers. It exports an image, not
a per-photon event stream or a set of per-photon energy estimates.
"""
function MKIDArrayDetector(; exposure_duration::Real=1.0,
    noise::NoiseModel=NoisePhoton(),
    sensor::MKIDArraySensor=MKIDArraySensor(),
    output_type::Union{Nothing,DataType}=nothing,
    gate_model::AbstractCountingGateModel=NullCountingGate(),
    thermal_model::AbstractDetectorThermalModel=NullDetectorThermalModel(),
    T::Type{<:AbstractFloat}=Float64,
    backend::AbstractArrayBackend=CPUBackend())
    backend = _resolve_array_backend(backend)
    return _build_mkid_array_detector(noise;
        exposure_duration=exposure_duration,
        gate_model=gate_model,
        thermal_model=thermal_model,
        sensor=sensor,
        output_type=output_type,
        T=T,
        backend=backend)
end

function convert_mkid_sensor(sensor::MKIDArraySensor{T},
    ::Type{T}) where {T<:AbstractFloat}
    return sensor
end

function convert_mkid_sensor(sensor::AbstractMKIDArraySensor,
    ::Type{T}) where {T<:AbstractFloat}
    return MKIDArraySensor(
        qe=sensor.qe,
        dark_count_rate=sensor.dark_count_rate,
        fill_factor=sensor.fill_factor,
        characteristics=sensor.characteristics,
        dead_time_model=sensor.dead_time_model,
        T=T,
    )
end
