struct InGaAsDetectorAcquisitionNode{T<:AbstractFloat} end

"""Immutable configuration for one complete-frame InGaAs acquisition."""
struct InGaAsDetectorAcquisitionNodeConfig{T,N,P,L}
    rows::Int
    columns::Int
    pixel_sampling_m::Tuple{T,T}
    spectral_channel::Symbol
    exposure_duration_s::T
    quantum_efficiency::T
    gain::T
    dark_current_e_per_pixel_s::T
    glow_e_per_pixel_s::T
    bits::Union{Nothing,Int}
    full_well_e::Union{Nothing,T}
    noise::N
    persistence_model::P
    nonlinearity_model::L
    rng_seed::UInt64
    photon_rate_schema::String
    frame_schema::String
end

function _ingaas_detector_acquisition_config(::Type{T};
        rows::Integer, columns::Integer, pixel_sampling_m::Tuple{Real,Real},
        spectral_channel::Symbol, exposure_duration_s::Real,
        quantum_efficiency::Real, gain::Real, dark_current_e_per_pixel_s::Real,
        glow_e_per_pixel_s::Real, bits::Union{Nothing,Integer},
        full_well_e::Union{Nothing,Real}, photon_noise::Bool,
        readout_noise::Bool, readout_noise_e::Real,
        persistence_model::Union{NullPersistence,ExponentialPersistence},
        nonlinearity_model::Union{NullFrameNonlinearity,SaturatingFrameNonlinearity},
        rng_seed::Integer, photon_rate_schema::AbstractString,
        frame_schema::AbstractString) where {T<:AbstractFloat}
    T in (Float32, Float64) || throw(AlgorithmGraphError(
        "InGaAs graph acquisition requires Float32 or Float64"))
    n_rows, n_columns = Int(rows), Int(columns)
    n_rows > 0 && n_columns > 0 || throw(AlgorithmGraphError(
        "InGaAs photon-rate dimensions must be positive"))
    sampling = (T(pixel_sampling_m[1]), T(pixel_sampling_m[2]))
    all(value -> isfinite(value) && value > zero(T), sampling) ||
        throw(AlgorithmGraphError("InGaAs pixel_sampling_m must be finite and positive"))
    isempty(string(spectral_channel)) && throw(AlgorithmGraphError(
        "InGaAs spectral_channel must not be empty"))
    exposure, qe, linear_gain = T(exposure_duration_s), T(quantum_efficiency), T(gain)
    dark, glow = T(dark_current_e_per_pixel_s), T(glow_e_per_pixel_s)
    isfinite(exposure) && exposure > zero(T) || throw(AlgorithmGraphError(
        "InGaAs exposure_duration_s must be finite and positive"))
    isfinite(qe) && zero(T) <= qe <= one(T) || throw(AlgorithmGraphError(
        "InGaAs quantum_efficiency must lie in [0, 1]"))
    isfinite(linear_gain) && linear_gain > zero(T) || throw(AlgorithmGraphError(
        "InGaAs gain must be finite and positive"))
    isfinite(dark) && dark >= zero(T) || throw(AlgorithmGraphError(
        "InGaAs dark_current_e_per_pixel_s must be finite and nonnegative"))
    isfinite(glow) && glow >= zero(T) || throw(AlgorithmGraphError(
        "InGaAs glow_e_per_pixel_s must be finite and nonnegative"))
    output_bits = isnothing(bits) ? nothing : Int(bits)
    isnothing(output_bits) || 1 <= output_bits <= 64 || throw(AlgorithmGraphError(
        "InGaAs bits must lie in 1:64"))
    isnothing(output_bits) || !isnothing(full_well_e) || throw(AlgorithmGraphError(
        "InGaAs bits requires full_well_e"))
    full_well = isnothing(full_well_e) ? nothing : T(full_well_e)
    isnothing(full_well) || isfinite(full_well) && full_well > zero(T) ||
        throw(AlgorithmGraphError("InGaAs full_well_e must be finite and positive"))
    0 <= rng_seed <= typemax(UInt64) || throw(AlgorithmGraphError(
        "InGaAs rng_seed must lie in the UInt64 range"))
    isempty(photon_rate_schema) && throw(AlgorithmGraphError(
        "InGaAs photon_rate_schema must not be empty"))
    isempty(frame_schema) && throw(AlgorithmGraphError(
        "InGaAs frame_schema must not be empty"))
    noise = _frame_detector_noise(T, "InGaAs", photon_noise,
        readout_noise, readout_noise_e)
    # The sensor validates and converts the immutable persistence coefficients;
    # evolving latent charge is created only during detector preparation.
    sensor = InGaAsSensor(; glow_rate=glow, persistence_model, T)
    nonlinearity = _ingaas_graph_nonlinearity(nonlinearity_model, T)
    return InGaAsDetectorAcquisitionNodeConfig(n_rows, n_columns, sampling,
        spectral_channel, exposure, qe, linear_gain, dark, glow, output_bits,
        full_well, noise, sensor.persistence_model, nonlinearity, UInt64(rng_seed),
        String(photon_rate_schema), String(frame_schema))
end

_ingaas_graph_nonlinearity(::NullFrameNonlinearity, ::Type{T}) where {T} =
    NullFrameNonlinearity()

function _ingaas_graph_nonlinearity(model::SaturatingFrameNonlinearity,
        ::Type{T}) where {T}
    coefficient = T(model.coefficient)
    isfinite(coefficient) && coefficient >= zero(T) || throw(AlgorithmGraphError(
        "InGaAs nonlinearity coefficient must be finite and nonnegative"))
    return SaturatingFrameNonlinearity{T}(coefficient)
end

"""
    ingaas_detector_acquisition_node(name; rows, columns, pixel_sampling_m,
        spectral_channel, exposure_duration_s, quantum_efficiency, rng_seed,
        photon_rate_schema, frame_schema, ...)

Declare one complete-frame, global-shutter InGaAs acquisition. The `:photon_rate`
input contains cell-integrated photons per second on a metric detector grid,
with `pixel_sampling_m=(row_pitch_m, column_pitch_m)` and one explicit integrated
`spectral_channel` identifier. The `:frame` output has the same shape and caller
schema. Array axes are preserved; no transpose, resampling, or camera profile is
inferred.

Scalar `quantum_efficiency` is an effective constant across the declared band;
this node infers no wavelength response. If upstream spectral synthesis already
applies QE, set this value to one to avoid applying it twice. Photon shot noise,
uniform dark current and integration-time glow, full-well saturation, read noise,
post-readout gain, and optional ADC quantization follow `InGaAsSensor` and
`Detector` ordering. Dark and glow are Poisson charge processes even when photon
shot noise is disabled. `glow_e_per_pixel_s` does not model read-cadence glow.

Non-null `persistence_model` or `nonlinearity_model` uses ordinary stream
execution and is rejected by captured execution. Capture admission otherwise
uses the existing exact prepared detector gate: unit sampling/binning,
full-frame global shutter, null response/coupling/defect/correction/thermal/
background/readout-product models and a prepared device counter RNG. A separate
upstream pixel MTF supplies a calibrated spatial response. Defect maps and richer
acquisition models remain available through the general detector API; they are
not implicit options of this scalar node. GPU capture requires target validation
of the complete composed graph.
"""
function ingaas_detector_acquisition_node(name::Symbol;
        rows::Integer, columns::Integer, pixel_sampling_m::Tuple{Real,Real},
        spectral_channel::Symbol, exposure_duration_s::Real,
        quantum_efficiency::Real, rng_seed::Integer,
        photon_rate_schema::AbstractString, frame_schema::AbstractString,
        gain::Real=1, dark_current_e_per_pixel_s::Real=0,
        glow_e_per_pixel_s::Real=0, bits::Union{Nothing,Integer}=nothing,
        full_well_e::Union{Nothing,Real}=nothing, photon_noise::Bool=true,
        readout_noise::Bool=false, readout_noise_e::Real=0,
        persistence_model::Union{NullPersistence,ExponentialPersistence}=NullPersistence(),
        nonlinearity_model::Union{NullFrameNonlinearity,SaturatingFrameNonlinearity}=
            NullFrameNonlinearity(), T::Type{<:AbstractFloat}=Float32)
    config = _ingaas_detector_acquisition_config(T; rows, columns,
        pixel_sampling_m, spectral_channel, exposure_duration_s, quantum_efficiency,
        gain, dark_current_e_per_pixel_s, glow_e_per_pixel_s, bits, full_well_e,
        photon_noise, readout_noise, readout_noise_e, persistence_model,
        nonlinearity_model, rng_seed, photon_rate_schema, frame_schema)
    return algorithm_node(name, InGaAsDetectorAcquisitionNode{T}, config;
        props=NamedTuple())
end

function graph_node_ports(::Type{InGaAsDetectorAcquisitionNode{T}},
        config::InGaAsDetectorAcquisitionNodeConfig{T}) where {T}
    shape = (config.rows, config.columns)
    return (graph_port_contract(:photon_rate, :input, :data, T, shape,
            config.photon_rate_schema, :column_major),
        graph_port_contract(:frame, :output, :data, T, shape,
            config.frame_schema, :column_major))
end

function prepare_graph_node(::Type{InGaAsDetectorAcquisitionNode{T}},
        config::InGaAsDetectorAcquisitionNodeConfig{T}, ::NamedTuple{()},
        inputs::NamedTuple{(:photon_rate,)}, outputs::NamedTuple{(:frame,)},
        ::NamedTuple{()}, target) where {T}
    metadata = OpticalPlaneMetadata(DetectorPlane(), inputs.photon_rate;
        coordinate_domain=MetricCoordinates(), sampling=config.pixel_sampling_m,
        spectral=IntegratedSpectralChannel(config.spectral_channel),
        normalization=PhotonRateNormalization(), spatial_measure=CellIntegratedMeasure(),
        coherence=IncoherentIntensityAddition())
    photon_rate = IntensityMap(metadata, inputs.photon_rate)
    sensor = InGaAsSensor(glow_rate=config.glow_e_per_pixel_s,
        persistence_model=config.persistence_model, T=T)
    detector = Detector(exposure_duration=config.exposure_duration_s,
        qe=config.quantum_efficiency, psf_sampling=1, binning=1,
        gain=config.gain, dark_current=config.dark_current_e_per_pixel_s,
        bits=config.bits, full_well=config.full_well_e,
        output_type=_graph_adc_output_type(config.bits), noise=config.noise,
        sensor=sensor, response_model=NullFrameResponse(),
        nonlinearity_model=config.nonlinearity_model, T=T,
        backend=compute_device_backend(target))
    prepared = prepare_detector_acquisition(detector, photon_rate)
    size(output_frame(detector)) == size(outputs.frame) || throw(AlgorithmGraphError(
        "prepared InGaAs output shape does not match its graph frame port"))
    state = detector_acquisition_state(prepared)
    return _FrameDetectorAcquisitionOwner(prepared,
        _prepare_graph_rng(target, config.rng_seed), outputs.frame,
        state.latent_buffer, config.rng_seed)
end
