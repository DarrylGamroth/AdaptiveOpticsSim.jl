"""Graph operation for one complete-frame HgCdTe linear-avalanche acquisition."""
struct HgCdTeAvalancheDetectorAcquisitionNode{T<:AbstractFloat} end

"""Immutable configuration for one complete-frame HgCdTe linear-avalanche acquisition."""
struct HgCdTeAvalancheDetectorAcquisitionNodeConfig{T,N,S,L}
    rows::Int
    columns::Int
    pixel_sampling_m::Tuple{T,T}
    spectral_channel::Symbol
    exposure_duration_s::T
    quantum_efficiency::T
    gain::T
    dark_current_e_per_pixel_s::T
    bits::Union{Nothing,Int}
    full_well_e::Union{Nothing,T}
    noise::N
    sensor::S
    nonlinearity_model::L
    rng_seed::UInt64
    photon_rate_schema::String
    frame_schema::String
end

function _hgcdte_avalanche_detector_acquisition_config(::Type{T};
        rows::Integer, columns::Integer, pixel_sampling_m::Tuple{Real,Real},
        spectral_channel::Symbol, exposure_duration_s::Real,
        quantum_efficiency::Real, gain::Real, dark_current_e_per_pixel_s::Real,
        glow_e_per_pixel_s::Real, avalanche_gain::Real, excess_noise_factor::Real,
        multiplication_model::Union{ConditionalGammaAvalancheMultiplication,
            ClippedGaussianAvalancheMultiplicationApproximation},
        bits::Union{Nothing,Integer},
        full_well_e::Union{Nothing,Real}, photon_noise::Bool,
        readout_noise::Bool, readout_noise_e::Real,
        persistence_model::Union{NullPersistence,ExponentialPersistence},
        nonlinearity_model::Union{NullFrameNonlinearity,SaturatingFrameNonlinearity},
        rng_seed::Integer, photon_rate_schema::AbstractString,
        frame_schema::AbstractString) where {T<:AbstractFloat}
    T in (Float32, Float64) || throw(AlgorithmGraphError(
        "HgCdTe linear-avalanche graph acquisition requires Float32 or Float64"))
    n_rows, n_columns = Int(rows), Int(columns)
    n_rows > 0 && n_columns > 0 || throw(AlgorithmGraphError(
        "HgCdTe linear-avalanche photon-rate dimensions must be positive"))
    sampling = (T(pixel_sampling_m[1]), T(pixel_sampling_m[2]))
    all(value -> isfinite(value) && value > zero(T), sampling) ||
        throw(AlgorithmGraphError("HgCdTe linear-avalanche pixel_sampling_m must be finite and positive"))
    isempty(string(spectral_channel)) && throw(AlgorithmGraphError(
        "HgCdTe linear-avalanche spectral_channel must not be empty"))
    exposure, qe, linear_gain = T(exposure_duration_s), T(quantum_efficiency), T(gain)
    dark, glow = T(dark_current_e_per_pixel_s), T(glow_e_per_pixel_s)
    isfinite(exposure) && exposure > zero(T) || throw(AlgorithmGraphError(
        "HgCdTe linear-avalanche exposure_duration_s must be finite and positive"))
    isfinite(qe) && zero(T) <= qe <= one(T) || throw(AlgorithmGraphError(
        "HgCdTe linear-avalanche quantum_efficiency must lie in [0, 1]"))
    isfinite(linear_gain) && linear_gain > zero(T) || throw(AlgorithmGraphError(
        "HgCdTe linear-avalanche gain must be finite and positive"))
    isfinite(dark) && dark >= zero(T) || throw(AlgorithmGraphError(
        "HgCdTe linear-avalanche dark_current_e_per_pixel_s must be finite and nonnegative"))
    isfinite(glow) && glow >= zero(T) || throw(AlgorithmGraphError(
        "HgCdTe linear-avalanche glow_e_per_pixel_s must be finite and nonnegative"))
    output_bits = isnothing(bits) ? nothing : Int(bits)
    isnothing(output_bits) || 1 <= output_bits <= 64 || throw(AlgorithmGraphError(
        "HgCdTe linear-avalanche bits must lie in 1:64"))
    isnothing(output_bits) || !isnothing(full_well_e) || throw(AlgorithmGraphError(
        "HgCdTe linear-avalanche bits requires full_well_e"))
    full_well = isnothing(full_well_e) ? nothing : T(full_well_e)
    isnothing(full_well) || isfinite(full_well) && full_well > zero(T) ||
        throw(AlgorithmGraphError("HgCdTe linear-avalanche full_well_e must be finite and positive"))
    0 <= rng_seed <= typemax(UInt64) || throw(AlgorithmGraphError(
        "HgCdTe linear-avalanche rng_seed must lie in the UInt64 range"))
    isempty(photon_rate_schema) && throw(AlgorithmGraphError(
        "HgCdTe linear-avalanche photon_rate_schema must not be empty"))
    isempty(frame_schema) && throw(AlgorithmGraphError(
        "HgCdTe linear-avalanche frame_schema must not be empty"))
    noise = _frame_detector_noise(T, "HgCdTeAvalanche", photon_noise,
        readout_noise, readout_noise_e)
    # The sensor validates and converts the immutable persistence coefficients;
    # evolving latent charge is created only during detector preparation.
    avalanche, factor = T(avalanche_gain), T(excess_noise_factor)
    isfinite(avalanche) && avalanche >= one(T) || throw(AlgorithmGraphError(
        "HgCdTe linear-avalanche avalanche_gain must be finite and >= 1"))
    isfinite(factor) && factor >= one(T) || throw(AlgorithmGraphError(
        "HgCdTe linear-avalanche excess_noise_factor must be finite and >= 1"))
    sensor = HgCdTeAvalancheArraySensor(; glow_rate=glow, persistence_model,
        avalanche_gain=avalanche, excess_noise_factor=factor, multiplication_model, T)
    nonlinearity = _hgcdte_avalanche_graph_nonlinearity(nonlinearity_model, T)
    return HgCdTeAvalancheDetectorAcquisitionNodeConfig(n_rows, n_columns, sampling,
        spectral_channel, exposure, qe, linear_gain, dark, output_bits,
        full_well, noise, sensor, nonlinearity, UInt64(rng_seed),
        String(photon_rate_schema), String(frame_schema))
end

_hgcdte_avalanche_graph_nonlinearity(::NullFrameNonlinearity, ::Type{T}) where {T} =
    NullFrameNonlinearity()

function _hgcdte_avalanche_graph_nonlinearity(model::SaturatingFrameNonlinearity,
        ::Type{T}) where {T}
    coefficient = T(model.coefficient)
    isfinite(coefficient) && coefficient >= zero(T) || throw(AlgorithmGraphError(
        "HgCdTe linear-avalanche nonlinearity coefficient must be finite and nonnegative"))
    return SaturatingFrameNonlinearity{T}(coefficient)
end

"""
    hgcdte_avalanche_detector_acquisition_node(name; rows, columns, pixel_sampling_m,
        spectral_channel, exposure_duration_s, quantum_efficiency, rng_seed,
        photon_rate_schema, frame_schema, ...)

Declare one complete-frame, global-shutter HgCdTe linear-avalanche array
acquisition, including SAPHIRA configurations supplied by the caller. The
`:photon_rate` input contains cell-integrated photons per second on a metric
detector grid, with `pixel_sampling_m=(row_pitch_m, column_pitch_m)` and an
explicit integrated `spectral_channel`. The `:frame` output preserves shape,
axes and caller schema. Scalar `quantum_efficiency` is effective across that
band; use one if upstream synthesis already applies QE.

This node prepares `HgCdTeAvalancheArraySensor` with single-read sampling and
uses the existing `Detector` acquisition pipeline. `avalanche_gain` is the
mean pre-readout multiplication (at least one); `excess_noise_factor` uses
F = E[M²]/E[M]² (at least one). `multiplication_model` selects either the CPU
`ConditionalGammaAvalancheMultiplication` reference or the numerical sensor's
default `ClippedGaussianAvalancheMultiplicationApproximation`, whose qualified
regime is q/(F−1) ≥ 25. F = 1 is deterministic. `readout_noise_e` is noise
before post-readout `gain`, after avalanche multiplication. `full_well_e` is
the multiplied-charge saturation limit; ADC quantization follows conversion
gain. Dark current and `glow_e_per_pixel_s` remain Poisson charge processes
even when photon shot noise is disabled. Glow scales with integration time.

Persistence and nonlinearity follow the numerical sensor's ordering. Spatial
response belongs to an explicit upstream operation. No camera parameters,
wavelength-dependent QE, read-cadence glow, or nondestructive-read timing are
inferred. CPU stream execution is supported; this adapter adds no accelerator
capture qualification.
"""
function hgcdte_avalanche_detector_acquisition_node(name::Symbol;
        rows::Integer, columns::Integer, pixel_sampling_m::Tuple{Real,Real},
        spectral_channel::Symbol, exposure_duration_s::Real,
        quantum_efficiency::Real, rng_seed::Integer,
        photon_rate_schema::AbstractString, frame_schema::AbstractString,
        gain::Real=1, dark_current_e_per_pixel_s::Real=0,
        glow_e_per_pixel_s::Real=0, avalanche_gain::Real=1,
        excess_noise_factor::Real=1,
        multiplication_model::Union{ConditionalGammaAvalancheMultiplication,
            ClippedGaussianAvalancheMultiplicationApproximation}=
            ClippedGaussianAvalancheMultiplicationApproximation(),
        bits::Union{Nothing,Integer}=nothing,
        full_well_e::Union{Nothing,Real}=nothing, photon_noise::Bool=true,
        readout_noise::Bool=false, readout_noise_e::Real=0,
        persistence_model::Union{NullPersistence,ExponentialPersistence}=NullPersistence(),
        nonlinearity_model::Union{NullFrameNonlinearity,SaturatingFrameNonlinearity}=
            NullFrameNonlinearity(), T::Type{<:AbstractFloat}=Float32)
    config = _hgcdte_avalanche_detector_acquisition_config(T; rows, columns,
        pixel_sampling_m, spectral_channel, exposure_duration_s, quantum_efficiency,
        gain, dark_current_e_per_pixel_s, glow_e_per_pixel_s, avalanche_gain,
        excess_noise_factor, multiplication_model, bits, full_well_e,
        photon_noise, readout_noise, readout_noise_e, persistence_model,
        nonlinearity_model, rng_seed, photon_rate_schema, frame_schema)
    return algorithm_node(name, HgCdTeAvalancheDetectorAcquisitionNode{T}, config;
        props=NamedTuple())
end

function graph_node_ports(::Type{HgCdTeAvalancheDetectorAcquisitionNode{T}},
        config::HgCdTeAvalancheDetectorAcquisitionNodeConfig{T}) where {T}
    shape = (config.rows, config.columns)
    return (graph_port_contract(:photon_rate, :input, :data, T, shape,
            config.photon_rate_schema, :column_major),
        graph_port_contract(:frame, :output, :data, T, shape,
            config.frame_schema, :column_major))
end

function prepare_graph_node(::Type{HgCdTeAvalancheDetectorAcquisitionNode{T}},
        config::HgCdTeAvalancheDetectorAcquisitionNodeConfig{T}, ::NamedTuple{()},
        inputs::NamedTuple{(:photon_rate,)}, outputs::NamedTuple{(:frame,)},
        ::NamedTuple{()}, target) where {T}
    metadata = OpticalPlaneMetadata(DetectorPlane(), inputs.photon_rate;
        coordinate_domain=MetricCoordinates(), sampling=config.pixel_sampling_m,
        spectral=IntegratedSpectralChannel(config.spectral_channel),
        normalization=PhotonRateNormalization(), spatial_measure=CellIntegratedMeasure(),
        coherence=IncoherentIntensityAddition())
    photon_rate = IntensityMap(metadata, inputs.photon_rate)
    detector = Detector(exposure_duration=config.exposure_duration_s,
        qe=config.quantum_efficiency, psf_sampling=1, binning=1,
        gain=config.gain, dark_current=config.dark_current_e_per_pixel_s,
        bits=config.bits, full_well=config.full_well_e,
        output_type=_graph_adc_output_type(config.bits), noise=config.noise,
        sensor=config.sensor, response_model=NullFrameResponse(),
        nonlinearity_model=config.nonlinearity_model, T=T,
        backend=compute_device_backend(target))
    prepared = prepare_detector_acquisition(detector, photon_rate)
    size(output_frame(detector)) == size(outputs.frame) || throw(AlgorithmGraphError(
        "prepared HgCdTe linear-avalanche output shape does not match its graph frame port"))
    state = detector_acquisition_state(prepared)
    return _FrameDetectorAcquisitionOwner(prepared,
        _prepare_graph_rng(target, config.rng_seed), outputs.frame,
        state.latent_buffer, config.rng_seed)
end
