abstract type AbstractCMOSOutputModel end
abstract type AbstractCMOSReadNoiseModel end

struct NullCMOSOutputModel <: AbstractCMOSOutputModel end
struct NullCMOSReadNoise <: AbstractCMOSReadNoiseModel end

struct CMOSReadNoiseMap{T<:AbstractFloat,A<:AbstractMatrix{T}} <:
    AbstractCMOSReadNoiseModel
    sigma::A
    function CMOSReadNoiseMap{T,A}(sigma::A) where
        {T<:AbstractFloat,A<:AbstractMatrix{T}}
        host_sigma = Array(sigma)
        isempty(host_sigma) && throw(InvalidConfiguration(
            "CMOSReadNoiseMap sigma map must not be empty"))
        minimum(host_sigma) >= zero(T) || throw(InvalidConfiguration(
            "CMOSReadNoiseMap sigma values must be >= 0"))
        all(isfinite, host_sigma) || throw(InvalidConfiguration(
            "CMOSReadNoiseMap sigma values must be finite"))
        return new{T,A}(sigma)
    end
end

function CMOSReadNoiseMap(sigma::AbstractMatrix;
    T::Type{<:AbstractFloat}=Float64,
    backend::AbstractArrayBackend=CPUBackend())
    converted = _to_backend_matrix(T.(Array(sigma)), backend)
    return validate_cmos_read_noise_model(
        CMOSReadNoiseMap{T,typeof(converted)}(converted))
end

struct StaticCMOSOutputPattern{T<:AbstractFloat,V1<:AbstractVector{T},V2<:AbstractVector{T}} <: AbstractCMOSOutputModel
    output_cols::Int
    gains::V1
    offsets::V2
end

function StaticCMOSOutputPattern(output_cols::Integer, gains::AbstractVector, offsets::AbstractVector;
    T::Type{<:AbstractFloat}=Float64, backend::AbstractArrayBackend=CPUBackend())
    backend = _resolve_array_backend(backend)
    backend_gains = _to_backend_vector(T.(gains), backend)
    backend_offsets = _to_backend_vector(T.(offsets), backend)
    return validate_cmos_output_model(
        StaticCMOSOutputPattern{T,typeof(backend_gains),typeof(backend_offsets)}(Int(output_cols), backend_gains, backend_offsets))
end

@kernel function apply_cmos_output_pattern_kernel!(frame, gains, offsets, output_cols::Int, n::Int, m::Int)
    i, j = @index(Global, NTuple)
    if i <= n && j <= m
        output_idx = cld(j, output_cols)
        @inbounds frame[i, j] = frame[i, j] * gains[output_idx] + offsets[output_idx]
    end
end

@kernel function add_row_noise_kernel!(frame, noise, sigma, n::Int, m::Int)
    i, j = @index(Global, NTuple)
    if i <= n && j <= m
        @inbounds frame[i, j] += sigma * noise[i, 1]
    end
end

struct CMOSSensor{T<:AbstractFloat,RN<:AbstractCMOSReadNoiseModel,
    O<:AbstractCMOSOutputModel,M<:AbstractFrameTimingModel} <: FrameSensorType
    column_readout_sigma::T
    row_readout_sigma::T
    readout_noise_model::RN
    output_model::O
    timing_model::M
end

function CMOSSensor(; column_readout_sigma::Real=0.0,
    row_readout_sigma::Real=0.0,
    readout_noise_model::AbstractCMOSReadNoiseModel=NullCMOSReadNoise(),
    output_model::AbstractCMOSOutputModel=NullCMOSOutputModel(),
    timing_model::AbstractFrameTimingModel=GlobalShutter(),
    T::Type{<:AbstractFloat}=Float64,
    backend::AbstractArrayBackend=CPUBackend())
    backend = _resolve_array_backend(backend)
    column_readout_sigma >= 0 || throw(InvalidConfiguration("CMOSSensor column_readout_sigma must be >= 0"))
    row_readout_sigma >= 0 || throw(InvalidConfiguration("CMOSSensor row_readout_sigma must be >= 0"))
    converted_noise = convert_cmos_read_noise_model(readout_noise_model, T, backend)
    validated_noise = validate_cmos_read_noise_model(converted_noise)
    converted_output = convert_cmos_output_model(output_model, T, backend)
    validated_output = validate_cmos_output_model(converted_output)
    validated_timing = validate_frame_timing_model(convert_frame_timing_model(timing_model, T))
    return CMOSSensor{T,typeof(validated_noise),typeof(validated_output),typeof(validated_timing)}(
        T(column_readout_sigma), T(row_readout_sigma), validated_noise,
        validated_output, validated_timing)
end

detector_sensor_symbol(::CMOSSensor) = :cmos
supports_column_readout_noise(::CMOSSensor) = true
supports_detector_defect_maps(::CMOSSensor) = true
supports_shutter_timing(::CMOSSensor) = true
default_frame_timing_model(sensor::CMOSSensor; T::Type{<:AbstractFloat}=Float64) = sensor.timing_model
default_response_model(::CMOSSensor; T::Type{<:AbstractFloat}=Float64, backend::AbstractArrayBackend=CPUBackend()) =
    NullFrameResponse()
is_null_cmos_output_model(::AbstractCMOSOutputModel) = false
is_null_cmos_output_model(::NullCMOSOutputModel) = true
is_null_cmos_read_noise_model(::AbstractCMOSReadNoiseModel) = false
is_null_cmos_read_noise_model(::NullCMOSReadNoise) = true

convert_cmos_read_noise_model(::NullCMOSReadNoise, ::Type{T}, backend) where {T<:AbstractFloat} =
    NullCMOSReadNoise()

function convert_cmos_read_noise_model(model::CMOSReadNoiseMap, ::Type{T}, backend) where {T<:AbstractFloat}
    sigma = _to_backend_matrix(T.(Array(model.sigma)), backend)
    return CMOSReadNoiseMap{T,typeof(sigma)}(sigma)
end

validate_cmos_read_noise_model(::NullCMOSReadNoise) = NullCMOSReadNoise()
validate_cmos_read_noise_model(model::CMOSReadNoiseMap) = model

convert_cmos_output_model(::NullCMOSOutputModel, ::Type{T}, backend) where {T<:AbstractFloat} = NullCMOSOutputModel()

function convert_cmos_output_model(model::StaticCMOSOutputPattern, ::Type{T}, backend) where {T<:AbstractFloat}
    gains = _to_backend_vector(T.(Array(model.gains)), backend)
    offsets = _to_backend_vector(T.(Array(model.offsets)), backend)
    return StaticCMOSOutputPattern{T,typeof(gains),typeof(offsets)}(model.output_cols, gains, offsets)
end

validate_cmos_output_model(::NullCMOSOutputModel) = NullCMOSOutputModel()

function validate_cmos_output_model(model::StaticCMOSOutputPattern)
    model.output_cols > 0 || throw(InvalidConfiguration("StaticCMOSOutputPattern output_cols must be > 0"))
    length(model.gains) > 0 || throw(InvalidConfiguration("StaticCMOSOutputPattern gains must not be empty"))
    length(model.gains) == length(model.offsets) ||
        throw(InvalidConfiguration("StaticCMOSOutputPattern gains and offsets must have matching length"))
    minimum(model.gains) >= zero(eltype(model.gains)) ||
        throw(InvalidConfiguration("StaticCMOSOutputPattern gains must be >= 0"))
    return model
end

function sampling_wallclock_time(sensor::CMOSSensor, integration_time, frame_size::Tuple{Int,Int},
    window::Union{Nothing,FrameWindow}, ::Type{T}) where {T<:AbstractFloat}
    is_global_shutter(sensor.timing_model) && return T(integration_time)
    active_rows = window === nothing ? frame_size[1] : length(window.rows)
    return T(integration_time) + T(cld(active_rows, sensor.timing_model.row_group_size)) * T(sensor.timing_model.line_time)
end

function apply_sensor_readout_noise!(sensor::CMOSSensor, det::Detector,
    rng::AbstractRNG)
    apply_cmos_column_noise!(sensor, det, rng)
    apply_cmos_row_noise!(sensor, det, rng)
    apply_cmos_pixel_noise!(sensor.readout_noise_model, det, rng)
    return det.state.frame
end

function apply_cmos_column_noise!(sensor::CMOSSensor, det::Detector,
    rng::AbstractRNG)
    sigma = sensor.column_readout_sigma
    sigma <= zero(sigma) && return det.state.frame
    randn_frame_noise!(det, rng, det.state.noise_buffer)
    return apply_column_noise!(execution_style(det.state.frame), det.state.frame, det.state.noise_buffer, sigma)
end

function apply_cmos_row_noise!(sensor::CMOSSensor, det::Detector,
    rng::AbstractRNG)
    sigma = sensor.row_readout_sigma
    sigma <= zero(sigma) && return det.state.frame
    randn_frame_noise!(det, rng, det.state.noise_buffer)
    return apply_row_noise!(execution_style(det.state.frame), det.state.frame,
        det.state.noise_buffer, sigma)
end

apply_cmos_pixel_noise!(::NullCMOSReadNoise, det::Detector,
    rng::AbstractRNG) = det.state.frame

function apply_cmos_pixel_noise!(model::CMOSReadNoiseMap, det::Detector,
    rng::AbstractRNG)
    size(model.sigma) == size(det.state.frame) ||
        throw(DimensionMismatchError(
            "CMOSReadNoiseMap sigma size must match detector frame size"))
    randn_frame_noise!(det, rng, det.state.noise_buffer)
    det.state.frame .+= model.sigma .* det.state.noise_buffer
    return det.state.frame
end

function apply_column_noise!(::ScalarCPUStyle, frame, noise, sigma)
    n, m = size(frame)
    @inbounds for j in 1:m
        offset = sigma * noise[1, j]
        for i in 1:n
            frame[i, j] += offset
        end
    end
    return frame
end

function apply_column_noise!(style::AcceleratorStyle, frame, noise, sigma)
    n, m = size(frame)
    launch_kernel!(style, add_column_noise_kernel!, frame, noise, sigma, n, m; ndrange=(n, m))
    return frame
end

function apply_row_noise!(::ScalarCPUStyle, frame, noise, sigma)
    @inbounds for j in axes(frame, 2), i in axes(frame, 1)
        frame[i, j] += sigma * noise[i, 1]
    end
    return frame
end

function apply_row_noise!(style::AcceleratorStyle, frame, noise, sigma)
    n, m = size(frame)
    launch_kernel!(style, add_row_noise_kernel!, frame, noise, sigma, n, m;
        ndrange=(n, m))
    return frame
end

function apply_post_readout_gain!(::CMOSSensor, det::Detector)
    det.state.frame .*= det.params.gain
    apply_output_model!(execution_style(det.state.frame), det.params.sensor.output_model, det.state.frame)
    return det.state.frame
end

apply_output_model!(::ExecutionStyle, ::NullCMOSOutputModel, frame) = frame
apply_output_model!(::ScalarCPUStyle, ::NullCMOSOutputModel, frame) = frame

function apply_output_model!(::ScalarCPUStyle, model::StaticCMOSOutputPattern, frame)
    n, m = size(frame)
    n_outputs = length(model.gains)
    n_outputs * model.output_cols >= m ||
        throw(DimensionMismatchError("StaticCMOSOutputPattern does not cover detector columns"))
    @inbounds for j in 1:m
        output_idx = cld(j, model.output_cols)
        gain = model.gains[output_idx]
        offset = model.offsets[output_idx]
        for i in 1:n
            frame[i, j] = frame[i, j] * gain + offset
        end
    end
    return frame
end

function apply_output_model!(style::AcceleratorStyle, model::StaticCMOSOutputPattern, frame)
    n, m = size(frame)
    length(model.gains) * model.output_cols >= m ||
        throw(DimensionMismatchError("StaticCMOSOutputPattern does not cover detector columns"))
    launch_kernel!(style, apply_cmos_output_pattern_kernel!, frame, model.gains, model.offsets, model.output_cols, n, m; ndrange=(n, m))
    return frame
end

function _require_batched_sensor_compat(sensor::CMOSSensor)
    sensor.column_readout_sigma <= zero(sensor.column_readout_sigma) ||
        throw(InvalidConfiguration("batched detector capture does not yet support CMOSSensor column_readout_sigma"))
    sensor.row_readout_sigma <= zero(sensor.row_readout_sigma) ||
        throw(InvalidConfiguration("batched detector capture does not yet support CMOSSensor row_readout_sigma"))
    is_null_cmos_read_noise_model(sensor.readout_noise_model) ||
        throw(InvalidConfiguration("batched detector capture does not yet support CMOSSensor readout-noise maps"))
    is_null_cmos_output_model(sensor.output_model) ||
        throw(InvalidConfiguration("batched detector capture does not yet support CMOSSensor output-group patterns"))
    is_global_shutter(sensor.timing_model) ||
        throw(InvalidConfiguration("batched detector capture does not yet support CMOSSensor rolling-shutter timing"))
    return nothing
end

function _batched_post_readout_gain!(::CMOSSensor, det::Detector, cube::AbstractArray)
    isone(det.params.gain) || (cube .*= det.params.gain)
    return cube
end
