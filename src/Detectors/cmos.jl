abstract type AbstractCMOSOutputModel end

struct NullCMOSOutputModel <: AbstractCMOSOutputModel end

struct StaticCMOSOutputPattern{T<:AbstractFloat,V1<:AbstractVector{T},V2<:AbstractVector{T}} <: AbstractCMOSOutputModel
    output_cols::Int
    gains::V1
    offsets::V2
end

function StaticCMOSOutputPattern(output_cols::Integer, gains::AbstractVector, offsets::AbstractVector;
    T::Type{<:AbstractFloat}=Float64, backend=CPUBackend())
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

struct CMOSSensor{T<:AbstractFloat,O<:AbstractCMOSOutputModel,M<:AbstractFrameTimingModel} <: FrameSensorType
    column_readout_sigma::T
    output_model::O
    timing_model::M
end

function CMOSSensor(; column_readout_sigma::Real=0.0, output_model::AbstractCMOSOutputModel=NullCMOSOutputModel(),
    timing_model::AbstractFrameTimingModel=GlobalShutter(), T::Type{<:AbstractFloat}=Float64, backend=CPUBackend())
    backend = _resolve_array_backend(backend)
    column_readout_sigma >= 0 || throw(InvalidConfiguration("CMOSSensor column_readout_sigma must be >= 0"))
    converted_output = convert_cmos_output_model(output_model, T, backend)
    validated_output = validate_cmos_output_model(converted_output)
    validated_timing = validate_frame_timing_model(convert_frame_timing_model(timing_model, T))
    return CMOSSensor{T,typeof(validated_output),typeof(validated_timing)}(T(column_readout_sigma), validated_output, validated_timing)
end

detector_sensor_symbol(::CMOSSensor) = :cmos
supports_column_readout_noise(::CMOSSensor) = true
supports_detector_defect_maps(::CMOSSensor) = true
supports_shutter_timing(::CMOSSensor) = true
default_frame_timing_model(sensor::CMOSSensor; T::Type{<:AbstractFloat}=Float64) = sensor.timing_model
is_null_cmos_output_model(::AbstractCMOSOutputModel) = false
is_null_cmos_output_model(::NullCMOSOutputModel) = true

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
    return T(integration_time) + T(active_rows) * T(sensor.timing_model.line_time)
end

function apply_sensor_statistics!(sensor::CMOSSensor, det::Detector, rng::AbstractRNG)
    sigma = sensor.column_readout_sigma
    sigma <= zero(sigma) && return det.state.frame
    randn_frame_noise!(det, rng, det.state.noise_buffer)
    return apply_column_noise!(execution_style(det.state.frame), det.state.frame, det.state.noise_buffer, sigma)
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
    is_null_cmos_output_model(sensor.output_model) ||
        throw(InvalidConfiguration("batched detector capture does not yet support CMOSSensor output-group patterns"))
    is_global_shutter(sensor.timing_model) ||
        throw(InvalidConfiguration("batched detector capture does not yet support CMOSSensor rolling-shutter timing"))
    return nothing
end

_batched_post_readout_gain!(::CMOSSensor, det::Detector, cube::AbstractArray) = (cube .*= det.params.gain; cube)
