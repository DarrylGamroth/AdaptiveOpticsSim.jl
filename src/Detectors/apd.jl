function ensure_channel_buffers!(det::APDDetector, dims::Tuple{Int,Int})
    if size(det.state.channels) != dims
        det.state.channels = similar(det.state.channels, dims...)
    end
    if size(det.state.noise_buffer) != dims
        det.state.noise_buffer = similar(det.state.noise_buffer, dims...)
    end
    if det.state.output_buffer !== nothing && size(det.state.output_buffer) != dims
        det.state.output_buffer = similar(det.state.output_buffer, dims...)
        fill!(det.state.output_buffer, zero(eltype(det.state.output_buffer)))
    end
    return det
end

function apply_apd_gain_map!(det::APDDetector)
    isnothing(det.channel_gain_map) && return det.state.channels
    size(det.channel_gain_map) == size(det.state.channels) ||
        throw(DimensionMismatchError("APDDetector channel_gain_map size must match counting-channel size"))
    det.state.channels .*= det.channel_gain_map
    return det.state.channels
end

function apply_apd_dark_counts!(det::APDDetector, exposure_time::Real)
    dark = det.params.dark_count_rate * exposure_time
    dark <= 0 && return det.state.channels
    det.state.channels .+= dark
    return det.state.channels
end

apply_apd_dead_time!(det::APDDetector) = apply_apd_dead_time!(det.params.dead_time_model, det)
apply_apd_dead_time!(::NoDeadTime, det::APDDetector) = det.state.channels

function apply_apd_dead_time!(model::NonParalyzableDeadTime, det::APDDetector)
    exposure_time = det.params.integration_time
    exposure_time > zero(exposure_time) || return det.state.channels
    scale = model.dead_time / exposure_time
    scale <= zero(scale) && return det.state.channels
    @. det.state.channels = det.state.channels / (1 + det.state.channels * scale)
    return det.state.channels
end

apply_apd_counting_noise!(det::APDDetector{NoiseNone}, rng::AbstractRNG) = det.state.channels

function apply_apd_counting_noise!(det::APDDetector{NoisePhoton}, rng::AbstractRNG)
    poisson_noise!(rng, det.state.channels)
    return det.state.channels
end

function apply_apd_gain!(det::APDDetector)
    det.state.channels .*= det.params.gain
    return det.state.channels
end

function write_channel_output!(det::APDDetector)
    output = det.state.output_buffer
    output === nothing && return det.state.channels
    output_eltype = eltype(output)
    if output_eltype <: Integer
        lo = typemin(output_eltype)
        hi = typemax(output_eltype)
        @. output = output_eltype(clamp(round(det.state.channels), lo, hi))
    else
        copyto!(output, det.state.channels)
    end
    return output
end

function capture!(det::APDDetector, channels::AbstractMatrix{T}; rng::AbstractRNG=Random.default_rng()) where {T<:AbstractFloat}
    ensure_channel_buffers!(det, size(channels))
    copyto!(det.state.channels, channels)
    det.state.channels .*= det.params.qe * det.params.integration_time
    apply_apd_gain_map!(det)
    apply_apd_dark_counts!(det, det.params.integration_time)
    apply_apd_counting_noise!(det, rng)
    apply_apd_dead_time!(det)
    apply_apd_gain!(det)
    return write_channel_output!(det)
end

function capture!(det::APDDetector, channels::AbstractMatrix{T}, rng::AbstractRNG) where {T<:AbstractFloat}
    return capture!(det, channels; rng=rng)
end
