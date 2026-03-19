struct TelemetryRow{T<:AbstractFloat}
    iter::Int
    t::T
    wfe_rms::T
    strehl::T
    loop_gain::T
end

function TelemetryRow(iter::Int, t::Real; wfe_rms::Real=NaN, strehl::Real=NaN, loop_gain::Real=NaN)
    T = promote_type(
        typeof(float(t)),
        typeof(float(wfe_rms)),
        typeof(float(strehl)),
        typeof(float(loop_gain)),
    )
    return TelemetryRow{T}(iter, T(t), T(wfe_rms), T(strehl), T(loop_gain))
end

struct ClosedLoopTraceRow{T<:AbstractFloat}
    iter::Int
    t::T
    forcing_rms_nm::T
    residual_rms_nm::T
    strehl::T
    slope_norm::T
    command_norm::T
end

struct GSCClosedLoopTraceRow{T<:AbstractFloat}
    iter::Int
    t::T
    forcing_rms_nm::T
    residual_rms_nm::T
    strehl::T
    slope_norm::T
    mean_optical_gain::T
    command_norm::T
end

struct GSCAtmosphereReplayTraceRow{T<:AbstractFloat}
    iter::Int
    t::T
    ngs_forcing_rms_nm::T
    ngs_residual_rms_nm::T
    sci_residual_rms_nm::T
    ngs_strehl::T
    sci_strehl::T
    slope_norm::T
    mean_optical_gain::T
end

struct Telemetry{R}
    rows::Vector{R}
end

Telemetry{R}() where {R} = Telemetry{R}(R[])
Telemetry() = Telemetry{TelemetryRow{Float64}}(TelemetryRow{Float64}[])

Base.length(t::Telemetry) = length(t.rows)
Base.iterate(t::Telemetry, state=1) = state > length(t.rows) ? nothing : (t.rows[state], state + 1)
Base.getindex(t::Telemetry, idx::Int) = t.rows[idx]

struct ClosedLoopTrace{R}
    rows::Vector{R}
end

struct GSCClosedLoopTrace{R}
    rows::Vector{R}
end

struct GSCAtmosphereReplayTrace{R}
    rows::Vector{R}
end

Base.length(t::ClosedLoopTrace) = length(t.rows)
Base.iterate(t::ClosedLoopTrace, state=1) = state > length(t.rows) ? nothing : (t.rows[state], state + 1)
Base.getindex(t::ClosedLoopTrace, idx::Int) = t.rows[idx]
Base.eltype(::Type{<:ClosedLoopTrace{R}}) where {R} = R

Base.length(t::GSCClosedLoopTrace) = length(t.rows)
Base.iterate(t::GSCClosedLoopTrace, state=1) = state > length(t.rows) ? nothing : (t.rows[state], state + 1)
Base.getindex(t::GSCClosedLoopTrace, idx::Int) = t.rows[idx]
Base.eltype(::Type{<:GSCClosedLoopTrace{R}}) where {R} = R

Base.length(t::GSCAtmosphereReplayTrace) = length(t.rows)
Base.iterate(t::GSCAtmosphereReplayTrace, state=1) = state > length(t.rows) ? nothing : (t.rows[state], state + 1)
Base.getindex(t::GSCAtmosphereReplayTrace, idx::Int) = t.rows[idx]
Base.eltype(::Type{<:GSCAtmosphereReplayTrace{R}}) where {R} = R

function record!(t::Telemetry{TelemetryRow{T}}, iter::Int, time::Real;
    wfe_rms::Real=NaN, strehl::Real=NaN, loop_gain::Real=NaN) where {T<:AbstractFloat}
    row = TelemetryRow{T}(iter, T(time), T(wfe_rms), T(strehl), T(loop_gain))
    push!(t.rows, row)
    return t
end

function record!(t::Telemetry{R}, row::R) where {R}
    push!(t.rows, row)
    return t
end

function write_telemetry_csv(::AbstractString, ::Telemetry; kwargs...)
    throw(InvalidConfiguration("CSV.jl not available; load CSV to enable write_telemetry_csv."))
end

function write_telemetry_csv(::AbstractString, ::Union{ClosedLoopTrace,GSCClosedLoopTrace,GSCAtmosphereReplayTrace}; kwargs...)
    throw(InvalidConfiguration("CSV.jl not available; load CSV to enable write_telemetry_csv."))
end

function _trace_value_type(trace::AbstractMatrix{<:Real}, dt::Real, t0::Real)
    return promote_type(typeof(float(zero(eltype(trace)))), typeof(float(dt)), typeof(float(t0)))
end

function ClosedLoopTrace(trace::AbstractMatrix{<:Real}; dt::Real=1.0, t0::Real=0.0)
    size(trace, 2) == 5 || throw(DimensionMismatchError("ClosedLoopTrace expects a matrix with 5 columns"))
    T = _trace_value_type(trace, dt, t0)
    rows = Vector{ClosedLoopTraceRow{T}}(undef, size(trace, 1))
    @inbounds for iter in axes(trace, 1)
        t = T(t0 + (iter - 1) * dt)
        rows[iter] = ClosedLoopTraceRow(
            iter,
            t,
            T(trace[iter, 1]),
            T(trace[iter, 2]),
            T(trace[iter, 3]),
            T(trace[iter, 4]),
            T(trace[iter, 5]),
        )
    end
    return ClosedLoopTrace(rows)
end

function GSCClosedLoopTrace(trace::AbstractMatrix{<:Real}; dt::Real=1.0, t0::Real=0.0)
    size(trace, 2) == 6 || throw(DimensionMismatchError("GSCClosedLoopTrace expects a matrix with 6 columns"))
    T = _trace_value_type(trace, dt, t0)
    rows = Vector{GSCClosedLoopTraceRow{T}}(undef, size(trace, 1))
    @inbounds for iter in axes(trace, 1)
        t = T(t0 + (iter - 1) * dt)
        rows[iter] = GSCClosedLoopTraceRow(
            iter,
            t,
            T(trace[iter, 1]),
            T(trace[iter, 2]),
            T(trace[iter, 3]),
            T(trace[iter, 4]),
            T(trace[iter, 5]),
            T(trace[iter, 6]),
        )
    end
    return GSCClosedLoopTrace(rows)
end

function GSCAtmosphereReplayTrace(trace::AbstractMatrix{<:Real}; dt::Real=1.0, t0::Real=0.0)
    size(trace, 2) == 7 || throw(DimensionMismatchError("GSCAtmosphereReplayTrace expects a matrix with 7 columns"))
    T = _trace_value_type(trace, dt, t0)
    rows = Vector{GSCAtmosphereReplayTraceRow{T}}(undef, size(trace, 1))
    @inbounds for iter in axes(trace, 1)
        t = T(t0 + (iter - 1) * dt)
        rows[iter] = GSCAtmosphereReplayTraceRow(
            iter,
            t,
            T(trace[iter, 1]),
            T(trace[iter, 2]),
            T(trace[iter, 3]),
            T(trace[iter, 4]),
            T(trace[iter, 5]),
            T(trace[iter, 6]),
            T(trace[iter, 7]),
        )
    end
    return GSCAtmosphereReplayTrace(rows)
end
