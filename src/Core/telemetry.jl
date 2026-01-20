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

struct Telemetry{R}
    rows::Vector{R}
end

Telemetry{R}() where {R} = Telemetry{R}(R[])
Telemetry() = Telemetry{TelemetryRow{Float64}}(TelemetryRow{Float64}[])

Base.length(t::Telemetry) = length(t.rows)
Base.iterate(t::Telemetry, state=1) = state > length(t.rows) ? nothing : (t.rows[state], state + 1)
Base.getindex(t::Telemetry, idx::Int) = t.rows[idx]

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
