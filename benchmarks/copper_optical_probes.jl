# Keep the Float64 readout and fixed probes from the October 2026 Copper
# convergence benchmark, recorded under benchmarks/results/platform.
function relative_rate_error(actual, reference)
    difference = norm(Float64.(actual) .- reference)
    scale = norm(Float64.(reference))
    return iszero(scale) ? (iszero(difference) ? 0.0 : Inf) : difference / scale
end

const POKE = 2e-8
const PROBES = ("flat", "piston", "tilt_x_20nm", "tilt_y_20nm", "focus_20nm",
    "astigmatism_20nm", "coma_20nm", "tilt_x_100nm", "focus_100nm",
    "fourier_4_20nm", "fourier_8_20nm", "fourier_12_20nm", "mixed_100nm")
# Identical detector support at every grid: geometric radius15 in a30x30 ROI.
const SUPPORT = [(i - 15.5)^2 + (j - 15.5)^2 <= 15^2 for i in 1:30, j in 1:30]
const OFFSETS = findall(SUPPORT)
const ORIGINS = ((1, 1), (33, 1), (33, 33), (1, 33)) # zero-based q1,q2,q3,q4

function opd_probe(n, name)
    axis = (2 .* collect(1:n) .- n .- 1) ./ n
    [Float32(probe_value(name, x, y)) for x in axis, y in axis]
end
function probe_value(name, x, y)
    name == "flat" && return 0.0
    name == "piston" && return 100e-9
    name == "tilt_x_20nm" && return 20e-9*x
    name == "tilt_y_20nm" && return 20e-9*y
    name == "focus_20nm" && return 20e-9*(2*(x*x+y*y)-1)
    name == "astigmatism_20nm" && return 20e-9*(x*x-y*y)
    name == "coma_20nm" && return 20e-9*(3*(x*x+y*y)-2)*x
    name == "tilt_x_100nm" && return 100e-9*x
    name == "focus_100nm" && return 100e-9*(2*(x*x+y*y)-1)
    name == "fourier_4_20nm" && return 20e-9*sinpi(4*x)*cospi(3*y)
    name == "fourier_8_20nm" && return 20e-9*sinpi(8*x)*cospi(5*y)
    name == "fourier_12_20nm" && return 20e-9*sinpi(12*x)*cospi(7*y)
    name == "mixed_100nm" && return 100e-9*(0.5sinpi(x+2y)+0.3cospi(3x-y)+0.2sinpi(7x+5y))
    error("unknown probe $name")
end

"""FGA I4Q differential after AOS (x,y) → image (row=y,column=x) transposition.
q1TL q2BL q3BR q4TR; fixed finest-flat incidence divisor.
"""
function raw_slopes(rate, divisor)
    result = Matrix{Float64}(undef, length(OFFSETS), 2)
    for (k, I) in enumerate(OFFSETS)
        i, j = Tuple(I)
        a, b, c, d = (Float64(rate[col+j, r+i]) for (r, col) in ORIGINS)
        result[k, 1] = (a-b+d-c)/divisor
        result[k, 2] = (a-d+b-c)/divisor
    end
    return vec(result)
end
function flat_divisor(rate)
    sum(sum(Float64(rate[c+I[2], r+I[1]]) for (r,c) in ORIGINS) for I in OFFSETS) / length(OFFSETS)
end

function probe_rates(system, n)
    input = graph_input(system.graph, :uncompensated_opd)
    fill!(graph_input(system.graph, :pdm_command), 0f0)
    rates = Array{Float32}(undef, 64, 64, length(PROBES))
    for (index, name) in enumerate(PROBES)
        copyto!(input, opd_probe(n, name))
        CUDA.synchronize()
        step_graph!(system.graph)
        rates[:,:,index] .= Array(graph_output(system.graph, :pwfs_photon_rate))
    end
    fill!(input, 0f0)
    CUDA.synchronize()
    @assert all(isfinite, rates) && minimum(rates) >= 0
    flat=rates[:,:,1]
    @assert relative_rate_error(rates[:,:,2],flat) < 5e-6 "piston must not change intensity"
    total=Float64(sum(Float64,flat))
    @assert maximum(abs(sum(Float64,rates[:,:,k])/total-1) for k in eachindex(PROBES)) < 5e-5
    @assert all(k->relative_rate_error(rates[:,:,k],flat)>1e-8,3:length(PROBES)) "capture must consume changed OPD"
    return rates
end
