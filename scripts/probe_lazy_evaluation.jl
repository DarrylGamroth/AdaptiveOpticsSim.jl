#!/usr/bin/env julia

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using AdaptiveOpticsSim
using BenchmarkTools
using LinearAlgebra

const AOS = AdaptiveOpticsSim

BenchmarkTools.DEFAULT_PARAMETERS.seconds = parse(Float64, get(ENV, "AOSIM_BENCH_SECONDS", "0.5"))
BenchmarkTools.DEFAULT_PARAMETERS.samples = parse(Int, get(ENV, "AOSIM_BENCH_SAMPLES", "8"))

struct LazyShiftedCoords{
    T<:AbstractFloat,
    X<:AbstractMatrix{T},
    Y<:AbstractMatrix{T},
    P<:AbstractVector{Int},
} <: AbstractVector{Complex{T}}
    x::X
    y::Y
    positions::P
    beta_x::T
    beta_y::T
    scale::T
end

Base.IndexStyle(::Type{<:LazyShiftedCoords}) = IndexLinear()
Base.size(coords::LazyShiftedCoords) = (length(coords.positions),)

@inline function Base.getindex(coords::LazyShiftedCoords{T}, i::Int) where {T}
    idx = @inbounds coords.positions[i]
    return @inbounds complex(coords.x[idx] * coords.scale + coords.beta_x,
        coords.y[idx] * coords.scale + coords.beta_y)
end

function fill_shifted_coords!(
    out::AbstractVector{Complex{T}},
    x::AbstractMatrix{T},
    y::AbstractMatrix{T},
    positions::AbstractVector{Int},
    beta_x::T,
    beta_y::T,
    scale::T,
) where {T<:AbstractFloat}
    @inbounds for k in eachindex(out, positions)
        idx = positions[k]
        out[k] = complex(x[idx] * scale + beta_x, y[idx] * scale + beta_y)
    end
    return out
end

function covariance_materialized!(cov, iz, jz, x1, y1, x2, y2, positions, cst, var_term, inv_L0, fraction)
    fill_shifted_coords!(iz, x1, y1, positions, 1.2, -0.4, 0.94)
    fill_shifted_coords!(jz, x2, y2, positions, -0.7, 0.2, 0.91)
    AOS._covariance_matrix!(cov, iz, jz, cst, var_term, inv_L0, fraction)
    return cov
end

function covariance_lazy!(cov, x1, y1, x2, y2, positions, cst, var_term, inv_L0, fraction)
    iz = LazyShiftedCoords(x1, y1, positions, 1.2, -0.4, 0.94)
    jz = LazyShiftedCoords(x2, y2, positions, -0.7, 0.2, 0.91)
    AOS._covariance_matrix!(cov, iz, jz, cst, var_term, inv_L0, fraction)
    return cov
end

function field_fill_package!(out, tel, src)
    AOS.fill_telescope_field!(out, tel, src; zero_padding=2)
    return out
end

function field_fill_workspace_amplitude!(out, amplitude, tel, src)
    n = tel.params.resolution
    n_pad = n * 2
    fill!(out, zero(eltype(out)))
    @. amplitude = sqrt(tel.state.pupil_reflectivity)
    opd_to_cycles = 2 / AOS.wavelength(src)
    amp_scale = sqrt(AOS.photon_flux(src) * tel.params.sampling_time * (tel.params.diameter / tel.params.resolution)^2)
    ox, oy = AOS.field_embedding_offsets(n, n_pad)
    @views @. out[ox+1:ox+n, oy+1:oy+n] = amp_scale * amplitude * cispi(opd_to_cycles * tel.state.opd)
    if iseven(n_pad)
        phase_shift = -pi * (n_pad + 1) / n_pad
        AOS.apply_centering_phase!(AOS.ScalarCPUStyle(), out, phase_shift)
    end
    return out
end

function summarize(name, trial)
    println(name)
    println("  median_ns: ", Float64(BenchmarkTools.median(trial).time))
    println("  memory_bytes: ", BenchmarkTools.memory(trial))
    println("  allocs: ", BenchmarkTools.allocs(trial))
end

function main()
    println("lazy_evaluation_probe")

    tel = AOS.Telescope(resolution=96, diameter=8.0, sampling_time=1e-3, pupil_reflectivity=0.72)
    src = AOS.Source(coordinates=(0.15, 30.0))
    out = Matrix{ComplexF64}(undef, 192, 192)
    amplitude = similar(tel.state.pupil_reflectivity)
    field_fill_package!(out, tel, src)
    field_fill_workspace_amplitude!(out, amplitude, tel, src)
    summarize("field_fill_package_fused_lazy_sqrt", @benchmark field_fill_package!($out, $tel, $src))
    summarize("field_fill_materialized_amplitude_workspace", @benchmark field_fill_workspace_amplitude!($out, $amplitude, $tel, $src))

    sampling = 12
    mask = trues(sampling, sampling)
    positions = findall(vec(mask))
    rotations = [0.0, 0.04]
    offsets_x = [0.0, 0.01]
    offsets_y = [0.0, -0.02]
    guide_x, guide_y = AOS._guide_star_grids(sampling, 8.0, rotations, offsets_x, offsets_y)
    target_x, target_y = AOS._guide_star_grid(sampling, 8.0, 0.0, 0.0, 0.0)
    iz = Vector{ComplexF64}(undef, length(positions))
    jz = similar(iz)
    cov = Matrix{Float64}(undef, length(positions), length(positions))
    cst, var_term, inv_L0 = AOS._covariance_constants(0.16, 25.0)

    covariance_materialized!(cov, iz, jz, @view(guide_x[:, :, 1]), @view(guide_y[:, :, 1]),
        target_x, target_y, positions, cst, var_term, inv_L0, 0.5)
    covariance_lazy!(cov, @view(guide_x[:, :, 1]), @view(guide_y[:, :, 1]),
        target_x, target_y, positions, cst, var_term, inv_L0, 0.5)

    summarize("covariance_materialized_shifted_coords", @benchmark covariance_materialized!($cov, $iz, $jz,
        $(view(guide_x, :, :, 1)), $(view(guide_y, :, :, 1)), $target_x, $target_y,
        $positions, $cst, $var_term, $inv_L0, 0.5))
    summarize("covariance_lazy_shifted_coords", @benchmark covariance_lazy!($cov,
        $(view(guide_x, :, :, 1)), $(view(guide_y, :, :, 1)), $target_x, $target_y,
        $positions, $cst, $var_term, $inv_L0, 0.5))
end

main()
