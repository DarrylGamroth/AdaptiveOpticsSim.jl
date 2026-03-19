using AdaptiveOpticsSim
using LinearAlgebra
using Random
using Statistics

function tiptilt_power(phs::AbstractMatrix)
    n = size(phs, 1)
    coords = collect(LinRange(-1.0, 1.0, n))
    x = repeat(reshape(coords, 1, n), n, 1)
    y = repeat(reshape(coords, n, 1), 1, n)
    A = hcat(vec(x), vec(y), ones(n * n))
    coeffs = A \ vec(phs)
    return coeffs[1]^2 + coeffs[2]^2
end

function structure_energy(phs::AbstractMatrix, shift::Int)
    diff = @views phs[:, 1:end-shift] .- phs[:, 1+shift:end]
    return mean(abs2, diff)
end

function subharmonic_metrics(atm::KolmogorovAtmosphere, D::Real;
    n::Int=32, nsamp::Int=24, kwargs...)
    delta = D / n
    shift = n ÷ 2
    ws = PhaseStatsWorkspace(n; T=Float64)
    tiptilt = Float64[]
    structure = Float64[]
    for seed in 1:nsamp
        phs = ft_sh_phase_screen(atm, n, delta; rng=MersenneTwister(seed), ws=ws, kwargs...)
        push!(tiptilt, tiptilt_power(phs))
        push!(structure, structure_energy(phs, shift))
    end
    return (; tiptilt=mean(tiptilt), structure=mean(structure))
end

function main()
    tel = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    println("L0\tmode\ttiptilt_power\tstructure_energy")
    for L0 in (25.0, 50.0, 100.0, 200.0)
        atm = KolmogorovAtmosphere(tel; r0=0.2, L0=L0)
        for (label, kwargs) in (
            ("none", (subharmonics=false,)),
            ("legacy", (subharmonics=true, n_levels=3, subharmonic_radius=1)),
            ("default", (subharmonics=true,)),
        )
            m = subharmonic_metrics(atm, tel.params.diameter; kwargs...)
            println("$(L0)\t$(label)\t$(m.tiptilt)\t$(m.structure)")
        end
    end
end

main()
