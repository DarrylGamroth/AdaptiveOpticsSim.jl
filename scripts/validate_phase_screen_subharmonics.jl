using AdaptiveOpticsSim
using LinearAlgebra
using Random
using Statistics

function subharmonic_theory_variance(r0::Real, L0::Real, D::Real, l0::Real;
    n_levels::Int, radius::Int)
    fm = 5.92 / l0 / (2 * pi)
    f0 = 1 / L0
    total = 0.0
    for p in 1:n_levels
        del_f = 1 / (D * 3^p)
        for fx in -radius:radius, fy in -radius:radius
            if fx == 0 && fy == 0
                continue
            end
            f = sqrt((fx * del_f)^2 + (fy * del_f)^2)
            psd = 0.023 * r0^(-5 / 3) * exp(-((f / fm)^2)) / ((f^2 + f0^2)^(11 / 6))
            total += psd * del_f^2
        end
    end
    return total
end

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
    added_variance = Float64[]
    for seed in 1:nsamp
        base = ft_sh_phase_screen(atm, n, delta;
            rng=MersenneTwister(seed), ws=ws, subharmonics=false)
        phs = ft_sh_phase_screen(atm, n, delta; rng=MersenneTwister(seed), ws=ws, kwargs...)
        sub = phs .- base
        push!(tiptilt, tiptilt_power(phs))
        push!(structure, structure_energy(phs, shift))
        push!(added_variance, mean(abs2, sub))
    end
    return (; tiptilt=mean(tiptilt), structure=mean(structure), added_variance=mean(added_variance))
end

function main()
    tel = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    println("L0\tmode\ttiptilt_power\tstructure_energy\tadded_variance\ttheory_variance\tratio")
    for L0 in (25.0, 50.0, 100.0, 200.0)
        atm = KolmogorovAtmosphere(tel; r0=0.2, L0=L0)
        for (label, kwargs, theory) in (
            ("none", (subharmonics=false,), nothing),
            ("fast", (subharmonics=true, mode=FastSubharmonics()),
                subharmonic_theory_variance(atm.params.r0, atm.params.L0, tel.params.diameter, 1e-10;
                    n_levels=3, radius=1)),
            ("fidelity", (subharmonics=true, mode=FidelitySubharmonics()),
                subharmonic_theory_variance(atm.params.r0, atm.params.L0, tel.params.diameter, 1e-10;
                    n_levels=AdaptiveOpticsSim.resolve_subharmonic_levels(atm.params.L0, tel.params.diameter), radius=2)),
        )
            m = subharmonic_metrics(atm, tel.params.diameter; kwargs...)
            ratio = theory === nothing ? "NA" : string(m.added_variance / theory)
            theory_str = theory === nothing ? "NA" : string(theory)
            println("$(L0)\t$(label)\t$(m.tiptilt)\t$(m.structure)\t$(m.added_variance)\t$(theory_str)\t$(ratio)")
        end
    end
end

main()
