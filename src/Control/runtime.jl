mutable struct ClosedLoopRuntime{SIM<:AOSimulation,R,V,WD,SD,RNG,T<:AbstractFloat}
    simulation::SIM
    reconstructor::R
    command::V
    wfs_detector::WD
    science_detector::SD
    rng::RNG
    control_sign::T
    science_zero_padding::Int
end

function ClosedLoopRuntime(simulation::AOSimulation, reconstructor;
    wfs_detector=nothing, science_detector=nothing, rng=MersenneTwister(0),
    control_sign::Real=-1.0, science_zero_padding::Int=2)
    T = eltype(simulation.dm.state.coefs)
    command = similar(simulation.dm.state.coefs)
    fill!(command, zero(T))
    return ClosedLoopRuntime{typeof(simulation), typeof(reconstructor), typeof(command),
        typeof(wfs_detector), typeof(science_detector), typeof(rng), T}(
        simulation,
        reconstructor,
        command,
        wfs_detector,
        science_detector,
        rng,
        T(control_sign),
        science_zero_padding,
    )
end

@inline function set_command!(runtime::ClosedLoopRuntime, command::AbstractVector)
    copyto!(runtime.command, command)
    copyto!(runtime.simulation.dm.state.coefs, command)
    return runtime.simulation.dm.state.coefs
end

function sense!(runtime::ClosedLoopRuntime)
    sim = runtime.simulation
    tel = sim.tel
    atm = sim.atm
    dm = sim.dm
    wfs = sim.wfs
    src = sim.src
    advance!(atm, tel, runtime.rng)
    propagate!(atm, tel)
    apply!(dm, tel, DMAdditive())
    if runtime.wfs_detector === nothing
        measure!(wfs, tel, src)
    else
        measure!(wfs, tel, src, runtime.wfs_detector; rng=runtime.rng)
    end
    if !(runtime.science_detector === nothing)
        psf = compute_psf!(tel, src; zero_padding=runtime.science_zero_padding)
        capture!(runtime.science_detector, psf; rng=runtime.rng)
    end
    return runtime
end

function step!(runtime::ClosedLoopRuntime)
    sim = runtime.simulation
    sense!(runtime)
    reconstruct!(runtime.command, runtime.reconstructor, sim.wfs.state.slopes)
    copyto!(sim.dm.state.coefs, runtime.command)
    runtime.control_sign == one(runtime.control_sign) || rmul!(sim.dm.state.coefs, runtime.control_sign)
    return runtime
end
