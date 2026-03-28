module SubaruAO3kSimulation

include(joinpath(dirname(@__FILE__), "subaru_ao188_simulation.jl"))

using AdaptiveOpticsSim
using .SubaruAO188Simulation

export AO3kNIRPyramidModel, AO3kSimulationParams, subaru_ao3k_simulation

struct AO3kNIRPyramidModel{T<:AbstractFloat} <: SubaruAO188Simulation.SubaruHighOrderWFSModel
    modulation::T
    modulation_points::Int
end

AO3kNIRPyramidModel(; modulation::Real=2.0, modulation_points::Int=8, T::Type{<:AbstractFloat}=Float32) =
    AO3kNIRPyramidModel{T}(T(modulation), modulation_points)

function AO3kSimulationParams(; kwargs...)
    nt = (; kwargs...)
    T0 = get(nt, :T, Float32)
    sampling = get(nt, :sampling_time, 1e-3)
    high_detector = get(nt, :high_detector, SubaruAO188Simulation.AO188WFSDetectorConfig(
        T=T0,
        integration_time=sampling,
        qe=0.9,
        psf_sampling=1,
        binning=1,
        gain=1.0,
        dark_current=0.0,
        noise=NoiseReadout(0.1),
        sensor=HgCdTeAvalancheArraySensor(read_time=sampling / 4, sampling_mode=CorrelatedDoubleSampling(), T=T0),
        correction_model=ReferencePixelCommonModeCorrection(4, 4),
    ))
    rest = Base.structdiff(nt, (; source_band=nothing, n_control_modes=nothing, n_subap=nothing,
        resolution=nothing, source_magnitude=nothing, high_order_sensor_model=nothing, high_detector=nothing))
    return SubaruAO188Simulation.AO188SimulationParams(
        ;
        source_band=:H,
        n_control_modes=get(nt, :n_control_modes, 1024),
        n_subap=get(nt, :n_subap, 32),
        resolution=get(nt, :resolution, 160),
        source_magnitude=get(nt, :source_magnitude, 10.0),
        high_order_sensor_model=AO3kNIRPyramidModel(T=T0),
        high_detector=high_detector,
        rest...,
    )
end

function SubaruAO188Simulation._build_high_order_wfs(model::AO3kNIRPyramidModel, tel::Telescope, params; backend=Array)
    T = eltype(tel.state.opd)
    return PyramidWFS(
        tel;
        n_subap=params.n_subap,
        mode=Diffractive(),
        modulation=model.modulation,
        modulation_points=model.modulation_points,
        T=T,
        backend=backend,
    )
end

function subaru_ao3k_simulation(;
    params::SubaruAO188Simulation.AO188SimulationParams=AO3kSimulationParams(),
    kwargs...,
)
    return SubaruAO188Simulation.subaru_ao188_simulation(; params=params, kwargs...)
end

end # module SubaruAO3kSimulation
