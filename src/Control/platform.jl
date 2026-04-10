#
# Julia-native platform orchestration layer
#
# This file introduces typed scenario/config objects for runtime composition
# without switching the package to config-file-driven control.
#

abstract type AbstractPlatformConfig end

struct ClosedLoopBranchConfig{SIM,REC,WD,SD,RNG,PR,CL}
    label::Symbol
    simulation::SIM
    reconstructor::REC
    wfs_detector::WD
    science_detector::SD
    rng::RNG
    products::PR
    command_layout::CL
end

function ClosedLoopBranchConfig(label::Symbol, simulation, reconstructor;
    wfs_detector=nothing,
    science_detector=nothing,
    rng::AbstractRNG=Random.default_rng(),
    products=nothing,
    command_layout=nothing)
    return ClosedLoopBranchConfig{
        typeof(simulation),
        typeof(reconstructor),
        typeof(wfs_detector),
        typeof(science_detector),
        typeof(rng),
        typeof(products),
        typeof(command_layout),
    }(
        label,
        simulation,
        reconstructor,
        wfs_detector,
        science_detector,
        rng,
        products,
        command_layout,
    )
end

struct SinglePlatformConfig{RP<:AbstractRuntimeProfile,PR<:RuntimeProductRequirements,T<:AbstractFloat} <: AbstractPlatformConfig
    name::Symbol
    branch_label::Symbol
    profile::RP
    products::PR
    latency::RuntimeLatencyModel
    control_sign::T
    science_zero_padding::Int
end

function SinglePlatformConfig(;
    name::Symbol=:platform_runtime,
    branch_label::Symbol=:main,
    profile::AbstractRuntimeProfile=default_runtime_profile(),
    products::RuntimeProductRequirements=default_runtime_products(),
    latency::RuntimeLatencyModel=default_runtime_latency(profile),
    control_sign::Real=1.0,
    science_zero_padding::Integer=runtime_science_zero_padding(profile),
)
    T = typeof(float(control_sign))
    return SinglePlatformConfig{
        typeof(profile),
        typeof(products),
        T,
    }(
        name,
        branch_label,
        profile,
        products,
        latency,
        T(control_sign),
        Int(science_zero_padding),
    )
end

struct GroupedPlatformConfig{N,RP<:AbstractRuntimeProfile,PR<:GroupedRuntimeProductRequirements,T<:AbstractFloat} <: AbstractPlatformConfig
    name::Symbol
    branch_labels::NTuple{N,Symbol}
    profile::RP
    products::PR
    latency::RuntimeLatencyModel
    control_sign::T
    science_zero_padding::Int
end

function GroupedPlatformConfig(branch_labels::NTuple{N,Symbol};
    name::Symbol=:grouped_platform_runtime,
    profile::AbstractRuntimeProfile=default_runtime_profile(),
    products::GroupedRuntimeProductRequirements=default_grouped_runtime_products(),
    latency::RuntimeLatencyModel=default_runtime_latency(profile),
    control_sign::Real=1.0,
    science_zero_padding::Integer=runtime_science_zero_padding(profile),
) where {N}
    T = typeof(float(control_sign))
    return GroupedPlatformConfig{
        N,
        typeof(profile),
        typeof(products),
        T,
    }(
        name,
        branch_labels,
        profile,
        products,
        latency,
        T(control_sign),
        Int(science_zero_padding),
    )
end

function GroupedPlatformConfig(branches::Vararg{ClosedLoopBranchConfig,N}; kwargs...) where {N}
    return GroupedPlatformConfig(ntuple(i -> branches[i].label, N); kwargs...)
end

struct PlatformScenario{CFG<:AbstractPlatformConfig,BR,BD} <: AbstractControlSimulation
    config::CFG
    branches::BR
    boundary::BD
end

const RuntimeBranch = ClosedLoopBranchConfig
const SingleRuntimeConfig = SinglePlatformConfig
const GroupedRuntimeConfig = GroupedPlatformConfig
const RuntimeScenario = PlatformScenario

@inline platform_config(scenario::PlatformScenario) = scenario.config
@inline platform_boundary(scenario::PlatformScenario) = scenario.boundary
@inline platform_name(config::SinglePlatformConfig) = config.name
@inline platform_name(config::GroupedPlatformConfig) = config.name
@inline platform_name(scenario::PlatformScenario) = platform_name(platform_config(scenario))
@inline platform_branch_labels(config::SinglePlatformConfig) = (config.branch_label,)
@inline platform_branch_labels(config::GroupedPlatformConfig) = config.branch_labels
@inline platform_branch_labels(scenario::PlatformScenario) = platform_branch_labels(platform_config(scenario))
@inline runtime_profile(scenario::PlatformScenario) = platform_config(scenario).profile
@inline runtime_latency(scenario::PlatformScenario) = platform_config(scenario).latency

@inline function _branch_runtime_products(config::SinglePlatformConfig, branch::ClosedLoopBranchConfig)
    isnothing(branch.products) || return branch.products
    return config.products
end

@inline function _branch_runtime_products(config::GroupedPlatformConfig, branch::ClosedLoopBranchConfig)
    isnothing(branch.products) || return branch.products
    return RuntimeProductRequirements(
        slopes=true,
        wfs_pixels=(config.products.wfs_frames || config.products.wfs_stack) && !isnothing(branch.wfs_detector),
        science_pixels=(config.products.science_frames || config.products.science_stack) && !isnothing(branch.science_detector),
    )
end

@inline function _build_platform_runtime(config::AbstractPlatformConfig, branch::ClosedLoopBranchConfig)
    return ClosedLoopRuntime(
        branch.simulation,
        branch.reconstructor;
        wfs_detector=branch.wfs_detector,
        science_detector=branch.science_detector,
        rng=branch.rng,
        profile=config.profile,
        products=_branch_runtime_products(config, branch),
        latency=config.latency,
        control_sign=config.control_sign,
        science_zero_padding=config.science_zero_padding,
        command_layout=branch.command_layout,
    )
end

function build_platform_scenario(config::SinglePlatformConfig, branch::ClosedLoopBranchConfig)
    branch.label == config.branch_label ||
        throw(InvalidConfiguration("single platform branch label $(branch.label) must match config branch_label $(config.branch_label)"))
    runtime = _build_platform_runtime(config, branch)
    boundary = SimulationInterface(runtime)
    return PlatformScenario{
        typeof(config),
        Tuple{typeof(branch)},
        typeof(boundary),
    }(
        config,
        (branch,),
        boundary,
    )
end

function build_platform_scenario(config::GroupedPlatformConfig{N}, branches::Vararg{ClosedLoopBranchConfig,N}) where {N}
    labels = ntuple(i -> branches[i].label, N)
    labels == config.branch_labels ||
        throw(InvalidConfiguration("grouped platform branch labels $labels must match config branch_labels $(config.branch_labels)"))
    runtimes = ntuple(i -> _build_platform_runtime(config, branches[i]), N)
    boundary = CompositeSimulationInterface(runtimes...; products=config.products)
    return PlatformScenario{
        typeof(config),
        typeof(branches),
        typeof(boundary),
    }(
        config,
        branches,
        boundary,
    )
end

build_runtime_scenario(config::SingleRuntimeConfig, branch::RuntimeBranch) = build_platform_scenario(config, branch)
build_runtime_scenario(config::GroupedRuntimeConfig{N}, branches::Vararg{RuntimeBranch,N}) where {N} = build_platform_scenario(config, branches...)

simulation_interface(scenario::PlatformScenario) = platform_boundary(scenario)

readout(sim::AbstractControlSimulation) = simulation_readout(sim)
readout(runtime::ClosedLoopRuntime) = simulation_readout(runtime)
readout(interface::SimulationInterface) = simulation_readout(interface)
readout(interface::CompositeSimulationInterface) = simulation_readout(interface)

function prepare!(scenario::PlatformScenario)
    prepare!(platform_boundary(scenario))
    return scenario
end

function sense!(scenario::PlatformScenario)
    sense!(platform_boundary(scenario))
    return scenario
end

function step!(scenario::PlatformScenario)
    step!(platform_boundary(scenario))
    return scenario
end

function snapshot_outputs!(scenario::PlatformScenario)
    snapshot_outputs!(platform_boundary(scenario))
    return scenario
end

@inline set_command!(scenario::PlatformScenario, command::AbstractVector) = set_command!(platform_boundary(scenario), command)
@inline runtime_phase_timing(scenario::PlatformScenario; kwargs...) = runtime_phase_timing(platform_boundary(scenario); kwargs...)

function command_layout(scenario::PlatformScenario)
    labels = platform_branch_labels(scenario)
    boundary = platform_boundary(scenario)
    if boundary isa SimulationInterface
        return command_layout(boundary)
    end
    segments = ntuple(i -> RuntimeCommandSegment(labels[i], i == 1 ? 1 : 1 + sum(length(boundary.interfaces[j].command) for j in 1:(i - 1)),
        length(boundary.interfaces[i].command)), length(labels))
    return RuntimeCommandLayout(segments)
end

@inline function branch_command_layout(scenario::PlatformScenario, label::Symbol)
    boundary = platform_boundary(scenario)
    labels = platform_branch_labels(scenario)
    idx = findfirst(==(label), labels)
    isnothing(idx) && throw(InvalidConfiguration("runtime scenario does not contain branch $(label)"))
    if boundary isa SimulationInterface
        return command_layout(boundary)
    end
    return command_layout(boundary.interfaces[idx])
end

@inline function branch_command_layouts(scenario::PlatformScenario)
    labels = platform_branch_labels(scenario)
    layouts = map(label -> branch_command_layout(scenario, label), labels)
    return NamedTuple{labels}(Tuple(layouts))
end

@inline function update_command!(scenario::PlatformScenario, command::NamedTuple)
    boundary = platform_boundary(scenario)
    labels = platform_branch_labels(scenario)
    if boundary isa SimulationInterface
        update_command!(boundary, command)
    else
        all(label -> label in labels, Tuple(keys(command))) ||
            throw(InvalidConfiguration("structured scenario command keys $(Tuple(keys(command))) must be a subset of branch labels $labels"))
        @inbounds for (idx, label) in pairs(labels)
            if hasproperty(command, label)
                branch_command = getproperty(command, label)
                if branch_command isa NamedTuple
                    update_command!(boundary.interfaces[idx], branch_command)
                elseif branch_command isa AbstractVector
                    set_command!(boundary.interfaces[idx], branch_command)
                else
                    throw(InvalidConfiguration("structured branch command $(label) must be an AbstractVector or NamedTuple"))
                end
            end
        end
    end
    snapshot_outputs!(scenario)
    return simulation_command(scenario)
end

@inline function set_command!(scenario::PlatformScenario, command::NamedTuple)
    boundary = platform_boundary(scenario)
    labels = platform_branch_labels(scenario)
    if boundary isa SimulationInterface
        set_command!(boundary, command)
    else
        Tuple(keys(command)) == labels ||
            throw(InvalidConfiguration("structured scenario command keys $(Tuple(keys(command))) must match branch labels $labels"))
        @inbounds for (idx, label) in pairs(labels)
            branch_command = getproperty(command, label)
            if branch_command isa NamedTuple
                set_command!(boundary.interfaces[idx], branch_command)
            elseif branch_command isa AbstractVector
                set_command!(boundary.interfaces[idx], branch_command)
            else
                throw(InvalidConfiguration("structured branch command $(label) must be an AbstractVector or NamedTuple"))
            end
        end
    end
    snapshot_outputs!(scenario)
    return simulation_command(scenario)
end

supports_prepared_runtime(scenario::PlatformScenario) =
    any(branch -> supports_prepared_runtime(branch.simulation.wfs, branch.simulation.src), scenario.branches)
supports_detector_output(scenario::PlatformScenario) = any(branch -> !isnothing(branch.wfs_detector) || !isnothing(branch.science_detector), scenario.branches)
supports_grouped_execution(scenario::PlatformScenario) = length(scenario.branches) > 1
