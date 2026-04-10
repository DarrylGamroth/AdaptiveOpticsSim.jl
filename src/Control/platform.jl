#
# Julia-native platform orchestration layer
#
# This file introduces typed scenario/config objects for runtime composition
# without switching the package to config-file-driven control.
#

abstract type AbstractPlatformConfig end

struct RuntimeBranch{SIM,REC,WD,SD,RNG,PR,CL,B<:AbstractArrayBackend}
    label::Symbol
    simulation::SIM
    reconstructor::REC
    wfs_detector::WD
    science_detector::SD
    rng::RNG
    products::PR
    command_layout::CL
end

function RuntimeBranch(label::Symbol, simulation, reconstructor;
    wfs_detector=nothing,
    science_detector=nothing,
    rng::AbstractRNG=Random.default_rng(),
    products=nothing,
    command_layout=nothing)
    selector = require_same_backend(simulation, wfs_detector, science_detector)
    return RuntimeBranch{
        typeof(simulation),
        typeof(reconstructor),
        typeof(wfs_detector),
        typeof(science_detector),
        typeof(rng),
        typeof(products),
        typeof(command_layout),
        typeof(selector),
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

@inline backend(::RuntimeBranch{<:Any,<:Any,<:Any,<:Any,<:Any,<:Any,<:Any,B}) where {B} = B()

struct SingleRuntimeConfig{RP<:AbstractRuntimeProfile,PR<:RuntimeProductRequirements,T<:AbstractFloat} <: AbstractPlatformConfig
    name::Symbol
    branch_label::Symbol
    profile::RP
    products::PR
    latency::RuntimeLatencyModel
    control_sign::T
    science_zero_padding::Int
end

function SingleRuntimeConfig(;
    name::Symbol=:platform_runtime,
    branch_label::Symbol=:main,
    profile::AbstractRuntimeProfile=default_runtime_profile(),
    products::RuntimeProductRequirements=default_runtime_products(),
    latency::RuntimeLatencyModel=default_runtime_latency(profile),
    control_sign::Real=1.0,
    science_zero_padding::Integer=runtime_science_zero_padding(profile),
)
    T = typeof(float(control_sign))
    return SingleRuntimeConfig{
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

struct GroupedRuntimeConfig{N,RP<:AbstractRuntimeProfile,PR<:GroupedRuntimeProductRequirements,T<:AbstractFloat} <: AbstractPlatformConfig
    name::Symbol
    branch_labels::NTuple{N,Symbol}
    profile::RP
    products::PR
    latency::RuntimeLatencyModel
    control_sign::T
    science_zero_padding::Int
end

function GroupedRuntimeConfig(branch_labels::NTuple{N,Symbol};
    name::Symbol=:grouped_platform_runtime,
    profile::AbstractRuntimeProfile=default_runtime_profile(),
    products::GroupedRuntimeProductRequirements=default_grouped_runtime_products(),
    latency::RuntimeLatencyModel=default_runtime_latency(profile),
    control_sign::Real=1.0,
    science_zero_padding::Integer=runtime_science_zero_padding(profile),
) where {N}
    T = typeof(float(control_sign))
    return GroupedRuntimeConfig{
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

function GroupedRuntimeConfig(branches::Vararg{RuntimeBranch,N}; kwargs...) where {N}
    return GroupedRuntimeConfig(ntuple(i -> branches[i].label, N); kwargs...)
end

struct RuntimeScenario{CFG<:AbstractPlatformConfig,BR,BD,B<:AbstractArrayBackend} <: AbstractControlSimulation
    config::CFG
    branches::BR
    boundary::BD
end


@inline platform_config(scenario::RuntimeScenario) = scenario.config
@inline platform_boundary(scenario::RuntimeScenario) = scenario.boundary
@inline platform_name(config::SingleRuntimeConfig) = config.name
@inline platform_name(config::GroupedRuntimeConfig) = config.name
@inline platform_name(scenario::RuntimeScenario) = platform_name(platform_config(scenario))
@inline platform_branch_labels(config::SingleRuntimeConfig) = (config.branch_label,)
@inline platform_branch_labels(config::GroupedRuntimeConfig) = config.branch_labels
@inline platform_branch_labels(scenario::RuntimeScenario) = platform_branch_labels(platform_config(scenario))
@inline runtime_profile(scenario::RuntimeScenario) = platform_config(scenario).profile
@inline runtime_latency(scenario::RuntimeScenario) = platform_config(scenario).latency

@inline function _branch_runtime_products(config::SingleRuntimeConfig, branch::RuntimeBranch)
    isnothing(branch.products) || return branch.products
    return config.products
end

@inline function _branch_runtime_products(config::GroupedRuntimeConfig, branch::RuntimeBranch)
    isnothing(branch.products) || return branch.products
    return RuntimeProductRequirements(
        slopes=true,
        wfs_pixels=(config.products.wfs_frames || config.products.wfs_stack) && !isnothing(branch.wfs_detector),
        science_pixels=(config.products.science_frames || config.products.science_stack) && !isnothing(branch.science_detector),
    )
end

@inline function _build_platform_runtime(config::AbstractPlatformConfig, branch::RuntimeBranch)
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

function build_runtime_scenario(config::SingleRuntimeConfig, branch::RuntimeBranch)
    branch.label == config.branch_label ||
        throw(InvalidConfiguration("single platform branch label $(branch.label) must match config branch_label $(config.branch_label)"))
    runtime = _build_platform_runtime(config, branch)
    boundary = SimulationInterface(runtime)
    selector = backend(branch)
    return RuntimeScenario{
        typeof(config),
        Tuple{typeof(branch)},
        typeof(boundary),
        typeof(selector),
    }(
        config,
        (branch,),
        boundary,
    )
end

function build_runtime_scenario(config::GroupedRuntimeConfig{N}, branches::Vararg{RuntimeBranch,N}) where {N}
    selector = require_same_backend(branches...)
    labels = ntuple(i -> branches[i].label, N)
    labels == config.branch_labels ||
        throw(InvalidConfiguration("grouped platform branch labels $labels must match config branch_labels $(config.branch_labels)"))
    runtimes = ntuple(i -> _build_platform_runtime(config, branches[i]), N)
    boundary = CompositeSimulationInterface(runtimes...; products=config.products)
    return RuntimeScenario{
        typeof(config),
        typeof(branches),
        typeof(boundary),
        typeof(selector),
    }(
        config,
        branches,
        boundary,
    )
end

@inline backend(::RuntimeScenario{<:Any,<:Any,<:Any,B}) where {B} = B()
simulation_interface(scenario::RuntimeScenario) = platform_boundary(scenario)

readout(sim::AbstractControlSimulation) = simulation_readout(sim)
readout(runtime::ClosedLoopRuntime) = simulation_readout(runtime)
readout(interface::SimulationInterface) = simulation_readout(interface)
readout(interface::CompositeSimulationInterface) = simulation_readout(interface)

function prepare!(scenario::RuntimeScenario)
    prepare!(platform_boundary(scenario))
    return scenario
end

function sense!(scenario::RuntimeScenario)
    sense!(platform_boundary(scenario))
    return scenario
end

function step!(scenario::RuntimeScenario)
    step!(platform_boundary(scenario))
    return scenario
end

function snapshot_outputs!(scenario::RuntimeScenario)
    snapshot_outputs!(platform_boundary(scenario))
    return scenario
end

@inline set_command!(scenario::RuntimeScenario, command::AbstractVector) = set_command!(platform_boundary(scenario), command)
@inline runtime_phase_timing(scenario::RuntimeScenario; kwargs...) = runtime_phase_timing(platform_boundary(scenario); kwargs...)

function command_layout(scenario::RuntimeScenario)
    labels = platform_branch_labels(scenario)
    boundary = platform_boundary(scenario)
    if boundary isa SimulationInterface
        return command_layout(boundary)
    end
    segments = ntuple(i -> RuntimeCommandSegment(labels[i], i == 1 ? 1 : 1 + sum(length(boundary.interfaces[j].command) for j in 1:(i - 1)),
        length(boundary.interfaces[i].command)), length(labels))
    return RuntimeCommandLayout(segments)
end

@inline function branch_command_layout(scenario::RuntimeScenario, label::Symbol)
    boundary = platform_boundary(scenario)
    labels = platform_branch_labels(scenario)
    idx = findfirst(==(label), labels)
    isnothing(idx) && throw(InvalidConfiguration("runtime scenario does not contain branch $(label)"))
    if boundary isa SimulationInterface
        return command_layout(boundary)
    end
    return command_layout(boundary.interfaces[idx])
end

@inline function branch_command_layouts(scenario::RuntimeScenario)
    labels = platform_branch_labels(scenario)
    layouts = map(label -> branch_command_layout(scenario, label), labels)
    return NamedTuple{labels}(Tuple(layouts))
end

@inline function update_command!(scenario::RuntimeScenario, command::NamedTuple)
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

@inline function set_command!(scenario::RuntimeScenario, command::NamedTuple)
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

supports_prepared_runtime(scenario::RuntimeScenario) =
    any(branch -> supports_prepared_runtime(branch.simulation.wfs, branch.simulation.src), scenario.branches)
supports_detector_output(scenario::RuntimeScenario) = any(branch -> !isnothing(branch.wfs_detector) || !isnothing(branch.science_detector), scenario.branches)
supports_grouped_execution(scenario::RuntimeScenario) = length(scenario.branches) > 1
