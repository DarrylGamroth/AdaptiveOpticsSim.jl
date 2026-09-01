module REVOLTHILGraphs

export graph_path
export supported_architectures
export supported_profiles

const _GRAPH_DIRECTORY = normpath(joinpath(@__DIR__, "..", "graphs"))
const _SUPPORTED_ARCHITECTURES = (:classic, :copper)
const _SUPPORTED_PROFILES = (:coordinate_gaussian, :grid_gaussian)

"""Return the REVOLT architectures with maintained external-RTC graphs."""
supported_architectures() = _SUPPORTED_ARCHITECTURES

"""Return the run-immutable profiles supported by the REVOLT graphs."""
supported_profiles() = _SUPPORTED_PROFILES

@inline _graph_filename(::Val{:classic}, ::Val{:coordinate_gaussian}) =
    "revolt_classic_hil_coordinate_gaussian.toml"
@inline _graph_filename(::Val{:classic}, ::Val{:grid_gaussian}) =
    "revolt_classic_hil_grid_gaussian.toml"
@inline _graph_filename(::Val{:copper}, ::Val{:coordinate_gaussian}) =
    "revolt_copper_hil_coordinate_gaussian.toml"
@inline _graph_filename(::Val{:copper}, ::Val{:grid_gaussian}) =
    "revolt_copper_hil_grid_gaussian.toml"

"""
    graph_path(architecture, profile=:grid_gaussian)

Return the maintained REVOLT graph file for one architecture and run-immutable
profile. `:coordinate_gaussian` selects the general actuator-coordinate
Gaussian PDM node. `:grid_gaussian` selects the mathematically equivalent
separable regular-grid Gaussian implementation; the remaining optical and
detector configuration is identical.
"""
function graph_path(
    architecture::Symbol,
    profile::Symbol=:grid_gaussian,
)
    architecture in _SUPPORTED_ARCHITECTURES || throw(ArgumentError(
        "unsupported REVOLT architecture '$architecture'; expected one of " *
        "$(join(_SUPPORTED_ARCHITECTURES, ", "))",
    ))
    profile in _SUPPORTED_PROFILES || throw(ArgumentError(
        "unsupported REVOLT profile '$profile'; expected one of " *
        "$(join(_SUPPORTED_PROFILES, ", "))",
    ))
    filename = _graph_filename(Val(architecture), Val(profile))
    return joinpath(_GRAPH_DIRECTORY, filename)
end

end # module REVOLTHILGraphs
