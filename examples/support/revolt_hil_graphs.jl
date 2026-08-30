module REVOLTHILGraphs

export graph_path
export supported_architectures
export supported_profiles

const _GRAPH_DIRECTORY = normpath(joinpath(@__DIR__, "..", "graphs"))
const _SUPPORTED_ARCHITECTURES = (:classic, :copper)
const _SUPPORTED_PROFILES = (:full_optical, :fast_dm)

"""Return the REVOLT architectures with maintained external-RTC graphs."""
supported_architectures() = _SUPPORTED_ARCHITECTURES

"""Return the run-immutable profiles supported by the REVOLT graphs."""
supported_profiles() = _SUPPORTED_PROFILES

@inline _graph_filename(::Val{:classic}, ::Val{:full_optical}) =
    "revolt_classic_hil.toml"
@inline _graph_filename(::Val{:classic}, ::Val{:fast_dm}) =
    "revolt_classic_hil_fast_dm.toml"
@inline _graph_filename(::Val{:copper}, ::Val{:full_optical}) =
    "revolt_copper_hil.toml"
@inline _graph_filename(::Val{:copper}, ::Val{:fast_dm}) =
    "revolt_copper_hil_fast_dm.toml"

"""
    graph_path(architecture, profile=:full_optical)

Return the maintained REVOLT graph file for one architecture and run-immutable
profile. `:full_optical` selects the general coordinate-sampled Gaussian PDM
node. `:fast_dm` changes only that node to the separable regular-grid Gaussian
implementation; the remaining optical and detector configuration is identical.
"""
function graph_path(
    architecture::Symbol,
    profile::Symbol=:full_optical,
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
