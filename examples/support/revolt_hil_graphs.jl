module REVOLTHILGraphs

export graph_path
export supported_architectures
export supported_fidelities

const _GRAPH_DIRECTORY = normpath(joinpath(@__DIR__, "..", "graphs"))
const _SUPPORTED_ARCHITECTURES = (:classic, :copper)
const _SUPPORTED_FIDELITIES = (:full_optical, :reduced_resolution)

"""Return the REVOLT architectures with maintained external-RTC graphs."""
supported_architectures() = _SUPPORTED_ARCHITECTURES

"""Return the immutable fidelity choices supported by the REVOLT graphs."""
supported_fidelities() = _SUPPORTED_FIDELITIES

@inline _graph_filename(::Val{:classic}, ::Val{:full_optical}) =
    "revolt_classic_hil.toml"
@inline _graph_filename(::Val{:classic}, ::Val{:reduced_resolution}) =
    "revolt_classic_hil_reduced_resolution.toml"
@inline _graph_filename(::Val{:copper}, ::Val{:full_optical}) =
    "revolt_copper_hil.toml"
@inline _graph_filename(::Val{:copper}, ::Val{:reduced_resolution}) =
    "revolt_copper_hil_reduced_resolution.toml"

"""
    graph_path(architecture, fidelity=:full_optical)

Return the maintained REVOLT graph file for one architecture and run-immutable
fidelity choice. `:full_optical` selects the primary instrument-scale graph.
`:reduced_resolution` retains atmosphere evolution, physical-DM command
response, diffractive sensor formation, output frame shape, and HIL sequencing
while reducing internal optical sampling and stochastic detector work.
"""
function graph_path(
    architecture::Symbol,
    fidelity::Symbol=:full_optical,
)
    architecture in _SUPPORTED_ARCHITECTURES || throw(ArgumentError(
        "unsupported REVOLT architecture '$architecture'; expected one of " *
        "$(join(_SUPPORTED_ARCHITECTURES, ", "))",
    ))
    fidelity in _SUPPORTED_FIDELITIES || throw(ArgumentError(
        "unsupported REVOLT fidelity '$fidelity'; expected one of " *
        "$(join(_SUPPORTED_FIDELITIES, ", "))",
    ))
    filename = _graph_filename(Val(architecture), Val(fidelity))
    return joinpath(_GRAPH_DIRECTORY, filename)
end

end # module REVOLTHILGraphs
