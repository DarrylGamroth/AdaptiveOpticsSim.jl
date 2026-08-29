module AdaptiveOpticsSimCalculonAlgorithmsExt

import AdaptiveOpticsSim: AlgorithmGraphs, Calibration, Control
import CalculonAlgorithms

const AG = AlgorithmGraphs
const CA = CalculonAlgorithms

"""Calculon execution slot retaining separate AOS state and scratch-workspace owners."""
struct _StateWorkspace{State,Workspace}
    state::State
    workspace::Workspace
end

include("calculon_algorithms/control.jl")
include("calculon_algorithms/calibration.jl")

@inline _direction(::CA.InputPort) = :input
@inline _direction(::CA.OutputPort) = :output

@inline _role(::CA.FrameData) = :data
@inline _role(::CA.SparseParameter) = :parameter

@inline _layout(::CA.ColumnMajorLayout) = :column_major
@inline _layout(::CA.RowMajorLayout) = :row_major

function AlgorithmGraphs.prepare_algorithm_instance(
    ::Type{Declaration},
    configuration,
) where {Declaration<:CA.AbstractAlgorithmDeclaration}
    return CA.prepare_algorithm(Declaration, configuration)
end

function _port_contract(port, format)
    return AG.graph_port_contract(
        port.name,
        _direction(port.direction),
        _role(port.role),
        format.element_type,
        format.shape,
        format.schema,
        _layout(format.layout),
    )
end

function AlgorithmGraphs.algorithm_port_contracts(
    ::Type{Declaration},
    prepared::CA.PreparedAlgorithm{Declaration},
) where {Declaration<:CA.AbstractAlgorithmDeclaration}
    ports = CA.algorithm_ports(Declaration)
    formats = CA.algorithm_formats(prepared)
    length(ports) == length(formats) || throw(AG.AlgorithmGraphError(
        "Calculon declaration $Declaration returned inconsistent port and format counts",
    ))
    return map(_port_contract, ports, formats)
end

@inline function AlgorithmGraphs.process_algorithm!(
    prepared::CA.PreparedAlgorithm,
    outputs::NamedTuple,
    inputs::NamedTuple,
)
    CA.process!(values(outputs)..., prepared, values(inputs)...)
    return nothing
end

@inline function AlgorithmGraphs.reset_algorithm!(prepared::CA.PreparedAlgorithm)
    CA.reset!(prepared)
    return nothing
end

function AlgorithmGraphs.replace_algorithm_parameter!(
    prepared::CA.PreparedAlgorithm,
    name::Val,
    values,
)
    CA.replace_parameter!(prepared, name, values)
    return nothing
end

end # module AdaptiveOpticsSimCalculonAlgorithmsExt
