module AdaptiveOpticsSimCalculonAlgorithmsExt

import AdaptiveOpticsSim: AlgorithmGraphs
import CalculonAlgorithms

const AG = AlgorithmGraphs
const CA = CalculonAlgorithms

@inline _direction(::CA.InputPort) = :input
@inline _direction(::CA.OutputPort) = :output

@inline _role(::CA.FrameData) = :data
@inline _role(::CA.SparseParameter) = :parameter

@inline _layout(::CA.ColumnMajorLayout) = :column_major
@inline _layout(::CA.RowMajorLayout) = :row_major

function AlgorithmGraphs._prepare_algorithm_instance(
    ::Type{Declaration},
    configuration,
) where {Declaration<:CA.AbstractAlgorithmDeclaration}
    return CA.prepare_algorithm(Declaration, configuration)
end

function _port_contract(port, format)
    return AG._graph_port_contract(
        port.name,
        _direction(port.direction),
        _role(port.role),
        format.element_type,
        format.shape,
        format.schema,
        _layout(format.layout),
    )
end

function AlgorithmGraphs._algorithm_port_contracts(
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

@inline function AlgorithmGraphs._process_algorithm!(
    prepared::CA.PreparedAlgorithm,
    outputs::NamedTuple,
    inputs::NamedTuple,
)
    CA.process!(values(outputs)..., prepared, values(inputs)...)
    return nothing
end

@inline function AlgorithmGraphs._reset_algorithm!(prepared::CA.PreparedAlgorithm)
    CA.reset!(prepared)
    return nothing
end

function AlgorithmGraphs._replace_algorithm_parameter!(
    prepared::CA.PreparedAlgorithm,
    name::Val,
    values,
)
    CA.replace_parameter!(prepared, name, values)
    return nothing
end

end # module AdaptiveOpticsSimCalculonAlgorithmsExt
