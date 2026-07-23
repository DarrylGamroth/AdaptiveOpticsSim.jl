#
# Prepared controller-output routing
#
# A controller may own one flat output buffer while the plant owns several
# independently timed command endpoints.  Routing binds named, caller-owned
# controller products (including views into a larger buffer) to exact prepared
# endpoints.  It deliberately owns no sequence, timestamp, admission,
# transaction, or transport semantics.
#

"""
    ControllerOutputRoute(product, endpoint)

Declare that the named controller output `product` supplies payload values for
one independently identified command `endpoint`.  The declaration contains no
storage selection, packed offset, command timing, or atomicity.
"""
struct ControllerOutputRoute
    product::Symbol
    endpoint::CommandEndpointID

    function ControllerOutputRoute(product::Symbol,
        endpoint::CommandEndpointID)
        isempty(String(product)) && throw(PlantPreparationError(
            :controller_output_routing, :empty_product,
            "controller-output product identity must not be empty"))
        return new(product, endpoint)
    end
end

ControllerOutputRoute(product::Symbol, endpoint) =
    ControllerOutputRoute(product, _as_command_endpoint_id(endpoint))

function ControllerOutputRoute(product, endpoint)
    throw(PlantPreparationError(:controller_output_routing,
        :invalid_product,
        "controller-output product identity must be a Symbol; got " *
        "$(typeof(product))"))
end

@inline controller_output_product(route::ControllerOutputRoute) =
    route.product
@inline controller_output_endpoint(route::ControllerOutputRoute) =
    route.endpoint

struct _PreparedControllerOutputRouteToken end
const _PREPARED_CONTROLLER_OUTPUT_ROUTE_TOKEN =
    _PreparedControllerOutputRouteToken()

"""
One validated, zero-copy controller-product binding.

The product is borrowed caller storage. Array products may be views into a
larger controller output, but their exact element type, shape, backend, and
physical device match the prepared endpoint. Scalar products use an assigned
`Ref{T}` and are host-resident.
"""
struct PreparedControllerOutputRoute{
    N,
    E<:PreparedCommandEndpoint,
    P,
}
    product::Symbol
    endpoint::E
    output::P

    function PreparedControllerOutputRoute(
        ::_PreparedControllerOutputRouteToken,
        product::Symbol,
        endpoint::E,
        output::P,
    ) where {E<:PreparedCommandEndpoint,P}
        endpoint_name = command_endpoint_id(endpoint).name
        return new{endpoint_name,E,P}(product, endpoint, output)
    end
end

struct _PreparedControllerOutputRoutingToken end
const _PREPARED_CONTROLLER_OUTPUT_ROUTING_TOKEN =
    _PreparedControllerOutputRoutingToken()

"""
    PreparedControllerOutputRouting

Canonical tuple of exact controller-product-to-endpoint bindings. The plan
borrows the declared controller products and performs no per-sample packing or
copy. Plant command admission remains the operation that copies a presented
payload into endpoint-owned bounded storage.
"""
struct PreparedControllerOutputRouting{R<:Tuple}
    routes::R

    function PreparedControllerOutputRouting(
        ::_PreparedControllerOutputRoutingToken,
        routes::R,
    ) where {R<:Tuple}
        return new{R}(routes)
    end
end

@inline controller_output_product(route::PreparedControllerOutputRoute) =
    route.product
@inline controller_output_endpoint(route::PreparedControllerOutputRoute) =
    command_endpoint_id(route.endpoint)
@inline controller_output_schema(route::PreparedControllerOutputRoute) =
    command_schema(route.endpoint)
@inline prepared_controller_output_routes(
    routing::PreparedControllerOutputRouting) = routing.routes

@inline controller_output_payload(
    route::PreparedControllerOutputRoute{N,E,P},
) where {N,E,P<:AbstractArray} =
    route.output
@inline controller_output_payload(
    route::PreparedControllerOutputRoute{N,E,P},
) where {N,E,P<:Base.RefValue} =
    route.output[]

function _prepared_command_endpoint_binding(plant::PreparedPlant,
    id::CommandEndpointID)
    @inbounds for binding in plant.command_endpoints
        command_endpoint_id(binding) == id && return binding
    end
    throw(PlantPreparationError(:controller_output_routing,
        :unknown_endpoint,
        "prepared plant has no command endpoint $id"))
end

function _require_controller_output_product(
    binding::_PreparedPlantCommandEndpoint{
        <:PreparedCommandEndpoint{<:PlantCommandSchema{T,0}}},
    product::Base.RefValue,
    name::Symbol,
) where {T}
    isassigned(product) || throw(PlantPreparationError(
        :controller_output_routing, :unassigned_scalar_product,
        "controller-output product $(repr(name)) must be assigned before " *
        "routing preparation"))
    typeof(product[]) === T || throw(PlantPreparationError(
        :controller_output_routing, :element_type,
        "controller-output product $(repr(name)) scalar type " *
        "$(typeof(product[])) does not match endpoint numeric type $T"))
    return product
end

function _require_controller_output_product(
    ::_PreparedPlantCommandEndpoint{
        <:PreparedCommandEndpoint{<:PlantCommandSchema{T,0}}},
    product::AbstractArray,
    name::Symbol,
) where {T}
    throw(PlantPreparationError(:controller_output_routing,
        :scalar_storage,
        "scalar controller-output product $(repr(name)) must use an " *
        "assigned Ref{$T}; got $(typeof(product))"))
end

function _require_controller_output_product(
    ::_PreparedPlantCommandEndpoint{
        <:PreparedCommandEndpoint{<:PlantCommandSchema{T,0}}},
    product,
    name::Symbol,
) where {T}
    throw(PlantPreparationError(:controller_output_routing,
        :scalar_storage,
        "scalar controller-output product $(repr(name)) must use an " *
        "assigned Ref{$T}; got $(typeof(product))"))
end

function _require_controller_output_product(
    binding::_PreparedPlantCommandEndpoint{
        <:PreparedCommandEndpoint{<:PlantCommandSchema{T,N}}},
    product::AbstractArray,
    name::Symbol,
) where {T,N}
    schema = command_schema(binding.endpoint)
    eltype(product) === T || throw(PlantPreparationError(
        :controller_output_routing, :element_type,
        "controller-output product $(repr(name)) element type " *
        "$(eltype(product)) does not match endpoint numeric type $T"))
    size(product) == command_dimensions(schema) || throw(
        PlantPreparationError(:controller_output_routing, :shape,
            "controller-output product $(repr(name)) shape " *
            "$(size(product)) does not match endpoint shape " *
            "$(command_dimensions(schema))"))
    typeof(backend(product)) === typeof(backend(binding.endpoint)) || throw(
        PlantPreparationError(:controller_output_routing, :backend,
            "controller-output product $(repr(name)) backend " *
            "$(typeof(backend(product))) does not match endpoint backend " *
            "$(typeof(backend(binding.endpoint)))"))
    expected_device = plane_device(binding.initial_command)
    actual_device = plane_device(product)
    actual_device == expected_device || throw(PlantPreparationError(
        :controller_output_routing, :physical_device,
        "controller-output product $(repr(name)) physical device " *
        "$actual_device does not match endpoint device $expected_device"))
    return product
end

function _require_controller_output_product(
    ::_PreparedPlantCommandEndpoint{
        <:PreparedCommandEndpoint{<:PlantCommandSchema{T,N}}},
    product,
    name::Symbol,
) where {T,N}
    throw(PlantPreparationError(:controller_output_routing,
        :array_storage,
        "array controller-output product $(repr(name)) must be an " *
        "AbstractArray{$T,$N}; got $(typeof(product))"))
end

@inline _require_controller_output_route(
    route::ControllerOutputRoute) = route

function _require_controller_output_route(route)
    throw(PlantPreparationError(
        :controller_output_routing, :invalid_route,
        "controller-output routes must contain ControllerOutputRoute " *
        "values; got $(typeof(route))"))
end

function _controller_output_route_tuple(routes::Tuple)
    foreach(_require_controller_output_route, routes)
    return routes
end

function _require_controller_output_route_set(products::NamedTuple,
    routes::Tuple)
    isempty(routes) && throw(PlantPreparationError(
        :controller_output_routing, :empty_routes,
        "controller-output routing requires at least one route"))
    length(routes) == length(products) || throw(PlantPreparationError(
        :controller_output_routing, :route_count,
        "controller-output routing declares $(length(products)) products " *
        "but $(length(routes)) routes"))
    product_names = Set{Symbol}()
    endpoint_ids = Set{CommandEndpointID}()
    for route in routes
        name = controller_output_product(route)
        hasproperty(products, name) || throw(PlantPreparationError(
            :controller_output_routing, :unknown_product,
            "controller-output route references undeclared product " *
            "$(repr(name))"))
        name in product_names && throw(PlantPreparationError(
            :controller_output_routing, :duplicate_product,
            "controller-output product $(repr(name)) has more than one route"))
        endpoint = controller_output_endpoint(route)
        endpoint in endpoint_ids && throw(PlantPreparationError(
            :controller_output_routing, :duplicate_endpoint,
            "command endpoint $endpoint has more than one controller-output " *
            "route"))
        push!(product_names, name)
        push!(endpoint_ids, endpoint)
    end
    for name in keys(products)
        name in product_names || throw(PlantPreparationError(
            :controller_output_routing, :unrouted_product,
            "declared controller-output product $(repr(name)) has no route"))
    end
    return nothing
end

"""
    prepare_controller_output_routing(plant, products, routes...)

Bind every named, caller-owned controller output product to one distinct
prepared command endpoint. Array products may be views into a flat controller
buffer; each view must exactly match its endpoint's element type, shape,
backend, and physical device. Scalar products use assigned `Ref` storage.

The returned routes are canonicalized by endpoint identity. They borrow the
controller products without copying and add no sequence, effective-time,
admission, transaction, optical-grouping, queue, or transport behavior.
"""
function prepare_controller_output_routing(
    plant::PreparedPlant,
    products::NamedTuple,
    routes::Tuple,
)
    checked_routes = _controller_output_route_tuple(routes)
    _require_controller_output_route_set(products, checked_routes)
    prepared = Any[]
    sizehint!(prepared, length(checked_routes))
    for route in checked_routes
        binding = _prepared_command_endpoint_binding(plant,
            controller_output_endpoint(route))
        product = getproperty(products, controller_output_product(route))
        _require_controller_output_product(binding, product,
            controller_output_product(route))
        push!(prepared, PreparedControllerOutputRoute(
            _PREPARED_CONTROLLER_OUTPUT_ROUTE_TOKEN,
            controller_output_product(route), binding.endpoint, product))
    end
    sort!(prepared; by=route ->
        String(controller_output_endpoint(route).name))
    return PreparedControllerOutputRouting(
        _PREPARED_CONTROLLER_OUTPUT_ROUTING_TOKEN, Tuple(prepared))
end

prepare_controller_output_routing(plant::PreparedPlant,
    products::NamedTuple, routes::ControllerOutputRoute...) =
    prepare_controller_output_routing(plant, products, routes)

function controller_output_route(routing::PreparedControllerOutputRouting,
    endpoint)
    resolved = _as_command_endpoint_id(endpoint)
    @inbounds for route in routing.routes
        controller_output_endpoint(route) == resolved && return route
    end
    throw(PlantPreparationError(:controller_output_routing,
        :unknown_endpoint,
        "prepared controller-output routing has no endpoint $resolved"))
end

@generated function controller_output_route(
    routing::PreparedControllerOutputRouting{R},
    ::Val{name},
) where {R<:Tuple,name}
    for (index, route_type) in enumerate(R.parameters)
        route_type.parameters[1] === name &&
            return :(getfield(routing.routes, $index))
    end
    return :(throw(PlantPreparationError(:controller_output_routing,
        :unknown_endpoint,
        "prepared controller-output routing has no endpoint " *
        string(CommandEndpointID($(QuoteNode(name)))))))
end

@inline controller_output_payload(
    routing::PreparedControllerOutputRouting, endpoint) =
    controller_output_payload(controller_output_route(routing, endpoint))
