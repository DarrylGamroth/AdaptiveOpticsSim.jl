struct ControllerRoutingTestAtmosphere <: AdaptiveOpticsSim.AbstractAtmosphere end
struct ControllerRoutingTestOpticModel end
struct PreparedControllerRoutingTestOptic end

Plant.plant_model_definition_style(
    ::Type{ControllerRoutingTestOpticModel}) = ColdPlantModelDefinition()

Plant.prepare_controllable_optic(
    ::ControllerRoutingTestOpticModel,
    ::ControllableOpticDefinition,
    ::Telescope,
    ::AdaptiveOpticsSim.AbstractAtmosphere,
) = PreparedControllerRoutingTestOptic()

function controller_routing_schema(::Type{T}, endpoint::Symbol,
    dimensions::Tuple) where {T<:AbstractFloat}
    return PlantCommandSchema(
        T,
        dimensions;
        id=Symbol(endpoint, :_schema),
        version=1,
        endpoint,
        units=:metre,
        sign_convention=:positive_surface_increases_opd,
        basis=CommandBasis(:test, Symbol(endpoint, :_basis)),
        basis_revision=1,
        semantics=AbsoluteCommand,
        bounds=UnboundedCommandValues(),
        value_policy=CommandValuePolicy(),
        sequence_policy=CommandSequencePolicy(),
        effective_time_policy=CommandEffectiveTimePolicy(),
        silence_policy=CommandSilencePolicy(),
    )
end

function controller_routing_plant(;
    woofer_dimensions=(2,),
    tweeter_dimensions=(3,),
    scalar=false,
)
    tel = Telescope(resolution=8, diameter=8.0,
        central_obstruction=0.0)
    atmosphere = ControllerRoutingTestAtmosphere()
    woofer_schema = controller_routing_schema(Float64, :woofer_command,
        scalar ? () : woofer_dimensions)
    tweeter_schema = controller_routing_schema(Float64, :tweeter_command,
        tweeter_dimensions)
    definition = PlantDefinition(
        telescope=tel,
        atmosphere=atmosphere,
        controllable_optics=(
            ControllableOpticDefinition(:woofer,
                ControllerRoutingTestOpticModel(),
                (woofer_schema,);
                placement=PupilPlanePlacement(),
                visibility=AllPathVisibility()),
            ControllableOpticDefinition(:tweeter,
                ControllerRoutingTestOpticModel(),
                (tweeter_schema,);
                placement=PupilPlanePlacement(),
                visibility=AllPathVisibility()),
        ),
    )
    initial_woofer = scalar ? 0.0 : zeros(Float64, woofer_dimensions)
    return prepare_plant(definition;
        run_seed=1,
        command_endpoints=(
            CommandEndpointConfiguration(:woofer_command,
                initial_woofer; capacity=2),
            CommandEndpointConfiguration(:tweeter_command,
                zeros(Float64, tweeter_dimensions); capacity=2),
        ))
end

function captured_controller_routing_error(f)
    try
        f()
    catch error
        return error
    end
    return nothing
end

function assert_controller_routing_error(f, reason::Symbol)
    error = captured_controller_routing_error(f)
    @test error isa PlantPreparationError
    if error isa PlantPreparationError
        @test error.component === :controller_output_routing
        @test error.reason === reason
    end
    return error
end

@testset "Prepared controller-output routing" begin
    plant = controller_routing_plant()
    flat_output = zeros(Float64, 5)
    products = (
        high_order=@view(flat_output[3:5]),
        low_order=@view(flat_output[1:2]),
    )
    routing = prepare_controller_output_routing(
        plant,
        products,
        ControllerOutputRoute(:high_order, :tweeter_command),
        ControllerOutputRoute(:low_order, :woofer_command),
    )

    routes = prepared_controller_output_routes(routing)
    @test map(controller_output_endpoint, routes) ==
        (CommandEndpointID(:tweeter_command),
            CommandEndpointID(:woofer_command))
    @test map(controller_output_product, routes) ==
        (:high_order, :low_order)
    @test controller_output_schema(first(routes)) ===
        command_schema(prepared_command_endpoint(
            plant, :tweeter_command))
    @test controller_output_payload(routing, :woofer_command) ===
        products.low_order
    @test controller_output_payload(routing,
        CommandEndpointID(:tweeter_command)) === products.high_order

    flat_output .= (1.0, 2.0, 3.0, 4.0, 5.0)
    @test controller_output_payload(routing, :woofer_command) ==
        [1.0, 2.0]
    @test controller_output_payload(routing, :tweeter_command) ==
        [3.0, 4.0, 5.0]

    reordered = prepare_controller_output_routing(
        plant,
        products,
        (
            ControllerOutputRoute(:low_order, :woofer_command),
            ControllerOutputRoute(:high_order, :tweeter_command),
        ),
    )
    @test map(controller_output_endpoint,
        prepared_controller_output_routes(reordered)) ==
        map(controller_output_endpoint, routes)

    if coverage_instrumented()
        @test_skip "routing allocation assertions are disabled under coverage instrumentation"
    else
        @test @inferred(controller_output_payload(
            routing, Val(:woofer_command))) === products.low_order
        @test @allocated(controller_output_payload(
            routing, Val(:woofer_command))) == 0
        @test @allocated(controller_output_payload(
            first(routes))) == 0
    end

    @test Base.isexported(Plant, :ControllerOutputRoute)
    @test Base.isexported(Plant, :prepare_controller_output_routing)
    @test !Base.isexported(AdaptiveOpticsSim, :ControllerOutputRoute)
    @test Base.ispublic(Plant, :PreparedControllerOutputRouting)
    @test Base.ispublic(Plant, :controller_output_payload)
end

@testset "Controller-output routing validation" begin
    plant = controller_routing_plant()
    outputs = zeros(Float64, 5)
    products = (
        low=@view(outputs[1:2]),
        high=@view(outputs[3:5]),
    )
    valid = (
        ControllerOutputRoute(:low, :woofer_command),
        ControllerOutputRoute(:high, :tweeter_command),
    )

    assert_controller_routing_error(:empty_product) do
        ControllerOutputRoute(Symbol(""), :woofer_command)
    end
    assert_controller_routing_error(:invalid_product) do
        ControllerOutputRoute("low", :woofer_command)
    end
    assert_controller_routing_error(:empty_routes) do
        prepare_controller_output_routing(plant, NamedTuple(), ())
    end
    assert_controller_routing_error(:route_count) do
        prepare_controller_output_routing(plant, products, (first(valid),))
    end
    assert_controller_routing_error(:unknown_product) do
        prepare_controller_output_routing(plant, products, (
            ControllerOutputRoute(:missing, :woofer_command),
            last(valid),
        ))
    end
    assert_controller_routing_error(:duplicate_product) do
        prepare_controller_output_routing(plant, products, (
            ControllerOutputRoute(:low, :woofer_command),
            ControllerOutputRoute(:low, :tweeter_command),
        ))
    end
    assert_controller_routing_error(:duplicate_endpoint) do
        prepare_controller_output_routing(plant, products, (
            ControllerOutputRoute(:low, :woofer_command),
            ControllerOutputRoute(:high, :woofer_command),
        ))
    end
    assert_controller_routing_error(:unknown_endpoint) do
        prepare_controller_output_routing(plant, products, (
            ControllerOutputRoute(:low, :missing_command),
            last(valid),
        ))
    end
    assert_controller_routing_error(:shape) do
        prepare_controller_output_routing(plant,
            (low=zeros(Float64, 3), high=products.high), valid)
    end
    assert_controller_routing_error(:element_type) do
        prepare_controller_output_routing(plant,
            (low=zeros(Float32, 2), high=products.high), valid)
    end
    assert_controller_routing_error(:array_storage) do
        prepare_controller_output_routing(plant,
            (low=1.0, high=products.high), valid)
    end
    assert_controller_routing_error(:unknown_endpoint) do
        routing = prepare_controller_output_routing(plant, products, valid)
        controller_output_route(routing, :missing_command)
    end
    assert_controller_routing_error(:unknown_endpoint) do
        routing = prepare_controller_output_routing(plant, products, valid)
        controller_output_route(routing, Val(:missing_command))
    end
    assert_controller_routing_error(:invalid_route) do
        prepare_controller_output_routing(plant, products,
            (first(valid), :not_a_route))
    end
end

@testset "Scalar controller-output routing" begin
    plant = controller_routing_plant(scalar=true)
    scalar_output = Ref(1.25)
    vector_output = zeros(Float64, 3)
    routing = prepare_controller_output_routing(
        plant,
        (low=scalar_output, high=vector_output),
        ControllerOutputRoute(:low, :woofer_command),
        ControllerOutputRoute(:high, :tweeter_command),
    )
    @test controller_output_payload(routing, :woofer_command) == 1.25
    scalar_output[] = -0.5
    @test controller_output_payload(routing, :woofer_command) == -0.5

    assert_controller_routing_error(:scalar_storage) do
        prepare_controller_output_routing(
            plant,
            (low=1.0, high=vector_output),
            ControllerOutputRoute(:low, :woofer_command),
            ControllerOutputRoute(:high, :tweeter_command),
        )
    end
    assert_controller_routing_error(:scalar_storage) do
        prepare_controller_output_routing(
            plant,
            (low=fill(1.0), high=vector_output),
            ControllerOutputRoute(:low, :woofer_command),
            ControllerOutputRoute(:high, :tweeter_command),
        )
    end
    assert_controller_routing_error(:element_type) do
        prepare_controller_output_routing(
            plant,
            (low=Ref(1.0f0), high=vector_output),
            ControllerOutputRoute(:low, :woofer_command),
            ControllerOutputRoute(:high, :tweeter_command),
        )
    end
    assert_controller_routing_error(:unassigned_scalar_product) do
        prepare_controller_output_routing(
            plant,
            (low=Ref{AbstractFloat}(), high=vector_output),
            ControllerOutputRoute(:low, :woofer_command),
            ControllerOutputRoute(:high, :tweeter_command),
        )
    end
end
