struct ResourceFactFakeBackend <: AbstractArrayBackend end
struct ResourceFactExactFixture end
struct ResourceFactWrongOwnerFixture end

function Plant.structural_resource_fact(::ResourceFactExactFixture,
    id::StructuralResourceOwnerID, target::AbstractComputeDevice)
    return KnownStructuralResourceFact(id, target, 7, 11,
        ResourceEstimateMethod(:fixture_bytes, 1))
end

function Plant.structural_resource_fact(::ResourceFactWrongOwnerFixture,
    ::StructuralResourceOwnerID, target::AbstractComputeDevice)
    return KnownStructuralResourceFact(
        StructuralResourceOwnerID(:fixture, :wrong), target, 7, 11,
        ResourceEstimateMethod(:fixture_bytes, 1))
end

function structural_resource_fact_totals(report)
    resident = UInt64(0)
    workspace = UInt64(0)
    for fact in structural_resource_facts(report)
        resident += structural_resident_bytes(fact)
        workspace += structural_workspace_bytes(fact)
    end
    return (; resident, workspace)
end

@testset "Structural resource fact values" begin
    target = HostComputeDevice()
    owner = StructuralResourceOwnerID(:path, :science)
    method = ResourceEstimateMethod(:fixture_bytes, 3)
    fact = KnownStructuralResourceFact(owner, target, 16, 32, method)
    unknown = UnknownStructuralResourceFact(
        owner, target, :unsupported_fixture, method)
    reserve = OpaqueResourceReserve(
        owner, target, :fft_provider, 64,
        ResourceEstimateMethod(:declared_fft_reserve, 1))

    @test structural_resource_owner_id(fact) == owner
    @test structural_resource_owner_id(reserve) == owner
    @test compute_device(fact) == target
    @test compute_device(unknown) == target
    @test structural_resource_known(fact)
    @test !structural_resource_known(unknown)
    @test structural_resource_unknown_reason(unknown) ==
        :unsupported_fixture
    @test structural_resident_bytes(fact) == UInt64(16)
    @test structural_workspace_bytes(fact) == UInt64(32)
    @test opaque_resource_reserve_bytes(reserve) == UInt64(64)
    @test resource_estimate_method(fact) == method
    @test sprint(show, owner) ==
        "StructuralResourceOwnerID(:path, :science)"
    @test sprint(show, method) ==
        "ResourceEstimateMethod(:fixture_bytes, 3)"

    @test_throws StructuralResourceError StructuralResourceOwnerID(
        Symbol(""), :science)
    @test_throws StructuralResourceError StructuralResourceOwnerID(
        :path, Symbol(""))
    @test_throws StructuralResourceError ResourceEstimateMethod(
        Symbol(""), 1)
    @test_throws StructuralResourceError ResourceEstimateMethod(
        :fixture_bytes, 0)
    @test_throws StructuralResourceError ResourceEstimateMethod(
        :fixture_bytes, true)
    @test_throws StructuralResourceError ResourceEstimateMethod(
        :fixture_bytes, big(typemax(UInt32)) + 1)
    @test_throws StructuralResourceError KnownStructuralResourceFact(
        owner, target, -1, 0, method)
    @test_throws StructuralResourceError KnownStructuralResourceFact(
        owner, target, true, 0, method)
    @test_throws StructuralResourceError KnownStructuralResourceFact(
        owner, target, big(typemax(UInt64)) + 1, 0, method)
    @test_throws StructuralResourceError KnownStructuralResourceFact(
        owner, target, 1.0, 0, method)
    @test_throws StructuralResourceError UnknownStructuralResourceFact(
        owner, target, Symbol(""), method)
    @test_throws StructuralResourceError OpaqueResourceReserve(
        owner, target, Symbol(""), 0, method)
    @test_throws StructuralResourceError OpaqueResourceReserve(
        owner, target, :fft_provider, true, method)
    @test_throws StructuralResourceError OpaqueResourceReserve(
        owner, target, :fft_provider, big(typemax(UInt64)) + 1, method)

    dense = zeros(Float32, 3, 5)
    mask = Matrix{Bool}(undef, 2, 3)
    @test structural_array_bytes(dense, target) == UInt64(60)
    @test structural_array_bytes(mask, target) == UInt64(6)
    @test structural_array_bytes(
        Memory{Any}(undef, 3), target) == UInt64(3 * sizeof(Ptr{Cvoid}))
    wrong_target = AcceleratorComputeDevice(
        ResourceFactFakeBackend(), UInt32(1))
    error = try
        structural_array_bytes(dense, wrong_target)
        nothing
    catch caught
        caught
    end
    @test error isa StructuralResourceError
    @test error.reason == :wrong_device
    @test_throws StructuralResourceError structural_array_bytes(
        @view(dense[:, 1:2]), target)
    @test_throws StructuralResourceError structural_array_bytes(
        falses(2, 3), target)
    @test_throws StructuralResourceError structural_array_bytes(
        String["unsupported"], target)
    @test_throws StructuralResourceError Plant._checked_resource_multiply(
        typemax(UInt64), UInt64(2), :array_storage)
end

@testset "Structural resource aggregation" begin
    target = HostComputeDevice()
    method = ResourceEstimateMethod(:fixture_bytes, 1)
    telescope = KnownStructuralResourceFact(
        StructuralResourceOwnerID(:telescope, :primary),
        target, 100, 10, method)
    path = KnownStructuralResourceFact(
        StructuralResourceOwnerID(:path, :science),
        target, 200, 20, method)
    reserve_method = ResourceEstimateMethod(
        :fixture_fft_reserve, 2)
    reserve = OpaqueResourceReserve(
        StructuralResourceOwnerID(:path, :science),
        target, :fft_provider, 300, reserve_method)

    report = aggregate_structural_resource_facts(
        (path, telescope), target; opaque_reserves=(reserve,))
    @test compute_device(report) == target
    @test collect(structural_resource_facts(report)) == [path, telescope]
    @test collect(opaque_resource_reserves(report)) == [reserve]
    @test structural_resident_bytes(report) == UInt64(300)
    @test structural_workspace_bytes(report) == UInt64(30)
    @test opaque_resource_reserve_bytes(report) == UInt64(300)
    @test resource_estimate_method(report) == method
    @test getfield(report, :facts) isa AbstractVector{
        KnownStructuralResourceFact{HostComputeDevice}}
    @test getfield(report, :reserves) isa AbstractVector{
        OpaqueResourceReserve{HostComputeDevice}}
    @test getfield(getfield(report, :facts), :_storage) isa Tuple
    @test getfield(getfield(report, :reserves), :_storage) isa Tuple
    @test @inferred(first(structural_resource_facts(report))) === path
    @test @inferred(structural_resource_fact_totals(report)) ==
        (resident=UInt64(300), workspace=UInt64(30))
    structural_resource_fact_totals(report)
    if !coverage_instrumented()
        @test @allocated(structural_resource_fact_totals(report)) == 0
    end
    reversed_report = aggregate_structural_resource_facts(
        (telescope, path), target; opaque_reserves=(reserve,))
    @test structural_resource_facts(reversed_report) ==
        structural_resource_facts(report)
    @test opaque_resource_reserves(reversed_report) ==
        opaque_resource_reserves(report)

    unknown = UnknownStructuralResourceFact(
        StructuralResourceOwnerID(:acquisition, :science_camera),
        target, :unsupported_model, method)
    error = try
        aggregate_structural_resource_facts((telescope, unknown), target)
        nothing
    catch caught
        caught
    end
    @test error isa StructuralResourceError
    @test error.reason == :unknown_resource

    duplicate = KnownStructuralResourceFact(
        structural_resource_owner_id(telescope), target, 1, 1, method)
    error = try
        aggregate_structural_resource_facts(
            (telescope, duplicate), target)
        nothing
    catch caught
        caught
    end
    @test error isa StructuralResourceError
    @test error.reason == :duplicate_owner

    inconsistent = KnownStructuralResourceFact(
        StructuralResourceOwnerID(:path, :other), target, 1, 1,
        ResourceEstimateMethod(:fixture_bytes, 2))
    error = try
        aggregate_structural_resource_facts(
            (telescope, inconsistent), target)
        nothing
    catch caught
        caught
    end
    @test error isa StructuralResourceError
    @test error.reason == :inconsistent_estimate_method

    wrong_target = AcceleratorComputeDevice(
        ResourceFactFakeBackend(), UInt32(1))
    wrong = KnownStructuralResourceFact(
        StructuralResourceOwnerID(:path, :wrong), wrong_target,
        1, 1, method)
    error = try
        aggregate_structural_resource_facts((telescope, wrong), target)
        nothing
    catch caught
        caught
    end
    @test error isa StructuralResourceError
    @test error.reason == :wrong_device

    overflow = KnownStructuralResourceFact(
        StructuralResourceOwnerID(:path, :overflow), target,
        typemax(UInt64), 0, method)
    one_byte = KnownStructuralResourceFact(
        StructuralResourceOwnerID(:path, :one_byte), target,
        1, 0, method)
    error = try
        aggregate_structural_resource_facts((overflow, one_byte), target)
        nothing
    catch caught
        caught
    end
    @test error isa StructuralResourceError
    @test error.reason == :arithmetic_overflow

    workspace_overflow = KnownStructuralResourceFact(
        StructuralResourceOwnerID(:path, :workspace_overflow), target,
        0, typemax(UInt64), method)
    workspace_one_byte = KnownStructuralResourceFact(
        StructuralResourceOwnerID(:path, :workspace_one_byte), target,
        0, 1, method)
    error = try
        aggregate_structural_resource_facts(
            (workspace_overflow, workspace_one_byte), target)
        nothing
    catch caught
        caught
    end
    @test error isa StructuralResourceError
    @test error.reason == :arithmetic_overflow

    reserve_overflow = OpaqueResourceReserve(
        structural_resource_owner_id(telescope), target,
        :fft_provider_a, typemax(UInt64), reserve_method)
    reserve_one_byte = OpaqueResourceReserve(
        structural_resource_owner_id(telescope), target,
        :fft_provider_b, 1, reserve_method)
    error = try
        aggregate_structural_resource_facts(
            (telescope,), target;
            opaque_reserves=(reserve_overflow, reserve_one_byte))
        nothing
    catch caught
        caught
    end
    @test error isa StructuralResourceError
    @test error.reason == :arithmetic_overflow

    duplicate_reserve = OpaqueResourceReserve(
        structural_resource_owner_id(reserve), target,
        reserve.provider, 1, reserve_method)
    error = try
        aggregate_structural_resource_facts(
            (path,), target;
            opaque_reserves=(reserve, duplicate_reserve))
        nothing
    catch caught
        caught
    end
    @test error isa StructuralResourceError
    @test error.reason == :duplicate_reserve

    wrong_reserve = OpaqueResourceReserve(
        StructuralResourceOwnerID(:path, :wrong_reserve), wrong_target,
        :fft_provider, 1, reserve_method)
    error = try
        aggregate_structural_resource_facts(
            (telescope,), target; opaque_reserves=(wrong_reserve,))
        nothing
    catch caught
        caught
    end
    @test error isa StructuralResourceError
    @test error.reason == :wrong_device

    unknown_owner_reserve = OpaqueResourceReserve(
        StructuralResourceOwnerID(:path, :unregistered), target,
        :fft_provider, 1, reserve_method)
    error = try
        aggregate_structural_resource_facts(
            (telescope,), target;
            opaque_reserves=(unknown_owner_reserve,))
        nothing
    catch caught
        caught
    end
    @test error isa StructuralResourceError
    @test error.reason == :unknown_owner

    @test_throws StructuralResourceError aggregate_structural_resource_facts(
        (), target)
    @test_throws StructuralResourceError aggregate_structural_resource_facts(
        (telescope, 1), target)
    @test_throws StructuralResourceError aggregate_structural_resource_facts(
        (telescope,), target; opaque_reserves=(1,))

    unsupported_owner = (fixture=:unsupported,)
    fact = structural_resource_fact(
        unsupported_owner,
        StructuralResourceOwnerID(:fixture, :unsupported),
        target,
    )
    @test !structural_resource_known(fact)
    @test_throws StructuralResourceError require_exact_structural_resource_facts(
        unsupported_owner,
        StructuralResourceOwnerID(:fixture, :unsupported),
        target,
    )

    fixture_id = StructuralResourceOwnerID(:fixture, :exact)
    fixture_report = require_exact_structural_resource_facts(
        ResourceFactExactFixture(), fixture_id, target)
    @test only(structural_resource_facts(fixture_report)) ==
        KnownStructuralResourceFact(fixture_id, target, 7, 11,
            ResourceEstimateMethod(:fixture_bytes, 1))

    error = try
        require_exact_structural_resource_facts(
            ResourceFactWrongOwnerFixture(), fixture_id, target)
        nothing
    catch caught
        caught
    end
    @test error isa StructuralResourceError
    @test error.reason == :wrong_owner

    wrong_owner_reserve = OpaqueResourceReserve(
        StructuralResourceOwnerID(:fixture, :wrong), target,
        :fft_provider, 1, reserve_method)
    error = try
        require_exact_structural_resource_facts(
            ResourceFactExactFixture(), fixture_id, target;
            opaque_reserves=(wrong_owner_reserve,))
        nothing
    catch caught
        caught
    end
    @test error isa StructuralResourceError
    @test error.reason == :wrong_owner
end
