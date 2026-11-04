@testset "external pupils preserve selected optical precision and Float64 exposure means" begin
    domain = ntuple(_ -> UInt8(1), 16)
    for T in (Float32, Float64)
        scc = external_test_graph(; T)
        llowfs = external_test_graph(; T)
        unfringed = external_test_graph(; T, steps=((op="multiply", amplitude=0.5),))
        owner = prepare_external_pupil((; scc, llowfs);
            exposure_durations=(scc=ModelDuration(4), llowfs=ModelDuration(2)),
            observations=(scc=(output=:intensity, sampling_node=:optics,
                chopping=(; unfringed, initial_phase=:fringed)),
                llowfs=(output=:intensity, sampling_node=:optics)))
        @test eltype(owner.pupil_opd) === eltype(owner.pupil_amplitude) === T
        @test eltype(owner.exposures.scc.integrated) === Float64
        products = NamedTuple[]
        emit(name, image, metadata) = push!(products, (; name, image=copy(image), metadata))
        for (sequence, (time, span)) in enumerate(((0,1), (1,2), (3,2), (5,2), (7,1)))
            @test step_external_pupil!(emit, owner, zeros(Float32,4,4), ones(Float32,4,4),
                ModelTimestamp(time), ModelDuration(span); domain,
                generation=UInt64(1), sequence=UInt64(sequence))
        end
        chopped = filter(p -> p.name == :scc, products)
        @test [p.metadata.chopped_image_kind for p in chopped] == [:fringed, :unfringed]
        @test [p.metadata.pts for p in chopped] == [0, 4]
        @test all(p -> p.image == fill(p.metadata.pts == 0 ? 1.0 : 0.25, 4,4), chopped)
        @test length(filter(p -> p.name == :llowfs, products)) == 4
        @test all(p -> eltype(p.image) === Float64, products)
        reset_external_pupil!(owner)
        @test owner.state.domain === nothing
        @test all(iszero, owner.exposures.scc.image)
    end
    @test_throws AlgorithmGraphError external_test_owner((
        slow=external_test_graph(; T=Float32), fast=external_test_graph()))
    @test_throws AlgorithmGraphError chopped_test_owner(
        fringed=external_test_graph(; T=Float32), unfringed=external_test_graph())

    owner = external_test_owner((slow=external_test_graph(; T=Float32),
        fast=external_test_graph(; T=Float32)))
    emit = (args...) -> nothing
    @test step_external_pupil!(emit, owner, zeros(4,4), ones(4,4),
        ModelTimestamp(0), ModelDuration(2); domain, generation=UInt64(1), sequence=UInt64(1))
    before = map(exposure_state, owner.exposures)
    state = (owner.state.domain, owner.state.generation, owner.state.sequence, owner.state.discontinuities)
    inputs = (copy(owner.pupil_opd), copy(owner.pupil_amplitude))
    for (opd, amplitude) in ((fill(Float64(floatmax(Float32))*2,4,4), ones(4,4)),
            (zeros(4,4), fill(Float64(floatmax(Float32))*2,4,4)))
        @test_throws AlgorithmGraphError step_external_pupil!(emit, owner, opd, amplitude,
            ModelTimestamp(2), ModelDuration(2); domain, generation=UInt64(2), sequence=UInt64(1))
        @test map(exposure_state, owner.exposures) == before
        @test (owner.state.domain, owner.state.generation, owner.state.sequence, owner.state.discontinuities) == state
        @test (owner.pupil_opd, owner.pupil_amplitude) == inputs
    end
    # Blocked OPD remains sanitized before conversion, even if it is nonfinite.
    @test step_external_pupil!(emit, owner, fill(NaN,4,4), zeros(4,4),
        ModelTimestamp(2), ModelDuration(2); domain, generation=UInt64(1), sequence=UInt64(2))
    @test all(iszero, owner.pupil_opd)

    # Multiply the optical value by nanoseconds in Float64, as well as accumulating
    # in Float64; a finite Float32 optical value must not overflow in that product.
    exposure = IntensityExposure((1,1), ModelDuration(10_000_000))
    value = 1f35
    images = Matrix{Float64}[]
    integrate_intensity!((image,args...) -> push!(images,copy(image)), exposure,
        fill(value,1,1), ModelTimestamp(0), ModelDuration(10_000_000))
    @test only(images)[1] ≈ Float64(value)
end

@testset "prepared exposure kernels preserve validation and ordered Float64 arithmetic" begin
    # Run the portable kernels on the host independently of graph capture. The
    # oracle is the original exposure broadcast, including conversion order.
    cpu = AlgorithmGraphs.KernelAbstractions.CPU()
    for T in (Float32,Float64), count in (1,257,1025)
        values = reshape(T[isodd(i) ? 1e-3 : 1e30 for i in 1:count],count,1)
        flags = fill(UInt8(255),cld(count,256))
        validate!() = AlgorithmGraphs._validate_external_intensity_kernel!(cpu)(
            flags,values,Float32,1.0;ndrange=256length(flags),workgroupsize=256)
        validate!()
        @test all(iszero,flags)
        for invalid in (T(-1),T(NaN),T(Inf))
            values[end] = invalid
            validate!()
            @test any(flag -> !iszero(flag & UInt8(1)),flags)
            values[end] = T(1)
            validate!()
            @test all(iszero,flags) # A replay overwrites earlier rejection flags.
        end
        values[end] = T(1e35)
        AlgorithmGraphs._validate_external_intensity_kernel!(cpu)(
            flags,values,Float32,1e10;ndrange=256length(flags),workgroupsize=256)
        @test all(flag -> iszero(flag & UInt8(1)),flags)
        @test any(flag -> !iszero(flag & UInt8(2)),flags)
        actual = zeros(Float64,size(values))
        expected = copy(actual)
        overlap = Int64[0]
        for span in (1,3,10_000_000,7)
            overlap[1] = span
            expected .+= Float64.(values) .* span
            AlgorithmGraphs._accumulate_external_intensity_kernel!(cpu)(
                actual,values,overlap;ndrange=length(values))
            @test isequal(actual,expected)
        end
    end
end

@testset "exposure callback failure retains publication and clear ordering" begin
    exposure = IntensityExposure((1,1),ModelDuration(4))
    failure = ErrorException("consumer rejected image")
    callback = function (image,pts,duration,sequence)
        @test image == fill(3.0,1,1)
        @test exposure.integrated == fill(12.0,1,1)
        @test exposure.covered_ns == 4 && exposure.next_ns == 4
        @test exposure.start_ns == 0 && exposure.sequence == sequence == 1
        throw(failure)
    end
    @test_throws typeof(failure) integrate_intensity!(callback,exposure,fill(3f0,1,1),
        ModelTimestamp(0),ModelDuration(4))
    @test exposure.integrated == fill(12.0,1,1)
    @test exposure.image == fill(3.0,1,1)
    @test exposure.covered_ns == 4 && exposure.next_ns == 4
    @test exposure.start_ns == 0 && exposure.sequence == 1
end
