function direction_batch_test_atmosphere(
    ::Type{MultiLayerAtmosphere},
    tel::Telescope;
    T::Type{<:AbstractFloat}=Float64,
)
    return MultiLayerAtmosphere(
        tel;
        r0=T(0.2),
        L0=T(25),
        fractional_cn2=T[0.65, 0.35],
        wind_speed=T[7, 3],
        wind_direction_deg=T[20, 115],
        altitude=T[0, 4_000],
        T=T,
        backend=backend(tel),
    )
end

function direction_batch_test_atmosphere(
    ::Type{InfiniteMultiLayerAtmosphere},
    tel::Telescope;
    T::Type{<:AbstractFloat}=Float64,
)
    return InfiniteMultiLayerAtmosphere(
        tel;
        r0=T(0.2),
        L0=T(25),
        fractional_cn2=T[0.65, 0.35],
        wind_speed=T[7, 3],
        wind_direction_deg=T[20, 115],
        altitude=T[0, 4_000],
        screen_resolution=25,
        stencil_size=27,
        T=T,
        backend=backend(tel),
    )
end

function direction_batch_test_sources(::Type{T}=Float64) where {
    T<:AbstractFloat,
}
    return Asterism(AdaptiveOpticsSim.Optics.AbstractSource[
        Source(
            band=:I,
            magnitude=zero(T),
            coordinates=(zero(T), zero(T)),
            T=T,
        ),
        Source(
            band=:I,
            magnitude=zero(T),
            coordinates=(T(8), T(35)),
            T=T,
        ),
        LGSSource(
            magnitude=zero(T),
            coordinates=(T(-12), T(70)),
            altitude=T(90_000),
            T=T,
        ),
        Source(
            band=:K,
            magnitude=one(T),
            coordinates=(T(5), T(120)),
            T=T,
        ),
    ])
end

function direction_batch_has_concrete_members(::Type{T}) where {T}
    members = T isa Union ? Base.uniontypes(T) : (T,)
    return all(isconcretetype, members)
end

mutable struct MutableDirectionSource{T<:AbstractFloat} <:
    AdaptiveOpticsSim.Optics.AbstractSource
    coordinates_xy_arcsec::NTuple{2,T}
    height_m::T
end

AdaptiveOpticsSim.Optics.coordinates_xy_arcsec(src::MutableDirectionSource) =
    src.coordinates_xy_arcsec
AdaptiveOpticsSim.Optics.source_height_m(src::MutableDirectionSource) = src.height_m

function serial_atmosphere_direction_stack(
    atm::AbstractTimedAtmosphere,
    tel::Telescope,
    sources,
    epoch::AtmosphereEpoch,
)
    T = atmosphere_numeric_type(atm)
    output = Array{T}(undef, tel.params.resolution, tel.params.resolution,
        length(sources))
    @inbounds for direction in eachindex(sources)
        pupil = PupilFunction(tel; T=T)
        renderer = prepare_atmosphere_renderer(
            atm,
            tel,
            sources[direction];
            T=T,
        )
        render_atmosphere!(pupil, renderer, atm, epoch)
        copyto!(@view(output[:, :, direction]), pupil.opd)
    end
    return output
end

@testset "Homogeneous atmosphere-direction batch parity" begin
    for AtmosphereType in (
        MultiLayerAtmosphere,
        InfiniteMultiLayerAtmosphere,
    )
        T = Float64
        n = 12
        tel = Telescope(
            resolution=n,
            diameter=T(8),
            central_obstruction=T(0.15),
            T=T,
        )
        sources = direction_batch_test_sources(T)
        atm = direction_batch_test_atmosphere(AtmosphereType, tel; T=T)
        output = fill(T(19), n, n, length(sources))

        # Preparation is independent of atmosphere publication.
        prepared = prepare_atmosphere_direction_batch(
            atm,
            tel,
            sources,
            output,
        )
        @test prepared isa PreparedAtmosphereDirectionBatch
        @test prepared.params.sources isa FixedSizeVector
        @test prepared.params.source_geometry_bindings isa FixedSizeVector
        @test prepared.params.source_bindings isa FixedSizeVector
        @test prepared.params.layer_bindings isa FixedSizeVector
        @test all(
            direction_batch_has_concrete_members(eltype(storage))
            for storage in (
                prepared.params.sources,
                prepared.params.source_geometry_bindings,
                prepared.params.source_bindings,
                prepared.params.layer_bindings,
            )
        )
        @test all(isconcretetype, fieldtypes(typeof(prepared.params)))
        @test atmosphere_direction_output(prepared) === output
        @test atmosphere_direction_count(prepared) == length(sources)
        @test atmosphere_direction_capacity(prepared) == length(sources)
        @test plane_metadata(prepared) ===
            atmosphere_direction_metadata(prepared)
        @test atmosphere_direction_metadata(prepared).dimensions == (n, n)
        @test atmosphere_direction_metadata(prepared).numeric_type === T
        @test atmosphere_direction_metadata(prepared).normalization isa
            UnspecifiedNormalization
        @test atmosphere_direction_metadata(prepared).coherence isa
            NonCombinableProduct
        @test backend(prepared.workspace.shift_x) isa CPUBackend
        @test compute_device(prepared.workspace.shift_x) ==
            HostComputeDevice()
        @test size(prepared.workspace.shift_x) ==
            (length(atm.layers), length(sources))
        @test prepared.workspace.footprint_scale[2, 3] < one(T)
        @test prepared.workspace.footprint_scale[2, 2] == one(T)

        epoch = advance_by!(
            atm,
            T(1e-3);
            rng=MersenneTwister(70),
        )
        rendered = @inferred render_atmosphere_directions!(
            prepared,
            atm,
            epoch,
        )
        reference = serial_atmosphere_direction_stack(
            atm,
            tel,
            prepared.params.sources,
            epoch,
        )
        @test rendered === output
        @test output == reference
        @test @view(output[:, :, 1]) != @view(output[:, :, 2])
        @test @view(output[:, :, 2]) != @view(output[:, :, 3])
        @test validate_atmosphere_direction_batch(
            prepared,
            atm,
            epoch,
        ) === prepared

        render_atmosphere_directions!(prepared, atm, epoch)
        if coverage_instrumented()
            @test_skip "allocation assertions are disabled under coverage instrumentation"
        else
            @test @allocated(render_atmosphere_directions!(
                prepared,
                atm,
                epoch,
            )) == 0
        end

        single_source = first(prepared.params.sources)
        single_output = fill(T(23), n, n, 1)
        single = prepare_atmosphere_direction_batch(
            atm,
            tel,
            single_source,
            single_output,
        )
        render_atmosphere_directions!(single, atm, epoch)
        single_pupil = PupilFunction(tel; T=T)
        single_renderer = prepare_atmosphere_renderer(
            atm,
            tel,
            single_source;
            T=T,
        )
        render_atmosphere!(
            single_pupil,
            single_renderer,
            atm,
            epoch,
        )
        @test @view(single_output[:, :, 1]) == single_pupil.opd
    end
end

@testset "Atmosphere-direction source expansion" begin
    T = Float64
    n = 10
    tel = Telescope(
        resolution=n,
        diameter=T(8),
        central_obstruction=zero(T),
        T=T,
    )
    atm = direction_batch_test_atmosphere(
        MultiLayerAtmosphere,
        tel;
        T=T,
    )
    source = Source(
        band=:I,
        magnitude=zero(T),
        coordinates=(T(3), T(40)),
        T=T,
    )
    extended = with_extended_source(
        source,
        GaussianDiskSourceModel(
            sigma_arcsec=T(0.2),
            n_side=3,
            T=T,
        ),
    )
    output = zeros(T, n, n, 9)
    prepared = prepare_atmosphere_direction_batch(
        atm,
        tel,
        extended,
        output,
    )
    epoch = advance_by!(atm, T(1e-3); rng=MersenneTwister(71))
    render_atmosphere_directions!(prepared, atm, epoch)
    @test atmosphere_direction_count(prepared) == 9
    @test output == serial_atmosphere_direction_stack(
        atm,
        tel,
        prepared.params.sources,
        epoch,
    )

    original_coordinates = map(
        coordinates_xy_arcsec,
        prepared.params.sources,
    )
    extended.model.offsets_xy_arcsec[1] = (T(99), T(99))
    @test map(coordinates_xy_arcsec, prepared.params.sources) ==
        original_coordinates
end

@testset "Atmosphere-direction batch preflight and fallback" begin
    T = Float64
    n = 10
    tel = Telescope(
        resolution=n,
        diameter=T(8),
        central_obstruction=zero(T),
        T=T,
    )
    sources = direction_batch_test_sources(T)
    atm = direction_batch_test_atmosphere(
        MultiLayerAtmosphere,
        tel;
        T=T,
    )
    output = fill(T(31), n, n, length(sources))
    prepared = prepare_atmosphere_direction_batch(
        atm,
        tel,
        sources,
        output,
    )
    epoch = advance_by!(atm, T(1e-3); rng=MersenneTwister(72))
    render_atmosphere_directions!(prepared, atm, epoch)

    fill!(output, T(31))
    stale_epoch = epoch
    current = advance_by!(atm, T(1e-3); rng=MersenneTwister(72))
    @test_throws AtmosphereEpochError render_atmosphere_directions!(
        prepared,
        atm,
        stale_epoch,
    )
    @test all(==(T(31)), output)

    other = direction_batch_test_atmosphere(
        MultiLayerAtmosphere,
        tel;
        T=T,
    )
    other_epoch = advance_by!(other, T(2e-3); rng=MersenneTwister(73))
    @test_throws AtmosphereEpochError render_atmosphere_directions!(
        prepared,
        other,
        other_epoch,
    )
    @test all(==(T(31)), output)

    prepared.params.sources[1], prepared.params.sources[2] =
        prepared.params.sources[2], prepared.params.sources[1]
    @test_throws AtmosphereEpochError render_atmosphere_directions!(
        prepared,
        atm,
        current,
    )
    @test all(==(T(31)), output)
    prepared.params.sources[1], prepared.params.sources[2] =
        prepared.params.sources[2], prepared.params.sources[1]

    mutable_source = MutableDirectionSource(
        (T(2), T(4)),
        T(Inf),
    )
    mutable_output = fill(T(31), n, n, 1)
    mutable_prepared = prepare_atmosphere_direction_batch(
        atm,
        tel,
        mutable_source,
        mutable_output,
    )
    mutable_source.coordinates_xy_arcsec = (T(3), T(4))
    @test_throws AtmosphereEpochError render_atmosphere_directions!(
        mutable_prepared,
        atm,
        current,
    )
    @test all(==(T(31)), mutable_output)

    atm.layers[1], atm.layers[2] = atm.layers[2], atm.layers[1]
    @test_throws AtmosphereEpochError render_atmosphere_directions!(
        prepared,
        atm,
        current,
    )
    @test all(==(T(31)), output)
    atm.layers[1], atm.layers[2] = atm.layers[2], atm.layers[1]

    saved_output = prepared.workspace.output
    wrong_capacity = fill(T(31), n, n, length(sources) + 1)
    prepared.workspace.output = wrong_capacity
    @test_throws DimensionMismatchError render_atmosphere_directions!(
        prepared,
        atm,
        current,
    )
    @test all(==(T(31)), wrong_capacity)
    prepared.workspace.output = saved_output

    bad_geometry = fill(
        T(31),
        size(prepared.workspace.shift_x, 1) + 1,
        size(prepared.workspace.shift_x, 2),
    )
    bad_workspace = AtmosphereDirectionBatchWorkspace(
        bad_geometry,
        prepared.workspace.shift_y,
        prepared.workspace.footprint_scale,
        prepared.workspace.pupil,
        prepared.workspace.output,
    )
    bad_geometry_prepared = PreparedAtmosphereDirectionBatch(
        prepared.params,
        bad_workspace,
    )
    @test_throws DimensionMismatchError render_atmosphere_directions!(
        bad_geometry_prepared,
        atm,
        current,
    )
    @test all(==(T(31)), output)

    wrong_type_output = fill(Float32(31), n, n, length(sources))
    wrong_type_workspace = AtmosphereDirectionBatchWorkspace(
        prepared.workspace.shift_x,
        prepared.workspace.shift_y,
        prepared.workspace.footprint_scale,
        prepared.workspace.pupil,
        wrong_type_output,
    )
    wrong_type_prepared = PreparedAtmosphereDirectionBatch(
        prepared.params,
        wrong_type_workspace,
    )
    @test_throws InvalidConfiguration render_atmosphere_directions!(
        wrong_type_prepared,
        atm,
        current,
    )
    @test all(==(Float32(31)), wrong_type_output)

    foreign_device_output = ContractDeviceArray(
        fill(T(31), n, n, length(sources)),
        ContractComputeDevice(2),
    )
    foreign_device_workspace = AtmosphereDirectionBatchWorkspace(
        prepared.workspace.shift_x,
        prepared.workspace.shift_y,
        prepared.workspace.footprint_scale,
        prepared.workspace.pupil,
        foreign_device_output,
    )
    foreign_device_prepared = PreparedAtmosphereDirectionBatch(
        prepared.params,
        foreign_device_workspace,
    )
    @test_throws InvalidConfiguration render_atmosphere_directions!(
        foreign_device_prepared,
        atm,
        current,
    )
    @test all(==(T(31)), foreign_device_output.storage)

    wrong_grid_telescope = Telescope(
        resolution=n,
        diameter=T(10),
        central_obstruction=zero(T),
        T=T,
    )
    wrong_grid_output = fill(T(31), n, n, length(sources))
    @test_throws InvalidConfiguration prepare_atmosphere_direction_batch(
        atm,
        wrong_grid_telescope,
        sources,
        wrong_grid_output,
    )
    @test all(==(T(31)), wrong_grid_output)

    @test_throws DimensionMismatchError prepare_atmosphere_direction_batch(
        atm,
        tel,
        sources,
        fill(T(31), n, n, length(sources) - 1),
    )
    @test_throws DimensionMismatchError prepare_atmosphere_direction_batch(
        atm,
        tel,
        sources,
        fill(T(31), n, n),
    )
    @test_throws InvalidConfiguration prepare_atmosphere_direction_batch(
        atm,
        tel,
        sources,
        fill(Float32(31), n, n, length(sources)),
    )
    @test_throws InvalidConfiguration prepare_atmosphere_direction_batch(
        atm,
        tel,
        Asterism(AdaptiveOpticsSim.Optics.AbstractSource[]),
        Array{T}(undef, n, n, 0),
    )

    spectral = with_spectrum(
        first(sources.sources),
        SpectralBundle(T[700e-9, 800e-9], T[0.5, 0.5]; T=T),
    )
    @test_throws UnsupportedAlgorithm prepare_atmosphere_direction_batch(
        atm,
        tel,
        spectral,
        Array{T}(undef, n, n, 1),
    )

    kolmogorov = KolmogorovAtmosphere(tel; r0=T(0.2), T=T)
    @test atmosphere_direction_batch_capability(typeof(kolmogorov)) isa
        UnsupportedAtmosphereDirectionBatchCapability
    @test_throws UnsupportedAlgorithm prepare_atmosphere_direction_batch(
        kolmogorov,
        tel,
        first(sources.sources),
        Array{T}(undef, n, n, 1),
    )
    @test atmosphere_direction_batch_capability(typeof(atm)) isa
        ExtractedScreenDirectionBatchCapability
end
