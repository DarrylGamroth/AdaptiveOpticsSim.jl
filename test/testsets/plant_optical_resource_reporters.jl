struct OpticalResourceFactDirectPathModel end

Plant.plant_model_definition_style(
    ::Type{OpticalResourceFactDirectPathModel},
) = Plant.ColdPlantModelDefinition()

@inline cpu_array_bytes(::Type{T}, dimensions::Vararg{Int}) where {T} =
    UInt64(sizeof(T) * prod(dimensions))

@inline cpu_pupil_bytes(::Type{T}, n::Int) where {T} =
    cpu_array_bytes(Bool, n, n) + 2 * cpu_array_bytes(T, n, n)

@inline cpu_kolmogorov_bytes(::Type{T}, n::Int) where {T} =
    8 * cpu_array_bytes(T, n, n) + cpu_array_bytes(T, n)

@testset "Dispatch-owned optical structural resource reporters" begin
    target = HostComputeDevice()
    T = Float64
    n = 4
    telescope = Telescope(resolution=n, diameter=T(2), T=T)
    telescope_fact = structural_resource_fact(telescope,
        StructuralResourceOwnerID(:telescope, :primary), target)
    @test structural_resident_bytes(telescope_fact) ==
        cpu_array_bytes(Bool, n, n) + cpu_array_bytes(T, n, n)
    @test structural_workspace_bytes(telescope_fact) == UInt64(0)

    pupil = PupilFunction(telescope; T=T)
    pupil_fact = structural_resource_fact(pupil,
        StructuralResourceOwnerID(:path_input, :science), target)
    @test structural_resident_bytes(pupil_fact) == cpu_pupil_bytes(T, n)
    @test structural_workspace_bytes(pupil_fact) == UInt64(0)

    atmosphere = KolmogorovAtmosphere(telescope; r0=T(0.2), L0=T(25),
        T=T)
    atmosphere_fact = structural_resource_fact(atmosphere,
        StructuralResourceOwnerID(:atmosphere, :primary), target)
    @test structural_resident_bytes(atmosphere_fact) ==
        cpu_kolmogorov_bytes(T, n)
    @test structural_workspace_bytes(atmosphere_fact) == UInt64(0)

    multilayer = MultiLayerAtmosphere(telescope; r0=T(0.2), L0=T(25),
        fractional_cn2=T[0.6, 0.4], wind_speed=T[3, 5],
        wind_direction=T[0, 90], altitude=T[0, 3_000],
        layer_ids=(:ground, :high), T=T)
    screen_n = 3 * n
    layer_bytes = cpu_array_bytes(Bool, screen_n, screen_n) +
        cpu_array_bytes(T, screen_n, screen_n) +
        cpu_kolmogorov_bytes(T, screen_n)
    multilayer_fact = structural_resource_fact(multilayer,
        StructuralResourceOwnerID(:atmosphere, :multilayer), target)
    @test structural_resident_bytes(multilayer_fact) == 2 * layer_bytes
    @test structural_workspace_bytes(multilayer_fact) == UInt64(0)

    source = Source(band=:resource_fact, wavelength=T(0.8e-6),
        photon_irradiance=one(T), T=T)
    materialization = prepare_pupil_opd_materialization(atmosphere,
        telescope, source, pupil)
    materialization_fact = structural_resource_fact(materialization,
        StructuralResourceOwnerID(:path_materialization, :science), target)
    @test structural_resident_bytes(materialization_fact) ==
        cpu_array_bytes(Bool, n, n)
    @test structural_workspace_bytes(materialization_fact) == UInt64(0)

    direct = prepare_direct_imaging(pupil, source; zero_padding=2)
    padded_n = 2 * n
    direct_fact = structural_resource_fact(direct,
        StructuralResourceOwnerID(:direct_science, :leaf), target)
    expected_direct_resident = cpu_pupil_bytes(T, n) +
        cpu_array_bytes(Complex{T}, padded_n, padded_n) +
        cpu_array_bytes(T, padded_n, padded_n)
    expected_direct_workspace =
        cpu_array_bytes(Complex{T}, padded_n, padded_n) +
        cpu_array_bytes(T, padded_n, padded_n)
    @test structural_resident_bytes(direct_fact) == expected_direct_resident
    @test structural_workspace_bytes(direct_fact) == expected_direct_workspace

    second_pupil = PupilFunction(telescope; T=T)
    second_source = Source(band=:resource_fact,
        wavelength=T(0.8e-6), photon_irradiance=T(2), T=T)
    batch = prepare_direct_imaging_batch([pupil, second_pupil],
        [source, second_source]; zero_padding=2)
    batch_fact = structural_resource_fact(batch,
        StructuralResourceOwnerID(:direct_science, :batch_workspace), target)
    expected_batch_workspace =
        cpu_array_bytes(Complex{T}, padded_n, padded_n, 2) +
        cpu_array_bytes(T, padded_n, padded_n, 2) +
        2 * cpu_array_bytes(Int, 2)
    @test structural_resident_bytes(batch_fact) == UInt64(0)
    @test structural_workspace_bytes(batch_fact) == expected_batch_workspace

    definition = OpticalPathDefinition(:science, source,
        OpticalResourceFactDirectPathModel())
    context = Backends._prepare_device_execution_context(direct.output.values)
    path = PreparedPathExecutor(definition, source, telescope, atmosphere,
        pupil, direct.output, direct;
        context=context,
        materialization=materialization,
        optical_model=:direct_science,
        propagation_model=:fraunhofer_fft)
    path_fact = structural_resource_fact(path,
        StructuralResourceOwnerID(:path, :science), target)
    @test structural_resident_bytes(path_fact) == expected_direct_resident +
        cpu_array_bytes(Bool, n, n)
    @test structural_workspace_bytes(path_fact) == expected_direct_workspace
end
