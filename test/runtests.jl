using Test
using AdaptiveOpticsSim
using LinearAlgebra
using Random
using SpecialFunctions
using Statistics
using Tables
using TOML

include(joinpath(dirname(@__DIR__), "examples", "support", "subaru_ao188_simulation.jl"))
include(joinpath(dirname(@__DIR__), "examples", "support", "subaru_ao3k_simulation.jl"))
using .SubaruAO188Simulation
using .SubaruAO3kSimulation

include("reference_harness.jl")
include("ka_cpu_matrix.jl")
include("tomography.jl")

function run_tutorial_example(name::AbstractString)
    path = joinpath(dirname(@__DIR__), "examples", "tutorials", name)
    mod = Module(Symbol("Tutorial_", replace(basename(name), "." => "_")))
    Core.eval(mod, :(include(path::AbstractString) = Base.include($mod, path)))
    Base.include(mod, path)
    main_fn = Core.eval(mod, :main)
    return Base.invokelatest(main_fn)
end

function assert_source_interface(src)
    @test hasmethod(wavelength, Tuple{typeof(src)})
end

function assert_atmosphere_interface(atm, tel)
    @test applicable(advance!, atm, tel)
    @test applicable(propagate!, atm, tel)
end

function assert_atmosphere_layer_interface(layer, tel, rng, src)
    @test layer isa AdaptiveOpticsSim.AbstractAtmosphereLayer
    @test applicable(AdaptiveOpticsSim.sample_layer!, similar(tel.state.opd), layer, tel, rng)
    @test applicable(AdaptiveOpticsSim.sample_layer_accumulate!, similar(tel.state.opd), layer, tel, rng)
    @test applicable(AdaptiveOpticsSim.render_layer!, similar(tel.state.opd), layer, zero(eltype(tel.state.opd)),
        zero(eltype(tel.state.opd)), one(eltype(tel.state.opd)))
    @test applicable(AdaptiveOpticsSim.render_layer_accumulate!, similar(tel.state.opd), layer, zero(eltype(tel.state.opd)),
        zero(eltype(tel.state.opd)), one(eltype(tel.state.opd)))
    altitude = AdaptiveOpticsSim.layer_altitude(layer)
    shift_x, shift_y, footprint_scale = AdaptiveOpticsSim.layer_source_geometry(src, altitude, tel, eltype(tel.state.opd))
    sample = similar(tel.state.opd)
    fill!(sample, zero(eltype(sample)))
    AdaptiveOpticsSim.sample_layer!(sample, layer, tel, rng)
    @test size(sample) == size(tel.state.opd)
    fill!(sample, zero(eltype(sample)))
    AdaptiveOpticsSim.sample_layer_accumulate!(sample, layer, tel, rng)
    @test size(sample) == size(tel.state.opd)
    fill!(sample, zero(eltype(sample)))
    AdaptiveOpticsSim.render_layer!(sample, layer, shift_x, shift_y, footprint_scale)
    @test size(sample) == size(tel.state.opd)
    fill!(sample, zero(eltype(sample)))
    AdaptiveOpticsSim.render_layer_accumulate!(sample, layer, shift_x, shift_y, footprint_scale)
    @test size(sample) == size(tel.state.opd)
end

function assert_wfs_interface(wfs, tel)
    @test applicable(update_valid_mask!, wfs, tel)
    @test applicable(measure!, wfs, tel)
end

function assert_detector_interface(det, psf)
    @test applicable(capture!, det, psf)
end

function assert_dm_interface(dm, tel)
    @test applicable(build_influence_functions!, dm, tel)
    @test applicable(apply!, dm, tel, DMAdditive())
end

function assert_optical_element_interface(element, tel)
    @test applicable(apply!, element, tel, DMAdditive())
    @test applicable(apply!, element, tel, DMReplace())
end

function assert_reconstructor_interface(recon, slopes, expected_length::Int)
    @test recon isa AbstractReconstructorOperator
    out = zeros(eltype(slopes), expected_length)
    @test applicable(reconstruct!, out, recon, slopes)
    reconstruct!(out, recon, slopes)
    @test length(out) == expected_length
    allocated = reconstruct(recon, slopes)
    @test length(allocated) == expected_length
end

function assert_controller_interface(ctrl, input, dt::Real)
    @test applicable(update!, ctrl, input, dt)
    output = update!(ctrl, input, dt)
    @test length(output) == length(input)
end

function assert_control_simulation_interface(sim)
    @test sim isa AbstractControlSimulation
    @test applicable(step!, sim)
    @test applicable(simulation_interface, sim)
    @test applicable(simulation_readout, sim)
    iface = simulation_interface(sim)
    readout = simulation_readout(sim)
    @test simulation_command(readout) === simulation_command(sim)
    @test simulation_slopes(readout) === simulation_slopes(sim)
    @test simulation_command(simulation_readout(iface)) === simulation_command(iface)
    return iface
end

function assert_interaction_matrix_contract(imat, expected_rows::Int, expected_cols::Int, amplitude::Real)
    @test imat isa InteractionMatrix
    @test size(imat.matrix) == (expected_rows, expected_cols)
    @test imat.amplitude ≈ amplitude
end

function assert_calibration_vault_contract(vault, forward::AbstractMatrix; inverted::Bool=true)
    @test vault isa CalibrationVault
    @test vault.D === forward
    if inverted
        @test !isnothing(vault.M)
        @test size(vault.M, 2) == size(forward, 1)
        @test length(vault.singular_values) == min(size(forward)...)
        @test vault.effective_rank >= 0
    else
        @test isnothing(vault.M)
        @test isempty(vault.singular_values)
    end
end

function assert_modal_basis_contract(basis::ModalBasis, n_commands::Int, n_modes::Int)
    @test size(basis.M2C) == (n_commands, n_modes)
    @test size(basis.basis, 2) == n_modes
    if !isnothing(basis.projector)
        @test size(basis.projector, 2) == size(basis.basis, 1)
        @test size(basis.projector, 1) == n_modes
    end
end

function assert_ao_calibration_contract(calib::AOCalibration, n_commands::Int, n_modes::Int)
    @test size(calib.M2C) == (n_commands, n_modes)
    @test size(calib.basis, 2) == n_modes
    @test calib.calibration isa CalibrationVault
end

function assert_meta_sensitivity_contract(meta::MetaSensitivity, n_fields::Int)
    @test meta.calib0 isa CalibrationVault
    @test meta.meta isa CalibrationVault
    @test length(meta.field_order) == n_fields
    @test size(meta.meta.D, 2) == n_fields
    @test meta.meta.M !== nothing
end

function subharmonic_tiptilt_power(phs::AbstractMatrix)
    n = size(phs, 1)
    coords = collect(LinRange(-1.0, 1.0, n))
    x = repeat(reshape(coords, 1, n), n, 1)
    y = repeat(reshape(coords, n, 1), 1, n)
    A = hcat(vec(x), vec(y), ones(n * n))
    coeffs = A \ vec(phs)
    return coeffs[1]^2 + coeffs[2]^2
end

function subharmonic_metrics(atm::KolmogorovAtmosphere, D::Real;
    n::Int=24, nsamp::Int=8, kwargs...)
    delta = D / n
    shift = n ÷ 2
    ws = PhaseStatsWorkspace(n; T=Float64)
    tiptilt = Float64[]
    structure = Float64[]
    for seed in 1:nsamp
        phs = ft_sh_phase_screen(atm, n, delta; rng=MersenneTwister(seed), ws=ws, kwargs...)
        push!(tiptilt, subharmonic_tiptilt_power(phs))
        diff = @views phs[:, 1:end-shift] .- phs[:, 1+shift:end]
        push!(structure, Statistics.mean(abs2, diff))
    end
    return (; tiptilt=Statistics.mean(tiptilt), structure=Statistics.mean(structure))
end

function moving_atmosphere_trace(;
    seed::Integer=1,
    steps::Integer=4,
    resolution::Int=16,
    diameter::Real=8.0,
    sampling_time::Real=1e-3,
    r0::Real=0.2,
    L0::Real=25.0,
    fractional_cn2::AbstractVector=[1.0],
    wind_speed::AbstractVector=[0.0],
    wind_direction::AbstractVector=[0.0],
    altitude::AbstractVector=[0.0],
)
    tel = Telescope(resolution=resolution, diameter=diameter, sampling_time=sampling_time, central_obstruction=0.0)
    atm = MultiLayerAtmosphere(tel;
        r0=r0,
        L0=L0,
        fractional_cn2=fractional_cn2,
        wind_speed=wind_speed,
        wind_direction=wind_direction,
        altitude=altitude,
    )
    rng = MersenneTwister(seed)
    trace = Matrix{Float64}[]
    for _ in 1:steps
        advance!(atm, tel; rng=rng)
        push!(trace, copy(atm.state.opd))
    end
    return trace
end

function moving_wfs_slope_trace(;
    seed::Integer=1,
    steps::Integer=4,
)
    tel = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    src = Source(band=:I, magnitude=0.0)
    delta = tel.params.diameter / tel.params.resolution
    atm = MultiLayerAtmosphere(tel;
        r0=0.2,
        L0=25.0,
        fractional_cn2=[0.7, 0.3],
        wind_speed=[delta / tel.params.sampling_time, 0.5 * delta / tel.params.sampling_time],
        wind_direction=[0.0, 90.0],
        altitude=[0.0, 5000.0],
    )
    wfs = ShackHartmann(tel; n_subap=4)
    rng = MersenneTwister(seed)
    trace = Vector{Vector{Float64}}(undef, steps)
    for i in 1:steps
        advance!(atm, tel; rng=rng)
        propagate!(atm, tel)
        measure!(wfs, tel, src)
        trace[i] = copy(wfs.state.slopes)
    end
    return trace
end

function moving_closed_loop_trace(;
    seed::Integer=1,
    steps::Integer=5,
)
    tel = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    src = Source(band=:I, magnitude=0.0)
    delta = tel.params.diameter / tel.params.resolution
    atm = MultiLayerAtmosphere(tel;
        r0=0.2,
        L0=25.0,
        fractional_cn2=[0.7, 0.3],
        wind_speed=[delta / tel.params.sampling_time, 0.5 * delta / tel.params.sampling_time],
        wind_direction=[0.0, 90.0],
        altitude=[0.0, 5000.0],
    )
    dm = DeformableMirror(tel; n_act=4, influence_width=0.3)
    wfs = ZernikeWFS(tel; n_subap=4, diffraction_padding=2)
    sim = AOSimulation(tel, atm, src, dm, wfs)
    imat = interaction_matrix(dm, wfs, tel, src; amplitude=1e-8)
    recon = ModalReconstructor(imat; gain=0.2)
    wfs_det = Detector(noise=NoiseNone(), integration_time=1.0, qe=1.0, binning=1)
    science_det = Detector(noise=NoiseNone(), integration_time=1.0, qe=1.0, binning=1)
    runtime = ClosedLoopRuntime(sim, recon; rng=MersenneTwister(seed), wfs_detector=wfs_det, science_detector=science_det)
    prepare!(runtime)
    delay = VectorDelayLine(runtime.command, 1)

    slope_norms = zeros(Float64, steps)
    command_norms = zeros(Float64, steps)
    wfs_energy = zeros(Float64, steps)
    science_energy = zeros(Float64, steps)

    for i in 1:steps
        sense!(runtime)
        reconstruct!(runtime)
        delayed = copy(shift_delay!(delay, runtime.command))
        copyto!(runtime.command, delayed)
        AdaptiveOpticsSim.apply_runtime_command!(runtime)
        slope_norms[i] = norm(simulation_slopes(runtime))
        command_norms[i] = norm(simulation_command(runtime))
        wfs_energy[i] = sum(abs, simulation_wfs_frame(runtime))
        science_energy[i] = sum(abs, simulation_science_frame(runtime))
    end

    return (; slope_norms, command_norms, wfs_energy, science_energy)
end

@test AdaptiveOpticsSim.PROJECT_STATUS == :in_development

@testset "GPU backend registry" begin
    @test !gpu_backend_loaded(CUDABackendTag)
    @test !gpu_backend_loaded(MetalBackendTag)
    @test !gpu_backend_loaded(AMDGPUBackendTag)
    @test gpu_backend_array_type(CUDABackendTag) === nothing
    @test gpu_backend_array_type(MetalBackendTag) === nothing
    @test gpu_backend_array_type(AMDGPUBackendTag) === nothing
    @test gpu_backend_name(Matrix{Float64}) === nothing
    @test available_gpu_backends() == ()
    @test GPUArrayBuildBackend(CUDABackendTag) isa GPUArrayBuildBackend{CUDABackendTag}
end

@testset "Telescope and PSF" begin
    tel = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3, central_obstruction=0.2)
    src = Source(band=:I, magnitude=0.0)
    field = ElectricField(tel, src; zero_padding=2)
    @test size(field.state.field) == (64, 64)
    @test field.params.wavelength == wavelength(src)
    centered_psf = AdaptiveOpticsSim.centered_psf_from_field!(similar(field.state.intensity), field)
    psf = compute_psf!(tel, src; zero_padding=2)
    @test size(psf) == (64, 64)
    @test maximum(psf) > 0
    @test isfinite(sum(psf))
    @test centered_psf ≈ psf
    @test tel.state.psf_workspace !== nothing
    cached_ws = tel.state.psf_workspace
    compute_psf!(tel, src; zero_padding=2)
    @test tel.state.psf_workspace === cached_ws

    tel_dim = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3, central_obstruction=0.2,
        pupil_reflectivity=0.25)
    psf_dim = compute_psf!(tel_dim, src; zero_padding=2)
    @test sum(psf_dim) ≈ 0.25 * sum(psf)
    fmap = flux_map(tel_dim, src)
    @test size(fmap) == size(tel_dim.state.pupil)
    @test maximum(fmap) > 0
    @test optical_path(src, tel_dim) == "source(I) -> telescope"

    tel_simple = Telescope(resolution=8, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    src_simple = Source(band=:I, magnitude=0.0)
    field_simple = ElectricField(tel_simple, src_simple; zero_padding=1)
    direct = similar(field_simple.state.field)
    fill_telescope_field!(direct, tel_simple, src_simple; zero_padding=1)
    @test direct == field_simple.state.field

    phase_map = fill(pi / 2, tel_simple.params.resolution, tel_simple.params.resolution)
    baseline = copy(field_simple.state.field)
    apply_phase!(field_simple, phase_map; units=:phase)
    @test field_simple.state.field ≈ baseline .* cis.(phase_map)

    field_simple = ElectricField(tel_simple, src_simple; zero_padding=1)
    opd_map = fill(eltype(tel_simple.state.opd)(wavelength(src_simple) / 4), tel_simple.params.resolution, tel_simple.params.resolution)
    baseline = copy(field_simple.state.field)
    apply_phase!(field_simple, opd_map; units=:opd)
    @test field_simple.state.field ≈ baseline .* cispi.(2 .* opd_map ./ wavelength(src_simple))

    field_simple = ElectricField(tel_simple, src_simple; zero_padding=1)
    amplitude_map = fill(0.5, tel_simple.params.resolution, tel_simple.params.resolution)
    baseline = copy(field_simple.state.field)
    apply_amplitude!(field_simple, amplitude_map)
    @test field_simple.state.field ≈ baseline .* amplitude_map

    intensity_buffer = similar(field_simple.state.intensity)
    intensity!(intensity_buffer, field_simple)
    @test intensity_buffer ≈ abs2.(field_simple.state.field)
    @test intensity!(field_simple) === field_simple.state.intensity

    fraunhofer = FraunhoferPropagation(field)
    propagated = similar(field.state.field)
    propagate_field!(propagated, field, fraunhofer)
    propagated_psf = similar(field.state.intensity)
    @. propagated_psf = abs2(propagated)
    @test propagated_psf ≈ centered_psf atol=1e-10 rtol=1e-10
    @test fraunhofer.params.output_sampling_rad ≈ wavelength(src) / (field.params.padded_resolution * field.params.sampling_m)

    fresnel = FresnelPropagation(field; distance_m=25.0)
    fresnel_out = similar(field.state.field)
    propagate_field!(fresnel_out, field, fresnel)
    @test size(fresnel_out) == size(field.state.field)
    @test isfinite(sum(abs, fresnel_out))
    reverse_fresnel = FresnelPropagation(field; distance_m=-25.0)
    reverse_input = ElectricField(tel, src; zero_padding=2)
    copyto!(reverse_input.state.field, fresnel_out)
    propagate_field!(reverse_input, reverse_fresnel)
    @test reverse_input.state.field ≈ field.state.field atol=1e-9 rtol=1e-9

    zero_fresnel = FresnelPropagation(field; distance_m=0.0)
    zero_out = similar(field.state.field)
    propagate_field!(zero_out, field, zero_fresnel)
    @test zero_out ≈ field.state.field atol=1e-10 rtol=1e-10

    atm_tel = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    atm_src = Source(band=:I, magnitude=0.0)
    atm = MultiLayerAtmosphere(atm_tel;
        r0=0.2,
        L0=25.0,
        fractional_cn2=[0.7, 0.3],
        wind_speed=[8.0, 4.0],
        wind_direction=[0.0, 90.0],
        altitude=[0.0, 5000.0],
    )
    advance!(atm, atm_tel; rng=MersenneTwister(1))
    geom_prop = AtmosphericFieldPropagation(atm, atm_tel, atm_src;
        model=GeometricAtmosphericPropagation(T=Float64),
        zero_padding=1,
        T=Float64)
    geom_field = propagate_atmosphere_field!(geom_prop, atm, atm_tel, atm_src)
    tel_geom = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    propagate!(atm, tel_geom, atm_src)
    collapsed = ElectricField(tel_geom, atm_src; zero_padding=1)
    @test geom_field.state.field ≈ collapsed.state.field atol=1e-8 rtol=1e-8

    fresnel_atm = MultiLayerAtmosphere(atm_tel;
        r0=0.2,
        L0=25.0,
        fractional_cn2=[1.0],
        wind_speed=[0.0],
        wind_direction=[0.0],
        altitude=[0.0],
    )
    advance!(fresnel_atm, atm_tel; rng=MersenneTwister(2))
    fresnel_prop = AtmosphericFieldPropagation(fresnel_atm, atm_tel, atm_src;
        model=LayeredFresnelAtmosphericPropagation(T=Float64),
        zero_padding=1,
        T=Float64)
    fresnel_field = propagate_atmosphere_field!(fresnel_prop, fresnel_atm, atm_tel, atm_src)
    geom_single = AtmosphericFieldPropagation(fresnel_atm, atm_tel, atm_src;
        model=GeometricAtmosphericPropagation(T=Float64),
        zero_padding=1,
        T=Float64)
    geom_single_field = propagate_atmosphere_field!(geom_single, fresnel_atm, atm_tel, atm_src)
    @test fresnel_field.state.field ≈ geom_single_field.state.field atol=1e-8 rtol=1e-8

    spectral = with_spectrum(atm_src, SpectralBundle([wavelength(atm_src), 1.1 * wavelength(atm_src)], [0.75, 0.25]))
    spectral_prop = AtmosphericFieldPropagation(atm, atm_tel, spectral;
        model=GeometricAtmosphericPropagation(T=Float64),
        zero_padding=2,
        T=Float64)
    spectral_intensity = atmospheric_intensity!(spectral_prop, atm, atm_tel, spectral)
    mono_prop = AtmosphericFieldPropagation(atm, atm_tel, atm_src;
        model=GeometricAtmosphericPropagation(T=Float64),
        zero_padding=2,
        T=Float64)
    mono_intensity = atmospheric_intensity!(mono_prop, atm, atm_tel, atm_src)
    single_spectral = with_spectrum(atm_src, SpectralBundle([wavelength(atm_src)], [1.0]))
    single_spectral_prop = AtmosphericFieldPropagation(atm, atm_tel, single_spectral;
        model=GeometricAtmosphericPropagation(T=Float64),
        zero_padding=2,
        T=Float64)
    @test atmospheric_intensity!(single_spectral_prop, atm, atm_tel, single_spectral) ≈ mono_intensity atol=1e-8 rtol=1e-8
    @test size(spectral_intensity) == size(mono_intensity)
    @test sum(spectral_intensity) > 0
end

@testset "Aperture masks" begin
    bool_mask = falses(32, 32)
    build_mask!(bool_mask, CircularAperture(radius=1.0))
    @test bool_mask[16, 16]
    @test !bool_mask[1, 1]

    annulus = falses(32, 32)
    build_mask!(annulus, AnnularAperture(inner_radius=0.25, outer_radius=1.0))
    @test !annulus[16, 16]
    @test count(annulus) < count(bool_mask)

    weighted = fill(-1.0, 16, 16)
    build_mask!(weighted, RectangularROI(5:8, 6:10); inside=2.0, outside=-1.0)
    @test all(weighted[5:8, 6:10] .== 2.0)
    @test weighted[4, 6] == -1.0

    spider = trues(32, 32)
    apply_mask!(spider, SpiderMask(thickness=0.08, angle_rad=pi / 2))
    @test count(spider) < length(spider)

    tel = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3, central_obstruction=0.25)
    expected = falses(32, 32)
    build_mask!(expected, AnnularAperture(inner_radius=0.25, outer_radius=1.0))
    @test tel.state.pupil == expected
    apply_spiders!(tel; thickness=0.4, angles=[0.0, 90.0])
    manual = copy(expected)
    apply_mask!(manual, SpiderMask(thickness=0.1, angle_rad=0.0))
    apply_mask!(manual, SpiderMask(thickness=0.1, angle_rad=pi / 2))
    @test tel.state.pupil == manual

    sf_tel = Telescope(resolution=8, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    sf = SpatialFilter(sf_tel; shape=SquareFilter(), diameter=4, zero_padding=2)
    @test count(x -> !iszero(x), sf.state.mask) > 0
    foucault = SpatialFilter(sf_tel; shape=FoucaultFilter(), diameter=4, zero_padding=2)
    foucault_count = count(x -> !iszero(x), foucault.state.mask)
    @test 0 < foucault_count < length(foucault.state.mask)

    pupil = falses(8, 8)
    pupil[1:4, 1:4] .= true
    pupil[5:8, 5:8] .= true
    valid = falses(2, 2)
    build_mask!(valid, SubapertureGridMask(threshold=0.5), pupil)
    @test valid == Bool[true false; false true]
    valid2 = similar(valid)
    AdaptiveOpticsSim.set_valid_subapertures!(valid2, pupil, 0.5)
    @test valid2 == valid
end

@testset "Zernike basis" begin
    tel = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    zb = ZernikeBasis(tel, 5)
    compute_zernike!(zb, tel)
    @test size(zb.modes) == (32, 32, 5)
    @test noll_to_nm(1) == (0, 0)
    @test noll_to_nm(2) == (1, -1)
    @test noll_to_nm(3) == (1, 1)
    @test noll_to_nm(4) == (2, -2)
end

@testset "Atmosphere propagation" begin
    tel = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    atm = KolmogorovAtmosphere(tel; r0=0.2, L0=25.0)
    advance!(atm, tel; rng=MersenneTwister(1))
    propagate!(atm, tel)
    @test size(atm.state.opd) == (32, 32)
    @test sum(abs.(tel.state.opd)) > 0

    delta = tel.params.diameter / tel.params.resolution
    ensure_psd!(atm, delta)
    psd_snapshot = copy(atm.state.psd)
    ensure_psd!(atm, delta)
    @test psd_snapshot == atm.state.psd

    atm = KolmogorovAtmosphere(tel; r0=0.1, L0=25.0)
    ensure_psd!(atm, delta)
    @test psd_snapshot != atm.state.psd
end

@testset "Infinite atmosphere boundary math" begin
    tel = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    @test AdaptiveOpticsSim.default_infinite_screen_resolution(tel.params.resolution) == 3 * tel.params.resolution
    @test AdaptiveOpticsSim.default_infinite_stencil_size(tel.params.resolution) >= 257

    stencil = AdaptiveOpticsSim.infinite_boundary_stencil(4, 0.25;
        stencil_size=9,
        orientation=:column,
        side=:positive,
        tail_stride=4,
        T=Float64,
    )
    @test size(stencil.boundary_coords) == (9, 2)
    @test size(stencil.stencil_positions, 2) == 2
    @test all(stencil.boundary_coords[:, 1] .== 0)

    row_stencil = AdaptiveOpticsSim.infinite_boundary_stencil(4, 0.25;
        stencil_size=9,
        orientation=:row,
        side=:negative,
        tail_stride=4,
        T=Float64,
    )
    @test all(row_stencil.boundary_coords[:, 2] .== 10)

    op = AdaptiveOpticsSim.boundary_injection_operator(stencil, 0.2, 25.0)
    @test size(op.predictor, 1) == size(stencil.boundary_coords, 1)
    @test size(op.predictor, 2) == size(stencil.stencil_coords, 1)
    @test size(op.residual_factor, 1) == size(stencil.boundary_coords, 1)
    @test isfinite(op.condition_ratio)
    @test op.condition_ratio > 0
    @test isapprox(op.residual_covariance, op.residual_factor * transpose(op.residual_factor); rtol=1e-10, atol=1e-10)

    rng = MersenneTwister(12)
    stencil_data = randn(rng, size(op.predictor, 2))
    nsamp = 4000
    samples = Matrix{Float64}(undef, size(op.predictor, 1), nsamp)
    expected_mean = op.predictor * stencil_data
    for i in 1:nsamp
        samples[:, i] = AdaptiveOpticsSim.sample_boundary_line(op, stencil_data, rng)
    end
    empirical_mean = vec(mean(samples; dims=2))
    centered = samples .- empirical_mean
    empirical_cov = centered * transpose(centered) / (nsamp - 1)
    @test isapprox(empirical_mean, expected_mean; rtol=0.08, atol=0.08)
    @test isapprox(empirical_cov, op.residual_covariance; rtol=0.12, atol=0.12)

    screen = InfinitePhaseScreen(tel; r0=0.2, L0=25.0, screen_resolution=33, stencil_size=35)
    @test size(screen.state.screen) == (35, 35)
    @test size(screen.state.extract_buffer) == (32, 32)
    @test size(screen.state.column_positive.operator.predictor, 1) == 35
    @test size(screen.state.row_negative.operator.predictor, 1) == 35

    atm = InfiniteMultiLayerAtmosphere(tel;
        r0=0.2,
        L0=25.0,
        fractional_cn2=[1.0],
        wind_speed=[0.0],
        wind_direction=[0.0],
        altitude=[0.0],
        screen_resolution=33,
        stencil_size=35,
    )
    @test length(atm.layers) == 1
    @test size(atm.state.opd) == (32, 32)
    advance!(atm, tel; rng=MersenneTwister(1))
    propagate!(atm, tel)
    @test size(atm.layers[1].screen.state.screen) == (35, 35)
    @test atm.layers[1].screen.state.initialized
    @test atm.layers[1].state.integer_shift_x == 0
    @test atm.layers[1].state.integer_shift_y == 0
    @test tel.state.opd == atm.state.opd .* tel.state.pupil

    @test_throws InvalidConfiguration AdaptiveOpticsSim.infinite_boundary_stencil(4, 0.25;
        stencil_size=8,
        orientation=:column,
        side=:positive,
    )
    @test_throws NumericalConditionError AdaptiveOpticsSim.boundary_injection_operator(stencil, 0.2, 25.0; conditioning_tol=1.0)
end

@testset "Infinite atmosphere stepping regressions" begin
    tel = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)

    normalized_correlation(a::AbstractMatrix, b::AbstractMatrix) =
        dot(vec(a), vec(b)) / sqrt(dot(vec(a), vec(a)) * dot(vec(b), vec(b)))

    stationary = InfiniteMultiLayerAtmosphere(tel;
        r0=0.2,
        L0=25.0,
        fractional_cn2=[1.0],
        wind_speed=[0.0],
        wind_direction=[0.0],
        altitude=[0.0],
        screen_resolution=33,
        stencil_size=35,
    )
    advance!(stationary, tel; rng=MersenneTwister(4))
    stationary_snapshot = copy(stationary.state.opd)
    advance!(stationary, tel; rng=MersenneTwister(5))
    @test stationary.state.opd == stationary_snapshot
    @test stationary.layers[1].state.integer_shift_x == 0
    @test stationary.layers[1].state.integer_shift_y == 0

    small_step = InfiniteMultiLayerAtmosphere(tel;
        r0=0.2,
        L0=25.0,
        fractional_cn2=[1.0],
        wind_speed=[50.0],
        wind_direction=[0.0],
        altitude=[0.0],
        screen_resolution=33,
        stencil_size=35,
    )
    rng = MersenneTwister(6)
    advance!(small_step, tel; rng=rng)
    opd_1 = copy(small_step.state.opd)
    advance!(small_step, tel; rng=rng)
    opd_2 = copy(small_step.state.opd)
    corr = normalized_correlation(opd_1, opd_2)
    @test corr > 0.95
    @test small_step.layers[1].state.integer_shift_x == 0

    row_step = InfiniteMultiLayerAtmosphere(tel;
        r0=0.2,
        L0=25.0,
        fractional_cn2=[1.0],
        wind_speed=[500.0],
        wind_direction=[90.0],
        altitude=[0.0],
        screen_resolution=33,
        stencil_size=35,
    )
    advance!(row_step, tel; rng=MersenneTwister(7))
    @test row_step.layers[1].state.integer_shift_y == 1
    @test row_step.layers[1].state.integer_shift_x == 0

    diagonal_speed = 250.0 * sqrt(2)
    diagonal_step = InfiniteMultiLayerAtmosphere(tel;
        r0=0.2,
        L0=25.0,
        fractional_cn2=[1.0],
        wind_speed=[diagonal_speed],
        wind_direction=[45.0],
        altitude=[0.0],
        screen_resolution=33,
        stencil_size=35,
    )
    diagonal_rng = MersenneTwister(71)
    advance!(diagonal_step, tel; rng=diagonal_rng)
    advance!(diagonal_step, tel; rng=diagonal_rng)
    @test diagonal_step.layers[1].state.integer_shift_x == 1
    @test diagonal_step.layers[1].state.integer_shift_y == 1
    @test abs(diagonal_step.layers[1].state.offset_x) < 1e-10
    @test abs(diagonal_step.layers[1].state.offset_y) < 1e-10

    function infinite_trace(; seed::Integer=73, steps::Int=10)
        tel_local = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
        src_local = Source(band=:I, magnitude=0.0, coordinates=(15.0, 30.0))
        atm = InfiniteMultiLayerAtmosphere(tel_local;
            r0=0.2,
            L0=25.0,
            fractional_cn2=[0.7, 0.3],
            wind_speed=[6.0, 3.0],
            wind_direction=[0.0, 120.0],
            altitude=[0.0, 5000.0],
            screen_resolution=33,
            stencil_size=35,
        )
        rng = MersenneTwister(seed)
        opd_trace = Matrix{Float64}[]
        for _ in 1:steps
            advance!(atm, tel_local; rng=rng)
            propagate!(atm, tel_local, src_local)
            push!(opd_trace, copy(Array(tel_local.state.opd)))
        end
        return (; opd_trace, screen=copy(Array(atm.layers[1].screen.state.screen)))
    end

    trace_a = infinite_trace()
    trace_b = infinite_trace()
    @test trace_a.opd_trace == trace_b.opd_trace
    @test trace_a.screen == trace_b.screen

    function trajectory_std_windows(; seed::Integer=79, steps::Int=24)
        tel_local = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
        atm = InfiniteMultiLayerAtmosphere(tel_local;
            r0=0.2,
            L0=25.0,
            fractional_cn2=[0.5, 0.5],
            wind_speed=[10.0, 6.0],
            wind_direction=[0.0, 90.0],
            altitude=[0.0, 5000.0],
            screen_resolution=33,
            stencil_size=35,
        )
        rng = MersenneTwister(seed)
        stds = Float64[]
        for _ in 1:steps
            advance!(atm, tel_local; rng=rng)
            push!(stds, std(vec(Array(atm.state.opd))))
        end
        midpoint = steps ÷ 2
        return mean(@view(stds[1:midpoint])), mean(@view(stds[midpoint + 1:end]))
    end

    early_std, late_std = trajectory_std_windows()
    @test isapprox(late_std, early_std; rtol=0.2)

    wind_speed_px = tel.params.diameter / tel.params.resolution / tel.params.sampling_time
    finite = MultiLayerAtmosphere(tel;
        r0=0.2,
        L0=25.0,
        fractional_cn2=[1.0],
        wind_speed=[wind_speed_px],
        wind_direction=[0.0],
        altitude=[0.0],
    )
    infinite = InfiniteMultiLayerAtmosphere(tel;
        r0=0.2,
        L0=25.0,
        fractional_cn2=[1.0],
        wind_speed=[wind_speed_px],
        wind_direction=[0.0],
        altitude=[0.0],
        screen_resolution=33,
        stencil_size=35,
    )
    finite_rng = MersenneTwister(8)
    infinite_rng = MersenneTwister(8)
    advance!(finite, tel; rng=finite_rng)
    advance!(infinite, tel; rng=infinite_rng)
    finite_snapshot = copy(finite.state.opd)
    infinite_snapshot = copy(infinite.state.opd)
    period = AdaptiveOpticsSim.moving_layer_screen_resolution(tel.params.resolution)
    for _ in 1:period
        advance!(finite, tel; rng=finite_rng)
        advance!(infinite, tel; rng=infinite_rng)
    end
    @test finite.state.opd == finite_snapshot
    @test !isapprox(infinite.state.opd, infinite_snapshot; rtol=1e-8, atol=1e-8)
    for _ in 1:(2 * period)
        advance!(finite, tel; rng=finite_rng)
        advance!(infinite, tel; rng=infinite_rng)
    end
    @test finite.state.opd == finite_snapshot
    @test normalized_correlation(infinite.state.opd, infinite_snapshot) < 0.9
end

@testset "Source-aware atmosphere extraction" begin
    tel = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    onaxis = Source(band=:I, magnitude=0.0, coordinates=(0.0, 0.0))
    offaxis = Source(band=:I, magnitude=0.0, coordinates=(100.0, 0.0))
    lgs = LGSSource(magnitude=0.0, coordinates=(100.0, 0.0), altitude=10_000.0)

    shift_x_ngs, shift_y_ngs, footprint_ngs = AdaptiveOpticsSim.layer_source_geometry(offaxis, 5000.0, tel, Float64)
    shift_x_lgs, shift_y_lgs, footprint_lgs = AdaptiveOpticsSim.layer_source_geometry(lgs, 5000.0, tel, Float64)
    @test shift_x_ngs ≈ shift_x_lgs
    @test shift_y_ngs ≈ shift_y_lgs
    @test shift_y_ngs ≈ 0.0
    @test footprint_ngs == 1.0
    @test footprint_lgs ≈ 0.5

    function fill_x_ramp!(screen::AbstractMatrix)
        @inbounds for j in axes(screen, 2), i in axes(screen, 1)
            screen[i, j] = j
        end
        return screen
    end

    function offaxis_shift_regression(constructor; seed::Integer=41, kwargs...)
        atm = constructor(tel;
            r0=0.2,
            L0=25.0,
            fractional_cn2=[1.0],
            wind_speed=[0.0],
            wind_direction=[0.0],
            altitude=[5000.0],
            kwargs...,
        )
        advance!(atm, tel; rng=MersenneTwister(seed))
        if atm isa MultiLayerAtmosphere
            fill_x_ramp!(atm.layers[1].generator.state.opd)
            atm.layers[1].state.offset_x = 0.0
            atm.layers[1].state.offset_y = 0.0
        else
            fill_x_ramp!(atm.layers[1].screen.state.screen)
            atm.layers[1].state.offset_x = 0.0
            atm.layers[1].state.offset_y = 0.0
        end
        propagate!(atm, tel, onaxis)
        onaxis_opd = copy(tel.state.opd)
        propagate!(atm, tel, offaxis)
        offaxis_opd = copy(tel.state.opd)
        propagate!(atm, tel, lgs)
        lgs_opd = copy(tel.state.opd)
        return (; onaxis_opd, offaxis_opd, lgs_opd)
    end

    finite_offaxis = offaxis_shift_regression(MultiLayerAtmosphere)
    infinite_offaxis = offaxis_shift_regression(InfiniteMultiLayerAtmosphere; screen_resolution=33, stencil_size=35)
    @test mean(finite_offaxis.offaxis_opd) > mean(finite_offaxis.onaxis_opd)
    @test mean(infinite_offaxis.offaxis_opd) > mean(infinite_offaxis.onaxis_opd)
    @test std(vec(finite_offaxis.lgs_opd)) < std(vec(finite_offaxis.offaxis_opd))
    @test std(vec(infinite_offaxis.lgs_opd)) < std(vec(infinite_offaxis.offaxis_opd))

    function onaxis_fast_path_regression(constructor; kwargs...)
        atm = constructor(tel;
            r0=0.2,
            L0=25.0,
            fractional_cn2=[0.7, 0.3],
            wind_speed=[6.0, 3.0],
            wind_direction=[0.0, 120.0],
            altitude=[0.0, 5000.0],
            kwargs...,
        )
        rng = MersenneTwister(52)
        advance!(atm, tel; rng=rng)
        propagate!(atm, tel)
        base_opd = copy(tel.state.opd)
        @test !atm.state.source_geometry.valid
        propagate!(atm, tel, onaxis)
        @test tel.state.opd == base_opd
        @test !atm.state.source_geometry.valid
        propagate!(atm, tel, offaxis)
        @test atm.state.source_geometry.valid
        @test atm.state.source_geometry.cached_x_arcsec ≈ coordinates_xy_arcsec(offaxis)[1]
        @test atm.state.source_geometry.cached_y_arcsec ≈ coordinates_xy_arcsec(offaxis)[2]
        return atm
    end

    onaxis_fast_path_regression(MultiLayerAtmosphere)
    onaxis_fast_path_regression(InfiniteMultiLayerAtmosphere; screen_resolution=33, stencil_size=35)

    function ensemble_std(constructor, fractions; nsamp::Int=12, kwargs...)
        acc = 0.0
        for s in 1:nsamp
            atm = constructor(tel;
                r0=0.2,
                L0=25.0,
                fractional_cn2=fractions,
                wind_speed=fill(0.0, length(fractions)),
                wind_direction=fill(0.0, length(fractions)),
                altitude=fill(0.0, length(fractions)),
                kwargs...,
            )
            advance!(atm, tel; rng=MersenneTwister(s))
            acc += std(vec(atm.state.opd))
        end
        return acc / nsamp
    end

    finite_single_std = ensemble_std(MultiLayerAtmosphere, [1.0])
    infinite_single_std = ensemble_std(InfiniteMultiLayerAtmosphere, [1.0]; screen_resolution=33, stencil_size=35)
    infinite_equal_std = ensemble_std(InfiniteMultiLayerAtmosphere, [0.5, 0.5]; screen_resolution=33, stencil_size=35)
    infinite_uneven_std = ensemble_std(InfiniteMultiLayerAtmosphere, [0.8, 0.2]; screen_resolution=33, stencil_size=35)
    @test isapprox(infinite_single_std, finite_single_std; rtol=0.2)
    @test isapprox(infinite_equal_std, infinite_single_std; rtol=0.2)
    @test isapprox(infinite_uneven_std, infinite_single_std; rtol=0.2)
end

@testset "Sub-harmonic phase screens" begin
    tel = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    atm = KolmogorovAtmosphere(tel; r0=0.2, L0=25.0)
    delta = tel.params.diameter / tel.params.resolution
    rng = MersenneTwister(11)
    ws = PhaseStatsWorkspace(32; T=Float64)
    phs_base = ft_sh_phase_screen(atm, 32, delta; rng=rng, ws=ws, subharmonics=false)
    rng = MersenneTwister(11)
    phs_sh = ft_sh_phase_screen(atm, 32, delta; rng=rng, ws=ws, subharmonics=true,
        mode=FidelitySubharmonics())
    @test sum(abs.(phs_sh .- phs_base)) > 0
    @test AdaptiveOpticsSim.resolve_subharmonic_levels(25.0, 8.0) == 4
    @test AdaptiveOpticsSim.resolve_subharmonic_levels(200.0, 8.0) > 4

    atm_large = KolmogorovAtmosphere(tel; r0=0.2, L0=200.0)
    rng = MersenneTwister(11)
    phs_default = ft_sh_phase_screen(atm_large, 32, delta; rng=rng, ws=ws,
        subharmonics=true, mode=FidelitySubharmonics())
    rng = MersenneTwister(11)
    phs_legacy = ft_sh_phase_screen(atm_large, 32, delta; rng=rng, ws=ws,
        subharmonics=true, mode=FastSubharmonics())
    @test sum(abs.(phs_default .- phs_legacy)) > 0

    metrics_25_none = subharmonic_metrics(atm, tel.params.diameter; subharmonics=false)
    metrics_25_legacy = subharmonic_metrics(atm, tel.params.diameter;
        subharmonics=true, mode=FastSubharmonics())
    metrics_25_default = subharmonic_metrics(atm, tel.params.diameter;
        subharmonics=true, mode=FidelitySubharmonics())
    @test metrics_25_none.tiptilt < metrics_25_legacy.tiptilt < metrics_25_default.tiptilt
    @test metrics_25_none.structure < metrics_25_legacy.structure < metrics_25_default.structure

    metrics_200_none = subharmonic_metrics(atm_large, tel.params.diameter; subharmonics=false)
    metrics_200_legacy = subharmonic_metrics(atm_large, tel.params.diameter;
        subharmonics=true, mode=FastSubharmonics())
    metrics_200_default = subharmonic_metrics(atm_large, tel.params.diameter;
        subharmonics=true, mode=FidelitySubharmonics())
@test metrics_200_none.tiptilt < metrics_200_legacy.tiptilt < metrics_200_default.tiptilt
@test metrics_200_none.structure < metrics_200_legacy.structure < metrics_200_default.structure
end

@testset "Fidelity profiles" begin
    @test default_fidelity_profile() isa ScientificProfile
    @test default_subharmonic_mode(ScientificProfile()) isa FidelitySubharmonics
    @test default_subharmonic_mode(FastProfile()) isa FastSubharmonics
    @test default_ncpa_basis(ScientificProfile()).method isa KLHHtPSD
    @test default_ncpa_basis(FastProfile()).method isa KLDMModes

    mixed = ProfileBundle(ScientificProfile(); lift=FastProfile(), tomography=FastProfile())
    @test atmosphere_profile(mixed) isa ScientificProfile
    @test calibration_profile(mixed) isa ScientificProfile
    @test detector_profile(mixed) isa ScientificProfile
    @test lift_profile(mixed) isa FastProfile
    @test tomography_profile(mixed) isa FastProfile
    @test default_subharmonic_mode(mixed) isa FidelitySubharmonics
    @test default_ncpa_basis(mixed).method isa KLHHtPSD

    fast_cal = ProfileBundle(ScientificProfile(); calibration=FastProfile())
    @test default_ncpa_basis(fast_cal).method isa KLDMModes

    tel = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    atm = KolmogorovAtmosphere(tel; r0=0.15, L0=200.0)
    scientific = ft_sh_phase_screen(atm, 16, 0.1; rng=MersenneTwister(3), profile=ScientificProfile())
    fast = ft_sh_phase_screen(atm, 16, 0.1; rng=MersenneTwister(3), profile=FastProfile())
    @test !all(scientific .== fast)

    dm = DeformableMirror(tel; n_act=4, influence_width=0.3)
    coeffs = zeros(4)
    ncpa_fast = NCPA(tel, dm, atm; profile=FastProfile(), coefficients=coeffs)
    ncpa_scientific = NCPA(tel, dm, atm; profile=ScientificProfile(), coefficients=coeffs)
    @test size(ncpa_fast.opd) == size(tel.state.opd)
    @test size(ncpa_scientific.opd) == size(tel.state.opd)
end

@testset "Multi-layer atmosphere" begin
    tel = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    delta = tel.params.diameter / tel.params.resolution
    one_pixel_speed = delta / tel.params.sampling_time
    quarter_pixel_speed = 0.25 * one_pixel_speed

    atm = MultiLayerAtmosphere(tel;
        r0=0.2,
        L0=25.0,
        fractional_cn2=[1.0],
        wind_speed=[one_pixel_speed],
        wind_direction=[0.0],
        altitude=[0.0],
    )
    advance!(atm, tel; rng=MersenneTwister(3))
    first = copy(atm.state.opd)
    advance!(atm, tel; rng=MersenneTwister(3))
    second = copy(atm.state.opd)
    propagate!(atm, tel)
    @test size(atm.state.opd) == (16, 16)
    @test sum(abs.(tel.state.opd)) > 0
    @test second[:, 2:end] ≈ first[:, 1:end-1]

    stationary = MultiLayerAtmosphere(tel;
        r0=0.2,
        L0=25.0,
        fractional_cn2=[1.0],
        wind_speed=[0.0],
        wind_direction=[0.0],
        altitude=[0.0],
    )
    advance!(stationary, tel; rng=MersenneTwister(4))
    static_first = copy(stationary.state.opd)
    advance!(stationary, tel; rng=MersenneTwister(4))
    @test stationary.state.opd ≈ static_first

    subpixel = MultiLayerAtmosphere(tel;
        r0=0.2,
        L0=25.0,
        fractional_cn2=[1.0],
        wind_speed=[quarter_pixel_speed],
        wind_direction=[0.0],
        altitude=[0.0],
    )
    advance!(subpixel, tel; rng=MersenneTwister(5))
    subpixel_first = copy(subpixel.state.opd)
    for _ in 1:4
        advance!(subpixel, tel; rng=MersenneTwister(5))
    end
    @test subpixel.state.opd[:, 2:end] ≈ subpixel_first[:, 1:end-1]

    function ensemble_std(fractions; nsamp::Int=16)
        acc = 0.0
        for s in 1:nsamp
            atm_local = MultiLayerAtmosphere(tel;
                r0=0.2,
                L0=25.0,
                fractional_cn2=fractions,
                wind_speed=fill(0.0, length(fractions)),
                wind_direction=fill(0.0, length(fractions)),
                altitude=fill(0.0, length(fractions)),
            )
            advance!(atm_local, tel; rng=MersenneTwister(s))
            acc += std(vec(atm_local.state.opd))
        end
        return acc / nsamp
    end

    single_std = ensemble_std([1.0])
    equal_std = ensemble_std([0.5, 0.5])
    uneven_std = ensemble_std([0.8, 0.2])
    @test isapprox(equal_std, single_std; rtol=0.15)
    @test isapprox(uneven_std, single_std; rtol=0.15)

    replay_a = moving_atmosphere_trace(
        seed=11,
        steps=5,
        resolution=16,
        fractional_cn2=[0.7, 0.3],
        wind_speed=[one_pixel_speed, quarter_pixel_speed],
        wind_direction=[0.0, 90.0],
        altitude=[0.0, 5000.0],
    )
    replay_b = moving_atmosphere_trace(
        seed=11,
        steps=5,
        resolution=16,
        fractional_cn2=[0.7, 0.3],
        wind_speed=[one_pixel_speed, quarter_pixel_speed],
        wind_direction=[0.0, 90.0],
        altitude=[0.0, 5000.0],
    )
    replay_c = moving_atmosphere_trace(
        seed=12,
        steps=5,
        resolution=16,
        fractional_cn2=[0.7, 0.3],
        wind_speed=[one_pixel_speed, quarter_pixel_speed],
        wind_direction=[0.0, 90.0],
        altitude=[0.0, 5000.0],
    )
    @test replay_a == replay_b
    @test any(!=(0.0), abs.(replay_a[1] .- replay_c[1]))
end

@testset "Moving-atmosphere regressions" begin
    slopes_a = moving_wfs_slope_trace(seed=21, steps=4)
    slopes_b = moving_wfs_slope_trace(seed=21, steps=4)
    slopes_c = moving_wfs_slope_trace(seed=22, steps=4)
    @test slopes_a == slopes_b
    @test any(norm(slopes_a[i + 1] - slopes_a[i]) > 0 for i in 1:length(slopes_a)-1)
    @test any(norm(slopes_a[i] - slopes_c[i]) > 0 for i in eachindex(slopes_a))

    loop_a = moving_closed_loop_trace(seed=31, steps=5)
    loop_b = moving_closed_loop_trace(seed=31, steps=5)
    loop_c = moving_closed_loop_trace(seed=32, steps=5)
    @test loop_a == loop_b
    @test any(diff(loop_a.slope_norms) .!= 0)
    @test all(isfinite, loop_a.command_norms)
    @test maximum(loop_a.command_norms) > 0
    @test all(>(0), loop_a.wfs_energy)
    @test all(>(0), loop_a.science_energy)
    @test any(abs.(loop_a.slope_norms .- loop_c.slope_norms) .> 0)
end

@testset "Deformable mirror and WFS" begin
    tel = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    dm = DeformableMirror(tel; n_act=4, influence_width=0.3)
    dm.state.coefs .= 0.1
    apply!(dm, tel, DMReplace())
    @test sum(abs.(tel.state.opd)) > 0

    for i in 1:tel.params.resolution, j in 1:tel.params.resolution
        tel.state.opd[i, j] = i
    end
    wfs = ShackHartmann(tel; n_subap=4)
    slopes = measure!(wfs, tel)
    @test length(slopes) == 2 * 4 * 4
    @test maximum(slopes) > 0
end

@testset "Calibration and control" begin
    tel = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    dm = DeformableMirror(tel; n_act=2, influence_width=0.4)
    wfs = ShackHartmann(tel; n_subap=2)
    imat = interaction_matrix(dm, wfs, tel; amplitude=0.1)
    @test size(imat.matrix) == (length(wfs.state.slopes), length(dm.state.coefs))

    recon = ModalReconstructor(imat; gain=1.0)
    cmd = reconstruct(recon, wfs.state.slopes)
    @test length(cmd) == length(dm.state.coefs)

    mapped = MappedReconstructor(Matrix{Float64}(I, length(dm.state.coefs), length(dm.state.coefs)), imat; gain=0.5)
    mapped_cmd = reconstruct(mapped, wfs.state.slopes)
    @test length(mapped_cmd) == length(dm.state.coefs)
    @test mapped.n_control_modes == length(dm.state.coefs)

    delay = VectorDelayLine(dm.state.coefs, 1)
    delayed0 = shift_delay!(delay, fill(1.0, length(dm.state.coefs)))
    @test all(iszero, delayed0)
    delayed1 = shift_delay!(delay, fill(2.0, length(dm.state.coefs)))
    @test all(==(1.0), delayed1)

    ctrl = DiscreteIntegratorController(length(wfs.state.slopes); gain=0.1, tau=0.02)
    dm_cmd = update!(ctrl, wfs.state.slopes, 0.01)
    @test length(dm_cmd) == length(wfs.state.slopes)

    vault = CalibrationVault(imat.matrix; build_backend=CPUBuildBackend())
    @test vault.M isa Matrix
    recon_cpu = ModalReconstructor(imat; build_backend=CPUBuildBackend())
    @test recon_cpu.reconstructor isa Matrix
end

@testset "AO188/3k simulation" begin
    default_params = AO188SimulationParams()
    @test default_params.n_act == 64
    @test default_params.n_active_actuators == 3228
    @test default_params.n_control_modes == 188
    @test default_params.n_low_order_subap == 2
    @test default_params.low_order_resolution == 28
    @test default_params.n_low_order_modes == 4
    @test default_params.latency.high_measurement_delay_frames == 1
    @test default_params.high_detector.noise isa NoisePhotonReadout
    @test default_params.branch_execution isa SequentialExecution
    @test default_params.replay_mode isa DirectReplayMode

    params = AO188SimulationParams(
        T=Float32,
        resolution=48,
        n_act=16,
        n_active_actuators=180,
        n_control_modes=24,
        control_grid_side=6,
        n_subap=4,
        n_low_order_subap=2,
        n_low_order_modes=3,
        source_magnitude=0.0,
    )
    @test params.low_order_resolution == 12
    surrogate = subaru_ao188_simulation(; params=params, rng=MersenneTwister(1))
    @test count(surrogate.active_mask) == params.n_active_actuators
    @test length(surrogate.active_indices) == params.n_active_actuators
    @test !isnothing(surrogate.dm.state.separable_x)
    @test !isnothing(surrogate.dm.state.separable_y_t)
    @test size(surrogate.high_M2C) == (params.n_act^2, params.n_control_modes)
    @test size(surrogate.low_M2C) == (params.n_act^2, params.n_low_order_modes)
    @test size(surrogate.high_reconstructor.command_basis, 1) == params.n_act^2
    @test size(surrogate.high_reconstructor.command_basis, 2) == params.n_control_modes
    @test size(surrogate.high_reconstructor.reconstructor, 1) == params.n_control_modes
    @test size(surrogate.high_reconstructor.reconstructor, 2) == length(surrogate.high_wfs.state.slopes)
    @test size(surrogate.low_reconstructor.command_basis, 2) == params.n_low_order_modes
    @test size(surrogate.low_reconstructor.reconstructor, 2) == length(surrogate.low_wfs.state.slopes)
    @test surrogate.high_reconstructor.n_control_modes == params.n_control_modes
    @test surrogate.low_reconstructor.n_control_modes == params.n_low_order_modes

    shifted = DeformableMirror(Telescope(resolution=32, diameter=8.0, sampling_time=1e-3, T=Float32);
        n_act=8, T=Float32, misregistration=Misregistration(shift_x=0.01, shift_y=-0.02, T=Float32))
    @test !isnothing(shifted.state.separable_x)

    rotated = DeformableMirror(Telescope(resolution=32, diameter=8.0, sampling_time=1e-3, T=Float32);
        n_act=8, T=Float32, misregistration=Misregistration(rotation_deg=1.0, T=Float32))
    @test isnothing(rotated.state.separable_x)

    step!(surrogate)
    step!(surrogate)
    step!(surrogate)
    @test length(surrogate.command) == params.n_act^2
    @test maximum(abs, surrogate.command) > 0
    @test supports_prepared_runtime(typeof(surrogate))
    @test supports_detector_output(typeof(surrogate))
    @test supports_grouped_execution(typeof(surrogate))
    simif = simulation_interface(surrogate)
    @test simulation_command(simif) === surrogate.command
    @test length(simulation_slopes(simif)) == 2
    timing = runtime_timing(surrogate; warmup=1, samples=2, gc_before=false)
    @test timing.samples == 2
    phase = subaru_ao188_phase_timing(surrogate; warmup=1, samples=2, gc_before=false)
    @test phase.samples == 2
    @test phase.delay_mean_ns >= 0

    experimental = AO188SimulationParams(
        T=Float32,
        resolution=32,
        n_act=12,
        n_active_actuators=96,
        n_control_modes=12,
        control_grid_side=4,
        n_subap=4,
        n_low_order_subap=2,
        n_low_order_modes=2,
        source_magnitude=0.0,
        branch_execution=ThreadedExecution(),
        replay_mode=PreparedReplayMode(),
    )
    @test experimental.low_order_resolution == 8
    surrogate_exp = subaru_ao188_simulation(; params=experimental, rng=MersenneTwister(2))
    @test surrogate_exp.replay_prepared
    step!(surrogate_exp)
    @test maximum(abs, surrogate_exp.command) >= 0

    stream_mode = AO188SimulationParams(
        T=Float32,
        resolution=32,
        n_act=12,
        n_active_actuators=96,
        n_control_modes=12,
        control_grid_side=4,
        n_subap=4,
        n_low_order_subap=2,
        n_low_order_modes=2,
        source_magnitude=0.0,
        branch_execution=BackendStreamExecution(),
    )
    surrogate_stream = subaru_ao188_simulation(; params=stream_mode, rng=MersenneTwister(3))
    step!(surrogate_stream)
    @test maximum(abs, surrogate_stream.command) >= 0

    curvature_params = AO188CurvatureSimulationParams(
        T=Float32,
        resolution=48,
        n_act=16,
        n_active_actuators=180,
        n_control_modes=24,
        control_grid_side=6,
        n_subap=4,
        n_low_order_subap=2,
        n_low_order_modes=3,
        source_magnitude=0.0,
    )
    curvature_sim = subaru_ao188_curvature_simulation(; params=curvature_params, rng=MersenneTwister(4))
    @test curvature_sim.high_wfs isa CurvatureWFS
    @test curvature_sim.high_wfs.params.readout_model isa CurvatureCountingReadout
    @test curvature_sim.high_wfs.params.readout_crop_resolution == 32
    @test curvature_sim.high_wfs.params.readout_pixels_per_subap == 1
    @test curvature_params.high_detector isa AO188APDDetectorConfig
    @test curvature_sim.high_detector isa APDDetector
    step!(curvature_sim)
    @test length(curvature_sim.command) == curvature_params.n_act^2
    curvature_readout = simulation_interface(curvature_sim)
    @test size(simulation_wfs_frame(curvature_readout)[1]) == (2, curvature_params.n_subap^2)
    @test simulation_wfs_metadata(curvature_readout)[1] isa CountingDetectorExportMetadata

    ao3k_params = AO3kSimulationParams(
        T=Float32,
        resolution=64,
        n_act=16,
        n_active_actuators=180,
        n_control_modes=32,
        control_grid_side=6,
        n_subap=8,
        n_low_order_subap=2,
        n_low_order_modes=3,
        source_magnitude=0.0,
    )
    ao3k_sim = subaru_ao3k_simulation(; params=ao3k_params, rng=MersenneTwister(5))
    @test ao3k_sim.high_wfs isa PyramidWFS
    @test ao3k_sim.high_detector isa Detector
    @test ao3k_sim.high_detector.params.sensor isa HgCdTeAvalancheArraySensor
    @test ao3k_params.high_detector.thermal_model isa FixedTemperature
    step!(ao3k_sim)
    @test length(ao3k_sim.command) == ao3k_params.n_act^2
    ao3k_iface = simulation_interface(ao3k_sim)
    ao3k_high_meta = simulation_wfs_metadata(ao3k_iface)[1]
    @test ao3k_high_meta.sensor == :hgcdte_avalanche_array
    @test ao3k_high_meta.thermal_model == :fixed_temperature
    @test ao3k_high_meta.detector_temperature_K == 80.0f0
    @test ao3k_high_meta.provides_combined_frame
    @test ao3k_high_meta.provides_reference_frame
    @test ao3k_high_meta.provides_signal_frame
end

function closed_loop_runtime_allocations()
    rng = MersenneTwister(0)
    tel = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    src = Source(band=:I, magnitude=0.0)
    atm = KolmogorovAtmosphere(tel; r0=0.2, L0=25.0)
    dm = DeformableMirror(tel; n_act=4, influence_width=0.3)
    wfs = ShackHartmann(tel; n_subap=4)
    sim = AOSimulation(tel, atm, src, dm, wfs)
    imat = interaction_matrix(dm, wfs, tel; amplitude=0.1)
    recon = ModalReconstructor(imat; gain=0.5)
    runtime = ClosedLoopRuntime(sim, recon; rng=rng)
    step!(runtime)
    step!(runtime)
    return @allocated step!(runtime)
end

@testset "Closed-loop runtime" begin
    rng = MersenneTwister(0)
    tel = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    src = Source(band=:I, magnitude=0.0)
    atm = KolmogorovAtmosphere(tel; r0=0.2, L0=25.0)
    dm = DeformableMirror(tel; n_act=4, influence_width=0.3)
    wfs = ShackHartmann(tel; n_subap=4)
    sim = AOSimulation(tel, atm, src, dm, wfs)
    imat = interaction_matrix(dm, wfs, tel; amplitude=0.1)
    recon = ModalReconstructor(imat; gain=0.5)
    det = Detector(noise=NoiseNone(), integration_time=1.0, qe=1.0, binning=1)
    runtime = ClosedLoopRuntime(sim, recon; rng=rng, science_detector=det)
    @test runtime isa AbstractControlSimulation
    @test !supports_prepared_runtime(runtime)
    @test supports_detector_output(runtime)
    @test runtime_profile(runtime) isa ScientificRuntimeProfile
    @test runtime_latency(runtime).measurement_delay_frames == 0
    @test runtime.science_zero_padding == 2
    @test !runtime.prepared

    step!(runtime)
    @test length(runtime.command) == length(dm.state.coefs)
    @test size(output_frame(det)) == (32, 32)
    @test closed_loop_runtime_allocations() == 0
    @test simulation_command(runtime) === runtime.command
    @test simulation_science_frame(runtime) === output_frame(det)

    boundary = SimulationInterface(runtime)
    readout = simulation_readout(boundary)
    @test length(simulation_slopes(boundary)) == length(wfs.state.slopes)
    @test simulation_slopes(readout) === simulation_slopes(boundary)
    @test size(simulation_science_frame(boundary)) == size(output_frame(det))
    science_metadata = simulation_science_metadata(boundary)
    @test science_metadata isa DetectorExportMetadata
    @test science_metadata.output_size == size(output_frame(det))
    @test science_metadata.frame_size == size(det.state.frame)
    step!(boundary)
    @test simulation_command(boundary) == runtime.command

    rng2 = MersenneTwister(2)
    tel2 = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    src2 = Source(band=:I, magnitude=0.0)
    atm2 = KolmogorovAtmosphere(tel2; r0=0.2, L0=25.0)
    dm2 = DeformableMirror(tel2; n_act=4, influence_width=0.3)
    wfs2 = ShackHartmann(tel2; n_subap=4, mode=Diffractive())
    sim2 = AOSimulation(tel2, atm2, src2, dm2, wfs2)
    imat2 = interaction_matrix(dm2, wfs2, tel2, src2; amplitude=0.1)
    recon2 = ModalReconstructor(imat2; gain=0.5)
    wfs_det = Detector(noise=NoiseNone(), integration_time=1.0, qe=1.0, binning=1)
    runtime2 = ClosedLoopRuntime(sim2, recon2; rng=rng2, wfs_detector=wfs_det)
    @test supports_prepared_runtime(runtime2)
    prepare!(runtime2)
    @test runtime2.wfs.state.calibrated
    @test supports_detector_output(runtime2)
    step!(runtime2)
    @test simulation_wfs_frame(runtime2) === wfs2.state.spot_cube
    boundary2 = SimulationInterface(runtime2)
    @test ndims(simulation_wfs_frame(boundary2)) == 3
    @test size(simulation_wfs_frame(boundary2), 1) == wfs2.params.n_subap^2

    composite = CompositeSimulationInterface(boundary, boundary2)
    @test supports_grouped_execution(composite)
    @test length(simulation_command(composite)) == length(simulation_command(boundary)) + length(simulation_command(boundary2))
    @test length(simulation_slopes(composite)) == length(simulation_slopes(boundary)) + length(simulation_slopes(boundary2))
    @test length(simulation_wfs_frame(composite)) == 2
    step!(composite)
    @test simulation_command(composite) == vcat(simulation_command(boundary), simulation_command(boundary2))

    timing = runtime_timing(runtime; warmup=1, samples=5, gc_before=false)
    @test timing.samples == 5
    @test timing.min_ns >= 0
    @test timing.max_ns >= timing.min_ns
    phase_runtime = runtime_phase_timing(runtime; warmup=1, samples=3, gc_before=false)
    @test phase_runtime.delay_mean_ns >= 0

    rng3 = MersenneTwister(3)
    tel3 = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    src3 = Source(band=:I, magnitude=0.0)
    atm3 = KolmogorovAtmosphere(tel3; r0=0.2, L0=25.0)
    dm3 = DeformableMirror(tel3; n_act=4, influence_width=0.3)
    wfs3 = ZernikeWFS(tel3; n_subap=4, diffraction_padding=2)
    sim3 = AOSimulation(tel3, atm3, src3, dm3, wfs3)
    imat3 = interaction_matrix(dm3, wfs3, tel3, src3; amplitude=1e-8)
    recon3 = ModalReconstructor(imat3; gain=0.5)
    wfs_det3 = Detector(noise=NoiseNone(), integration_time=1.0, qe=1.0, binning=1)
    runtime3 = ClosedLoopRuntime(sim3, recon3; rng=rng3, wfs_detector=wfs_det3)
    @test supports_prepared_runtime(runtime3)
    @test supports_detector_output(runtime3)
    prepare!(runtime3)
    @test runtime3.prepared
    @test runtime3.wfs.state.calibrated
    step!(runtime3)
    @test length(runtime3.command) == length(dm3.state.coefs)
    @test simulation_wfs_frame(runtime3) === wfs3.state.camera_frame
    @test simulation_wfs_frame(runtime3) == output_frame(wfs_det3)
    @test size(simulation_wfs_frame(runtime3)) == size(wfs3.state.camera_frame)
    @test all(isfinite, simulation_slopes(runtime3))

    rng4 = MersenneTwister(4)
    tel4 = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    src4 = Source(band=:I, magnitude=0.0)
    atm4 = KolmogorovAtmosphere(tel4; r0=0.2, L0=25.0)
    dm4 = DeformableMirror(tel4; n_act=4, influence_width=0.3)
    wfs4 = ShackHartmann(tel4; n_subap=4)
    sim4 = AOSimulation(tel4, atm4, src4, dm4, wfs4)
    imat4 = interaction_matrix(dm4, wfs4, tel4; amplitude=0.1)
    recon4 = ModalReconstructor(imat4; gain=0.5)
    runtime4 = ClosedLoopRuntime(sim4, recon4;
        rng=rng4,
        profile=HILRuntimeProfile(),
        latency=RuntimeLatencyModel(
            measurement_delay_frames=1,
            readout_delay_frames=1,
            reconstruction_delay_frames=1,
            dm_delay_frames=1,
        ),
    )
    @test runtime_profile(runtime4) isa HILRuntimeProfile
    @test runtime4.science_zero_padding == 0
    slope_norms = Float64[]
    command_norms = Float64[]
    dm_norms = Float64[]
    for _ in 1:5
        step!(runtime4)
        push!(slope_norms, norm(simulation_slopes(runtime4)))
        push!(command_norms, norm(simulation_command(runtime4)))
        push!(dm_norms, norm(runtime4.dm.state.coefs))
    end
    @test slope_norms[1] == 0
    @test slope_norms[2] == 0
    @test slope_norms[3] > 0
    @test command_norms[1] == 0
    @test command_norms[2] == 0
    @test command_norms[3] == 0
    @test command_norms[4] > 0
    @test dm_norms[1] == 0
    @test dm_norms[2] == 0
    @test dm_norms[3] == 0
    @test dm_norms[4] == 0
    @test dm_norms[5] > 0
end

@testset "Detector" begin
    psf = fill(1.0, 8, 8)
    det = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=2)
    frame = capture!(det, psf; rng=MersenneTwister(2))
    @test size(frame) == (4, 4)
    @test sum(frame) == sum(psf)

    det_sampling = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, psf_sampling=2, binning=2)
    frame_sampling = capture!(det_sampling, psf; rng=MersenneTwister(2))
    @test size(frame_sampling) == (2, 2)
    @test sum(frame_sampling) == sum(psf)

    det_tuple = Detector(integration_time=1.0, noise=(NoisePhoton(), NoiseReadout(0.5)),
        qe=1.0, binning=1)
    @test det_tuple.noise isa NoisePhotonReadout

    det_sat = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1, full_well=5.0)
    frame_sat = capture!(det_sat, fill(10.0, 4, 4); rng=MersenneTwister(2))
    @test maximum(frame_sat) == 5.0

    det_adc = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1, bits=8, full_well=10.0)
    frame_adc = capture!(det_adc, fill(10.0, 4, 4); rng=MersenneTwister(2))
    @test frame_adc isa Matrix{UInt8}
    @test output_frame(det_adc) === frame_adc
    @test maximum(frame_adc) == 0xff
    @test minimum(frame_adc) >= 0x00
    @test eltype(det_adc.state.frame) == Float64
    metadata_adc = detector_export_metadata(det_adc)
    @test metadata_adc.noise == :none
    @test metadata_adc.sensor == :ccd
    @test metadata_adc.output_precision == UInt8
    @test metadata_adc.frame_size == (4, 4)
    @test metadata_adc.output_size == (4, 4)

    det_adc_float = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        bits=8, full_well=10.0, output_precision=Float32)
    frame_adc_float = capture!(det_adc_float, fill(10.0, 4, 4); rng=MersenneTwister(2))
    @test frame_adc_float isa Matrix{Float32}
    @test maximum(frame_adc_float) == Float32(255.0)

    det_window = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        readout_window=FrameWindow(2:3, 2:4))
    psf_window = reshape(collect(1.0:16.0), 4, 4)
    frame_window = copy(capture!(det_window, psf_window; rng=MersenneTwister(2)))
    @test size(frame_window) == (2, 3)
    @test frame_window == psf_window[2:3, 2:4]
    meta_window = detector_export_metadata(det_window)
    @test meta_window.window_rows == (2, 3)
    @test meta_window.window_cols == (2, 4)
    det_window_oob = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        readout_window=FrameWindow(2:5, 1:2))
    @test_throws DimensionMismatchError capture!(det_window_oob, psf_window; rng=MersenneTwister(2))
    @test_throws InvalidConfiguration FrameWindow(0:1, 1:2)

    det_dark = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1, dark_current=100.0)
    frame_dark = capture!(det_dark, zeros(4, 4); rng=MersenneTwister(2))
    @test sum(frame_dark) > 0

    zero_psf = zeros(4, 4)
    rng_ccd = MersenneTwister(7)
    rng_emccd = MersenneTwister(7)
    det_ccd = Detector(integration_time=1.0, noise=NoiseReadout(1.0), qe=1.0, binning=1,
        gain=10.0, sensor=CCDSensor())
    det_emccd = Detector(integration_time=1.0, noise=NoiseReadout(1.0), qe=1.0, binning=1,
        gain=10.0, sensor=EMCCDSensor())
    frame_ccd = copy(capture!(det_ccd, zero_psf; rng=rng_ccd))
    frame_emccd = copy(capture!(det_emccd, zero_psf; rng=rng_emccd))
    @test frame_ccd ≈ 10 .* frame_emccd
    @test_throws InvalidConfiguration EMCCDSensor(excess_noise_factor=0.5)

    uniform_signal = fill(50.0, 8, 8)
    det_emccd_base = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        gain=1.0, sensor=EMCCDSensor())
    det_emccd_excess = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        gain=1.0, sensor=EMCCDSensor(excess_noise_factor=sqrt(2.0)))
    frame_emccd_base = copy(capture!(det_emccd_base, uniform_signal; rng=MersenneTwister(8)))
    frame_emccd_excess = copy(capture!(det_emccd_excess, uniform_signal; rng=MersenneTwister(8)))
    @test frame_emccd_base == uniform_signal
    @test std(vec(frame_emccd_excess)) > 0
    det_emccd_stochastic = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        gain=5.0, sensor=EMCCDSensor(excess_noise_factor=1.4,
            multiplication_model=StochasticMultiplicationRegister(0.6)))
    frame_emccd_stochastic = copy(capture!(det_emccd_stochastic, uniform_signal; rng=MersenneTwister(124)))
    @test std(vec(frame_emccd_stochastic)) > 0

    det_ccd_cic = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        sensor=CCDSensor(clock_induced_charge_rate=5.0))
    frame_ccd_cic = capture!(det_ccd_cic, zero_psf; rng=MersenneTwister(11))
    @test sum(frame_ccd_cic) > 0
    @test supports_clock_induced_charge(det_ccd_cic.params.sensor)
    @test_throws InvalidConfiguration CCDSensor(clock_induced_charge_rate=-1.0)
    det_emccd_cic = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        sensor=EMCCDSensor(cic_rate=3.0))
    @test sum(capture!(det_emccd_cic, zero_psf; rng=MersenneTwister(125))) > 0
    det_emccd_sat = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        gain=5.0, sensor=EMCCDSensor(register_full_well=100.0))
    @test maximum(capture!(det_emccd_sat, fill(50.0, 4, 4); rng=MersenneTwister(126))) == 100.0
    @test_throws InvalidConfiguration EMCCDSensor(cic_rate=-1.0)
    @test_throws InvalidConfiguration EMCCDSensor(register_full_well=0.0)
    @test_throws InvalidConfiguration EMCCDSensor(multiplication_model=StochasticMultiplicationRegister(-1.0))

    det_cmos = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        sensor=CMOSSensor(column_readout_sigma=1.0))
    frame_cmos = copy(capture!(det_cmos, zeros(8, 8); rng=MersenneTwister(12)))
    @test !all(iszero, frame_cmos)
    @test all(j -> isapprox(std(frame_cmos[:, j]), 0.0; atol=1e-8), axes(frame_cmos, 2))
    @test std(vec(frame_cmos[1, :])) > 0
    @test supports_column_readout_noise(det_cmos.params.sensor)
    @test detector_export_metadata(det_cmos).frame_response == :gaussian
    @test_throws InvalidConfiguration CMOSSensor(column_readout_sigma=-1.0)
    prnu_map = [1.0 0.5 1.0 0.5; 1.0 0.5 1.0 0.5; 1.0 0.5 1.0 0.5; 1.0 0.5 1.0 0.5]
    dsnu_map = fill(0.25, 4, 4)
    bad_mask = falses(4, 4)
    bad_mask[2, 3] = true
    det_cmos_structured = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        sensor=CMOSSensor(output_model=StaticCMOSOutputPattern(2, [1.0, 2.0], [0.0, 10.0]),
            timing_model=RollingShutter(1e-3)),
        response_model=NullFrameResponse(),
        defect_model=CompositeDetectorDefectModel(
            PixelResponseNonuniformity(prnu_map),
            DarkSignalNonuniformity(dsnu_map),
            BadPixelMask(bad_mask; throughput=0.0)))
    structured_frame = capture!(det_cmos_structured, fill(2.0, 4, 4); rng=MersenneTwister(120))
    @test structured_frame[1, 1] ≈ 2.25
    @test structured_frame[1, 2] ≈ 1.25
    @test structured_frame[1, 3] ≈ 14.5
    @test structured_frame[2, 3] ≈ 10.5
    structured_meta = detector_export_metadata(det_cmos_structured)
    @test structured_meta.detector_defects == :composite
    @test structured_meta.has_prnu
    @test structured_meta.has_dsnu
    @test structured_meta.has_bad_pixels
    @test structured_meta.timing_model == :rolling_shutter
    @test structured_meta.timing_line_time == 1e-3
    @test structured_meta.sampling_wallclock_time == 1.004
    @test supports_detector_defect_maps(det_cmos_structured.params.sensor)
    @test supports_shutter_timing(det_cmos_structured.params.sensor)
    @test_throws InvalidConfiguration CMOSSensor(timing_model=RollingShutter(-1.0))
    @test_throws InvalidConfiguration CMOSSensor(output_model=StaticCMOSOutputPattern(2, [1.0], [0.0, 1.0]))

    det_ingaas = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        sensor=InGaAsSensor(glow_rate=3.0))
    frame_ingaas = capture!(det_ingaas, zero_psf; rng=MersenneTwister(13))
    @test sum(frame_ingaas) > 0
    @test supports_sensor_glow(det_ingaas.params.sensor)
    @test detector_export_metadata(det_ingaas).frame_response == :gaussian
    @test_throws InvalidConfiguration InGaAsSensor(glow_rate=-1.0)
    det_ingaas_persist = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        response_model=NullFrameResponse(),
        sensor=InGaAsSensor(persistence_model=ExponentialPersistence(0.5, 0.0)))
    capture!(det_ingaas_persist, fill(4.0, 4, 4); rng=MersenneTwister(121))
    persisted = capture!(det_ingaas_persist, zeros(4, 4); rng=MersenneTwister(122))
    @test sum(persisted) ≈ 32.0
    persist_meta = detector_export_metadata(det_ingaas_persist)
    @test persist_meta.persistence_model == :exponential
    det_ingaas_nonlinear = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        response_model=NullFrameResponse(),
        nonlinearity_model=SaturatingFrameNonlinearity(0.1),
        sensor=InGaAsSensor())
    nonlinear_frame = capture!(det_ingaas_nonlinear, fill(10.0, 2, 2); rng=MersenneTwister(123))
    @test nonlinear_frame == fill(5.0, 2, 2)
    nonlinear_meta = detector_export_metadata(det_ingaas_nonlinear)
    @test nonlinear_meta.nonlinearity_model == :saturating
    @test supports_detector_persistence(det_ingaas_persist.params.sensor)
    @test supports_detector_nonlinearity(det_ingaas_nonlinear.params.sensor)
    @test_throws InvalidConfiguration InGaAsSensor(persistence_model=ExponentialPersistence(1.1, 0.0))
    @test_throws InvalidConfiguration Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        nonlinearity_model=SaturatingFrameNonlinearity(-0.1), sensor=InGaAsSensor())

    arrhenius = ArrheniusRateLaw(300.0, 6000.0)
    linear = LinearTemperatureLaw(300.0, 0.01)
    exp_law = ExponentialTemperatureLaw(300.0, 0.01)
    @test evaluate_temperature_law(arrhenius, 10.0, 80.0) < 10.0
    @test evaluate_temperature_law(linear, 2.0, 250.0) ≈ 1.0
    @test evaluate_temperature_law(exp_law, 2.0, 250.0) < 2.0

    thermal_det = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        dark_current=10.0,
        response_model=NullFrameResponse(),
        thermal_model=FixedTemperature(temperature_K=80.0, dark_current_law=arrhenius),
        sensor=CCDSensor())
    thermal_meta = detector_export_metadata(thermal_det)
    @test supports_detector_thermal_model(thermal_det)
    @test !supports_dynamic_thermal_state(thermal_det.params.thermal_model)
    @test supports_temperature_dependent_dark_current(thermal_det)
    @test detector_temperature(thermal_det) == 80.0
    @test thermal_meta.thermal_model == :fixed_temperature
    @test thermal_meta.detector_temperature_K == 80.0
    @test thermal_meta.cooling_setpoint_K == 80.0
    @test thermal_meta.dark_current_law == :arrhenius
    @test effective_dark_current(thermal_det) < thermal_det.params.dark_current
    @test thermal_state(thermal_det) isa NoThermalState
    @test advance_thermal!(thermal_det, 1.0) === thermal_det

    dynamic_model = FirstOrderThermalModel(
        ambient_temperature_K=295.0,
        setpoint_temperature_K=120.0,
        initial_temperature_K=300.0,
        time_constant_s=2.0,
        min_temperature_K=80.0,
        max_temperature_K=320.0,
        dark_current_law=arrhenius,
    )
    dynamic_det = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        dark_current=10.0,
        response_model=NullFrameResponse(),
        thermal_model=dynamic_model,
        sensor=CCDSensor())
    dynamic_meta = detector_export_metadata(dynamic_det)
    @test supports_detector_thermal_model(dynamic_det)
    @test supports_dynamic_thermal_state(dynamic_det.params.thermal_model)
    @test thermal_state(dynamic_det) isa DetectorThermalState
    @test detector_temperature(dynamic_det) == 300.0
    @test dynamic_meta.thermal_model == :first_order
    @test dynamic_meta.detector_temperature_K == 300.0
    @test dynamic_meta.ambient_temperature_K == 295.0
    @test dynamic_meta.cooling_setpoint_K == 120.0
    @test dynamic_meta.thermal_time_constant_s == 2.0
    dark_current_initial = effective_dark_current(dynamic_det)
    @test advance_thermal!(dynamic_det, 2.0) === dynamic_det
    @test detector_temperature(dynamic_det) ≈ 120.0 + 180.0 * exp(-1.0)
    @test effective_dark_current(dynamic_det) < dark_current_initial
    reset_integration!(dynamic_det)
    capture!(dynamic_det, fill(1.0f0, 4, 4); rng=MersenneTwister(24))
    @test detector_temperature(dynamic_det) < 120.0 + 180.0 * exp(-1.0)
    @test_throws InvalidConfiguration advance_thermal!(dynamic_det, -1.0)
    @test_throws InvalidConfiguration FirstOrderThermalModel(
        ambient_temperature_K=295.0,
        setpoint_temperature_K=120.0,
        initial_temperature_K=60.0,
        time_constant_s=2.0,
        min_temperature_K=80.0,
        max_temperature_K=320.0,
    )

    thermal_ingaas = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        response_model=NullFrameResponse(),
        thermal_model=FixedTemperature(temperature_K=250.0, glow_rate_law=linear),
        sensor=InGaAsSensor(glow_rate=2.0))
    @test supports_temperature_dependent_glow(thermal_ingaas)
    @test effective_glow_rate(thermal_ingaas) ≈ 1.0

    thermal_emccd = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        response_model=NullFrameResponse(),
        thermal_model=FixedTemperature(temperature_K=250.0, cic_rate_law=linear),
        sensor=EMCCDSensor(cic_rate=2.0))
    @test effective_cic_rate(thermal_emccd) ≈ 1.0

    det_saphira = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        gain=1.0, sensor=HgCdTeAvalancheArraySensor(avalanche_gain=5.0))
    frame_saphira = copy(capture!(det_saphira, uniform_signal; rng=MersenneTwister(14)))
    @test frame_saphira == 5.0 .* uniform_signal
    @test supports_avalanche_gain(det_saphira.params.sensor)
    @test supports_sensor_glow(det_saphira.params.sensor)
    @test detector_export_metadata(det_saphira).frame_response == :sampled
    det_saphira_excess = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        gain=1.0, sensor=HgCdTeAvalancheArraySensor(avalanche_gain=1.0, excess_noise_factor=sqrt(2.0)))
    frame_saphira_excess = copy(capture!(det_saphira_excess, uniform_signal; rng=MersenneTwister(14)))
    @test std(vec(frame_saphira_excess)) > 0
    det_saphira_sat = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        gain=1.0, full_well=100.0, sensor=HgCdTeAvalancheArraySensor(avalanche_gain=5.0))
    frame_saphira_sat = copy(capture!(det_saphira_sat, uniform_signal; rng=MersenneTwister(15)))
    @test maximum(frame_saphira_sat) == 100.0
    det_saphira_single = Detector(integration_time=1.0, noise=NoiseReadout(4.0), qe=1.0, binning=1,
        gain=1.0, sensor=HgCdTeAvalancheArraySensor())
    frame_saphira_single = copy(capture!(det_saphira_single, zero_psf; rng=MersenneTwister(16)))
    single_products = readout_products(det_saphira_single)
    @test single_products isa HgCdTeReadoutProducts
    @test detector_reference_frame(det_saphira_single) === nothing
    @test detector_signal_frame(det_saphira_single) !== nothing
    @test detector_combined_frame(det_saphira_single) == frame_saphira_single
    @test detector_reference_cube(det_saphira_single) === nothing
    @test detector_signal_cube(det_saphira_single) !== nothing
    @test detector_read_cube(det_saphira_single) === nothing
    @test detector_read_times(det_saphira_single) === nothing
    det_saphira_ndr = Detector(integration_time=1.0, noise=NoiseReadout(4.0), qe=1.0, binning=1,
        gain=1.0, sensor=HgCdTeAvalancheArraySensor(sampling_mode=AveragedNonDestructiveReads(4)))
    frame_saphira_ndr = copy(capture!(det_saphira_ndr, zero_psf; rng=MersenneTwister(16)))
    @test std(vec(frame_saphira_ndr)) < std(vec(frame_saphira_single))
    @test supports_nondestructive_reads(det_saphira_ndr.params.sensor)
    @test supports_readout_correction(det_saphira_ndr.params.sensor)
    @test supports_read_cube(det_saphira_ndr.params.sensor)
    saphira_meta = detector_export_metadata(det_saphira_ndr)
    @test saphira_meta.sampling_mode == :averaged_non_destructive_reads
    @test saphira_meta.sampling_reads == 4
    @test saphira_meta.sampling_reference_reads == 0
    @test saphira_meta.sampling_signal_reads == 4
    @test saphira_meta.readout_sigma == 2.0
    @test saphira_meta.provides_signal_frame
    @test !saphira_meta.provides_reference_frame
    @test saphira_meta.provides_combined_frame
    @test !saphira_meta.provides_reference_cube
    @test saphira_meta.provides_signal_cube
    @test saphira_meta.provides_read_cube
    @test saphira_meta.signal_cube_reads == 4
    @test saphira_meta.read_cube_reads == 4
    @test detector_combined_frame(det_saphira_ndr) == frame_saphira_ndr
    @test size(detector_signal_cube(det_saphira_ndr)) == (4, 4, 4)
    @test length(detector_read_times(det_saphira_ndr)) == 4
    det_saphira_cds = Detector(integration_time=1.0, noise=NoiseReadout(4.0), qe=1.0, binning=1,
        gain=1.0, sensor=HgCdTeAvalancheArraySensor(sampling_mode=CorrelatedDoubleSampling()))
    frame_saphira_cds = copy(capture!(det_saphira_cds, zero_psf; rng=MersenneTwister(16)))
    @test std(vec(frame_saphira_cds)) > std(vec(frame_saphira_single))
    cds_meta = detector_export_metadata(det_saphira_cds)
    @test cds_meta.sampling_mode == :correlated_double_sampling
    @test cds_meta.sampling_reference_reads == 1
    @test cds_meta.sampling_signal_reads == 1
    @test cds_meta.readout_sigma ≈ 4.0 * sqrt(2.0)
    @test cds_meta.provides_reference_cube
    @test cds_meta.provides_signal_cube
    @test cds_meta.reference_cube_reads == 1
    @test cds_meta.signal_cube_reads == 1
    @test size(detector_reference_cube(det_saphira_cds)) == (4, 4, 1)
    @test size(detector_signal_cube(det_saphira_cds)) == (4, 4, 1)
    @test detector_combined_frame(det_saphira_cds) ≈
        detector_signal_frame(det_saphira_cds) .- detector_reference_frame(det_saphira_cds)
    det_saphira_fowler = Detector(integration_time=1.0, noise=NoiseReadout(4.0), qe=1.0, binning=1,
        gain=1.0, sensor=HgCdTeAvalancheArraySensor(sampling_mode=FowlerSampling(8)))
    frame_saphira_fowler = copy(capture!(det_saphira_fowler, zero_psf; rng=MersenneTwister(16)))
    fowler_meta = detector_export_metadata(det_saphira_fowler)
    @test fowler_meta.sampling_mode == :fowler_sampling
    @test fowler_meta.sampling_reads == 16
    @test fowler_meta.sampling_reference_reads == 8
    @test fowler_meta.sampling_signal_reads == 8
    @test fowler_meta.readout_sigma == 2.0
    @test fowler_meta.reference_cube_reads == 8
    @test fowler_meta.signal_cube_reads == 8
    @test detector_combined_frame(det_saphira_fowler) ≈
        detector_signal_frame(det_saphira_fowler) .- detector_reference_frame(det_saphira_fowler)
    det_saphira_timed_single = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        dark_current=1000.0, gain=1.0, sensor=HgCdTeAvalancheArraySensor(read_time=1.0))
    frame_saphira_timed_single = copy(capture!(det_saphira_timed_single, zero_psf; rng=MersenneTwister(17)))
    det_saphira_timed_cds = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        dark_current=1000.0, gain=1.0,
        sensor=HgCdTeAvalancheArraySensor(read_time=1.0, sampling_mode=CorrelatedDoubleSampling()))
    frame_saphira_timed_cds = copy(capture!(det_saphira_timed_cds, zero_psf; rng=MersenneTwister(17)))
    @test sum(frame_saphira_timed_cds) > sum(frame_saphira_timed_single)
    timed_meta = detector_export_metadata(det_saphira_timed_cds)
    @test timed_meta.sampling_read_time == 1.0
    @test timed_meta.sampling_wallclock_time == 3.0
    det_saphira_windowed = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        gain=1.0, sensor=HgCdTeAvalancheArraySensor(read_time=1.0, sampling_mode=CorrelatedDoubleSampling()),
        readout_window=FrameWindow(2:3, 2:3))
    capture!(det_saphira_windowed, fill(10.0, 4, 4); rng=MersenneTwister(18))
    windowed_meta = detector_export_metadata(det_saphira_windowed)
    @test windowed_meta.sampling_read_time == 0.5
    @test windowed_meta.sampling_wallclock_time == 2.0
    @test detector_combined_frame(det_saphira_windowed) !== nothing
    @test size(detector_reference_cube(det_saphira_windowed)) == (2, 2, 1)
    @test size(detector_signal_cube(det_saphira_windowed)) == (2, 2, 1)
    @test size(detector_signal_frame(det_saphira_windowed)) == (2, 2)
    @test size(detector_read_cube(det_saphira_windowed)) == (2, 2, 2)
    @test detector_read_times(det_saphira_windowed) == [0.5, 1.0]
    det_saphira_corrected = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        gain=1.0, sensor=HgCdTeAvalancheArraySensor(),
        response_model=NullFrameResponse(),
        correction_model=ReferencePixelCommonModeCorrection(1, 1))
    corrected_frame = capture!(det_saphira_corrected, fill(5.0, 4, 4); rng=MersenneTwister(19))
    corrected_meta = detector_export_metadata(det_saphira_corrected)
    @test corrected_meta.readout_correction == :reference_pixel_common_mode
    @test corrected_meta.correction_edge_rows == 1
    @test corrected_meta.correction_edge_cols == 1
    @test abs(mean(corrected_frame)) < 1e-6
    det_saphira_row_corrected = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        gain=1.0, sensor=HgCdTeAvalancheArraySensor(),
        response_model=NullFrameResponse(),
        correction_model=ReferenceRowCommonModeCorrection(1))
    row_pattern = repeat(reshape([1.0, 2.0, 3.0, 4.0], :, 1), 1, 4)
    row_corrected = capture!(det_saphira_row_corrected, row_pattern; rng=MersenneTwister(20))
    @test maximum(abs, row_corrected) < 1e-6
    det_saphira_col_corrected = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        gain=1.0, sensor=HgCdTeAvalancheArraySensor(),
        response_model=NullFrameResponse(),
        correction_model=ReferenceColumnCommonModeCorrection(1))
    col_pattern = repeat(reshape([1.0, 2.0, 3.0, 4.0], 1, :), 4, 1)
    col_corrected = capture!(det_saphira_col_corrected, col_pattern; rng=MersenneTwister(21))
    @test maximum(abs, col_corrected) < 1e-6
    det_saphira_output_corrected = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        gain=1.0, sensor=HgCdTeAvalancheArraySensor(),
        response_model=NullFrameResponse(),
        correction_model=ReferenceOutputCommonModeCorrection(2; edge_rows=1, edge_cols=1))
    output_pattern = hcat(fill(5.0, 4, 2), fill(10.0, 4, 2))
    output_corrected = capture!(det_saphira_output_corrected, output_pattern; rng=MersenneTwister(22))
    output_meta = detector_export_metadata(det_saphira_output_corrected)
    @test output_meta.readout_correction == :reference_output_common_mode
    @test output_meta.correction_group_cols == 2
    @test maximum(abs, output_corrected) < 1e-6
    det_saphira_composite = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        gain=1.0, sensor=HgCdTeAvalancheArraySensor(),
        response_model=NullFrameResponse(),
        correction_model=CompositeFrameReadoutCorrection((
            ReferenceRowCommonModeCorrection(1),
            ReferenceColumnCommonModeCorrection(1))))
    composite_pattern = row_pattern .+ col_pattern
    composite_corrected = capture!(det_saphira_composite, composite_pattern; rng=MersenneTwister(23))
    composite_meta = detector_export_metadata(det_saphira_composite)
    @test composite_meta.readout_correction == :composite
    @test composite_meta.correction_stage_count == 2
    @test maximum(abs, composite_corrected) < 1e-6
    @test_throws InvalidConfiguration HgCdTeAvalancheArraySensor(avalanche_gain=0.5)
    @test_throws InvalidConfiguration HgCdTeAvalancheArraySensor(excess_noise_factor=0.5)
    @test_throws InvalidConfiguration HgCdTeAvalancheArraySensor(glow_rate=-1.0)
    @test_throws InvalidConfiguration HgCdTeAvalancheArraySensor(read_time=-1.0)
    @test_throws InvalidConfiguration HgCdTeAvalancheArraySensor(sampling_mode=AveragedNonDestructiveReads(0))
    @test_throws InvalidConfiguration HgCdTeAvalancheArraySensor(sampling_mode=FowlerSampling(0))
    @test_throws InvalidConfiguration ReferencePixelCommonModeCorrection(0, 0)
    @test_throws InvalidConfiguration ReferenceRowCommonModeCorrection(0)
    @test_throws InvalidConfiguration ReferenceColumnCommonModeCorrection(0)
    @test_throws InvalidConfiguration ReferenceOutputCommonModeCorrection(0)
    @test_throws InvalidConfiguration CompositeFrameReadoutCorrection(())

    @test_throws InvalidConfiguration Detector(integration_time=1.0, noise=NoisePhoton(), qe=1.0, binning=1,
        sensor=APDSensor())

    det_default_ccd = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1, sensor=CCDSensor())
    @test detector_export_metadata(det_default_ccd).frame_response == :none

    apd = APDDetector(integration_time=1.0, qe=0.5, gain=2.0, dark_count_rate=0.0, noise=NoiseNone())
    channels = fill(4.0, 2, 8)
    apd_out = capture!(apd, channels; rng=MersenneTwister(9))
    @test apd_out == fill(4.0, 2, 8)
    @test channel_output(apd) === apd_out
    apd_meta = detector_export_metadata(apd)
    @test apd_meta isa CountingDetectorExportMetadata
    @test apd_meta.sensor == :apd
    @test apd_meta.readout.output_size == (2, 8)
    @test apd_meta.readout.n_channels == 16
    @test apd_meta.dead_time_model == :none
    @test apd_meta.dead_time === nothing
    @test supports_counting_noise(apd)
    @test !supports_dead_time(apd)
    @test !supports_channel_gain_map(apd)

    apd_gain_map = APDDetector(integration_time=1.0, qe=1.0, gain=1.0, dark_count_rate=0.0,
        noise=NoiseNone(), channel_gain_map=fill(0.5, 2, 8))
    @test capture!(apd_gain_map, fill(2.0, 2, 8); rng=MersenneTwister(9)) == fill(1.0, 2, 8)
    @test supports_channel_gain_map(apd_gain_map)

    det_mtf = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        response_model=GaussianPixelResponse(response_width_px=0.75))
    impulse = zeros(9, 9)
    impulse[5, 5] = 1.0
    frame_mtf = capture!(det_mtf, impulse; rng=MersenneTwister(3))
    @test sum(frame_mtf) ≈ 1.0 atol=1e-6
    @test frame_mtf[5, 5] < 1.0
    @test frame_mtf[5, 4] > 0
    @test supports_detector_mtf(det_mtf)
    mtf_meta = detector_export_metadata(det_mtf)
    @test mtf_meta.frame_response == :gaussian
    @test mtf_meta.response_width_px == 0.75
    @test mtf_meta.response_application_domain == :image
    @test mtf_meta.response_is_separable
    @test mtf_meta.response_is_shift_invariant
    @test mtf_meta.response_support_rows == mtf_meta.response_support_cols
    @test mtf_meta.pitch_x_px === nothing
    @test mtf_meta.aperture_shape === nothing
    @test_throws InvalidConfiguration GaussianPixelResponse(response_width_px=0.0)

    sampled_kernel = [0.0 0.1 0.0; 0.1 0.6 0.1; 0.0 0.1 0.0]
    sampled_det = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        response_model=SampledFrameResponse(sampled_kernel))
    sampled_frame = capture!(sampled_det, impulse; rng=MersenneTwister(6))
    @test sum(sampled_frame) ≈ 1.0 atol=1e-6
    @test sampled_frame[5, 5] ≈ 0.6 atol=1e-6
    @test sampled_frame[5, 4] ≈ 0.1 atol=1e-6
    sampled_meta = detector_export_metadata(sampled_det)
    @test sampled_meta.frame_response == :sampled
    @test sampled_meta.response_application_domain == :image
    @test !sampled_meta.response_is_separable
    @test sampled_meta.response_support_rows == 3
    @test sampled_meta.response_support_cols == 3
    @test sampled_meta.aperture_shape == :sampled
    @test supports_detector_mtf(sampled_det)
    @test_throws InvalidConfiguration SampledFrameResponse(zeros(3, 3))
    @test_throws InvalidConfiguration SampledFrameResponse(ones(2, 3))

    cube_mtf = Array{Float64}(undef, 2, size(impulse, 1), size(impulse, 2))
    cube_mtf[1, :, :] .= impulse
    cube_mtf[2, :, :] .= impulse
    scratch_mtf = similar(cube_mtf)
    stack_mtf = AdaptiveOpticsSim.capture_stack!(det_mtf, cube_mtf, scratch_mtf; rng=MersenneTwister(10))
    @test size(stack_mtf) == size(cube_mtf)
    @test all(isfinite, stack_mtf)
    cube_sampled = Array{Float64}(undef, 2, size(impulse, 1), size(impulse, 2))
    cube_sampled[1, :, :] .= impulse
    cube_sampled[2, :, :] .= impulse
    scratch_sampled = similar(cube_sampled)
    stack_sampled = AdaptiveOpticsSim.capture_stack!(sampled_det, cube_sampled, scratch_sampled; rng=MersenneTwister(10))
    @test size(stack_sampled) == size(cube_sampled)
    @test all(isfinite, stack_sampled)
    @test stack_sampled[1, :, :] ≈ sampled_frame atol=1e-6
    @test stack_sampled[2, :, :] ≈ sampled_frame atol=1e-6

    rect_det = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        response_model=RectangularPixelAperture(pitch_x_px=2.0, pitch_y_px=2.0,
            fill_factor_x=0.6, fill_factor_y=0.8))
    rect_frame = capture!(rect_det, impulse; rng=MersenneTwister(4))
    @test sum(rect_frame) ≈ 1.0 atol=1e-6
    @test rect_frame[5, 5] < 1.0
    rect_meta = detector_export_metadata(rect_det)
    @test rect_meta.frame_response == :rectangular_aperture
    @test rect_meta.pitch_x_px == 2.0
    @test rect_meta.pitch_y_px == 2.0
    @test rect_meta.fill_factor_x == 0.6
    @test rect_meta.fill_factor_y == 0.8
    @test rect_meta.aperture_shape == :rectangular
    @test rect_meta.response_application_domain == :image
    @test supports_detector_mtf(rect_det)

    mtf_det = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        response_model=SeparablePixelMTF(pitch_x_px=1.0, pitch_y_px=1.0,
            fill_factor_x=0.7, fill_factor_y=0.7))
    mtf_frame = capture!(mtf_det, impulse; rng=MersenneTwister(5))
    @test sum(mtf_frame) ≈ 1.0 atol=1e-6
    mtf_meta2 = detector_export_metadata(mtf_det)
    @test mtf_meta2.frame_response == :separable_mtf
    @test mtf_meta2.response_is_separable
    @test mtf_meta2.response_application_domain == :image
    @test mtf_meta2.aperture_shape == :rectangular
    @test_throws InvalidConfiguration RectangularPixelAperture(fill_factor_x=0.0)
    @test_throws InvalidConfiguration SeparablePixelMTF(fill_factor_y=1.5)

    det_window_stack = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        readout_window=FrameWindow(2:8, 2:8))
    cube_window = Array{Float64}(undef, 2, size(impulse, 1), size(impulse, 2))
    cube_window[1, :, :] .= impulse
    cube_window[2, :, :] .= impulse
    scratch_window = similar(cube_window)
    @test_throws InvalidConfiguration AdaptiveOpticsSim.capture_stack!(det_window_stack, cube_window, scratch_window; rng=MersenneTwister(10))
    input_window_stack = copy(cube_window)
    output_window_stack = Array{Float64}(undef, 2, 7, 7)
    generalized_window = AdaptiveOpticsSim.capture_stack!(det_window_stack, output_window_stack, input_window_stack; rng=MersenneTwister(10))
    @test size(generalized_window) == (2, 7, 7)
    @test generalized_window[1, :, :] ≈ capture!(det_window_stack, impulse; rng=MersenneTwister(10))

    corrected_stack_det = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        sensor=HgCdTeAvalancheArraySensor(sampling_mode=SingleRead()),
        correction_model=ReferencePixelCommonModeCorrection(1, 1))
    corrected_stack_in = fill(10.0, 2, 5, 5)
    corrected_stack = AdaptiveOpticsSim.capture_stack!(corrected_stack_det, corrected_stack_in,
        similar(corrected_stack_in); rng=MersenneTwister(10))
    @test size(corrected_stack) == size(corrected_stack_in)
    @test maximum(abs, corrected_stack) ≤ 1e-6

    det_cmos_batched = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        sensor=CMOSSensor(column_readout_sigma=1.0))
    @test_throws InvalidConfiguration AdaptiveOpticsSim.capture_stack!(det_cmos_batched, cube_mtf, scratch_mtf; rng=MersenneTwister(10))

    det_generalized = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, psf_sampling=2, binning=2,
        bits=8, full_well=10.0, output_precision=UInt8)
    input_generalized = zeros(Float64, 2, 8, 8)
    input_generalized[1, 4, 4] = 10.0
    input_generalized[2, 5, 5] = 10.0
    output_generalized = Array{UInt8}(undef, 2, 2, 2)
    generalized_stack = AdaptiveOpticsSim.capture_stack!(det_generalized, output_generalized, input_generalized; rng=MersenneTwister(10))
    @test size(generalized_stack) == (2, 2, 2)
    @test generalized_stack[1, :, :] == capture!(det_generalized, @view(input_generalized[1, :, :]); rng=MersenneTwister(10))
    @test generalized_stack[2, :, :] == capture!(det_generalized, @view(input_generalized[2, :, :]); rng=MersenneTwister(10))

    apd_dead_time = APDDetector(integration_time=1.0, qe=1.0, gain=1.0, dark_count_rate=0.0,
        noise=NoiseNone(), dead_time_model=NonParalyzableDeadTime(0.5))
    dead_time_out = capture!(apd_dead_time, fill(4.0, 2, 8); rng=MersenneTwister(9))
    @test dead_time_out ≈ fill(4.0 / 3.0, 2, 8)
    dead_time_meta = detector_export_metadata(apd_dead_time)
    @test dead_time_meta.dead_time_model == :nonparalyzable
    @test dead_time_meta.dead_time == 0.5
    @test supports_dead_time(apd_dead_time)
    apd_paralyzable = APDDetector(integration_time=1.0, qe=1.0, gain=1.0, dark_count_rate=0.0,
        noise=NoiseNone(), dead_time_model=ParalyzableDeadTime(0.5))
    @test capture!(apd_paralyzable, fill(4.0, 2, 8); rng=MersenneTwister(9)) ≈ fill(4.0 * exp(-2.0), 2, 8)
    @test supports_paralyzable_dead_time(apd_paralyzable)
    apd_gated = APDDetector(integration_time=1.0, qe=1.0, gain=1.0, dark_count_rate=0.0,
        noise=NoiseNone(), gate_model=DutyCycleGate(0.5))
    @test capture!(apd_gated, fill(4.0, 2, 8); rng=MersenneTwister(9)) == fill(2.0, 2, 8)
    gated_meta = detector_export_metadata(apd_gated)
    @test gated_meta.gate_model == :duty_cycle
    @test gated_meta.duty_cycle == 0.5
    @test supports_counting_gating(apd_gated)
    apd_afterpulse = APDDetector(integration_time=1.0, qe=1.0, gain=1.0, dark_count_rate=0.0,
        noise=NoiseNone(), correlation_model=AfterpulsingModel(0.25))
    @test capture!(apd_afterpulse, fill(4.0, 2, 8); rng=MersenneTwister(9)) == fill(5.0, 2, 8)
    @test supports_afterpulsing(apd_afterpulse)
    apd_crosstalk = APDDetector(integration_time=1.0, qe=1.0, gain=1.0, dark_count_rate=0.0,
        noise=NoiseNone(), correlation_model=ChannelCrosstalkModel(0.4))
    crosstalk_in = zeros(3, 3)
    crosstalk_in[2, 2] = 10.0
    crosstalk_out = capture!(apd_crosstalk, crosstalk_in; rng=MersenneTwister(9))
    @test crosstalk_out[2, 2] ≈ 6.0
    @test crosstalk_out[1, 2] ≈ 1.0
    @test crosstalk_out[2, 1] ≈ 1.0
    @test crosstalk_out[2, 3] ≈ 1.0
    @test crosstalk_out[3, 2] ≈ 1.0
    @test supports_channel_crosstalk(apd_crosstalk)
    apd_composite_corr = APDDetector(integration_time=1.0, qe=1.0, gain=1.0, dark_count_rate=0.0,
        noise=NoiseNone(), correlation_model=CompositeCountingCorrelation(AfterpulsingModel(0.1), ChannelCrosstalkModel(0.2)))
    composite_meta_apd = detector_export_metadata(apd_composite_corr)
    @test composite_meta_apd.correlation_model == :composite
    @test composite_meta_apd.afterpulse_probability == 0.1
    @test composite_meta_apd.crosstalk == 0.2

    @test_throws InvalidConfiguration APDDetector(noise=NoiseReadout(1.0))
    @test_throws InvalidConfiguration APDDetector(dead_time_model=NonParalyzableDeadTime(-1.0))
    @test_throws InvalidConfiguration APDDetector(dead_time_model=ParalyzableDeadTime(-1.0))
    @test_throws InvalidConfiguration APDDetector(gate_model=DutyCycleGate(0.0))
    @test_throws InvalidConfiguration APDDetector(correlation_model=AfterpulsingModel(1.5))
    @test_throws InvalidConfiguration APDDetector(correlation_model=ChannelCrosstalkModel(-0.1))
    apd_thermal = APDDetector(integration_time=1.0, qe=1.0, gain=1.0, dark_count_rate=10.0,
        noise=NoiseNone(), thermal_model=FixedTemperature(temperature_K=80.0, dark_count_law=arrhenius))
    apd_thermal_meta = detector_export_metadata(apd_thermal)
    @test supports_detector_thermal_model(apd_thermal)
    @test supports_temperature_dependent_dark_counts(apd_thermal)
    @test detector_temperature(apd_thermal) == 80.0
    @test apd_thermal_meta.thermal_model == :fixed_temperature
    @test apd_thermal_meta.dark_count_law == :arrhenius
    @test effective_dark_count_rate(apd_thermal) < apd_thermal.params.dark_count_rate
    apd_dynamic = APDDetector(integration_time=1.0, qe=1.0, gain=1.0, dark_count_rate=10.0,
        noise=NoiseNone(), thermal_model=FirstOrderThermalModel(
            ambient_temperature_K=295.0,
            setpoint_temperature_K=120.0,
            initial_temperature_K=300.0,
            time_constant_s=2.0,
            min_temperature_K=80.0,
            max_temperature_K=320.0,
            dark_count_law=arrhenius))
    apd_dynamic_initial = effective_dark_count_rate(apd_dynamic)
    @test supports_detector_thermal_model(apd_dynamic)
    @test supports_dynamic_thermal_state(apd_dynamic.params.thermal_model)
    @test thermal_state(apd_dynamic) isa DetectorThermalState
    @test advance_thermal!(apd_dynamic, 2.0) === apd_dynamic
    @test detector_temperature(apd_dynamic) ≈ 120.0 + 180.0 * exp(-1.0)
    @test effective_dark_count_rate(apd_dynamic) < apd_dynamic_initial
    capture!(apd_dynamic, fill(1.0, 2, 2); rng=MersenneTwister(25))
    @test detector_temperature(apd_dynamic) < 120.0 + 180.0 * exp(-1.0)

    det_buffered = Detector(integration_time=2.0, noise=NoiseNone(), qe=1.0, binning=1)
    frame_partial = copy(capture!(det_buffered, fill(1.0, 4, 4); rng=MersenneTwister(2), sample_time=1.0))
    @test !readout_ready(det_buffered)
    @test sum(frame_partial) == 16.0
    frame_buffered = copy(capture!(det_buffered, fill(1.0, 4, 4); rng=MersenneTwister(2), sample_time=1.0))
    @test readout_ready(det_buffered)
    @test sum(frame_buffered) == 32.0
    reset_integration!(det_buffered)
    @test readout_ready(det_buffered)
    @test det_buffered.state.integrated_time == 0.0
    metadata_buffered = detector_export_metadata(det_buffered)
    @test metadata_buffered.output_size == size(output_frame(det_buffered))
    @test metadata_buffered.psf_sampling == 1
    @test metadata_buffered.binning == 1

    det_background_flux = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        background_flux=2.0)
    frame_background_flux = capture!(det_background_flux, zeros(4, 4); rng=MersenneTwister(2))
    @test sum(frame_background_flux) > 0

    det_background_map = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        background_map=fill(1.0, 4, 4))
    frame_background_map = capture!(det_background_map, zeros(4, 4); rng=MersenneTwister(2))
    @test frame_background_map == fill(-1.0, 4, 4)

    cube = Array{Float64}(undef, 2, 4, 4)
    cube[1, :, :] .= fill(1.0, 4, 4)
    cube[2, :, :] .= fill(2.0, 4, 4)
    scratch = similar(cube)
    det_stack = Detector(integration_time=1.0, noise=NoiseNone(), qe=0.5, binning=1)
    AdaptiveOpticsSim.capture_stack!(det_stack, cube, scratch; rng=MersenneTwister(10))
    @test cube[1, :, :] ≈ fill(0.5, 4, 4)
    @test cube[2, :, :] ≈ fill(1.0, 4, 4)

    psf = reshape(Float64.(1:256), 16, 16)
    det_fused = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, psf_sampling=2, binning=2)
    frame_fused = copy(AdaptiveOpticsSim.fill_frame!(det_fused, psf, 1.0))
    manual_mid = zeros(Float64, 8, 8)
    manual_out = zeros(Float64, 4, 4)
    AdaptiveOpticsSim.bin2d!(manual_mid, psf, 2)
    AdaptiveOpticsSim.bin2d!(manual_out, manual_mid, 2)
    @test frame_fused == manual_out
end

@testset "Asterism PSF" begin
    tel = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    src1 = Source(band=:I, magnitude=0.0, coordinates=(0.0, 0.0))
    src2 = Source(band=:I, magnitude=0.0, coordinates=(1.0, 90.0))
    @test coordinates_xy_arcsec(src1) == (0.0, 0.0)
    @test coordinates_xy_arcsec(src2)[1] ≈ 0.0 atol=1e-12
    @test coordinates_xy_arcsec(src2)[2] ≈ 1.0
    ast = Asterism([src1, src2])
    psf = compute_psf!(tel, ast; zero_padding=2)
    @test size(tel.state.psf_stack, 3) == 2
    @test size(psf) == (32, 32)
    @test sum(psf) >= sum(@view tel.state.psf_stack[:, :, 1])
end

@testset "Polychromatic WFS" begin
    tel = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    src = Source(band=:I, magnitude=0.0)
    λ0 = wavelength(src)
    bundle_single = SpectralBundle([SpectralSample(λ0, 1.0)])
    bundle_broad = SpectralBundle([SpectralSample(0.9 * λ0, 0.4), SpectralSample(1.1 * λ0, 0.6)])
    poly_single = with_spectrum(src, bundle_single)
    poly_broad = with_spectrum(src, bundle_broad)

    @test sum(sample.weight for sample in bundle_broad) ≈ 1.0
    @test weighted_wavelength(bundle_broad) ≈ (0.9 * λ0 * 0.4 + 1.1 * λ0 * 0.6)
    @test has_spectral_bundle(poly_broad)
    @test is_polychromatic(poly_broad)
    @test !is_polychromatic(poly_single)
    @test spectral_reference_source(poly_broad) === src

    zb = ZernikeBasis(tel, 5)
    compute_zernike!(zb, tel)
    focus = @view zb.modes[:, :, 5]
    @. tel.state.opd = 5e-8 * focus

    sh_mono = ShackHartmann(tel; n_subap=8, mode=Diffractive())
    sh_single = ShackHartmann(tel; n_subap=8, mode=Diffractive())
    sh_broad = ShackHartmann(tel; n_subap=8, mode=Diffractive())

    mono_slopes = copy(measure!(sh_mono, tel, src))
    single_slopes = copy(measure!(sh_single, tel, poly_single))
    broad_slopes_1 = copy(measure!(sh_broad, tel, poly_broad))
    broad_slopes_2 = copy(measure!(sh_broad, tel, poly_broad))

    @test single_slopes ≈ mono_slopes atol=1e-10 rtol=1e-10
    @test broad_slopes_1 ≈ broad_slopes_2 atol=1e-10 rtol=1e-10
    @test norm(broad_slopes_1 - mono_slopes) > 1e-8
    @test supports_stacked_sources(sh_broad, poly_broad)

    pyr_mono = PyramidWFS(tel; n_subap=8, mode=Diffractive(), modulation=1.0)
    pyr_single = PyramidWFS(tel; n_subap=8, mode=Diffractive(), modulation=1.0)
    pyr_broad = PyramidWFS(tel; n_subap=8, mode=Diffractive(), modulation=1.0)

    mono_pyr = copy(measure!(pyr_mono, tel, src))
    single_pyr = copy(measure!(pyr_single, tel, poly_single))
    broad_pyr_1 = copy(measure!(pyr_broad, tel, poly_broad))
    broad_pyr_2 = copy(measure!(pyr_broad, tel, poly_broad))

    @test single_pyr ≈ mono_pyr atol=1e-10 rtol=1e-10
    @test broad_pyr_1 ≈ broad_pyr_2 atol=1e-10 rtol=1e-10
    @test norm(broad_pyr_1 - mono_pyr) > 1e-8
    @test supports_stacked_sources(pyr_broad, poly_broad)

    det = Detector(noise=NoiseNone(), binning=1)
    spectral_frame = measure!(sh_broad, tel, poly_broad, det)
    @test size(sh_broad.state.detector_noise_cube) == size(sh_broad.state.spot_cube)
    @test spectral_frame ≈ broad_slopes_1 atol=1e-10 rtol=1e-10

    pyr_det = Detector(noise=NoiseNone(), binning=1)
    pyr_det_slopes = measure!(pyr_broad, tel, poly_broad, pyr_det)
    @test size(output_frame(pyr_det)) == size(pyr_broad.state.camera_frame)
    @test pyr_det_slopes ≈ broad_pyr_1 atol=1e-10 rtol=1e-10

    fill!(tel.state.opd, 0.0)
end

@testset "Extended-source WFS" begin
    tel = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    src = Source(band=:I, magnitude=0.0)
    point_model = PointCloudSourceModel([(0.0, 0.0)], [1.0])
    gaussian_model = GaussianDiskSourceModel(sigma_arcsec=0.35, n_side=5)
    image_model = SampledImageSourceModel([0.0 1.0 0.0; 1.0 2.0 1.0; 0.0 1.0 0.0], pixel_scale_arcsec=0.2)
    ext_point = with_extended_source(src, point_model)
    ext_gauss = with_extended_source(src, gaussian_model)
    ext_image = with_extended_source(src, image_model)

    point_ast = extended_source_asterism(ext_point)
    image_ast = extended_source_asterism(ext_image)
    @test has_extended_source_model(ext_gauss)
    @test !has_extended_source_model(src)
    @test length(point_ast) == 1
    @test length(image_ast) == 5
    @test AdaptiveOpticsSim.photon_flux(point_ast.sources[1]) ≈ AdaptiveOpticsSim.photon_flux(src)
    @test sum(AdaptiveOpticsSim.photon_flux(sample) for sample in image_ast.sources) ≈ AdaptiveOpticsSim.photon_flux(src)

    zb = ZernikeBasis(tel, 5)
    compute_zernike!(zb, tel)
    focus = @view zb.modes[:, :, 5]
    @. tel.state.opd = 5e-8 * focus

    sh_point = ShackHartmann(tel; n_subap=8, mode=Diffractive())
    sh_ext_point = ShackHartmann(tel; n_subap=8, mode=Diffractive())
    sh_ext = ShackHartmann(tel; n_subap=8, mode=Diffractive())
    point_slopes = copy(measure!(sh_point, tel, src))
    ext_point_slopes = copy(measure!(sh_ext_point, tel, ext_point))
    point_peak = AdaptiveOpticsSim.sampled_spots_peak!(sh_point, tel, src)
    ext_peak = AdaptiveOpticsSim.sampled_spots_peak!(sh_ext, tel, ext_gauss)
    ext_slopes_1 = copy(measure!(sh_ext, tel, ext_gauss))
    ext_slopes_2 = copy(measure!(sh_ext, tel, ext_gauss))

    @test ext_point_slopes ≈ point_slopes atol=1e-10 rtol=1e-10
    @test ext_slopes_1 ≈ ext_slopes_2 atol=1e-10 rtol=1e-10
    @test ext_peak < point_peak
    @test norm(sh_ext.state.spot_cube - sh_point.state.spot_cube) > 1e-8
    @test supports_stacked_sources(sh_ext, ext_gauss)

    pyr_point = PyramidWFS(tel; n_subap=8, mode=Diffractive(), modulation=1.0)
    pyr_ext_point = PyramidWFS(tel; n_subap=8, mode=Diffractive(), modulation=1.0)
    pyr_ext = PyramidWFS(tel; n_subap=8, mode=Diffractive(), modulation=1.0)
    pyr_point_slopes = copy(measure!(pyr_point, tel, src))
    pyr_ext_point_slopes = copy(measure!(pyr_ext_point, tel, ext_point))
    pyr_ext_slopes_1 = copy(measure!(pyr_ext, tel, ext_gauss))
    pyr_ext_slopes_2 = copy(measure!(pyr_ext, tel, ext_gauss))

    @test pyr_ext_point_slopes ≈ pyr_point_slopes atol=1e-10 rtol=1e-10
    @test pyr_ext_slopes_1 ≈ pyr_ext_slopes_2 atol=1e-10 rtol=1e-10
    @test norm(pyr_ext.state.intensity - pyr_point.state.intensity) > 1e-10
    @test supports_stacked_sources(pyr_ext, ext_gauss)

    det = Detector(noise=NoiseNone(), binning=1)
    sh_det_slopes = measure!(sh_ext, tel, ext_gauss, det)
    @test size(sh_ext.state.detector_noise_cube) == size(sh_ext.state.spot_cube)
    @test sh_det_slopes ≈ ext_slopes_1 atol=1e-10 rtol=1e-10

    pyr_det = Detector(noise=NoiseNone(), binning=1)
    pyr_det_slopes = measure!(pyr_ext, tel, ext_gauss, pyr_det)
    @test size(output_frame(pyr_det)) == size(pyr_ext.state.camera_frame)
    @test pyr_det_slopes ≈ pyr_ext_slopes_1 atol=1e-10 rtol=1e-10

    fill!(tel.state.opd, 0.0)
end

@testset "Pupil masks and misregistration" begin
    tel = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    base_sum = sum(tel.state.pupil)
    apply_spiders!(tel; thickness=0.5, angles=[0.0, 90.0])
    @test sum(tel.state.pupil) < base_sum

    custom = trues(16, 16)
    custom[:, 9:end] .= false
    set_pupil!(tel, custom)
    @test sum(tel.state.pupil) == sum(custom)
    @test tel.state.pupil_reflectivity == Float64.(custom)

    reflectivity = fill(0.5, 16, 16)
    set_pupil_reflectivity!(tel, reflectivity)
    @test tel.state.pupil_reflectivity[:, 1:8] == fill(0.5, 16, 8)
    @test tel.state.pupil_reflectivity[:, 9:end] == fill(0.0, 16, 8)

    tel2 = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    dm1 = DeformableMirror(tel2; n_act=2, influence_width=0.3)
    mis = Misregistration(shift_x=0.1, shift_y=0.0, rotation_deg=5.0, T=Float64)
    @test rotation_deg(mis) ≈ 5.0
    @test rotation_rad(mis) ≈ deg2rad(5.0)
    dm2 = DeformableMirror(tel2; n_act=2, influence_width=0.3, misregistration=mis)
    @test dm1.state.modes != dm2.state.modes
end

@testset "Pyramid, BioEdge, and LGS" begin
    tel = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    for i in 1:tel.params.resolution, j in 1:tel.params.resolution
        tel.state.opd[i, j] = i
    end

    pyr = PyramidWFS(tel; n_subap=4, modulation=1.0)
    pyr_slopes = measure!(pyr, tel)
    @test length(pyr_slopes) == 2 * 4 * 4

    bio = BioEdgeWFS(tel; n_subap=4)
    bio_slopes = measure!(bio, tel)
    @test length(bio_slopes) == 2 * 4 * 4

    sh = ShackHartmann(tel; n_subap=4)
    ngs = Source(band=:I, magnitude=0.0)
    lgs = LGSSource(elongation_factor=2.0)
    slopes_ngs = measure!(sh, tel, ngs)
    slopes_lgs = measure!(sh, tel, lgs)
    n = sh.params.n_subap * sh.params.n_subap
    @test slopes_lgs[n+1:end] ≈ slopes_ngs[n+1:end] .* 2.0
end

@testset "Diffractive WFS" begin
    tel = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    for i in 1:tel.params.resolution, j in 1:tel.params.resolution
        tel.state.opd[i, j] = i
    end
    ngs = Source(band=:I, magnitude=0.0)
    lgs = LGSSource(elongation_factor=1.5)

    sh = ShackHartmann(tel; n_subap=4, mode=Diffractive())
    @test_throws InvalidConfiguration measure!(sh, tel)
    sh_slopes = measure!(sh, tel, ngs)
    @test length(sh_slopes) == 2 * 4 * 4
    @test all(isfinite, sh_slopes)
    sh_lgs = measure!(sh, tel, lgs)
    @test all(isfinite, sh_lgs)

    na_profile = [80000.0 90000.0 100000.0; 0.2 0.6 0.2]
    lgs_profile = LGSSource(elongation_factor=1.2, na_profile=na_profile, fwhm_spot_up=1.0)
    sh_profile = ShackHartmann(tel; n_subap=4, mode=Diffractive())
    sh_profile_slopes = measure!(sh_profile, tel, lgs_profile)
    @test all(isfinite, sh_profile_slopes)

    sh_sampled = ShackHartmann(tel; n_subap=4, mode=Diffractive(), pixel_scale=0.06, n_pix_subap=8)
    sh_sampled_slopes = measure!(sh_sampled, tel, ngs)
    @test length(sh_sampled_slopes) == 2 * 4 * 4

    pyr_sampled = PyramidWFS(tel; n_subap=4, mode=Diffractive(), n_pix_separation=4, binning=2)
    pyr_sampled_slopes = measure!(pyr_sampled, tel, ngs)
    @test length(pyr_sampled_slopes) == 2 * count(pyr_sampled.state.valid_i4q)
    pyr_intensity = reshape(Float64.(1:size(tel.state.opd, 1)^2), size(tel.state.opd))
    pyr_frame = copy(AdaptiveOpticsSim.sample_pyramid_intensity!(pyr_sampled, tel, pyr_intensity))
    pyr_camera = zeros(Float64, 4, 4)
    pyr_manual = zeros(Float64, 2, 2)
    AdaptiveOpticsSim.bin2d!(pyr_camera, pyr_intensity, 8)
    AdaptiveOpticsSim.bin2d!(pyr_manual, pyr_camera, 2)
    @test pyr_frame == pyr_manual

    bio_sampled = BioEdgeWFS(tel; n_subap=4, mode=Diffractive(), binning=2)
    bio_sampled_slopes = measure!(bio_sampled, tel, ngs)
    @test length(bio_sampled_slopes) == 2 * count(bio_sampled.state.valid_i4q)
    bio_intensity = reshape(Float64.(1:size(tel.state.opd, 1)^2), size(tel.state.opd))
    bio_frame = copy(AdaptiveOpticsSim.sample_bioedge_intensity!(bio_sampled, tel, bio_intensity))
    bio_camera = zeros(Float64, 4, 4)
    bio_manual = similar(bio_frame)
    AdaptiveOpticsSim.bin2d!(bio_camera, bio_intensity, 8)
    AdaptiveOpticsSim.bin2d!(bio_manual, bio_camera, div(size(bio_camera, 1), size(bio_frame, 1)))
    @test bio_frame == bio_manual

    pyr_profile = PyramidWFS(tel; n_subap=4, mode=Diffractive())
    pyr_profile_slopes = measure!(pyr_profile, tel, lgs_profile)
    @test all(isfinite, pyr_profile_slopes)

    bio_profile = BioEdgeWFS(tel; n_subap=4, mode=Diffractive())
    bio_profile_slopes = measure!(bio_profile, tel, lgs_profile)
    @test all(isfinite, bio_profile_slopes)

    pyr = PyramidWFS(tel; n_subap=4, mode=Diffractive())
    @test_throws InvalidConfiguration measure!(pyr, tel)
    pyr_slopes = measure!(pyr, tel, ngs)
    @test length(pyr_slopes) == 2 * 4 * 4

    bio = BioEdgeWFS(tel; n_subap=4, mode=Diffractive())
    @test_throws InvalidConfiguration measure!(bio, tel)
    bio_slopes = measure!(bio, tel, ngs)
    @test length(bio_slopes) == 2 * 4 * 4

    det = Detector(noise=NoiseNone(), binning=1)
    sh_det = ShackHartmann(tel; n_subap=4, mode=Diffractive())
    sh_det_slopes = measure!(sh_det, tel, ngs, det)
    @test length(sh_det_slopes) == 2 * 4 * 4
    pyr_det = PyramidWFS(tel; n_subap=4, mode=Diffractive())
    pyr_det_slopes = measure!(pyr_det, tel, ngs, det)
    @test length(pyr_det_slopes) == 2 * 4 * 4

    ast = Asterism([ngs, Source(band=:I, magnitude=0.0, coordinates=(0.0, 0.0))])
    sh_ast = ShackHartmann(tel; n_subap=4, mode=Diffractive())
    sh_ast_slopes = copy(measure!(sh_ast, tel, ast))
    @test length(sh_ast_slopes) == 2 * 4 * 4
    sh_ast_serial = ShackHartmann(tel; n_subap=4, mode=Diffractive())
    AdaptiveOpticsSim.prepare_sampling!(sh_ast_serial, tel, ast.sources[1])
    AdaptiveOpticsSim.ensure_sh_calibration!(sh_ast_serial, tel, ast.sources[1])
    fill!(sh_ast_serial.state.detector_noise_cube, zero(eltype(sh_ast_serial.state.detector_noise_cube)))
    for src in ast.sources
        AdaptiveOpticsSim.sampled_spots_peak!(sh_ast_serial, tel, src)
        sh_ast_serial.state.detector_noise_cube .+= sh_ast_serial.state.spot_cube
    end
    copyto!(sh_ast_serial.state.spot_cube, sh_ast_serial.state.detector_noise_cube)
    sh_ast_serial_peak = maximum(sh_ast_serial.state.spot_cube)
    AdaptiveOpticsSim.sh_signal_from_spots!(sh_ast_serial, sh_ast_serial_peak, slope_extraction_model(sh_ast_serial))
    AdaptiveOpticsSim.subtract_reference_and_scale!(sh_ast_serial)
    sh_ast_serial_slopes = copy(sh_ast_serial.state.slopes)
    @test sh_ast_slopes ≈ sh_ast_serial_slopes
    mixed_ngs = Source(wavelength=wavelength(lgs_profile), magnitude=0.0, coordinates=(0.0, 0.0))
    mixed_ast = Asterism([mixed_ngs, lgs_profile])
    sh_mixed_det = ShackHartmann(tel; n_subap=4, mode=Diffractive())
    sh_mixed_det_slopes = copy(measure!(sh_mixed_det, tel, mixed_ast, det; rng=MersenneTwister(14)))
    sh_mixed_det_frame = copy(sh_mixed_det.state.spot_cube)
    sh_mixed_det_manual = ShackHartmann(tel; n_subap=4, mode=Diffractive())
    AdaptiveOpticsSim.prepare_sampling!(sh_mixed_det_manual, tel, mixed_ast.sources[1])
    AdaptiveOpticsSim.ensure_sh_calibration!(sh_mixed_det_manual, tel, mixed_ast.sources[1])
    fill!(sh_mixed_det_manual.state.detector_noise_cube, zero(eltype(sh_mixed_det_manual.state.detector_noise_cube)))
    for src in mixed_ast.sources
        AdaptiveOpticsSim.sampled_spots_peak!(sh_mixed_det_manual, tel, src, det, MersenneTwister(14))
        sh_mixed_det_manual.state.detector_noise_cube .+= sh_mixed_det_manual.state.spot_cube
    end
    copyto!(sh_mixed_det_manual.state.spot_cube, sh_mixed_det_manual.state.detector_noise_cube)
    @test sh_mixed_det_frame ≈ sh_mixed_det_manual.state.spot_cube
    @test length(sh_mixed_det_slopes) == 2 * 4 * 4
    pyr_ast = PyramidWFS(tel; n_subap=4, mode=Diffractive())
    pyr_ast_slopes = copy(measure!(pyr_ast, tel, ast))
    @test length(pyr_ast_slopes) == 2 * 4 * 4
    pyr_ast_serial = PyramidWFS(tel; n_subap=4, mode=Diffractive())
    AdaptiveOpticsSim.ensure_pyramid_calibration!(pyr_ast_serial, tel, ast.sources[1])
    pyr_ast_stack = @view AdaptiveOpticsSim.ensure_pyramid_asterism_stack!(pyr_ast_serial, length(ast.sources))[:, :, 1:length(ast.sources)]
    fill!(pyr_ast_serial.state.intensity, zero(eltype(pyr_ast_serial.state.intensity)))
    for (src_idx, src) in pairs(ast.sources)
        AdaptiveOpticsSim.pyramid_intensity!(@view(pyr_ast_stack[:, :, src_idx]), pyr_ast_serial, tel, src)
        pyr_ast_serial.state.intensity .+= @view(pyr_ast_stack[:, :, src_idx])
    end
    pyr_ast_intensity = AdaptiveOpticsSim.sample_pyramid_intensity!(pyr_ast_serial, tel, pyr_ast_serial.state.intensity)
    AdaptiveOpticsSim.pyramid_signal!(pyr_ast_serial, tel, pyr_ast_intensity)
    @. pyr_ast_serial.state.slopes *= pyr_ast_serial.state.optical_gain
    @test pyr_ast_slopes ≈ pyr_ast_serial.state.slopes
    bio_ast = BioEdgeWFS(tel; n_subap=4, mode=Diffractive())
    bio_ast_slopes = copy(measure!(bio_ast, tel, ast))
    @test length(bio_ast_slopes) == 2 * 4 * 4
    bio_ast_serial = BioEdgeWFS(tel; n_subap=4, mode=Diffractive())
    AdaptiveOpticsSim.ensure_bioedge_calibration!(bio_ast_serial, tel, ast.sources[1])
    fill!(bio_ast_serial.state.binned_intensity, zero(eltype(bio_ast_serial.state.binned_intensity)))
    for src in ast.sources
        AdaptiveOpticsSim.bioedge_intensity!(bio_ast_serial.state.intensity, bio_ast_serial, tel, src)
        bio_ast_serial.state.binned_intensity .+= bio_ast_serial.state.intensity
    end
    bio_ast_intensity = AdaptiveOpticsSim.sample_bioedge_intensity!(bio_ast_serial, tel, bio_ast_serial.state.binned_intensity)
    AdaptiveOpticsSim.bioedge_signal!(bio_ast_serial, tel, bio_ast_intensity)
    @test bio_ast_slopes ≈ bio_ast_serial.state.slopes

    pyr_ast_det = PyramidWFS(tel; n_subap=4, mode=Diffractive())
    pyr_ast_det_slopes = copy(measure!(pyr_ast_det, tel, ast, det))
    pyr_ast_det_serial = PyramidWFS(tel; n_subap=4, mode=Diffractive())
    AdaptiveOpticsSim.ensure_pyramid_calibration!(pyr_ast_det_serial, tel, ast.sources[1])
    pyr_ast_det_stack = @view AdaptiveOpticsSim.ensure_pyramid_asterism_stack!(pyr_ast_det_serial, length(ast.sources))[:, :, 1:length(ast.sources)]
    fill!(pyr_ast_det_serial.state.intensity, zero(eltype(pyr_ast_det_serial.state.intensity)))
    for (src_idx, src) in pairs(ast.sources)
        AdaptiveOpticsSim.pyramid_intensity!(@view(pyr_ast_det_stack[:, :, src_idx]), pyr_ast_det_serial, tel, src)
        pyr_ast_det_serial.state.intensity .+= @view(pyr_ast_det_stack[:, :, src_idx])
    end
    pyr_ast_det_intensity = AdaptiveOpticsSim.sample_pyramid_intensity!(pyr_ast_det_serial, tel, pyr_ast_det_serial.state.intensity)
    pyr_ast_det_frame = capture!(det, pyr_ast_det_intensity; rng=MersenneTwister(12))
    AdaptiveOpticsSim.resize_pyramid_signal_buffers!(pyr_ast_det_serial, size(pyr_ast_det_frame, 1))
    AdaptiveOpticsSim.pyramid_signal!(pyr_ast_det_serial, tel, pyr_ast_det_frame)
    @. pyr_ast_det_serial.state.slopes *= pyr_ast_det_serial.state.optical_gain
    @test pyr_ast_det_slopes ≈ pyr_ast_det_serial.state.slopes

    bio_ast_det = BioEdgeWFS(tel; n_subap=4, mode=Diffractive())
    bio_ast_det_slopes = copy(measure!(bio_ast_det, tel, ast, det; rng=MersenneTwister(13)))
    bio_ast_det_serial = BioEdgeWFS(tel; n_subap=4, mode=Diffractive())
    AdaptiveOpticsSim.ensure_bioedge_calibration!(bio_ast_det_serial, tel, ast.sources[1])
    fill!(bio_ast_det_serial.state.binned_intensity, zero(eltype(bio_ast_det_serial.state.binned_intensity)))
    for src in ast.sources
        AdaptiveOpticsSim.bioedge_intensity!(bio_ast_det_serial.state.intensity, bio_ast_det_serial, tel, src)
        bio_ast_det_serial.state.binned_intensity .+= bio_ast_det_serial.state.intensity
    end
    bio_ast_det_intensity = AdaptiveOpticsSim.sample_bioedge_intensity!(bio_ast_det_serial, tel, bio_ast_det_serial.state.binned_intensity)
    bio_ast_det_frame = capture!(det, bio_ast_det_intensity; rng=MersenneTwister(13))
    AdaptiveOpticsSim.resize_bioedge_signal_buffers!(bio_ast_det_serial, size(bio_ast_det_frame, 1))
    AdaptiveOpticsSim.bioedge_signal!(bio_ast_det_serial, tel, bio_ast_det_frame)
    @test bio_ast_det_slopes ≈ bio_ast_det_serial.state.slopes
end

@testset "Shack-Hartmann subapertures" begin
    tel = Telescope(resolution=24, diameter=8.0, sampling_time=1e-3, central_obstruction=0.1)
    src = Source(band=:I, magnitude=0.0)
    sh = ShackHartmann(tel; n_subap=6, mode=Diffractive(), pixel_scale=0.06, n_pix_subap=8, threshold_cog=0.02)

    layout = subaperture_layout(sh)
    calibration = subaperture_calibration(sh)
    @test layout isa SubapertureLayout
    @test calibration isa SubapertureCalibration
    @test layout.n_subap == 6
    @test layout.subap_pixels == 4
    @test layout.pitch_m ≈ tel.params.diameter / 6
    @test !calibration.calibrated
    @test slope_extraction_model(sh) isa CenterOfGravityExtraction
    @test slope_extraction_model(sh).threshold ≈ 0.02
    @test n_valid_subapertures(layout) == count(layout.valid_mask_host)

    prepare_runtime_wfs!(sh, tel, src)
    @test calibration.calibrated
    @test calibration.slopes_units == sh.state.slopes_units
    @test calibration.wavelength == sh.state.calibration_wavelength
    @test calibration.signature == sh.state.calibration_signature
    @test calibration.reference_signal_2d === sh.state.reference_signal_2d
    @test calibration.reference_signal_host === sh.state.reference_signal_host
    @test length(valid_subaperture_indices(layout)) == n_valid_subapertures(layout)

    slopes = measure!(sh, tel, src)
    @test all(isfinite, slopes)
    meta = AdaptiveOpticsSim.wfs_output_metadata(sh)
    @test meta.n_valid_subap == n_valid_subapertures(layout)
    @test meta.subap_pixels == layout.subap_pixels
    @test meta.calibrated

    dm = DeformableMirror(tel; n_act=5)
    imat = interaction_matrix(dm, sh, tel, src; amplitude=1e-8)
    @test size(imat.matrix, 1) == length(sh.state.slopes)
    @test size(imat.matrix, 2) == length(dm.state.coefs)
end

@testset "Zernike WFS" begin
    tel = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    src = Source(band=:I, magnitude=0.0)
    wfs = ZernikeWFS(tel; n_subap=8, diffraction_padding=2)

    @test size(wfs.state.camera_frame) == (8, 8)
    @test length(wfs.state.slopes) == count(wfs.state.valid_mask)
    @test_throws InvalidConfiguration measure!(wfs, tel)
    @test_throws InvalidConfiguration measure!(wfs, tel, Asterism([src, Source(band=:I, magnitude=0.0)]))

    flat_slopes = copy(measure!(wfs, tel, src))
    @test wfs.state.calibrated
    @test all(isfinite, flat_slopes)
    @test all(>=(0.0), wfs.state.camera_frame)
    @test flat_slopes ≈ zero.(flat_slopes) atol=1e-10

    det = Detector(noise=NoiseNone(), binning=1)
    det_slopes = copy(measure!(wfs, tel, src, det))
    @test det_slopes ≈ flat_slopes atol=1e-10
    @test size(output_frame(det)) == size(wfs.state.camera_frame)

    zb = ZernikeBasis(tel, 5)
    compute_zernike!(zb, tel)
    focus = @view zb.modes[:, :, 5]
    @. tel.state.opd = 5e-8 * focus
    slopes_plus = copy(measure!(wfs, tel, src))
    @. tel.state.opd = -5e-8 * focus
    slopes_minus = copy(measure!(wfs, tel, src))
    fill!(tel.state.opd, 0.0)

    @test norm(slopes_plus) > 1e-6
    @test norm(slopes_minus) > 1e-6
    @test dot(slopes_plus, slopes_minus) < 0
end

@testset "Curvature WFS" begin
    tel = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    src = Source(band=:I, magnitude=0.0)
    wfs = CurvatureWFS(tel; n_subap=8, defocus_rms_nm=500.0)

    @test size(wfs.state.camera_frame) == (16, 8)
    @test length(wfs.state.slopes) == 64
    @test_throws InvalidConfiguration measure!(wfs, tel)
    @test_throws InvalidConfiguration measure!(wfs, tel, Asterism([src, Source(band=:I, magnitude=0.0)]))

    flat_slopes = copy(measure!(wfs, tel, src))
    @test wfs.state.calibrated
    @test all(isfinite, flat_slopes)
    @test all(>=(0.0), wfs.state.camera_frame)
    @test flat_slopes ≈ zero.(flat_slopes) atol=1e-10

    det = Detector(noise=NoiseNone(), binning=1)
    det_slopes = copy(measure!(wfs, tel, src, det))
    @test det_slopes ≈ flat_slopes atol=1e-10
    @test size(output_frame(det)) == size(wfs.state.camera_frame)

    zb = ZernikeBasis(tel, 5)
    compute_zernike!(zb, tel)
    focus = @view zb.modes[:, :, 5]
    @. tel.state.opd = 5e-8 * focus
    slopes_plus = copy(measure!(wfs, tel, src))
    @. tel.state.opd = -5e-8 * focus
    slopes_minus = copy(measure!(wfs, tel, src))
    fill!(tel.state.opd, 0.0)

    @test norm(slopes_plus) > 1e-6
    @test norm(slopes_minus) > 1e-6
    @test dot(slopes_plus, slopes_minus) < 0

    counting = CurvatureWFS(tel; n_subap=8, defocus_rms_nm=500.0, readout_model=CurvatureCountingReadout())
    counting_flat = copy(measure!(counting, tel, src))
    @test size(counting.state.camera_frame) == (2, 64)
    @test counting_flat ≈ zero.(counting_flat) atol=1e-10
    @test_throws InvalidConfiguration measure!(counting, tel, src, det)
    apd = APDDetector(integration_time=1.0, qe=1.0, gain=1.0, dark_count_rate=0.0, noise=NoiseNone())
    counting_apd = copy(measure!(counting, tel, src, apd))
    @test counting_apd ≈ counting_flat atol=1e-10
    @test detector_export_metadata(apd).readout.output_size == size(counting.state.camera_frame)
    apd_dead = APDDetector(integration_time=1.0, qe=1.0, gain=1.0, dark_count_rate=0.0,
        noise=NoiseNone(), dead_time_model=NonParalyzableDeadTime(0.25))
    counting_dead = copy(measure!(counting, tel, src, apd_dead))
    @test counting_dead ≈ counting_flat atol=1e-10
    @test_throws InvalidConfiguration CurvatureWFS(tel; n_subap=8, readout_model=CurvatureCountingReadout(),
        readout_pixels_per_subap=2)

    response = CurvatureBranchResponse(T=Float64, plus_throughput=1.2, minus_throughput=0.8,
        plus_background=5.0, minus_background=1.0)
    imbalanced = CurvatureWFS(tel; n_subap=8, defocus_rms_nm=500.0, branch_response=response)
    imbalanced_flat = copy(measure!(imbalanced, tel, src))
    @test imbalanced_flat ≈ zero.(imbalanced_flat) atol=1e-10
    plus_mean = mean(@view imbalanced.state.camera_frame[1:imbalanced.params.n_subap, :])
    minus_mean = mean(@view imbalanced.state.camera_frame[imbalanced.params.n_subap+1:end, :])
    @test plus_mean > minus_mean
    @test_throws InvalidConfiguration CurvatureBranchResponse(plus_throughput=-1.0)

    oversampled = CurvatureWFS(tel; n_subap=8, readout_crop_resolution=16, readout_pixels_per_subap=2)
    oversampled_flat = copy(measure!(oversampled, tel, src))
    @test size(oversampled.state.camera_frame) == (32, 16)
    @test size(oversampled.state.frame_plus) == (16, 16)
    @test size(oversampled.state.reduced_plus) == (8, 8)
    @test oversampled_flat ≈ zero.(oversampled_flat) atol=1e-10
    @test_throws InvalidConfiguration CurvatureWFS(tel; n_subap=8, readout_crop_resolution=18, readout_pixels_per_subap=2)

    atm = MultiLayerAtmosphere(tel;
        r0=0.2,
        L0=25.0,
        fractional_cn2=[0.7, 0.3],
        wind_speed=[8.0, 4.0],
        wind_direction=[0.0, 90.0],
        altitude=[0.0, 5000.0],
    )
    advance!(atm, tel; rng=MersenneTwister(3))
    atm_slopes = copy(measure!(wfs, tel, src, atm))
    @test all(isfinite, atm_slopes)
    @test norm(atm_slopes) > 0

    det_atm = Detector(noise=NoiseNone(), binning=1)
    det_atm_slopes = copy(measure!(wfs, tel, src, atm, det_atm))
    @test all(isfinite, det_atm_slopes)

    ast = Asterism([src, Source(band=:I, magnitude=0.0, coordinates=(1.0, 90.0))])
    ast_slopes = copy(measure!(wfs, tel, ast, atm))
    @test length(ast_slopes) == length(wfs.state.slopes)
    @test all(isfinite, ast_slopes)
    @test norm(ast_slopes) > 0
end

@testset "OOPAO parity knobs" begin
    tel = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    src = Source(band=:I, magnitude=0.0)
    for i in 1:tel.params.resolution, j in 1:tel.params.resolution
        tel.state.opd[i, j] = i + j / 10
    end

    sh_plain = ShackHartmann(tel; n_subap=4, mode=Diffractive(), pixel_scale=0.06, n_pix_subap=8)
    sh_shift = ShackHartmann(tel; n_subap=4, mode=Diffractive(), pixel_scale=0.06, n_pix_subap=8,
        half_pixel_shift=true)
    sh_thresh = ShackHartmann(tel; n_subap=4, mode=Diffractive(), pixel_scale=0.06, n_pix_subap=8,
        threshold_cog=0.2)
    @test measure!(sh_plain, tel, src) != measure!(sh_shift, tel, src)
    @test measure!(sh_plain, tel, src) != measure!(sh_thresh, tel, src)

    pyr_auto = PyramidWFS(tel; n_subap=4, mode=Diffractive(), modulation=1.0)
    @test size(pyr_auto.state.modulation_phases, 3) == 8

    pyr_path = PyramidWFS(tel; n_subap=4, mode=Diffractive(), modulation=0.0,
        user_modulation_path=((1.0, 0.0), (0.0, 1.0)))
    @test size(pyr_path.state.modulation_phases, 3) == 2

    pyr_default = PyramidWFS(tel; n_subap=4, mode=Diffractive(), modulation=1.0)
    pyr_rooftop = PyramidWFS(tel; n_subap=4, mode=Diffractive(), modulation=1.0,
        rooftop=0.5, theta_rotation=0.2)
    pyr_old = PyramidWFS(tel; n_subap=4, mode=Diffractive(), modulation=1.0, old_mask=true)
    @test pyr_default.state.pyramid_mask != pyr_rooftop.state.pyramid_mask
    @test pyr_default.state.pyramid_mask != pyr_old.state.pyramid_mask

    bio_plain = BioEdgeWFS(tel; n_subap=4, mode=Diffractive(), modulation=1.0)
    bio_gray = BioEdgeWFS(tel; n_subap=4, mode=Diffractive(), modulation=1.0,
        grey_width=0.5, grey_length=1.0)
    amps = real.(bio_gray.state.bioedge_masks[:, :, 1])
    @test any(x -> 0 < x < 1, amps)
    @test bio_plain.state.bioedge_masks != bio_gray.state.bioedge_masks

    bio_gray_slopes = measure!(bio_gray, tel, src)
    @test length(bio_gray_slopes) == 2 * 4 * 4
    @test all(isfinite, bio_gray_slopes)
end
@testset "Calibration vault and modal basis" begin
    D = rand(4, 3)
    vault = CalibrationVault(D)
    @test size(vault.M) == (3, 4)
    @test vault.cond > 0
    @test vault.effective_rank == 3
    vault_trunc = with_truncation(vault, 1)
    @test vault_trunc.n_trunc == 1

    D_sing = [1.0 0.0; 0.0 1e-12]
    vault_exact = CalibrationVault(D_sing; policy=ExactPseudoInverse())
    vault_tsvd = CalibrationVault(D_sing; policy=TSVDInverse(rtol=1e-9))
    @test vault_exact.effective_rank == 2
    @test vault_tsvd.effective_rank == 1
    @test maximum(abs, vault_exact.M .- vault_tsvd.M) > 0
    @test CalibrationVault(D_sing).policy isa TSVDInverse

    imat = InteractionMatrix(D_sing, 0.1)
    recon_exact = ModalReconstructor(imat; policy=ExactPseudoInverse())
    recon_tsvd = ModalReconstructor(imat; policy=TSVDInverse(rtol=1e-9))
    @test recon_exact.effective_rank == 2
    @test recon_tsvd.effective_rank == 1
    @test ModalReconstructor(imat).policy isa TSVDInverse

    tel = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    dm = DeformableMirror(tel; n_act=2, influence_width=0.4)
    basis = modal_basis(dm, tel; n_modes=2)
    @test size(basis.M2C, 2) == 2
    opd = rand(16, 16)
    fit, corr, turb = fitting_error(opd, basis.projector, basis.basis)
    @test size(fit) == size(opd)
    @test size(corr) == size(opd)
    @test size(turb) == size(opd)

    atm = KolmogorovAtmosphere(tel; r0=0.2, L0=25.0)
    M2C, basis_hht = kl_modal_basis(KLHHtPSD(), dm, tel, atm; n_modes=2)
    @test size(M2C, 2) == 2
    @test size(basis_hht, 3) == 2
    basis2 = modal_basis(dm, tel; n_modes=2, method=KLHHtPSD(), atm=atm)
    @test size(basis2.M2C, 2) == 2
end

@testset "OPD maps and NCPA" begin
    tel = Telescope(resolution=8, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    dm = DeformableMirror(tel; n_act=2, influence_width=0.4)
    atm = KolmogorovAtmosphere(tel; r0=0.2, L0=25.0)
    map = OPDMap(fill(1.0, 8, 8))
    apply!(map, tel, DMReplace())
    @test sum(tel.state.opd) ≈ 64.0

    basis_default = AdaptiveOpticsSim.ncpa_basis(KLBasis(), tel, dm, atm; n_modes=2)
    basis_hht = AdaptiveOpticsSim.ncpa_basis(KLBasis(KLHHtPSD()), tel, dm, atm; n_modes=2)
    basis_dm = AdaptiveOpticsSim.ncpa_basis(KLBasis(KLDMModes()), tel, dm, atm; n_modes=2)
    @test basis_default ≈ basis_hht
    @test sum(abs.(basis_default .- basis_dm)) > 0

    coeffs = [1e-9, 2e-9]
    ncpa_default_kl = NCPA(tel, dm, atm; basis=KLBasis(), coefficients=coeffs)
    ncpa_hht = NCPA(tel, dm, atm; basis=KLBasis(KLHHtPSD()), coefficients=coeffs)
    ncpa_dm = NCPA(tel, dm, atm; basis=KLBasis(KLDMModes()), coefficients=coeffs)
    @test ncpa_default_kl.opd ≈ ncpa_hht.opd
    @test sum(abs.(ncpa_default_kl.opd .- ncpa_dm.opd)) > 0

    ncpa = NCPA(tel, dm, atm; basis=ZernikeModalBasis(), coefficients=[0.0, 1e-9, 2e-9])
    @test size(ncpa.opd) == (8, 8)
end

@testset "Spatial filter" begin
    tel = Telescope(resolution=8, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    src = Source(band=:I, magnitude=0.0)
    sf = SpatialFilter(tel; shape=CircularFilter(), diameter=4, zero_padding=2)
    phase, amp = filter!(sf, tel, src)
    @test size(phase) == (8, 8)
    @test size(amp) == (8, 8)
end

@testset "Gain sensing camera" begin
    mask = ones(8, 8)
    basis = rand(8, 8, 3)
    gsc = GainSensingCamera(mask, basis)
    frame = abs.(randn(8, 8))
    calibrate!(gsc, frame)
    og = compute_optical_gains!(gsc, frame)
    @test length(og) == 3
    @test length(weak_mode_mask(gsc)) == 3
    @test all(isfinite, og)
    @test detector_metadata(gsc) === nothing

    weak_gsc = GainSensingCamera(mask, zeros(8, 8, 2); sensitivity_floor=1e-6)
    calibrate!(weak_gsc, frame)
    weak_og = compute_optical_gains!(weak_gsc, frame)
    @test all(weak_mode_mask(weak_gsc))
    @test weak_og == ones(2)

    det = Detector(noise=NoiseReadout(1e-3), integration_time=2.0, qe=0.8, psf_sampling=2, binning=4)
    gsc_with_det = GainSensingCamera(mask, basis; detector=det)
    metadata = detector_metadata(gsc_with_det)
    @test metadata isa GSCDetectorMetadata
    @test metadata.integration_time == 2.0
    @test metadata.qe == 0.8
    @test metadata.psf_sampling == 2
    @test metadata.binning == 4
    @test metadata.noise == :readout
    @test metadata.readout_sigma == 1e-3
    @test occursin("psf_sampling=2", sprint(show, MIME"text/plain"(), gsc_with_det))

    detach_detector!(gsc_with_det)
    @test detector_metadata(gsc_with_det) === nothing
    attach_detector!(gsc_with_det, det)
    @test detector_metadata(gsc_with_det) isa GSCDetectorMetadata
end

@testset "LiFT" begin
    tel = Telescope(resolution=8, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    src = Source(band=:I, magnitude=0.0)
    det = Detector(noise=NoiseNone(), psf_sampling=1)
    basis = rand(8, 8, 3)
    diversity = zeros(8, 8)
    lift = LiFT(tel, src, basis, det; diversity_opd=diversity, iterations=2, img_resolution=8, numerical=true)
    H = lift_interaction_matrix(lift, zeros(3), [1, 2])
    @test size(H) == (64, 2)
    kernel = [0.0 1.0 0.0; 1.0 4.0 1.0; 0.0 1.0 0.0]
    lift_analytic = LiFT(tel, src, basis, det; diversity_opd=diversity, iterations=2,
        img_resolution=8, numerical=false, object_kernel=kernel)
    H_analytic = lift_interaction_matrix(lift_analytic, zeros(3), [1, 2])
    @test size(H_analytic) == (64, 2)
    psf = compute_psf!(tel, src; zero_padding=1)
    coeffs = reconstruct(lift, psf, [1, 2])
    @test length(coeffs) == 2
    @test AdaptiveOpticsSim.effective_solve_mode(AdaptiveOpticsSim.ScalarCPUStyle(), LiFTSolveAuto()) isa LiFTSolveQR
    diag = diagnostics(lift)
    @test diag.used_qr isa Bool
    @test isfinite(diag.residual_norm)
    @test isfinite(diag.weighted_residual_norm)
    @test isfinite(diag.update_norm)
    @test isfinite(diag.condition_ratio) || isinf(diag.condition_ratio)
    lift_normal = LiFT(tel, src, basis, det; diversity_opd=diversity, iterations=2,
        img_resolution=8, numerical=true, solve_mode=LiFTSolveNormalEquations())
    coeffs_normal = reconstruct(lift_normal, psf, [1, 2])
    @test length(coeffs_normal) == 2
    @test all(isfinite, coeffs_normal)
    @test !diagnostics(lift_normal).used_qr
    lift_damped = LiFT(tel, src, basis, det; diversity_opd=diversity, iterations=2,
        img_resolution=8, numerical=true, damping=LiFTLevenbergMarquardt())
    coeffs_damped = reconstruct(lift_damped, psf, [1, 2])
    @test length(coeffs_damped) == 2
    @test all(isfinite, coeffs_damped)
    @test diagnostics(lift_damped).regularization >= 0
    lift_adaptive = LiFT(tel, src, basis, det; diversity_opd=diversity, iterations=2,
        img_resolution=8, numerical=true, damping=LiFTAdaptiveLevenbergMarquardt())
    coeffs_adaptive = reconstruct(lift_adaptive, psf, [1, 2])
    @test length(coeffs_adaptive) == 2
    @test all(isfinite, coeffs_adaptive)
    @test diagnostics(lift_adaptive).regularization >= 0
    det_binned = Detector(noise=NoiseNone(), psf_sampling=1, binning=2)
    frame_binned = capture!(det_binned, psf; rng=MersenneTwister(3))
    lift_binned = LiFT(tel, src, basis, det_binned; diversity_opd=diversity, iterations=2,
        img_resolution=size(frame_binned, 1), numerical=true)
    coeffs_binned = reconstruct(lift_binned, frame_binned, [1, 2])
    @test length(coeffs_binned) == 2
    @test all(isfinite, coeffs_binned)
    det_readout = Detector(noise=NoiseReadout(1e-3), psf_sampling=1)
    lift_readout = LiFT(tel, src, basis, det_readout; diversity_opd=diversity, iterations=2,
        img_resolution=8, numerical=true)
    coeffs_readout = reconstruct(lift_readout, psf, [1, 2])
    @test length(coeffs_readout) == 2
    @test all(isfinite, coeffs_readout)

    sep_kernel = [1.0, 2.0, 1.0] * transpose([1.0, 0.5, 1.0])
    lift_sep = LiFT(tel, src, basis, det; diversity_opd=diversity, iterations=2,
        img_resolution=8, numerical=false, object_kernel=sep_kernel)
    @test lift_sep.params.object_kernel isa AdaptiveOpticsSim.LiFTSeparableObjectKernel
    dense_conv = similar(psf)
    sep_conv = similar(psf)
    tmp_conv = similar(psf)
    AdaptiveOpticsSim.conv2d_same!(dense_conv, psf, sep_kernel)
    AdaptiveOpticsSim.conv2d_same_separable!(
        sep_conv,
        tmp_conv,
        psf,
        lift_sep.params.object_kernel.row,
        lift_sep.params.object_kernel.col,
    )
    @test isapprox(sep_conv, dense_conv; rtol=1e-6, atol=1e-6)
end

@testset "Phase statistics" begin
    tel = Telescope(resolution=8, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    atm = KolmogorovAtmosphere(tel; r0=0.2, L0=25.0)
    rho = [0.0, 1e-6, 0.1, 1.0]
    cov = phase_covariance(rho, atm)
    @test length(cov) == length(rho)
    @test all(isfinite, cov)
    @test cov[1] >= cov[2]
    @test cov[2] <= cov[1]
    @test cov[3] <= cov[2]
    @test cov[4] <= cov[3]
    @test abs(cov[1] - cov[2]) / cov[1] < 1e-3
    var = phase_variance(atm)
    @test var > 0
    psd = phase_spectrum([0.1], atm)
    @test length(psd) == 1
    @test psd[1] > 0
    delta = tel.params.diameter / tel.params.resolution
    screen = ft_phase_screen(atm, 8, 0.1; rng=MersenneTwister(1))
    @test size(screen) == (8, 8)
    ensure_psd!(atm, delta)
    runtime_screen_rng = MersenneTwister(7)
    helper_screen_rng = MersenneTwister(7)
    advance!(atm, tel; rng=runtime_screen_rng)
    helper_screen, helper_psd = ft_phase_screen(atm, tel.params.resolution, delta; rng=helper_screen_rng, return_psd=true)
    @test helper_screen ≈ atm.state.opd
    @test helper_psd ≈ atm.state.psd

    for z in (1e-6, 1e-3, 0.1, 1.0, 4.0, 10.0, 40.0, 140.0)
        ref = SpecialFunctions.besselk(5 / 6, z)
        approx = AdaptiveOpticsSim._kv56_scalar(z)
        scaled_ref = z^(5 / 6) * ref
        scaled_approx = AdaptiveOpticsSim._scaled_kv56_scalar(z)
        @test isapprox(real(approx), ref; rtol=2e-4, atol=1e-10)
        @test isapprox(scaled_approx, scaled_ref; rtol=1e-7, atol=1e-12)
    end
end

@testset "Mis-registration identification" begin
    tel = Telescope(resolution=8, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    dm = DeformableMirror(tel; n_act=2, influence_width=0.4)
    wfs = ShackHartmann(tel; n_subap=2)
    basis = modal_basis(dm, tel; n_modes=2)
    meta = compute_meta_sensitivity_matrix(tel, dm, wfs, basis.M2C[:, 1:2]; n_mis_reg=2)
    est = estimate_misregistration(meta, meta.calib0.D; misregistration_zero=Misregistration())
    @test est.shift_x ≈ 0.0
    @test est.shift_y ≈ 0.0
end

@testset "Interface conformance" begin
    tel = Telescope(resolution=8, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    src = Source(band=:I, magnitude=0.0)
    lgs = LGSSource()
    atm = KolmogorovAtmosphere(tel; r0=0.2, L0=25.0)
    wfs = ShackHartmann(tel; n_subap=2)
    dm = DeformableMirror(tel; n_act=2, influence_width=0.4)
    det = Detector(noise=NoiseNone())
    psf = fill(1.0, 8, 8)
    opd_map = OPDMap(fill(0.1, size(tel.state.opd)))
    ncpa = NCPA(tel, dm, atm; coefficients=[0.01, -0.02])
    imat = interaction_matrix(dm, wfs, tel; amplitude=0.1)
    modal = ModalReconstructor(imat; gain=1.0)
    mapped = MappedReconstructor(Matrix{Float64}(I, length(dm.state.coefs), length(dm.state.coefs)), imat; gain=0.5)
    ctrl = DiscreteIntegratorController(length(wfs.state.slopes); gain=0.1, tau=0.02)
    sim = AOSimulation(tel, atm, src, dm, wfs)
    runtime = ClosedLoopRuntime(sim, modal; rng=MersenneTwister(9))
    wfs_diffractive = ShackHartmann(tel; n_subap=2, mode=Diffractive())
    zwfs = ZernikeWFS(tel; n_subap=2)
    curv_count = CurvatureWFS(tel; n_subap=2, readout_model=CurvatureCountingReadout())
    ast = Asterism([src, Source(band=:I, magnitude=1.0, coordinates=(1.0, -45.0))])
    moving_atm = MultiLayerAtmosphere(tel; r0=0.2, L0=25.0, fractional_cn2=[1.0],
        wind_speed=[0.0], wind_direction=[0.0], altitude=[0.0])
    infinite_atm = InfiniteMultiLayerAtmosphere(tel; r0=0.2, L0=25.0, fractional_cn2=[1.0],
        wind_speed=[0.0], wind_direction=[0.0], altitude=[0.0], screen_resolution=33, stencil_size=35)
    @test CCDSensor <: FrameSensorType
    @test CMOSSensor <: FrameSensorType
    @test AvalancheFrameSensorType <: FrameSensorType
    @test HgCdTeAvalancheArraySensorType <: AvalancheFrameSensorType
    @test EMCCDSensor <: AvalancheFrameSensorType
    @test InGaAsSensor <: FrameSensorType
    @test HgCdTeAvalancheArraySensor <: HgCdTeAvalancheArraySensorType
    @test !supports_avalanche_gain(CCDSensor())
    @test !supports_sensor_glow(CMOSSensor())
    @test supports_detector_defect_maps(CMOSSensor())
    @test supports_detector_defect_maps(InGaAsSensor())
    @test supports_shutter_timing(CMOSSensor())
    @test !supports_shutter_timing(CCDSensor())
    @test !supports_detector_persistence(CMOSSensor())
    @test supports_detector_persistence(InGaAsSensor(persistence_model=ExponentialPersistence(0.1, 0.9)))
    @test !supports_detector_nonlinearity(CMOSSensor())
    @test supports_detector_nonlinearity(InGaAsSensor())
    @test !supports_nondestructive_reads(CCDSensor())
    @test supports_nondestructive_reads(HgCdTeAvalancheArraySensor())
    @test !supports_reference_read_subtraction(EMCCDSensor())
    @test supports_reference_read_subtraction(HgCdTeAvalancheArraySensor())
    @test !supports_readout_correction(EMCCDSensor())
    @test supports_readout_correction(HgCdTeAvalancheArraySensor())
    @test supports_read_cube(HgCdTeAvalancheArraySensor())
    @test AdaptiveOpticsSim.readout_correction_symbol(ReferenceRowCommonModeCorrection()) == :reference_row_common_mode
    @test AdaptiveOpticsSim.readout_correction_symbol(ReferenceColumnCommonModeCorrection()) == :reference_column_common_mode
    @test AdaptiveOpticsSim.readout_correction_symbol(ReferenceOutputCommonModeCorrection(4)) == :reference_output_common_mode
    @test APDSensor <: CountingSensorType
    @test curv_count.params.readout_model isa CurvatureCountingReadout

    # IF-SRC
    assert_source_interface(src)
    assert_source_interface(lgs)
    # IF-ATM
    assert_atmosphere_interface(atm, tel)
    @test applicable(propagate!, moving_atm, tel, src)
    @test applicable(propagate!, infinite_atm, tel, src)
    assert_atmosphere_layer_interface(moving_atm.layers[1], tel, MersenneTwister(11), src)
    assert_atmosphere_layer_interface(infinite_atm.layers[1], tel, MersenneTwister(12), src)
    # IF-WFS
    assert_wfs_interface(wfs, tel)
    # IF-DM
    assert_dm_interface(dm, tel)
    # IF-DET
    assert_detector_interface(det, psf)
    # IF-OPT
    assert_optical_element_interface(opd_map, tel)
    assert_optical_element_interface(ncpa, tel)
    # IF-REC
    assert_reconstructor_interface(modal, wfs.state.slopes, length(dm.state.coefs))
    assert_reconstructor_interface(mapped, wfs.state.slopes, length(dm.state.coefs))
    # IF-CTRL
    assert_controller_interface(ctrl, wfs.state.slopes, 0.01)
    # IF-SIM
    iface = assert_control_simulation_interface(runtime)
    @test !supports_prepared_runtime(runtime)
    @test !supports_detector_output(runtime)
    @test !supports_stacked_sources(runtime)
    @test !supports_grouped_execution(runtime)
    @test simulation_interface(iface) === iface
    @test !supports_prepared_runtime(wfs, src)
    @test supports_prepared_runtime(wfs_diffractive, src)
    @test supports_prepared_runtime(wfs_diffractive, ast)
    @test supports_prepared_runtime(zwfs, src)
    @test supports_prepared_runtime(CurvatureWFS(tel; n_subap=2), src)
    @test supports_prepared_runtime(PyramidWFS(tel; n_subap=2, mode=Diffractive()), src)
    @test !supports_stacked_sources(wfs, src)
    @test supports_stacked_sources(wfs, ast)
    @test supports_stacked_sources(wfs_diffractive, ast)
    prepare_runtime_wfs!(wfs_diffractive, tel, src)
    @test wfs_diffractive.state.calibrated
    prepare_runtime_wfs!(zwfs, tel, src)
    @test zwfs.state.calibrated
    curv = CurvatureWFS(tel; n_subap=2)
    prepare_runtime_wfs!(curv, tel, src)
    @test curv.state.calibrated
end

@testset "Calibration workflow contracts" begin
    tel = Telescope(resolution=8, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    src = Source(band=:I, magnitude=0.0)
    atm = KolmogorovAtmosphere(tel; r0=0.2, L0=25.0)
    dm = DeformableMirror(tel; n_act=2, influence_width=0.4)
    wfs = ShackHartmann(tel; n_subap=2)
    det = Detector(noise=NoiseNone(), integration_time=1.0, qe=1.0, binning=1)

    basis = modal_basis(dm, tel; n_modes=2)
    assert_modal_basis_contract(basis, length(dm.state.coefs), 2)

    imat = interaction_matrix(dm, wfs, tel; amplitude=0.1)
    assert_interaction_matrix_contract(imat, length(wfs.state.slopes), length(dm.state.coefs), 0.1)

    imat_basis = interaction_matrix(dm, wfs, tel, basis.M2C; amplitude=0.1)
    assert_interaction_matrix_contract(imat_basis, length(wfs.state.slopes), size(basis.M2C, 2), 0.1)

    vault = CalibrationVault(imat.matrix)
    assert_calibration_vault_contract(vault, imat.matrix)
    vault_noinv = CalibrationVault(imat.matrix; invert=false)
    assert_calibration_vault_contract(vault_noinv, imat.matrix; inverted=false)
    vault_trunc = with_truncation(vault, 0)
    assert_calibration_vault_contract(vault_trunc, imat.matrix)
    @test vault_trunc.n_trunc == 0

    calib = ao_calibration(tel, dm, wfs; n_modes=2, amplitude=0.1, basis=basis)
    assert_ao_calibration_contract(calib, length(dm.state.coefs), 2)
    @test calib.calibration.D == imat_basis.matrix

    meta = compute_meta_sensitivity_matrix(tel, dm, wfs, basis.M2C[:, 1:2]; n_mis_reg=2)
    assert_meta_sensitivity_contract(meta, 2)

    sprint = SPRINT(tel, dm, wfs, basis.M2C[:, 1:2]; n_mis_reg=2)
    @test sprint.meta isa MetaSensitivity
    est = estimate!(sprint, meta.calib0.D)
    @test est isa Misregistration

    diversity = fill(eltype(tel.state.opd)(1e-9), size(tel.state.opd))
    lift_basis = basis_from_m2c(dm, tel, basis.M2C)
    lift = LiFT(tel, src, lift_basis, det; diversity_opd=diversity, iterations=2, img_resolution=8, numerical=true)
    psf_in = compute_psf!(tel, src; zero_padding=1)
    coeffs = zeros(eltype(psf_in), 2)
    reconstruct!(coeffs, lift, psf_in, [1, 2]; check_convergence=false)
    @test length(coeffs) == 2
    @test diagnostics(lift).residual_norm >= 0
end

@testset "Telemetry and config" begin
    telemetry = Telemetry()
    record!(telemetry, 1, 0.0; wfe_rms=1.0, strehl=0.5, loop_gain=0.2)
    @test length(telemetry) == 1
    row = telemetry[1]
    @test row.iter == 1
    @test row.strehl ≈ 0.5

    @test Tables.istable(Telemetry)
    @test Tables.rowaccess(Telemetry)
    rows = collect(Tables.rows(telemetry))
    @test length(rows) == 1
    @test rows[1].loop_gain ≈ 0.2

    closed_loop = ClosedLoopTrace(Float32[
        100 80 0.50 3 4
        120 90 0.45 5 6
    ]; dt=0.002f0, t0=0.1f0)
    @test length(closed_loop) == 2
    @test eltype(typeof(closed_loop)) == ClosedLoopTraceRow{Float32}
    @test closed_loop[2].iter == 2
    @test closed_loop[2].t ≈ 0.102f0
    @test closed_loop[2].residual_rms_nm ≈ 90.0f0
    @test Tables.istable(typeof(closed_loop))
    @test Tables.rowaccess(typeof(closed_loop))
    closed_rows = collect(Tables.rows(closed_loop))
    @test length(closed_rows) == 2
    @test closed_rows[1].command_norm ≈ 4.0f0

    gsc_closed_loop = GSCClosedLoopTrace(Float32[
        100 80 0.50 3 0.9 4
        120 90 0.45 5 0.8 6
    ]; dt=0.002f0, t0=0.1f0)
    @test length(gsc_closed_loop) == 2
    @test eltype(typeof(gsc_closed_loop)) == GSCClosedLoopTraceRow{Float32}
    @test gsc_closed_loop[2].mean_optical_gain ≈ 0.8f0
    @test Tables.istable(typeof(gsc_closed_loop))
    @test Tables.rowaccess(typeof(gsc_closed_loop))
    gsc_closed_rows = collect(Tables.rows(gsc_closed_loop))
    @test length(gsc_closed_rows) == 2
    @test gsc_closed_rows[1].slope_norm ≈ 3.0f0

    replay = GSCAtmosphereReplayTrace(Float32[
        140 100 110 0.50 0.45 3 0.9
        150 105 115 0.48 0.42 5 0.8
    ]; dt=0.002f0, t0=0.1f0)
    @test length(replay) == 2
    @test eltype(typeof(replay)) == GSCAtmosphereReplayTraceRow{Float32}
    @test replay[2].sci_strehl ≈ 0.42f0
    @test replay[2].slope_norm ≈ 5.0f0
    @test Tables.istable(typeof(replay))
    @test Tables.rowaccess(typeof(replay))
    replay_rows = collect(Tables.rows(replay))
    @test length(replay_rows) == 2
    @test replay_rows[1].ngs_forcing_rms_nm ≈ 140.0f0

    @test_throws DimensionMismatchError ClosedLoopTrace(zeros(2, 4))
    @test_throws DimensionMismatchError GSCClosedLoopTrace(zeros(2, 5))
    @test_throws DimensionMismatchError GSCAtmosphereReplayTrace(zeros(2, 6))

    tel = Telescope(resolution=8, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    src = Source(band=:I, magnitude=0.0)
    cfg = snapshot_config(tel=tel, src=src)
    @test haskey(cfg, "tel")
    @test haskey(cfg, "src")
    path, io = mktemp()
    close(io)
    write_config_toml(path, cfg)
    parsed = TOML.parsefile(path)
    @test parsed["tel"]["resolution"] == tel.params.resolution
end

@testset "Reference compare conventions" begin
    convention = parse_reference_compare_convention(Dict(
        "convention" => "oopao_geometric_sh_signal_2d",
    ))
    raw = Float64[1, 2, 3, 4]
    adapted = adapt_compare_convention(convention, raw)
    @test adapted ≈ OOPAO_GEOMETRIC_SH_SLOPE_SCALE .* Float64[3, 4, 1, 2]

    legacy = parse_reference_compare_convention(Dict(
        "swap_halves" => true,
        "scale" => OOPAO_GEOMETRIC_SH_SLOPE_SCALE,
    ))
    @test adapt_compare_convention(legacy, raw) == adapted
    @test parse_reference_compare_convention(nothing) isa IdentityCompareConvention
    @test adapt_compare_convention(IdentityCompareConvention(), raw) === raw
    @test_throws InvalidConfiguration parse_reference_compare_convention(Dict(
        "swap_halves" => true,
    ))
end

@testset "Reference storage and detector conventions" begin
    col_data = Float64[1, 2, 3, 4]
    row_data = Float64[1, 2, 3, 4, 5, 6]
    @test reshape_reference_data(col_data, (2, 2), JuliaColumnMajorStorage()) == [1 3; 2 4]
    @test reshape_reference_data(row_data, (2, 3), NumPyRowMajorStorage()) == [1 2 3; 4 5 6]
    @test parse_reference_storage_convention("F") isa JuliaColumnMajorStorage
    @test parse_reference_storage_convention("numpy_row_major") isa NumPyRowMajorStorage
    @test_throws InvalidConfiguration parse_reference_storage_convention("weird")

    tel = Telescope(resolution=8, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    det = Detector(noise=NoiseNone(), psf_sampling=2, binning=1)
    @test reference_lift_img_resolution(tel, det, Dict{String,Any}()) == 16
    @test reference_lift_img_resolution(tel, det, Dict{String,Any}("img_resolution" => 12)) == 12
end

@testset "Reference harness fixture" begin
    root = mktempdir()
    create_reference_fixture(root)
    bundle = load_reference_bundle(root)
    @test length(bundle.cases) == 6
    for case in bundle.cases
        result = validate_reference_case(case)
        @test size(result.actual) == size(result.expected)
        @test result.ok
    end
end

@testset "OOPAO reference regression" begin
    root = default_reference_root()
    if has_reference_bundle(root)
        bundle = load_reference_bundle(root)
        @test !isempty(bundle.cases)
        for case in bundle.cases
            result = validate_reference_case(case)
            @test size(result.actual) == size(result.expected)
            @test result.ok
        end
    else
        @info "Skipping OOPAO reference regression; no manifest found" root=root
        @test true
    end
end

@testset "Tutorial examples" begin
    image = run_tutorial_example("image_formation.jl")
    @test size(image.psf_nominal) == size(image.psf_aberrated)
    @test maximum(image.psf_nominal) > maximum(image.psf_aberrated)

    detector = run_tutorial_example("detector.jl")
    @test sum(detector.frame_native) ≈ sum(detector.psf)
    @test size(detector.frame_sampled, 1) < size(detector.frame_native, 1)

    asterism = run_tutorial_example("asterism.jl")
    @test size(asterism.per_source_psf, 3) == 4
    @test sum(asterism.combined_psf) >= sum(@view asterism.per_source_psf[:, :, 1])

    extended = run_tutorial_example("extended_source_sensing.jl")
    @test extended.n_samples == 25
    @test norm(extended.sh_spot_delta) > 1e-12
    @test norm(extended.pyramid_intensity_delta) > 1e-12

    sh_subaps = run_tutorial_example("shack_hartmann_subapertures.jl")
    @test sh_subaps.n_valid > 0
    @test sh_subaps.calibrated
    @test sh_subaps.metadata.n_valid_subap == sh_subaps.n_valid
    @test all(isfinite, sh_subaps.slopes)

    spatial = run_tutorial_example("spatial_filter.jl")
    @test size(spatial.filtered_phase) == size(spatial.filtered_amplitude)
    @test all(isfinite, spatial.filtered_phase)

    ncpa = run_tutorial_example("ncpa.jl")
    @test size(ncpa.psf, 1) == 2 * size(ncpa.ncpa_opd, 1)
    @test size(ncpa.psf, 2) == 2 * size(ncpa.ncpa_opd, 2)
    @test maximum(ncpa.psf) > 0

    lift = run_tutorial_example("lift.jl")
    @test length(lift.coeffs_true) == length(lift.coeffs_fit)
    @test all(isfinite, lift.coeffs_fit)

    sprint = run_tutorial_example("sprint.jl")
    @test isfinite(sprint.estimate.shift_x)
    @test isfinite(sprint.estimate.shift_y)

    gsc = run_tutorial_example("gain_sensing_camera.jl")
    @test length(gsc.optical_gains) == 4
    @test all(isfinite, gsc.optical_gains)
    @test sum(gsc.calibration_frame) > 0
    @test sum(gsc.frame) > 0
    @test size(gsc.atmosphere_trace, 2) == 7
    @test all(isfinite, gsc.atmosphere_trace)

    transfer = run_tutorial_example("transfer_function.jl")
    @test size(transfer.rejection_db, 2) == length(transfer.loop_gains)
    @test all(isfinite, transfer.rejection_db)
    @test all(isfinite, transfer.closed_loop_db)

    tomography = run_tutorial_example("tomography.jl")
    @test all(isfinite, tomography.wavefront[.!isnan.(tomography.wavefront)])
    @test length(tomography.commands) == 4
    @test all(isfinite, tomography.commands)

    for name in ("closed_loop_shack_hartmann.jl", "closed_loop_pyramid.jl", "closed_loop_bioedge.jl", "closed_loop_zernike.jl")
        loop = run_tutorial_example(name)
        @test length(loop.residual_before) == length(loop.residual_after)
        @test all(isfinite, loop.residual_before)
        @test all(isfinite, loop.residual_after)
        if hasproperty(loop, :final_psf)
            @test maximum(loop.final_psf) > 0
        else
            @test maximum(loop.final_frame) > 0
            @test all(isfinite, loop.final_slopes)
            @test all(isfinite, loop.final_command)
        end
    end
end
