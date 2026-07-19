using Test
using Aqua
using AdaptiveOpticsSim
using LinearAlgebra
using Random
using SpecialFunctions
using Statistics
using Tables
using TOML

BLAS.set_num_threads(1)
AdaptiveOpticsSim.set_fft_provider_threads!(1)

const TEST_ATMOSPHERE_STEP = 1e-3

# `Pkg.test(coverage=true)` normally sets Julia's coverage option in the test
# process. The explicit workflow flag keeps allocation gating deterministic if
# that propagation changes while ordinary CPU CI continues to enforce it.
@inline coverage_runner_flag(value::AbstractString) =
    lowercase(strip(value)) in ("1", "true", "yes", "on")

@inline coverage_instrumented() =
    coverage_runner_flag(get(ENV, "ADAPTIVEOPTICS_TEST_COVERAGE", "false")) ||
    Base.JLOptions().code_coverage != 0

# The package exports only the user-facing API. The tests intentionally exercise
# internal extension seams too, so make those names available without expanding
# the public export list.
for name in names(AdaptiveOpticsSim; all=true)
    s = String(name)
    if Base.isidentifier(s) && !startswith(s, "#") && !isdefined(@__MODULE__, name)
        @eval const $(name) = getfield(AdaptiveOpticsSim, $(QuoteNode(name)))
    end
end

include(joinpath(dirname(@__DIR__), "examples", "support", "subaru_ao188_simulation.jl"))
include(joinpath(dirname(@__DIR__), "examples", "support", "subaru_ao3k_simulation.jl"))
using .SubaruAO188Simulation
using .SubaruAO3kSimulation

include("reference_harness.jl")
include("gate0_characterization_harness.jl")
include("ka_cpu_matrix.jl")
include("tomography.jl")

_rng_family(::Xoshiro) = :xoshiro
_rng_family(::MersenneTwister) = :mersenne_twister
_rng_family(::AbstractRNG) = :other

@testset "RNG policy helpers" begin
    @test _rng_family(runtime_rng(1)) == :xoshiro
    @test _rng_family(deterministic_reference_rng(1)) == :mersenne_twister
    @test rand(runtime_rng(42), UInt64) == rand(runtime_rng(42), UInt64)
    @test rand(deterministic_reference_rng(42), UInt64) == rand(deterministic_reference_rng(42), UInt64)
    @test coverage_runner_flag("true")
    @test coverage_runner_flag(" YES ")
    @test !coverage_runner_flag("false")
    @test coverage_instrumented() ==
        (coverage_runner_flag(get(ENV, "ADAPTIVEOPTICS_TEST_COVERAGE",
            "false")) || Base.JLOptions().code_coverage != 0)
end

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
    @test hasmethod(photon_irradiance, Tuple{typeof(src)})
end

function assert_atmosphere_interface(atm, tel)
    pupil = PupilFunction(tel)
    test_time = one(eltype(pupil.opd))
    @test applicable(advance_by!, atm, test_time)
    @test applicable(advance_to!, atm, test_time)
    @test applicable(prepare_atmosphere_renderer, atm, tel)
    renderer = prepare_atmosphere_renderer(atm, tel)
    epoch = advance_by!(atm, zero(test_time); rng=MersenneTwister(0))
    @test applicable(render_atmosphere!, pupil, renderer, atm,
        epoch)
end

function assert_atmosphere_layer_interface(layer, tel, rng, src)
    pupil = PupilFunction(tel)
    @test layer isa AdaptiveOpticsSim.AbstractAtmosphereLayer
    @test applicable(AdaptiveOpticsSim.render_layer!, similar(pupil.opd),
        layer, zero(eltype(pupil.opd)), zero(eltype(pupil.opd)),
        one(eltype(pupil.opd)))
    @test applicable(AdaptiveOpticsSim.render_layer_accumulate!,
        similar(pupil.opd), layer, zero(eltype(pupil.opd)),
        zero(eltype(pupil.opd)), one(eltype(pupil.opd)))
    altitude = AdaptiveOpticsSim.layer_altitude(layer)
    shift_x, shift_y, footprint_scale =
        AdaptiveOpticsSim.layer_source_geometry(src, altitude, tel,
            eltype(pupil.opd))
    sample = similar(pupil.opd)
    fill!(sample, zero(eltype(sample)))
    AdaptiveOpticsSim.render_layer!(sample, layer, shift_x, shift_y, footprint_scale)
    @test size(sample) == size(pupil.opd)
    fill!(sample, zero(eltype(sample)))
    AdaptiveOpticsSim.render_layer_accumulate!(sample, layer, shift_x, shift_y, footprint_scale)
    @test size(sample) == size(pupil.opd)
end

function assert_wfs_interface(wfs, tel)
    pupil = PupilFunction(tel)
    @test applicable(update_valid_mask!, wfs, pupil)
    @test applicable(measure!, wfs, pupil)
    @test slopes(wfs) isa AbstractVector
    @test supports_valid_subaperture_mask(wfs) == !isnothing(valid_subaperture_mask(wfs))
    @test supports_reference_signal(wfs) == !isnothing(reference_signal(wfs))
    @test supports_camera_frame(wfs) == !isnothing(camera_frame(wfs))
end

function assert_detector_interface(det, psf)
    @test applicable(capture!, det, psf)
end

function assert_dm_interface(dm, tel)
    pupil = PupilFunction(tel)
    @test applicable(build_influence_functions!, dm, tel)
    @test applicable(update_surface!, dm)
    @test applicable(apply_surface!, pupil, dm, DMAdditive())
    @test topology_command_count(topology(dm)) == length(dm.state.coefs)
    @test size(actuator_coordinates(dm), 2) == length(dm.state.coefs)
    @test length(valid_actuator_mask(dm)) >= length(dm.state.coefs)
end

function assert_optical_element_interface(element, tel)
    pupil = PupilFunction(tel)
    @test element isa AbstractOpticalElement
    @test applicable(surface_opd, element)
    @test applicable(apply_surface!, pupil, element, DMAdditive())
    @test applicable(apply_surface!, pupil, element, DMReplace())
end

function assert_reconstructor_interface(recon, slopes, expected_length::Int)
    @test recon isa AbstractReconstructorOperator
    out = zeros(eltype(slopes), expected_length)
    @test applicable(reconstruct!, out, recon, slopes)
    reconstruct!(out, recon, slopes)
    @test length(out) == expected_length
    policy = inverse_policy(recon)
    spectrum = singular_values(recon)
    cond = condition_number(recon)
    rank = effective_rank(recon)
    @test policy isa InversePolicy
    @test spectrum isa AbstractVector
    @test cond isa Real
    @test rank isa Integer
    @test 0 <= rank <= length(spectrum)
    allocated = reconstruct(recon, slopes)
    @test length(allocated) == expected_length
end

function assert_controller_interface(ctrl, input, dt::Real)
    @test applicable(update!, ctrl, input, dt)
    output = update!(ctrl, input, dt)
    @test length(output) == length(input)
    @test controller_output(ctrl) === output
    if supports_controller_reset(ctrl)
        reset_controller!(ctrl)
        @test all(iszero, controller_output(ctrl))
    end
end

function assert_control_simulation_interface(sim)
    @test sim isa AbstractControlSimulation
    @test applicable(step!, sim)
    @test applicable(AdaptiveOpticsSim.simulation_interface, sim)
    @test applicable(readout, sim)
    iface = AdaptiveOpticsSim.simulation_interface(sim)
    sim_readout = readout(sim)
    iface_readout = readout(iface)
    @test command(sim_readout) === command(sim)
    @test slopes(sim_readout) === slopes(sim)
    @test command(iface_readout) === command(iface)
    return iface
end

function assert_interaction_matrix_contract(imat, expected_rows::Int, expected_cols::Int, amplitude::Real)
    @test imat isa InteractionMatrix
    @test size(imat.matrix) == (expected_rows, expected_cols)
    @test forward_operator(imat) === imat.matrix
    @test calibration_amplitude(imat) ≈ amplitude
end

function assert_control_matrix_contract(control_matrix, forward::AbstractMatrix; inverted::Bool=true)
    @test control_matrix isa ControlMatrix
    @test control_matrix.D === forward
    @test forward_operator(control_matrix) === forward
    @test inverse_policy(control_matrix) === control_matrix.policy
    @test singular_values(control_matrix) === control_matrix.singular_values
    @test isequal(condition_number(control_matrix), control_matrix.cond)
    @test effective_rank(control_matrix) == control_matrix.effective_rank
    @test truncation_count(control_matrix) == control_matrix.n_trunc
    if inverted
        @test inverse_operator_matrix(control_matrix) === control_matrix.M
        @test !isnothing(control_matrix.M)
        @test size(control_matrix.M, 2) == size(forward, 1)
        @test length(singular_values(control_matrix)) == min(size(forward)...)
        @test effective_rank(control_matrix) >= 0
    else
        @test isnothing(inverse_operator_matrix(control_matrix))
        @test isnothing(control_matrix.M)
        @test isempty(singular_values(control_matrix))
    end
end

function assert_modal_basis_contract(basis::ModalBasis, n_commands::Int, n_modes::Int)
    @test modal_to_command(basis) === basis.M2C
    @test sampled_basis(basis) === basis.basis
    @test modal_projector(basis) === basis.projector
    @test size(basis.M2C) == (n_commands, n_modes)
    @test size(basis.basis, 2) == n_modes
    if !isnothing(basis.projector)
        @test size(basis.projector, 2) == size(basis.basis, 1)
        @test size(basis.projector, 1) == n_modes
    end
end

function assert_ao_calibration_contract(calib::AOCalibration, n_commands::Int, n_modes::Int)
    @test modal_to_command(calib) === calib.M2C
    @test sampled_basis(calib) === calib.basis
    @test modal_projector(calib) === calib.projector
    @test control_matrix(calib) === calib.calibration
    @test size(calib.M2C) == (n_commands, n_modes)
    @test size(calib.basis, 2) == n_modes
    @test calib.calibration isa ControlMatrix
end

function assert_meta_sensitivity_contract(meta::AdaptiveOpticsSim.MetaSensitivity, n_fields::Int)
    @test meta.calib0 isa ControlMatrix
    @test meta.meta isa ControlMatrix
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
    atmosphere_step::Real=1e-3,
    r0::Real=0.2,
    L0::Real=25.0,
    fractional_cn2::AbstractVector=[1.0],
    wind_speed::AbstractVector=[0.0],
    wind_direction::AbstractVector=[0.0],
    altitude::AbstractVector=[0.0],
)
    tel = Telescope(resolution=resolution, diameter=diameter, central_obstruction=0.0)
    atm = MultiLayerAtmosphere(tel;
        r0=r0,
        L0=L0,
        fractional_cn2=fractional_cn2,
        wind_speed=wind_speed,
        wind_direction=wind_direction,
        altitude=altitude,
    )
    rng = MersenneTwister(seed)
    pupil = PupilFunction(tel)
    renderer = prepare_atmosphere_renderer(atm, tel)
    trace = Matrix{Float64}[]
    for _ in 1:steps
        advance_by!(atm, atmosphere_step; rng=rng)
        render_atmosphere!(pupil, renderer, atm, current_epoch(atm))
        push!(trace, copy(pupil.opd))
    end
    return trace
end

function moving_wfs_slope_trace(;
    seed::Integer=1,
    steps::Integer=4,
)
    tel = Telescope(resolution=16, diameter=8.0, central_obstruction=0.0)
    src = Source(band=:I, magnitude=0.0)
    atmosphere_step = 1e-3
    delta = tel.params.diameter / tel.params.resolution
    atm = MultiLayerAtmosphere(tel;
        r0=0.2,
        L0=25.0,
        fractional_cn2=[0.7, 0.3],
        wind_speed=[delta / atmosphere_step, 0.5 * delta / atmosphere_step],
        wind_direction=[0.0, 90.0],
        altitude=[0.0, 5000.0],
    )
    wfs = ShackHartmannWFS(tel; n_lenslets=4)
    pupil = PupilFunction(tel)
    renderer = prepare_atmosphere_renderer(atm, tel, src)
    rng = MersenneTwister(seed)
    trace = Vector{Vector{Float64}}(undef, steps)
    for i in 1:steps
        advance_by!(atm, atmosphere_step; rng=rng)
        render_atmosphere!(pupil, renderer, atm, current_epoch(atm))
        measure!(wfs, pupil, src)
        trace[i] = copy(slopes(wfs))
    end
    return trace
end

function moving_closed_loop_trace(;
    seed::Integer=1,
    steps::Integer=5,
)
    tel = Telescope(resolution=16, diameter=8.0, central_obstruction=0.0)
    src = Source(band=:I, magnitude=0.0)
    atmosphere_step = 1e-3
    delta = tel.params.diameter / tel.params.resolution
    atm = MultiLayerAtmosphere(tel;
        r0=0.2,
        L0=25.0,
        fractional_cn2=[0.7, 0.3],
        wind_speed=[delta / atmosphere_step, 0.5 * delta / atmosphere_step],
        wind_direction=[0.0, 90.0],
        altitude=[0.0, 5000.0],
    )
    dm = DeformableMirror(tel; n_act=4, influence_width=0.3)
    wfs = ZernikeWFS(tel; pupil_samples=4, diffraction_padding=2)
    sim = AOSimulation(tel, src, atm, dm, wfs)
    calibration_pupil = PupilFunction(tel)
    imat = interaction_matrix(dm, wfs, calibration_pupil, src;
        amplitude=1e-8)
    recon = ModalReconstructor(imat; gain=0.2)
    wfs_det = Detector(noise=NoiseNone(), integration_time=1.0, qe=1.0, binning=1)
    science_det = Detector(noise=NoiseNone(), integration_time=1.0, qe=1.0, binning=1)
    runtime = ClosedLoopRuntime(sim, recon; atmosphere_step=1e-3, rng=MersenneTwister(seed), wfs_detector=wfs_det, science_detector=science_det)
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
        slope_norms[i] = norm(slopes(runtime))
        command_norms[i] = norm(command(runtime))
        wfs_energy[i] = sum(abs, wfs_frame(runtime))
        science_energy[i] = sum(abs, science_frame(runtime))
    end

    return (; slope_norms, command_norms, wfs_energy, science_energy)
end

@test AdaptiveOpticsSim.PROJECT_STATUS == :in_development

@testset "Aqua" begin
    Aqua.test_all(
        AdaptiveOpticsSim;
        undocumented_names=false,
    )
end
