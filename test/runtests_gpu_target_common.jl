using Test
using AdaptiveOpticsSim
using AdaptiveOpticsSim.Optics
import AdaptiveOpticsSim.Optics: filter!
using AdaptiveOpticsSim.Backends
using AdaptiveOpticsSim.Detectors
using AdaptiveOpticsSim.Atmospheres
using AdaptiveOpticsSim.WavefrontSensors
using AdaptiveOpticsSim.Calibration
using AdaptiveOpticsSim.Control
using AdaptiveOpticsSim.Tomography
using AdaptiveOpticsSim.Ensembles
using AdaptiveOpticsSim.AlgorithmGraphs
using FixedSizeArrays: FixedSizeVector
using LinearAlgebra
using Random
using Statistics

# Hardware targets exercise qualified-public and internal backend contracts in
# addition to the routine exported workflow.
for name in names(Backends; all=true)
    s = String(name)
    if Base.isidentifier(s) && !startswith(s, "#") &&
            !isdefined(@__MODULE__, name)
        @eval const $(name) = getfield(Backends, $(QuoteNode(name)))
    end
end

for name in names(Detectors; all=true)
    s = String(name)
    if Base.isidentifier(s) && !startswith(s, "#") &&
            !isdefined(@__MODULE__, name)
        @eval const $(name) = getfield(Detectors, $(QuoteNode(name)))
    end
end

for name in names(Atmospheres; all=true)
    s = String(name)
    if Base.isidentifier(s) && !startswith(s, "#") &&
            !isdefined(@__MODULE__, name)
        @eval const $(name) = getfield(Atmospheres, $(QuoteNode(name)))
    end
end

for name in names(WavefrontSensors; all=true)
    s = String(name)
    if Base.isidentifier(s) && !startswith(s, "#") &&
            !isdefined(@__MODULE__, name)
        @eval const $(name) =
            getfield(WavefrontSensors, $(QuoteNode(name)))
    end
end

include("backend_optional_common.jl")
include(normpath(joinpath(@__DIR__, "..", "benchmarks", "support", "revolt_like_hil_common.jl")))
include(normpath(joinpath(@__DIR__, "..", "scripts", "gpu_builder_contract.jl")))

BLAS.set_num_threads(1)
Backends.set_fft_provider_threads!(1)

backend_target_branch_mode(::Type{Backends.CUDABackendTag}) =
    BackendStreamExecution()
backend_target_branch_mode(::Type{Backends.AMDGPUBackendTag}) =
    SequentialExecution()

graph_rng_device_resident(::Type{Backends.CUDABackendTag}) = true
graph_rng_device_resident(::Type{Backends.AMDGPUBackendTag}) = true

function run_graph_rng_capture_replay(
    ::Type{B},
) where {B<:Backends.GPUBackendTag}
    BackendArray = Backends.gpu_backend_array_type(B)
    output = BackendArray(zeros(Float32, 256))
    poisson_output = BackendArray(fill(8.0f0, 256))
    target = compute_device(output)
    context = Backends._prepare_device_execution_context(target)
    rng = Backends._prepare_graph_rng(target, UInt64(0x5eed))
    style = Backends.execution_style(output)

    Backends._with_prepared_device_execution_context(context) do
        # Compile every stochastic kernel before native command-graph capture
        # begins. Both kernels consume the same device-resident draw sequence.
        Backends.randn_backend_async!(style, rng, output)
        Backends.poisson_noise_async!(style, rng, poisson_output)
        Backends._synchronize_prepared_device_execution_context!(context)
        Backends._reset_graph_rng!(rng, UInt64(0x5eed))
        Backends._synchronize_prepared_device_execution_context!(context)

        captured = Backends._capture_prepared_device_graph(context) do
            Backends.randn_backend_async!(style, rng, output)
            Backends.poisson_noise_async!(style, rng, poisson_output)
        end

        fill!(poisson_output, 8.0f0)
        Backends._launch_prepared_device_graph!(captured, context)
        Backends._synchronize_prepared_device_execution_context!(context)
        first_draw = Array(output)
        first_poisson_draw = Array(poisson_output)
        @test all(x -> x >= 0 && isinteger(x), first_poisson_draw)

        fill!(poisson_output, 8.0f0)
        Backends._launch_prepared_device_graph!(captured, context)
        Backends._synchronize_prepared_device_execution_context!(context)
        second_draw = Array(output)
        second_poisson_draw = Array(poisson_output)
        @test second_draw != first_draw
        @test second_poisson_draw != first_poisson_draw

        Backends._reset_graph_rng!(rng, UInt64(0x5eed))
        Backends._synchronize_prepared_device_execution_context!(context)
        fill!(poisson_output, 8.0f0)
        Backends._launch_prepared_device_graph!(captured, context)
        Backends._synchronize_prepared_device_execution_context!(context)
        @test Array(output) == first_draw
        @test Array(poisson_output) == first_poisson_draw
    end
    return nothing
end

function require_backend_target!(::Type{B}) where {B<:Backends.GPUBackendTag}
    pkg = backend_package_name(B)
    pkg_path = Base.find_package(pkg)
    pkg_path === nothing && error("$(backend_label(B)) target requires $(pkg).jl in the active environment")
    import_backend_package!(B)
    backend_functional(B) || error("$(backend_label(B)) target requires a functional backend/device on this host")
    Backends.disable_scalar_backend!(B)
    return nothing
end

function run_gpu_backend_target(::Type{B}) where {B<:Backends.GPUBackendTag}
    require_backend_target!(B)
    @testset "$(backend_label(B)) hardware target" begin
        run_captured_graph_execution_smoke(B)
        run_algorithm_graph_backend_smoke(B)
        run_graph_rng_capture_replay(B)
        run_optional_backend_smoke(B)
        run_gpu_builder_smoke(B)
        run_revolt_like_hil_backend_smoke(B)
    end
    return nothing
end

function run_captured_graph_execution_smoke(
    ::Type{B},
) where {B<:Backends.GPUBackendTag}
    BackendArray = Backends.gpu_backend_array_type(B)
    command = BackendArray(Float32[2, -1, 0.5, 1.5, -0.25] .* 1.0f-8)
    actuator_grid_indices = BackendArray(Int32[1, 3, 5, 7, 9])
    uncompensated = BackendArray(zeros(Float32, 16, 16))
    command_schema = "test.graph.captured-command.f32/1"
    surface_schema = "test.graph.captured-surface.f32/1"
    pupil_schema = "test.graph.captured-pupil-opd.f32/1"
    definition = algorithm_graph(
        (
            grid_gaussian_deformable_mirror_surface_node(
                :dm;
                resolution=16,
                telescope_diameter_m=1.22,
                actuator_count=length(command),
                actuator_axis_count=3,
                actuator_pitch=0.4f0,
                influence_width=0.2f0,
                pdm_command_schema=command_schema,
                surface_opd_schema=surface_schema,
                actuator_grid_indices_schema=
                    "test.graph.captured-grid-indices.i32/1",
            ),
            pupil_opd_composition_node(
                :compose;
                resolution=16,
                uncompensated_opd_schema=
                    "test.graph.captured-uncompensated.f32/1",
                surface_opd_schema=surface_schema,
                pupil_opd_schema=pupil_schema,
            ),
        );
        name=:captured_gpu_dm,
        inputs=(
            graph_input(:command, :dm => :pdm_command, command),
            graph_input(
                :uncompensated,
                :compose => :uncompensated_opd,
                uncompensated,
            ),
        ),
        outputs=(
            graph_output(:surface, :dm => :surface_opd),
            graph_output(:pupil_opd, :compose => :pupil_opd),
        ),
        links=(link(:dm => :surface_opd, :compose => :surface_opd),),
        parameters=(sparse_parameter(
            :dm => :actuator_grid_indices,
            actuator_grid_indices,
        ),),
    )
    graph = prepare_algorithm_graph(
        definition;
        target=compute_device(command),
        execution=CapturedGraphExecution(),
    )
    @test graph_execution_policy(graph) isa CapturedGraphExecution
    @test captured_graph_node_count(graph) == 2
    @test all(iszero, Array(graph_output(graph, Val(:surface))))

    ticket = step_graph_async!(graph)
    @test graph_step_pending(graph)
    wait_graph_step!(ticket)
    first_surface = Array(graph_output(graph, Val(:surface)))
    @test !all(iszero, first_surface)
    @test Array(graph_output(graph, Val(:pupil_opd))) == first_surface

    copyto!(command, Float32[-1, 2, -0.5, 0.25, 1] .* 1.0f-8)
    Backends.synchronize_backend!(Backends.execution_style(command))
    dm_owner = AlgorithmGraphs.prepared_graph_node(graph, Val(:dm))
    AlgorithmGraphs.step_graph_node!(dm_owner)
    direct_surface = Array(graph_output(graph, Val(:surface)))
    step_graph!(graph)
    second_surface = Array(graph_output(graph, Val(:surface)))
    @test second_surface != first_surface
    @test second_surface ≈ direct_surface rtol = 2.0f-5 atol = 1.0f-12
    @test Array(graph_output(graph, Val(:pupil_opd))) == second_surface

    reset_graph!(graph)
    @test captured_graph_node_count(graph) == 2
    @test all(iszero, Array(graph_output(graph, Val(:surface))))
    return nothing
end

function run_algorithm_graph_backend_smoke(
    ::Type{B},
) where {B<:Backends.GPUBackendTag}
    BackendArray = Backends.gpu_backend_array_type(B)
    residual = BackendArray(Float32[1, 2])
    basis_host = zeros(Float32, 2, 2, 2)
    fill!(@view(basis_host[:, :, 1]), 1.0f0)
    fill!(@view(basis_host[:, :, 2]), 2.0f0)
    basis = BackendArray(basis_host)
    pupil_support = BackendArray(Bool[true false; true true])
    coefficient_schema = "test.graph.modal-coefficients.f32/1"
    target = compute_device(residual)

    definition = algorithm_graph(
        (
            discrete_integrator_node(
                :controller;
                extent=2,
                sample_period_s=0.1f0,
                input_schema="test.graph.residual.f32/1",
                output_schema=coefficient_schema,
                gain=2.0f0,
                tau_s=0.2f0,
            ),
            modal_opd_expansion_node(
                :modal_opd;
                pupil_rows=2,
                pupil_columns=2,
                mode_count=2,
                coefficients_schema=coefficient_schema,
                opd_schema="test.graph.opd.f32/1",
                basis_schema="test.graph.modal-basis.f32/1",
                pupil_support_schema="test.graph.pupil-support.bool/1",
            ),
        );
        name=:gpu_native_modal_control,
        inputs=(graph_input(:residual, :controller => :input, residual),),
        outputs=(graph_output(:opd, :modal_opd => :opd),),
        links=(link(:controller => :output, :modal_opd => :coefficients),),
        parameters=(
            sparse_parameter(:modal_opd => :basis, basis),
            sparse_parameter(:modal_opd => :pupil_support, pupil_support),
        ),
    )
    graph = prepare_algorithm_graph(definition; target)
    ticket = @inferred step_graph_async!(graph)
    @test graph_step_pending(graph)
    @test iszero(graph_step_sequence(graph))
    @test_throws AlgorithmGraphError step_graph_async!(graph)
    @test @inferred(wait_graph_step!(ticket)) === graph
    @test !graph_step_pending(graph)
    @test graph_step_sequence(graph) == UInt64(1)

    @test compute_device(graph) == target
    @test compute_device(graph_output(graph, Val(:opd))) == target
    @test Array(graph_output(graph, Val(:opd))) ≈ Float32[0.5 0.0; 0.5 0.5]

    atmosphere_definition = algorithm_graph(
        (
            multilayer_atmosphere_opd_node(
                :atmosphere;
                resolution=8,
                telescope_diameter_m=1.22,
                r0=0.15,
                reference_wavelength_m=500e-9,
                L0=30.0,
                fractional_cn2=(1.0,),
                wind_speed=(5.0,),
                wind_direction_deg=(0.0,),
                altitude=(0.0,),
                layer_ids=(:ground,),
                atmosphere_step=0.001,
                rng_seed=17,
                atmosphere_opd_schema=
                    "test.graph.atmosphere-opd-m.f32/1",
            ),
        );
        name=:gpu_multilayer_atmosphere_opd,
        outputs=(
            graph_output(
                :atmosphere_opd,
                :atmosphere => :atmosphere_opd,
            ),
        ),
    )
    atmosphere_graph = prepare_algorithm_graph(
        atmosphere_definition;
        target,
    )
    step_graph!(atmosphere_graph)
    atmosphere_opd = graph_output(
        atmosphere_graph,
        Val(:atmosphere_opd),
    )
    first_atmosphere_opd = Array(atmosphere_opd)
    @test compute_device(atmosphere_graph) == target
    @test compute_device(atmosphere_opd) == target
    @test all(isfinite, first_atmosphere_opd)
    @test !all(iszero, first_atmosphere_opd)
    step_graph!(atmosphere_graph)
    @test Array(atmosphere_opd) != first_atmosphere_opd
    reset_graph!(atmosphere_graph)
    step_graph!(atmosphere_graph)
    @test Array(atmosphere_opd) == first_atmosphere_opd

    pdm_command = BackendArray(Float32[2, 3])
    actuator_coordinates = BackendArray(Float32[-0.5 0.5; 0 0])
    influence_functions = BackendArray(Float32[
        1 0
        0 1
        0 1
        1 0
    ])
    uncompensated_opd = BackendArray(fill(4.0f0, 2, 2))
    dm_definition = algorithm_graph(
        (
            deformable_mirror_surface_node(
                :dm;
                resolution=2,
                telescope_diameter_m=1.22,
                actuator_count=2,
                pdm_command_schema="test.graph.pdm-command.f32/1",
                surface_opd_schema="test.graph.dm-surface-opd.f32/1",
                actuator_coordinates_schema=
                    "test.graph.dm-actuator-coordinates.f32/1",
                influence_functions_schema=
                    "test.graph.dm-influence-functions.f32/1",
            ),
            pupil_opd_composition_node(
                :composition;
                resolution=2,
                uncompensated_opd_schema=
                    "test.graph.uncompensated-opd.f32/1",
                surface_opd_schema="test.graph.dm-surface-opd.f32/1",
                pupil_opd_schema="test.graph.pupil-opd.f32/1",
            ),
        );
        name=:gpu_measured_deformable_mirror,
        inputs=(
            graph_input(:pdm_command, :dm => :pdm_command, pdm_command),
            graph_input(
                :uncompensated_opd,
                :composition => :uncompensated_opd,
                uncompensated_opd,
            ),
        ),
        outputs=(
            graph_output(:surface_opd, :dm => :surface_opd),
            graph_output(:pupil_opd, :composition => :pupil_opd),
        ),
        links=(
            link(:dm => :surface_opd, :composition => :surface_opd),
        ),
        parameters=(
            sparse_parameter(
                :dm => :actuator_coordinates,
                actuator_coordinates,
            ),
            sparse_parameter(
                :dm => :influence_functions,
                influence_functions,
            ),
        ),
    )
    dm_graph = prepare_algorithm_graph(dm_definition; target)
    dm_boundary = prepare_graph_hil_boundary(
        dm_graph;
        command_input=:pdm_command,
        frame_output=:pupil_opd,
    )
    @test step_hil_frame!(dm_boundary) == UInt64(1)
    dm_owner = AlgorithmGraphs.prepared_graph_node(dm_graph, Val(:dm))
    dm_surface = graph_output(dm_graph, Val(:surface_opd))
    pupil_opd = graph_output(dm_graph, Val(:pupil_opd))
    @test compute_device(
        influence_model(dm_owner.deformable_mirror).modes,
    ) == target
    @test influence_model(dm_owner.deformable_mirror).modes !==
        influence_functions
    @test compute_device(dm_surface) == target
    @test compute_device(pupil_opd) == target
    @test Array(dm_surface) == Float32[2 3; 3 2]
    @test Array(pupil_opd) == Float32[6 7; 7 6]
    @test hil_frame_buffer(dm_boundary) == Float32[6 7; 7 6]
    hil_command_buffer(dm_boundary) .= Float32[1, 1]
    @test adopt_hil_command!(dm_boundary, UInt64(1)) == UInt64(1)
    @test Array(graph_input(dm_graph, Val(:pdm_command))) == Float32[1, 1]

    gaussian_command = BackendArray(Float32[2.0f-8, -1.0f-8])
    gaussian_coordinates = BackendArray(Float32[-0.25 0.25; -0.25 0.25])
    gaussian_definition = algorithm_graph(
        (
            gaussian_deformable_mirror_surface_node(
                :dm;
                resolution=4,
                telescope_diameter_m=1.22,
                actuator_count=2,
                influence_width=0.2,
                pdm_command_schema=
                    "test.graph.pdm-actuator-opd-m.f32/1",
                surface_opd_schema=
                    "test.graph.dm-surface-opd-m.f32/1",
                actuator_coordinates_schema=
                    "test.graph.dm-normalized-pupil-coordinates.f32/1",
            ),
        );
        name=:gpu_analytic_gaussian_deformable_mirror,
        inputs=(
            graph_input(
                :pdm_command,
                :dm => :pdm_command,
                gaussian_command,
            ),
        ),
        outputs=(
            graph_output(:surface_opd, :dm => :surface_opd),
        ),
        parameters=(
            sparse_parameter(
                :dm => :actuator_coordinates,
                gaussian_coordinates,
            ),
        ),
    )
    gaussian_graph = prepare_algorithm_graph(gaussian_definition; target)
    step_graph!(gaussian_graph)
    gaussian_surface = graph_output(
        gaussian_graph,
        Val(:surface_opd),
    )
    gaussian_surface_host = Array(gaussian_surface)
    diagonal_response = exp(-0.5f0 / (2 * 0.2f0^2))
    cross_response = exp(-0.25f0 / (2 * 0.2f0^2))
    @test compute_device(gaussian_surface) == target
    @test gaussian_surface_host[2, 2] ≈
        2.0f-8 - 1.0f-8 * diagonal_response
    @test gaussian_surface_host[3, 3] ≈
        2.0f-8 * diagonal_response - 1.0f-8
    @test gaussian_surface_host[2, 3] ≈ 1.0f-8 * cross_response

    cmos_photon_rate = BackendArray(fill(1_000.0f0, 4, 4))
    cmos_target = compute_device(cmos_photon_rate)
    cmos_definition = algorithm_graph(
        (
            cmos_detector_acquisition_node(
                :cmos_detector;
                rows=4,
                columns=4,
                pixel_scale_arcsec=0.1,
                wavelength_m=0.8e-6,
                exposure_duration_s=0.01,
                quantum_efficiency=0.5,
                gain=2,
                bits=12,
                full_well_e=1_000,
                photon_noise=false,
                rng_seed=0x7ff,
                photon_rate_schema="test.graph.shwfs-photon-rate.f32/1",
                frame_schema="test.graph.shwfs-frame.f32/1",
            ),
        );
        name=:gpu_cmos_detector_acquisition,
        inputs=(
            graph_input(
                :photon_rate,
                :cmos_detector => :photon_rate,
                cmos_photon_rate,
            ),
        ),
        outputs=(
            graph_output(:frame, :cmos_detector => :frame),
        ),
    )
    cmos_graph = prepare_algorithm_graph(cmos_definition; target=cmos_target)
    step_graph!(cmos_graph)
    cmos_frame = graph_output(cmos_graph, Val(:frame))
    @test compute_device(cmos_frame) == cmos_target
    @test Array(cmos_frame) == fill(41.0f0, 4, 4)

    stochastic_rate = BackendArray(fill(1_000.0f0, 8, 8))
    stochastic_definition = algorithm_graph(
        (
            ccd_detector_acquisition_node(
                :detector;
                rows=8,
                columns=8,
                pixel_scale_arcsec=0.1,
                wavelength_m=0.75e-6,
                exposure_duration_s=0.01,
                quantum_efficiency=0.5,
                readout_noise_e=2.0,
                photon_noise=true,
                readout_noise=true,
                rng_seed=0x810,
                photon_rate_schema=
                    "test.graph.shwfs-photon-rate.f32/1",
                frame_schema="test.graph.shwfs-frame.f32/1",
            ),
        );
        name=:gpu_stochastic_ccd,
        inputs=(
            graph_input(
                :photon_rate,
                :detector => :photon_rate,
                stochastic_rate,
            ),
        ),
        outputs=(graph_output(:frame, :detector => :frame),),
    )
    stochastic_graph = prepare_algorithm_graph(
        stochastic_definition;
        target=compute_device(stochastic_rate),
    )
    stochastic_owner = AlgorithmGraphs.prepared_graph_node(
        stochastic_graph,
        Val(:detector),
    )
    if graph_rng_device_resident(B)
        @test stochastic_owner.rng isa Backends._PreparedCounterRNG
        @test compute_device(
            stochastic_owner.rng.state.draw_sequence,
        ) == compute_device(stochastic_rate)
    else
        @test stochastic_owner.rng isa SplitMix64RNG
    end
    step_graph!(stochastic_graph)
    stochastic_frame = graph_output(stochastic_graph, Val(:frame))
    first_stochastic_frame = Array(stochastic_frame)
    step_graph!(stochastic_graph)
    @test Array(stochastic_frame) != first_stochastic_frame
    reset_graph!(stochastic_graph)
    step_graph!(stochastic_graph)
    @test Array(stochastic_frame) == first_stochastic_frame

    pyramid_opd = BackendArray(zeros(Float32, 8, 8))
    pyramid_target = compute_device(pyramid_opd)
    pyramid_definition = algorithm_graph(
        (
            pyramid_rate_node(
                :pwfs;
                resolution=8,
                telescope_diameter_m=1.22,
                pupil_samples=4,
                modulation=2,
                modulation_points=4,
                n_pix_separation=2,
                n_pix_edge=1,
                source_wavelength_m=0.55e-6,
                source_photon_irradiance_m2_s=1.0,
                opd_schema="test.graph.pupil-opd.f32/1",
                photon_rate_schema=
                    "test.graph.pwfs-photon-rate.f32/1",
            ),
            emccd_detector_acquisition_node(
                :pwfs_detector;
                rows=12,
                columns=12,
                normalized_pupil_sampling=0.25,
                wavelength_m=0.55e-6,
                exposure_duration_s=0.5,
                quantum_efficiency=0.5,
                rng_seed=0x800,
                photon_noise=false,
                photon_rate_schema=
                    "test.graph.pwfs-photon-rate.f32/1",
                frame_schema="test.graph.pwfs-frame.f32/1",
            ),
        );
        name=:gpu_pwfs,
        inputs=(graph_input(:pupil_opd, :pwfs => :opd, pyramid_opd),),
        outputs=(
            graph_output(:pwfs_photon_rate, :pwfs => :photon_rate),
            graph_output(:pwfs_frame, :pwfs_detector => :frame),
        ),
        links=(
            link(
                :pwfs => :photon_rate,
                :pwfs_detector => :photon_rate,
            ),
        ),
    )
    pyramid_graph = prepare_algorithm_graph(
        pyramid_definition;
        target=pyramid_target,
    )
    step_graph!(pyramid_graph)
    pyramid_rate = graph_output(pyramid_graph, Val(:pwfs_photon_rate))
    pyramid_frame = graph_output(pyramid_graph, Val(:pwfs_frame))
    @test compute_device(pyramid_rate) == pyramid_target
    @test compute_device(pyramid_frame) == pyramid_target
    @test size(pyramid_rate) == (12, 12)
    @test size(pyramid_frame) == (12, 12)
    @test sum(Array(pyramid_rate)) > 0
    @test sum(Array(pyramid_frame)) ≈ sum(Array(pyramid_rate)) * 0.25f0

    pupil_opd = BackendArray(zeros(Float32, 16, 16))
    shwfs_target = compute_device(pupil_opd)
    reconstruction_matrix = BackendArray(Float32[
        1 0 0 0 0 0
        0 0 0 0 1 0
    ])
    controller_constraint_feedback = BackendArray(zeros(Float32, 2))
    shwfs_definition = algorithm_graph(
        (
            shack_hartmann_rate_node(
                :shwfs;
                resolution=16,
                telescope_diameter_m=8.0,
                n_lenslets=4,
                n_pix_subap=6,
                pixel_scale_arcsec=0.1,
                source_wavelength_m=0.75e-6,
                source_photon_irradiance_m2_s=1.0,
                opd_schema="test.graph.pupil-opd.f32/1",
                photon_rate_schema="test.graph.shwfs-photon-rate.f32/1",
            ),
            ccd_detector_acquisition_node(
                :detector;
                rows=24,
                columns=24,
                pixel_scale_arcsec=0.1,
                wavelength_m=0.75e-6,
                exposure_duration_s=0.5,
                quantum_efficiency=0.5,
                rng_seed=0x801,
                photon_rate_schema="test.graph.shwfs-photon-rate.f32/1",
                frame_schema="test.graph.shwfs-frame.f32/1",
                photon_noise=false,
            ),
            shack_hartmann_centroid_node(
                :centroid;
                resolution=16,
                telescope_diameter_m=8.0,
                n_lenslets=4,
                n_pix_subap=6,
                centroid_cutoff_fraction=0.1,
                centroid_response=1.0,
                calibration_wavelength_m=0.75e-6,
                calibration_signature=0x802,
                frame_schema="test.graph.shwfs-frame.f32/1",
                slopes_schema=
                    "test.graph.shwfs-centroid-slopes.f32/1",
                valid_subapertures_schema=
                    "test.graph.shwfs-valid-subapertures.bool/1",
                reference_signal_schema=
                    "test.graph.shwfs-centroid-reference.f32/1",
            ),
            shack_hartmann_slope_selection_node(
                :slope_selection;
                n_lenslets=4,
                selected_lenslet_count=3,
                full_slopes_schema=
                    "test.graph.shwfs-centroid-slopes.f32/1",
                selected_slopes_schema=
                    "test.graph.shwfs-selected-slopes.f32/1",
                lenslet_order_schema=
                    "test.graph.shwfs-lenslet-order.u32/1",
            ),
            control_matrix_reconstruction_node(
                :reconstruction;
                slope_count=6,
                reconstructed_count=2,
                slopes_schema=
                    "test.graph.shwfs-selected-slopes.f32/1",
                reconstructed_schema=
                    "test.graph.controller-error.f32/1",
                control_matrix_schema=
                    "test.graph.shwfs-control-matrix.f32/1",
            ),
            closed_loop_correction_node(
                :controller;
                extent=2,
                residual_error_schema=
                    "test.graph.controller-error.f32/1",
                constraint_feedback_schema=
                    "test.graph.controller-constraint-feedback.f32/1",
                correction_schema=
                    "test.graph.controller-command.f32/1",
                controller_state_schema=
                    "test.graph.controller-state.f32/1",
                gain=-0.3f0,
                pole=0.99f0,
                anti_windup_gain=1.0f0,
            ),
        );
        name=:gpu_shwfs_centroid,
        inputs=(
            graph_input(:pupil_opd, :shwfs => :opd, pupil_opd),
            graph_input(
                :controller_constraint_feedback,
                :controller => :constraint_feedback,
                controller_constraint_feedback,
            ),
        ),
        outputs=(
            graph_output(:shwfs_photon_rate, :shwfs => :photon_rate),
            graph_output(:shwfs_frame, :detector => :frame),
            graph_output(:shwfs_full_slopes, :centroid => :slopes),
            graph_output(
                :shwfs_selected_slopes,
                :slope_selection => :selected_slopes,
            ),
            graph_output(
                :controller_residual_error,
                :reconstruction => :reconstructed,
            ),
            graph_output(
                :controller_correction,
                :controller => :correction,
            ),
            graph_output(
                :controller_state,
                :controller => :controller_state,
            ),
        ),
        links=(
            link(:shwfs => :photon_rate, :detector => :photon_rate),
            link(:detector => :frame, :centroid => :frame),
            link(
                :centroid => :slopes,
                :slope_selection => :full_slopes,
            ),
            link(
                :slope_selection => :selected_slopes,
                :reconstruction => :slopes,
            ),
            link(
                :reconstruction => :reconstructed,
                :controller => :residual_error,
            ),
        ),
        parameters=(
            sparse_parameter(
                :centroid => :valid_subapertures,
                BackendArray(fill(true, 4, 4)),
            ),
            sparse_parameter(
                :centroid => :reference_signal,
                BackendArray(zeros(Float32, 16, 2)),
            ),
            sparse_parameter(
                :slope_selection => :lenslet_order,
                BackendArray(UInt32[16, 1, 6]),
            ),
            sparse_parameter(
                :reconstruction => :control_matrix,
                reconstruction_matrix,
            ),
        ),
    )
    shwfs_graph = prepare_algorithm_graph(
        shwfs_definition;
        target=shwfs_target,
    )
    step_graph!(shwfs_graph)
    photon_rate = graph_output(shwfs_graph, Val(:shwfs_photon_rate))
    frame = graph_output(shwfs_graph, Val(:shwfs_frame))
    full_slopes = graph_output(shwfs_graph, Val(:shwfs_full_slopes))
    selected_slopes = graph_output(
        shwfs_graph,
        Val(:shwfs_selected_slopes),
    )
    controller_residual_error = graph_output(
        shwfs_graph,
        Val(:controller_residual_error),
    )
    controller_correction = graph_output(
        shwfs_graph,
        Val(:controller_correction),
    )
    controller_state = graph_output(shwfs_graph, Val(:controller_state))
    @test compute_device(photon_rate) == shwfs_target
    @test compute_device(frame) == shwfs_target
    @test compute_device(full_slopes) == shwfs_target
    @test compute_device(selected_slopes) == shwfs_target
    @test compute_device(controller_residual_error) == shwfs_target
    @test compute_device(controller_correction) == shwfs_target
    @test compute_device(controller_state) == shwfs_target
    @test size(photon_rate) == (24, 24)
    @test size(frame) == (24, 24)
    @test size(full_slopes) == (32,)
    @test size(selected_slopes) == (6,)
    @test size(controller_residual_error) == (2,)
    @test size(controller_correction) == (2,)
    @test size(controller_state) == (2,)
    @test sum(Array(photon_rate)) > 0
    @test sum(Array(frame)) ≈ sum(Array(photon_rate)) * 0.25f0
    @test all(isfinite, Array(full_slopes))
    full_slopes_host = Array(full_slopes)
    @test Array(selected_slopes) == Float32[
        full_slopes_host[16],
        full_slopes_host[32],
        full_slopes_host[1],
        full_slopes_host[17],
        full_slopes_host[6],
        full_slopes_host[22],
    ]
    @test Array(controller_residual_error) ==
        Array(selected_slopes)[[1, 5]]
    @test Array(controller_correction) ≈
        -0.3f0 .* Array(controller_residual_error)
    @test all(iszero, Array(controller_state))
    step_graph!(shwfs_graph)
    @test all(isfinite, Array(full_slopes))
    @test all(isfinite, Array(selected_slopes))
    @test all(isfinite, Array(controller_residual_error))
    @test all(isfinite, Array(controller_correction))
    @test all(isfinite, Array(controller_state))
    return nothing
end

"""
    run_gpu_detector_target(::Type{<:Backends.GPUBackendTag})

Run the maintained detector qualification surface on a required hardware
backend. This target deliberately excludes unrelated optics, WFS, control, and
orchestration checks so a failure in another subsystem cannot hide detector
qualification evidence. The complete hardware target remains the release
integration gate.
"""
function run_gpu_detector_target(::Type{B}) where {B<:Backends.GPUBackendTag}
    require_backend_target!(B)
    BackendArray = Backends.gpu_backend_array_type(B)
    @test BackendArray !== nothing

    @testset "$(backend_label(B)) detector hardware target" begin
        run_optional_cmos_family_checks(B, BackendArray)
        run_optional_shared_detector_ipc_checks(B, BackendArray)
        run_optional_avalanche_detector_parity(B, BackendArray)
        run_optional_skipper_ccd_checks(B, BackendArray)
        run_optional_ingaas_checks(B, BackendArray)
        run_optional_spad_qualification_checks(B, BackendArray)
        run_optional_counting_detector_parity(B, BackendArray)
    end
    return nothing
end

function run_revolt_like_hil_backend_smoke(
    ::Type{B},
) where {B<:Backends.GPUBackendTag}
    config_dir = normpath(joinpath(@__DIR__, "..", "benchmarks", "assets", "revolt_like"))
    backend_name = lowercase(backend_label(B))
    cpu_ctx = build_revolt_like_hil_context(;
        backend_name="cpu",
        config_dir,
        sensor=CMOSSensor(T=Float32),
        T=Float32,
        rng=runtime_rng(20260713),
    )
    ctx = build_revolt_like_hil_context(;
        backend_name,
        config_dir,
        sensor=CMOSSensor(T=Float32),
        T=Float32,
        rng=runtime_rng(20260713),
    )
    revolt_like_step!(cpu_ctx)
    revolt_like_step!(ctx)
    @test size(ctx.tiled_frame) == (352, 352)
    @test all(isfinite, Array(ctx.tiled_frame))
    @test sum(Array(ctx.tiled_frame)) > 0
    @test isapprox(Array(ctx.tiled_frame), cpu_ctx.tiled_frame;
        rtol=1f-5, atol=1f-3)
    return nothing
end
