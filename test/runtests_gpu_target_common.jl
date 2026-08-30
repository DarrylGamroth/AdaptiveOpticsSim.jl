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
using AdaptiveOpticsSim: Plant
using AdaptiveOpticsSim.Plant
using FixedSizeArrays: FixedSizeVector
using LinearAlgebra
using Random
using Statistics

# Hardware targets exercise qualified-public and internal plant contracts in
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

for name in names(Plant; all=true)
    s = String(name)
    if Base.isidentifier(s) && !startswith(s, "#") && !isdefined(@__MODULE__, name)
        @eval const $(name) = getfield(Plant, $(QuoteNode(name)))
    end
end

include("plant_device_batching_fixtures.jl")
include("plant_device_model_matrix_fixtures.jl")
include("backend_optional_common.jl")
include(normpath(joinpath(@__DIR__, "..", "benchmarks", "support", "revolt_like_hil_common.jl")))
include(normpath(joinpath(@__DIR__, "..", "scripts", "gpu_builder_contract.jl")))

BLAS.set_num_threads(1)
Backends.set_fft_provider_threads!(1)

backend_target_branch_mode(::Type{Backends.CUDABackendTag}) =
    BackendStreamExecution()
backend_target_branch_mode(::Type{Backends.AMDGPUBackendTag}) =
    SequentialExecution()

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
        run_algorithm_graph_backend_smoke(B)
        run_optional_backend_smoke(B)
        run_gpu_builder_smoke(B)
        run_revolt_like_hil_backend_smoke(B)
    end
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
    step_graph!(graph)

    @test compute_device(graph) == target
    @test compute_device(graph_output(graph, Val(:opd))) == target
    @test Array(graph_output(graph, Val(:opd))) ≈ Float32[0.5 0.0; 0.5 0.5]

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
        run_optional_detector_device_model_matrix_checks(B, BackendArray)
        run_optional_cmos_family_checks(B, BackendArray)
        run_optional_shared_detector_ipc_checks(B, BackendArray)
        run_optional_detector_event_checks(B, BackendArray)
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
