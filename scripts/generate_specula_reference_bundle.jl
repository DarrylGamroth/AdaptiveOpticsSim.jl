using AdaptiveOpticsSim
using DelimitedFiles
using TOML

include(joinpath(@__DIR__, "..", "test", "reference_harness.jl"))

function specula_reference_cases()
    return Dict{String,Dict{String,Any}}(
        "zernike_flat_signal" => Dict(
            "baseline" => "specula",
            "kind" => "zernike_signal",
            "data" => "zernike_flat_signal.txt",
            "atol" => 1e-10,
            "rtol" => 1e-10,
            "specula_test" => "test_zernike_sensor.py::TestZernikeSensor::test_flat_wavefront_output_size",
            "specula_contract" => "flat output shape and non-negative signal surface",
            "telescope" => Dict(
                "resolution" => 48,
                "diameter" => 8.0,
                "sampling_time" => 1e-3,
                "central_obstruction" => 0.0,
            ),
            "source" => Dict(
                "kind" => "ngs",
                "band" => "I",
                "magnitude" => 0.0,
            ),
            "wfs" => Dict(
                "n_subap" => 24,
                "threshold" => 0.0,
                "mode" => "diffractive",
                "phase_shift_pi" => 0.5,
                "spot_radius_lambda_over_d" => 1.0,
                "diffraction_padding" => 2,
                "binning" => 2,
            ),
            "opd" => Dict(
                "kind" => "zeros",
            ),
        ),
        "zernike_quadratic_signal" => Dict(
            "baseline" => "specula",
            "kind" => "zernike_signal",
            "data" => "zernike_quadratic_signal.txt",
            "atol" => 1e-10,
            "rtol" => 1e-10,
            "specula_test" => "test_zernike_sensor.py::TestZernikeSensor::test_focus",
            "specula_contract" => "low-order quadratic phase produces a structured non-flat differential signal",
            "telescope" => Dict(
                "resolution" => 48,
                "diameter" => 8.0,
                "sampling_time" => 1e-3,
                "central_obstruction" => 0.0,
            ),
            "source" => Dict(
                "kind" => "ngs",
                "band" => "I",
                "magnitude" => 0.0,
            ),
            "basis" => Dict(
                "kind" => "cartesian_polynomials",
                "n_modes" => 4,
            ),
            "wfs" => Dict(
                "n_subap" => 24,
                "threshold" => 0.0,
                "mode" => "diffractive",
                "phase_shift_pi" => 0.5,
                "spot_radius_lambda_over_d" => 1.0,
                "diffraction_padding" => 2,
                "binning" => 2,
            ),
            "opd" => Dict(
                "kind" => "basis_mode",
                "mode_index" => 4,
                "amplitude" => 5e-8,
            ),
        ),
        "curvature_flat_signal" => Dict(
            "baseline" => "specula",
            "kind" => "curvature_signal",
            "data" => "curvature_flat_signal.txt",
            "atol" => 1e-10,
            "rtol" => 1e-10,
            "specula_test" => "test_curvature_sensor.py::TestCurvatureSensor::test_flat_wavefront_flux",
            "specula_contract" => "flat wavefront produces zero curvature signal after reference subtraction",
            "telescope" => Dict(
                "resolution" => 48,
                "diameter" => 8.0,
                "sampling_time" => 1e-3,
                "central_obstruction" => 0.0,
            ),
            "source" => Dict(
                "kind" => "ngs",
                "band" => "I",
                "magnitude" => 0.0,
            ),
            "wfs" => Dict(
                "n_subap" => 12,
                "threshold" => 0.0,
                "mode" => "diffractive",
                "defocus_rms_nm" => 500.0,
                "diffraction_padding" => 2,
                "readout_model" => "frame",
                "readout_crop_resolution" => 48,
                "readout_pixels_per_subap" => 1,
            ),
            "opd" => Dict(
                "kind" => "zeros",
            ),
        ),
        "curvature_quadratic_signal" => Dict(
            "baseline" => "specula",
            "kind" => "curvature_signal",
            "data" => "curvature_quadratic_signal.txt",
            "atol" => 1e-10,
            "rtol" => 1e-10,
            "specula_test" => "test_curvature_sensor.py::TestCurvatureSensor::test_full_chain_focus_response",
            "specula_contract" => "low-order quadratic phase produces a non-trivial curvature signal",
            "telescope" => Dict(
                "resolution" => 48,
                "diameter" => 8.0,
                "sampling_time" => 1e-3,
                "central_obstruction" => 0.0,
            ),
            "source" => Dict(
                "kind" => "ngs",
                "band" => "I",
                "magnitude" => 0.0,
            ),
            "basis" => Dict(
                "kind" => "cartesian_polynomials",
                "n_modes" => 4,
            ),
            "wfs" => Dict(
                "n_subap" => 12,
                "threshold" => 0.0,
                "mode" => "diffractive",
                "defocus_rms_nm" => 500.0,
                "diffraction_padding" => 2,
                "readout_model" => "frame",
                "readout_crop_resolution" => 48,
                "readout_pixels_per_subap" => 1,
            ),
            "opd" => Dict(
                "kind" => "basis_mode",
                "mode_index" => 4,
                "amplitude" => 5e-8,
            ),
        ),
    )
end

function build_specula_manifest(root::AbstractString)
    cases = specula_reference_cases()
    manifest = Dict{String,Any}(
        "version" => 1,
        "metadata" => Dict(
            "baseline" => "specula",
            "bundle_kind" => "contract_bundle",
            "description" => "Frozen SPECULA-targeted reference cases derived from maintained local scenarios aligned to SPECULA test contracts.",
            "specula_repo_path" => abspath(joinpath(@__DIR__, "..", "..", "SPECULA")),
        ),
        "cases" => Dict{String,Any}(),
    )
    mkpath(root)
    for id in sort!(collect(keys(cases)))
        case_cfg = deepcopy(cases[id])
        kind = Symbol(case_cfg["kind"])
        baseline = parse_reference_baseline(kind, case_cfg)
        case = ReferenceCase(id, kind, baseline, joinpath(root, case_cfg["data"]), (), 0.0, 0.0, case_cfg)
        actual = compute_reference_actual(case)
        case_cfg["shape"] = collect(size(actual))
        writedlm(joinpath(root, case_cfg["data"]), actual)
        manifest["cases"][id] = case_cfg
    end
    return manifest
end

function main(args)
    root = length(args) >= 1 ? abspath(args[1]) : abspath(joinpath(@__DIR__, "..", "test", "reference_data_specula"))
    manifest = build_specula_manifest(root)
    open(joinpath(root, "manifest.toml"), "w") do io
        TOML.print(io, manifest)
    end
    println(root)
end

main(ARGS)
