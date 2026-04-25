using AdaptiveOpticsSim
using Logging
using Random
using TOML

const OUTDIR = joinpath(@__DIR__, "..", "benchmarks", "results", "workflows")
const OUTFILE = joinpath(OUTDIR, "2026-04-01-phase1-pvp05.toml")
const MANIFEST = joinpath(OUTDIR, "manifest.toml")

function combine_modes!(out::AbstractMatrix{T}, basis::AbstractArray{<:Real,3},
    coeffs::AbstractVector{<:Real}) where {T<:AbstractFloat}
    fill!(out, zero(T))
    n_modes = min(size(basis, 3), length(coeffs))
    @inbounds for k in 1:n_modes
        coeff = T(coeffs[k])
        @views @. out += coeff * basis[:, :, k]
    end
    return out
end

function cartesian_basis(tel::Telescope, n_modes::Int)
    n = tel.params.resolution
    basis = zeros(Float64, n, n, n_modes)
    x = collect(range(-1.0, 1.0; length=n + 1))[1:n]
    y = collect(range(-1.0, 1.0; length=n + 1))[1:n]
    @inbounds for j in 1:n, i in 1:n
        px = x[i]
        py = y[j]
        pupil = tel.state.pupil[i, j] ? 1.0 : 0.0
        if n_modes >= 1
            basis[i, j, 1] = pupil * px
        end
        if n_modes >= 2
            basis[i, j, 2] = pupil * py
        end
        if n_modes >= 3
            basis[i, j, 3] = pupil * px * py
        end
        if n_modes >= 4
            basis[i, j, 4] = pupil * (px * px - py * py)
        end
    end
    return basis
end

function _alloc_bytes(f!::Function)
    GC.gc()
    return @allocated f!()
end

function lift_profile()
    tel = Telescope(resolution=24, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    src = Source(band=:I, magnitude=0.0)
    zb = ZernikeBasis(tel, 6)
    compute_zernike!(zb, tel)
    basis = zb.modes[:, :, 1:4]
    coeffs_true = [25e-9, -10e-9, 5e-9, 0.0]
    opd = zeros(Float64, size(tel.state.opd))
    combine_modes!(opd, basis, coeffs_true)
    apply_opd!(tel, opd)
    psf = copy(compute_psf!(tel, src; zero_padding=2))
    det = Detector(noise=NoiseNone(), integration_time=1.0, qe=1.0, psf_sampling=2, binning=1)
    diversity = zeros(Float64, size(tel.state.opd))
    t0 = time_ns()
    lift = LiFT(tel, src, basis, det; diversity_opd=diversity, iterations=3, numerical=false)
    build_time_ns = Int(time_ns() - t0)

    mode_ids = collect(1:length(coeffs_true))
    coeffs_fit = zeros(Float64, length(coeffs_true))
    reconstruct!(coeffs_fit, lift, psf, mode_ids; coeffs0=nothing, check_convergence=true)
    timing = runtime_timing(() -> reconstruct!(coeffs_fit, lift, psf, mode_ids;
            coeffs0=nothing, check_convergence=true);
        warmup=3, samples=20, gc_before=false)
    alloc_bytes = _alloc_bytes(() -> reconstruct!(coeffs_fit, lift, psf, mode_ids;
        coeffs0=nothing, check_convergence=true))
    coeff_error = maximum(abs.(coeffs_fit .- coeffs_true))

    return Dict(
        "scenario" => "tutorial_like_modal_psf_fit",
        "build_time_ns" => build_time_ns,
        "reconstruct_mean_ns" => timing.mean_ns,
        "reconstruct_p95_ns" => timing.p95_ns,
        "reconstruct_alloc_bytes" => alloc_bytes,
        "max_coeff_error" => coeff_error,
        "n_modes" => length(coeffs_true),
        "img_resolution" => size(psf, 1),
    )
end

function gsc_profile()
    tel = Telescope(resolution=24, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    src = Source(band=:R, magnitude=8.0)
    wfs = PyramidWFS(tel; pupil_samples=4, mode=Diffractive(), threshold=0.5, modulation=3.0,
        normalization=IncidenceFluxNormalization(),
        modulation_points=8, diffraction_padding=2, n_pix_separation=2, n_pix_edge=1)
    basis = cartesian_basis(tel, 4)
    t0 = time_ns()
    gsc = GainSensingCamera(wfs.state.pyramid_mask, basis)
    build_time_ns = Int(time_ns() - t0)
    calibration_frame = similar(wfs.state.intensity)
    reset_opd!(tel)
    pyramid_modulation_frame!(calibration_frame, wfs, tel, src)
    calib_time_ns = time_ns()
    with_logger(NullLogger()) do
        calibrate!(gsc, calibration_frame)
    end
    calib_build_ns = Int(time_ns() - calib_time_ns)

    opd = zeros(Float64, size(tel.state.opd))
    combine_modes!(opd, basis, [20e-9, -12e-9, 8e-9, 0.0])
    apply_opd!(tel, opd)
    frame = similar(calibration_frame)
    pyramid_modulation_frame!(frame, wfs, tel, src)
    og = copy(compute_optical_gains!(gsc, frame))

    calib_timing = runtime_timing(() -> with_logger(NullLogger()) do
            calibrate!(gsc, calibration_frame)
        end;
        warmup=2, samples=10, gc_before=false)
    calib_alloc = _alloc_bytes(() -> with_logger(NullLogger()) do
        calibrate!(gsc, calibration_frame)
    end)

    measure_timing = runtime_timing(() -> compute_optical_gains!(gsc, frame);
        warmup=3, samples=20, gc_before=false)
    measure_alloc = _alloc_bytes(() -> compute_optical_gains!(gsc, frame))

    return Dict(
        "scenario" => "tutorial_like_pyramid_gain_sensing",
        "build_time_ns" => build_time_ns,
        "first_calibration_time_ns" => calib_build_ns,
        "calibration_mean_ns" => calib_timing.mean_ns,
        "calibration_p95_ns" => calib_timing.p95_ns,
        "calibration_alloc_bytes" => calib_alloc,
        "measurement_mean_ns" => measure_timing.mean_ns,
        "measurement_p95_ns" => measure_timing.p95_ns,
        "measurement_alloc_bytes" => measure_alloc,
        "mean_optical_gain" => sum(og) / length(og),
        "n_modes" => length(og),
        "frame_size" => collect(size(frame)),
    )
end

function build_report()
    lift = lift_profile()
    gsc = gsc_profile()
    return Dict(
        "artifact_id" => "WORKFLOW-VAL-2026-04-01",
        "generated_on" => "2026-04-01",
        "scope" => Dict(
            "backend" => "cpu",
            "artifact_kind" => "workflow_profile_validation",
            "families" => ["lift", "gain_sensing_camera"],
        ),
        "cases" => Dict(
            "lift_reconstruct" => lift,
            "gain_sensing_camera" => gsc,
        ),
    )
end

function update_manifest!(artifact_path::AbstractString)
    mkpath(dirname(MANIFEST))
    manifest = isfile(MANIFEST) ? TOML.parsefile(MANIFEST) : Dict{String,Any}()
    artifacts = get!(manifest, "artifacts", Any[])
    kept = Any[item for item in artifacts if get(item, "id", "") != "WORKFLOW-VAL-2026-04-01"]
    push!(kept, Dict(
        "purpose" => "LiFT and gain-sensing workflow profile artifact for PVP-05",
        "id" => "WORKFLOW-VAL-2026-04-01",
        "path" => basename(artifact_path),
    ))
    manifest["artifacts"] = kept
    open(MANIFEST, "w") do io
        TOML.print(io, manifest)
    end
end

function main()
    mkpath(OUTDIR)
    report = build_report()
    open(OUTFILE, "w") do io
        TOML.print(io, report)
    end
    update_manifest!(OUTFILE)
    println(OUTFILE)
end

main()
