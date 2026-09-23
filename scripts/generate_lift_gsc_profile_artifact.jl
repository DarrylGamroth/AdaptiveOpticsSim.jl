using AdaptiveOpticsSim
using AdaptiveOpticsSim.Optics
using AdaptiveOpticsSim.WavefrontSensors
using AdaptiveOpticsSim.Calibration
import AdaptiveOpticsCalibration
using Random
using TOML

const ProfilePhaseRetrieval = AdaptiveOpticsCalibration.PhaseRetrieval
const ProfileOpticalGains = AdaptiveOpticsCalibration.OpticalGains

const OUTDIR = joinpath(@__DIR__, "..", "benchmarks", "results", "workflows")
const OUTFILE = joinpath(OUTDIR, "2026-09-23-lift-aoc-gsc-profile.toml")
const MANIFEST = joinpath(OUTDIR, "manifest.toml")
const ARTIFACT_ID = "WORKFLOW-VAL-2026-09-23-AOS-AOC-GSC"

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
        pupil = pupil_mask(tel)[i, j] ? 1.0 : 0.0
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

function centered_focal_basis(pupil_basis::AbstractArray{T,3},
    mask::AbstractMatrix{Complex{T}}) where {T<:AbstractFloat}
    side = size(mask, 1)
    size(mask, 2) == side || throw(ArgumentError(
        "Pyramid focal mask must be square for complete-image gain sensing",
    ))
    size(pupil_basis, 1) == size(pupil_basis, 2) || throw(ArgumentError(
        "modal pupil basis must be square for complete-image gain sensing",
    ))
    side >= size(pupil_basis, 1) || throw(ArgumentError(
        "Pyramid focal mask must not be smaller than the modal pupil basis",
    ))
    iseven(side - size(pupil_basis, 1)) || throw(ArgumentError(
        "Pyramid focal mask and modal pupil basis must have aligned centers",
    ))

    basis = zeros(T, side, side, size(pupil_basis, 3))
    offset = div(side - size(pupil_basis, 1), 2)
    @views basis[offset+1:offset+size(pupil_basis, 1),
        offset+1:offset+size(pupil_basis, 2), :] .= pupil_basis
    return basis
end

function _alloc_bytes(f!::Function)
    GC.gc()
    return @allocated f!()
end

function lift_profile()
    tel = Telescope(resolution=24, diameter=8.0, central_obstruction=0.0)
    src = Source(band=:I, magnitude=0.0)
    zb = ZernikeBasis(tel, 6)
    compute_zernike!(zb, tel)
    basis = zb.modes[:, :, 1:4]
    coeffs_true = [25e-9, -10e-9, 5e-9, 0.0]
    pupil = PupilFunction(tel)
    opd = zeros(Float64, size(pupil.opd))
    combine_modes!(opd, basis, coeffs_true)
    apply_opd!(pupil, opd)
    imaging = prepare_direct_imaging(pupil, src; zero_padding=2)
    psf = copy(intensity_values(form_direct_image!(imaging)))
    diversity = zeros(Float64, size(pupil.opd))
    model_opd = zeros(Float64, size(pupil.opd))
    t0 = time_ns()
    forward = prepare_lift_forward_model(tel, src, basis, model_opd;
        diversity_opd=diversity, zero_padding=2)
    observation = LiFTObservation(forward, psf)
    specification = ProfilePhaseRetrieval.LiFTSpecification(forward,
        observation)
    method = ProfilePhaseRetrieval.LiFT(iterations=3,
        mode_indices=1:length(coeffs_true), check_convergence=true)
    plan = AdaptiveOpticsCalibration.prepare(method, specification)
    result = AdaptiveOpticsCalibration.allocate_result(plan)
    workspace = AdaptiveOpticsCalibration.allocate_workspace(plan)
    inputs = ProfilePhaseRetrieval.LiFTInputs(observation.values)
    build_time_ns = Int(time_ns() - t0)

    reconstruct!() = AdaptiveOpticsCalibration.process!(result, workspace,
        plan, inputs)
    reconstruct!()
    timing = runtime_timing(reconstruct!;
        warmup=3, samples=20, gc_before=false)
    alloc_bytes = _alloc_bytes(reconstruct!)
    coeff_error = maximum(abs.(
        ProfilePhaseRetrieval.lift_coefficients(result) .- coeffs_true))

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
    tel = Telescope(resolution=24, diameter=8.0, central_obstruction=0.0)
    src = Source(band=:R, magnitude=8.0)
    wfs = PyramidWFS(tel; pupil_samples=4, modulation=3.0,
        modulation_points=8, diffraction_padding=2, n_pix_separation=2, n_pix_edge=1)
    pupil_basis = cartesian_basis(tel, 4)
    pupil = PupilFunction(tel)
    reference_frame_time_ns = time_ns()
    calibration_frame = pyramid_modulation_frame(wfs, pupil, src)
    reference_frame_time_ns = Int(time_ns() - reference_frame_time_ns)
    focal_mask = pyramid_focal_mask(wfs)
    focal_basis = centered_focal_basis(pupil_basis, focal_mask)
    plan_time_ns = time_ns()
    specification = ProfileOpticalGains.GainSensingSpecification(
        focal_mask, focal_basis, calibration_frame,
    )
    plan = AdaptiveOpticsCalibration.prepare(
        ProfileOpticalGains.GainSensing(), specification,
    )
    product = AdaptiveOpticsCalibration.allocate_result(plan)
    workspace = AdaptiveOpticsCalibration.allocate_workspace(plan)
    plan_build_time_ns = Int(time_ns() - plan_time_ns)

    opd = zeros(Float64, size(pupil.opd))
    combine_modes!(opd, pupil_basis, [20e-9, -12e-9, 8e-9, 0.0])
    apply_opd!(pupil, opd)
    frame = similar(calibration_frame)
    pyramid_modulation_frame!(frame, wfs, pupil, src)
    estimate!() = AdaptiveOpticsCalibration.process!(product, workspace,
        plan, ProfileOpticalGains.GainSensingInputs(frame))
    estimate!()
    og = copy(ProfileOpticalGains.optical_gains(product))

    measure_timing = runtime_timing(estimate!;
        warmup=3, samples=20, gc_before=false)
    measure_alloc = _alloc_bytes(estimate!)

    return Dict(
        "scenario" => "aos_pyramid_modulation_frame_plus_aoc_gain_sensing",
        "measurement_scope" =>
            "AOC process! on a preformed complete frame; excludes AOS frame formation",
        "reference_frame_build_time_ns" => reference_frame_time_ns,
        "aoc_plan_build_time_ns" => plan_build_time_ns,
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
        "artifact_id" => ARTIFACT_ID,
        "generated_on" => "2026-09-23",
        "scope" => Dict(
            "backend" => "cpu",
            "artifact_kind" => "workflow_profile_validation",
            "families" => ["lift", "aos_pyramid_aoc_gain_sensing"],
            "gain_sensing_ownership" =>
                "AOS forms the physical Pyramid modulation frame; AOC owns the prepared complete-image estimator",
        ),
        "cases" => Dict(
            "lift_reconstruct" => lift,
            "aos_pyramid_aoc_gain_sensing" => gsc,
        ),
    )
end

function update_manifest!(artifact_path::AbstractString)
    mkpath(dirname(MANIFEST))
    manifest = isfile(MANIFEST) ? TOML.parsefile(MANIFEST) : Dict{String,Any}()
    artifacts = get!(manifest, "artifacts", Any[])
    kept = Any[item for item in artifacts if get(item, "id", "") != ARTIFACT_ID]
    push!(kept, Dict(
        "purpose" => "LiFT and AOS-Pyramid/AOC-gain-sensing workflow profile artifact",
        "id" => ARTIFACT_ID,
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

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
