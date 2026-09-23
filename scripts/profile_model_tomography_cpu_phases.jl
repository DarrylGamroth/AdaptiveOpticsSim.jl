"""
    julia --project=. scripts/profile_model_tomography_cpu_phases.jl

Opt-in, single-process CPU profile of the model-tomography cold build.
The input values deliberately match `gpu_profile_model_tomography_phases_contract.jl`;
the setup is repeated here so the GPU contract remains unchanged. The CPU solve
phase follows the production AdaptiveOpticsCalibration path, not the GPU solve.
"First use" includes type-specific compilation after package loading, not
process startup or package-load time. Phase measurements are warmed wall time.
Set `AOS_CPU_PROFILE_REPEATS` to override the adaptive profile repeat count.
Set `AOS_TOMO_PROFILE_LENSLETS` to change the lenslet count (default 3).
Pass `--float32-only` to omit the optional Float64 case.
"""

using AdaptiveOpticsSim
using AdaptiveOpticsSim.Tomography
using LinearAlgebra
using Profile
using SparseArrays

const AOS = AdaptiveOpticsSim
const PROFILE_TARGET_SECONDS = 3.0

function profile_case(::Type{T}) where {T<:AbstractFloat}
    n_lenslets = parse(Int, get(ENV, "AOS_TOMO_PROFILE_LENSLETS", "3"))
    n_lenslets >= 1 || error("AOS_TOMO_PROFILE_LENSLETS must be positive")
    n_lgs = 2
    n_dm = 2
    grid_side = 2 * n_lenslets
    atmosphere = TomographyAtmosphereParams(
        zenith_angle_deg=T(0),
        layer_altitudes_m=T[0, 10_000],
        L0=T(25),
        r0_zenith=T(0.2),
        fractional_cn2=T[0.6, 0.4],
        reference_wavelength_m=T(500e-9),
        wind_direction_deg=T[0, 45],
        wind_speed=T[10, 20],
    )
    asterism = LGSAsterismParams(
        radius_arcsec=T(7.6),
        wavelength_m=T(589e-9),
        base_height_m=T(90_000),
        n_lgs=n_lgs,
    )
    wfs = LGSWFSParams(
        pupil_diameter_m=T(8),
        n_lenslets=n_lenslets,
        n_px=8,
        field_stop_size_arcsec=T(2),
        valid_lenslet_map=trues(n_lenslets, n_lenslets),
        lenslet_grid_rotations_rad=zeros(T, n_lgs),
        lenslet_grid_offsets_fraction=zeros(T, 2, n_lgs),
    )
    tomography = TomographyParams(
        n_fit_src=2,
        fov_optimization_arcsec=T(15),
        fit_src_height_m=T(Inf),
    )
    dm = TomographyDMParams(
        heights_m=collect(T, range(T(0), length=n_dm, step=T(6000))),
        pitch_m=fill(T(8 / n_lenslets), n_dm),
        cross_coupling=T(0.2),
        n_actuators=fill(grid_side, n_dm),
        valid_actuators=trues(grid_side, grid_side),
    )
    return (; atmosphere, asterism, wfs, tomography, dm,
        noise_model=AOS.Tomography.RelativeSignalNoise(T(0.1)),
        build_backend=AOS.Calibration.CPUBuildBackend())
end

function build_model(case)
    return build_reconstructor(
        ModelBasedTomography(), case.atmosphere, case.asterism, case.wfs,
        case.tomography, case.dm;
        noise_model=case.noise_model, build_backend=case.build_backend)
end

function measured_phase!(f, rows, name::Symbol)
    result = @timed f()
    push!(rows, (name, result.time * 1e3, result.bytes, result.gctime * 1e3))
    return result.value
end

function build_phases(case)
    T = eltype(case.atmosphere.layer_altitudes_m)
    tomo = AOS.Tomography
    calibration = AOS.Calibration
    rows = Tuple{Symbol,Float64,Int,Float64}[]
    support = tomo.valid_lenslet_support(case.wfs)
    gamma_single, grid_mask = measured_phase!(rows, :gamma_single) do
        tomo.sparse_gradient_matrix(support; over_sampling=2)
    end
    gamma_t = measured_phase!(rows, :gamma_convert) do
        SparseMatrixCSC{T,Int}(gamma_single)
    end
    gamma = measured_phase!(rows, :gamma_blockdiag) do
        blockdiag(ntuple(_ -> gamma_t, case.asterism.n_lgs)...)
    end
    cxx = measured_phase!(rows, :auto_correlation) do
        tomo.auto_correlation(case.build_backend, case.atmosphere,
            case.asterism, case.wfs, grid_mask)
    end
    cross = measured_phase!(rows, :cross_correlation) do
        tomo.cross_correlation(case.build_backend, case.atmosphere,
            case.asterism, case.wfs, case.tomography; grid_mask=grid_mask)
    end
    weights = tomo._equal_fit_source_weights(case.tomography)
    cox = measured_phase!(rows, :fit_source_average) do
        tomo._fit_source_average(cross, weights)
    end
    gamma_native = measured_phase!(rows, :materialize_gamma) do
        calibration.materialize_build(case.build_backend, gamma, gamma)
    end
    cxx_native = measured_phase!(rows, :materialize_cxx) do
        calibration.materialize_build(case.build_backend, gamma_native, cxx)
    end
    cox_native = measured_phase!(rows, :materialize_cox) do
        calibration.materialize_build(case.build_backend, gamma_native, cox)
    end
    measured_phase!(rows, :materialize_mask) do
        calibration.materialize_build(case.build_backend, gamma_native, grid_mask)
    end
    css_signal = measured_phase!(rows, :css_signal) do
        tomo.backend_symmetric_product(gamma_native, cxx_native)
    end
    reference_diag = measured_phase!(rows, :reference_diag) do
        tomo.tomography_reference_diagonal(case.build_backend, css_signal)
    end
    cnz = measured_phase!(rows, :cnz) do
        tomo.tomography_noise_covariance(case.build_backend, case.noise_model,
            reference_diag)
    end
    recstat = measured_phase!(rows, :aoc_covariance_solve) do
        tomo._tomographic_covariance_reconstructor(
            tomo.execution_style(cxx_native), case.build_backend, gamma_native,
            cxx_native, cox_native, cnz)
    end
    d = tomo.lenslet_grid_support_diameter_m(case.wfs) / size(support, 1)
    wavefront_to_meter = case.asterism.wavelength_m / d / 2
    recon = measured_phase!(rows, :recon_scale) do
        d * wavefront_to_meter .* recstat
    end
    return recon, rows
end

function print_measurement(label, timed)
    println("  ", label, ": elapsed_ms=", round(timed.time * 1e3; digits=3),
        " bytes=", timed.bytes,
        " gc_ms=", round(timed.gctime * 1e3; digits=3))
end

function profile_text(format; sortedby=:count)
    return sprint() do io
        Profile.print(io; format=format, sortedby=sortedby, C=true,
            mincount=5, noisefloor=2.0, maxdepth=45)
    end
end

function print_profile_summary()
    tree = profile_text(:tree)
    snapshot_match = match(r"Total snapshots: (\d+)", tree)
    println("  profile_snapshots: ", isnothing(snapshot_match) ? "unavailable" :
        snapshot_match.captures[1])
    domain_lines = filter(split(tree, '\n')) do line
        occursin("AdaptiveOpticsSim/", line) ||
            occursin("AdaptiveOpticsCalibration/", line) ||
            occursin("SpecialFunctions/", line) ||
            occursin("libopenspecfun", line)
    end
    println("  profile_tree_domain_frames (first ", min(30, length(domain_lines)), " lines):")
    for line in Iterators.take(domain_lines, 30)
        println("    ", line)
    end

    flat = profile_text(:flat; sortedby=:overhead)
    lines = filter(line -> occursin(r"^\s*\d+\s+\d+", line), split(flat, '\n'))
    println("  profile_top_self_frames (highest overhead, up to 25 lines):")
    for line in Iterators.take(Iterators.reverse(lines), 25)
        println("    ", line)
    end
end

function run_profile(::Type{T}) where {T<:AbstractFloat}
    case = profile_case(T)
    println("CPU model tomography cold profile")
    println("  Julia: ", VERSION, "; type: ", T,
        "; Julia threads: ", Threads.nthreads(),
        "; BLAS threads: ", BLAS.get_num_threads())
    println("  case: ", case.wfs.n_lenslets,
        " lenslets, 2 LGS, 2 fit sources, 2 layers, 2 DMs, CPU build")

    first = @timed build_model(case)
    print_measurement("first_use_full_build", first)
    println("  shapes: reconstructor=", size(first.value.reconstructor),
        " grid_mask=", size(first.value.grid_mask),
        " Cxx=", size(first.value.operators.cxx),
        " Cox=", size(first.value.operators.cox))
    warm = [@timed build_model(case) for _ in 1:3]
    for (i, measurement) in enumerate(warm)
        print_measurement("warmed_full_build_$i", measurement)
    end

    # Compile the measurement wrapper before reporting phase-level steady state.
    phase_recon, _ = build_phases(case)
    @assert isapprox(phase_recon, first.value.reconstructor; rtol=1e-5, atol=0)
    for i in 1:3
        _, rows = build_phases(case)
        println("  warmed_phases_$i (elapsed_ms, bytes, gc_ms):")
        for (name, elapsed_ms, bytes, gc_ms) in rows
            println("    ", name, ": ", round(elapsed_ms; digits=3),
                ", ", bytes, ", ", round(gc_ms; digits=3))
        end
        println("    timed_phase_sum: ", round(sum(row[2] for row in rows); digits=3),
            ", ", sum(row[3] for row in rows),
            ", ", round(sum(row[4] for row in rows); digits=3))
    end

    warmed_median = sort([measurement.time for measurement in warm])[2]
    default_repeats = clamp(ceil(Int, PROFILE_TARGET_SECONDS / max(warmed_median, 1e-3)), 10, 200)
    repeats = parse(Int, get(ENV, "AOS_CPU_PROFILE_REPEATS", string(default_repeats)))
    repeats >= 1 || error("AOS_CPU_PROFILE_REPEATS must be positive")
    Profile.init(n=10^7, delay=0.001)
    Profile.clear()
    Profile.@profile for _ in 1:repeats
        build_model(case)
    end
    raw = Profile.fetch()
    println("  profile_repeats: ", repeats,
        "; sample_buffer_words: ", length(raw),
        "; delay_s: 0.001")
    print_profile_summary()
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_profile(Float32)
    if !("--float32-only" in ARGS)
        run_profile(Float64)
    end
end
