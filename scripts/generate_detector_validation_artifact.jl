using AdaptiveOpticsSim
using Random
using Statistics
using TOML

const OUTDIR = joinpath(@__DIR__, "..", "benchmarks", "results", "detectors")
const OUTFILE = joinpath(OUTDIR, "2026-04-01-phase1-pvp04.toml")
const MANIFEST = joinpath(OUTDIR, "manifest.toml")

function ccd_emccd_case()
    zero_psf = zeros(4, 4)
    rng_ccd = MersenneTwister(7)
    rng_emccd = MersenneTwister(7)
    det_ccd = Detector(
        integration_time=1.0,
        noise=NoiseReadout(1.0),
        qe=1.0,
        binning=1,
        gain=10.0,
        sensor=CCDSensor(),
    )
    det_emccd = Detector(
        integration_time=1.0,
        noise=NoiseReadout(1.0),
        qe=1.0,
        binning=1,
        gain=10.0,
        sensor=EMCCDSensor(),
    )
    frame_ccd = copy(capture!(det_ccd, zero_psf; rng=rng_ccd))
    frame_emccd = copy(capture!(det_emccd, zero_psf; rng=rng_emccd))
    gain_error = maximum(abs.(frame_ccd .- 10.0 .* frame_emccd))
    return Dict(
        "ccd_mean" => mean(frame_ccd),
        "emccd_mean" => mean(frame_emccd),
        "gain_error_max" => gain_error,
        "contract_holds" => gain_error ≤ 1e-10,
    )
end

function cmos_structured_case()
    prnu_map = [1.0 0.5 1.0 0.5; 1.0 0.5 1.0 0.5; 1.0 0.5 1.0 0.5; 1.0 0.5 1.0 0.5]
    dsnu_map = fill(0.25, 4, 4)
    bad_mask = falses(4, 4)
    bad_mask[2, 3] = true
    det = Detector(
        integration_time=1.0,
        noise=NoiseNone(),
        qe=1.0,
        binning=1,
        sensor=CMOSSensor(
            output_model=StaticCMOSOutputPattern(2, [1.0, 2.0], [0.0, 10.0]),
            timing_model=RollingShutter(1e-3),
        ),
        response_model=NullFrameResponse(),
        defect_model=CompositeDetectorDefectModel(
            PixelResponseNonuniformity(prnu_map),
            DarkSignalNonuniformity(dsnu_map),
            BadPixelMask(bad_mask; throughput=0.0),
        ),
    )
    frame = capture!(det, fill(2.0, 4, 4); rng=MersenneTwister(120))
    return Dict(
        "pixel_1_1" => frame[1, 1],
        "pixel_1_2" => frame[1, 2],
        "pixel_1_3" => frame[1, 3],
        "pixel_2_3" => frame[2, 3],
        "mean" => mean(frame),
        "contract_holds" => frame[1, 1] == 2.25 &&
            frame[1, 2] == 1.25 &&
            frame[1, 3] == 14.5 &&
            frame[2, 3] == 10.5,
    )
end

function ingaas_case()
    det = Detector(
        integration_time=1.0,
        noise=NoiseNone(),
        qe=1.0,
        binning=1,
        response_model=NullFrameResponse(),
        sensor=InGaAsSensor(persistence_model=ExponentialPersistence(0.5, 0.0)),
    )
    capture!(det, fill(4.0, 4, 4); rng=MersenneTwister(121))
    persisted = capture!(det, zeros(4, 4); rng=MersenneTwister(122))
    return Dict(
        "persisted_sum" => sum(persisted),
        "persisted_mean" => mean(persisted),
        "contract_holds" => isapprox(sum(persisted), 32.0; atol=1e-10, rtol=1e-10),
    )
end

function hgcdte_case()
    uniform_signal = fill(50.0, 8, 8)
    det_gain = Detector(
        integration_time=1.0,
        noise=NoiseNone(),
        qe=1.0,
        binning=1,
        gain=1.0,
        sensor=HgCdTeAvalancheArraySensor(avalanche_gain=5.0),
    )
    gain_frame = copy(capture!(det_gain, uniform_signal; rng=MersenneTwister(14)))
    det_excess = Detector(
        integration_time=1.0,
        noise=NoiseNone(),
        qe=1.0,
        binning=1,
        gain=1.0,
        sensor=HgCdTeAvalancheArraySensor(avalanche_gain=1.0, excess_noise_factor=sqrt(2.0)),
    )
    excess_frame = copy(capture!(det_excess, uniform_signal; rng=MersenneTwister(14)))
    gain_error = maximum(abs.(gain_frame .- 5.0 .* uniform_signal))
    return Dict(
        "gain_error_max" => gain_error,
        "excess_std" => std(vec(excess_frame)),
        "contract_holds" => gain_error ≤ 1e-10 && std(vec(excess_frame)) > 0,
    )
end

function apd_case()
    channels = fill(2.0, 4, 4)
    det = APDDetector(
        integration_time=1.0,
        qe=1.0,
        gain=1.0,
        dark_count_rate=0.0,
        noise=NoiseNone(),
        gate_model=DutyCycleGate(0.5),
        dead_time_model=NonParalyzableDeadTime(0.25),
        correlation_model=AfterpulsingModel(0.1),
    )
    frame = copy(capture!(det, channels; rng=MersenneTwister(31)))
    gated = 2.0 * 0.5
    dead_time_scale = 0.25 / 0.5
    expected_value = (gated / (1 + gated * dead_time_scale)) * 1.1
    expected = fill(expected_value, 4, 4)
    return Dict(
        "mean" => mean(frame),
        "expected_mean" => mean(expected),
        "error_max" => maximum(abs.(frame .- expected)),
        "contract_holds" => isapprox(mean(frame), mean(expected); atol=1e-10, rtol=1e-10),
    )
end

function build_report()
    ccd_emccd = ccd_emccd_case()
    cmos = cmos_structured_case()
    ingaas = ingaas_case()
    hgcdte = hgcdte_case()
    apd = apd_case()
    return Dict(
        "artifact_id" => "DET-VAL-2026-04-01",
        "generated_on" => "2026-04-01",
        "scope" => Dict(
            "families" => [
                "ccd",
                "emccd",
                "cmos",
                "ingaas",
                "hgcdte_avalanche_array",
                "apd",
            ],
            "artifact_kind" => "detector_fixture_validation",
            "rng_policy" => "fixed_seed_mersenne_twister",
        ),
        "cases" => Dict(
            "ccd_emccd_gain" => ccd_emccd,
            "cmos_structured_output" => cmos,
            "ingaas_persistence" => ingaas,
            "hgcdte_avalanche" => hgcdte,
            "apd_counting_chain" => apd,
        ),
        "interpretation" => Dict(
            "ccd_emccd_gain_ok" => ccd_emccd["contract_holds"],
            "cmos_structured_output_ok" => cmos["contract_holds"],
            "ingaas_persistence_ok" => ingaas["contract_holds"],
            "hgcdte_avalanche_ok" => hgcdte["contract_holds"],
            "apd_counting_chain_ok" => apd["contract_holds"],
        ),
    )
end

function update_manifest!(artifact_path::AbstractString)
    mkpath(dirname(MANIFEST))
    manifest = isfile(MANIFEST) ? TOML.parsefile(MANIFEST) : Dict{String,Any}()
    artifacts = get!(manifest, "artifacts", Any[])
    kept = Any[item for item in artifacts if get(item, "id", "") != "DET-VAL-2026-04-01"]
    push!(kept, Dict(
        "purpose" => "detector-family validation artifact for PVP-04",
        "id" => "DET-VAL-2026-04-01",
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
