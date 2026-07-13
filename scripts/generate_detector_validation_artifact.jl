using AdaptiveOpticsSim
using Random
using Statistics
using TOML

import AdaptiveOpticsSim: AveragedNonDestructiveReads, FrameWindow,
    StaticCMOSOutputPattern, detector_combined_frame, detector_read_cube,
    detector_read_times, detector_reference_cube, detector_reference_frame,
    detector_signal_cube, detector_signal_frame, readout_products

const OUTDIR = joinpath(@__DIR__, "..", "benchmarks", "results", "detectors")
const OUTFILE = joinpath(OUTDIR, "2026-07-12-detector-mkid-validation.toml")
const MANIFEST = joinpath(OUTDIR, "manifest.toml")
const ARTIFACT_ID = "DET-VAL-2026-07-12"

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

function hgcdte_avalanche_case()
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

function hgcdte_multiread_case()
    zero_psf = zeros(64, 64)
    raw_sigma = 4.0
    det_ndr = Detector(
        integration_time=1.0,
        noise=NoiseReadout(raw_sigma),
        qe=1.0,
        binning=1,
        gain=1.0,
        response_model=NullFrameResponse(),
        sensor=HgCdTeAvalancheArraySensor(read_time=1.0,
            sampling_mode=AveragedNonDestructiveReads(4)),
    )
    frame_ndr = copy(capture!(det_ndr, zero_psf; rng=MersenneTwister(41)))
    det_cds = Detector(
        integration_time=1.0,
        noise=NoiseReadout(raw_sigma),
        qe=1.0,
        binning=1,
        gain=1.0,
        response_model=NullFrameResponse(),
        sensor=HgCdTeAvalancheArraySensor(read_time=1.0,
            sampling_mode=CorrelatedDoubleSampling()),
    )
    frame_cds = copy(capture!(det_cds, zero_psf; rng=MersenneTwister(41)))
    det_fowler = Detector(
        integration_time=1.0,
        noise=NoiseReadout(raw_sigma),
        qe=1.0,
        binning=1,
        gain=1.0,
        response_model=NullFrameResponse(),
        sensor=HgCdTeAvalancheArraySensor(read_time=1.0,
            sampling_mode=FowlerSampling(8)),
    )
    frame_fowler = copy(capture!(det_fowler, zero_psf; rng=MersenneTwister(41)))

    ndr_products = readout_products(det_ndr)
    cds_products = readout_products(det_cds)
    fowler_products = readout_products(det_fowler)
    ndr_expected_std = raw_sigma / sqrt(4)
    cds_expected_std = raw_sigma * sqrt(2)
    fowler_expected_std = raw_sigma * sqrt(2 / 8)

    return Dict(
        "ndr_std" => std(vec(frame_ndr)),
        "cds_std" => std(vec(frame_cds)),
        "fowler_std" => std(vec(frame_fowler)),
        "ndr_expected_std" => ndr_expected_std,
        "cds_expected_std" => cds_expected_std,
        "fowler_expected_std" => fowler_expected_std,
        "ndr_read_cube_reads" => size(detector_read_cube(ndr_products), 3),
        "cds_reference_cube_reads" => size(detector_reference_cube(cds_products), 3),
        "fowler_signal_cube_reads" => size(detector_signal_cube(fowler_products), 3),
        "ndr_read_times_last" => last(detector_read_times(ndr_products)),
        "combined_consistent_cds" =>
            isapprox(maximum(abs, detector_combined_frame(cds_products) .-
                (detector_signal_frame(cds_products) .- detector_reference_frame(cds_products))), 0.0; atol=1e-10),
        "combined_consistent_fowler" =>
            isapprox(maximum(abs, detector_combined_frame(fowler_products) .-
                (detector_signal_frame(fowler_products) .- detector_reference_frame(fowler_products))), 0.0; atol=1e-10),
        "contract_holds" => std(vec(frame_ndr)) < std(vec(frame_cds)) &&
            std(vec(frame_fowler)) < std(vec(frame_cds)) &&
            isapprox(std(vec(frame_ndr)), ndr_expected_std; rtol=0.15) &&
            isapprox(std(vec(frame_cds)), cds_expected_std; rtol=0.15) &&
            isapprox(std(vec(frame_fowler)), fowler_expected_std; rtol=0.15) &&
            size(detector_read_cube(ndr_products), 3) == 4 &&
            size(detector_reference_cube(cds_products), 3) == 1 &&
            size(detector_signal_cube(fowler_products), 3) == 8 &&
            last(detector_read_times(ndr_products)) == 4.0 &&
            maximum(abs, detector_combined_frame(cds_products) .-
                (detector_signal_frame(cds_products) .- detector_reference_frame(cds_products))) ≤ 1e-10 &&
            maximum(abs, detector_combined_frame(fowler_products) .-
                (detector_signal_frame(fowler_products) .- detector_reference_frame(fowler_products))) ≤ 1e-10,
    )
end

function hgcdte_timing_and_correction_case()
    zero_psf = zeros(4, 4)
    det_timed_single = Detector(
        integration_time=1.0,
        noise=NoiseNone(),
        qe=1.0,
        binning=1,
        dark_current=3.0,
        gain=1.0,
        response_model=NullFrameResponse(),
        sensor=HgCdTeAvalancheArraySensor(glow_rate=2.0, read_time=1.0),
    )
    frame_single = copy(capture!(det_timed_single, zero_psf; rng=MersenneTwister(42)))
    det_timed_cds = Detector(
        integration_time=1.0,
        noise=NoiseNone(),
        qe=1.0,
        binning=1,
        dark_current=3.0,
        gain=1.0,
        response_model=NullFrameResponse(),
        sensor=HgCdTeAvalancheArraySensor(glow_rate=2.0, read_time=1.0, sampling_mode=CorrelatedDoubleSampling()),
    )
    frame_cds = copy(capture!(det_timed_cds, zero_psf; rng=MersenneTwister(42)))
    timed_meta = detector_export_metadata(det_timed_cds)

    det_windowed = Detector(
        integration_time=1.0,
        noise=NoiseNone(),
        qe=1.0,
        binning=1,
        gain=1.0,
        response_model=NullFrameResponse(),
        sensor=HgCdTeAvalancheArraySensor(read_time=1.0, sampling_mode=CorrelatedDoubleSampling()),
        readout_window=FrameWindow(2:3, 2:3),
    )
    capture!(det_windowed, fill(10.0, 4, 4); rng=MersenneTwister(43))
    windowed_meta = detector_export_metadata(det_windowed)

    composite_pattern = repeat(reshape([1.0, 2.0, 3.0, 4.0], :, 1), 1, 4) .+
        repeat(reshape([1.0, 2.0, 3.0, 4.0], 1, :), 4, 1)
    det_corrected = Detector(
        integration_time=1.0,
        noise=NoiseNone(),
        qe=1.0,
        binning=1,
        gain=1.0,
        response_model=NullFrameResponse(),
        sensor=HgCdTeAvalancheArraySensor(sampling_mode=CorrelatedDoubleSampling()),
        correction_model=CompositeFrameReadoutCorrection((
            ReferenceRowCommonModeCorrection(1),
            ReferenceColumnCommonModeCorrection(1),
        )),
    )
    corrected = capture!(det_corrected, composite_pattern; rng=MersenneTwister(44))
    corrected_products = readout_products(det_corrected)

    return Dict(
        "single_sum" => sum(frame_single),
        "cds_sum" => sum(frame_cds),
        "sampling_read_time" => timed_meta.sampling_read_time,
        "sampling_wallclock_time" => timed_meta.sampling_wallclock_time,
        "windowed_sampling_read_time" => windowed_meta.sampling_read_time,
        "windowed_sampling_wallclock_time" => windowed_meta.sampling_wallclock_time,
        "windowed_read_times" => detector_read_times(det_windowed),
        "correction_residual_max" => maximum(abs, corrected),
        "combined_consistent_corrected" =>
            maximum(abs, detector_combined_frame(corrected_products) .-
                (detector_signal_frame(corrected_products) .- detector_reference_frame(corrected_products))),
        "contract_holds" => sum(frame_cds) > sum(frame_single) &&
            timed_meta.sampling_read_time == 1.0 &&
            timed_meta.sampling_wallclock_time == 3.0 &&
            windowed_meta.sampling_read_time == 0.5 &&
            windowed_meta.sampling_wallclock_time == 2.0 &&
            detector_read_times(det_windowed) == [0.5, 1.0] &&
            maximum(abs, corrected) ≤ 1e-6 &&
            maximum(abs, detector_combined_frame(corrected_products) .-
                (detector_signal_frame(corrected_products) .- detector_reference_frame(corrected_products))) ≤ 1e-10,
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

function spad_case()
    det = SPADArrayDetector(
        integration_time=1.0,
        noise=NoiseNone(),
        sensor=SPADArraySensor(pde=0.5, dark_count_rate=0.0, fill_factor=0.8),
    )
    frame = copy(capture!(det, fill(10.0, 8, 8), MersenneTwister(32)))
    metadata = detector_export_metadata(det)
    return Dict(
        "mean" => mean(frame),
        "expected_mean" => 4.0,
        "sensor" => String(metadata.sensor),
        "contract_holds" => frame == fill(4.0, 8, 8) &&
            metadata.sensor == :spad_array && metadata.fill_factor == 0.8,
    )
end

function mkid_case()
    det = MKIDArrayDetector(
        integration_time=2.0,
        noise=NoiseNone(),
        sensor=MKIDArraySensor(
            qe=0.75,
            dark_count_rate=0.0,
            fill_factor=0.8,
            energy_resolution=12.0,
            timing_jitter_s=2e-6,
            wavelength_range_m=(0.8e-6, 1.4e-6),
        ),
    )
    input = fill(10.0, 8, 8)
    matrix_frame = copy(capture!(det, input, MersenneTwister(33)))
    inside = Source(band=:custom, wavelength=1.0e-6)
    outside = Source(band=:custom, wavelength=0.55e-6)
    inside_frame = copy(capture!(det, input, inside, MersenneTwister(33)))
    outside_frame = copy(capture!(det, input, outside, MersenneTwister(33)))
    spectral = with_spectrum(
        inside,
        SpectralBundle([0.6e-6, 1.0e-6, 1.6e-6], [0.2, 0.3, 0.5]),
    )
    spectral_frame = copy(capture!(det, input, spectral, MersenneTwister(33)))
    metadata = detector_export_metadata(det)
    return Dict(
        "matrix_mean" => mean(matrix_frame),
        "inside_band_mean" => mean(inside_frame),
        "outside_band_max" => maximum(abs, outside_frame),
        "spectral_bundle_mean" => mean(spectral_frame),
        "energy_resolution" => metadata.energy_resolution,
        "timing_jitter_s" => metadata.timing_jitter_s,
        "wavelength_min_m" => metadata.wavelength_min_m,
        "wavelength_max_m" => metadata.wavelength_max_m,
        "contract_holds" => isapprox(matrix_frame, fill(12.0, 8, 8); atol=1e-12, rtol=1e-12) &&
            isapprox(inside_frame, matrix_frame; atol=1e-12, rtol=1e-12) && all(iszero, outside_frame) &&
            isapprox(mean(spectral_frame), 3.6; atol=1e-12, rtol=1e-12) &&
            metadata.energy_resolution == 12.0 && metadata.timing_jitter_s == 2e-6 &&
            metadata.wavelength_min_m == 0.8e-6 && metadata.wavelength_max_m == 1.4e-6,
    )
end

function build_report()
    ccd_emccd = ccd_emccd_case()
    cmos = cmos_structured_case()
    ingaas = ingaas_case()
    hgcdte_avalanche = hgcdte_avalanche_case()
    hgcdte_multiread = hgcdte_multiread_case()
    hgcdte_timing = hgcdte_timing_and_correction_case()
    apd = apd_case()
    spad = spad_case()
    mkid = mkid_case()
    return Dict(
        "artifact_id" => ARTIFACT_ID,
        "generated_on" => "2026-07-12",
        "scope" => Dict(
            "families" => [
                "ccd",
                "emccd",
                "cmos",
                "ingaas",
                "hgcdte_avalanche_array",
                "apd",
                "spad_array",
                "mkid_array",
            ],
            "artifact_kind" => "detector_fixture_validation",
            "rng_policy" => "fixed_seed_mersenne_twister",
        ),
        "cases" => Dict(
            "ccd_emccd_gain" => ccd_emccd,
            "cmos_structured_output" => cmos,
            "ingaas_persistence" => ingaas,
            "hgcdte_avalanche" => hgcdte_avalanche,
            "hgcdte_multiread" => hgcdte_multiread,
            "hgcdte_timing_and_correction" => hgcdte_timing,
            "apd_counting_chain" => apd,
            "spad_counting_array" => spad,
            "mkid_counting_array" => mkid,
        ),
        "interpretation" => Dict(
            "ccd_emccd_gain_ok" => ccd_emccd["contract_holds"],
            "cmos_structured_output_ok" => cmos["contract_holds"],
            "ingaas_persistence_ok" => ingaas["contract_holds"],
            "hgcdte_avalanche_ok" => hgcdte_avalanche["contract_holds"],
            "hgcdte_multiread_ok" => hgcdte_multiread["contract_holds"],
            "hgcdte_timing_and_correction_ok" => hgcdte_timing["contract_holds"],
            "apd_counting_chain_ok" => apd["contract_holds"],
            "spad_counting_array_ok" => spad["contract_holds"],
            "mkid_counting_array_ok" => mkid["contract_holds"],
        ),
    )
end

function validate_report(report)
    failures = sort!([name for (name, case) in report["cases"] if !case["contract_holds"]])
    isempty(failures) || throw(NumericalConditionError(
        "detector validation contracts failed: $(join(failures, ", "))"))
    return report
end

function update_manifest!(artifact_path::AbstractString)
    mkpath(dirname(MANIFEST))
    manifest = isfile(MANIFEST) ? TOML.parsefile(MANIFEST) : Dict{String,Any}()
    artifacts = get!(manifest, "artifacts", Any[])
    for item in artifacts
        if get(item, "id", "") == "DET-VAL-2026-04-23"
            item["purpose"] = "superseded detector validation; retained as the historical failed HgCdTe multi-read result"
            item["status"] = "superseded"
            item["superseded_by"] = ARTIFACT_ID
        end
    end
    kept = Any[item for item in artifacts if get(item, "id", "") != ARTIFACT_ID]
    push!(kept, Dict(
        "purpose" => "detector-family validation with HgCdTe multi-read and MKID counting coverage",
        "id" => ARTIFACT_ID,
        "path" => basename(artifact_path),
        "status" => "active",
    ))
    manifest["artifacts"] = kept
    open(MANIFEST, "w") do io
        TOML.print(io, manifest)
    end
end

function main()
    mkpath(OUTDIR)
    report = validate_report(build_report())
    open(OUTFILE, "w") do io
        TOML.print(io, report)
    end
    update_manifest!(OUTFILE)
    println(OUTFILE)
end

main()
