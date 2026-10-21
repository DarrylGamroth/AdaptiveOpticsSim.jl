# Independent instrument science review — 2026-10-21

## Scope and evidence

Read-only source review under `review-julia-package-design`, focusing on physical units, numerical approximations, cadence and validation oracles. No package environment was instantiated or modified; no Julia tests, GPU work, transport services or hardware experiments were run. A few closed-form calculations were evaluated with Python 3. All claims below are source observations or explicit derivations, not newly measured instrument performance.

Reviewed checkouts: REVOLTClassicSim `603e58ba8ab278369d3b1ac138952a8fc6234db8`; REVOLTCopperSim `1ae815314582f666a4c9de3d87a931c591e7a422`; GPI2Sim `1d9b170a5fe1cae4cbdbb2ac0cd9cd0fc8270929`. Instrument working trees were clean when inspected. AOS development source is the supplied `4bb78685`; the relevant detector/atmosphere/Pyramid paths were also inspected in registered AOS 0.2.1 at `/home/dgamroth/.julia/packages/AdaptiveOpticsSim/uvrvC`. This does not establish which dependency version a future newly instantiated application will resolve. Application policy prohibits AOS development path overrides.

Use package roots `/home/dgamroth/workspaces/codex/REVOLTClassicSim.jl`, `REVOLTCopperSim.jl`, `GPI2Sim.jl`, and `AdaptiveOpticsSim.jl` for relative references below.

## Reconstructed model

All commands are additive OPD coefficients in metres. The controller supplies the corrective sign. No implicit factor of two for mirror reflection is applied. The REVOLT 277-element command order is reconstructed from explicitly identified source-map payloads. GPI2 has 97 woofer + 2304 tweeter + 2 OPD-ramp entries, distinct from a modal reconstructor vector or hardware wiring.

The Gaussian mirror surface is w(x,y) = Σⱼ cⱼ exp(−((x−xⱼ)²+(y−yⱼ)²)/(2σ²)), with σ = pitch / √(−2 log(coupling)). Thus the neighboring-center response is the selected coupling. REVOLT pitch 0.125 and coupling 0.35 produce σ = 0.08626550214; grid and coordinate profiles are implementations of the same provisional model. GPI2 pitch/coupling pairs are (2/10, 0.30) and (2/43, 0.15). Tip/tilt OPD coefficients multiply normalized pupil x/y; an angular tilt θ corresponds to c = θD/2. Source: instrument `src/geometry.jl`, GPI2 `src/graphs.jl:140,163,210`.

AOS atmosphere layers are weighted in OPD amplitude by √fⱼ, not fⱼ (`src/atmosphere/multilayer.jl:264,286`), so independent layers add the intended variance fractions. Phase screens use the specified r₀ reference wavelength and convert to metres using λref/(2π) (`src/atmosphere/kolmogorov.jl:57,81`). Subsequent WFS phase is 2πw/λWFS. All three are pupil-OPD plants: their current graph prescriptions do not perform atmospheric Fresnel propagation/scintillation. With one on-axis source, multilayer altitudes do not establish anisoplanatic validation.

| Model | Pupil grid / aperture | WFS | Model period | Detector exposure |
| --- | --- | --- | --- | --- |
| Classic | 240², 1.22 m, obstruction diameter ratio 0.204918 | 16² SH lenslets, 22² detector pixels each, 800 nm, 0.697 arcsec/pixel | 2 ms / 500 Hz | 1.896 ms |
| Copper | 480², 1.22 m, no obstruction | 30-pixel pupils, 64² frame, 550 nm, radius 5 λ/D, 32 modulation points | 2 ms / 500 Hz | 2 ms |
| GPI2 | 240², 7.77 m, 1.024 m obstruction | 60-pixel pupils, 240² optical frame, 750 nm, radius 3 λ/D, 24 modulation points | 0.5 ms / 2 kHz by default | 0.5 ms |

The exposure is a photon-count multiplier, and the model period advances atmosphere once for each graph step. Wall-clock execution rate is a separate question. Each complete modulation cycle uses the same atmospheric/DM state, not evolving samples over one exposure. The initial published atmosphere is after one advance, not the t=0 screen (AOS `src/algorithm_graphs/native_nodes.jl:385`).

Photometry uses F = F₀,band × 10^(−0.4m), aperture intensity throughput, then μ = rate × QE × exposure. GPI2 throughput is intensity transmission, so 0.47808 × QE 0.9 = 0.430272; the code does not square this throughput or apply QE twice. Its ideal continuous-annulus default electron budget is approximately 46,232 electrons per exposure before sampled-pupil/FOV effects, derived from the local AOS I-band zero point. Copper's magnitude 5.752574989 gives precisely half the flux of magnitude 5. Classic's 800 nm diffraction and R-band magnitude are a monochromatic scenario, not a broadband sensor response.

## Findings

### SCI-001 — GPI2 low-flux EM multiplication has a known mean and distribution bias

**Severity:** medium for scientific use, high if used to establish faint-star sensitivity. **Confidence:** high. **Disposition:** confirmed approximation; already disclosed in model provenance, not a newly discovered undocumented defect.

**Evidence:** GPI2 `src/graphs.jl:237–249` selects excess-noise factor √2 when noise is enabled; AOS 0.2.1 `src/algorithm_graphs/native_nodes.jl:3192` builds the default EM sensor without a selectable multiplication model. AOS `src/detectors/emccd.jl:368–382` implements Y/g = max(N + √N Z, 0) at this excess-noise factor. Conditional on N=1, Python evaluation of Φ(1)+φ(1) gives mean 1.0833154706; P(Y=0 | N=1) = Φ(−1) = 0.1586552539. The ideal high-gain one-carrier exponential has mean 1 and no point mass at zero. This distortion also applies to sparse CIC events. It does not affect Copper's configured gain=1, excess-noise factor=1 branch.

**Impact:** matching read-noise variance and deterministic flux bookkeeping does not validate signal-plus-EM statistics. `GPI2Sim/test/runtests.jl:100` deliberately removes photoelectrons/CIC to test read noise, so it cannot detect this bias. Frame-wide bias depends on the actual spatial flux distribution; no claimed numeric whole-frame bias has been established here.

**Remediation/validation:** retain the current approximation as an explicitly selected throughput profile. Provide a CPU reference profile selecting AOS's existing `ConditionalGammaMultiplication` (already in released 0.2.1, CPU-only), which first requires exposing that choice through the graph node because the current graph API does not pass it through. Validate fixed N=1,2 and Poisson μ=0.01,0.1,1,10 mean/variance/CDF, then compare noisy centroid/normalized Pyramid signals for instrument scenes. Do not substitute a new accelerator approximation solely for speed without the same distribution tests.

### SCI-002 — GPI2 relabels the PWFS frame with twice the normalized spatial sampling

**Severity:** low today. **Confidence:** high. **Disposition:** confirmed metadata defect, no demonstrated current pixel-value defect.

**Evidence:** GPI2 `src/graphs.jl:238` supplies `normalized_pupil_sampling=2f0/60`. AOS 0.2.1 `src/wfs/pyramid/stages.jl:478–487` assigns sampling = output sampling factor / pupil resolution, hence 1/60 for a 60-pixel pupil. AOS graph detector preparation (`src/algorithm_graphs/native_nodes.jl:3180–3190`) discards producer metadata and reconstructs metadata from the supplied value. Both stages claim `NormalizedPupilCoordinates` but differ by two. Copper correctly supplies 1/30.

**Impact:** AOS cell-integrated acquisition uses spatial scale 1 (`src/detectors/detector_acquisition.jl:141`), so the existing frame values and total counts are unaffected. Consumers of the detector optical contract can infer a field size twice as large. Future spatial-density or response operations could turn this into a numerical error.

**Remediation/validation:** use the established detector-plane convention, 1/60; assert identical sampling/domain/spectral metadata across the producer-consumer binding, while confirming unchanged frame values. Keep normalized mirror coordinates −1…+1 separate from AOS detector-plane normalization. Prefer propagating optical-plane metadata across graph bindings over silently rebuilding it when feasible.

### SCI-003 — Full-frame sampling does not establish finite-exposure or servo-lag fidelity

**Severity:** medium validity limit for all three. **Confidence:** high for omission; unquantified scientific impact. **Disposition:** explicit model approximation requiring a convergence study, not a proven unacceptable error.

**Evidence:** AOS `src/algorithm_graphs/native_nodes.jl:385–397` advances once and renders once; the WFS nodes consume that OPD once; detector `step_graph_node!` captures one rate map (`native_nodes.jl:2511` in development; same operation in 0.2.1). `src/detectors/frame_response.jl:204–225` multiplies rate by exposure. No instrument graph repeats propagation inside the exposure. GPI2 cadence comes from `src/graphs.jl:190`; REVOLT periods/exposures from their `src/parameters.jl`.

**Derived scale:** at the 10 m/s REVOLT layers, one model step shifts turbulence by 3.93 Classic pupil pixels / 7.87 Copper pupil pixels. GPI2 default shifts by 0.154 pupil pixels. These numbers are an indication that temporal integration deserves measurement; pixel displacement alone is not an error bound. The missing exposure averaging cannot be repaired by averaging OPD before nonlinear WFS propagation.

**Validation:** hold one common atmospheric history and command timing fixed, integrate detector photon rates from 1,2,4,8 time samples per exposure, then compare frame signals, calibrated interaction responses, residual PSD and closed-loop Strehl. State exposure-center/end timestamps and command application delay explicitly. Separate this study from real-time service latency. AOS already has generic exposure machinery; instrument adoption must preserve total exposure/QE and apply detector read noise once per completed readout.

### SCI-004 — Fast Pyramid propagation and present numerical settings need instrument-specific accuracy gates

**Severity:** medium validation gap for Copper and GPI2. **Confidence:** high that current reviewed package tests do not supply the gate; no asserted numerical failure. **Disposition:** investigation before changing resolution, quadrature or strategy.

**Evidence:** Copper `src/graphs.jl:108–115` and GPI2 `src/graphs.jl:228–231` unconditionally select `PyramidShiftedMaskStrategy`. AOS explicitly calls this approximate (`src/wfs/pyramid/setup.jl:210`). Its parity fixture (`test/testsets/wfs_common_and_parity.jl:176–192,284–288`) tests a small four-pixel-pupil example with one off-grid modulation point, conserving total flux but accepting a nonzero image relative error below 2%. That is useful algorithm evidence; it is not a bound for either instrument's complete modulation ring, aberration range or pupil registration. The instrument tests verify finite output, command response, replay and local calibration, without comparing these physical settings against the pupil-tilt reference or increasing quadrature/resolution.

**Validation:** reference-vs-fast frames and central-difference interaction columns on flat pupil, low-order modes, representative high-order DM probes, and seeded atmosphere. Compare total flux, detector-image norm, normalized Pyramid signals, singular spectra and optical gains. Refine modulation 24→48→96 for GPI2 and 32→64→128 for Copper; refine pupil sampling separately. Stop when error is below a declared application tolerance. GPI2 currently has no public propagation-strategy setting, so a small explicit numerical configuration seam would improve reviewability. Preserve the reference path as an oracle. Do not assume the small-fixture 2% ceiling applies to instruments.

### SCI-005 — Instrument calibration and prediction remain provisional despite self-consistency tests

**Severity:** high limit on claims of physical instrument accuracy, not a code defect. **Confidence:** high. **Disposition:** missing measured inputs/qualification.

All three use analytic Gaussian influence functions and provisional pupil/command registration; none models measured DM temporal response, stroke/nonlinearity or hardware voltage conversion. Copper has symmetric pupil placement rather than the measured ROI registration and an unobscured pupil. Classic uses typical CMOS characterization and an assumed 12-bit mode/full-well ADC mapping, not unit calibration. GPI2 uses a legacy annular pupil without spiders, synthetic single-layer turbulence, homogeneous reduction of eight output statistics, provisional 240→242-row packing, and no measured full-well/charge-diffusion/gain/bias map. These limits are explicitly documented in instrument source comments and GPI2 `docs/model-provenance.md`; they should not be mislabeled software bugs.

Classic/Copper optional pyRTC tests construct and invert a simulation-local push-pull matrix and close the same simulated plant. Their synthetic AOC tests validate ordering/scaling/pseudoinverse construction. These are useful controls/orchestration regressions, not independent measured optical or detector oracles. GPI2's electron bookkeeping test is independent of the gain/QE conversion but uses the simulation's own photon-rate result, so it is not an absolute radiometric oracle.

**Next evidence:** prioritize measured DM influence/sign/order calibration and flat/±probe WFS frames under matched pupil/registration/wavelength; independently measured absolute photon/electron flux and dark/read/EM statistics; then closed-loop residual PSD/Strehl. For simulation calibration itself, check push-pull amplitude convergence because Classic/Copper noiseless calibration still retains ADC quantization and clipping. Classic's inferred input-referred ADC increment is about 1.448 electrons (94,000/(4095×15.8489)); 'noise=false' does not imply an unquantized derivative. No claim of failed calibration is made without amplitude-sweep evidence.

### SCI-006 — Prepared accelerator photon noise uses a declared high-mean Gaussian approximation

**Severity:** medium validation gap when photon-count tails or thresholded measurements matter. **Confidence:** high for the source behavior; scientific impact unmeasured. **Disposition:** known numerical approximation, not a demonstrated instrument accuracy failure.

**Evidence:** registered AOS 0.2.1 `src/backends/random.jl:5–6,60–103` sets `_POISSON_NORMAL_APPROXIMATION_MEAN = 30` and `_POISSON_INVERSION_DRAWS = 96`. The prepared counter-RNG sampler returns bounded product-inversion Poisson draws for 0 < λ < 30 unless its bound is exhausted. For λ ≥ 30, or the bounded fallback, it returns max(0, floor(λ + √λ Z + 0.5)). The source explicitly documents the normal approximation. This applies to the prepared accelerator random path and is separate from SCI-001: it approximates the incident photon/electron count distribution before EM multiplication. No claim is made that all CPU RNG paths use this sampler. The noiseless CPU/GPU parity tests cannot validate it.

**Impact:** matching approximate first moments, reproducible seeds and finite output does not establish Poisson skewness/tail/CDF accuracy or continuity of output statistics across λ=30. No unacceptable bias/tail error for the current instrument scenes has been measured here.

**Validation/remediation:** compare independent exact/reference Poisson CDFs, moments, quantiles and relevant tail/threshold probabilities at λ=29.9,30,30.1 and representative larger means from the actual photon-rate maps. Use sufficient samples and statistical confidence intervals, with analytic PMF/CDF values as the oracle. If the application tolerance requires exact Poisson statistics, qualify an exact/rejection reference and any accelerator implementation, including capture, deterministic stream ownership and bounded-resource behavior. Keep changing the distribution separate from optimizing the present algorithm.

**Performance hypothesis only:** the low-mean source loop continues through 96 iterations using an `active` flag after a draw terminates. It is plausible that generated predication/control flow costs work, but source alone cannot establish actual device instructions or a bottleneck. Inspect compiled code and selected-kernel measurements before proposing early-exit/warp/control-flow changes; preserve per-element key/domain/sub-draw mapping and distribution. No performance conclusion follows merely from the loop's written trip count.

## Accuracy-preserving priorities

1. Fix SCI-002 with a metadata contract test; it does not change the present detector pixels.
2. Establish SCI-004 reference/gate cases before tuning FFT size, quadrature or numerical resolution. Same-model separable Gaussian evaluation and precomputed fixed optical masks are stronger speed candidates than silently reducing scientific sampling.
3. Expose an explicit CPU gamma multiplication reference profile for SCI-001 and keep measured noise/calibration identities attached to profiles.
4. Quantify SCI-003 before making claims about exposure transfer functions or controller bandwidth.
5. Gate the prepared accelerator photon sampler against independent Poisson distributions (SCI-006), especially around its λ=30 branch.
6. Acquire the measured inputs in SCI-005 before promoting simulated regression values to instrument predictions.

No production files changed. No high-severity previously unknown numerical implementation defect was established within this bounded review.
