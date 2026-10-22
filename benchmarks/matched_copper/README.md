# Matched Copper optical benchmark

This compares AOS, pinned current OOPAO and pinned current SPECULA on the
Copper 240/96 optical workload. It measures completed photon-rate production;
detector noise, atmospheric propagation, DM runtime and RTC are excluded.
Slopes and the complete 277-column optical calibration response use a shared
readout on the saved physical photon rates. They do not validate each external
package's native slope or calibration APIs.

## Sources and isolation

`python prepare_sources.py --root /path/to/isolated/sources` fetches the exact
commits in `sources.json` into missing directories and verifies existing ones.
`--verify-only` works offline. It never resets an existing checkout or follows
an unpinned branch. Existing checkouts must be at the pinned revision and
completely clean, including untracked files. Source acquisition and dependency installation are setup;
they do not occur in a timing sample. Existing user checkouts are not modified.
Record package locks, source revisions, dirty state and hashes with each run.

A benchmark should recreate pinned checkouts when needed, then reuse verified
ones. A normal run should not silently clone or update latest upstream.

## Physical contract

`fixture.jl` generates a shared binary pupil, seventeen OPD probes, and the
provisional Copper Gaussian277 OPD basis. Every engine receives those exact
arrays. Parameters are 1.22 m, 550 nm, nominal5lambda/D circular modulation,
96 points, pupil240, pupil-image30, separation2, edge1, FFT512 and output64.
All runtime optical arrays and transforms use Float32/ComplexF32. The common
endpoint tilt convention samples nominal5*240/239lambda/D on a D/N physical
pixel pitch; this is reported rather than silently correcting one engine.

SPECULA's native padding guard, circular stop and prism orientation differ
from this prescription. `specula_adapter.py` explicitly changes these cold
choices while retaining native propagation and lifecycle. Preserve an
unmodified native diagnostic separately. This comparison is of engines under
a common prescription, not proof of identical native defaults.

`protocol.toml` freezes accuracy thresholds and timing boundaries before
numerical results. No fitted gain, flux, registration or response scaling is
allowed. Investigate and retain failed comparisons before accepting timings.

The timed input is resident OPD in the engine's native units: metres in
AOS/OOPAO, nanometres in SPECULA. Input staging is excluded; native phase
formation remains included. SPECULA's native auxiliary PSF work remains part
of its propagation, while OOPAO's optional focal-plane diagnostic is disabled.
Report that output-work difference when interpreting performance.

## Reproduce the captured campaign

The archived fixture is the comparison input. Reuse its `fixture/` directory
unchanged for validation and performance runs; `common.py` checks the saved
array metadata and checksums. Regenerate it only when intentionally creating a
new campaign from the exact AOS and Copper source revisions recorded in
`fixture.toml`.

Use one campaign root with this layout:

```text
ROOT/
  fixture/
  julia-env/Project.toml
  julia-env/Manifest.toml
  sources/oopao/
  sources/specula/
  venv/
```

Prepare and verify the pinned Python sources. The bootstrap creates only
missing checkouts; it does not fetch or alter existing ones when they exist.

```sh
ROOT=/path/to/campaign
python3.13 benchmarks/matched_copper/prepare_sources.py --root "$ROOT/sources"
python3.13 benchmarks/matched_copper/prepare_sources.py --root "$ROOT/sources" --verify-only
```

Create the isolated Python environment from the checked-in full lock. It pins
direct and transitive Python dependencies, including CuPy and the CUDA 13.4
cuBLAS runtime used by SPECULA's CUDA path. The benchmark commands set CPU
thread limits and affinity; no system-wide settings are needed.

```sh
python3.13 -m venv "$ROOT/venv"
"$ROOT/venv/bin/python" -m pip install -r benchmarks/matched_copper/requirements-frozen.txt
"$ROOT/venv/bin/python" -m pip check
```

The captured environment used Python 3.13.5, NumPy 2.2.4, SciPy 1.18.1,
CuPy 14.2.0 and CUDA runtime/NVRTC packages 13.4.92. CuPy links CUDA
13.2 and discovers the locally installed CUDA 13.4; Julia uses CUDA 13.4. Both
GPU paths use cuFFT 12.4.0. Full Python package versions,
NumPy and SciPy BLAS build information, runtime FFT dtypes and the final
dependency check are in the campaign's `env-info.txt`.

The Julia environment was captured with Julia 1.13.2. Its Manifest pins
AdaptiveOpticsCalibration 0.17.0, CUDA 6.4.2, FixedSizeArrays 1.4.2,
KernelAbstractions 0.9.45, AdaptiveOpticsSim 0.2.2 and REVOLTCopperSim 0.2.0.
The captured AOS and Copper source commits are recorded in `fixture.toml`; the
Python repositories are pinned in `sources.json`. AOS and Copper are local path sources. Prepare separate clean library
checkouts at the commits in `fixture.toml` (AOS `b590aa26`, Copper `781f41d3`).
The benchmark harness lives in the later results commit; using that commit as
the library source against this archived fixture intentionally fails the exact
revision check. A new source revision needs a new fixture campaign. The archived `Project.toml` and
`Manifest.toml` must travel together. Their local source paths are absolute;
when moving the campaign, update the `[sources]` entries in `Project.toml` and
the matching `path` entries in `Manifest.toml` to the new checkout locations,
then instantiate without updating or re-resolving package versions. Keep the
manifest's other package versions unchanged.

Run the complete numerical validation before measuring service time. The default
configuration list includes CPU and CUDA runs for all three engines, and the
collector expects that complete set.

```sh
"$ROOT/venv/bin/python" benchmarks/matched_copper/campaign.py --root "$ROOT" --mode validate
"$ROOT/venv/bin/python" benchmarks/matched_copper/campaign.py --root "$ROOT" --mode benchmark
"$ROOT/venv/bin/python" benchmarks/matched_copper/collect.py --root "$ROOT" --output "$ROOT/collected"
```

Validation and performance outputs are separate per-configuration directories
under `ROOT`; campaign logs and a JSONL execution journal are also retained.
The benchmark mode uses the warmups, sample count and round count frozen in
`protocol.toml`. Interpret timing only alongside the completed numerical checks. Failed
accuracy gates remain failures; timing collection does not turn them into a
qualification pass.

The source OPD is staged in native units before each sample: metres for AOS and
OOPAO, nanometres for SPECULA. Repository updates, input staging/unit
conversion and diagnostic copies are outside the timed section. Each native
OPD setter, phase formation and optical propagation remain inside the measured
step when required by that engine. The shared post-readout slope and
calibration evaluator is a common comparison tool. It does not qualify
OOPAO's or SPECULA's native slope estimators, native DM implementations, or a
closed-loop controller.

The captured evidence and finalized host-specific interpretation will be
linked from the [campaign results guide](../results/platform/2026-10-22-matched-copper/README.md).
