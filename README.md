# AdaptiveOpticsSim.jl

[![CPU Validation](https://github.com/DarrylGamroth/AdaptiveOpticsSim.jl/actions/workflows/cpu-validation.yml/badge.svg)](https://github.com/DarrylGamroth/AdaptiveOpticsSim.jl/actions/workflows/cpu-validation.yml)
[![Coverage](https://codecov.io/gh/DarrylGamroth/AdaptiveOpticsSim.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/DarrylGamroth/AdaptiveOpticsSim.jl)

Julia adaptive optics simulation toolkit. The package keeps the OOPAO modeling
surface recognizable, but uses idiomatic Julia design for performance,
reproducibility, and backend portability.

## Start Here

If you are a normal user, read these in order:

- [docs/user-guide.md](docs/user-guide.md)
- [docs/model-cookbook.md](docs/model-cookbook.md)
- [docs/api-reference.md](docs/api-reference.md)
- `examples/tutorials/`

If you are working on validation, production support, or backend work, then use:

- [docs/supported-production-surfaces.md](docs/supported-production-surfaces.md)
- [docs/release-validation-runbook.md](docs/release-validation-runbook.md)
- [docs/documentation-map.md](docs/documentation-map.md)

## First Model

### 1. Optics-only PSF

```julia
using AdaptiveOpticsSim

tel = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3, central_obstruction=0.1)
src = Source(band=:I, magnitude=8.0)
psf = compute_psf!(tel, src; zero_padding=2)
```

### 2. Add atmosphere and a wavefront sensor

```julia
atm = MultiLayerAtmosphere(
    tel;
    r0=0.15,
    L0=25.0,
    fractional_cn2=(0.6, 0.4),
    wind_speed=(8.0, 12.0),
    wind_direction=(0.0, 90.0),
    altitude=(0.0, 5000.0),
)

wfs = ShackHartmann(tel; n_lenslets=4, mode=Diffractive(), pixel_scale=0.1, n_pix_subap=6)

advance!(atm, tel)
propagate!(atm, tel)
slopes = measure!(wfs, tel, src)
```

### 3. Build a closed-loop AO model

```julia
using Random

rng = runtime_rng(0)
dm = DeformableMirror(tel; n_act=4, influence_width=0.3)
sim = AOSimulation(tel, src, atm, dm, wfs)

imat = interaction_matrix(dm, wfs, tel, src; amplitude=0.1)
recon = ModalReconstructor(imat; gain=0.5)

branch = RuntimeBranch(:main, sim, recon; rng=rng)
cfg = SingleRuntimeConfig(name=:demo, branch_label=:main)
scenario = build_runtime_scenario(cfg, branch)
prepare!(scenario)

for _ in 1:5
    step!(scenario)
end

rt = readout(scenario)
cmd = command(rt)
frame = wfs_frame(rt)
```

`prepare!(...)` performs one-time runtime/WFS precomputation. `step!(...)`
runs the full closed-loop update: sense + reconstruct + apply.

The main modeling objects are:

- `Telescope` and `Source` for optical geometry and illumination
- `MultiLayerAtmosphere` or `KolmogorovAtmosphere` for turbulence
- `ShackHartmann`, `PyramidWFS`, `BioEdgeWFS`, `CurvatureWFS`, `ZernikeWFS` for sensing
- `DeformableMirror` plus a reconstructor for control
- `RuntimeScenario` when you want the maintained step-wise AO simulation surface

For external-control / HIL paths, use `NullReconstructor()` plus
`set_command!(scenario, cmd)` and `sense!(scenario)`.

For advanced controllable-optic and DM modeling, see:

- `ModalControllableOptic(...)` with basis specs such as
  `CartesianTiltBasis(...)` and `ZernikeOpticBasis(...)`
- `DeformableMirror(...; mechanical_coupling=...)`
- `DeformableMirror(...; influence_model=...)` for explicit DM influence models

The user-facing details for those surfaces live in:

- [docs/model-cookbook.md](docs/model-cookbook.md)
- [docs/user-guide.md](docs/user-guide.md)
- [docs/api-reference.md](docs/api-reference.md)

## Tutorials

Runnable tutorial ports live in `examples/tutorials/`. Start with:

```bash
julia --project examples/tutorials/image_formation.jl
julia --project examples/tutorials/detector.jl
julia --project examples/tutorials/closed_loop_shack_hartmann.jl
julia --project examples/tutorials/closed_loop_pyramid.jl
```

Use `examples/closed_loop_demo.jl` for the smallest direct closed-loop script.

## Detector ADU Output

Detector quantization is controlled by `bits` and `full_well`. The Julia array
element type returned to a HIL/RTC boundary is controlled separately by
`output_type`:

```julia
det = Detector(
    noise=NoiseNone(),
    full_well=30_000.0,
    bits=12,
    output_type=UInt16,
)

measure!(wfs, tel, src, det; rng=rng)
adu = wfs_detector_image(wfs, det)
```

In this example, `adu` is a `UInt16` detector image with 12-bit quantized
values. Use `output_type=nothing` when you want the floating-point internal
readout instead of a typed digital export.

For stored references and regression fixtures, prefer
`deterministic_reference_rng(seed)` to preserve the historical
`MersenneTwister` stream. For new RTC/HIL runtime simulations and benchmarks,
prefer `runtime_rng(seed)`, which uses `Xoshiro`.

## Documentation

For users:

- [docs/user-guide.md](docs/user-guide.md)
- [docs/model-cookbook.md](docs/model-cookbook.md)
- [docs/api-reference.md](docs/api-reference.md)
- [docs/julia-tutorial-mappings.md](docs/julia-tutorial-mappings.md)

For validation and supported scope:

- [docs/supported-production-surfaces.md](docs/supported-production-surfaces.md)
- [docs/production-readiness-checklist.md](docs/production-readiness-checklist.md)
- [docs/release-validation-runbook.md](docs/release-validation-runbook.md)

For deeper developer reference:

- [docs/documentation-map.md](docs/documentation-map.md)

For maintainers extending subsystem families:

- [docs/adding-detectors.md](docs/adding-detectors.md)
- [docs/adding-deformable-mirrors.md](docs/adding-deformable-mirrors.md)
- [docs/adding-wavefront-sensors.md](docs/adding-wavefront-sensors.md)
- [docs/adding-controllable-optics.md](docs/adding-controllable-optics.md)
- [docs/adding-controllers-and-reconstructors.md](docs/adding-controllers-and-reconstructors.md)

For full plotted examples, use the companion plotting package in
`../AdaptiveOpticsSimPlots.jl`. The examples in this repo remain plotting-free
by design.

## Supported Production Surface

The package has a documented maintained production-validation surface. Use:

- [docs/supported-production-surfaces.md](docs/supported-production-surfaces.md)
- [docs/production-readiness-checklist.md](docs/production-readiness-checklist.md)
- [docs/release-validation-runbook.md](docs/release-validation-runbook.md)

Release validation entry point:

```bash
./scripts/run_release_validation.sh
```
