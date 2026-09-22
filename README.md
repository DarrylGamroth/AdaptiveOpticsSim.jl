# AdaptiveOpticsSim.jl

[![CPU Validation](https://github.com/DarrylGamroth/AdaptiveOpticsSim.jl/actions/workflows/cpu-validation.yml/badge.svg)](https://github.com/DarrylGamroth/AdaptiveOpticsSim.jl/actions/workflows/cpu-validation.yml)
[![Coverage](https://codecov.io/gh/DarrylGamroth/AdaptiveOpticsSim.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/DarrylGamroth/AdaptiveOpticsSim.jl)

AdaptiveOpticsSim.jl (AOS) is an independent Julia adaptive-optics plant
simulator for external-RTC HIL development, deterministic validation, and
offline CPU/GPU studies. The maintained design uses idiomatic Julia dispatch,
explicit optical products, prepared workspaces, and backend-portable
algorithms.

Requires Julia 1.12 or newer. The package relies on current Julia atomics and
backend behavior for maintained CPU/GPU execution paths.

## Start Here

If you are a normal user, read these in order:

- [docs/user-guide.md](docs/user-guide.md)
- [Use AOS with pyRTC](docs/user-guide.md#use-aos-with-pyrtc) for the external-
  RTC reference workflow
- [docs/model-cookbook.md](docs/model-cookbook.md)
- [docs/api-reference.md](docs/api-reference.md)
- `examples/tutorials/`

If you are working on validation, production support, or backend work, then use:

- [docs/supported-production-surfaces.md](docs/supported-production-surfaces.md)
- [docs/release-validation-runbook.md](docs/release-validation-runbook.md)
- [docs/documentation-map.md](docs/documentation-map.md)

## First Model

### 1. Optics-only direct image

```julia
using AdaptiveOpticsSim
using AdaptiveOpticsSim.Optics
using AdaptiveOpticsSim.WavefrontSensors

tel = Telescope(resolution=32, diameter=8.0, central_obstruction=0.1)
src = Source(band=:I, magnitude=8.0)
pupil = PupilFunction(tel)
imaging = prepare_direct_imaging(pupil, src; zero_padding=2)
form_direct_image!(imaging)
photon_rate_image = intensity_values(direct_imaging_output(imaging))
```

The caller-owned result is a source-brightness-scaled, cell-integrated photon-
arrival rate on focal-plane angular coordinates before detector exposure; it
is not inherently a unit-normalized PSF.

### 2. Add atmosphere and a wavefront sensor

```julia
atm = MultiLayerAtmosphere(
    tel;
    r0=0.15,
    reference_wavelength_m=500e-9,
    L0=25.0,
    fractional_cn2=[0.6, 0.4],
    wind_speed=[8.0, 12.0],
    wind_direction_deg=[0.0, 90.0],
    altitude=[0.0, 5000.0],
)

wfs = ShackHartmannWFS(
    tel;
    n_lenslets=4,
    pixel_scale_arcsec=0.1,
    n_pix_subap=6,
)

rng = runtime_rng(0)
renderer = prepare_atmosphere_renderer(atm, tel, src)
pupil = PupilFunction(tel)
epoch = advance_by!(atm, 1e-3; rng=rng)
render_atmosphere!(pupil, renderer, atm, epoch)
photon_rate = shack_hartmann_rate_map(wfs, pupil, src)
optics = shack_hartmann_optics(wfs, src)
prepared_optics = prepare_wfs_optics(optics, pupil, photon_rate)
form_wfs_optical_products!(photon_rate, pupil, prepared_optics)
```

The resulting mosaic is the physical detector-plane photon rate. A detector
acquisition produces the complete frame consumed by an external RTC. AOS does
not convert Shack–Hartmann frames to operational slopes; the maintained FGA/JFG
composition below owns that estimation contract.

### 3. Compose an AO plant with an RTC

AOS owns the physical plant: atmosphere, optics, WFS products, detector
acquisition, and application of complete PDM commands. The maintained
in-process RTC reference is
[`examples/integrations/filter_graph_algorithms/`](examples/integrations/filter_graph_algorithms/):
FilterGraphAlgorithms/JuliaFilterGraph owns reconstruction, controller state,
frame delay, and VDM/PDM routing. The fixture contains maintained
Shack–Hartmann and complete-frame Pyramid plant/RTC references.

The legacy AOS closed-loop scripts and tutorials were retired. Bi-O-edge RTC
composition belongs in a downstream package with its calibration and acceptance
evidence: AOS supplies the physical optics and complete detector acquisition,
while registered FilterGraphAlgorithms v0.3.0 and JuliaFilterGraph v0.2.1 own
signal estimation and RTC processing. AOS retains the Zernike and Curvature
estimators until their corresponding FGA targets are released and adopted. AOS
retains the LiFT forward model and retains its inverse workflow
until the AdaptiveOpticsCalibration target is released and adopted. The
combined Subaru AO188/AO3k model likewise moves to a downstream instrument
package.

The main modeling objects are:

- `Telescope` and `Source` for optical geometry and illumination
- `MultiLayerAtmosphere` or `KolmogorovAtmosphere` for turbulence
- `ShackHartmannWFS`, `PyramidWFS`, `BiOEdgeWFS`, `CurvatureWFS`, `ZernikeWFS` for sensing
- `DeformableMirror` for the physical response to a complete PDM command
- `AdaptiveOpticsSim.AlgorithmGraphs` for static, single-rate, complete-frame
  composition on CPU, CUDA, or AMDGPU arrays
- direct Julia composition for generated topology, multiple rates, conditional
  execution, or sub-frame optical sampling

For a lockstep external RTC, `prepare_graph_hil_boundary` binds one graph input
as the complete DM command and one graph output as the complete detector frame.
Transport and wall-clock pacing remain outside the simulation graph.

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

Runnable tutorials live in `examples/tutorials/`. Start with:

```bash
julia --project=. examples/tutorials/image_formation.jl
julia --project=. examples/tutorials/detector.jl
```

For the maintained Shack–Hartmann, Pyramid, and Bi-O-edge plant/RTC references,
use `examples/integrations/filter_graph_algorithms/` with its documented
cross-package environment. Bi-O-edge and Zernike RTC composition is not an AOS
tutorial surface. Curvature RTC composition is
likewise downstream; LiFT remains a distinct AOS forward and phase-retrieval
tutorial while its calibration-package inverse target is pending.
To verify the maintained core examples as a group, run:

```bash
./scripts/run_core_examples.sh
```

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

rng = runtime_rng(0)
observation = WFSObservation(
    similar(intensity_values(photon_rate), UInt16);
    units=:adu,
    layout=:lenslet_mosaic,
)
acquisition = prepare_wfs_acquisition(
    det,
    photon_rate,
    observation;
    source=src,
)
acquire_wfs_observation!(observation, photon_rate, acquisition, rng)
adu = observation_storage(observation)
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
- [docs/release-validation-runbook.md](docs/release-validation-runbook.md)
- [docs/backend-validation-guide.md](docs/backend-validation-guide.md)

For deeper developer reference:

- [docs/documentation-map.md](docs/documentation-map.md)

For maintainers extending subsystem families:

- [docs/extension-guide.md](docs/extension-guide.md)

For full plotted examples, use the companion plotting package in
`../AdaptiveOpticsSimPlots.jl`. The examples in this repo remain plotting-free
by design. In the plotting package, run visual examples with:

```bash
GKSwstype=100 julia --project=. examples/image_formation_visual.jl
GKSwstype=100 julia --project=. examples/wfs_detector_comparison_visual.jl
```

## Supported Production Surface

The package has a documented maintained production-validation surface. Use:

- [docs/supported-production-surfaces.md](docs/supported-production-surfaces.md)
- [docs/release-validation-runbook.md](docs/release-validation-runbook.md)

Release validation entry point:

```bash
./scripts/run_release_validation.sh
```

Optional release-validation tracks are enabled with environment flags, for
example:

```bash
ADAPTIVEOPTICS_VALIDATE_EXAMPLES=1 ./scripts/run_release_validation.sh
ADAPTIVEOPTICS_VALIDATE_AMDGPU=1 ./scripts/run_release_validation.sh
ADAPTIVEOPTICS_VALIDATE_CUDA=1 ./scripts/run_release_validation.sh
```

AMDGPU and CUDA both have current hardware validation through their dedicated
targets. CUDA is the primary GPU performance-optimization target; AMDGPU is the
secondary portability and qualification target. AMDGPU remains the automated
release gate because a continuously available CUDA CI runner has not been
established, while CUDA is validated manually on the WSL RTX host. Both targets
retain the same correctness and steady-state allocation requirements.
