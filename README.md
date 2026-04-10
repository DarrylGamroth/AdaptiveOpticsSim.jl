# AdaptiveOpticsSim.jl

Julia adaptive optics simulation toolkit. The package keeps the OOPAO modeling
surface recognizable, but uses idiomatic Julia design for performance,
reproducibility, and backend portability.

## Start Here

If you are a normal user, read these in order:

- [docs/user-guide.md](docs/user-guide.md)
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

wfs = ShackHartmann(tel; n_subap=4, mode=Diffractive(), pixel_scale=0.1, n_pix_subap=6)

advance!(atm, tel)
propagate!(atm, tel)
slopes = measure!(wfs, tel, src)
```

### 3. Build a closed-loop AO model

```julia
using Random

rng = MersenneTwister(0)
dm = DeformableMirror(tel; n_act=4, influence_width=0.3)
sim = AdaptiveOpticsSim.AOSimulation(tel, atm, src, dm, wfs)

imat = interaction_matrix(dm, wfs, tel, src; amplitude=0.1)
recon = ModalReconstructor(imat; gain=0.5)
runtime = ClosedLoopRuntime(sim, recon; rng=rng)
interface = simulation_interface(runtime)

for _ in 1:5
    step!(interface)
end

cmd = simulation_command(interface)
frame = simulation_wfs_frame(interface)
```

The main modeling objects are:

- `Telescope` and `Source` for optical geometry and illumination
- `MultiLayerAtmosphere` or `KolmogorovAtmosphere` for turbulence
- `ShackHartmann`, `PyramidWFS`, `BioEdgeWFS`, `CurvatureWFS`, `ZernikeWFS` for sensing
- `DeformableMirror` plus a reconstructor for control
- `ClosedLoopRuntime` when you want a step-wise AO simulation surface

## Tutorials

Runnable tutorial ports live in `examples/tutorials/`. Start with:

```bash
julia --project examples/tutorials/image_formation.jl
julia --project examples/tutorials/detector.jl
julia --project examples/tutorials/closed_loop_shack_hartmann.jl
julia --project examples/tutorials/closed_loop_pyramid.jl
```

Use `examples/closed_loop_demo.jl` for the smallest direct closed-loop script.

## Documentation

For users:

- [docs/user-guide.md](docs/user-guide.md)
- [docs/api-reference.md](docs/api-reference.md)
- [docs/julia-tutorial-mappings.md](docs/julia-tutorial-mappings.md)

For validation and supported scope:

- [docs/supported-production-surfaces.md](docs/supported-production-surfaces.md)
- [docs/production-readiness-checklist.md](docs/production-readiness-checklist.md)
- [docs/release-validation-runbook.md](docs/release-validation-runbook.md)

For deeper developer reference:

- [docs/documentation-map.md](docs/documentation-map.md)

## Supported Production Surface

The package has a documented maintained production-validation surface. Use:

- [docs/supported-production-surfaces.md](docs/supported-production-surfaces.md)
- [docs/production-readiness-checklist.md](docs/production-readiness-checklist.md)
- [docs/release-validation-runbook.md](docs/release-validation-runbook.md)

Release validation entry point:

```bash
./scripts/run_release_validation.sh
```
