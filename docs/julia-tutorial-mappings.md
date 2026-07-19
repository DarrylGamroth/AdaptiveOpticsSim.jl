# Tutorial Mapping Guide

Status: active

This document maps high-value OOPAO tutorials to deterministic Julia examples in
`examples/tutorials/`. The Julia versions favor multiple dispatch, explicit
state transitions, and small runnable scripts over notebook-style mutation.

## How to run

From the package root:

```bash
julia --project examples/tutorials/image_formation.jl
julia --project examples/tutorials/closed_loop_pyramid.jl
```

Each script exposes a `main()` function and logs a short completion summary with
`Logging.jl`.

## Mapping table

| OOPAO tutorial | Julia example | Coverage |
| --- | --- | --- |
| `tutorials/image_formation.py` | `examples/tutorials/image_formation.jl` | Telescope, source-scaled direct image, Zernike aberrations |
| `tutorials/how_to_detector.py` | `examples/tutorials/detector.jl` | Detector sampling, binning, noise model wiring |
| `tutorials/how_to_multi_sources.py` | `examples/tutorials/asterism.jl` | Multiple sources combined through `Asterism` |
| `tutorials/how_to_asterism.py` | `examples/tutorials/asterism.jl` | Per-source photon-arrival-rate images and their compatible incoherent sum |
| `tutorials/how_to_spatial_filter.py` | `examples/tutorials/spatial_filter.jl` | Spatial filtering without baking optics into the telescope type |
| `tutorials/how_to_NCPA.py` | `examples/tutorials/ncpa.jl` | Basis-driven NCPA synthesis and application |
| `tutorials/how_to_LIFT.ipynb` | `examples/tutorials/lift.jl` | LiFT setup and coefficient recovery |
| `tutorials/how_to_SPRINT.py` | `examples/tutorials/sprint.jl` | Mis-registration sensitivity and estimation |
| `tutorials/AO_transfer_function.py` | `examples/tutorials/transfer_function.jl` | Closed-loop rejection and closed-loop transfer functions |
| `tutorials/AO_closed_loop_ShackHartmannWFS_WFS.py` | `examples/tutorials/closed_loop_shack_hartmann.jl` | Deterministic diffractive SH loop; compact OOPAO regression trace committed |
| `tutorials/AO_closed_loop_Pyramid_WFS.py` | `examples/tutorials/closed_loop_pyramid.jl` | Deterministic diffractive Pyramid loop; compact OOPAO regression trace committed |
| `tutorials/AO_closed_loop_BioEdge_WFS.py` | `examples/tutorials/closed_loop_bioedge.jl` | Deterministic diffractive BioEdge loop; compact OOPAO regression trace committed |
| `tutorials/AO_closed_loop_Pyramid_WFS_GSC.py` | `examples/tutorials/gain_sensing_camera.jl` | Pyramid modulation-frame, optical-gain estimation, and compact GSC closed-loop regression |
| `tutorials/how_to_tomography.py` | `examples/tutorials/tomography.jl` | Compact model-based tomography workflow plus committed pyTomoAO KAPA regression for wavefront and DM-command reconstruction |

## Julia patterns behind the mapping

- OOPAO’s `ngs*tel*wfs` chain becomes explicit preparation and execution such
  as `prepare_direct_imaging` plus `form_direct_image!`, or
  `measure!(wfs, pupil, src)`.
- WFS sensing mode is encoded in the WFS type parameter via
  `mode=Geometric()` or `mode=Diffractive()`, not a mutable string flag.
- Detector noise is encoded by the detector’s `noise` type, for example
  `Detector(noise=(NoisePhoton(), NoiseReadout(0.5)))`.
- Closed-loop examples preallocate their work buffers and use `reconstruct!`
  for the hot path.

## Representative translations

### Image formation

```julia
tel = Telescope(resolution=32, diameter=8.0, central_obstruction=0.1)
src = Source(band=:I, magnitude=10.0)
pupil = PupilFunction(tel)
imaging = prepare_direct_imaging(pupil, src; zero_padding=2)
form_direct_image!(imaging)
photon_rate_image = intensity_values(direct_imaging_output(imaging))
```

### Diffractive WFS measurement

```julia
wfs = PyramidWFS(tel; pupil_samples=4, mode=Diffractive(), modulation=1.0, modulation_points=4)
pupil = PupilFunction(tel)
slopes = measure!(wfs, pupil, src)
```

### Closed loop

```julia
calibration_pupil = PupilFunction(tel)
imat = interaction_matrix(dm, wfs, calibration_pupil, src; amplitude=1e-9)
recon = ModalReconstructor(imat; gain=0.4)
cmd = similar(dm.state.coefs)
renderer = prepare_atmosphere_renderer(atm, tel, src)
pupil = PupilFunction(tel)

epoch = advance_by!(atm, 1e-3; rng=rng)
render_atmosphere!(pupil, renderer, atm, epoch)
measure!(wfs, pupil, src)
reconstruct!(cmd, recon, slopes(wfs))
dm.state.coefs .= -cmd
update_surface!(dm)
apply_surface!(pupil, dm, DMAdditive())
```

## Logging in examples

Use structured logging instead of print statements:

```julia
using Logging

@info "Closed-loop Pyramid tutorial complete" final_residual=result.residual_after[end]
```

That keeps examples composable in scripts, tests, and notebooks.
