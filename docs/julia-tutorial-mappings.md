# Tutorial Mapping Guide

Status: active

This provenance map records high-value OOPAO source tutorials and external
regression traces alongside AOS-native Julia examples in
`examples/tutorials/`. The examples favor multiple dispatch, explicit state
transitions, and small runnable scripts over notebook-style mutation; the
mapping does not impose an OOPAO API or class layout.

## How to run

From the package root:

```bash
julia --project examples/tutorials/image_formation.jl
```

Each script exposes a `main()` function and logs a short completion summary with
`Logging.jl`.

## Mapping table

| OOPAO source tutorial (provenance) | AOS Julia example | Coverage |
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
| `tutorials/AO_closed_loop_ShackHartmannWFS_WFS.py` | [`examples/integrations/filter_graph_algorithms/`](../examples/integrations/filter_graph_algorithms/) | Maintained Shack–Hartmann plant/RTC fixture; AOS owns the plant and FGA/JFG owns the RTC |
| `tutorials/AO_closed_loop_Pyramid_WFS.py` | Downstream package | AOS Pyramid sensing and the FGA Pyramid estimator remain available; their complete RTC composition and acceptance are not an AOS tutorial surface |
| `tutorials/AO_closed_loop_BioEdge_WFS.py` | Downstream package | AOS Bi-O-edge sensing remains available; its RTC composition and acceptance are not an AOS tutorial surface |
| `tutorials/AO_closed_loop_Pyramid_WFS_GSC.py` | `examples/tutorials/gain_sensing_camera.jl` | Pyramid modulation-frame and optical-gain estimation; no maintained AOS RTC composition |
| `tutorials/how_to_tomography.py` | `examples/tutorials/tomography.jl` | Compact model-based tomography workflow plus committed pyTomoAO KAPA regression for wavefront and DM-command reconstruction |

## Julia patterns behind the mapping

- Where an OOPAO tutorial is cited, its `ngs*tel*wfs` source expression maps to
  explicit AOS preparation and execution. Shack–Hartmann composition uses
  `prepare_wfs_optics` plus `form_wfs_optical_products!`, followed by explicit
  detector acquisition; operational estimation belongs to FGA.
- Sensor families that retain alternate sensing models encode them in a type
  parameter, not a mutable string flag. The AOS Shack–Hartmann surface is the
  physical diffractive front end; direct OPD-gradient truth uses the explicitly
  named `geometric_wavefront_slopes!` calculation.
- Detector noise is encoded by the detector’s `noise` type, for example
  `Detector(noise=(NoisePhoton(), NoiseReadout(0.5)))`.
- The maintained closed-loop fixtures keep the AOS plant and the FGA/JFG RTC
  in separate packages with explicit cross-package products.

## Representative AOS mappings

The snippets below assume:

```julia
using AdaptiveOpticsSim
using AdaptiveOpticsSim.Optics
using AdaptiveOpticsSim.WavefrontSensors
using AdaptiveOpticsSim.Calibration
using AdaptiveOpticsSim.Tomography
```

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

## Logging in examples

Use structured logging instead of print statements:

```julia
using Logging

@info "Pyramid sensing tutorial complete" measurement_length=length(slopes(wfs))
```

That keeps examples composable in scripts, tests, and notebooks.
