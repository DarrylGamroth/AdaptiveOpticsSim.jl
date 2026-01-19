# Tutorial Mapping Sketches (Julia)

These examples map existing OOPAO tutorials to an idiomatic Julia API in
AdaptiveOptics.jl. They are not complete implementations, but show how the
object graph and propagation steps can translate to multiple dispatch and
explicit state.

## Shack-Hartmann WFS (from `tutorials/AO_closed_loop_ShackHartmann_WFS.py`)

```julia
using AdaptiveOptics

nsub = 20
tel = Telescope(resolution=6nsub, diameter=8.0, sampling_time=1e-3,
                central_obstruction=0.1, fov_arcsec=0.0)

ngs = Source(band=:I, magnitude=8.0, coordinates=(0.0, 0.0))
src = Source(band=:K, magnitude=8.0, coordinates=(0.0, 0.0))

atm = KolmogorovAtmosphere(tel; r0=0.15, L0=25.0,
                           fractional_r0=[0.45,0.1,0.1,0.25,0.1],
                           wind_speed=[10,12,11,15,20],
                           wind_direction=[0,72,144,216,288],
                           altitude=[0,1000,5000,10000,12000])
initialize!(atm, tel)
advance!(atm, tel)

wfs = ShackHartmann(tel; n_subap=20, mode=Geometric())
cam = Detector(integration_time=tel.params.sampling_time, photon_noise=true,
               readout_noise=0.5, qe=1.0, psf_sampling=2, binning=1)
# Detector noise model is encoded in the type: Detector{NoisePhotonReadout, ...}

ws = Workspace(tel.params.resolution)

propagate!(ws, ngs, atm, tel, wfs)  # slopes, valid mask, spot cubes
propagate!(ws, src, tel, cam)       # PSF image on detector
```

Trait-based sensing mode switch (diffractive vs. geometric):
```julia
wfs = ShackHartmann(tel; n_subap=20, mode=Geometric())
# sensing_mode(wfs) returns Geometric() from the type parameter.
```

## Asterism and PSF Grid (from `tutorials/how_to_asterism.py`)

```julia
sources = Source[]
for x in range(-6.0, 6.0; length=5), y in range(-6.0, 6.0; length=5)
    r = hypot(x, y)
    θ = rad2deg(atan(y, x))
    push!(sources, Source(band=:I, magnitude=0.0, coordinates=(r, θ)))
end

ast = Asterism(sources)
dm = DeformableMirror(tel; n_act=20)

propagate!(ws, ast, atm, tel, dm, cam)

# Results:
# tel.state.psf -> combined PSF
# tel.state.psf_list[i] -> per-source PSF
```

## Image Formation + Zernike (from `tutorials/image_formation.py`)

```julia
tel = Telescope(resolution=120, diameter=8.0, sampling_time=1e-3,
                central_obstruction=0.1)
ngs = Source(band=:I, magnitude=10.0)

ws = Workspace(tel)
propagate!(ws, ngs, tel)
psf = compute_psf!(ws, tel; zero_padding=4)

Z = ZernikeBasis(tel, n_modes=100)
modes = compute!(Z, tel)

amp = 100e-9
n = 80
apply_opd!(tel, modes[:, :, n] .* amp)

psf_aberr = compute_psf!(ws, tel; zero_padding=8)
```

## Closed-loop sketch (minimal)

```julia
recon = ModalReconstructor(interaction_matrix)
for k in 1:n_iter
    advance!(atm, tel)
    propagate!(ws, ngs, atm, tel, wfs)
    dm_command = reconstruct!(recon, wfs.slopes)
    dm.state.coefs .= dm_command
    apply!(dm, tel, DMAdditive())
    propagate!(ws, ngs, atm, tel, dm, cam)
end
```

## Logging example

```julia
using Logging

global_logger(ConsoleLogger(stderr, Logging.Info))

@info "Starting AO loop" n_iter=100
for k in 1:n_iter
    advance!(atm, tel)
    propagate!(ws, ngs, atm, tel, wfs)
    dm_command = reconstruct!(recon, wfs.slopes)
    dm.state.coefs .= dm_command
    apply!(dm, tel, DMAdditive())
    @debug "Loop step complete" iter=k wfe_rms=compute_wfe(tel)
end
@info "AO loop complete"
```
