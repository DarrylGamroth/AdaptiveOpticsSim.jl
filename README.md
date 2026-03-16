# AdaptiveOpticsSim.jl

Julia adaptive optics simulation toolkit (in development). This package is an
idiomatic Julia port of OOPAO with a focus on performance, reproducibility, and
extensible modeling.

## Quick start

```julia
using AdaptiveOpticsSim

tel = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3, central_obstruction=0.1)
src = Source(band=:I, magnitude=8.0)
psf = compute_psf!(tel, src; zero_padding=2)
```

Runnable tutorial ports live in `examples/tutorials/`:

```bash
julia --project examples/tutorials/image_formation.jl
julia --project examples/tutorials/closed_loop_pyramid.jl
```

Reference-bundle regeneration can be pinned directly to OOPAO upstream:

```bash
python3 scripts/generate_oopao_reference_bundle.py /tmp/oopao-bundle \
  --oopao-repo https://github.com/cheritier/OOPAO.git \
  --oopao-ref 085d5e50ace0d20fe13cc2da20129d5400166973
```

## Documentation

- `docs/user-guide.md`
- `docs/api-reference.md`
- `docs/julia-port-design.md`
- `docs/julia-tutorial-mappings.md`
- `docs/roadmap.md`
- `docs/oopao-reference-datasets.md`
- `docs/deterministic-simulation.md`
- `docs/phase7-traceability.md`
