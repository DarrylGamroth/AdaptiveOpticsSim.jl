# Phase Statistics Accuracy Note

Status: active

Plan traceability:

- [post-review-platform-plan.md](./post-review-platform-plan.md) `PVP-01`
- [model-validity-matrix.md](./model-validity-matrix.md) `MV-02`

## Purpose

This note summarizes the maintained accuracy evidence for the shared
`K_{5/6}` helper used by:

- [`phase_covariance`](../src/Atmosphere/phase_stats.jl)
- tomography covariance paths in
  [`reconstructors.jl`](../src/Tomography/reconstructors.jl)

The implementation lives in [kv56.jl](../src/Core/kv56.jl) and evaluates the
fixed-order quantity needed by the package rather than a generic `besselk(ν,x)`
surface.

## Implemented Approximation

The helper is split into three regimes:

- small `x`: truncated series
- intermediate `x`: fixed-range numerical quadrature over `t ∈ [0, 9]`
- large `x`: asymptotic expansion

The maintained branch points are:

- small cutoff: `x < 0.1`
- large cutoff: `x >= 20`
- intermediate bins: `64`

## Reference Baseline

The reference function is `SpecialFunctions.besselk(5/6, x)` on CPU.

The maintained helper is validated in two forms:

- plain `K_{5/6}(x)` via `AdaptiveOpticsSim._kv56_scalar(x)`
- scaled `x^(5/6) K_{5/6}(x)` via `AdaptiveOpticsSim._scaled_kv56_scalar(x)`

The scaled form is the numerically natural quantity for the current covariance
surfaces, but both are checked because the unscaled helper is still used in
tests and diagnostics.

## Validation Range

Maintained sweep range:

- `x ∈ [1e-6, 140]`

Sampling used for the maintained summary below:

- `4000` log-spaced points over `[1e-6, 140]`
- explicit inclusion of `{0.1, 20, 40, 80, 140}`

This range covers the intended operational envelope for the current atmosphere
and tomography covariance use.

## Summary Results

Measured against `SpecialFunctions.besselk(5/6, x)`:

### Plain `K_{5/6}(x)`

- max relative error: `5.3624489134393e-8`
- mean relative error: `6.488853084596498e-10`
- 95th percentile relative error: `1.6474627526975024e-9`
- worst sampled point: `x = 20.0`

### Scaled `x^(5/6) K_{5/6}(x)`

- max relative error: `5.362448921425661e-8`
- mean relative error: `6.48885307067207e-10`
- 95th percentile relative error: `1.6474626809472649e-9`
- worst sampled point: `x = 20.0`

## Interpretation

- The helper is accurate enough for the maintained covariance surfaces.
- The worst point appears at the regime boundary around `x = 20`, but the
  error remains below `6e-8`.
- There is no practical difference between the plain and scaled summaries at
  this sweep density.

## Maintained Evidence Surfaces

- direct regression checks in
  [calibration_and_analysis.jl](../test/testsets/calibration_and_analysis.jl)
- GPU contract coverage through
  [gpu_smoke_contract.jl](../scripts/gpu_smoke_contract.jl)
- this note as the maintained summary artifact

## Reproduction

The summary above was generated from a direct Julia sweep against
`SpecialFunctions` using the project environment:

```bash
julia --project=. --startup-file=no -e '
using AdaptiveOpticsSim, SpecialFunctions, Statistics
xs = sort!(vcat(collect(10.0 .^ range(-6, log10(140.0), length=4000)), [0.1, 20.0, 40.0, 80.0, 140.0]))
rel_k = Float64[]
rel_scaled = Float64[]
worst_k = (x=0.0, err=-1.0)
worst_scaled = (x=0.0, err=-1.0)
for x in xs
    ref = SpecialFunctions.besselk(5/6, x)
    approx = AdaptiveOpticsSim._kv56_scalar(x)
    scaled_ref = x^(5/6) * ref
    scaled_approx = AdaptiveOpticsSim._scaled_kv56_scalar(x)
    err_k = abs((approx - ref) / ref)
    err_scaled = abs((scaled_approx - scaled_ref) / scaled_ref)
    push!(rel_k, err_k)
    push!(rel_scaled, err_scaled)
    if err_k > worst_k.err
        global worst_k = (x=x, err=err_k)
    end
    if err_scaled > worst_scaled.err
        global worst_scaled = (x=x, err=err_scaled)
    end
end
println("plain_max_rel=", maximum(rel_k))
println("plain_mean_rel=", mean(rel_k))
println("plain_p95_rel=", quantile(rel_k, 0.95))
println("plain_worst_x=", worst_k.x)
println("scaled_max_rel=", maximum(rel_scaled))
println("scaled_mean_rel=", mean(rel_scaled))
println("scaled_p95_rel=", quantile(rel_scaled, 0.95))
println("scaled_worst_x=", worst_scaled.x)
'
```
