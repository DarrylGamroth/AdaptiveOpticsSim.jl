# Independent assessment of shifted-mask equivalence

Date: 2026-10-22. Reviewed source: `AdaptiveOpticsSim-batches`, based on `e7b4fbb4`, before removal. No production changes or GPU jobs were performed by this reviewer.

## Decision

Remove the public shifted-mask propagation strategy and migrate callers to physical pupil tilt. The observed errors are consistent with a structural difference between the finite discrete operators. A bounded source and CPU analysis identified no isolated sign, scale, or axis correction capable of restoring general equivalence. A periodic-wrap correction would address only one smaller error and would leave fractional modulation errors.

This recommendation implements the user's stated preference to remove the approximation unless it can be improved defensibly. It does not claim that an alternative accurate fast algorithm is impossible; such an algorithm would need its own design and validation.

## Source observations

- `_prepare_pyramid_shifted_focal_field!` prepares one unmodulated, centered, padded pupil FFT.
- `_build_pyramid_shifted_masks!` evaluates the analytic pyramid mask at translated coordinates. The translated angles are `2π × offset / (resolution − 1)`, with axis swapping matching the pupil modulation convention. Both dense and separable implementations encode this same operation.
- Subsequent propagation multiplies that FFT pointwise by each translated mask and takes an inverse FFT.
- Physical pupil tilt multiplies the pupil by the modulation phase before its forward FFT, then applies the fixed sampled mask.
- The effective shift in FFT-bin units is `δ = offset × padded_side / (resolution − 1)`. Circular modulation generally produces fractional bin shifts.

## Derived discrete-operator difference

Let `F` be the finite unitary DFT, `Dδ` the diagonal pupil tilt, and `M` the diagonal sampled pyramid mask. Physical propagation is:

```text
y = F⁻¹ M F Dδ u
```

Removing an output tilt preserves intensity. With `Sδ = F Dδ F⁻¹`, the equivalent unmodulated-input focal operator is:

```text
Dδ⁻¹ y = F⁻¹ [Sδ⁻¹ M Sδ] F u
```

For integer bin shifts, `Sδ` is a cyclic permutation. Its conjugation of `M` is a periodically shifted diagonal mask, and the mask-shift shortcut is exact when periodic wrapping is respected.

For fractional bin shifts, `Sδ` is a dense finite Fourier interpolation operator. Its conjugation of a nonconstant pyramid mask is generally non-diagonal. Replacing that operator with pointwise samples of the translated analytic mask changes the finite discrete model. The continuous Fourier shift theorem does not make this replacement exact at fixed discrete sampling.

## Observed CPU discrimination

`shifted-mask-discrete-check.py` reconstructs the source equations using NumPy on the CPU. It uses pupil resolution 120, 30 pupil samples, separation 2, padded side 256, a circular pupil, and deterministic nonzero OPD. It compares raw propagated intensity, before detector binning. Raw results are in `shifted-mask-discrete-check.json`.

| Shift in FFT bins | Analytic translated mask versus physical pupil tilt, relative L2 | Periodically shifted sampled mask versus physical pupil tilt |
|---:|---:|---:|
| 0 | 0 | 0 |
| 0.25 | 0.1927 | Not an integer permutation |
| 0.5 | 0.4025 | Not an integer permutation |
| 1 | 0.00424 | 5.47 × 10⁻¹⁶ |
| 1.5 | 0.1898 | Not an integer permutation |
| 2 | 0.00717 | 5.38 × 10⁻¹⁶ |
| 3 | 0.00831 | 5.82 × 10⁻¹⁶ |
| 6 | 0.01206 | 7.26 × 10⁻¹⁶ |

The integer-permutation result validates the shift sign and scale in the reconstructed identity. The smaller integer-shift discrepancy in the current analytic mask also identifies an edge-wrapping difference. Correcting that edge behavior alone cannot address the much larger fractional-bin mismatch.

A separate 32-point explicit matrix experiment evaluates the exact conjugated mask. The relative Frobenius norm of its off-diagonal entries is approximately 0.394 at a quarter-bin shift and 0.557 at a half-bin shift. It is below 6.2 × 10⁻¹⁶ at the tested integer shifts. This directly demonstrates why a diagonal translated mask cannot represent that conjugated finite operator.

These results complement the actual package/OOPAO comparison: the existing shifted strategy has roughly 2.30% flat-image, 4.29% ramp-image, and 27.85% tip-derivative relative errors on the recorded fixture, while the physical pupil-tilt implementation agrees much more closely. The mathematical experiment establishes a structural mismatch; it does not attribute every measured package error exclusively to one mechanism.

## Stop criterion and required follow-through

The bounded investigation has reached its stop criterion: fractional-bin mismatch remains intrinsic to the current pointwise-mask formulation, and the commensurate control validates the basic conventions. Do not try coefficient tuning, coordinate changes, or interpolation substitutions as speculative production repairs.

Remove the public strategy, its implementation-only owners and helpers, and active callers across AOS, GPI, and COPPER. Preserve historical benchmark and external-comparison evidence with explicit superseded-model labeling. Keep physical pupil tilt and validate representative photons, calibration, captured execution, and closed-loop behavior after migration. Reassess timing claims using that retained model; the earlier shifted-mask speedups are historical implementation measurements, not qualification of an accepted physical model.
