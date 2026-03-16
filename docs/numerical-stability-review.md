# Numerical Stability Review

Date: 2026-03-15
Commit reviewed: `ca8a97d`

## Scope

This review focuses on numerical stability risks in the current `AdaptiveOptics.jl`
implementation, not on line-by-line equivalence with OOPAO. The target is feature
parity with equal or better numerical behavior where practical.

The review was based on static inspection of the current implementation plus the
existing deterministic test suite. It is not yet backed by full condition-number
sweeps or floating-point error characterization across parameter ranges.

## Summary

The main stability concerns are concentrated in inverse problems and covariance-based
reconstructors:

1. Modal/control reconstructors use pseudoinverses or reciprocal singular values
   with weak conditioning policy.
2. Tomography reconstructors use heuristic diagonal loading and dense right-division
   on covariance systems that should be treated as structured linear solves.
3. LiFT currently solves weighted least squares through normal equations, which
   squares the condition number.
4. Gain Sensing Camera optical gains divide by calibration sensitivities without
   an internal sensitivity floor.
5. The custom covariance/Bessel approximations likely need regime validation and
   guardrails at extreme parameter values.

## Findings

### High: Unconstrained pseudoinverse use in modal reconstruction

Files:
- [src/Calibration/reconstructor.jl](/home/dgamroth/workspaces/codex/AdaptiveOptics.jl/src/Calibration/reconstructor.jl#L8)
- [src/Calibration/calibration_vault.jl](/home/dgamroth/workspaces/codex/AdaptiveOptics.jl/src/Calibration/calibration_vault.jl#L31)
- [src/Calibration/modal_basis.jl](/home/dgamroth/workspaces/codex/AdaptiveOptics.jl/src/Calibration/modal_basis.jl#L29)

Current behavior:
- `ModalReconstructor` builds `pinv(imat.matrix)`.
- `CalibrationVault` inverts all nonzero singular values.
- `basis_projector` falls back to `pinv(basis)` when the basis is not sufficiently
  diagonal under a heuristic criterion.

Risk:
- Small singular values are amplified aggressively.
- No explicit rank cutoff, relative tolerance, or regularization policy is exposed.
- No diagnostic indicates when a matrix is effectively rank-deficient.

Impact:
- Noise amplification in closed-loop control.
- Actuator chatter or unstable controller behavior for nearly dependent modes.
- Poor reproducibility across BLAS/LAPACK implementations when systems are ill-conditioned.

Recommendation:
- Introduce explicit inverse policies: TSVD, Tikhonov, and exact pseudoinverse only
  when requested.
- Report effective rank, condition number, and truncation threshold.
- Make the conditioning policy part of calibration/reconstruction configuration.

### High: Tomography solve path relies on heuristic diagonal loading

Files:
- [src/Tomography/reconstructors.jl](/home/dgamroth/workspaces/codex/AdaptiveOptics.jl/src/Tomography/reconstructors.jl#L441)
- [src/Tomography/reconstructors.jl](/home/dgamroth/workspaces/codex/AdaptiveOptics.jl/src/Tomography/reconstructors.jl#L443)
- [src/Tomography/reconstructors.jl](/home/dgamroth/workspaces/codex/AdaptiveOptics.jl/src/Tomography/reconstructors.jl#L671)
- [src/Tomography/reconstructors.jl](/home/dgamroth/workspaces/codex/AdaptiveOptics.jl/src/Tomography/reconstructors.jl#L673)

Current behavior:
- `Cnz` is synthesized with simple heuristics like `1e-3 * α` or `mean_diag / 10`.
- `css` is formed densely and solved through right-division.

Risk:
- The regularization scale is not tied to a physical or statistical noise model.
- Different WFS/channel scalings can change the effective conditioning substantially.
- Dense right-division hides whether the system is close to singular.

Impact:
- Tomographic reconstructor quality can vary non-physically with scaling choices.
- Poor robustness when moving from compact parity cases to larger systems.

Recommendation:
- Treat `css` as a Hermitian positive-definite solve when possible.
- Use Cholesky or LDLt with explicit fallback behavior.
- Derive `Cnz` from detector/noise assumptions or expose a structured user policy.

### Medium: LiFT still needs explicit damping / conditioning policy

File:
- [src/WFS/lift.jl](/home/dgamroth/workspaces/codex/AdaptiveOptics.jl/src/WFS/lift.jl#L237)

Current behavior:
- LiFT now defaults to `LiFTSolveAuto()`.
- `LiFTSolveAuto()` resolves to QR on scalar CPU arrays and to the
  normal-equation path on accelerator backends.
- An explicit `LiFTSolveNormalEquations()` compatibility mode still exists.
- Rank-deficient QR cases currently fall back to a damped normal-equation path.

Risk:
- The default path is materially better than pure normal equations, but there is
  still no explicit damping or conditioning policy.
- The rank-deficient fallback still uses the weaker normal-equation solve.

Impact:
- Mode estimates may become unstable for weakly observed modes or noisy PSFs.
- Iterative convergence can degrade or stall in high-dynamic-range cases.

Recommendation:
- Keep backend-aware `Auto` solve selection, with QR preferred on scalar CPU.
- Add optional Levenberg-Marquardt damping for the nonlinear iteration.
- Track residual norm, update norm, and Jacobian conditioning per iteration.

### Medium: Gain Sensing Camera optical gains lack a calibration-sensitivity floor

File:
- [src/Calibration/gain_sensing_camera.jl](/home/dgamroth/workspaces/codex/AdaptiveOptics.jl/src/Calibration/gain_sensing_camera.jl#L120)

Current behavior:
- Optical gains are computed as `real(sensi_sky / sensi_calib)`.

Risk:
- Modes with very small `sensi_calib` produce unstable or misleading gains.
- Some workflows floor gains later, but the primitive itself does not.

Impact:
- Large or sign-flipping gains on weak modes.
- Harder to distinguish physical gain loss from numerical instability.

Recommendation:
- Add a calibration sensitivity floor or explicit invalid-mode mask.
- Report which modes are weakly observed at calibration time.

### Medium: Custom covariance / Bessel approximations need validation at regime boundaries

Files:
- [src/Atmosphere/phase_stats.jl](/home/dgamroth/workspaces/codex/AdaptiveOptics.jl/src/Atmosphere/phase_stats.jl#L47)
- [src/Tomography/reconstructors.jl](/home/dgamroth/workspaces/codex/AdaptiveOptics.jl/src/Tomography/reconstructors.jl#L74)

Current behavior:
- `phase_covariance` uses a direct `besselk` expression with a special case only at exact zero.
- `_kv56_scalar` uses a custom small-argument series and large-argument asymptotic expansion.

Risk:
- Accuracy near the switching boundary and at very small `rho` / `z` is not yet characterized.
- Extreme parameter combinations may cause avoidable cancellation or underflow.

Impact:
- Small but systematic errors in covariance assembly.
- Hard-to-diagnose tomography drift when scaling up system size or changing atmospheric parameters.

Recommendation:
- Sweep these kernels against trusted references across operating ranges.
- Add regime-specific guards and error-tested approximations if needed.

## Lower-Risk Notes

- Many hot paths already protect against divide-by-zero with `total <= 0` or `max(..., eps(T))`.
- WFS signal extraction generally handles empty/zero-flux cases defensively.
- The current architecture is compatible with introducing stronger factorization-based
  solvers without breaking the public API.

## Recommended Next Steps

1. Add a configurable inverse-policy layer for reconstructors and calibration objects.
2. Add damping and conditioning diagnostics to the LiFT QR solve.
3. Refactor tomography solve paths to use structured factorizations and explicit noise models.
4. Add GSC sensitivity flooring and weak-mode diagnostics.
5. Add numerical sweep tests for:
   - singular-value truncation sensitivity,
   - reconstructor perturbation sensitivity,
   - covariance-kernel approximation error,
   - LiFT convergence under noise and low observability.

## Status

This is a findings document. The LiFT item has been partially addressed by moving
the scalar CPU default away from normal equations, adding backend-aware solve
selection, and exposing damping/conditioning diagnostics. More advanced damping
policy and stronger conditioning diagnostics are still open.
