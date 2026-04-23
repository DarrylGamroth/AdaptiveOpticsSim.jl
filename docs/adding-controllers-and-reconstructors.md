# Adding Controllers And Reconstructors

Status: active

## Purpose

This guide is the maintained authoring path for adding a new controller or
reconstructor family to `AdaptiveOpticsSim.jl`.

Use it when you need to:

- add a new `AbstractController`
- add a new `AbstractReconstructorOperator`
- decide whether a new control behavior belongs in the controller layer or the
  reconstructor layer

Use together with:

- [api-reference.md](./api-reference.md)
- [runtime-dataflow.md](./runtime-dataflow.md)
- [interface-contract-closeout-plan-2026-04.md](./interface-contract-closeout-plan-2026-04.md)

## Design Rule

The control stack is organized around this split:

- reconstructors map measured WFS signals to command-like vectors
- controllers own temporal state and update laws

In practice:

- if the behavior is a linear or mapped slopes-to-command operator, it belongs
  in the reconstructor layer
- if the behavior adds dynamic state, lag, or temporal filtering, it belongs in
  the controller layer
- do not blur those two surfaces into one monolithic type

## Controllers

### Maintained Contract

Controllers are subtypes of `AbstractController`.

The maintained controller contract is:

- `update!(ctrl, slopes, dt)`
- `controller_output(ctrl)`
- optionally `reset_controller!(ctrl)`
- optionally `supports_controller_reset(ctrl)`

The key rule is:

- `update!` advances the internal controller state
- `controller_output(ctrl)` returns the current command-like output buffer

Current maintained example:

- [controller.jl](../src/Control/controller.jl)

### When To Add A Controller

Add a new controller when you need:

- temporal integration
- lag dynamics
- stateful filtering
- predictive or history-dependent control law behavior

Do not add a controller when a pure static linear slopes-to-command map is
enough. That belongs in a reconstructor.

### Validation Checklist

Every new controller should add evidence in
[runtests_head.jl](../test/runtests_head.jl) and runtime-facing tests when the
controller is part of a maintained runtime surface.

Minimum expectations:

1. `update!` returns or preserves a stable output buffer
2. `controller_output(ctrl)` has the expected shape
3. `reset_controller!` works when supported
4. no hidden allocations are introduced into the maintained hot path

## Reconstructors

### Maintained Contract

Reconstructors are subtypes of `AbstractReconstructorOperator`.

The maintained reconstructor contract is:

- `reconstruct!(out, recon, slopes)`
- `reconstruct(recon, slopes)`

And the maintained diagnostic accessors are:

- `inverse_policy(recon)`
- `singular_values(recon)`
- `condition_number(recon)`
- `effective_rank(recon)`

Current maintained families:

- `NullReconstructor`
- `ModalReconstructor`
- `MappedReconstructor`

### When To Add A Reconstructor

Add a new reconstructor when you need:

- a different slopes-to-command inverse operator
- a different command-basis mapping
- a different diagnostic or inverse-policy surface

Do not add controller-like temporal state here. Keep reconstructors static and
structural.

### Validation Checklist

Every new reconstructor should add evidence in
[calibration_and_analysis.jl](../test/testsets/calibration_and_analysis.jl)
and runtime tests when the operator participates in closed-loop runtime
control.

Minimum expectations:

1. `reconstruct!` validates input/output dimensions correctly
2. `reconstruct!` and `reconstruct` agree
3. diagnostic accessors are defined or rejected explicitly
4. runtime integration works through the maintained closed-loop path

## Reuse Rule

This is the main architectural test:

- if the behavior is static and operator-like, keep it in the reconstructor
  layer
- if the behavior is dynamic and stateful, keep it in the controller layer
- if a second family needs the same temporal or operator substep, factor it
  into a shared helper instead of duplicating code

## What Not To Do

- do not put temporal controller state into reconstructors
- do not hide controller output behind undocumented fields
- do not make runtime code depend on family-local controller or reconstructor
  storage details
- do not add new control behavior through `isa` branches in the runtime path
