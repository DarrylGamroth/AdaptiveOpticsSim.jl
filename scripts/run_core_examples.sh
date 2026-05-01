#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

examples=(
    "examples/closed_loop_demo.jl"
    "examples/closed_loop/control_loop_single_runtime.jl"
    "examples/closed_loop/control_loop_grouped_runtime.jl"
    "examples/closed_loop/run_cl.jl"
    "examples/closed_loop/run_cl_first_stage.jl"
    "examples/closed_loop/run_cl_from_phase_screens.jl"
    "examples/closed_loop/run_cl_long_push_pull.jl"
    "examples/closed_loop/run_cl_sinusoidal_modulation.jl"
    "examples/closed_loop/run_cl_two_stages.jl"
    "examples/closed_loop/run_cl_two_stages_atm_change.jl"
    "examples/closed_loop/run_cl_zernike.jl"
    "examples/tutorials/asterism.jl"
    "examples/tutorials/closed_loop_bioedge.jl"
    "examples/tutorials/closed_loop_pyramid.jl"
    "examples/tutorials/closed_loop_shack_hartmann.jl"
    "examples/tutorials/closed_loop_zernike.jl"
    "examples/tutorials/detector.jl"
    "examples/tutorials/extended_source_sensing.jl"
    "examples/tutorials/gain_sensing_camera.jl"
    "examples/tutorials/image_formation.jl"
    "examples/tutorials/lift.jl"
    "examples/tutorials/ncpa.jl"
    "examples/tutorials/shack_hartmann_subapertures.jl"
    "examples/tutorials/spatial_filter.jl"
    "examples/tutorials/sprint.jl"
    "examples/tutorials/tomography.jl"
    "examples/tutorials/transfer_function.jl"
)

cd "${ROOT_DIR}"

for example in "${examples[@]}"; do
    echo "==> ${example}"
    julia --project=. --startup-file=no "${example}"
done

echo "==> Core examples completed"
