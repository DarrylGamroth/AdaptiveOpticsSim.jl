#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
COMPARISON_ROOT="${ROOT_DIR%/AdaptiveOpticsSim.jl}/AdaptiveOpticsComparisons"
STAMP="$(date +%F)"

run_step() {
    local label="$1"
    shift
    echo "==> ${label}"
    "$@"
}

cd "${ROOT_DIR}"

if [[ "${ADAPTIVEOPTICS_SKIP_CPU_FULL_TESTS:-0}" != "1" ]]; then
    run_step "CPU full test suite" \
        julia --project=. --startup-file=no -e 'using Pkg; Pkg.test()'
else
    echo "==> Skipping CPU full test suite (ADAPTIVEOPTICS_SKIP_CPU_FULL_TESTS=1)"
fi

if [[ "${ADAPTIVEOPTICS_VALIDATE_CUDA:-0}" == "1" ]]; then
    run_step "CUDA test project instantiate" \
        julia --project=test/cuda --startup-file=no -e 'using Pkg; Pkg.instantiate()'
    run_step "CUDA hardware target" \
        julia --project=test/cuda --startup-file=no test/runtests_cuda.jl
fi

if [[ "${ADAPTIVEOPTICS_VALIDATE_AMDGPU:-0}" == "1" ]]; then
    run_step "AMDGPU test project instantiate" \
        julia --project=test/amdgpu --startup-file=no -e 'using Pkg; Pkg.instantiate()'
    run_step "AMDGPU hardware target" \
        julia --project=test/amdgpu --startup-file=no test/runtests_amdgpu.jl
fi

if [[ "${ADAPTIVEOPTICS_VALIDATE_EXAMPLES:-0}" == "1" ]]; then
    run_step "Core example scripts" \
        ./scripts/run_core_examples.sh
fi

if [[ "${ADAPTIVEOPTICS_VALIDATE_COMPARISONS:-0}" == "1" ]]; then
    if [[ -d "${COMPARISON_ROOT}" ]]; then
        run_step "HEART all-package runtime ladder" \
            bash -lc "cd '${COMPARISON_ROOT}' && julia --project=. --startup-file=no julia/runners/run_cross_package_benchmarks.jl contracts/heart_hil.toml results/archived/${STAMP}-heart-hil-all-packages.toml"
    else
        echo "==> Skipping cross-package comparisons: ${COMPARISON_ROOT} not found"
    fi
fi


if [[ "${ADAPTIVEOPTICS_VALIDATE_TRUTH:-0}" == "1" ]]; then
    if [[ -d "${ROOT_DIR%/AdaptiveOpticsSim.jl}/REVOLT" ]]; then
        run_step "HEART boundary truth artifact" \
            python3 scripts/generate_heart_boundary_truth_artifact.py
    else
        echo "==> Skipping HEART truth artifact: ${ROOT_DIR%/AdaptiveOpticsSim.jl}/REVOLT not found"
    fi
fi
echo "==> Release validation completed"
