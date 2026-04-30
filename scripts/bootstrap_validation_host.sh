#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

need_cmd() {
    local cmd="$1"
    if ! command -v "${cmd}" >/dev/null 2>&1; then
        echo "missing required command: ${cmd}" >&2
        exit 1
    fi
}

bootstrap_project() {
    julia --project="${ROOT_DIR}" --startup-file=no -e '
        using Pkg
        Pkg.instantiate()
        Pkg.precompile()
    '
}

bootstrap_cuda() {
    need_cmd nvidia-smi
    julia --project="${ROOT_DIR}/test/cuda" --startup-file=no -e '
        using Pkg
        Pkg.instantiate()
        Pkg.precompile()
    '
    nvidia-smi >/dev/null
}

bootstrap_amdgpu() {
    need_cmd rocminfo
    julia --project="${ROOT_DIR}/test/amdgpu" --startup-file=no -e '
        using Pkg
        Pkg.instantiate()
        Pkg.precompile()
    '
    rocminfo >/dev/null
}

main() {
    need_cmd julia
    need_cmd git

    bootstrap_project

    for backend in "$@"; do
        case "${backend}" in
            cuda)
                bootstrap_cuda
                ;;
            amdgpu)
                bootstrap_amdgpu
                ;;
            *)
                echo "unknown backend: ${backend}" >&2
                echo "usage: $0 [cuda] [amdgpu]" >&2
                exit 1
                ;;
        esac
    done

    echo "validation host bootstrap completed"
}

main "$@"
