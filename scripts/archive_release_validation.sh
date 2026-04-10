#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
    echo "usage: $0 <cpu|cuda|amdgpu|all> [host_label]" >&2
    exit 2
fi

TRACK="$1"
HOST_LABEL="${2:-$(hostname -s)}"
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
STAMP="$(date +%F)"
OUTDIR="${ROOT_DIR}/benchmarks/results/validation_runs"
mkdir -p "${OUTDIR}"
LOGFILE="${OUTDIR}/${STAMP}-${HOST_LABEL}-${TRACK}.log"
METAFILE="${OUTDIR}/${STAMP}-${HOST_LABEL}-${TRACK}.toml"

DEFAULT_SKIP_CPU_FULL_TESTS=0

case "${TRACK}" in
    cpu)
        export ADAPTIVEOPTICS_VALIDATE_CUDA=0
        export ADAPTIVEOPTICS_VALIDATE_AMDGPU=0
        ;;
    cuda)
        export ADAPTIVEOPTICS_VALIDATE_CUDA=1
        export ADAPTIVEOPTICS_VALIDATE_AMDGPU=0
        DEFAULT_SKIP_CPU_FULL_TESTS=1
        ;;
    amdgpu)
        export ADAPTIVEOPTICS_VALIDATE_CUDA=0
        export ADAPTIVEOPTICS_VALIDATE_AMDGPU=1
        DEFAULT_SKIP_CPU_FULL_TESTS=1
        ;;
    all)
        export ADAPTIVEOPTICS_VALIDATE_CUDA=1
        export ADAPTIVEOPTICS_VALIDATE_AMDGPU=1
        ;;
    *)
        echo "invalid track: ${TRACK}" >&2
        exit 2
        ;;
esac

export ADAPTIVEOPTICS_SKIP_CPU_FULL_TESTS="${ADAPTIVEOPTICS_SKIP_CPU_FULL_TESTS:-${DEFAULT_SKIP_CPU_FULL_TESTS}}"

cd "${ROOT_DIR}"
GIT_COMMIT="$(git rev-parse HEAD)"
START_TS="$(date --iso-8601=seconds)"
CMD="./scripts/run_release_validation.sh"
STATUS="passed"

{
    echo "==> host: ${HOST_LABEL}"
    echo "==> track: ${TRACK}"
    echo "==> commit: ${GIT_COMMIT}"
    echo "==> started: ${START_TS}"
    echo "==> command: ${CMD}"
} | tee "${LOGFILE}"

if ! "${ROOT_DIR}/scripts/run_release_validation.sh" 2>&1 | tee -a "${LOGFILE}"; then
    STATUS="failed"
fi

END_TS="$(date --iso-8601=seconds)"
cat > "${METAFILE}" <<EOF
artifact_kind = "release_validation_run"
host_label = "${HOST_LABEL}"
track = "${TRACK}"
status = "${STATUS}"
started_at = "${START_TS}"
finished_at = "${END_TS}"
git_commit = "${GIT_COMMIT}"
log_path = "$(basename "${LOGFILE}")"
command = "${CMD}"
validate_cuda = ${ADAPTIVEOPTICS_VALIDATE_CUDA}
validate_amdgpu = ${ADAPTIVEOPTICS_VALIDATE_AMDGPU}
skip_cpu_full_tests = ${ADAPTIVEOPTICS_SKIP_CPU_FULL_TESTS}
EOF

if [[ "${STATUS}" != "passed" ]]; then
    exit 1
fi
