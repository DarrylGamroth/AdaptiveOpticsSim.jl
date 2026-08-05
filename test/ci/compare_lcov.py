#!/usr/bin/env python3
"""Compare the four authoritative shard LCOV reports with a transition baseline."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
import re
import subprocess
import sys
from typing import Iterable, Mapping, Sequence


class CoverageError(RuntimeError):
    """Raised when coverage evidence is incomplete, malformed, or different."""


LineKey = tuple[str, int]
Coverage = dict[LineKey, bool]


@dataclass(frozen=True)
class CoverageTotals:
    attributed: int
    hit: int

    @property
    def missed(self) -> int:
        return self.attributed - self.hit


def _source_path(source: str, repo_root: Path, report: Path) -> str:
    raw_path = Path(source)
    absolute = raw_path if raw_path.is_absolute() else repo_root / raw_path
    try:
        relative = absolute.resolve().relative_to(repo_root)
    except ValueError as error:
        raise CoverageError(
            f"{report}: source path is outside the repository: {source}"
        ) from error
    if not relative.parts or relative.parts[0] not in {"src", "ext"}:
        raise CoverageError(
            f"{report}: source path is outside src/ and ext/: {relative.as_posix()}"
        )
    return relative.as_posix()


def parse_lcov(report: Path, repo_root: Path) -> Coverage:
    """Parse attributed line hit/miss state from one nonempty LCOV report."""
    report = report.resolve()
    if not report.is_file():
        raise CoverageError(f"missing LCOV report: {report}")
    if report.stat().st_size == 0:
        raise CoverageError(f"empty LCOV report: {report}")

    coverage: Coverage = {}
    source: str | None = None
    record_lines: dict[int, bool] = {}
    declared_lines: int | None = None
    declared_hits: int | None = None
    record_count = 0

    def finish_record(line_number: int) -> None:
        nonlocal source, record_lines, declared_lines, declared_hits, record_count
        if source is None:
            raise CoverageError(
                f"{report}:{line_number}: end_of_record without an SF record"
            )
        if not record_lines:
            if declared_lines not in {None, 0} or declared_hits not in {None, 0}:
                raise CoverageError(
                    f"{report}:{line_number}: {source} declares coverage for "
                    "a record with no attributed DA lines"
                )
            source = None
            record_lines = {}
            declared_lines = None
            declared_hits = None
            record_count += 1
            return
        hits = sum(record_lines.values())
        if declared_lines is not None and declared_lines != len(record_lines):
            raise CoverageError(
                f"{report}:{line_number}: LF={declared_lines} but parsed "
                f"{len(record_lines)} DA lines for {source}"
            )
        if declared_hits is not None and declared_hits != hits:
            raise CoverageError(
                f"{report}:{line_number}: LH={declared_hits} but parsed "
                f"{hits} hit lines for {source}"
            )
        for number, hit in record_lines.items():
            key = (source, number)
            coverage[key] = coverage.get(key, False) or hit
        source = None
        record_lines = {}
        declared_lines = None
        declared_hits = None
        record_count += 1

    try:
        lines = report.read_text(encoding="utf-8").splitlines()
    except UnicodeDecodeError as error:
        raise CoverageError(f"{report}: LCOV report is not UTF-8") from error

    for line_number, line in enumerate(lines, 1):
        if line.startswith("SF:"):
            if source is not None:
                raise CoverageError(
                    f"{report}:{line_number}: SF record started before end_of_record"
                )
            source_text = line[3:]
            if not source_text:
                raise CoverageError(f"{report}:{line_number}: empty SF path")
            source = _source_path(source_text, repo_root, report)
        elif line.startswith("DA:"):
            if source is None:
                raise CoverageError(f"{report}:{line_number}: DA outside an SF record")
            fields = line[3:].split(",")
            if len(fields) not in {2, 3}:
                raise CoverageError(f"{report}:{line_number}: malformed DA record")
            try:
                number = int(fields[0])
                count = int(fields[1])
            except ValueError as error:
                raise CoverageError(
                    f"{report}:{line_number}: non-integer DA line or count"
                ) from error
            if number <= 0 or count < 0:
                raise CoverageError(
                    f"{report}:{line_number}: DA line must be positive and count nonnegative"
                )
            if number in record_lines:
                raise CoverageError(
                    f"{report}:{line_number}: duplicate DA line {number} for {source}"
                )
            record_lines[number] = count > 0
        elif line.startswith("LF:"):
            if source is None or declared_lines is not None:
                raise CoverageError(f"{report}:{line_number}: misplaced or duplicate LF")
            try:
                declared_lines = int(line[3:])
            except ValueError as error:
                raise CoverageError(f"{report}:{line_number}: malformed LF") from error
        elif line.startswith("LH:"):
            if source is None or declared_hits is not None:
                raise CoverageError(f"{report}:{line_number}: misplaced or duplicate LH")
            try:
                declared_hits = int(line[3:])
            except ValueError as error:
                raise CoverageError(f"{report}:{line_number}: malformed LH") from error
        elif line == "end_of_record":
            finish_record(line_number)
        elif not line or line.startswith(
            ("TN:", "FN:", "FNDA:", "FNF:", "FNH:", "BRDA:", "BRF:", "BRH:")
        ):
            continue
        else:
            raise CoverageError(f"{report}:{line_number}: unrecognized LCOV record: {line}")

    if source is not None:
        raise CoverageError(f"{report}: unterminated SF record for {source}")
    if record_count == 0 or not coverage:
        raise CoverageError(f"{report}: no attributed src/ or ext/ lines")
    return coverage


def union_coverage(reports: Iterable[Coverage]) -> Coverage:
    combined: Coverage = {}
    for report in reports:
        for key, hit in report.items():
            combined[key] = combined.get(key, False) or hit
    return combined


def changed_lines(repo_root: Path, base_sha: str) -> set[LineKey]:
    """Return src line keys added or modified from base_sha through HEAD."""
    command = [
        "git",
        "-c",
        "core.quotePath=false",
        "diff",
        "--unified=0",
        "--no-color",
        "--no-ext-diff",
        "--diff-filter=ACMR",
        f"{base_sha}...HEAD",
        "--",
        "src",
    ]
    try:
        result = subprocess.run(
            command,
            cwd=repo_root,
            check=True,
            capture_output=True,
            text=True,
        )
    except subprocess.CalledProcessError as error:
        detail = error.stderr.strip() or error.stdout.strip()
        raise CoverageError(f"cannot diff supplied base SHA {base_sha}: {detail}") from error

    changed: set[LineKey] = set()
    current_path: str | None = None
    hunk_pattern = re.compile(r"^@@ -\d+(?:,\d+)? \+(\d+)(?:,(\d+))? @@")
    for line in result.stdout.splitlines():
        if line.startswith("+++ "):
            path = line[4:]
            current_path = None if path == "/dev/null" else path.removeprefix("b/")
            if current_path is not None and not current_path.startswith("src/"):
                raise CoverageError(f"git diff returned an unexpected path: {current_path}")
        elif line.startswith("@@ "):
            if current_path is None:
                continue
            match = hunk_pattern.match(line)
            if match is None:
                raise CoverageError(f"cannot parse git diff hunk: {line}")
            start = int(match.group(1))
            count = int(match.group(2) or "1")
            changed.update((current_path, number) for number in range(start, start + count))
    return changed


def coverage_totals(coverage: Mapping[LineKey, bool], keys: set[LineKey] | None = None) -> CoverageTotals:
    selected = coverage.keys() if keys is None else coverage.keys() & keys
    states = [coverage[key] for key in selected]
    return CoverageTotals(len(states), sum(states))


def _format_totals(label: str, totals: CoverageTotals) -> str:
    percent = "n/a" if totals.attributed == 0 else f"{100 * totals.hit / totals.attributed:.3f}%"
    return (
        f"{label}: attributed={totals.attributed}, hit={totals.hit}, "
        f"missed={totals.missed}, coverage={percent}"
    )


def _format_keys(keys: Iterable[LineKey], limit: int = 20) -> str:
    ordered = sorted(keys)
    rendered = [f"{path}:{line}" for path, line in ordered[:limit]]
    if len(ordered) > limit:
        rendered.append(f"... and {len(ordered) - limit} more")
    return ", ".join(rendered) or "none"


def validate_artifact_tree(artifact_root: Path, expected_reports: Sequence[Path]) -> None:
    artifact_root = artifact_root.resolve()
    if not artifact_root.is_dir():
        raise CoverageError(f"missing coverage artifact directory: {artifact_root}")
    expected = {path.resolve() for path in expected_reports}
    if len(expected) != len(expected_reports):
        raise CoverageError("coverage report paths must be distinct")
    actual = {path.resolve() for path in artifact_root.rglob("*") if path.is_file()}
    missing = expected - actual
    extra = actual - expected
    if missing or extra:
        raise CoverageError(
            "coverage artifact tree does not contain exactly the expected reports; "
            f"missing={_format_keys((str(path), 0) for path in missing)}, "
            f"extra={_format_keys((str(path), 0) for path in extra)}"
        )


def compare_reports(
    repo_root: Path,
    artifact_root: Path,
    shard_paths: Sequence[Path],
    baseline_path: Path,
    base_sha: str,
) -> tuple[CoverageTotals, CoverageTotals]:
    if len(shard_paths) != 4:
        raise CoverageError(f"expected exactly four shard reports, received {len(shard_paths)}")
    reports = [*shard_paths, baseline_path]
    validate_artifact_tree(artifact_root, reports)
    shard_union = union_coverage(parse_lcov(path, repo_root) for path in shard_paths)
    baseline = parse_lcov(baseline_path, repo_root)
    changed = changed_lines(repo_root, base_sha)
    shard_core_keys = {key for key in shard_union if key[0].startswith("src/")}
    baseline_core_keys = {key for key in baseline if key[0].startswith("src/")}

    missing_keys = baseline.keys() - shard_union.keys()
    extra_keys = shard_union.keys() - baseline.keys()
    state_differences = {
        key for key in baseline.keys() & shard_union.keys()
        if baseline[key] != shard_union[key]
    }
    shard_project = coverage_totals(shard_union, shard_core_keys)
    baseline_project = coverage_totals(baseline, baseline_core_keys)
    shard_patch = coverage_totals(shard_union, changed)
    baseline_patch = coverage_totals(baseline, changed)

    print(_format_totals("shard project", shard_project))
    print(_format_totals("baseline project", baseline_project))
    print(_format_totals("shard patch", shard_patch))
    print(_format_totals("baseline patch", baseline_patch))

    if missing_keys or extra_keys or state_differences:
        raise CoverageError(
            "shard coverage differs from the transition baseline; "
            f"missing lines: {_format_keys(missing_keys)}; "
            f"extra lines: {_format_keys(extra_keys)}; "
            f"hit/miss differences: {_format_keys(state_differences)}"
        )
    if shard_project != baseline_project or shard_patch != baseline_patch:
        raise CoverageError("project or patch totals differ from the transition baseline")
    return shard_project, shard_patch


def _arguments(arguments: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, required=True)
    parser.add_argument("--artifact-root", type=Path, required=True)
    parser.add_argument("--base-sha", required=True)
    parser.add_argument("--shard", action="append", type=Path, required=True)
    parser.add_argument("--baseline", type=Path, required=True)
    return parser.parse_args(arguments)


def main(arguments: Sequence[str] | None = None) -> int:
    args = _arguments(arguments)
    try:
        compare_reports(
            args.repo_root.resolve(),
            args.artifact_root,
            args.shard,
            args.baseline,
            args.base_sha,
        )
    except CoverageError as error:
        print(f"coverage comparison failed: {error}", file=sys.stderr)
        return 1
    print("coverage comparison passed: shard and baseline line evidence are identical")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
