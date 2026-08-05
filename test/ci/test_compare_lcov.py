#!/usr/bin/env python3

from __future__ import annotations

import importlib.util
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


SCRIPT = Path(__file__).with_name("compare_lcov.py")
SPEC = importlib.util.spec_from_file_location("compare_lcov", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
compare_lcov = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = compare_lcov
SPEC.loader.exec_module(compare_lcov)


def write_lcov(path: Path, records: dict[str, dict[int, int]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines: list[str] = []
    for source, counts in records.items():
        lines.append(f"SF:{source}")
        lines.extend(f"DA:{line},{count}" for line, count in counts.items())
        lines.append(f"LF:{len(counts)}")
        lines.append(f"LH:{sum(count > 0 for count in counts.values())}")
        lines.append("end_of_record")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


class CompareLcovTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary_directory.name)
        (self.root / "src").mkdir()
        (self.root / "ext").mkdir()
        (self.root / "src" / "core.jl").write_text("a\nb\nc\n", encoding="utf-8")
        (self.root / "ext" / "gpu.jl").write_text("x\ny\n", encoding="utf-8")
        subprocess.run(["git", "init", "-q"], cwd=self.root, check=True)
        subprocess.run(["git", "config", "user.email", "ci@example.invalid"], cwd=self.root, check=True)
        subprocess.run(["git", "config", "user.name", "CI Test"], cwd=self.root, check=True)
        subprocess.run(["git", "add", "src", "ext"], cwd=self.root, check=True)
        subprocess.run(["git", "commit", "-qm", "baseline"], cwd=self.root, check=True)
        self.base_sha = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=self.root, check=True,
            capture_output=True, text=True,
        ).stdout.strip()
        (self.root / "src" / "core.jl").write_text("a\nchanged\nc\n", encoding="utf-8")
        subprocess.run(["git", "add", "src/core.jl"], cwd=self.root, check=True)
        subprocess.run(["git", "commit", "-qm", "change"], cwd=self.root, check=True)
        self.artifacts = self.root / "coverage"
        self.shards = [self.artifacts / f"shard-{index}" / "lcov.info" for index in range(4)]
        self.baseline = self.artifacts / "baseline" / "lcov.info"

    def tearDown(self) -> None:
        self.temporary_directory.cleanup()

    def write_equivalent_reports(self) -> None:
        write_lcov(self.shards[0], {"src/core.jl": {1: 1, 2: 0}})
        write_lcov(self.shards[1], {"src/core.jl": {2: 2, 3: 0}})
        write_lcov(self.shards[2], {"ext/gpu.jl": {1: 0}})
        write_lcov(self.shards[3], {"ext/gpu.jl": {2: 3}})
        write_lcov(self.baseline, {
            "src/core.jl": {1: 4, 2: 1, 3: 0},
            "ext/gpu.jl": {1: 0, 2: 1},
        })

    def test_equivalent_union_and_patch_totals_pass(self) -> None:
        self.write_equivalent_reports()
        project, patch = compare_lcov.compare_reports(
            self.root, self.artifacts, self.shards, self.baseline, self.base_sha
        )
        self.assertEqual(project, compare_lcov.CoverageTotals(3, 2))
        self.assertEqual(patch, compare_lcov.CoverageTotals(1, 1))

    def test_hit_miss_difference_fails(self) -> None:
        self.write_equivalent_reports()
        write_lcov(self.shards[1], {"src/core.jl": {2: 0, 3: 0}})
        with self.assertRaisesRegex(compare_lcov.CoverageError, "hit/miss differences"):
            compare_lcov.compare_reports(
                self.root, self.artifacts, self.shards, self.baseline, self.base_sha
            )

    def test_extra_artifact_file_fails_closed(self) -> None:
        self.write_equivalent_reports()
        (self.artifacts / "unexpected.txt").write_text("unexpected", encoding="utf-8")
        with self.assertRaisesRegex(compare_lcov.CoverageError, "exactly the expected"):
            compare_lcov.compare_reports(
                self.root, self.artifacts, self.shards, self.baseline, self.base_sha
            )

    def test_missing_artifact_fails_closed(self) -> None:
        self.write_equivalent_reports()
        self.shards[3].unlink()
        with self.assertRaisesRegex(compare_lcov.CoverageError, "exactly the expected"):
            compare_lcov.compare_reports(
                self.root, self.artifacts, self.shards, self.baseline, self.base_sha
            )

    def test_empty_report_fails_closed(self) -> None:
        self.write_equivalent_reports()
        self.shards[0].write_text("", encoding="utf-8")
        with self.assertRaisesRegex(compare_lcov.CoverageError, "empty LCOV report"):
            compare_lcov.compare_reports(
                self.root, self.artifacts, self.shards, self.baseline, self.base_sha
            )

    def test_zero_line_source_record_is_ignored(self) -> None:
        self.write_equivalent_reports()
        existing = self.shards[0].read_text(encoding="utf-8")
        self.shards[0].write_text(
            "SF:src/empty.jl\nLF:0\nLH:0\nend_of_record\n" + existing,
            encoding="utf-8",
        )
        project, patch = compare_lcov.compare_reports(
            self.root, self.artifacts, self.shards, self.baseline, self.base_sha
        )
        self.assertEqual(project, compare_lcov.CoverageTotals(3, 2))
        self.assertEqual(patch, compare_lcov.CoverageTotals(1, 1))

    def test_extension_change_does_not_enter_core_patch_totals(self) -> None:
        self.write_equivalent_reports()
        base_sha = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=self.root, check=True,
            capture_output=True, text=True,
        ).stdout.strip()
        (self.root / "ext" / "gpu.jl").write_text("changed\ny\n", encoding="utf-8")
        subprocess.run(["git", "add", "ext/gpu.jl"], cwd=self.root, check=True)
        subprocess.run(["git", "commit", "-qm", "extension change"], cwd=self.root, check=True)
        project, patch = compare_lcov.compare_reports(
            self.root, self.artifacts, self.shards, self.baseline, base_sha
        )
        self.assertEqual(project, compare_lcov.CoverageTotals(3, 2))
        self.assertEqual(patch, compare_lcov.CoverageTotals(0, 0))

    def test_malformed_report_fails_closed(self) -> None:
        self.write_equivalent_reports()
        self.shards[0].write_text("SF:src/core.jl\nDA:not-a-line,1\n", encoding="utf-8")
        with self.assertRaisesRegex(compare_lcov.CoverageError, "non-integer DA"):
            compare_lcov.compare_reports(
                self.root, self.artifacts, self.shards, self.baseline, self.base_sha
            )


if __name__ == "__main__":
    unittest.main()
