#!/usr/bin/env python3
"""Prepare or verify the pinned external source checkouts for matched runs."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path
from typing import Any


MANIFEST = Path(__file__).with_name("sources.json")


class SourceError(RuntimeError):
    """A source checkout does not match the benchmark manifest safely."""


def git(path: Path, *args: str) -> str:
    result = subprocess.run(
        ["git", "-C", str(path), *args],
        check=False,
        capture_output=True,
        text=True,
    )
    if result.returncode:
        detail = result.stderr.strip() or result.stdout.strip()
        raise SourceError(f"git {' '.join(args)} failed in {path}: {detail}")
    return result.stdout.strip()


def load_sources() -> list[dict[str, Any]]:
    try:
        manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise SourceError(f"cannot read source manifest {MANIFEST}: {exc}") from exc
    if manifest.get("schema_version") != 1 or not isinstance(manifest.get("sources"), list):
        raise SourceError(f"unsupported or invalid source manifest: {MANIFEST}")
    return manifest["sources"]


def checkout_status(path: Path) -> list[str]:
    raw = git(path, "status", "--porcelain=v1", "--untracked-files=all")
    lines = raw.splitlines()
    if lines:
        raise SourceError(
            f"checkout {path} has changes that prevent safe use:\n  "
            + "\n  ".join(lines)
        )
    return lines


def prepare_one(root: Path, source: dict[str, Any], verify_only: bool) -> None:
    name, url, revision, relative = (source.get(key) for key in ("name", "url", "revision", "path"))
    if not all(isinstance(value, str) and value for value in (name, url, revision, relative)):
        raise SourceError(f"invalid source entry: {source!r}")
    target = (root / relative).resolve()
    if not target.is_relative_to(root):
        raise SourceError(f"source path escapes --root: {relative}")

    if not target.exists():
        if verify_only:
            raise SourceError(f"{name}: checkout is missing at {target} (verify-only performs no network access)")
        target.parent.mkdir(parents=True, exist_ok=True)
        target.mkdir()
        git(target, "init")
        git(target, "remote", "add", "origin", url)
        git(target, "fetch", "--depth=1", "origin", revision)
        git(target, "checkout", "--detach", "FETCH_HEAD")
    elif not (target / ".git").exists():
        raise SourceError(f"{name}: target exists but is not a Git checkout: {target}")

    actual = git(target, "rev-parse", "HEAD")
    if actual != revision:
        raise SourceError(f"{name}: expected HEAD {revision}, found {actual} at {target}; checkout was not changed")
    checkout_status(target)
    print(f"{name}: verified {actual} at {target}; clean worktree")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True, help="directory under which source checkouts are stored")
    parser.add_argument("--verify-only", action="store_true", help="verify existing pinned checkouts without network access")
    args = parser.parse_args()
    try:
        root = args.root.expanduser().resolve()
        for source in load_sources():
            prepare_one(root, source, args.verify_only)
    except (OSError, SourceError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
