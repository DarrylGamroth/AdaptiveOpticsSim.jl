"""Create compact, checksummed archives for a completed matched-copper run."""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tomllib

HERE = Path(__file__).resolve().parent
CONFIGS = (
    "aos-cpu1", "aos-cpu8", "oopao-cpu1", "oopao-cpu8",
    "specula-cpu1", "aos-cuda", "oopao-cuda", "specula-cuda",
)
ENGINE_FOR_CONFIG = {key: key.split("-", 1)[0] for key in CONFIGS}
VALIDATION_FALLBACK = {
    "aos-cpu1": "aos-cpu-validate", "aos-cpu8": "validate-aos-cpu8-r1",
    "oopao-cpu1": "oopao-cpu-validate", "oopao-cpu8": "validate-oopao-cpu8-r1",
    "specula-cpu1": "specula-cpu-validate", "aos-cuda": "aos-cuda-validate",
    "oopao-cuda": "oopao-cuda-validate", "specula-cuda": "specula-cuda-validate",
}
ALLOWED_SOURCE_SUFFIXES = {".py", ".jl", ".toml", ".json", ".md", ".txt"}
CORE_REQUIRED = {
    "env-info.txt", "requirements-matched-copper-frozen.txt", "julia-cuda-version.txt",
    "source-state.json", "upstream-verification.json", "source-verification.txt",
    "system-before.json", "system-after.json", "telemetry.jsonl",
    "qualification-commands.json", "validation8-campaign.log", "benchmark-supervisor.log",
    "campaign-supervisor.log", "review.md", "specula-alignment.md",
    "evaluator-self-check-final.json",
}


def safe_file(path: Path) -> None:
    """Reject symlinks and non-regular files, including symlinked ancestors."""
    if path.is_symlink():
        raise ValueError(f"Refusing symlink: {path}")
    current = path.parent
    while current != current.parent:
        if current.is_symlink():
            raise ValueError(f"Refusing path through symlink: {current}")
        current = current.parent
    if not path.is_file():
        raise ValueError(f"Expected regular file: {path}")


def files_below(root: Path, relative: str) -> list[Path]:
    directory = root / relative
    if directory.is_symlink():
        raise ValueError(f"Refusing symlink directory: {directory}")
    if not directory.is_dir():
        raise ValueError(f"Missing required directory: {directory}")
    result = []
    for path in sorted(directory.rglob("*")):
        if path.is_symlink():
            raise ValueError(f"Refusing symlink: {path}")
        if path.is_file():
            safe_file(path)
            result.append(path)
    return result


def require_complete(root: Path, run_name: str) -> list[Path]:
    directory = root / run_name
    record_paths = [directory / "record.json", directory / "record.toml"]
    present = [path for path in record_paths if path.exists()]
    if len(present) != 1:
        raise ValueError(f"Expected exactly one JSON or TOML record in {directory}")
    record_path = present[0]
    safe_file(record_path)
    if record_path.suffix == ".json":
        record = json.loads(record_path.read_text())
    else:
        record = tomllib.loads(record_path.read_text())
    if record.get("complete") is not True:
        raise ValueError(f"Run is not marked complete: {record_path}")
    return files_below(root, run_name)


def add_core_matches(root: Path) -> list[Path]:
    selected = set()
    for name in CORE_REQUIRED:
        path = root / name
        safe_file(path)
        selected.add(path)
    for pattern in ("*-campaign.jsonl", "pip*.log", "independent-*.py", "independent-*.json",
                    "specula-native-cpu-preflight.*", "specula-cpu-preflight.*"):
        for path in root.glob(pattern):
            safe_file(path)
            selected.add(path)
    selected.update(files_below(root, "cold-factor-diagnostic"))
    env = root / "julia-env"
    for name in ("Project.toml", "Manifest.toml"):
        path = env / name
        safe_file(path)
        selected.add(path)
    return sorted(selected)


def validate_inputs(root: Path, output: Path) -> dict[str, list[Path]]:
    if not root.is_dir() or root.is_symlink():
        raise ValueError(f"Invalid evidence root: {root}")
    summary = root / "collected" / "summary.json"
    safe_file(summary)
    groups: dict[str, list[Path]] = {"analysis": files_below(root, "collected"),
                                    "fixture": files_below(root, "fixture")}
    validations = {}
    for config in CONFIGS:
        fresh = root / f"validate-{config}-r1"
        fallback = root / VALIDATION_FALLBACK[config]
        choice = fresh.name if fresh.exists() else fallback.name
        validations[config] = choice
        require_complete(root, choice)
    for engine in ("aos", "oopao", "specula"):
        selected: set[Path] = set()
        for config in CONFIGS:
            if ENGINE_FOR_CONFIG[config] != engine:
                continue
            validation = validations[config]
            selected.update(require_complete(root, validation))
            validation_log = root / f"{validation}.log"
            safe_file(validation_log)
            selected.add(validation_log)
            for round_number in range(1, 4):
                run = f"benchmark-{config}-r{round_number}"
                selected.update(require_complete(root, run))
                log = root / f"{run}.log"
                safe_file(log)
                selected.add(log)
        groups[f"raw-{engine}"] = sorted(selected)

    core = add_core_matches(root)
    harness = []
    for path in sorted(HERE.rglob("*")):
        if path.is_symlink():
            raise ValueError(f"Refusing harness symlink: {path}")
        if path.is_file() and path.suffix in ALLOWED_SOURCE_SUFFIXES:
            safe_file(path)
            harness.append(path)
    if not harness:
        raise ValueError("No harness source files found")
    groups["core"] = core + harness

    output_resolved = output.resolve()
    for group, paths in groups.items():
        for path in paths:
            resolved = path.resolve(strict=True)
            if path.is_symlink() or (resolved != root and root not in resolved.parents
                                     and resolved != HERE and HERE not in resolved.parents):
                raise ValueError(f"Selected path escapes its permitted root: {path}")
            if resolved == output_resolved or output_resolved in resolved.parents:
                raise ValueError(f"Output directory cannot be archived: {path}")
    return groups


def sha256_file(path: Path) -> tuple[str, int]:
    digest = hashlib.sha256()
    size = 0
    with path.open("rb") as stream:
        while chunk := stream.read(1024 * 1024):
            digest.update(chunk)
            size += len(chunk)
    return digest.hexdigest(), size


def archive_group(name: str, paths: list[Path], root: Path, output: Path) -> Path:
    target = output / f"{name}.tar.zst"
    if target.exists():
        raise FileExistsError(f"Refusing to overwrite archive: {target}")
    command = ["tar", "--zstd", "-cf", str(target)]
    command.extend(["--transform=s,^benchmarks/matched_copper/,source-harness/benchmarks/matched_copper/,"])
    for path in paths:
        safe_file(path)
        base = root if root == path or root in path.parents else HERE.parent.parent
        relative = path.relative_to(base)
        command.extend(["-C", str(base), str(relative)])
    try:
        subprocess.run(command, check=True)
    except Exception:
        target.unlink(missing_ok=True)
        raise
    size = target.stat().st_size
    if size >= 90 * 1024 * 1024:
        target.unlink()
        raise ValueError(f"Archive exceeds 90 MiB: {target} ({size} bytes)")
    return target


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path, help="Completed evidence root")
    parser.add_argument("--output", required=True, type=Path, help="Results directory")
    args = parser.parse_args()
    root = args.root.resolve()
    output = args.output.resolve()
    groups = validate_inputs(root, output)
    output.mkdir(parents=True, exist_ok=True)
    archives = [output / f"{name}.tar.zst" for name in
                ("raw-aos", "raw-oopao", "raw-specula", "analysis", "fixture", "core")]
    for index in (output / "raw-sha256.json", output / "archives-sha256.json"):
        if index.exists():
            raise FileExistsError(f"Refusing to overwrite index: {index}")
    if any(path.exists() for path in archives):
        existing = next(path for path in archives if path.exists())
        raise FileExistsError(f"Refusing to overwrite archive: {existing}")

    raw_index = []
    for group, paths in groups.items():
        for path in paths:
            digest, size = sha256_file(path)
            origin = "evidence" if root == path or root in path.parents else "harness"
            relative = path.relative_to(root if origin == "evidence" else HERE.parent.parent)
            archive_path = (f"source-harness/{relative}" if origin == "harness" else str(relative))
            raw_index.append(dict(path=f"{origin}/{relative}", archive_path=archive_path,
                                  size_bytes=size,
                                  sha256=digest, archive_group=group))
    raw_index.sort(key=lambda item: (item["archive_group"], item["path"]))
    (output / "raw-sha256.json").write_text(json.dumps(raw_index, indent=2) + "\n")

    created = []
    try:
        for name in ("raw-aos", "raw-oopao", "raw-specula", "analysis", "fixture", "core"):
            created.append(archive_group(name, groups[name], root, output))
        archive_index = []
        for path in created:
            digest, size = sha256_file(path)
            archive_index.append(dict(path=path.name, size_bytes=size, sha256=digest))
        (output / "archives-sha256.json").write_text(json.dumps(archive_index, indent=2) + "\n")
    except Exception:
        for path in created:
            path.unlink(missing_ok=True)
        (output / "raw-sha256.json").unlink(missing_ok=True)
        raise
    print(f"Created {len(created)} archives in {output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
