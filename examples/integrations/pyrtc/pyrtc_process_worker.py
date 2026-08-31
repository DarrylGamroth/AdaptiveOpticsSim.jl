"""Run pyRTC signal processing and control outside the Julia process."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np


CONTROL_PREFIX = "AOS_PYRTC_WORKER "


def control(message: str) -> None:
    print(CONTROL_PREFIX + message, flush=True)


def python_config(**entries):
    return entries


def shack_hartmann_valid_subapertures() -> np.ndarray:
    count = 8
    radius = count / 2
    centre = (count + 1) / 2
    lenslet_axis = (np.arange(count, dtype=np.float64) + 1 - centre) / radius
    axis_1, axis_2 = np.meshgrid(lenslet_axis, lenslet_axis, indexing="ij")
    valid = axis_1 * axis_1 + axis_2 * axis_2 <= 1
    return np.concatenate((valid, valid), axis=0)


def slopes_config(sensor: str, temporary_directory: Path):
    if sensor == "shack_hartmann":
        valid_path = temporary_directory / "valid_subapertures.npy"
        np.save(valid_path, shack_hartmann_valid_subapertures())
        return python_config(
            type="SHWFS",
            signalType="slopes",
            subApSpacing=8,
            subApOffsetX=0,
            subApOffsetY=0,
            imageNoise=0.0,
            contrast=0,
            validSubApsFile=str(valid_path),
        )
    if sensor == "pyramid":
        return python_config(
            type="PYWFS",
            signalType="slopes",
            imageNoise=0.0,
            centralObscurationRatio=0.0,
            flatNorm=True,
            pupils=["9,9", "27,9", "9,27", "27,27"],
            pupilsRadius=8,
        )
    raise ValueError(f"unsupported sensor {sensor!r}")


def close_handle(handle) -> None:
    handle.shm.close()
    if hasattr(handle, "metadataShm"):
        handle.metadataShm.close()


def close_components(slopes, loop) -> None:
    handles = (
        loop.signalShm,
        loop.wfcShm,
        slopes.wfsShm,
        slopes.signal,
        slopes.signal2D,
    )
    for handle in handles:
        try:
            close_handle(handle)
        except Exception:
            pass


def configure_loop(loop, matrix_path: str, gain: float) -> None:
    matrix = np.fromfile(matrix_path, dtype=np.float32)
    expected = loop.signalSize * loop.numModes
    if matrix.size != expected:
        raise ValueError(
            f"interaction matrix contains {matrix.size} values; expected {expected}"
        )
    loop.IM = matrix.reshape((loop.signalSize, loop.numModes), order="F")
    loop.computeCM()
    loop.setGain(gain)


def serve(sensor: str, temporary_directory: Path) -> None:
    from pyRTC.Loop import Loop
    from pyRTC.SlopesProcess import SlopesProcess

    slopes = SlopesProcess(slopes_config(sensor, temporary_directory))
    loop = Loop(python_config(gain=0.4, numDroppedModes=0))
    try:
        control(f"READY {loop.signalSize} {loop.numModes}")
        for line in sys.stdin:
            command = line.strip().split()
            if not command:
                continue
            if command[0] == "PROCESS" and len(command) == 1:
                slopes.computeSignal()
                control("PROCESSED")
            elif command[0] == "SET_REF" and len(command) == 1:
                slopes.setRefSlopes(slopes.signal2D.read_noblock())
                control("REF_SET")
            elif command[0] == "CONFIGURE" and len(command) == 3:
                configure_loop(loop, command[1], float(command[2]))
                control("CONFIGURED")
            elif command[0] == "FLATTEN" and len(command) == 1:
                loop.flatten()
                control("FLATTENED")
            elif command[0] == "STEP" and len(command) == 1:
                slopes.computeSignal()
                loop.standardIntegrator()
                control("STEPPED")
            elif command[0] == "STOP" and len(command) == 1:
                control("STOPPED")
                return
            else:
                raise ValueError(f"invalid worker command: {line.strip()!r}")
    finally:
        close_components(slopes, loop)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--sensor",
        required=True,
        choices=("shack_hartmann", "pyramid"),
    )
    parser.add_argument("--temporary-directory", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    temporary_directory = Path(args.temporary_directory).resolve()
    temporary_directory.mkdir(parents=True, exist_ok=True)
    serve(args.sensor, temporary_directory)


if __name__ == "__main__":
    main()
