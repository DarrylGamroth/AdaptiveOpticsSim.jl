"""Cross-language test peer for the native Julia pyRTC SHM adapter."""

from __future__ import annotations

import argparse
import sys
from multiprocessing import shared_memory
from pathlib import Path

import numpy as np


CONTROL_PREFIX = "AOS_PYRTC_PEER "


def control(message: str) -> None:
    print(CONTROL_PREFIX + message, flush=True)


def expected_values(shape: tuple[int, ...]) -> np.ndarray:
    if len(shape) == 1:
        return np.arange(shape[0], dtype=np.float32) + np.float32(0.25)
    if len(shape) == 2:
        rows, columns = shape
        row_values = np.arange(rows, dtype=np.float32).reshape(rows, 1)
        column_values = np.arange(columns, dtype=np.float32).reshape(1, columns)
        return row_values * np.float32(100) + column_values + np.float32(0.25)
    raise ValueError("the interoperability peer supports vectors and matrices")


def close_stream(stream, *, unlink: bool) -> None:
    payload_name = stream.shm._name
    metadata_name = (
        stream.metadataShm._name if hasattr(stream, "metadataShm") else None
    )
    stream.shm.close()
    if hasattr(stream, "metadataShm"):
        stream.metadataShm.close()
    if unlink:
        for name in (payload_name, metadata_name):
            if name is None:
                continue
            try:
                shared_memory._posixshmem.shm_unlink(name)
            except FileNotFoundError:
                pass


def consume(name: str, shape: tuple[int, ...]) -> None:
    from pyRTC.Pipeline import initExistingShm

    stream, discovered_shape, dtype = initExistingShm(name)
    try:
        if tuple(discovered_shape) != shape:
            raise ValueError(
                f"discovered shape {tuple(discovered_shape)} does not match {shape}"
            )
        if np.dtype(dtype) != np.dtype(np.float32):
            raise TypeError(f"discovered dtype {dtype} is not float32")
        control("READY")
        actual = stream.read()
        np.testing.assert_array_equal(actual, expected_values(shape))
        control("DONE")
    finally:
        close_stream(stream, unlink=False)


def produce(name: str, shape: tuple[int, ...]) -> None:
    from pyRTC.Pipeline import ImageSHM

    stream = ImageSHM(name, shape, np.float32, consumer=False)
    try:
        control("READY")
        command = sys.stdin.readline().strip()
        if command != "PUBLISH":
            raise RuntimeError(f"expected PUBLISH, received {command!r}")
        if stream.write(expected_values(shape)) != 1:
            raise RuntimeError("pyRTC rejected the test publication")
        control("PUBLISHED")
        command = sys.stdin.readline().strip()
        if command != "STOP":
            raise RuntimeError(f"expected STOP, received {command!r}")
        control("DONE")
    finally:
        close_stream(stream, unlink=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("consume", "produce"))
    parser.add_argument("name")
    parser.add_argument("shape", nargs="+", type=int)
    parser.add_argument("--pyrtc-root", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    pyrtc_root = Path(args.pyrtc_root).resolve()
    if not (pyrtc_root / "pyRTC").is_dir():
        raise FileNotFoundError(f"{pyrtc_root} does not contain pyRTC/")
    sys.path.insert(0, str(pyrtc_root))
    shape = tuple(args.shape)
    if args.mode == "consume":
        consume(args.name, shape)
    else:
        produce(args.name, shape)


if __name__ == "__main__":
    main()
