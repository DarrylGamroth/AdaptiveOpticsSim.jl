#!/usr/bin/env python3
"""Import a typed YAML/FITS model or frozen TOML/FITS export as portable YAML.

The cold tool auto-detects model.yaml or historical manifest.toml. Runtime readers use typed YAML
array descriptors. --link preserves frozen generated array bytes without copying
large data; neither source metadata nor array files are modified.
--preserve-metadata retains canonical prepared YAML bytes and their provenance.
"""
import argparse
import copy
import hashlib
import os
import pathlib
import shutil
import tomllib


def digest(path):
    with open(path, "rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def contained(root, relative):
    path = pathlib.Path(relative)
    if str(path) != relative or path.is_absolute() or any(part in (".", "..") for part in path.parts):
        raise ValueError("Array references must stay inside the model directory")
    result = (root / path).resolve(strict=True)
    result.relative_to(root)
    return result


def model_document(metadata, output, source_root=None, link=False, include_references=False):
    import numpy as np
    from astropy.io import fits
    result = copy.deepcopy(metadata)
    if not include_references:
        result.pop("references", None)
    provenance = result.setdefault("provenance", {})
    provenance.pop("source_fit", None)
    provenance.pop("model_path", None)
    typed = metadata.get("arrays")
    if typed is not None and (not isinstance(typed, dict) or not typed):
        raise ValueError("Compiled YAML requires typed array descriptors")
    references = typed or metadata["array_sha256"]
    prepared = []
    descriptors = {}
    # Validate every referenced payload before staging any file.
    for reference in references:
        declared = typed[reference] if typed else None
        if declared is not None:
            required = {"path", "sha256", "shape", "type", "axes", "real_hdu", "imag_hdu"}
            if set(declared) != required:
                raise ValueError("Unknown or missing typed array descriptor fields")
            relative, expected = declared["path"], declared["sha256"]
            axes = declared["axes"]
            if axes not in ([1, 2], [2, 1]) or any(type(axis) is not int for axis in axes):
                raise ValueError("Optical array axes must explicitly permute[1,2]")
            shape = declared["shape"]
            if not isinstance(shape, list) or len(shape) != 2 or any(type(n) is not int or n <= 0 for n in shape):
                raise ValueError("Optical array shape must contain two positive integers")
            real_hdu, imag_hdu = declared["real_hdu"], declared["imag_hdu"]
            if type(real_hdu) is not int or type(imag_hdu) is not int or real_hdu < 1 or imag_hdu < 0:
                raise ValueError("Invalid optical array HDU indices")
        else:
            relative, expected = reference, metadata["array_sha256"][reference]
            axes, real_hdu, imag_hdu = [2, 1], 1, None
        source = contained(source_root or output, relative)
        if digest(source) != expected:
            raise ValueError(f"Frozen array SHA-256 mismatch: {reference}")
        with fits.open(source, memmap=True) as hdus:
            if imag_hdu is None:
                if len(hdus) not in (1, 2):
                    raise ValueError("Legacy optical arrays require one or two HDUs")
                imag_hdu = 2 if len(hdus) == 2 else 0
            if real_hdu > len(hdus) or imag_hdu > len(hdus):
                raise ValueError("Optical array HDU is missing")
            real = hdus[real_hdu - 1].data
            if real.ndim != 2:
                raise ValueError("Optical arrays must be matrices")
            shape = list(real.shape if axes == [2, 1] else real.T.shape)
            bits = int(hdus[real_hdu - 1].header["BITPIX"])
            if bits not in (-32, -64) or not np.isfinite(real).all():
                raise ValueError("Optical arrays require finite floating-point FITS images")
            if imag_hdu:
                imaginary = hdus[imag_hdu - 1].data
                if imaginary.shape != real.shape or hdus[imag_hdu - 1].header["BITPIX"] != bits or not np.isfinite(imaginary).all():
                    raise ValueError("Complex optical array components differ")
            element = ("ComplexF" if imag_hdu else "Float") + str(-bits)
        if declared is not None and (declared["shape"] != shape or declared["type"] != element):
            raise ValueError("Typed optical array shape or element type differs")
        destination_relative = relative if typed or source_root is None else "arrays/" + relative
        descriptor = dict(path=destination_relative, sha256=expected, shape=shape,
                          type=element, axes=axes, real_hdu=real_hdu, imag_hdu=imag_hdu)
        descriptors[reference] = descriptor
        prepared.append((source, destination_relative, expected))
    if source_root is not None:
        for source, relative, expected in prepared:
            destination = output / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            if destination.exists():
                if digest(destination) != expected:
                    raise ValueError(f"Destination array differs: {relative}")
            elif link:
                os.link(source, destination)
            else:
                shutil.copyfile(source, destination)
    result["arrays"] = descriptors
    for path in result["propagations"]:
        for section in ("prefix", "closed", "only"):
            for operation in path[section]:
                if operation["kind"] == "live_screen" and "name" in operation:
                    name = operation.pop("name")
                    operation["input"] = {"DM": "dm_opd", "Subaru_Pupil": "pupil_opd"}[name]
    return result


def read_model_source(source):
    if source.is_file():
        manifest = source
    elif (source / "model.yaml").is_file():
        manifest = source / "model.yaml"
    else:
        manifest = source / "manifest.toml"
    if manifest.suffix in (".yaml", ".yml"):
        import yaml
        metadata = yaml.safe_load(manifest.read_text())
        if not isinstance(metadata, dict) or metadata.get("schema_version") != 1 or "arrays" not in metadata:
            raise ValueError("Unsupported compiled YAML model")
    elif manifest.suffix == ".toml":
        metadata = tomllib.loads(manifest.read_text())
    else:
        raise ValueError("Compiled source must be YAML metadata or a historical TOML export")
    return metadata, manifest


def write_model(document, destination):
    import yaml
    with open(destination, "w", encoding="utf8") as stream:
        stream.write("# Generated by AdaptiveOpticsSim tools/export_optical_paths.py or tools/import_compiled_optical_model.py.\n")
        stream.write("# Source identities/options are recorded in provenance; see docs/proper-integration-guide.md.\n")
        yaml.safe_dump(document, stream, sort_keys=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=pathlib.Path)
    parser.add_argument("output", type=pathlib.Path)
    parser.add_argument("--link", action="store_true")
    parser.add_argument("--name", default="model.yaml")
    parser.add_argument("--preserve-metadata", action="store_true",
                        help="Validate and copy canonical prepared YAML without changing metadata bytes")
    args = parser.parse_args()
    name = pathlib.Path(args.name)
    if (name.name != args.name or name.suffix not in (".yaml", ".yml")
            or "\\" in args.name):
        parser.error("Metadata name must be a single relative YAML filename")
    source = args.source.resolve(strict=True)
    output = args.output.resolve()
    metadata, manifest = read_model_source(source)
    if args.preserve_metadata and manifest.suffix not in (".yaml", ".yml"):
        raise ValueError("Metadata preservation requires a typed canonical YAML model")
    source_root = manifest.parent
    if source_root == output:
        parser.error("Import output must differ from the frozen source directory")
    destination = output / args.name
    if args.preserve_metadata and destination.is_symlink():
        raise ValueError("Metadata destination must not be a symbolic link")
    output.mkdir(parents=True, exist_ok=True)
    document = model_document(metadata, output, source_root=source_root, link=args.link)
    if args.preserve_metadata:
        if document != metadata:
            raise ValueError("Metadata preservation requires canonical YAML without normalization")
        if destination.exists():
            if digest(destination) != digest(manifest):
                raise ValueError("Destination metadata differs from the frozen source")
        else:
            shutil.copyfile(manifest, destination)
    else:
        identity = "imported_manifest_sha256" if manifest.suffix == ".toml" else "imported_model_sha256"
        document["provenance"][identity] = digest(manifest)
        document["provenance"]["imported_metadata_sha256"] = digest(manifest)
        document["provenance"]["importer_sha256"] = digest(pathlib.Path(__file__))
        document["provenance"]["array_preservation"] = "source generated FITS bytes retained exactly"
        write_model(document, destination)
    print(f"Imported {len(document['arrays'])} typed optical arrays into {destination}")


if __name__ == "__main__":
    main()
