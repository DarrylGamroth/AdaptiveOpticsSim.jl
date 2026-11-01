"""Small cold I/O checks; no JAX, optical propagation or device imports."""
import copy
import contextlib
import io
import pathlib
import sys
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
from astropy.io import fits
import import_compiled_optical_model as loader
import export_optical_paths as compiler


class Artifacts(unittest.TestCase):
    def run_import(self, source, output, *options):
        with patch.object(sys, "argv", ["import_compiled_optical_model.py", str(source), str(output), *options]), contextlib.redirect_stdout(io.StringIO()):
            loader.main()

    def prepared_model(self, root):
        arrays = compiler.Arrays(root / "arrays")
        arrays.save(np.arange(6.0).reshape(2, 3))
        metadata = dict(schema_version=1, array_sha256=arrays.sha256, propagations=[],
                        provenance={"imported_metadata_sha256": "a" * 64, "importer_sha256": "b" * 64})
        metadata = loader.model_document(metadata, root)
        loader.write_model(metadata, root / "model.yaml")
        # Include authored formatting that a YAML rewrite would lose.
        with (root / "model.yaml").open("a") as stream:
            stream.write("\n# Frozen prepared model: retain these exact bytes.\n")
        return metadata

    def test_preserve_prepared_metadata_and_arrays(self):
        with tempfile.TemporaryDirectory() as directory:
            root = pathlib.Path(directory) / "source"
            root.mkdir()
            metadata = self.prepared_model(root)
            before = {path.relative_to(root): path.read_bytes() for path in root.rglob("*") if path.is_file()}
            for link in (False, True):
                target = pathlib.Path(directory) / str(link)
                options = ["--preserve-metadata", "--name", "prepared.yml"] + (["--link"] if link else [])
                self.run_import(root, target, *options)
                self.run_import(root / "model.yaml", target, *options)
                self.assertEqual((target / "prepared.yml").read_bytes(), before[pathlib.Path("model.yaml")])
                for descriptor in metadata["arrays"].values():
                    relative = pathlib.Path(descriptor["path"])
                    self.assertEqual((target / relative).read_bytes(), before[relative])
                    if link:
                        self.assertEqual((target / relative).stat().st_ino, (root / relative).stat().st_ino)
                self.assertEqual(before, {path.relative_to(root): path.read_bytes() for path in root.rglob("*") if path.is_file()})

    def test_preserve_rejects_invalid_or_noncanonical_models(self):
        with tempfile.TemporaryDirectory() as source, tempfile.TemporaryDirectory() as output:
            root, target = pathlib.Path(source), pathlib.Path(output)
            metadata = self.prepared_model(root)
            changes = [lambda m: m["arrays"][next(iter(m["arrays"]))].update(sha256="0" * 64),
                       lambda m: m["arrays"][next(iter(m["arrays"]))].update(shape=[3, 2]),
                       lambda m: m["arrays"][next(iter(m["arrays"]))].update(axes=[1, 1]),
                       lambda m: m["arrays"][next(iter(m["arrays"]))].update(path="../escape.fits"),
                       lambda m: m.update(arrays={}),
                       lambda m: m["provenance"].update(source_fit="/origin/fit"),
                       lambda m: m.update(references={"unused": "reference"}),
                       lambda m: m.update(propagations=[dict(prefix=[dict(kind="live_screen", name="DM")], closed=[], only=[])])]
            for change in changes:
                invalid = copy.deepcopy(metadata)
                change(invalid)
                loader.write_model(invalid, root / "model.yaml")
                with self.assertRaises(ValueError):
                    self.run_import(root, target, "--preserve-metadata")
                self.assertFalse((target / "model.yaml").exists())
            (root / "manifest.toml").write_text("schema_version = 1\npropagations = []\n[array_sha256]\n")
            with self.assertRaises(ValueError):
                self.run_import(root / "manifest.toml", target, "--preserve-metadata")
            (root / "model.yaml").write_text("schema_version: 1\npropagations: []\n")
            with self.assertRaises(ValueError):
                self.run_import(root, target, "--preserve-metadata")

    def test_preserve_rejects_unsafe_metadata_destinations(self):
        with tempfile.TemporaryDirectory() as directory:
            root = pathlib.Path(directory) / "source"
            root.mkdir()
            self.prepared_model(root)
            before = (root / "model.yaml").read_bytes()
            target = pathlib.Path(directory) / "output"
            for name in ("../escape.yaml", "/absolute.yaml", "nested/model.yaml", "model.toml", ""):
                with self.subTest(name=name), contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                    self.run_import(root, target, "--preserve-metadata", "--name", name)
            target.mkdir(exist_ok=True)
            (target / "model.yaml").symlink_to(root / "model.yaml")
            for options in ([], ["--preserve-metadata"]):
                with self.assertRaises(ValueError):
                    self.run_import(root, target, *options)
            self.assertEqual((root / "model.yaml").read_bytes(), before)
            (target / "model.yaml").unlink()
            payload = next((root / "arrays").iterdir())
            payload_before = payload.read_bytes()
            (target / "model.yaml").hardlink_to(payload)
            for options in ([], ["--preserve-metadata"]):
                with self.assertRaises(ValueError):
                    self.run_import(root, target, *options)
            self.assertEqual(payload.read_bytes(), payload_before)

    def test_preserve_rejects_duplicate_and_nonstring_yaml_keys_before_staging(self):
        with tempfile.TemporaryDirectory() as directory:
            root = pathlib.Path(directory) / "source"
            root.mkdir()
            self.prepared_model(root)
            manifest = root / "model.yaml"
            canonical = manifest.read_text()
            invalid_documents = [canonical + "schema_version: 1\n",
                                 canonical + "extra:\n  duplicate: 1\n  duplicate: 1\n",
                                 canonical + "1: nonstring\n",
                                 canonical + "extra:\n  true: nonstring\n",
                                 canonical + "extra:\n  <<: {duplicate: 1, duplicate: 1}\n",
                                 canonical + "extra:\n  <<: {1: nonstring}\n"]
            for index, text in enumerate(invalid_documents):
                with self.subTest(index=index):
                    manifest.write_text(text)
                    target = pathlib.Path(directory) / str(index)
                    with self.assertRaisesRegex(ValueError, "duplicate YAML mapping key|YAML mapping keys must be strings"):
                        self.run_import(root, target, "--preserve-metadata")
                    self.assertFalse(target.exists())
                    self.assertEqual(manifest.read_text(), text)
            # The ordinary importer keeps its existing YAML parsing behavior.
            manifest.write_text(invalid_documents[0])
            self.assertEqual(loader.read_model_source(root)[0]["schema_version"], 1)

    def test_preserve_allows_yaml_merge_defaults_and_overrides(self):
        with tempfile.TemporaryDirectory() as source, tempfile.TemporaryDirectory() as output:
            root, target = pathlib.Path(source), pathlib.Path(output)
            self.prepared_model(root)
            manifest = root / "model.yaml"
            with manifest.open("a") as stream:
                stream.write("extra:\n  base: &base\n    <<: {value: default}\n    value: override\n  reused:\n    <<: *base\n")
            self.run_import(root, target, "--preserve-metadata")
            self.assertEqual(manifest.read_bytes(), (target / "model.yaml").read_bytes())
            parsed, _ = loader.read_model_source(root, strict=True)
            self.assertEqual(parsed["extra"], {"base": {"value": "override"}, "reused": {"value": "override"}})

    def test_default_legacy_import_still_normalizes(self):
        with tempfile.TemporaryDirectory() as source, tempfile.TemporaryDirectory() as output:
            root, target = pathlib.Path(source), pathlib.Path(output)
            fits.PrimaryHDU(np.ones((2, 3))).writeto(root / "a.fits")
            before = loader.digest(root / "a.fits")
            (root / "manifest.toml").write_text(f'schema_version = 1\npropagations = []\n[array_sha256]\n"a.fits" = "{before}"\n')
            self.run_import(root, target)
            metadata, _ = loader.read_model_source(target)
            self.assertEqual(metadata["arrays"]["a.fits"]["path"], "arrays/a.fits")
            self.assertEqual(metadata["provenance"]["imported_manifest_sha256"], loader.digest(root / "manifest.toml"))
            self.assertEqual(loader.digest(target / "arrays/a.fits"), before)

    def test_frozen_link_descriptor_and_sanitization(self):
        with tempfile.TemporaryDirectory() as source, tempfile.TemporaryDirectory() as output:
            root, target = pathlib.Path(source), pathlib.Path(output)
            values = np.arange(6.0).reshape(2, 3) + 1j * np.arange(6.0, 12.0).reshape(2, 3)
            path = root / "asymmetric.fits"
            fits.HDUList([fits.PrimaryHDU(values.real), fits.ImageHDU(values.imag)]).writeto(path)
            before = loader.digest(path)
            metadata = dict(array_sha256={path.name: before}, provenance={"source_fit": "/origin/fit", "model_path": "/origin/model"},
                            propagations=[dict(prefix=[dict(kind="live_screen", name="DM")], closed=[], only=[])])
            result = loader.model_document(metadata, target, root, link=True)
            descriptor = result["arrays"][path.name]
            self.assertEqual(descriptor["shape"], [2, 3])
            self.assertEqual(descriptor["type"], "ComplexF64")
            self.assertEqual(descriptor["axes"], [2, 1])
            self.assertEqual(descriptor["sha256"], before)
            self.assertEqual(path.stat().st_ino, (target / descriptor["path"]).stat().st_ino)
            self.assertNotIn("source_fit", result["provenance"])
            self.assertNotIn("model_path", result["provenance"])
            self.assertEqual(result["propagations"][0]["prefix"][0]["input"], "dm_opd")
            self.assertEqual(loader.digest(path), before)
            self.assertEqual(metadata["propagations"][0]["prefix"][0]["name"], "DM")

    def test_hash_and_path_rejection(self):
        with tempfile.TemporaryDirectory() as source, tempfile.TemporaryDirectory() as output:
            root, target = pathlib.Path(source), pathlib.Path(output)
            fits.PrimaryHDU(np.ones((2, 2))).writeto(root / "a.fits")
            metadata = dict(array_sha256={"a.fits": "0" * 64}, propagations=[])
            with self.assertRaises(ValueError):
                loader.model_document(metadata, target, root)
            for path in ("../escape.fits", "/absolute.fits"):
                with self.assertRaises(ValueError):
                    loader.contained(root, path)

    def test_legacy_manifest_detection(self):
        with tempfile.TemporaryDirectory() as source:
            root = pathlib.Path(source)
            (root / "manifest.toml").write_text("schema_version = 1\npropagations = []\n[array_sha256]\n")
            metadata, path = loader.read_model_source(root)
            self.assertEqual(path.name, "manifest.toml")
            self.assertEqual(metadata["schema_version"], 1)

    def test_typed_yaml_roundtrip_and_bad_descriptor(self):
        with tempfile.TemporaryDirectory() as source, tempfile.TemporaryDirectory() as output:
            root, target = pathlib.Path(source), pathlib.Path(output)
            arrays = compiler.Arrays(root / "arrays")
            reference = arrays.save(np.arange(6.0).reshape(2, 3))
            metadata = dict(schema_version=1, array_sha256=arrays.sha256, propagations=[], provenance={})
            typed = loader.model_document(metadata, root)
            loader.write_model(typed, root / "model.yaml")
            parsed, manifest = loader.read_model_source(root)
            self.assertEqual(manifest.name, "model.yaml")
            imported = loader.model_document(parsed, target, root)
            descriptor = imported["arrays"][reference]
            self.assertEqual(loader.digest(target / descriptor["path"]), descriptor["sha256"])
            self.assertEqual((root / descriptor["path"]).read_bytes(), (target / descriptor["path"]).read_bytes())
            for key, value in (("sha256", "0" * 64), ("path", "../escape.fits"), ("shape", [3, 2])):
                original = parsed["arrays"][reference][key]
                parsed["arrays"][reference][key] = value
                with self.assertRaises(ValueError):
                    loader.model_document(parsed, target, root)
                parsed["arrays"][reference][key] = original

    def test_fresh_export_headers_and_deduplication(self):
        with tempfile.TemporaryDirectory() as output:
            root = pathlib.Path(output)
            arrays = compiler.Arrays(root / "arrays")
            values = np.array([[1.0 + 2j, 3.0 + 4j], [5.0 + 6j, 7.0 + 8j]])
            reference = arrays.save(values)
            self.assertEqual(arrays.save(values.copy()), reference)
            self.assertEqual(len(arrays.sha256), 1)
            with fits.open(root / reference) as hdus:
                self.assertIn("GENERATE", hdus[0].header)
                self.assertIn("SOURCE", hdus[0].header)
                np.testing.assert_array_equal(hdus[0].data + 1j * hdus[1].data, values)
            before = loader.digest(root / reference)
            # Reuse a previously frozen array without rewriting its header/data.
            second = compiler.Arrays(root / "arrays")
            self.assertEqual(second.save(values), reference)
            self.assertEqual(loader.digest(root / reference), before)


if __name__ == "__main__":
    unittest.main()
