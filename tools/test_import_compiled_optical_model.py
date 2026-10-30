"""Small cold I/O checks; no JAX, optical propagation or device imports."""
import hashlib
import pathlib
import tempfile
import unittest
import numpy as np
from astropy.io import fits
import import_compiled_optical_model as loader
import export_optical_paths as compiler


class Artifacts(unittest.TestCase):
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
