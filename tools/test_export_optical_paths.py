# Adapted from SPIDERSSim b479871 cold compiler tests.
"""Cheap exporter checks: synthetic inputs only, no JAX or propagation imports."""

import ast
import contextlib
import dataclasses
import io
import pathlib
import types
import unittest

import numpy as np

import export_optical_paths as compiler


def options(*flags):
    parser = compiler.argument_parser()
    return compiler.resolve_options(parser.parse_args(["--output", "/unused", *flags]), parser)


def synthetic_manifest():
    phantoms = [dict(label="stack0.724", source="main", opd=0.000724),
                dict(label="custom", source="main", opd=0.0009),
                dict(label="ghost0.100", source=None, opd=0.0001)]
    return dict(phantoms=phantoms, labels=["main", "lens", *[p["label"] for p in phantoms]],
                screens=["DM", "L2B", "Apodizer"],
                fast=dict(phantom_opd_um=np.array([1., 2., 3.]),
                          phantom_defocus_mm=np.array([4., 5., 6.]),
                          lowphase=np.arange(25).reshape(5, 5)))


class Arguments(unittest.TestCase):
    def test_defaults_and_presets(self):
        for preset, expected in compiler.PRESET_OPTIONS.items():
            args = options("--preset", preset)
            for key, value in expected.items():
                self.assertEqual(args.resolved[key], value)
            self.assertEqual(args.resolved["channels"], 8)
            self.assertEqual(args.configuration, "reference")
            self.assertEqual(args.mode, "sky")
            self.assertEqual(args.resolved["fixed_ghost"], "fitted")

    def test_checked_out_source_preset_parity(self):
        source = compiler.argument_parser().get_default("spiders_sim") / "src/spiders_sim/sim.py"
        if not source.exists():
            self.skipTest("No sibling spiders-sim checkout")
        tree = ast.parse(source.read_text())
        assignment = next(node for node in tree.body if isinstance(node, ast.Assign)
                          and any(isinstance(target, ast.Name) and target.id == "PRESETS"
                                  for target in node.targets))
        # Evaluate only the literal preset assignment with a detector sentinel.
        detector = object()
        actual = eval(compile(ast.Expression(assignment.value), str(source), "eval"),
                      {"dict": dict, "Detector": lambda: detector, "__builtins__": {}})
        for name, settings in actual.items():
            self.assertEqual(compiler.PRESET_OPTIONS[name],
                             dict(grid=settings["grid"], band=settings["band"],
                                  ghosts=settings["ghosts"], mtf=settings["detector"] is not None))

    def test_independent_knobs(self):
        args = options("--preset", "medium", "--monochromatic", "--no-ghosts", "--mtf")
        self.assertIsNone(args.resolved["band"])
        self.assertEqual(args.resolved["construction_preset"], "fast")
        self.assertTrue(args.resolved["mtf"])
        self.assertFalse(args.resolved["ghosts"])
        args = options("--ghosts", "--band-width", "0.03", "--band-centre-m", "1.6e-6",
                       "--anchors", "5", "--channels", "11", "--no-mtf")
        self.assertEqual(args.resolved["band"], (1.6e-6, 0.03, 5))
        self.assertEqual(args.resolved["channels"], 11)
        self.assertTrue(args.resolved["ghosts"])
        self.assertFalse(args.resolved["mtf"])
        self.assertEqual(options("--preset", "full", "--monochromatic").resolved["grid"], 1024)

    def test_configuration_modes_and_fixed_field(self):
        for name, mode, fixed in (("reference", "sky", "fitted"),
                                  ("scc-bench", "source", "fitted"), ("scc-sky", "sky", "off")):
            args = options("--configuration", name)
            self.assertEqual(args.mode, mode)
            self.assertEqual(args.resolved["fixed_ghost"], fixed)
        self.assertEqual(options("--configuration", "scc-sky", "--fixed-ghost", "fitted")
                         .resolved["fixed_ghost"], "fitted")

    def test_invalid_arguments_fail_before_import(self):
        invalid = [("--grid", "0"), ("--anchors", "-1"), ("--channels", "0"),
                   ("--band-width", "nan"), ("--band-width", "inf"), ("--band-width", "2"),
                   ("--band-width", "1e-100"),
                   ("--band-width", "1.9", "--band-centre-m", "1e308"),
                   ("--band-centre-m", "0"), ("--band-centre-m", "-1"),
                   ("--ghosts", "--no-ghosts"), ("--mtf", "--no-mtf"),
                   ("--monochromatic", "--band-width", "0.01"),
                   ("--monochromatic", "--anchors", "3"), ("--channels", "4"),
                   ("--configuration", "scc-bench", "--mode", "sky"),
                   ("--configuration", "scc-sky", "--mode", "source"),
                   ("--configuration", "llowfs")]
        for flags in invalid:
            with self.subTest(flags=flags), contextlib.redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit) as error:
                    options(*flags)
                self.assertEqual(error.exception.code, 2)


class Configuration(unittest.TestCase):
    def test_copy_model_removes_only_branch_fit_before_build(self):
        @dataclasses.dataclass
        class Model:
            screens: dict
            adjust: dict
            prescription: str

        model = Model(screens={"DM": object(), "L2B": object(), "Apodizer": object()},
                      adjust={"wedge:L2B.i2": np.array(0.1), "gamma:Apodizer": np.array(0.2)},
                      prescription="contains the actual L2B lens")
        manifest = dict(screens=list(model.screens), slots=list(model.adjust))
        for configuration in ("scc-bench", "scc-sky"):
            selected, recorded, screens, slots = compiler.configured_model(model, manifest, configuration)
            self.assertEqual(screens, ["L2B"])
            self.assertEqual(slots, ["wedge:L2B.i2"])
            self.assertEqual(list(selected.screens), ["DM", "Apodizer"])
            self.assertEqual(list(selected.adjust), ["gamma:Apodizer"])
            self.assertEqual(recorded["screens"], ["DM", "Apodizer"])
            self.assertEqual(selected.prescription, model.prescription)
            self.assertIs(selected.screens["DM"], model.screens["DM"])
        self.assertIn("L2B", model.screens)
        self.assertIn("wedge:L2B.i2", model.adjust)
        selected, recorded, screens, slots = compiler.configured_model(model, manifest, "reference")
        self.assertIs(selected, model)
        self.assertIs(recorded, manifest)
        self.assertEqual(screens, [])
        self.assertEqual(slots, [])

    def test_manifest_filter_preserves_original_and_unclassified_paths(self):
        original = synthetic_manifest()
        filtered, removed = compiler.configured_manifest(original, "scc-sky", "off")
        self.assertEqual(removed, ["stack0.724", "ghost0.100"])
        self.assertEqual([p["label"] for p in filtered["phantoms"]], ["custom"])
        np.testing.assert_array_equal(filtered["fast"]["phantom_opd_um"], [2.])
        np.testing.assert_array_equal(filtered["fast"]["phantom_defocus_mm"], [5.])
        self.assertEqual(filtered["screens"], original["screens"])
        self.assertEqual(filtered["labels"], original["labels"])
        filtered["fast"]["lowphase"][0, 0] = -50
        self.assertEqual(original["fast"]["lowphase"][0, 0], 0)
        self.assertEqual(len(original["phantoms"]), 3)
        np.testing.assert_array_equal(original["fast"]["phantom_opd_um"], [1., 2., 3.])

    def test_reference_retains_every_fit_component(self):
        manifest = synthetic_manifest()
        filtered, removed = compiler.configured_manifest(manifest, "reference", "fitted")
        self.assertEqual(removed, [])
        self.assertEqual(filtered["phantoms"], manifest["phantoms"])
        for key in manifest["fast"]:
            np.testing.assert_array_equal(filtered["fast"][key], manifest["fast"][key])

    def test_replace_tier_before_outputs_and_zero_only_gamma(self):
        args = options("--preset", "medium", "--configuration", "scc-sky", "--apodizer-gamma", "nominal")
        manifest = synthetic_manifest()
        sentinel = object()
        sim = types.SimpleNamespace(manifest=manifest, model=types.SimpleNamespace(band=(1.5e-6, 1.6e-6)),
                                    fwd=sentinel, wavelengths=np.array([1.54e-6, 1.55e-6, 1.56e-6]),
                                    channels=np.array([1.54e-6, 1.56e-6]),
                                    adjust={"gamma:Apodizer": np.array(0.2), "other": np.array(5.)})
        captured = {}

        def factory(m, fwd, wavelengths, channels, dsigma, ghosts):
            captured.update(manifest=m, fwd=fwd, wavelengths=wavelengths,
                            channels=channels, dsigma=dsigma, ghosts=ghosts)
            return types.SimpleNamespace(labels=("main", "lens", "custom"))

        compiler.configure_simulation(sim, args, factory)
        self.assertIs(captured["fwd"], sentinel)
        self.assertIs(captured["wavelengths"], sim.wavelengths)
        self.assertIs(captured["channels"], sim.channels)
        edges = (1.55e-6 * (1 - 0.015 / 2), 1.55e-6 * (1 + 0.015 / 2))
        sig = np.linspace(1 / edges[1], 1 / edges[0], 9)
        self.assertEqual(captured["dsigma"], float(sig[1] - sig[0]))
        self.assertTrue(captured["ghosts"])
        self.assertEqual(sim.labels, ("main", "lens", "custom"))
        self.assertEqual(float(sim.adjust["gamma:Apodizer"]), 0.)
        self.assertEqual(float(sim.adjust["other"]), 5.)
        self.assertEqual(len(manifest["phantoms"]), 3)


if __name__ == "__main__":
    unittest.main()
