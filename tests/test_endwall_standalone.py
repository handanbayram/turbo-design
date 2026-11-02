"""Run with Python alone, or with repository dependencies for live comparisons.

    python -S tests/test_endwall_standalone.py -v
    .venv/bin/python tests/test_endwall_standalone.py -v

Only the reference comparison requires NumPy/pyturbo. It loads the unchanged
endwall.py directly, without mocking or replacing any of its operations.
"""

import csv
import importlib.util
import math
from pathlib import Path
import random
import runpy
import shutil
import subprocess
import sys
import tempfile
import unittest


EXAMPLE_DIR = Path(__file__).resolve().parents[1] / "examples" / "radial-turbine"
SCRIPT = EXAMPLE_DIR / "endwall_standalone.py"
standalone = runpy.run_path(str(SCRIPT))
build_endwalls = standalone["build_endwalls"]
EXAMPLE_PARAMETERS = dict(
    radius=0.04, hub_outlet_radius_scale=0.95,
    shroud_inlet_offset=0.0045, shroud_outlet_offset_ratio=3,
    rhub_out=0.005, x_stretch_factor=1,
    inlet_ext_percent=0.3, outlet_ext_percent=0.4,
)


class StandaloneEndwallTests(unittest.TestCase):
    def test_geometry_and_joins(self):
        hi, hub, ho, si, shroud, so = build_endwalls(**EXAMPLE_PARAMETERS)
        self.assertEqual([len(c) for c in (hi, hub, ho, si, shroud, so)],
                         [500, 1000, 500, 500, 1000, 500])
        for inlet, curve, outlet in ((hi, hub, ho), (si, shroud, so)):
            for a, b in ((inlet[-1], curve[0]), (curve[-1], outlet[0])):
                for left, right in zip(a, b):
                    self.assertAlmostEqual(left, right, places=14)
            self.assertTrue(all(p[0] == inlet[0][0] for p in inlet))
            self.assertTrue(all(p[1] == outlet[0][1] for p in outlet))
            self.assertAlmostEqual(inlet[0][1] - curve[0][1], 0.04 * 0.3)
            self.assertAlmostEqual(outlet[-1][0] - curve[-1][0], 0.04 * 0.4)
        self.assertAlmostEqual(min(r for _, r in hub), 0.005)
        self.assertAlmostEqual(hub[0][0], -0.04)  # Preserve original x frame.
        self.assertAlmostEqual(hub[-1][0], 0.0)
        self.assertAlmostEqual(math.dist(hub[0], shroud[0]), 0.0045)
        self.assertAlmostEqual(math.dist(hub[-1], shroud[-1]), 0.0135)

    def test_zero_offsets_and_extensions(self):
        hi, hub, ho, si, shroud, so = build_endwalls(
            0.04, shroud_inlet_offset=0, inlet_ext_percent=0, outlet_ext_percent=0)
        self.assertEqual(hub, shroud)
        self.assertEqual(hi, [hub[0]] * 500)
        self.assertEqual(ho, [hub[-1]] * 500)
        self.assertEqual(si, hi)
        self.assertEqual(so, ho)

    def test_normals_and_offsets(self):
        normals = standalone["compute_normals"]([0, 1, 3], [0, 1, 1])
        expected = [(-1 / math.sqrt(2), 1 / math.sqrt(2)),
                    (-1 / math.sqrt(10), 3 / math.sqrt(10)), (0, 1)]
        for actual, reference in zip(zip(*normals), expected):
            for a, b in zip(actual, reference):
                self.assertAlmostEqual(a, b)
        offset = standalone["offset_curve"]
        self.assertEqual(offset([0, 1], [0, 0], 2), [(0, 2), (1, 2)])
        self.assertEqual(offset([0, 1], [0, 0], [2, 3]), [(0, 2), (1, 3)])
        for x, y, distance in (([0], [0], 1), ([0, 0], [0, 0], 1),
                               ([0, 1], [0], 1), ([0, 1], [0, 0], [1])):
            with self.subTest(x=x, y=y, distance=distance), self.assertRaises(ValueError):
                offset(x, y, distance)

    def test_invalid_parameters(self):
        for name in EXAMPLE_PARAMETERS:
            for value in (-1, math.inf, -math.inf, math.nan):
                with self.subTest(name=name, value=value), self.assertRaises(ValueError):
                    build_endwalls(**{**EXAMPLE_PARAMETERS, name: value})
        for name in ("radius", "hub_outlet_radius_scale", "x_stretch_factor"):
            with self.subTest(name=name), self.assertRaises(ValueError):
                build_endwalls(**{**EXAMPLE_PARAMETERS, name: 0})

    def test_copied_cli_without_site_packages(self):
        # Copy only this single file outside the repo and disable all site packages.
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            script = root / SCRIPT.name
            shutil.copyfile(SCRIPT, script)
            command = [sys.executable, "-I", "-S", str(script)]
            args = [arg for name, value in EXAMPLE_PARAMETERS.items()
                    for arg in ("--" + name.replace("_", "-"), str(value))]
            result = subprocess.run(command + args + ["--output-dir", "geometry"],
                                    cwd=root, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            expected = build_endwalls(**EXAMPLE_PARAMETERS)
            self.assertEqual(len(list((root / "geometry").glob("*.csv"))), 6)
            for name, curve in zip(standalone["CURVE_NAMES"], expected):
                with (root / "geometry" / f"{name}.csv").open(newline="") as stream:
                    reader = csv.reader(stream)
                    self.assertEqual(next(reader), ["x", "r"])
                    self.assertEqual([tuple(map(float, row)) for row in reader], curve)
            invalid = subprocess.run(command + ["--radius", "nan", "--output-dir", "bad"],
                                     cwd=root, capture_output=True, text=True)
            self.assertEqual(invalid.returncode, 2)
            self.assertIn("radius must be finite and positive", invalid.stderr)
            self.assertFalse((root / "bad").exists())


class OriginalComparisonTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if importlib.util.find_spec("pyturbo") is None:
            raise unittest.SkipTest("live comparison requires the original pyturbo dependencies")
        cls.original = staticmethod(runpy.run_path(str(EXAMPLE_DIR / "endwall.py"))["build_endwalls"])

    def test_all_coordinates_against_original(self):
        cases = [dict(radius=r) for r in (0.001, 0.04, 2.0)]
        cases += [EXAMPLE_PARAMETERS,
                  dict(radius=0.04, shroud_inlet_offset=0,
                       inlet_ext_percent=0, outlet_ext_percent=0)]
        rng = random.Random(20261102)
        for _ in range(20):
            radius = 10 ** rng.uniform(-3, 0.3)
            cases.append(dict(
                radius=radius, hub_outlet_radius_scale=rng.uniform(0.5, 1.8),
                shroud_inlet_offset=radius * rng.uniform(0.01, 0.5),
                shroud_outlet_offset_ratio=rng.uniform(0, 4),
                rhub_out=radius * rng.uniform(0, 0.3),
                x_stretch_factor=rng.uniform(0.4, 2.5),
                inlet_ext_percent=rng.uniform(0, 0.6),
                outlet_ext_percent=rng.uniform(0, 0.8),
            ))
        max_error = 0.0
        coordinates = 0
        for index, parameters in enumerate(cases):
            with self.subTest(case=index, parameters=parameters):
                reference = self.original(**parameters)
                result = build_endwalls(**parameters)
                self.assertEqual(len(result), len(reference))
                for name, actual, expected in zip(standalone["CURVE_NAMES"], result, reference):
                    self.assertEqual((len(actual), 2), expected.shape)
                    for point_index, (a, b) in enumerate(zip(actual, expected)):
                        for coordinate, (value, target) in enumerate(zip(a, b)):
                            error = abs(value - target)
                            max_error = max(max_error, error)
                            coordinates += 1
                            self.assertLessEqual(
                                error, 1e-12 + 1e-10 * abs(target),
                                f"{name}[{point_index}][{coordinate}]: {value} != {target}")
        print(f"\nOriginal comparison: {len(cases)} cases, {coordinates} coordinates; "
              f"maximum absolute error {max_error:.6g}")


if __name__ == "__main__":
    unittest.main()
