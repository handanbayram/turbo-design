"""GUI integration tests: python -S tests/test_endwall_gui.py -v.

Requires Tkinter and a display (for example, xvfb-run on headless Linux).
The geometry-only tests remain runnable without either in test_endwall_standalone.py.
"""

from pathlib import Path
import sys
import unittest


EXAMPLE_DIR = Path(__file__).resolve().parents[1] / "examples" / "radial-turbine"


class EndwallGuiTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        try:
            import tkinter as tk
        except ImportError as error:
            raise unittest.SkipTest(f"Tkinter unavailable: {error}")
        try:
            cls.root = tk.Tk()
        except tk.TclError as error:
            raise unittest.SkipTest(f"Display unavailable: {error}")
        cls.addClassCleanup(cls.root.destroy)
        sys.path.insert(0, str(EXAMPLE_DIR))
        try:
            import endwall_gui
        finally:
            sys.path.pop(0)
        cls.gui = endwall_gui
        cls.app = endwall_gui.EndwallApp(cls.root)

    def setUp(self):
        self.root.geometry("1100x720")
        self.app.reset_button.invoke()
        self.root.update()

    def assert_plot_matches_geometry(self):
        lines = self.app.canvas.find_withtag("contour")
        self.assertEqual(len(lines), 6)
        # Recover the affine transform from the hub's endpoints, then verify every
        # rendered vertex against the generator, including the four extensions.
        hub = self.app.curves[1]
        pixels = self.app.canvas.coords(lines[1])
        sx = (pixels[-2] - pixels[0]) / (hub[-1][0] - hub[0][0])
        sr = (pixels[-1] - pixels[1]) / (hub[-1][1] - hub[0][1])
        self.assertGreater(sx, 0)
        self.assertAlmostEqual(sx / -sr, 1.0)
        for line, curve in zip(lines, self.app.curves):
            actual = self.app.canvas.coords(line)
            self.assertEqual(len(actual), 2 * len(curve))
            for index, (x, r) in enumerate(curve):
                self.assertAlmostEqual(actual[2 * index], pixels[0] + (x - hub[0][0]) * sx)
                self.assertAlmostEqual(actual[2 * index + 1], pixels[1] + (r - hub[0][1]) * sr)
                self.assertGreaterEqual(actual[2 * index], 80)
                self.assertLessEqual(actual[2 * index], self.app.canvas.winfo_width() - 24)
                self.assertGreaterEqual(actual[2 * index + 1], 28)
                self.assertLessEqual(actual[2 * index + 1], self.app.canvas.winfo_height() - 65)

    def test_initial_plot_uses_example_geometry(self):
        expected = self.gui.build_endwalls(**{name: default for name, _, default in self.gui.FIELDS})
        self.assertEqual(self.app.curves, expected)
        self.assert_plot_matches_geometry()

    def test_edit_plot_and_reset(self):
        self.app.values["radius"].set("0.06")
        self.app.values["x_stretch_factor"].set("1.5")
        self.app.plot_button.invoke()
        self.assertAlmostEqual(self.app.curves[1][0][0], -0.09)
        self.assert_plot_matches_geometry()
        self.app.reset_button.invoke()
        for name, _, default in self.gui.FIELDS:
            self.assertEqual(float(self.app.values[name].get()), default)
        self.assertAlmostEqual(self.app.curves[1][0][0], -0.04)
        self.assert_plot_matches_geometry()

    def test_invalid_inputs_keep_last_plot_and_recover(self):
        curves = self.app.curves
        lines = self.app.canvas.find_withtag("contour")
        for name, value in (("radius", ""), ("radius", "abc"), ("radius", "nan"),
                            ("radius", "inf"), ("radius", "0"), ("radius", "-1"),
                            ("shroud_inlet_offset", "-0.1"),
                            ("x_stretch_factor", "1e309")):
            with self.subTest(name=name, value=value):
                original = self.app.values[name].get()
                self.app.values[name].set(value)
                self.app.plot_button.invoke()
                self.assertTrue(self.app.curves is curves, "Invalid inputs replaced the last plot")
                self.assertEqual(self.app.canvas.find_withtag("contour"), lines)
                self.assertIn("Plot unchanged", self.app.status.cget("text"))
                self.app.values[name].set(original)
        self.app.plot_button.invoke()
        self.assertNotIn("Plot unchanged", self.app.status.cget("text"))
        self.assert_plot_matches_geometry()

    def test_geometry_overflow_keeps_last_plot(self):
        curves = self.app.curves
        self.app.values["radius"].set("1e308")
        self.app.values["x_stretch_factor"].set("100")
        self.app.plot_button.invoke()
        self.assertTrue(self.app.curves is curves)
        self.assertIn("Plot unchanged", self.app.status.cget("text"))

    def test_resize_preserves_equal_scale_and_all_curves(self):
        for size in ("900x680", "1200x800"):
            with self.subTest(size=size):
                self.root.geometry(size)
                self.root.update()
                self.assert_plot_matches_geometry()

    def test_zero_extensions_and_offset(self):
        for name in ("inlet_ext_percent", "outlet_ext_percent", "shroud_inlet_offset"):
            self.app.values[name].set("0")
        self.app.plot_button.invoke()
        self.assertEqual(self.app.curves[1], self.app.curves[4])
        self.assert_plot_matches_geometry()


if __name__ == "__main__":
    unittest.main()
