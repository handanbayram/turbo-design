# Endwall contour GUI

From the repository root, run:

```sh
python examples/radial-turbine/endwall_gui.py
```

The dark-mode desktop window provides all eight inputs to
[`endwall_standalone.py`](endwall_standalone.py). Edit the values and click
**Plot contour** (or press Enter). **Reset** restores the geometry from the
`radial_turbine-1D.py` example. The plot shows the hub and shroud in the
meridional `(x, r)` plane, with dashed inlet/outlet extensions and equal axis
scales. Invalid inputs show an inline error and keep the last valid plot.

Use one consistent unit for all lengths. Scales and ratios are dimensionless;
extension fractions are relative to the hub arc radius (`0.1` means 10%). These
are passage/endwall contours, not a blade's pressure/suction surface profile.

Only Python 3.9+ and Tkinter are needed; no pip packages or project installation
are required. Tkinter is included with many Python installations. On Debian or
Ubuntu, install `python3-tk` if needed. A graphical desktop/display is required.
You can also copy `endwall_gui.py` and `endwall_standalone.py` together into any
directory and run the GUI there. The original standalone CLI and CSV export
remain available; see [its documentation](endwall_standalone.md).

Run GUI integration tests with `python -S tests/test_endwall_gui.py -v`, or
`xvfb-run -a python -S tests/test_endwall_gui.py -v` on headless Linux with Xvfb
installed. GUI tests skip when Tkinter or a display is unavailable.

# How to create a new gas for cantera

[Cantera Instructions](https://cantera.org/tutorials/thermobuild.html)
Create a new gas with thermobuild
[NASA Thermobuild](https://cearun.grc.nasa.gov/ThermoBuild/index_ds.html)
