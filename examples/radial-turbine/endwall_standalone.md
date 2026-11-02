# Standalone radial-turbine endwalls

`endwall_standalone.py` reproduces the geometry from `endwall.py` with **zero
third-party dependencies**. Copy that one file anywhere and run it with Python
3.9 or later. It uses `math` for cubic Bezier evaluation, the variable-radius
quarter arc, finite-difference normals, and straight extensions; no operation
imports or calls pyturbo, NumPy, or turbo-design.

From the repository root, generate the geometry used by `radial_turbine-1D.py`:

```sh
python examples/radial-turbine/endwall_standalone.py \
  --radius 0.04 \
  --hub-outlet-radius-scale 0.95 \
  --shroud-inlet-offset 0.0045 \
  --shroud-outlet-offset-ratio 3 \
  --rhub-out 0.005 \
  --x-stretch-factor 1 \
  --inlet-ext-percent 0.3 \
  --outlet-ext-percent 0.4 \
  --output-dir endwalls
```

The output directory contains six CSV files with an `x,r` header:

| File | Points |
| --- | ---: |
| `hub_inlet.csv` | 500 |
| `hub.csv` | 1000 |
| `hub_outlet.csv` | 500 |
| `shroud_inlet.csv` | 500 |
| `shroud.csv` | 1000 |
| `shroud_outlet.csv` | 500 |

Existing files with these names are overwritten. Coordinates use the same
length unit as the inputs. Extension percentages are fractions of `radius`,
so `0.3` means 30%. The flags mirror the original function's parameters:

| Parameter | Default | Meaning |
| --- | ---: | --- |
| `radius` | required | Hub arc radius |
| `hub_outlet_radius_scale` | 1.1 | Outlet-to-inlet arc radius ratio |
| `shroud_inlet_offset` | 0.1 | Initial normal distance from hub to shroud |
| `shroud_outlet_offset_ratio` | 1.2 | Final-to-initial normal offset ratio |
| `rhub_out` | 0.01 | Minimum hub radius after radial translation |
| `x_stretch_factor` | 1 | Axial coordinate multiplier |
| `inlet_ext_percent` | 0.1 | Inlet extension length divided by radius |
| `outlet_ext_percent` | 0.2 | Outlet extension length divided by radius |

All parameters must be finite. Radius, hub radius scale, and axial stretch
must be positive; the remaining parameters must be nonnegative. Zero offsets
and zero extension lengths are supported. These are geometric input checks,
not a guarantee of a physically suitable turbine passage.

## Using the function

With the script on your Python path:

```python
from endwall_standalone import build_endwalls

hub_inlet, hub, hub_outlet, shroud_inlet, shroud, shroud_outlet = build_endwalls(
    radius=0.04, shroud_inlet_offset=0.0045,
)
hub_all = hub_inlet[:-1] + hub + hub_outlet[1:]
shroud_all = shroud_inlet[:-1] + shroud + shroud_outlet[1:]
```

The parameters, defaults, return order, and sample counts match `endwall.py`.
To avoid requiring NumPy, each returned curve is a list of `(x, r)` tuples.
Use `hub[i][0]` to access x, or `[x for x, r in hub]` for all x coordinates.
Consumers already using NumPy can explicitly convert with `np.asarray(hub)`.
Each extension includes its join with the main curve; the slicing above removes
the duplicate join points when concatenating curves.

Two numerical details intentionally preserve the original result:

- Bezier values are evaluated at uniformly spaced **t** values, not uniformly
  spaced Bezier x coordinates. The Bezier x control points are unused.
- The hub x coordinate is **not** translated to start at zero. The original
  translates `xhub` after `np.vstack` has copied it into `hub`, leaving the
  returned geometry unchanged. Thus the hub starts at
  `x = -radius * x_stretch_factor`.

## Validation against the original

The test loads the unchanged `endwall.py` and runs its actual pyturbo operations.
It compares every coordinate of all six curves for 25 parameter sets: three
radius scales using defaults, the radial-turbine example above, zero offsets
and extensions, and 20 reproducible random combinations of all eight inputs.

On Python 3.12.14, NumPy 2.5.1, and pyturbo-aero 1.3.3, all **200,000 coordinates**
passed, with maximum absolute error **1.6925e-13** in the input length unit.
The per-coordinate tolerance is `1e-12 + 1e-10 * abs(reference)`, allowing for
floating-point differences between scalar math and NumPy, especially in normals.

Run the dependency-free tests from the repository root:

```sh
python -S tests/test_endwall_standalone.py -v
```

This runs five tests and skips the live reference comparison because site
packages are disabled. It also copies just the standalone script to a temporary
directory and executes its CLI with `-I -S`, then reads back all six CSV files.
Other checks cover normals, joins, extension lengths, offsets, invalid inputs,
and the original coordinate frame.

To run all six tests, including the original comparison, use an environment with
the repository's existing dependencies installed:

```sh
.venv/bin/python tests/test_endwall_standalone.py -v
# Or run through the repository's pytest configuration:
.venv/bin/python -m pytest tests/test_endwall_standalone.py -v
```

Those dependencies are needed only to run the original as a reference. They are
not needed to run or import `endwall_standalone.py`. The original script and its
existing consumers are unchanged.
