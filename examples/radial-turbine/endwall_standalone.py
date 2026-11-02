#!/usr/bin/env python3
"""Generate radial-turbine endwalls using only the Python standard library.

Run ``python endwall_standalone.py --radius 0.04 --output-dir endwalls``.
The six CSV files contain x,r coordinates in the same units as the inputs.
``build_endwalls`` accepts the parameters of endwall.py and returns the same
six curves, as lists of (x, r) tuples instead of NumPy arrays. No installation
of turbo-design, pyturbo, or other third-party packages is needed.
"""

import argparse
import csv
import math
from pathlib import Path


CURVE_NAMES = (
    "hub_inlet", "hub", "hub_outlet",
    "shroud_inlet", "shroud", "shroud_outlet",
)


def _linspace(count):
    """Include both endpoints, matching the original uniform parameter grid."""
    step = 1.0 / (count - 1)
    return [i * step for i in range(count - 1)] + [1.0]


def _bezier_value(start, stop, t):
    """Cubic Bernstein polynomial with control values [start, start, stop, stop]."""
    # The original Bezier x control points are unused: y is evaluated at t,
    # not interpolated as a function of the Bezier's x coordinate.
    return ((1.0 - t) ** 3 * start
            + 3.0 * t * (1.0 - t) ** 2 * start
            + 3.0 * t ** 2 * (1.0 - t) * stop
            + t ** 3 * stop)


def compute_normals(x, y):
    """Return unit left normals with NumPy gradient's default differences.

    Use centered differences internally and one-sided differences at the ends.
    Coordinates must have equal lengths of at least two and nonzero tangents.
    """
    if len(x) != len(y) or len(x) < 2:
        raise ValueError("x and y must have equal lengths of at least two")
    nx, ny = [], []
    for i in range(len(x)):
        lo, hi = max(0, i - 1), min(len(x) - 1, i + 1)
        divisor = hi - lo
        dx = (x[hi] - x[lo]) / divisor
        dy = (y[hi] - y[lo]) / divisor
        length = math.hypot(dx, dy)
        if not math.isfinite(length) or length == 0:
            raise ValueError("curve must have finite, nonzero tangents")
        nx.append(-dy / length)
        ny.append(dx / length)
    return nx, ny


def offset_curve(x, y, offset_distance):
    """Offset a curve along its left normals by a scalar or per-point distances."""
    nx, ny = compute_normals(x, y)
    if isinstance(offset_distance, (int, float)):
        distances = [offset_distance] * len(x)
    else:
        distances = list(offset_distance)
    if len(distances) != len(x):
        raise ValueError("offset_distance must be a scalar or one value per point")
    if not all(math.isfinite(d) for d in distances):
        raise ValueError("offset distances must be finite")
    return [(px + d * normal_x, py + d * normal_y)
            for px, py, d, normal_x, normal_y in zip(x, y, distances, nx, ny)]


def _line(start, stop):
    dx, dy = stop[0] - start[0], stop[1] - start[1]
    return [(dx * t + start[0], dy * t + start[1]) for t in _linspace(500)]


def build_endwalls(radius: float,
                   hub_outlet_radius_scale: float = 1.1,
                   shroud_inlet_offset: float = 0.1,
                   shroud_outlet_offset_ratio: float = 1.2,
                   rhub_out: float = 0.01,
                   x_stretch_factor: float = 1,
                   inlet_ext_percent: float = 0.1,
                   outlet_ext_percent: float = 0.2):
    """Return hub inlet, hub, hub outlet, shroud inlet, shroud, shroud outlet.

    Each main curve has 1000 points; each extension has 500, including shared
    endpoints. Length inputs share one unit; scales and extension fractions
    are dimensionless (0.1 means 10%). ``rhub_out`` sets the minimum hub radius.
    All parameters must be finite. Radius, hub radius scale and x stretch must
    be positive; offsets, minimum radius and extension fractions nonnegative.

    Preserve the original coordinate frame: x is NOT translated to start at
    zero. In endwall.py the xhub translation happens after vstack has copied
    the data into hub, so that translation does not affect the returned curves.
    """
    positive = {"radius": radius,
                "hub_outlet_radius_scale": hub_outlet_radius_scale,
                "x_stretch_factor": x_stretch_factor}
    nonnegative = {"shroud_inlet_offset": shroud_inlet_offset,
                   "shroud_outlet_offset_ratio": shroud_outlet_offset_ratio,
                   "rhub_out": rhub_out,
                   "inlet_ext_percent": inlet_ext_percent,
                   "outlet_ext_percent": outlet_ext_percent}
    for name, value in positive.items():
        if not math.isfinite(value) or value <= 0:
            raise ValueError(f"{name} must be finite and positive")
    for name, value in nonnegative.items():
        if not math.isfinite(value) or value < 0:
            raise ValueError(f"{name} must be finite and nonnegative")

    parameters = _linspace(1000)
    x, r = [], []
    for t in parameters:
        local_radius = radius * _bezier_value(1.0, hub_outlet_radius_scale, t)
        angle = math.radians(180.0 + 90.0 * t)
        x.append(local_radius * math.cos(angle) * x_stretch_factor)
        r.append(local_radius * math.sin(angle))
    radial_shift = -min(r) + rhub_out
    r = [value + radial_shift for value in r]
    hub = list(zip(x, r))

    offsets = [_bezier_value(shroud_inlet_offset,
                             shroud_outlet_offset_ratio * shroud_inlet_offset, t)
               for t in parameters]
    shroud = offset_curve(x, r, offsets)

    inlet_length = radius * inlet_ext_percent
    outlet_length = radius * outlet_ext_percent
    hub_inlet = _line((hub[0][0], hub[0][1] + inlet_length), hub[0])
    shroud_inlet = _line((shroud[0][0], shroud[0][1] + inlet_length), shroud[0])
    hub_outlet = _line(hub[-1], (hub[-1][0] + outlet_length, hub[-1][1]))
    shroud_outlet = _line(shroud[-1],
                         (shroud[-1][0] + outlet_length, shroud[-1][1]))
    return hub_inlet, hub, hub_outlet, shroud_inlet, shroud, shroud_outlet


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--radius", type=float, required=True,
                        help="hub arc radius (positive; all lengths use this unit)")
    parser.add_argument("--hub-outlet-radius-scale", type=float, default=1.1)
    parser.add_argument("--shroud-inlet-offset", type=float, default=0.1)
    parser.add_argument("--shroud-outlet-offset-ratio", type=float, default=1.2)
    parser.add_argument("--rhub-out", type=float, default=0.01,
                        help="minimum hub radius")
    parser.add_argument("--x-stretch-factor", type=float, default=1.0)
    parser.add_argument("--inlet-ext-percent", type=float, default=0.1,
                        help="inlet extension as a fraction of radius")
    parser.add_argument("--outlet-ext-percent", type=float, default=0.2,
                        help="outlet extension as a fraction of radius")
    parser.add_argument("--output-dir", type=Path, default=Path("endwalls"),
                        help="directory for six x,r CSV files (existing files overwritten)")
    args = vars(parser.parse_args(argv))
    output_dir = args.pop("output_dir")
    try:
        curves = build_endwalls(**args)
    except ValueError as error:
        parser.error(str(error))
    output_dir.mkdir(parents=True, exist_ok=True)
    for name, curve in zip(CURVE_NAMES, curves):
        path = output_dir / f"{name}.csv"
        with path.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream)
            writer.writerow(("x", "r"))
            writer.writerows(curve)
        print(f"Wrote {len(curve)} points to {path}")


if __name__ == "__main__":
    main()
