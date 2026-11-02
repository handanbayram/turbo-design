#!/usr/bin/env python3
"""Simple dark-mode viewer for endwall_standalone.py (Python + Tkinter only)."""

import math
import tkinter as tk
from tkinter import font

from endwall_standalone import build_endwalls


# Start with the geometry used in radial_turbine-1D.py. Lengths share one unit;
# extension values are fractions of the hub arc radius, not whole percentages.
FIELDS = (
    ("radius", "Hub arc radius", 0.04),
    ("hub_outlet_radius_scale", "Hub outlet radius scale", 0.95),
    ("shroud_inlet_offset", "Shroud inlet offset", 0.0045),
    ("shroud_outlet_offset_ratio", "Shroud outlet offset ratio", 3.0),
    ("rhub_out", "Minimum hub radius", 0.005),
    ("x_stretch_factor", "Axial stretch factor", 1.0),
    ("inlet_ext_percent", "Inlet extension fraction", 0.3),
    ("outlet_ext_percent", "Outlet extension fraction", 0.4),
)

BG = "#11151c"
PANEL = "#1a202b"
TEXT = "#e7edf5"
MUTED = "#a5b2c5"
GRID = "#2c3544"
HUB = "#67dce5"
SHROUD = "#f4bb72"
ERROR = "#ff9eaa"


class EndwallApp:
    def __init__(self, root):
        self.root = root
        self.curves = None
        self.values = {}
        root.title("Radial turbine · Endwall contour")
        root.configure(bg=BG)
        root.geometry("1100x720")
        root.minsize(900, 680)
        root.columnconfigure(1, weight=1)
        root.rowconfigure(0, weight=1)
        self.heading_font = font.nametofont("TkDefaultFont").copy()
        self.heading_font.configure(size=18, weight="bold")
        self.plot_font = font.nametofont("TkDefaultFont").copy()
        self.plot_font.configure(size=14, weight="bold")

        controls = tk.Frame(root, bg=PANEL, padx=24, pady=24)
        controls.grid(row=0, column=0, sticky="nsew", padx=(16, 0), pady=16)
        controls.columnconfigure(1, weight=1)
        tk.Label(controls, text="Endwall contour", font=self.heading_font,
                 bg=PANEL, fg=TEXT).grid(row=0, column=0, columnspan=2, sticky="w")
        tk.Label(controls, text="Radial turbine / meridional view", bg=PANEL,
                 fg=MUTED).grid(row=1, column=0, columnspan=2, sticky="w", pady=(6, 24))

        for row, (name, label, default) in enumerate(FIELDS, start=2):
            tk.Label(controls, text=label, bg=PANEL, fg=TEXT).grid(
                row=row, column=0, sticky="w", padx=(0, 12), pady=9)
            value = tk.StringVar(root, value=str(default))
            self.values[name] = value
            tk.Entry(controls, textvariable=value, width=10, bg=BG, fg=TEXT,
                     insertbackground=TEXT, selectbackground=GRID, selectforeground=TEXT,
                     relief="flat", highlightthickness=1, highlightbackground=GRID,
                     highlightcolor=HUB).grid(row=row, column=1, sticky="ew", ipady=6)

        tk.Label(controls, text="Use one length unit throughout.\nExtension fractions: 0.1 = 10% of radius.",
                 justify="left", bg=PANEL, fg=MUTED).grid(
                     row=10, column=0, columnspan=2, sticky="w", pady=(16, 20))
        actions = tk.Frame(controls, bg=PANEL)
        actions.grid(row=11, column=0, columnspan=2, sticky="ew")
        self.plot_button = tk.Button(
            actions, text="Plot contour", command=self.plot, bg=HUB, fg=BG,
            activebackground="#9cebf0", activeforeground=BG, relief="flat",
            padx=18, pady=9, cursor="hand2")
        self.plot_button.pack(side="left", padx=(0, 10))
        self.reset_button = tk.Button(
            actions, text="Reset", command=self.reset, bg=GRID, fg=TEXT,
            activebackground="#3c485b", activeforeground=TEXT, relief="flat",
            padx=16, pady=9, cursor="hand2")
        self.reset_button.pack(side="left")
        self.status = tk.Label(controls, text="", bg=PANEL, fg=MUTED,
                               wraplength=290, justify="left", anchor="nw")
        self.status.grid(row=12, column=0, columnspan=2, sticky="ew", pady=(18, 0))

        plot_panel = tk.Frame(root, bg=BG)
        plot_panel.grid(row=0, column=1, sticky="nsew", padx=16, pady=24)
        plot_panel.rowconfigure(1, weight=1)
        plot_panel.columnconfigure(0, weight=1)
        tk.Label(plot_panel, text="Meridional contour", bg=BG, fg=TEXT,
                 font=self.plot_font).grid(row=0, sticky="w", padx=12)
        self.canvas = tk.Canvas(plot_panel, bg=BG, highlightthickness=0)
        self.canvas.grid(row=1, sticky="nsew")
        legend = tk.Frame(plot_panel, bg=BG)
        legend.grid(row=2, sticky="w", padx=12)
        for label, color in (("━  Hub", HUB), ("━  Shroud", SHROUD),
                             ("┄  Inlet / outlet extensions", MUTED)):
            tk.Label(legend, text=label, bg=BG, fg=color).pack(side="left", padx=(0, 18))
        self.canvas.bind("<Configure>", self.draw_plot)
        root.bind("<Return>", lambda event: self.plot())
        self.plot()

    def plot(self):
        """Keep the last successful plot if an entry or generated curve is invalid."""
        try:
            parameters = {}
            for name, label, _ in FIELDS:
                try:
                    parameters[name] = float(self.values[name].get())
                except ValueError:
                    raise ValueError(f"{label} must be a number.") from None
            curves = build_endwalls(**parameters)
            points = [point for curve in curves for point in curve]
            if not all(math.isfinite(value) for point in points for value in point):
                raise ValueError("These inputs produce non-finite coordinates.")
            xs, rs = zip(*points)
            bounds = min(xs), max(xs), min(rs), max(rs)
            spans = bounds[1] - bounds[0], bounds[3] - bounds[2]
            if not all(math.isfinite(span) and span > 0 for span in spans):
                raise ValueError("These inputs produce an unplottable coordinate range.")
        except (ValueError, OverflowError) as error:
            self.status.configure(text=f"{error}\nPlot unchanged.", fg=ERROR)
            return
        self.curves = curves
        self.bounds = bounds
        self.draw_plot()
        self.status.configure(text="Contour updated. Equal scale on both axes.", fg=MUTED)

    def reset(self):
        for name, _, default in FIELDS:
            self.values[name].set(str(default))
        self.plot()

    def draw_plot(self, event=None):
        if self.curves is None:
            return
        canvas = self.canvas
        width, height = canvas.winfo_width(), canvas.winfo_height()
        left, top, right, bottom = 80, 28, width - 24, height - 65
        if right <= left or bottom <= top:
            return
        xmin, xmax, rmin, rmax = self.bounds
        dx, dr = xmax - xmin, rmax - rmin
        # Normalize first to avoid overflowing the scale for very small lengths.
        span = max(dx, dr)
        nx, nr = dx / span, dr / span
        scale = 0.88 * min((right - left) / nx if nx else math.inf,
                          (bottom - top) / nr if nr else math.inf)
        # Center using differences so large offsets do not overflow a midpoint sum.
        xpad = ((right - left) - nx * scale) / 2
        rpad = ((bottom - top) - nr * scale) / 2

        def project(x, r):
            return (left + xpad + (x - xmin) / span * scale,
                    bottom - rpad - (r - rmin) / span * scale)

        canvas.delete("all")
        for i in range(5):
            x, r = xmin + dx * (i / 4), rmin + dr * (i / 4)
            px, py = project(x, r)
            canvas.create_line(px, top, px, bottom, fill=GRID)
            canvas.create_line(left, py, right, py, fill=GRID)
            canvas.create_text(px, bottom + 18, text=f"{x:.3g}", fill=MUTED)
            canvas.create_text(left - 10, py, text=f"{r:.3g}", fill=MUTED, anchor="e")
        canvas.create_rectangle(left, top, right, bottom, outline=GRID)
        canvas.create_text((left + right) / 2, height - 16,
                           text="Axial position, x (input units)", fill=MUTED)
        canvas.create_text(16, (top + bottom) / 2, angle=90,
                           text="Radius, r (input units)", fill=MUTED)
        for index, curve in enumerate(self.curves):
            coords = [coordinate for point in curve for coordinate in project(*point)]
            canvas.create_line(*coords, fill=HUB if index < 3 else SHROUD,
                               width=2, dash=() if index in (1, 4) else (5, 4),
                               tags=("contour",))


def main():
    root = tk.Tk()
    EndwallApp(root)
    root.mainloop()


if __name__ == "__main__":
    main()
