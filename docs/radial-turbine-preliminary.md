# Minimum-diameter radial turbine design

This revision keeps mass flow fixed, meets the requested **aerodynamic rotor power** exactly, and chooses RPM at or below a configurable maximum to minimize the rotor inlet diameter. It preserves the original standalone, standard-library-only model and legacy `DesignInputs` API/JSON format.

## Run

Python 3.9 or newer; no packages to install.

```sh
python radial_turbine_preliminary.py --example
python radial_turbine_preliminary.py --input minimum_diameter_input.json --output results.json
python radial_turbine_preliminary.py --self-test
```

Edit `required_power_W`, `mdot`, and `rpm_max` in the supplied JSON for your actual operating requirement. The example's 100000 W, 0.5 kg/s, and 60000 RPM are illustrative. Keep the supplied pressures and temperature only if they match your operating point. To generate another editable example:

```sh
python radial_turbine_preliminary.py --write-example design.json
```

New Python API:

```python
from radial_turbine_preliminary import MinimumDiameterInputs, design_minimum_diameter

result = design_minimum_diameter(MinimumDiameterInputs(
    p0_in=400000.0, T0_in=1000.0, p_exit=100000.0,
    mdot=0.5, required_power_W=100000.0, rpm_max=60000.0,
))
print(result.geometry.inlet_diameter)
print(result.inputs.rpm)
print(result.sizing["active_constraints"])
```

## What is minimized

With fixed mass flow, inlet flow coefficient `phi`, inlet relative angle `beta_1_deg`, and zero absolute exit swirl, Euler work gives:

```text
P = mdot * U1^2 * (1 + phi*tan(beta1))
U1 = sqrt(P / (mdot * (1 + phi*tan(beta1))))
speed_ratio = U1 / spouting_speed
```

This fixes the flow velocities and thermodynamic states independently of RPM. The solver then finds the smallest inlet radius permitted by:

- the maximum RPM;
- the minimum exit hub-to-tip radius ratio;
- the maximum exit-tip-to-inlet radius ratio.

Let `q = exit_radius_ratio`, `h = min_exit_hub_tip_ratio`, `t = max_exit_tip_inlet_radius_ratio`, and `a = exit_area_geometric/(2*pi)`. The exact radius lower bounds within these fixed closures are:

```text
r_RPM = U1 / (rpm_max*pi/30)
r_hub = sqrt(a*(1+h^2)/(q^2*(1-h^2)))
r_tip = sqrt(a/(t^2-q^2))
r_min = max(r_RPM, r_hub, r_tip)
D_min = 2*r_min
rpm = U1/r_min*30/pi
```

The implementation adds a 1e-12 relative outward radius allowance to avoid roundoff crossing an active constraint. It reports the lower bounds, controlling constraints, chosen RPM, power error, and geometric margins in `sizing`.

The illustrative geometry defaults are `min_exit_hub_tip_ratio=0.20` and `max_exit_tip_inlet_radius_ratio=0.95`. These are adjustable assumptions, not validated correlations. They keep the optimum away from zero hub radius and an exit tip touching the inlet radius. Geometry may force RPM below its maximum; raising maximum RPM further then does not reduce diameter.

The minimum is for the **rotor inlet diameter at exact target power within the specified closure family**. It does not optimize inlet angle, meridional coefficient, exit mean radius ratio, losses, mass flow, or complete turbine envelope. A target incompatible with the available enthalpy drop or fixed velocity/loss closures raises `DesignError`; changing RPM cannot repair those energy/velocity inconsistencies. Power above the target is not used to trade against geometry.

## Example and validation

For the supplied example, the result is approximately:

- 100.000 kW aerodynamic power;
- 60000 RPM, controlled by the RPM cap;
- 0.142352509 m inlet diameter (142.35 mm).

All 20 embedded tests pass: 12 original physics/validation tests plus 8 new tests covering target power, each limiting constraint, geometry reconstruction at a smaller radius, RPM independence of flow states, other gas/angle choices, infeasible inputs, and JSON/report compatibility. Standalone CLI, JSON round-trip, legacy compatibility, and silent import were also checked.

Mechanical losses, structural stress limits, throat choking, blade counts/thicknesses, and manufacturing design remain outside the original model. Use `rpm_max` from an independently established mechanical limit; the sizing model does not calculate an allowable shaft speed. The original source audit and NASA Open Source Agreement notice are preserved in the Python file.
