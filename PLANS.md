# Standalone preliminary radial turbine design program

## Instructions to the implementing coding agent

Implement this plan. Deliver one self-contained Python source file named
`radial_turbine_preliminary.py` in the repository root. The file must perform a
useful preliminary design of a single-stage radial inflow turbine operating with
an ideal gas. Reuse and adapt the mathematical relationships in this repository,
with explicit equation provenance in the Python comments.

The user wants the preliminary aerodynamic design portion only. Do not build a
website, GUI, CFD interface, 3D blade generator, compressor model, performance-map
system, or structural analysis package. Preserve the existing library and examples.
Do not simply wrap or import `turbodesign`.

## 1. Reference repository and source audit

Reference checkout: `https://github.com/handanbayram/turbo-design`, version 1.4.4,
commit `4be23d16d82d7d7264821f038937ba233cac7204`.

Before implementation, inspect the actual source. If the checkout has changed,
record the new commit and relevant differences. Function names and equations are
the primary locators; do not rely exclusively on line numbers.

| Original file | Relevant symbols | Mathematical role |
| --- | --- | --- |
| `turbodesign/turbine_math.py` | `inlet_calc`, `stator_calc`, `rotor_calc` | Total/static states, stationary and rotating velocity triangles, rothalpy, relative total pressure, and loss application. |
| `turbodesign/flow_math.py` | `compute_streamline_areas`, `compute_massflow`, `compute_power`, `radii_for_area` | Flow area, continuity, blockage, energy power, Euler power, and geometric sizing relationships. Inspect geometric assumptions before reusing `radii_for_area` for a radial station. |
| `turbodesign/isentropic.py` | `IsenP`, `IsenT`, `FindMachP0P`, `mass_flow_function`, `Massflow`, `area_for_massflow`, `min_area_for_massflow` | Ideal-gas compressible flow and area/flow-capacity relationships. |
| `turbodesign/bladerow.py` | `compute_gas_constants`, `row_entropy_rise`, `entropy_change`, angle/RPM properties | State definitions, density, entropy, and angle conventions. Replace Cantera queries with constant ideal-gas properties. |
| `turbodesign/loss/fixedpressureloss.py` | `FixedPressureLoss` | Prescribed pressure-loss coefficient. The coefficient's physical meaning comes from its application in `stator_calc` and `rotor_calc`. |
| `turbodesign/turbine_spool.py` | `initialize`, `solve`, `_balance_pressure`, `_angle_match`, `step_pressures` | Existing iteration strategy and boundary-condition handling; use only what the selected meanline formulation actually needs. |
| `turbodesign/passage.py` | `get_area`, `get_cutting_line`, `streamline_curvature` | Geometry reference. A full meridional curve implementation is unnecessary for the baseline meanline program. |
| `turbodesign/radeq.py`, `turbodesign/solve_radeq.py` | `radeq`, `adjust_streamlines` | Understand the spanwise physics being omitted from the baseline. These modules are not required in a strictly meanline implementation. |
| `turbodesign/shaft_match.py` | `ShaftMatch.size` | Reference for sizing algebra such as temperature drop, loading, blade speed, and area. The class itself depends on a compressor and must not be imported. |
| `examples/radial-turbine/radial_turbine-1D.py` | Setup and `build_endwalls` | Reference operating point and conventions. Its geometry, losses, and tuned angles are case-specific. |
| `mass_flow_function.md`, `radial-equilibrium-derivation.md`, `entropy_based_efficiency.md` | Derivations | Context to review critically against the implementation and thermodynamic identities. |

`examples/radial-turbine/single_execution.py` is unfinished and is not a working
standalone solver. The existing radial example uses five streamlines; do not
describe a newly simplified meanline calculation as an exact reproduction of it.

## 2. Required model scope and assumptions

- Steady, adiabatic, single-stage radial inflow turbine: stationary nozzle/stator
  followed by one rotor with an approximately axial exit.
- One representative meanline per station; no spanwise radial-equilibrium
  integration or automatic streamline relocation.
- Calorically perfect ideal gas: constant `R` and `gamma`; derive
  `cp = gamma * R / (gamma - 1)` and `cv = R / (gamma - 1)`.
- Air defaults: `R = 287.15 J/(kg K)`, `gamma = 1.4`. Both must be configurable.
  Do not silently reuse the hot-air example's `gamma = 1.35` for every case.
- Density `rho = p / (R*T)` and sound speed `a = sqrt(gamma*R*T)`.
- Explicit prescribed nozzle and rotor pressure-loss coefficients. Zero is
  allowed for loss-free verification. Document the chosen coefficient definitions.
- No cooling, leakage, disk friction, bearing loss, heat transfer, chemistry,
  slip correlation, or empirical off-design model in the baseline. Report
  aerodynamic rotor power; document that mechanical losses are omitted.
- Angles are flow angles. Do not call them finished blade metal angles unless a
  documented deviation assumption is applied. Default to zero deviation and
  identify any resulting suggested metal angles as approximate.
- Use SI internally: Pa absolute, K, kg/s, m, m/s, J/kg, W, rad/s, and radians.
  Clearly label degree and RPM conversions at the interface.

Use station numbers consistently:

- Station 0: nozzle/stator supply, specified inlet total pressure and temperature.
- Station 1: nozzle exit / rotor inlet, predominantly radial meridional flow.
- Station 2: rotor exit, approximately axial meridional flow.

Use names such as `p0_1` for total pressure at station 1 and `p_1` for static
pressure. Explain that `0` within a property name denotes total/stagnation state.
Define `C` as absolute velocity, `W` as relative velocity, `C_m` as meridional
speed, `C_theta` as swirl, and `U = omega*r` as signed blade speed. Establish one
rotation/sign convention and use it throughout. Prefer positive design RPM with
positive rotor-inlet swirl in that convention; document the mapping to the
original example's negative RPM.

## 3. Equation-source comments are mandatory

For every governing equation or tightly related equation block, place a source
comment immediately above the implementation. Cover thermodynamics, velocity
triangles, loss application, continuity, sizing, geometry, power, and efficiency.
One traceable comment may cover a short coherent block; do not annotate only the
module header and leave the actual calculations unattributed.

Use this format, replacing the example with the actual relationship:

```python
# Equation RT-07 [ADAPTED]
# Source: turbodesign/turbine_math.py::rotor_calc
# Reference commit: 4be23d16d82d7d7264821f038937ba233cac7204
# Original relation: I = cp*T0_relative - U**2/2.
# Here: constant cp, one meanline, no coolant; U = omega*r.
# This retains the rotor-inlet to rotor-exit radius change.
I = cp * T0_relative_1 - 0.5 * U_1**2
T0_relative_2 = (I + 0.5 * U_2**2) / cp
```

Classify each relationship as:

- `REUSED`: the mathematical relationship is unchanged.
- `ADAPTED`: changed notation, averaging, gas-property treatment, or scope.
- `DERIVED`: algebraically derived from cited original relationships; show the
  short derivation or reasoning.
- `ADDED`: a new design closure or relationship not implemented in the source;
  state the assumption and derivation or a real external reference.
- `CORRECTED`: an original inconsistency was corrected; explain the defect and
  the independent physical check supporting the replacement.

Do not invent an original-file citation for a newly added design assumption.
For example, a chosen speed ratio or zero-exit-swirl closure must be identified
as a design choice even if related velocity equations exist in the repository.

Include an equation index in the file's module docstring mapping equation IDs
to source files/functions, classification, and significant changes. Keep the
index consistent with the inline comments. Preserve applicable original notices.

## 4. Inputs and a mathematically closed preliminary sizing calculation

The program must size a turbine from operating requirements and explicit design
choices, rather than requiring the user to supply every finished dimension.
At minimum, accept:

- Inlet total pressure `p0_in` and total temperature `T0_in`.
- Rotor-exit static pressure `p_exit`.
- Target mass flow `mdot`.
- Positive shaft speed `rpm`.
- `R`, `gamma`, nozzle/rotor loss coefficients, and inlet/exit blockage fractions.
- Explicit meanline design choices needed to close the model. A suitable baseline
  is blade-speed/spouting-speed ratio, inlet meridional flow coefficient,
  rotor-inlet relative flow angle, and exit-mean-radius/inlet-radius ratio.
- Zero absolute exit swirl as the baseline exit-flow closure.

Provide practical example values for the design choices, describe them as
assumptions, and allow the user to change them. Do not claim that they are
optimized or validated radial-turbine correlations. Treat target power as a
computed result in the baseline; do not quietly impose it as another independent
constraint on an already closed calculation.

The following sequence is one suitable closed formulation. Verify its algebra,
definitions, and feasibility while implementing it:

1. Calculate gas properties and the ideal total-to-static enthalpy drop from
   `p0_in`, `T0_in`, and `p_exit`. Define spouting speed
   `c0 = sqrt(2*delta_h_is_ts)`.
2. From the selected speed ratio determine `U_1`, then `r_1 = U_1/omega`.
3. From the inlet flow coefficient determine `C_m1`. With the chosen relative
   inlet angle, use `W_theta1 = C_m1*tan(beta_1)` and
   `C_theta1 = W_theta1 + U_1`. Recover absolute speed and inlet flow angle.
4. Set `T0_1 = T0_in` for the adiabatic nozzle. Obtain `T_1` from its absolute
   kinetic energy. Solve the nozzle loss relation together with the ideal-gas
   total/static pressure relation to obtain `p_1` and `p0_1`.
5. Set the exit mean radius from the chosen radius ratio; compute `U_2`.
6. Compute the rotor-inlet relative state and conserve rothalpy to determine
   exit relative total temperature. Carry ideal relative total pressure through
   the radius change, then apply the rotor loss coefficient and prescribed
   `p_exit` to determine exit static temperature and relative speed.
7. With zero absolute exit swirl, `W_theta2 = -U_2`. Recover
   `C_m2 = sqrt(W_2**2 - W_theta2**2)` and the exit flow angles. Reject an
   impossible radicand with a useful input-specific error.
8. Apply `mdot = rho*C_m*A_effective` independently at inlet and exit. At the
   radial rotor inlet, `A_geometric = 2*pi*r_1*b_1`; solve for inlet blade height
   `b_1`. At the axial rotor exit, solve for the required annular geometric area.
9. Define exit mean radius explicitly. For the root-mean-square convention,
   `r_mean2**2 = (r_tip2**2 + r_hub2**2)/2`; combining this with
   `A_geometric2 = pi*(r_tip2**2 - r_hub2**2)` determines both radii. Validate
   positive hub radius and an exit tip radius smaller than the rotor inlet
   radius. Do not also hard-fix an independent hub/tip ratio; report it as an
   output for this closure.
10. Calculate absolute/relative total states, mass-flow residuals, power,
    efficiencies, and physical consistency checks.

The sizing choices and annular reconstruction above include derived/added
relations. Attribute them honestly under section 3. An alternative formulation
is acceptable if it is equally explicit about inputs, unknowns, equations, and
closure, and actually produces dimensions from the requested operating point.

## 5. Thermodynamic and numerical requirements

- Preserve the active `rotor_calc` relationship for ideal relative pressure:
  `p0rel_2_is = p0rel_1 * (T0rel_2/T0rel_1)**(gamma/(gamma-1))` for an uncooled,
  constant-property rotor, before applying losses. Relative total temperature
  and pressure generally change when the rotor radius changes.
- Follow the active nozzle and rotor pressure-loss definitions, distinguishing
  loss coefficients from efficiency and from an arbitrary percentage pressure
  reduction. For the rotor, the reference is the ideal exit relative total
  pressure, not simply the inlet relative total pressure.
- Use one consistent definition of geometric area, effective area, and blockage:
  `A_effective = A_geometric*(1-blockage)`. Avoid applying blockage twice.
- For meridional cross-sectional areas, continuity uses meridional velocity.
  Applying a mass-flow parameter with absolute total states and meridional Mach
  without accounting for swirl is not automatically valid. Any mass-flow-function
  use must state its velocity, total-state, and area-normal conventions.
- Report absolute, relative, and meridional Mach numbers with clear labels.
  A station Mach number alone does not establish blade-throat choking. Without
  throat geometry, report throat choking as unevaluated; do not claim a resolved
  choking limit or silently clamp a Mach number to one.
- Calculate power from both the total-enthalpy drop and Euler work. Require
  their agreement for the stated adiabatic meanline assumptions.
- Report conventional total-to-static and total-to-total isentropic efficiencies
  with explicit inlet/exit reference states. Inspect `compute_power` and the
  spool efficiency helpers critically: the source contains different efficiency
  definitions. Do not silently label an entropy-based diagnostic as the standard
  isentropic efficiency. Document any corrected/replaced calculation.
- Calculate physical entropy changes using
  `delta_s = cp*log(T_b/T_a) - R*log(p_b/p_a)`. With changing relative total
  temperature, a relative-pressure ratio alone is insufficient for entropy.
- Reject nonfinite/nonphysical inputs and states: nonpositive pressure,
  temperature, mass flow, speed, or dimensions; `gamma <= 1`; invalid expansion;
  loss coefficients outside the documented admissible range; blockage outside
  `[0, 1)`; singular angle choices; and nonpositive recovered energy/area terms.
- Validate computed states, not just input ranges. Do not hide failures with
  `abs`, arbitrary clipping, or replacement of NaN values. Use tolerances only
  for numerical roundoff, and explain them.
- If iteration is needed, use bounded/bracketed methods with explicit residuals,
  iteration limits, and failure reporting. An optimizer success flag alone does
  not establish a physically converged solution.

## 6. Standalone file structure and interface

Prefer Python's standard library only: `math`, `dataclasses`, `argparse`, `json`,
and `unittest` if self-tests are embedded. This meanline formulation can be
implemented algebraically. If a third-party dependency proves necessary, justify
it and keep the dependency list explicit and minimal.

The delivered file must have no runtime imports from this repository, Cantera,
PyTurbo, plot3d, findiff, pandas, plotting libraries, or CFD software. It must not
download correlations or read caches, pickles, or other source-tree data.

Suggested internal organization:

1. Module docstring: purpose, assumptions, station/sign conventions, source
   commit, equation index, and executable usage examples.
2. Input, gas, station-state, geometry, and result dataclasses.
3. Small named thermodynamic and velocity-triangle functions with source comments.
4. Loss application and area/sizing functions with source comments.
5. `design_turbine(inputs) -> DesignResult`, independent of CLI and filesystem.
6. Human-readable reporting and JSON serialization.
7. Argument parsing and `main()` under `if __name__ == '__main__':`.
8. Embedded meaningful self-tests, invocable separately.

Support these interfaces, with equivalent option spelling acceptable:

```sh
python radial_turbine_preliminary.py --help
python radial_turbine_preliminary.py --example
python radial_turbine_preliminary.py --input design.json --output results.json
python radial_turbine_preliminary.py --self-test
```

JSON files are user input/output, not additional files required to distribute
the program. Keep an example input inside the Python source. Importing the file
must not run the design, print output, or write files. Bad inputs should produce
a clear CLI error and nonzero exit status.

Report inputs and assumptions, gas properties, station states, velocity triangles,
flow angles, rotor inlet radius/diameter and blade height, exit hub/tip radii and
height, effective/geometric areas, radius ratios, RPM, power, pressure ratios,
efficiencies, losses, entropy changes, Mach numbers, and consistency residuals.
Use ordinary JSON numbers with units documented in names or metadata. Distinguish
computed flowpath dimensions from designer-selected closure parameters. The
result does not define blade thickness, blade count, complete camber surfaces,
volute geometry, or a manufactured rotor.

## 7. Verification and acceptance criteria

Write focused tests for the mathematical identities and failure modes below.
They must be runnable without importing the original repository. Avoid tests
that merely repeat the same implementation expression as their expected value.

- **Ideal-gas consistency:** `cp-cv = R`, `cp/cv = gamma`; static/total conversions
  round-trip for independently chosen states.
- **Velocity triangles:** `C_theta = W_theta + U`, and both frames share the
  meridional component. Verify angle conversions and the chosen sign convention.
- **Rothalpy:** inlet and exit values agree for an uncooled rotor with materially
  different inlet/exit radii.
- **Loss-free entropy:** zero nozzle and rotor losses give zero entropy rise
  within numerical tolerance, including the changing-radius rotor case.
- **Positive loss:** admissible positive pressure losses generate nonnegative
  entropy. Apply this check to each component, not only the whole stage.
- **Continuity:** reconstruct target mass flow independently from each station's
  density, meridional velocity, and effective geometric area.
- **Energy:** Euler power and total-enthalpy-drop power agree; extracted power is
  positive for the bundled example. Standard efficiencies are physically bounded.
- **Geometry:** reconstructed exit annulus matches both its required area and its
  declared mean-radius convention. Check invalid radius/area combinations.
- **Scale check:** at fixed operating state, RPM, and closure parameters, changing
  target mass flow changes required areas and power consistently until a stated
  geometric feasibility limit is reached.
- **Invalid cases:** bad gas constants, non-expansion, invalid blockage/losses,
  nonpositive static temperature, and impossible exit velocity are rejected with
  intelligible errors.
- **Independent reference:** compare selected pure mathematical relationships
  against the cited original functions or hand calculations during development.
  Run a complete representative air design and record its residuals. Any original
  code used as an oracle is a development-only dependency.
- **Scope difference:** the five-streamline radial example may inform a comparison,
  but matching its results exactly is not an acceptance condition for a different
  meanline geometry and closure. Explain any such comparison honestly.
- **Isolation:** copy the single Python file to a temporary directory outside the
  repository and run `--example` and `--self-test` with a standard-library-only
  interpreter/environment. Verify import-time silence and JSON output.
- **Provenance audit:** every governing equation block has a valid source comment
  or an explicit added/derived/corrected explanation, and every equation ID occurs
  in the docstring index.

The repository's known `adjust_streamlines` missing-argument defect is outside
the selected meanline scope. Do not import that path or carry its automatic
relocation call into the standalone program. Also inspect area and efficiency
helpers before copying them; related legacy functions do not all implement the
same corrected formulas as the active solver.

## 8. Implementation sequence and final handoff

1. Audit the sources and define the station, angle, loss, geometry, and efficiency
   conventions. Write an input/unknown/closure table in the module docstring.
2. Implement constant-property thermodynamics, triangles, losses, and rothalpy
   with equation-source comments and focused invariant tests.
3. Implement the closed sizing calculation and its feasibility checks. Ensure it
   creates useful inlet/exit dimensions from operating requirements.
4. Add reporting, CLI, embedded example, JSON support, and self-tests.
5. Complete the provenance audit and isolation checks. Run only the validation
   needed for this new program; do not modify unrelated library behavior.
6. Deliver `radial_turbine_preliminary.py`. In the final response, state how to run
   it, summarize the example outputs and validation, and identify the main
   approximations and any departures from the original equations. Describe it as
   a preliminary aerodynamic design tool with prescribed losses, not as a
   validated efficiency predictor or finished turbine design.

Proceed autonomously within this scope. Keep the complete runtime implementation,
documentation, example configuration, and self-tests in the single Python file.
