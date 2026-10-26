#!/usr/bin/env python3
"""Standalone minimum-diameter preliminary design of an ideal-gas radial turbine.

Run with Python >= 3.9; only the standard library is required::

    python radial_turbine_preliminary.py --example
    python radial_turbine_preliminary.py --example --output results.json
    python radial_turbine_preliminary.py --write-example design.json
    python radial_turbine_preliminary.py --input design.json --output results.json
    python radial_turbine_preliminary.py --self-test

The API is design_turbine(DesignInputs(...)) -> DesignResult. Import is silent.
JSON input is a flat object using the DesignInputs field names below. Operating
requirements are mandatory; omitted gas/loss/closure fields use documented
defaults. Unknown fields, booleans, strings and nonfinite numbers are rejected.
The alternative sizing fields required_power_W and speed_ratio accept null.

Model and conventions
---------------------
Steady, adiabatic, uncooled nozzle and single radial-inflow rotor; constant R,
gamma and cp; one meanline at each station. Station 0 is the nozzle supply
(only total states specified), 1 the radial nozzle exit/rotor inlet, and 2 the
axial rotor exit. A '0' within p0/T0 denotes stagnation, not a station index;
'rel' denotes isentropic stagnation in the local rotor frame. C is absolute
velocity, W relative velocity, C_m the shared meridional component, C_theta
absolute swirl, and U = omega*r. C_m is positive along the flow: inward at 1,
axial at 2. Positive RPM and swirl follow the rotor rotation. alpha and beta
are atan2(tangential, meridional), in radians internally. Input beta_1_deg and
reported *_deg values are explicitly degrees; RPM is converted to rad/s.
Reversing U, swirl, alpha and beta maps this convention to the original radial
example's negative RPM. Its five-streamline, hot-gas, tuned-geometry calculation
is not reproduced by this different geometry and meanline closure.

Input / unknown / closure table
-------------------------------
Inputs                        Unknowns                      Closure
p0_in, T0_in, p_exit           ideal expansion, spouting C0  ideal-gas isentrope
rpm, required_power_W          U_1, inlet radius r_1         Euler power, U_1 = omega*r_1
speed_ratio (alternative)      U_1, inlet radius r_1         U_1 = speed_ratio*C0
inlet_flow_coefficient         C_m1                          C_m1 = coefficient*U_1
beta_1_deg                    W_theta1, C_theta1            inlet relative angle
nozzle_loss                   p_1, p0_1, T_1                adiabatic nozzle loss
exit_radius_ratio             RMS r_2, U_2                  r_2 = ratio*r_1
rotor_loss, p_exit             T_2, W_2, C_m2                rothalpy, loss, C_theta2=0
mdot, inlet/exit_blockage      areas, b_1, exit hub/tip       continuity and RMS annulus
R=287.15, gamma=1.4            cp, cv, density, sound speed  calorically perfect gas

The bundled example requires 100000 W at .5 kg/s and 60000 RPM. It solves for
the minimum rotor inlet diameter with those operating requirements and the
prescribed inlet triangle held fixed. With phi=inlet_flow_coefficient,
k=1+phi*tan(beta_1), and zero absolute exit swirl, Euler power is
P=mdot*k*U_1^2=mdot*k*(omega*D_1/2)^2. For k>0 this increases strictly with
D_1, so the exact power match is the diameter lower bound for P>=required_power_W.
The recovered thermodynamics and annulus must also be feasible; otherwise the
exact power requirement is rejected with an error. No overproduction is assumed.
This is a constrained algebraic minimum, not a search over RPM, mass flow,
flow coefficient, inlet angle, exit radius ratio, losses or blockage. Without
fixed RPM and inlet-triangle choices there is no such diameter guarantee.

Specify either required_power_W or speed_ratio, never both. Omitting both
preserves the legacy speed_ratio=.65 forward calculation, where power is an
output. The other bundled choices (inlet_flow_coefficient=.30, beta_1_deg=0,
exit_radius_ratio=.50, losses=.05/.10, blockage=.03/.03) are adjustable
assumptions, not optimized or validated correlations. At fixed operating
pressures, inlet triangle and speed ratio, Euler work (and therefore eta_ts)
is fixed. Increasing loss consumes exit kinetic energy and changes eta_tt,
geometry and feasibility;
it does not independently change the already prescribed work. Loss coefficients
are in [0,1):
Yn = (p0_in-p0_1)/(p0_in-p_1);
Yr = (p0rel_2_is-p0rel_2)/(p0rel_2_is-p_2).
They are neither efficiencies nor percentages of inlet total pressure.
Effective area = geometric area*(1-blockage), applied exactly once. Continuity
uses C_m normal to the meridional cross-section; no swirl-blind mass-flow
parameter or Mach clamp is used. Station Mach numbers do not resolve throat
choking: nozzle and rotor throat choking remain unevaluated.

eta_ts = cp*(T0_in-T0_2) / {cp*T0_in*[1-(p_2/p0_in)^((gamma-1)/gamma)]};
eta_tt uses p0_2 in the denominator instead. Both include nozzle and rotor
losses. Entropy changes use complete (T,p) pairs in one frame. Nozzle entropy
uses supply total to station-1 static; rotor uses static 1 to static 2.
The relative-total entropy difference is independently checked, including the
relative temperature change. Roundoff checks allow 1e-9 in normalized balances,
1e-8 J/(kg K) entropy, and 1e-10 efficiency overshoot; no state is clipped.

No cooling, leakage, disk friction, bearings, chemistry, heat transfer, slip or
off-design correlation is included. Reported power is aerodynamic rotor power;
mechanical losses are omitted. Flow angles are not finished blade metal angles:
zero deviation is assumed only for approximate suggested metal angles. Computed
flowpath dimensions do not define blade thickness/count, camber surfaces,
volute, structural suitability or a manufactured rotor. Spanwise gradients,
meridional curvature and streamline relocation are omitted.

Source audit and notices
------------------------
Audited turbo-design 1.4.4, reference commit
4be23d16d82d7d7264821f038937ba233cac7204 (the local checkout).
Upstream main at audit: 3addd1fa14307e347324b6bc717cd3164a7ebe6c; comparison
shows only PLANS.md added, with no changes to the cited mathematical sources.
Source: https://github.com/handanbayram/turbo-design . Original project by
Paht Juangphanich and contributors; README identifies the NASA Open Source
Agreement: https://opensource.org/licenses/NASA-1.3 . Preserve that notice
when redistributing these adapted relationships.

Audited turbine_math::{inlet_calc,stator_calc,rotor_calc}, flow_math::{
compute_streamline_areas,compute_massflow,compute_power,radii_for_area},
isentropic::{IsenP,IsenT,FindMachP0P,mass_flow_function,Massflow,
area_for_massflow,min_area_for_massflow}, bladerow::{compute_gas_constants,
row_entropy_rise,entropy_change,angle/RPM properties}, FixedPressureLoss,
ShaftMatch.size, TurbineSpool::{initialize,solve,_balance_pressure,_angle_match},
step_pressures, Passage::{get_area,get_cutting_line,streamline_curvature},
radeq and adjust_streamlines. The spool's iterative pressure split and angle
matching are unnecessary for this algebraically closed sizing problem.
single_execution.py is unfinished; adjust_streamlines requires a third argument
omitted by legacy spool callers. Neither path is imported or reproduced.
mass_flow_function.md, radial-equilibrium-derivation.md and
entropy_based_efficiency.md were reviewed against the active code. In particular,
a vanishing static-temperature factor is not a Mach-1 throat criterion, and
w/(w+T_exit*ds) is not the conventional isentropic efficiency. A relative total
pressure drop across a radius change alone is not entropy generation. Those
diagnostics are not substituted for the identities used here.

Equation index (all source paths below are relative to the audited repository)
---------------------------------------------------------------------------
RT-01 ADAPTED  bladerow.py::compute_gas_constants: invert cp-cv=R, cp/cv=gamma.
RT-02 ADAPTED  turbine_math.py::inlet_calc,rotor_calc; isentropic.py::IsenP,
               IsenT,FindMachP0P: constant-cp static/total conversions.
RT-03 ADAPTED  bladerow.py::compute_gas_constants; turbine_math.py::rotor_calc:
               ideal-gas density, sound speed and three labeled Mach numbers.
RT-04 ADAPTED  bladerow.py::row_entropy_rise,entropy_change: constant-cp entropy.
RT-05 ADAPTED  turbine_math.py::rotor_calc: triangles and signed flow angles.
RT-06 DERIVED  flow_math.py::compute_power (T_is); isentropic.py::IsenP:
               ideal drop and spouting speed; stable expm1 evaluation.
RT-07 ADDED    explicit speed, inlet flow and exit radius design closures;
               rpm conversion/U=omega*r adapted from rotor_calc/ShaftMatch.size.
RT-08 DERIVED  turbine_math.py::stator_calc; fixedpressureloss.py::
               FixedPressureLoss: solve nozzle loss and stagnation simultaneously.
RT-09 REUSED   turbine_math.py::rotor_calc: rothalpy through radius change.
RT-10 REUSED   turbine_math.py::rotor_calc: ideal relative pressure radius carry.
RT-11 REUSED   turbine_math.py::rotor_calc: rotor loss against ideal exit state.
RT-12 ADDED    zero absolute exit swirl; derived C_m2 from rotor_calc triangles.
RT-13 DERIVED  flow_math.py::compute_massflow: effective/geometric area sizing.
RT-14 DERIVED  flow_math.py::compute_streamline_areas: radial cylindrical inlet.
RT-15 ADDED    RMS exit mean radius plus axial annulus area; differs from
               flow_math.py::radii_for_area's arithmetic mean convention.
RT-16 ADAPTED  flow_math.py::compute_power: enthalpy and Euler rotor power.
RT-17 CORRECTED flow_math.py::compute_power; turbine_spool.py::
               overall_entropy_efficiency: conventional stage isentropic eta_tt
               replaces entropy diagnostic; zero-loss identity is independent check.
RT-18 DERIVED  RT-02/04/08/09/13/15/16: explicit independent balance residuals.
RT-19 DERIVED  RT-15/13: strict geometric upper mass-flow bound (not choking).
RT-20 DERIVED  RT-06/16 and turbine_spool.py::overall_pressure_ratio:
               report expansion ratios, ideal drops and power in SI units.
RT-21 DERIVED  RT-14/15: diameter, heights and geometric radius ratios.
RT-22 DERIVED  RT-05/07/16: invert Euler power for the minimum inlet diameter
               at fixed mdot, RPM and inlet triangle; validate exact power match.

Development verification snapshot
---------------------------------
Bundled air design: 400000 Pa / 1000 K supply, 100000 Pa exit, .5 kg/s,
60000 RPM, required aerodynamic power 100000 W. Solved speed ratio .551574976;
eta_ts=.608469908, eta_tt=.927257062; minimum inlet diameter .142353 m,
b_1=.008361 m; exit hub/tip radii .030492/.040041 m. Representative residuals:
inlet mass flow 0 kg/s, rothalpy 0 J/kg, enthalpy minus Euler power 4.37e-11 W,
power-target error 4.37e-11 W. Last digits can vary by interpreter/platform.
Original forward-mode development-only comparisons
matched 16 nozzle/rotor state properties against the actual stator_calc and
rotor_calc functions, plus IsenP, IsenT, FindMachP0P and the zero-swirl exit
area_for_massflow. Those original modules are never imported by this file.
Copying this file outside the repository and running Python with -I -S passed
the example, self-tests, JSON round-trip, CLI errors and silent-import checks.
"""

import argparse
from dataclasses import asdict, dataclass, field, fields, replace
import json
import math
import sys
from typing import Optional
import unittest


REFERENCE_COMMIT = "4be23d16d82d7d7264821f038937ba233cac7204"


class DesignError(ValueError):
    """An input or recovered state is incompatible with the stated model."""


def finite(name, value):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise DesignError(f"{name} must be a finite JSON number, got {value!r}")
    try:
        valid = math.isfinite(value)
    except OverflowError:
        valid = False
    if not valid:
        raise DesignError(f"{name} must be finite, got {value!r}")
    return value


def positive(name, value):
    finite(name, value)
    if value <= 0:
        raise DesignError(f"{name} must be positive, got {value:g}")
    return value


def fraction(name, value):
    finite(name, value)
    if not 0 <= value < 1:
        raise DesignError(f"{name} must be in [0, 1), got {value:g}")
    return value


@dataclass(frozen=True)
class Gas:
    R: float = 287.15
    gamma: float = 1.4
    cp: float = field(init=False)
    cv: float = field(init=False)

    def __post_init__(self):
        positive("R", self.R)
        finite("gamma", self.gamma)
        if self.gamma <= 1:
            raise DesignError("gamma must be greater than 1")
        # Equation RT-01 [ADAPTED]
        # Source: turbodesign/bladerow.py::compute_gas_constants
        # Reference commit: 4be23d16d82d7d7264821f038937ba233cac7204
        # Original: R=cp-cv, gamma=cp/cv; invert for constant ideal-gas cp,cv.
        cv = positive("cv", self.R / (self.gamma - 1))
        cp = positive("cp", self.gamma * cv)
        object.__setattr__(self, "cv", cv)
        object.__setattr__(self, "cp", cp)


@dataclass(frozen=True)
class DesignInputs:
    p0_in: float
    T0_in: float
    p_exit: float
    mdot: float
    rpm: float
    R: float = 287.15
    gamma: float = 1.4
    nozzle_loss: float = 0.05
    rotor_loss: float = 0.10
    inlet_blockage: float = 0.03
    exit_blockage: float = 0.03
    speed_ratio: Optional[float] = None
    inlet_flow_coefficient: float = 0.30
    beta_1_deg: float = 0.0
    exit_radius_ratio: float = 0.50
    required_power_W: Optional[float] = None

    @classmethod
    def from_dict(cls, data):
        if not isinstance(data, dict):
            raise DesignError("input JSON must be an object")
        unknown = set(data) - {f.name for f in fields(cls)}
        if unknown:
            raise DesignError(f"unknown input field(s): {', '.join(sorted(unknown))}")
        try:
            return cls(**data)
        except TypeError as exc:
            raise DesignError(f"missing required operating input: {exc}") from exc

    def validate(self):
        for name, value in asdict(self).items():
            if name in ("speed_ratio", "required_power_W") and value is None:
                continue
            finite(name, value)
        for name in ("p0_in", "T0_in", "p_exit", "mdot", "rpm",
                     "inlet_flow_coefficient", "exit_radius_ratio"):
            positive(name, getattr(self, name))
        for name in ("speed_ratio", "required_power_W"):
            if getattr(self, name) is not None:
                positive(name, getattr(self, name))
        if self.required_power_W is not None and self.speed_ratio is not None:
            raise DesignError("specify either required_power_W or speed_ratio, not both")
        for name in ("nozzle_loss", "rotor_loss", "inlet_blockage", "exit_blockage"):
            fraction(name, getattr(self, name))
        if self.p_exit >= self.p0_in:
            raise DesignError("expansion requires p_exit < p0_in")
        if self.exit_radius_ratio >= 1:
            raise DesignError("exit_radius_ratio must be strictly between 0 and 1")
        if not -90 < self.beta_1_deg < 90:
            raise DesignError("beta_1_deg must be strictly between -90 and 90 degrees")


EXAMPLE_INPUT = DesignInputs(p0_in=400000.0, T0_in=1000.0, p_exit=100000.0,
                             mdot=0.5, rpm=60000.0, required_power_W=100000.0)


@dataclass(frozen=True)
class SupplyState:
    p0_0: float
    T0_0: float


@dataclass(frozen=True)
class StationState:
    radius: float
    p: float
    T: float
    p0: float
    T0: float
    p0rel: float
    T0rel: float
    rho: float
    sound_speed: float
    C_m: float
    C_theta: float
    W_theta: float
    C: float
    W: float
    U: float
    alpha_rad: float
    beta_rad: float
    mach_absolute: float
    mach_relative: float
    mach_meridional: float


@dataclass(frozen=True)
class Geometry:
    inlet_radius: float
    inlet_diameter: float
    inlet_blade_height: float
    inlet_area_geometric: float
    inlet_area_effective: float
    exit_mean_radius_rms: float
    exit_hub_radius: float
    exit_tip_radius: float
    exit_blade_height: float
    exit_area_geometric: float
    exit_area_effective: float
    exit_hub_tip_ratio: float
    exit_tip_inlet_radius_ratio: float
    exit_mean_inlet_radius_ratio: float
    geometric_mdot_upper_bound: float


@dataclass(frozen=True)
class DesignResult:
    inputs: DesignInputs
    gas: Gas
    station_0: SupplyState
    station_1: StationState
    station_2: StationState
    geometry: Geometry
    performance: dict
    losses: dict
    residuals: dict

    def to_dict(self):
        data = asdict(self)
        data["metadata"] = {
            "model": "preliminary aerodynamic radial-inflow turbine meanline",
            "reference_commit": REFERENCE_COMMIT,
            "units": {"pressure": "Pa absolute", "temperature": "K", "radius_height": "m",
                      "area": "m^2", "velocity": "m/s", "rho": "kg/m^3",
                      "sound_speed": "m/s", "R_cp_cv_entropy": "J/(kg K)",
                      "enthalpy_rothalpy_specific_work": "J/kg", "power": "W",
                      "mdot": "kg/s", "omega": "rad/s", "rpm": "rev/min",
                      "angles": "radians unless suffixed _deg",
                      "ratios_mach_efficiency_loss_blockage": "dimensionless"},
            "station_0": "nozzle supply; static and relative states unspecified",
            "station_1": "radial nozzle exit / rotor inlet",
            "station_2": "axial rotor exit; zero absolute swirl",
            "total_state_notation": "0 in p0/T0 means stagnation; rel is the local rotor frame",
            "sizing_mode": ("minimum inlet diameter for required aerodynamic power"
                            if self.inputs.required_power_W is not None
                            else "prescribed speed ratio; power is an output"),
            "diameter_minimum_scope": ("fixed mdot, rpm and all other design closures; exact power match"
                                       if self.inputs.required_power_W is not None else "not applied"),
            "selected_closures": [("required_power_W" if self.inputs.required_power_W is not None
                                   else "speed_ratio (default .65 if omitted)"),
                                  "inlet_flow_coefficient", "beta_1_deg",
                                  "exit_radius_ratio", "zero absolute exit swirl"],
            "computed_dimensions": "geometry; hub/tip ratio is an output",
            "assumptions": ["steady adiabatic uncooled single stage", "constant ideal-gas R and gamma",
                            "one meanline; no spanwise equilibrium or relocation",
                            "prescribed losses, no slip or off-design correlation",
                            "aerodynamic rotor power; mechanical losses omitted",
                            "no leakage, disk friction, bearings, chemistry or heat transfer",
                            "zero deviation for approximate suggested metal angles only",
                            "no blade count, thickness, camber surfaces, volute or structural design"],
            "nozzle_loss_definition": "(p0_in-p0_1)/(p0_in-p_1)",
            "rotor_loss_definition": "(p0rel_2_is-p0rel_2)/(p0rel_2_is-p_2)",
            "throat_choking": "unevaluated for both rows: throat geometry is unspecified",
            "efficiency_reference": "supply total to exit static (ts) or exit absolute total (tt)",
        }
        # RT-05 interface conversion only; the solved angles remain radians.
        for name in ("station_1", "station_2"):
            data[name]["alpha_deg"] = math.degrees(data[name]["alpha_rad"])
            data[name]["beta_deg"] = math.degrees(data[name]["beta_rad"])
        data["approximate_metal_angles_deg_zero_deviation"] = {
            "nozzle_exit_alpha": data["station_1"]["alpha_deg"],
            "rotor_inlet_beta": data["station_1"]["beta_deg"],
            "rotor_exit_beta": data["station_2"]["beta_deg"],
        }
        return data


def total_from_static(gas, T, p, speed):
    positive("static temperature", T)
    positive("static pressure", p)
    finite("speed", speed)
    if speed < 0:
        raise DesignError("speed magnitude cannot be negative")
    # Equation RT-02 [ADAPTED]
    # Source: turbodesign/turbine_math.py::inlet_calc,rotor_calc;
    # turbodesign/isentropic.py::IsenP,IsenT,FindMachP0P
    # Reference commit: 4be23d16d82d7d7264821f038937ba233cac7204
    # Original: T0=T+V^2/(2cp), p0=p*(T0/T)^(gamma/(gamma-1)).
    # Here scalar constant-cp conversions; inverse uses the same isentrope.
    T0 = positive("total temperature", T + speed**2 / (2 * gas.cp))
    p0 = positive("total pressure", p * (T0 / T)**(gas.gamma / (gas.gamma - 1)))
    return T0, p0


def static_from_total(gas, T0, p0, speed):
    positive("total temperature", T0)
    positive("total pressure", p0)
    finite("speed", speed)
    if speed < 0:
        raise DesignError("speed magnitude cannot be negative")
    # Equation RT-02 [ADAPTED]
    # Source: turbodesign/turbine_math.py::inlet_calc,rotor_calc
    # Reference commit: 4be23d16d82d7d7264821f038937ba233cac7204
    # Inverse: rearrange total kinetic enthalpy, then the pressure isentrope.
    T = positive("static temperature (reduce speed or increase T0)", T0 - speed**2 / (2 * gas.cp))
    p = positive("static pressure", p0 * (T / T0)**(gas.gamma / (gas.gamma - 1)))
    return T, p


def entropy_change(gas, T_a, p_a, T_b, p_b):
    for name, value in (("T_a", T_a), ("p_a", p_a), ("T_b", T_b), ("p_b", p_b)):
        positive(name, value)
    # Equation RT-04 [ADAPTED]
    # Source: turbodesign/bladerow.py::row_entropy_rise,entropy_change
    # Reference commit: 4be23d16d82d7d7264821f038937ba233cac7204
    # Original fallback: cp_avg*ln(T_b/T_a)-R*ln(p_b/p_a).
    # Here constant cp; retain BOTH terms for changing relative total states.
    return finite("entropy change", gas.cp * math.log(T_b / T_a) - gas.R * math.log(p_b / p_a))


def make_station(gas, radius, T, p, C_m, C_theta, U):
    for name, value in (("radius", radius), ("T", T), ("p", p), ("C_m", C_m), ("U", U)):
        positive(name, value)
    finite("C_theta", C_theta)
    # Equation RT-05 [ADAPTED]
    # Source: turbodesign/turbine_math.py::rotor_calc
    # Reference commit: 4be23d16d82d7d7264821f038937ba233cac7204
    # Original: Wt=Vt-U, W=sqrt(Vm^2+Wt^2), beta=atan2(Wt,Vm).
    # Here C replaces V, one positive along-flow C_m in both velocity frames.
    W_theta = C_theta - U
    C, W = math.hypot(C_m, C_theta), math.hypot(C_m, W_theta)
    alpha, beta = math.atan2(C_theta, C_m), math.atan2(W_theta, C_m)
    T0, p0 = total_from_static(gas, T, p, C)
    T0rel, p0rel = total_from_static(gas, T, p, W)
    # Equation RT-03 [ADAPTED]
    # Source: turbodesign/bladerow.py::compute_gas_constants;
    # turbodesign/turbine_math.py::rotor_calc
    # Reference commit: 4be23d16d82d7d7264821f038937ba233cac7204
    # Original: rho=P/(R*T), a=sqrt(gamma*R*T), M=V/a, M_rel=W/a.
    # Also expose area-normal meridional Mach; none imply throat choking.
    rho = positive("density", p / (gas.R * T))
    a = positive("sound speed", math.sqrt(gas.gamma * gas.R * T))
    return StationState(radius, p, T, p0, T0, p0rel, T0rel, rho, a,
                        C_m, C_theta, W_theta, C, W, U, alpha, beta, C/a, W/a, C_m/a)


def ideal_enthalpy_drop(gas, T0, p0, p_end):
    for name, value in (("T0", T0), ("p0", p0), ("p_end", p_end)):
        positive(name, value)
    if p_end >= p0:
        raise DesignError("isentropic expansion requires end pressure below supply total pressure")
    # Equation RT-06 [DERIVED]
    # Source: turbodesign/flow_math.py::compute_power (T_is);
    # turbodesign/isentropic.py::IsenP
    # Reference commit: 4be23d16d82d7d7264821f038937ba233cac7204
    # T_is=T0*(p_end/p0)^q; dh=cp*(T0-T_is), C0=sqrt(2*dh).
    # -expm1(q*ln(p_end/p0)) avoids subtracting near-equal temperatures.
    return positive("ideal enthalpy drop", -gas.cp * T0 * math.expm1(
        ((gas.gamma - 1) / gas.gamma) * math.log(p_end / p0)))


def nozzle_state(gas, p0_in, T0_in, C, loss):
    fraction("nozzle_loss", loss)
    T, _ = static_from_total(gas, T0_in, p0_in, C)
    # Equation RT-08 [DERIVED]
    # Source: turbodesign/turbine_math.py::stator_calc;
    # turbodesign/loss/fixedpressureloss.py::FixedPressureLoss
    # Reference commit: 4be23d16d82d7d7264821f038937ba233cac7204
    # Original: p0_1=p0_in-Yn*(p0_in-p_1), T0_1=T0_in (no cooling).
    # Let K=(T0_in/T_1)^n=p0_1/p_1. Then (K-Yn)*p_1=(1-Yn)*p0_in.
    K = (T0_in / T)**(gas.gamma / (gas.gamma - 1))
    p = positive("p_1", p0_in * (1 - loss) / (K - loss))
    return T, p


def rotor_relative_exit(gas, inlet, U_2, p_exit, loss):
    positive("U_2", U_2)
    positive("p_exit", p_exit)
    fraction("rotor_loss", loss)
    # Equation RT-09 [REUSED]
    # Source: turbodesign/turbine_math.py::rotor_calc
    # Reference commit: 4be23d16d82d7d7264821f038937ba233cac7204
    # Original: I=cp*T0R_1-U_1^2/2; T0R_2=(I+U_2^2/2)/cp, no coolant.
    I = gas.cp * inlet.T0rel - 0.5 * inlet.U**2
    T0rel_2 = positive("T0rel_2 (reduce speed_ratio or inlet swirl)",
                      (I + 0.5 * U_2**2) / gas.cp)
    # Equation RT-10 [REUSED]
    # Source: turbodesign/turbine_math.py::rotor_calc (active P0R_is relation)
    # Reference commit: 4be23d16d82d7d7264821f038937ba233cac7204
    # Original: P0R_is=P0R_1*(T0R_2/T0R_1)^(gamma/(gamma-1)).
    # A radius change alters BOTH relative stagnation T and ideal p.
    p0rel_2_is = positive("p0rel_2_is", inlet.p0rel *
                          (T0rel_2 / inlet.T0rel)**(gas.gamma / (gas.gamma - 1)))
    if p0rel_2_is <= p_exit:
        raise DesignError(f"p_exit={p_exit:g} Pa must be below ideal exit relative total "
                          f"pressure {p0rel_2_is:g} Pa; reduce speed_ratio/inlet swirl or p_exit")
    # Equation RT-11 [REUSED]
    # Source: turbodesign/turbine_math.py::rotor_calc
    # Reference commit: 4be23d16d82d7d7264821f038937ba233cac7204
    # Original: P0R_2=P0R_is-Yr*(P0R_is-P_2), T_2=T0R_2*(P_2/P0R_2)^q.
    p0rel_2 = positive("p0rel_2", p0rel_2_is - loss * (p0rel_2_is - p_exit))
    T_2 = positive("T_2", T0rel_2 * (p_exit / p0rel_2)**((gas.gamma - 1) / gas.gamma))
    W_2_sq = positive("exit relative kinetic energy term", 2 * gas.cp * (T0rel_2 - T_2))
    return T_2, W_2_sq, p0rel_2_is, p0rel_2, T0rel_2, I


def area_for_flow(mdot, rho, C_m, blockage):
    for name, value in (("mdot", mdot), ("density", rho), ("C_m", C_m)):
        positive(name, value)
    fraction("blockage", blockage)
    # Equation RT-13 [DERIVED]
    # Source: turbodesign/flow_math.py::compute_massflow
    # Reference commit: 4be23d16d82d7d7264821f038937ba233cac7204
    # Original: mdot=rho*Vm*A_geometric*(1-blockage).
    # Solve with meridional area-normal velocity, not full C in swirling flow.
    effective = positive("effective area", mdot / (rho * C_m))
    geometric = positive("geometric area", effective / (1 - blockage))
    return geometric, effective


def exit_annulus(area, r_mean, r_inlet):
    for name, value in (("exit geometric area", area), ("exit RMS radius", r_mean),
                        ("inlet radius", r_inlet)):
        positive(name, value)
    # Equation RT-15 [ADDED]
    # Source: explicit RMS design convention; area identity from
    # turbodesign/flow_math.py::compute_streamline_areas (axial branch).
    # Reference commit: 4be23d16d82d7d7264821f038937ba233cac7204
    # A=pi*(rt^2-rh^2), rm^2=(rt^2+rh^2)/2 => rh,t^2=rm^2 -,+ A/(2*pi).
    # radii_for_area uses arithmetic mean and is NOT the RMS reconstruction.
    hub_sq, tip_sq = r_mean**2 - area/(2*math.pi), r_mean**2 + area/(2*math.pi)
    if hub_sq <= 0:
        raise DesignError("exit annulus has nonpositive hub radius squared; reduce mdot/exit_blockage "
                          "or increase exit_radius_ratio/inlet radius")
    hub, tip = math.sqrt(hub_sq), math.sqrt(tip_sq)
    if tip >= r_inlet:
        raise DesignError("exit tip radius must be smaller than inlet radius; reduce mdot/exit_blockage "
                          "or exit_radius_ratio, or increase inlet radius")
    positive("exit blade height", tip - hub)
    return hub, tip


def design_turbine(inputs: DesignInputs) -> DesignResult:
    """Size one stage, minimizing inlet diameter when required_power_W is given.

    Raises DesignError for an invalid input, an infeasible closure/geometry, or
    an unresolved physical balance. Dimensions are sized for the requested mdot;
    changing mdot designs a new flowpath rather than predicting off-design flow.
    Power sizing holds RPM and the other closures fixed and replaces speed_ratio
    with an exact Euler-power match. No filesystem, imports, or CLI side effects.
    """
    if not isinstance(inputs, DesignInputs):
        raise DesignError("inputs must be a DesignInputs instance")
    try:
        return _design_turbine(inputs)
    except DesignError as exc:
        if inputs.required_power_W is not None:
            raise DesignError(f"required-power design is infeasible with the specified operating "
                              f"point and fixed design choices: {exc}") from exc
        raise
    except (OverflowError, ZeroDivisionError, ValueError) as exc:
        raise DesignError("inputs produce an unrepresentable floating-point state; "
                          "check gas constants, operating point and design choices") from exc


def _design_turbine(x):
    x.validate()
    gas = Gas(x.R, x.gamma)
    dh_is_ts = ideal_enthalpy_drop(gas, x.T0_in, x.p0_in, x.p_exit)
    # Equation RT-06 [DERIVED]
    # Source: turbodesign/flow_math.py::compute_power (T_is), ideal_enthalpy_drop above
    # Reference commit: 4be23d16d82d7d7264821f038937ba233cac7204
    # Definition of spouting speed: C0^2/2 = ideal total-to-static enthalpy drop.
    c0 = positive("spouting speed", math.sqrt(2 * dh_is_ts))
    # Equation RT-07 [ADDED]
    # Source: explicit design assumptions: U_1/C0, C_m1/U_1, beta_1, r_2/r_1.
    # U=omega*r and omega=rpm*pi/30 adapted from turbine_math.py::rotor_calc
    # and shaft_match.py::ShaftMatch.size.
    # Reference commit: 4be23d16d82d7d7264821f038937ba233cac7204
    # No original correlation is claimed for the selected closure values.
    omega = positive("omega", x.rpm * math.pi / 30)
    if x.required_power_W is None:
        speed_ratio = 0.65 if x.speed_ratio is None else x.speed_ratio
        U_1 = positive("U_1", speed_ratio * c0)
    else:
        # Equation RT-22 [DERIVED]
        # Source: RT-05/07 inlet triangle and RT-16 Euler power above;
        # turbodesign/flow_math.py::compute_power, turbine_math.py::rotor_calc.
        # Reference commit: 4be23d16d82d7d7264821f038937ba233cac7204
        # C_theta1=U1*(1+phi*tan(beta1)), C_theta2=0 => P=mdot*k*U1^2.
        # For fixed positive mdot, RPM and k, P increases strictly with D1^2.
        # Thus U1=sqrt(P_required/(mdot*k)) gives the smallest diameter that
        # can supply the required power. Retain all downstream feasibility checks.
        k = positive("inlet loading factor 1+phi*tan(beta_1)",
                     1 + x.inlet_flow_coefficient * math.tan(math.radians(x.beta_1_deg)))
        specific_work = positive("required specific work", x.required_power_W / x.mdot)
        if specific_work >= dh_is_ts:
            raise DesignError("required_power_W must be below mdot times the ideal "
                              "total-to-static enthalpy drop (finite exit kinetic energy and losses)")
        U_1 = positive("U_1 for required power", math.sqrt(specific_work / k))
        speed_ratio = positive("solved speed_ratio", U_1 / c0)
    r_1 = positive("r_1", U_1 / omega)
    C_m1 = positive("C_m1", x.inlet_flow_coefficient * U_1)
    beta_1 = math.radians(x.beta_1_deg)
    C_theta1 = positive("C_theta1 (inlet swirl must follow rotation; change beta_1_deg)",
                        U_1 + C_m1 * math.tan(beta_1))
    r_2 = positive("r_2", x.exit_radius_ratio * r_1)
    U_2 = positive("U_2", omega * r_2)
    T_1, p_1 = nozzle_state(gas, x.p0_in, x.T0_in, math.hypot(C_m1, C_theta1), x.nozzle_loss)
    s1 = make_station(gas, r_1, T_1, p_1, C_m1, C_theta1, U_1)
    T_2, W_2_sq, p0rel_is, p0rel_2, T0rel_2, I_1 = rotor_relative_exit(
        gas, s1, U_2, x.p_exit, x.rotor_loss)
    # Equation RT-12 [ADDED]
    # Source: explicit zero absolute exit swirl closure, C_theta2=0.
    # Derivation from turbine_math.py::rotor_calc triangles: W_theta2=-U_2,
    # W_2^2=C_m2^2+W_theta2^2 => C_m2^2=W_2^2-U_2^2.
    # Reference commit: 4be23d16d82d7d7264821f038937ba233cac7204
    # This is a design choice, not a repository correlation or optimized angle.
    C_m2_sq = W_2_sq - U_2**2
    if not math.isfinite(C_m2_sq) or C_m2_sq <= 0:
        raise DesignError(f"impossible exit velocity: C_m2^2={C_m2_sq:g} m^2/s^2 "
                          f"with speed_ratio={speed_ratio:g}, rotor_loss={x.rotor_loss:g}, "
                          f"exit_radius_ratio={x.exit_radius_ratio:g}; reduce speed_ratio/rotor_loss "
                          "or revise the inlet triangle/exit pressure")
    s2 = make_station(gas, r_2, T_2, x.p_exit, math.sqrt(C_m2_sq), 0.0, U_2)
    A1, Ae1 = area_for_flow(x.mdot, s1.rho, s1.C_m, x.inlet_blockage)
    A2, Ae2 = area_for_flow(x.mdot, s2.rho, s2.C_m, x.exit_blockage)
    # Equation RT-14 [DERIVED]
    # Source: turbodesign/flow_math.py::compute_streamline_areas (radial branch)
    # Reference commit: 4be23d16d82d7d7264821f038937ba233cac7204
    # Constant radius cut: dr=0, dl=b_1 => A_1=pi*(r_1+r_1)*b_1.
    # Solve b_1=A_1/(2*pi*r_1); inlet diameter is 2*r_1.
    b_1 = positive("inlet blade height", A1 / (2 * math.pi * r_1))
    r_hub, r_tip = exit_annulus(A2, r_2, r_1)
    # Equation RT-19 [DERIVED]
    # Source: RT-15 and RT-13 above (RMS annulus and continuity).
    # hub>0 => A2<2*pi*r2^2; tip<r1 => A2<2*pi*(r1^2-r2^2).
    # Multiply the smaller bound by rho*C_m*(1-blockage). This strict bound
    # is geometric feasibility only, NOT blade-throat choking or a capacity map.
    mdot_bound = positive("geometric mass-flow upper bound", 2 * math.pi *
                          min(r_2**2, r_1**2 - r_2**2) * s2.rho * s2.C_m * (1-x.exit_blockage))
    # Equation RT-21 [DERIVED]
    # Source: RT-14 radial inlet and RT-15 RMS annulus above.
    # Reference commit: 4be23d16d82d7d7264821f038937ba233cac7204
    # D1=2*r1, b2=rt-rh; report computed hub/tip and exit/inlet radius ratios.
    geometry = Geometry(r_1, 2*r_1, b_1, A1, Ae1, r_2, r_hub, r_tip,
                        r_tip-r_hub, A2, Ae2, r_hub/r_tip, r_tip/r_1, r_2/r_1, mdot_bound)
    # Equation RT-16 [ADAPTED]
    # Source: turbodesign/flow_math.py::compute_power
    # Reference commit: 4be23d16d82d7d7264821f038937ba233cac7204
    # Original: power=mdot*cp*(T0_in-T0_out), Euler=mdot*(U1*Vt1-U2*Vt2).
    # Scalar meanline; adiabatic nozzle T0_1=T0_in, mechanical losses omitted.
    work_h = positive("extracted specific work", gas.cp * (x.T0_in - s2.T0))
    work_e = positive("Euler specific work", s1.U*s1.C_theta - s2.U*s2.C_theta)
    power_h, power_e = x.mdot*work_h, x.mdot*work_e
    # Equation RT-17 [CORRECTED]
    # Source: turbodesign/flow_math.py::compute_power (T_is,T0_is,eta_static);
    # turbodesign/turbine_spool.py::overall_entropy_efficiency
    # Reference commit: 4be23d16d82d7d7264821f038937ba233cac7204
    # Original active eta_total uses w/(w+T_exit*R*ln(P0R1/P0R2)), which counts
    # radius-related pressure change as loss. Replace that diagnostic with
    # eta_tt=w/[cp*T0_in*(1-(p0_2/p0_in)^q)]; eta_ts uses p_2 instead.
    # Physical check: isentropic changing-radius rotor gives eta_tt=1,
    # while eta_ts<1 with finite exit kinetic energy; entropy uses RT-04.
    dh_is_tt = ideal_enthalpy_drop(gas, x.T0_in, x.p0_in, s2.p0)
    eta_ts, eta_tt = work_h/dh_is_ts, work_h/dh_is_tt
    for name, eta in (("eta_ts", eta_ts), ("eta_tt", eta_tt)):
        if not math.isfinite(eta) or not 0 < eta <= 1 + 1e-10:
            raise DesignError(f"{name}={eta:g} is outside physical turbine efficiency bounds")
    ds_n = entropy_change(gas, x.T0_in, x.p0_in, s1.T, s1.p)
    ds_r = entropy_change(gas, s1.T, s1.p, s2.T, s2.p)
    ds_stage = entropy_change(gas, x.T0_in, x.p0_in, s2.T, s2.p)
    ds_relative = entropy_change(gas, s1.T0rel, s1.p0rel, s2.T0rel, s2.p0rel)
    for name, ds in (("nozzle", ds_n), ("rotor", ds_r), ("stage", ds_stage)):
        if ds < -1e-8:
            raise DesignError(f"{name} entropy rise is negative: {ds:g} J/(kg K)")
    # Equation RT-18 [DERIVED]
    # Source: RT-02/04/08/09/13/15/16; independently reconstruct balances.
    # Mass uses recovered dimensions (not stored effective area); loss residuals
    # use their defining equations; entropy checks static/relative frame identity.
    # Normalize by specified positive design scales, with tolerance 1e-9 for
    # floating-point roundoff only. Never repair a failed state by clipping it.
    mass_1 = s1.rho*s1.C_m*(2*math.pi*r_1*b_1)*(1-x.inlet_blockage)
    annulus_area = math.pi*(r_tip**2-r_hub**2)
    mass_2 = s2.rho*s2.C_m*annulus_area*(1-x.exit_blockage)
    I_2 = gas.cp*s2.T0rel - 0.5*s2.U**2
    residuals = {
        "massflow_1_kg_s": mass_1-x.mdot,
        "massflow_2_kg_s": mass_2-x.mdot,
        "massflow_1_relative": (mass_1-x.mdot)/x.mdot,
        "massflow_2_relative": (mass_2-x.mdot)/x.mdot,
        "rothalpy_J_kg": I_2-I_1,
        "rothalpy_normalized": (I_2-I_1)/(gas.cp*x.T0_in),
        "energy_euler_W": power_h-power_e,
        "energy_euler_relative": (work_h-work_e)/work_e,
        "energy_euler_normalized": (work_h-work_e)/(gas.cp*x.T0_in),
        "nozzle_loss_Pa": s1.p0-(x.p0_in-x.nozzle_loss*(x.p0_in-s1.p)),
        "rotor_loss_Pa": s2.p0rel-p0rel_2,
        "nozzle_T0_K": s1.T0-x.T0_in,
        "rotor_T0rel_K": s2.T0rel-T0rel_2,
        "exit_annulus_area_m2": annulus_area-A2,
        "exit_rms_radius_squared_m2": (r_tip**2+r_hub**2)/2-r_2**2,
        "entropy_additivity_J_kg_K": ds_stage-ds_n-ds_r,
        "rotor_entropy_frames_J_kg_K": ds_relative-ds_r,
    }
    checks = {
        "inlet continuity": residuals["massflow_1_relative"],
        "exit continuity": residuals["massflow_2_relative"],
        "rothalpy": residuals["rothalpy_normalized"],
        "energy/Euler": residuals["energy_euler_normalized"],
        "nozzle loss": residuals["nozzle_loss_Pa"]/x.p0_in,
        "rotor loss": residuals["rotor_loss_Pa"]/p0rel_2,
        "nozzle total temperature": residuals["nozzle_T0_K"]/x.T0_in,
        "rotor relative total temperature": residuals["rotor_T0rel_K"]/T0rel_2,
        "exit area": residuals["exit_annulus_area_m2"]/A2,
        "exit RMS radius": residuals["exit_rms_radius_squared_m2"]/r_2**2,
        "entropy additivity": residuals["entropy_additivity_J_kg_K"]/gas.cp,
        "entropy frame identity": residuals["rotor_entropy_frames_J_kg_K"]/gas.cp,
    }
    if x.required_power_W is not None:
        # RT-22: independently check the target using the recovered enthalpy drop,
        # rather than only the Euler expression used to select the diameter.
        residuals["power_target_W"] = power_h - x.required_power_W
        residuals["power_target_relative"] = residuals["power_target_W"] / x.required_power_W
        checks["required power"] = residuals["power_target_relative"]
    for name, residual in checks.items():
        if not math.isfinite(residual) or abs(residual) > 1e-9:
            raise DesignError(f"{name} consistency check failed: normalized residual={residual:g}")
    # Equation RT-20 [DERIVED]
    # Source: RT-06/16 above; turbodesign/turbine_spool.py::overall_pressure_ratio
    # Reference commit: 4be23d16d82d7d7264821f038937ba233cac7204
    # Expansion ratios are supply total / exit static or exit absolute total.
    # Report existing ideal drops and aerodynamic work/power without mechanical factors.
    performance = {
        "omega_rad_s": omega, "spouting_speed_m_s": c0,
        "speed_ratio": speed_ratio,
        "ideal_drop_ts_J_kg": dh_is_ts, "ideal_drop_tt_J_kg": dh_is_tt,
        "specific_work_enthalpy_J_kg": work_h, "specific_work_euler_J_kg": work_e,
        "aerodynamic_power_W": power_h, "euler_power_W": power_e,
        "eta_ts": eta_ts, "eta_tt": eta_tt,
        "expansion_ratio_total_static": x.p0_in/s2.p,
        "expansion_ratio_total_total": x.p0_in/s2.p0,
        "rotor_rothalpy_1_J_kg": I_1, "rotor_rothalpy_2_J_kg": I_2,
    }
    losses = {
        "nozzle_coefficient": x.nozzle_loss, "rotor_coefficient": x.rotor_loss,
        "nozzle_total_pressure_loss_Pa": x.p0_in-s1.p0,
        "rotor_p0rel_2_is_Pa": p0rel_is,
        "rotor_relative_total_pressure_loss_Pa": p0rel_is-s2.p0rel,
        "nozzle_entropy_J_kg_K": ds_n, "rotor_entropy_J_kg_K": ds_r,
        "rotor_entropy_relative_totals_J_kg_K": ds_relative, "stage_entropy_J_kg_K": ds_stage,
    }
    result = DesignResult(x, gas, SupplyState(x.p0_in, x.T0_in), s1, s2, geometry,
                          performance, losses, residuals)
    _validate_numbers(asdict(result))
    return result


def _validate_numbers(data, path="result"):
    if isinstance(data, dict):
        for key, value in data.items():
            _validate_numbers(value, f"{path}.{key}")
    elif isinstance(data, (int, float)):
        finite(path, data)


def format_report(result):
    """Show all solved numbers, with the same unit metadata as the JSON output."""
    d = result.to_dict()
    p, g = result.performance, result.geometry
    lines = ["Preliminary radial-inflow turbine aerodynamic design",
             f"Aerodynamic rotor power: {p['aerodynamic_power_W']/1000:.4f} kW",
             f"Isentropic efficiencies: total-to-static={p['eta_ts']:.6f}, total-to-total={p['eta_tt']:.6f}",
             f"Inlet diameter: {g.inlet_diameter:.6f} m; blade height: {g.inlet_blade_height:.6f} m",
             f"Exit hub/tip radii: {g.exit_hub_radius:.6f} / {g.exit_tip_radius:.6f} m",
             "Throat choking: unevaluated (no throat geometry). Mechanical losses omitted.",
             "Prescribed-loss meanline model; flowpath sizing and approximate zero-deviation angles only."]
    if result.inputs.required_power_W is not None:
        lines.insert(1, f"Required aerodynamic power: {result.inputs.required_power_W/1000:.4f} kW; "
                     "minimum inlet diameter at fixed mass flow, RPM and inlet triangle.")
        lines.insert(3, f"Solved speed ratio U_1/C0: {p['speed_ratio']:.6f}")
    for section, values in d.items():
        lines.append(f"\n{section}:")
        for name, value in values.items():
            if isinstance(value, float):
                value = f"{value:.10g}"
            elif isinstance(value, (dict, list)):
                value = json.dumps(value, ensure_ascii=False)
            lines.append(f"  {name}: {value}")
    return "\n".join(lines)


class SelfTests(unittest.TestCase):
    """Independent identities and failure modes, runnable outside the repository."""

    def setUp(self):
        self.x = replace(EXAMPLE_INPUT, required_power_W=None, speed_ratio=.65)
        self.r = design_turbine(self.x)

    def test_required_power_and_minimum_diameter(self):
        for beta in (-10, 0, 10):
            with self.subTest(beta=beta):
                x = replace(EXAMPLE_INPUT, beta_1_deg=beta)
                r = design_turbine(x)
                self.assertAlmostEqual(r.performance["aerodynamic_power_W"], x.required_power_W)
                self.assertAlmostEqual(r.performance["euler_power_W"], x.required_power_W)
                self.assertLess(abs(r.residuals["power_target_relative"]), 1e-10)
                # Independently integrate meanline torque: P=omega*mdot*r1*C_theta1.
                torque = x.mdot*r.geometry.inlet_radius*r.station_1.C_theta
                self.assertAlmostEqual(torque*x.rpm*math.pi/30, x.required_power_W)
                # Re-solve both sides of the candidate with the forward model.
                # Every smaller D lowers power for the fixed inlet triangle.
                for factor in (.99, 1.01):
                    other = design_turbine(replace(
                        x, required_power_W=None, speed_ratio=r.performance["speed_ratio"]*factor))
                    self.assertAlmostEqual(other.geometry.inlet_diameter/r.geometry.inlet_diameter,
                                           factor)
                    self.assertEqual(other.performance["aerodynamic_power_W"] < x.required_power_W,
                                     factor < 1)

    def test_required_power_scaling_and_fixed_rpm(self):
        r = design_turbine(EXAMPLE_INPUT)
        high = design_turbine(replace(EXAMPLE_INPUT, required_power_W=120000))
        fast = design_turbine(replace(EXAMPLE_INPUT, rpm=72000))
        self.assertAlmostEqual(high.geometry.inlet_diameter/r.geometry.inlet_diameter, math.sqrt(1.2))
        self.assertAlmostEqual(fast.geometry.inlet_diameter/r.geometry.inlet_diameter, 1/1.2)
        self.assertEqual(fast.inputs.rpm, 72000)
        self.assertAlmostEqual(fast.performance["aerodynamic_power_W"], EXAMPLE_INPUT.required_power_W)
        self.assertAlmostEqual(fast.station_1.U, r.station_1.U)
        self.assertAlmostEqual(fast.station_2.T, r.station_2.T)

    def test_infeasible_required_power(self):
        ideal_power = self.x.mdot*ideal_enthalpy_drop(
            Gas(), self.x.T0_in, self.x.p0_in, self.x.p_exit)
        for power, message in ((ideal_power, "ideal"), (2*ideal_power, "ideal"),
                               (160000, "impossible exit"), (1000, "hub")):
            with self.subTest(power=power), self.assertRaisesRegex(DesignError, message):
                design_turbine(replace(EXAMPLE_INPUT, required_power_W=power))
        with self.assertRaisesRegex(DesignError, "inlet loading factor"):
            design_turbine(replace(EXAMPLE_INPUT, beta_1_deg=-89))

    def test_required_power_input_validation(self):
        for power in (0, -1, True, "100000", float("nan"), float("inf")):
            with self.subTest(power=power), self.assertRaisesRegex(DesignError, "required_power_W"):
                design_turbine(replace(EXAMPLE_INPUT, required_power_W=power))
        with self.assertRaisesRegex(DesignError, "either required_power_W or speed_ratio"):
            design_turbine(replace(EXAMPLE_INPUT, speed_ratio=.65))

    def test_required_power_json_and_legacy_default(self):
        x = DesignInputs.from_dict(json.loads(json.dumps(asdict(EXAMPLE_INPUT))))
        r = design_turbine(x)
        data = json.loads(json.dumps(r.to_dict(), allow_nan=False))
        self.assertEqual(data["inputs"]["required_power_W"], 100000)
        self.assertIsNone(data["inputs"]["speed_ratio"])
        self.assertAlmostEqual(data["performance"]["speed_ratio"],
                               r.station_1.U/r.performance["spouting_speed_m_s"])
        self.assertIn("minimum inlet diameter", data["metadata"]["sizing_mode"])
        legacy_data = asdict(self.x)
        del legacy_data["speed_ratio"], legacy_data["required_power_W"]
        legacy = design_turbine(DesignInputs.from_dict(legacy_data))
        self.assertEqual(legacy.geometry, self.r.geometry)
        self.assertEqual(legacy.performance, self.r.performance)

    def test_gas_and_independent_total_static_reference(self):
        gas = Gas()
        self.assertAlmostEqual(gas.cp-gas.cv, gas.R)
        self.assertAlmostEqual(gas.cp/gas.cv, gas.gamma)
        # Hand case M=1, T=300 K: T0/T=1.2 and p0/p=1.2^3.5.
        speed = math.sqrt(1.4*287.15*300)
        T0, p0 = total_from_static(gas, 300, 100000, speed)
        self.assertAlmostEqual(T0, 360)
        self.assertAlmostEqual(p0, 189292.9158737854)
        T, p = static_from_total(gas, T0, p0, speed)
        self.assertAlmostEqual(T, 300)
        self.assertAlmostEqual(p, 100000)
        # Independent stationary zero-speed limit, and a non-air gas.
        self.assertEqual(total_from_static(gas, 410, 230000, 0), (410, 230000))
        other = design_turbine(replace(self.x, R=208, gamma=1.66))
        self.assertAlmostEqual(other.gas.cp/other.gas.cv, 1.66)

    def test_triangles_and_sign_convention(self):
        s = make_station(Gas(), .1, 300, 100000, 100, 100, 200)
        self.assertAlmostEqual(math.degrees(s.alpha_rad), 45)
        self.assertAlmostEqual(math.degrees(s.beta_rad), -45)
        self.assertAlmostEqual(s.C_theta, s.W_theta+s.U)
        self.assertAlmostEqual(s.W*math.cos(s.beta_rad), s.C_m)
        self.assertAlmostEqual(s.C*math.cos(s.alpha_rad), s.C_m)
        varied = design_turbine(replace(self.x, beta_1_deg=-10))
        self.assertAlmostEqual(math.degrees(varied.station_1.beta_rad), -10)
        self.assertGreater(varied.station_1.C_theta, 0)
        self.assertEqual(varied.station_2.alpha_rad, 0)
        self.assertLess(varied.station_2.beta_rad, 0)

    def test_rothalpy_and_changing_radius(self):
        a, b, gas = self.r.station_1, self.r.station_2, self.r.gas
        self.assertAlmostEqual(a.radius/b.radius, 2)
        # Alternative static-state expression I=h+W^2/2-U^2/2.
        I_a = gas.cp*a.T + (a.W*a.W-a.U*a.U)/2
        I_b = gas.cp*b.T + (b.W*b.W-b.U*b.U)/2
        self.assertAlmostEqual(I_a, I_b, places=7)
        self.assertLess(b.T0rel, .95*a.T0rel)

    def test_loss_free_entropy_and_standard_efficiency(self):
        r = design_turbine(replace(self.x, nozzle_loss=0, rotor_loss=0))
        for key in ("nozzle_entropy_J_kg_K", "rotor_entropy_J_kg_K", "stage_entropy_J_kg_K"):
            self.assertAlmostEqual(r.losses[key], 0, places=9)
        self.assertAlmostEqual(r.performance["eta_tt"], 1, places=10)
        self.assertLess(r.performance["eta_ts"], 1)
        # Absolute stagnation isentropic invariant, independent of relative-pressure carry.
        a, b = r.station_1, r.station_2
        self.assertAlmostEqual(b.p0/a.p0, (b.T0/a.T0)**3.5, places=10)
        self.assertLess(b.p0rel, .95*a.p0rel)

    def test_component_losses_and_entropy_are_independent(self):
        for yn, yr in ((0, .1), (.05, 0), (.05, .1)):
            with self.subTest(nozzle_loss=yn, rotor_loss=yr):
                r = design_turbine(replace(self.x, nozzle_loss=yn, rotor_loss=yr))
                self.assertGreaterEqual(r.losses["nozzle_entropy_J_kg_K"], -1e-9)
                self.assertGreaterEqual(r.losses["rotor_entropy_J_kg_K"], -1e-9)
                if yn:
                    self.assertGreater(r.losses["nozzle_entropy_J_kg_K"], 0)
                if yr:
                    self.assertGreater(r.losses["rotor_entropy_J_kg_K"], 0)
                self.assertAlmostEqual(r.losses["rotor_entropy_J_kg_K"],
                                       r.losses["rotor_entropy_relative_totals_J_kg_K"], places=9)
                # Independently recover the prescribed coefficients from pressure states.
                a, b = r.station_1, r.station_2
                self.assertAlmostEqual((self.x.p0_in-a.p0)/(self.x.p0_in-a.p), yn)
                pis = r.losses["rotor_p0rel_2_is_Pa"]
                self.assertAlmostEqual((pis-b.p0rel)/(pis-b.p), yr)

    def test_continuity_geometry_and_energy(self):
        r, gas = self.r, self.r.gas
        g, a, b = r.geometry, r.station_1, r.station_2
        A1 = 2*math.pi*g.inlet_radius*g.inlet_blade_height
        A2 = math.pi*(g.exit_tip_radius**2-g.exit_hub_radius**2)
        for s, area, blockage in ((a, A1, self.x.inlet_blockage), (b, A2, self.x.exit_blockage)):
            self.assertAlmostEqual(s.p/(gas.R*s.T)*s.C_m*area*(1-blockage), self.x.mdot)
        self.assertAlmostEqual(A2, g.exit_area_geometric)
        self.assertAlmostEqual(math.hypot(g.exit_hub_radius, g.exit_tip_radius)/math.sqrt(2),
                               g.exit_mean_radius_rms)
        self.assertGreater(g.exit_hub_radius, 0)
        self.assertLess(g.exit_tip_radius, g.inlet_radius)
        self.assertGreater(r.performance["aerodynamic_power_W"], 0)
        self.assertAlmostEqual(r.performance["aerodynamic_power_W"],
                               r.performance["euler_power_W"], places=7)
        self.assertTrue(0 < r.performance["eta_ts"] <= r.performance["eta_tt"] <= 1)
        for key in ("massflow_1_relative", "massflow_2_relative", "energy_euler_relative"):
            self.assertLess(abs(r.residuals[key]), 1e-10)

    def test_massflow_scaling_and_geometric_bound(self):
        r2 = design_turbine(replace(self.x, mdot=2*self.x.mdot))
        self.assertEqual(self.r.station_1, r2.station_1)
        self.assertEqual(self.r.station_2, r2.station_2)
        for name in ("inlet_area_geometric", "exit_area_geometric", "inlet_blade_height"):
            self.assertAlmostEqual(getattr(r2.geometry, name), 2*getattr(self.r.geometry, name))
        self.assertAlmostEqual(r2.performance["aerodynamic_power_W"],
                               2*self.r.performance["aerodynamic_power_W"])
        bound = self.r.geometry.geometric_mdot_upper_bound
        design_turbine(replace(self.x, mdot=.99*bound))
        with self.assertRaisesRegex(DesignError, "hub|tip"):
            design_turbine(replace(self.x, mdot=1.01*bound))

    def test_invalid_annuli(self):
        for args in ((.02, .04, .1), (.002, .09, .091), (0, .04, .1)):
            with self.subTest(args=args), self.assertRaisesRegex(DesignError, "hub|tip|area"):
                exit_annulus(*args)

    def test_invalid_inputs(self):
        cases = [("R", 0), ("R", -1), ("gamma", 1), ("gamma", .9),
                 ("p0_in", 0), ("T0_in", 0), ("p_exit", 0), ("p_exit", 400000),
                 ("mdot", 0), ("rpm", -1), ("speed_ratio", 0),
                 ("inlet_flow_coefficient", 0), ("exit_radius_ratio", 1),
                 ("beta_1_deg", 90), ("beta_1_deg", -90), ("beta_1_deg", -89),
                 ("nozzle_loss", -0.1), ("nozzle_loss", 1), ("rotor_loss", 1),
                 ("inlet_blockage", 1), ("exit_blockage", -0.1),
                 ("T0_in", float("nan")), ("rpm", float("inf")),
                 ("mdot", True), ("R", "287.15"), ("R", None)]
        for name, value in cases:
            with self.subTest(name=name, value=value), self.assertRaises(DesignError):
                design_turbine(replace(self.x, **{name: value}))
        for R, gamma in ((0, 1.4), (287.15, 1), (287.15, float("nan"))):
            with self.subTest(R=R, gamma=gamma), self.assertRaises(DesignError):
                Gas(R, gamma)

    def test_invalid_recovered_states(self):
        with self.assertRaisesRegex(DesignError, "static temperature"):
            design_turbine(replace(self.x, inlet_flow_coefficient=10))
        with self.assertRaisesRegex(DesignError, "impossible exit velocity"):
            design_turbine(replace(self.x, rotor_loss=.99))
        with self.assertRaisesRegex(DesignError, "relative total pressure|T0rel_2|impossible exit"):
            design_turbine(replace(self.x, speed_ratio=.9))
        with self.assertRaisesRegex(DesignError, "unrepresentable|positive|finite"):
            design_turbine(replace(self.x, rpm=1e-320))

    def test_swirl_and_mach_labels(self):
        s = self.r.station_1
        self.assertGreater(s.mach_absolute, s.mach_meridional)
        self.assertGreater(s.mach_absolute, s.mach_relative)
        self.assertAlmostEqual(s.mach_absolute*math.cos(s.alpha_rad), s.mach_meridional)
        self.assertIn("unevaluated", self.r.to_dict()["metadata"]["throat_choking"])
        supersonic = design_turbine(replace(self.x, inlet_flow_coefficient=.6))
        self.assertGreater(supersonic.station_1.mach_absolute, 1)
        self.assertIn("unevaluated", supersonic.to_dict()["metadata"]["throat_choking"])

    def test_json_contract(self):
        x = DesignInputs.from_dict(json.loads(json.dumps(asdict(self.x))))
        data = json.loads(json.dumps(design_turbine(x).to_dict(), allow_nan=False))
        self.assertEqual(data["inputs"]["mdot"], .5)
        self.assertIn("station_0", data)
        self.assertIn("units", data["metadata"])
        self.assertAlmostEqual(data["station_1"]["beta_deg"], 0)
        for bad in ({"p0_typo": 1}, {}, [], {"power": 1}):
            with self.subTest(data=bad), self.assertRaises(DesignError):
                DesignInputs.from_dict(bad)


def _unique_object(pairs):
    result = {}
    for name, value in pairs:
        if name in result:
            raise DesignError(f"duplicate JSON input field: {name}")
        result[name] = value
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--example", action="store_true", help="run the bundled air design")
    mode.add_argument("--input", metavar="PATH", help="read DesignInputs from a JSON file")
    mode.add_argument("--self-test", action="store_true", help="run embedded standard-library tests")
    mode.add_argument("--write-example", metavar="PATH", help="write the editable bundled JSON input")
    parser.add_argument("--output", metavar="PATH", help="write complete results as JSON")
    args = parser.parse_args(argv)
    if args.output and (args.self_test or args.write_example):
        parser.error("--output requires --example or --input")
    if args.self_test:
        suite = unittest.defaultTestLoader.loadTestsFromTestCase(SelfTests)
        return 0 if unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful() else 1
    try:
        if args.write_example:
            with open(args.write_example, "w", encoding="utf-8") as stream:
                json.dump(asdict(EXAMPLE_INPUT), stream, indent=2, allow_nan=False)
                stream.write("\n")
            print(f"Wrote example input to {args.write_example}")
            return 0
        inputs = EXAMPLE_INPUT
        if args.input:
            with open(args.input, encoding="utf-8") as stream:
                inputs = DesignInputs.from_dict(json.load(stream, object_pairs_hook=_unique_object))
        result = design_turbine(inputs)
        if args.output:
            with open(args.output, "w", encoding="utf-8") as stream:
                json.dump(result.to_dict(), stream, indent=2, allow_nan=False)
                stream.write("\n")
        print(format_report(result))
        return 0
    except (DesignError, OSError, ValueError, OverflowError) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    sys.exit(main())
