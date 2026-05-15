"""Step 7: Pharmacokinetics and safety simulation.

Implements one-compartment PK model for adapter molecule (tagged nanobody)
concentration dynamics. Simulates ON/OFF control kinetics for modular CAR-T.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field


@dataclass
class PKParameters:
    """Pharmacokinetic parameters for simulation."""

    half_life_min: float = 15.0
    dose_mg_kg: float = 1.0
    dosing_interval_hr: float = 4.0
    num_doses: int = 6
    volume_distribution_L: float = 5.0
    bioavailability: float = 1.0
    body_weight_kg: float = 70.0


@dataclass
class PKResult:
    """Results from PK simulation."""

    time_points: list[float] = field(default_factory=list)  # hours
    concentration: list[float] = field(default_factory=list)  # ng/mL
    effective_window: list[tuple[float, float]] = field(default_factory=list)
    cmax: float = 0.0
    cmin: float = 0.0
    trough: float = 0.0
    auc: float = 0.0
    on_time_fraction: float = 0.0
    off_time_to_threshold: float = 0.0  # minutes to drop below threshold


def simulate_one_compartment(params: PKParameters) -> PKResult:
    """Simulate one-compartment IV bolus PK model with repeated dosing.

    Args:
        params: PKParameters instance.

    Returns:
        PKResult with time course and summary statistics.
    """
    ke = math.log(2) / (params.half_life_min / 60)  # elimination rate (1/hr)
    dose_mg = params.dose_mg_kg * params.body_weight_kg * params.bioavailability
    c0_per_dose = (dose_mg * 1e6) / (params.volume_distribution_L * 1e3)  # ng/mL

    # Total simulation time
    total_time_hr = params.num_doses * params.dosing_interval_hr + 4.0
    dt = 0.01  # time step in hours

    time_points = []
    concentrations = []

    t = 0.0
    dose_times = [i * params.dosing_interval_hr for i in range(params.num_doses)]

    while t <= total_time_hr:
        conc = 0.0
        for dose_t in dose_times:
            if t >= dose_t:
                elapsed = t - dose_t
                conc += c0_per_dose * math.exp(-ke * elapsed)
        time_points.append(round(t, 4))
        concentrations.append(round(conc, 4))
        t += dt

    # Summary statistics
    cmax = max(concentrations)
    cmin = min(c for c in concentrations if c > 0) if any(c > 0 for c in concentrations) else 0
    auc = _calculate_auc(time_points, concentrations)

    # Trough (concentration just before next dose)
    trough_indices = []
    for dose_t in dose_times[1:]:
        idx = int(dose_t / dt) - 1
        if 0 <= idx < len(concentrations):
            trough_indices.append(concentrations[idx])
    trough = sum(trough_indices) / len(trough_indices) if trough_indices else 0

    return PKResult(
        time_points=time_points,
        concentration=concentrations,
        cmax=round(cmax, 2),
        cmin=round(cmin, 4),
        trough=round(trough, 2),
        auc=round(auc, 2),
    )


def simulate_on_off_dynamics(
    params: PKParameters,
    ec50_ng_ml: float = 100.0,
    safety_threshold_ng_ml: float = 1000.0,
) -> dict:
    """Simulate ON/OFF dynamics with therapeutic window.

    Args:
        params: PKParameters for the adapter molecule.
        ec50_ng_ml: Concentration for 50% CAR-T activation.
        safety_threshold_ng_ml: Upper safety limit.

    Returns:
        Dict with time courses, therapeutic window, and ON/OFF metrics.
    """
    pk_result = simulate_one_compartment(params)

    # Calculate activation fraction at each time point (Hill equation)
    hill_n = 2.0  # Hill coefficient
    activation = []
    on_time = 0
    above_safety = 0
    total_points = len(pk_result.time_points)

    for conc in pk_result.concentration:
        if conc > 0:
            act = (conc ** hill_n) / (ec50_ng_ml ** hill_n + conc ** hill_n)
        else:
            act = 0.0
        activation.append(round(act, 4))

        if act >= 0.5:
            on_time += 1
        if conc > safety_threshold_ng_ml:
            above_safety += 1

    on_fraction = on_time / total_points if total_points > 0 else 0
    safety_violation_fraction = above_safety / total_points if total_points > 0 else 0

    # Time to drop below EC50 after last dose
    ke = math.log(2) / (params.half_life_min / 60)
    last_dose_c0 = pk_result.cmax
    if last_dose_c0 > ec50_ng_ml:
        time_to_off_hr = math.log(last_dose_c0 / ec50_ng_ml) / ke
        time_to_off_min = time_to_off_hr * 60
    else:
        time_to_off_min = 0

    # Therapeutic window assessment
    therapeutic_window = []
    in_window = False
    window_start = 0.0
    for i, (t, conc) in enumerate(zip(pk_result.time_points, pk_result.concentration)):
        is_therapeutic = ec50_ng_ml <= conc <= safety_threshold_ng_ml
        if is_therapeutic and not in_window:
            window_start = t
            in_window = True
        elif not is_therapeutic and in_window:
            therapeutic_window.append((window_start, t))
            in_window = False
    if in_window:
        therapeutic_window.append((window_start, pk_result.time_points[-1]))

    return {
        "pk_result": pk_result,
        "activation": activation,
        "on_time_fraction": round(on_fraction, 3),
        "safety_violation_fraction": round(safety_violation_fraction, 3),
        "time_to_off_min": round(time_to_off_min, 1),
        "therapeutic_window": therapeutic_window,
        "ec50": ec50_ng_ml,
        "safety_threshold": safety_threshold_ng_ml,
    }


def _calculate_auc(time_points: list[float], concentrations: list[float]) -> float:
    """Calculate area under the curve using trapezoidal rule.

    Args:
        time_points: Time values in hours.
        concentrations: Concentration values in ng/mL.

    Returns:
        AUC in ng·hr/mL.
    """
    auc = 0.0
    for i in range(1, len(time_points)):
        dt = time_points[i] - time_points[i - 1]
        avg_conc = (concentrations[i] + concentrations[i - 1]) / 2
        auc += dt * avg_conc
    return auc


def generate_pk_report(result: dict) -> dict:
    """Generate a summary report from ON/OFF simulation results.

    Args:
        result: Output from simulate_on_off_dynamics.

    Returns:
        Dict with formatted summary metrics.
    """
    pk = result["pk_result"]
    return {
        "Cmax (ng/mL)": pk.cmax,
        "Trough (ng/mL)": pk.trough,
        "AUC (ng·hr/mL)": pk.auc,
        "ON time fraction": f"{result['on_time_fraction'] * 100:.1f}%",
        "Time to OFF after last dose (min)": result["time_to_off_min"],
        "Safety violation fraction": f"{result['safety_violation_fraction'] * 100:.1f}%",
        "EC50 (ng/mL)": result["ec50"],
        "Safety threshold (ng/mL)": result["safety_threshold"],
        "Therapeutic windows": len(result["therapeutic_window"]),
    }
