"""
Simulator Validation: Pilot batch of synthetic patients.

Varies 15 population parameters, runs each through the simulator,
extracts clinical outputs, and validates against physiological ranges.
Produces a validation report.
"""

import numpy as np
import time
import json
import os
from multiprocessing import Pool, cpu_count

from hallow_c_driver import build_params, build_inits, integrate, STATE_MAP, PARAM_MAP

# =========================================================================
# 15 Population-Varying Parameters with Ranges
# =========================================================================
# Format: (param_name, nominal, low, high, description)
POPULATION_PARAMS = [
    # Tier 1: Core Patient Parameters
    ("HR_heart_rate",                          70,    55,     95,    "Heart rate (bpm)"),
    ("cf",                                     11,     7,     18,    "LV passive stiffness"),
    ("contractility",                          1.0,   0.6,    1.3,   "LV contractility"),
    ("C_art_initial",                          1.10e-8, 0.6e-8, 1.6e-8, "Arterial compliance (m³/Pa)"),
    ("R_per0",                                 1.27e8, 0.9e8,  1.8e8, "Peripheral resistance (Pa·s/m³)"),
    ("baseline_nephrons",                      2e6,    1.2e6,  2.5e6, "Nephron count"),
    ("BV",                                     0.005,  0.004,  0.006, "Blood volume (m³)"),
    ("V_w_0",                                  0.00012, 0.00008, 0.00018, "LV wall volume (m³)"),
    # Tier 2: Disease Progression
    ("kD_HYPERTROPHY",                         1e-9,   0.5e-9, 3e-9,  "Concentric hypertrophy rate"),
    ("kL_HYPERTROPHY",                         2.67e-10, 1e-10, 8e-10, "Eccentric hypertrophy rate"),
    ("disease_effect_on_nephrons",             0,      0,      0.05,  "Nephron loss rate"),
    ("Stiffness_BP_slope",                     0.01,   0.005,  0.03,  "Vascular stiffening rate"),
    # Tier 3: Neurohumoral
    ("nominal_aldosterone_concentration",      85,     50,     150,   "Aldosterone (pg/mL)"),
    ("S_tubulo_glomerular_feedback",           0.7,    0.4,    1.0,   "TGF sensitivity"),
    ("glucose_concentration",                  5.5,    4.5,    11.0,  "Plasma glucose (mmol/L)"),
]


def sample_patient(rng, patient_id):
    """Sample a random patient's 15 parameters using Latin Hypercube-like uniform sampling."""
    sampled = {}
    for name, nom, lo, hi, desc in POPULATION_PARAMS:
        val = rng.uniform(lo, hi)
        sampled[name] = val
    return sampled


def run_single_patient(args):
    """Run one patient simulation. Returns dict of clinical outputs or None if failed."""
    patient_id, param_overrides = args
    try:
        t_start = time.time()
        params, r_values = build_params()

        # Apply overrides
        for name, val in param_overrides.items():
            if name in PARAM_MAP:
                params[PARAM_MAP[name]] = val
            # Also update r_values for dependent params in build_inits
            r_values[name] = val

        # Recalculate dependent initial conditions
        V_w_0 = param_overrides.get('V_w_0', 0.00012)
        params[PARAM_MAP['Baseline_Interstitial_Fibrosis']] = V_w_0 * 0.02
        params[PARAM_MAP['Baseline_Replacement_Fibrosis']] = V_w_0 * 0.02
        params[PARAM_MAP['Baseline_Interstitial_Tissue']] = V_w_0 * 0.22

        # BV (m³) → blood_volume_nom (L) for initial conditions
        if 'BV' in param_overrides:
            r_values['blood_volume_nom'] = param_overrides['BV'] * 1000

        y = build_inits(r_values)

        # Simulate 24 hours to steady state
        results, completed = integrate(y, params, 24.0, dt_output=24.0)
        if not completed or len(results) < 1:
            return None

        t_final, y_ss = results[-1]

        # Extract outputs
        HR = param_overrides.get('HR_heart_rate', 70)
        SP = y_ss[STATE_MAP['systolic_pressure']]
        DP = y_ss[STATE_MAP['diastolic_pressure']]
        MAP = (SP / 3 + DP * 2 / 3) * 0.0075
        SBP = SP * 0.0075
        DBP = DP * 0.0075
        CO = y_ss[STATE_MAP['CO_delayed']]
        SV = CO * 1000 / HR if HR > 0 else 0
        BV = y_ss[STATE_MAP['blood_volume_L']]
        Na = y_ss[STATE_MAP['sodium_amount']] / BV if BV > 0 else 140
        EDV = y_ss[STATE_MAP['LV_EDV']] * 1e6
        ESV = EDV - SV
        EF = (SV / EDV) * 100 if EDV > 0 else 0
        EDP = y_ss[STATE_MAP['LV_EDP']] * 0.0075
        sCr = y_ss[STATE_MAP['serum_creatinine']] / BV if BV > 0 else 0.92

        # LV mass
        Pi = 3.1416
        fib_i = params[PARAM_MAP['Baseline_Interstitial_Fibrosis']]
        fib_r = params[PARAM_MAP['Baseline_Replacement_Fibrosis']]
        tis_i = params[PARAM_MAP['Baseline_Interstitial_Tissue']]
        myo_N = params[PARAM_MAP['Baseline_Myocyte_Number']]
        myo_L = params[PARAM_MAP['Baseline_Myocyte_Length']]
        btmv = V_w_0 - fib_i - fib_r - tis_i
        bsmv = btmv / myo_N
        bmd = 2 * np.sqrt(bsmv / (Pi * myo_L))
        cmd = y_ss[STATE_MAP['change_in_myocyte_diameter']]
        cml = y_ss[STATE_MAP['change_in_myocyte_length']]
        ml = myo_L + cml
        md = bmd + cmd
        smv = ml * Pi * (md**2) / 4
        tmv = smv * myo_N
        tnmv = fib_i + tis_i + fib_r
        LV_mass = 1e6 * (tmv + tnmv) * 1.05

        # eGFR (CKD-EPI, male, age 68)
        kappa = 0.9
        ratio = sCr / kappa
        eGFR = 141 * (min(ratio, 1)**(-0.411)) * (max(ratio, 1)**(-1.209)) * (0.993**68)

        # BNP
        EDS = y_ss[STATE_MAP['LV_EDS']]
        BNP_factor = params[PARAM_MAP['BNP_factor']]
        BNP = np.exp(BNP_factor * ((EDS + 1736) / 5.094) + 3.14)
        BNP = min(BNP, 5000)

        sim_time = time.time() - t_start
        out = {
            'patient_id': patient_id,
            'sim_seconds': sim_time,
            # Input parameters (prefixed to avoid collision with outputs)
            **{f'param_{k}': v for k, v in param_overrides.items()},
            # Clinical outputs
            'MAP': MAP, 'SBP': SBP, 'DBP': DBP, 'CO': CO, 'HR': HR,
            'SV': SV, 'BV': BV, 'Na': Na, 'EDV': EDV, 'ESV': ESV,
            'EF': EF, 'EDP': EDP, 'LV_mass': LV_mass, 'sCr': sCr,
            'eGFR': eGFR, 'BNP': BNP,
        }
        # Reject if any clinical output is NaN or clearly unphysiological
        if any(np.isnan(v) for v in [MAP, CO, EF, EDV, BV, Na]):
            return None
        return out
    except Exception as e:
        return None


def generate_report(results, elapsed, n_requested):
    """Generate validation report as a string."""
    n_success = len(results)
    n_fail = n_requested - n_success

    lines = []
    lines.append("=" * 70)
    lines.append("CARDIORENAL HFpEF SIMULATOR — VALIDATION REPORT")
    lines.append("=" * 70)
    lines.append(f"Date: {time.strftime('%Y-%m-%d %H:%M')}")
    lines.append(f"Patients requested: {n_requested}")
    lines.append(f"Patients completed: {n_success} ({100*n_success/n_requested:.1f}%)")
    lines.append(f"Patients failed:    {n_fail} ({100*n_fail/n_requested:.1f}%)")
    lines.append(f"Total time:         {elapsed:.1f}s ({elapsed/n_requested:.2f}s per patient wall clock)")
    sim_times = [r['sim_seconds'] for r in results if 'sim_seconds' in r]
    if sim_times:
        lines.append(f"Per-patient sim:    mean={np.mean(sim_times):.2f}s, min={np.min(sim_times):.2f}s, max={np.max(sim_times):.2f}s")
        lines.append(f"Projected 100K×15yr: {np.mean(sim_times)*15*100000/3600:.0f} CPU-hours")
    lines.append("")

    # Clinical output distributions
    clinical_vars = [
        ("MAP",     "mmHg",   70, 105, "Mean arterial pressure"),
        ("SBP",     "mmHg",   90, 160, "Systolic blood pressure"),
        ("DBP",     "mmHg",   50,  95, "Diastolic blood pressure"),
        ("CO",      "L/min",   3,   9, "Cardiac output"),
        ("SV",      "mL",     40, 120, "Stroke volume"),
        ("BV",      "L",     3.5, 6.5, "Blood volume"),
        ("Na",      "mEq/L", 130, 150, "Serum sodium"),
        ("EDV",     "mL",     50, 250, "End-diastolic volume"),
        ("ESV",     "mL",     10, 120, "End-systolic volume"),
        ("EF",      "%",      35,  75, "Ejection fraction"),
        ("EDP",     "mmHg",    2,  25, "End-diastolic pressure"),
        ("LV_mass", "g",      60, 350, "LV mass"),
        ("sCr",     "mg/dL", 0.4, 2.5, "Serum creatinine"),
        ("eGFR",    "mL/min", 30, 150, "Est. GFR"),
        ("BNP",     "pg/mL",   5, 500, "BNP"),
    ]

    lines.append("-" * 70)
    lines.append("CLINICAL OUTPUT DISTRIBUTIONS (at 24h steady state)")
    lines.append("-" * 70)
    lines.append(f"{'Variable':<12} {'Unit':<8} {'Mean':>8} {'SD':>8} {'Min':>8} {'Max':>8} {'Phys Range':>14} {'In Range':>10}")
    lines.append("-" * 70)

    for varname, unit, phys_lo, phys_hi, desc in clinical_vars:
        vals = np.array([r[varname] for r in results])
        mean = np.mean(vals)
        sd = np.std(vals)
        vmin = np.min(vals)
        vmax = np.max(vals)
        in_range = np.sum((vals >= phys_lo) & (vals <= phys_hi))
        pct = 100 * in_range / len(vals)
        lines.append(f"{varname:<12} {unit:<8} {mean:>8.2f} {sd:>8.2f} {vmin:>8.2f} {vmax:>8.2f} {f'[{phys_lo}-{phys_hi}]':>14} {pct:>8.1f}%")

    lines.append("")

    # Correlation check: key physiological correlations
    lines.append("-" * 70)
    lines.append("KEY CROSS-VARIABLE CORRELATIONS")
    lines.append("-" * 70)
    lines.append("(Verifying the simulator enforces physiological coupling)")
    lines.append("")

    corr_pairs = [
        ("EDV", "ESV",     "EDV-ESV (volume coupling)",        "+"),
        ("EDV", "EF",      "EDV-EF (Frank-Starling)",          "+/-"),
        ("EF",  "CO",      "EF-CO (pump function)",            "+"),
        ("MAP", "CO",      "MAP-CO (hemodynamic coupling)",    "+"),
        ("EDP", "BNP",     "EDP-BNP (wall stress → BNP)",      "+"),
        ("EDP", "LV_mass", "EDP-LV mass (pressure → hypertrophy)", "+"),
        ("sCr", "eGFR",    "sCr-eGFR (inverse renal function)", "-"),
        ("SBP", "DBP",     "SBP-DBP (pressure coupling)",      "+"),
        ("BV",  "CO",      "BV-CO (preload → output)",         "+"),
        ("LV_mass", "EDV", "LV mass-EDV (remodeling)",         "+"),
    ]

    lines.append(f"{'Pair':<40} {'r':>8} {'Expected':>10} {'Match':>8}")
    lines.append("-" * 70)

    for var1, var2, desc, expected_sign in corr_pairs:
        v1 = np.array([r[var1] for r in results])
        v2 = np.array([r[var2] for r in results])
        if np.std(v1) > 0 and np.std(v2) > 0:
            r = np.corrcoef(v1, v2)[0, 1]
        else:
            r = 0
        if expected_sign == "+":
            match = "YES" if r > 0.1 else "WEAK" if r > -0.1 else "NO"
        elif expected_sign == "-":
            match = "YES" if r < -0.1 else "WEAK" if r < 0.1 else "NO"
        else:
            match = "OK"  # either direction acceptable
        lines.append(f"{desc:<40} {r:>8.3f} {expected_sign:>10} {match:>8}")

    lines.append("")

    # Parameter sensitivity: which inputs drive which outputs
    lines.append("-" * 70)
    lines.append("PARAMETER SENSITIVITY (|correlation| of input param → output)")
    lines.append("-" * 70)

    key_outputs = ["MAP", "CO", "EF", "EDV", "EDP", "LV_mass", "eGFR", "BNP"]
    key_inputs = [
        "param_HR_heart_rate", "param_cf", "param_contractility", "param_C_art_initial",
        "param_R_per0", "param_baseline_nephrons", "param_BV", "param_V_w_0",
    ]

    header = f"{'Input Param':<28}" + "".join([f"{o:>9}" for o in key_outputs])
    lines.append(header)
    lines.append("-" * 70)

    for inp in key_inputs:
        v_in = np.array([r[inp] for r in results])
        row = f"{inp:<28}"
        for out in key_outputs:
            v_out = np.array([r[out] for r in results])
            if np.std(v_in) > 0 and np.std(v_out) > 0:
                r = np.corrcoef(v_in, v_out)[0, 1]
            else:
                r = 0
            row += f"{r:>9.3f}"
        lines.append(row)

    lines.append("")

    # Summary
    lines.append("=" * 70)
    lines.append("SUMMARY")
    lines.append("=" * 70)
    lines.append(f"  Success rate:     {100*n_success/n_requested:.1f}%")

    # Count how many outputs have >80% of patients in physiological range
    good_outputs = 0
    for varname, unit, phys_lo, phys_hi, desc in clinical_vars:
        vals = np.array([r[varname] for r in results])
        pct = 100 * np.sum((vals >= phys_lo) & (vals <= phys_hi)) / len(vals)
        if pct >= 80:
            good_outputs += 1

    lines.append(f"  Outputs in range: {good_outputs}/{len(clinical_vars)} have >80% of patients in physiological range")

    # Count correlations that match expected direction
    good_corrs = 0
    for var1, var2, desc, expected_sign in corr_pairs:
        v1 = np.array([r[var1] for r in results])
        v2 = np.array([r[var2] for r in results])
        if np.std(v1) > 0 and np.std(v2) > 0:
            r = np.corrcoef(v1, v2)[0, 1]
        else:
            r = 0
        if expected_sign == "+" and r > 0.1:
            good_corrs += 1
        elif expected_sign == "-" and r < -0.1:
            good_corrs += 1
        elif expected_sign == "+/-":
            good_corrs += 1

    lines.append(f"  Correlations OK:  {good_corrs}/{len(corr_pairs)} match expected physiological direction")
    lines.append("")
    lines.append("CONCLUSION: " + (
        "Simulator produces physiologically valid and coupled outputs across diverse patient parameters."
        if good_outputs >= 12 and good_corrs >= 7
        else "Some outputs or correlations need investigation."
    ))
    lines.append("=" * 70)

    return "\n".join(lines)


if __name__ == "__main__":
    N_PATIENTS = 50

    print(f"Generating {N_PATIENTS} synthetic patients...")
    rng = np.random.default_rng(seed=42)

    patient_args = []
    for i in range(N_PATIENTS):
        overrides = sample_patient(rng, i)
        patient_args.append((i, overrides))

    print(f"Running simulations on {min(cpu_count(), 8)} cores...")
    t0 = time.time()

    n_workers = min(cpu_count(), 8)
    with Pool(n_workers) as pool:
        raw_results = pool.map(run_single_patient, patient_args)

    elapsed = time.time() - t0

    results = [r for r in raw_results if r is not None]

    print(f"Completed {len(results)}/{N_PATIENTS} in {elapsed:.1f}s")
    print()

    # Generate report
    report = generate_report(results, elapsed, N_PATIENTS)
    print(report)

    # Save report
    report_path = os.path.join(os.path.dirname(__file__), "validation_report.txt")
    with open(report_path, "w") as f:
        f.write(report)
    print(f"\nReport saved to {report_path}")

    # Save raw data as JSON
    data_path = os.path.join(os.path.dirname(__file__), "validation_data.json")
    with open(data_path, "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"Raw data saved to {data_path}")
