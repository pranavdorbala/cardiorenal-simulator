"""
Morris Elementary Effects screening for the Hallow cardiorenal model.

Identifies which of the 430 parameters are influential for each clinical
output using One-At-a-Time (OAT) perturbations across r trajectories.

Cost: r * (k+1) simulations, where k = number of perturbable parameters.
With k~350 params and r=20 trajectories: ~7,020 simulations.
"""

import numpy as np
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from hallow_c_driver import (
    build_params, build_inits, batch_integrate,
    STATE_MAP, PARAM_MAP, N_STATE, N_PARAM
)
from neural_surrogate.generate_training_data import extract_clinical_from_state

# Clinical output names (must match extract_clinical_from_state order)
CLINICAL_KEYS = [
    'MAP', 'SBP', 'DBP', 'CO', 'blood_volume_L', 'Na',
    'LV_EDV_mL', 'LV_EDP_mmHg', 'LV_EDS', 'LV_active_stress_peak',
    'LV_mass', 'BNP', 'serum_creatinine', 'pct_diameter', 'pct_length',
]

# Parameter indices to skip: unit conversions and mathematical constants
SKIP_INDICES = set(range(17))  # p[0]-p[16]: nL_mL, dl_ml, ..., Pi


def identify_perturbable_params(param_map, base_params, r_values=None):
    """Filter 430 params to those worth perturbing.

    Excludes: unit conversions (indices 0-16), zero-valued params with
    no physiological meaning, boolean on/off drug flags.

    Returns: list of dicts with keys: name, index, nominal, lo, hi
    """
    # Boolean/flag params that should be 0 or 1
    flag_params = {
        'heart_renal_link', 'aortic_valve_stenosis', 'mitral_regurgitation',
        'ANP_infusion', 'aortic_regurgitation',
    }

    # Rate constants that should be perturbed on log scale
    rate_params = {
        'C_renal_CV_timescale', 'C_cycle', 'C_cycle2', 'C_cycle3',
        'C_co', 'C_co_delay', 'C_map', 'C_co_error',
        'C_vasopressin_delay', 'C_Na_error', 'C_aldo_secretion',
        'C_tgf_reset', 'C_md_flow', 'C_tgf', 'C_rbf',
        'C_serum_creatinine', 'C_pt_water', 'C_rsna',
        'C_postglomerular_pressure', 'C_sglt2_delay', 'C_ruge',
    }

    perturbable = []
    for name, idx in sorted(param_map.items(), key=lambda x: x[1]):
        if idx in SKIP_INDICES:
            continue
        if name in flag_params:
            continue

        nominal = base_params[idx]

        # Skip zero-valued params that are likely unused drug effects
        if nominal == 0.0:
            # Some zero-valued params are meaningful (e.g., disease effects start at 0)
            meaningful_zeros = {
                'disease_effects_decreasing_Kf', 'disease_effect_on_nephrons',
                'SGLT2_inhibition', 'SGLT1_inhibition', 'loop_diuretic_effect',
                'CD_PN_loss_rate',
            }
            if name not in meaningful_zeros:
                continue

        # Set perturbation bounds
        if name in rate_params:
            # Log-uniform: 0.1x to 10x
            lo = nominal * 0.1 if nominal > 0 else 0.0
            hi = nominal * 10.0 if nominal > 0 else 1.0
        elif nominal == 0.0:
            # Meaningful zero params: small positive range
            lo = 0.0
            hi = 0.5
        elif nominal > 0:
            lo = nominal * 0.5
            hi = nominal * 1.5
        else:
            lo = nominal * 1.5  # negative: 1.5x is more negative
            hi = nominal * 0.5

        perturbable.append({
            'name': name,
            'index': idx,
            'nominal': nominal,
            'lo': lo,
            'hi': hi,
        })

    return perturbable


def generate_morris_trajectories(n_params, r=20, levels=4, seed=42):
    """Generate r Morris trajectories in [0,1]^k space.

    Each trajectory has k+1 points (base + one perturbation per param).
    Uses the simplified one-at-a-time design.

    Returns: array of shape (r, k+1, k) in unit hypercube.
    """
    rng = np.random.default_rng(seed)
    delta = levels / (2 * (levels - 1))  # step size

    trajectories = np.zeros((r, n_params + 1, n_params))

    for traj_idx in range(r):
        # Random base point on grid
        base = rng.choice(np.linspace(0, 1 - delta, levels), size=n_params)
        trajectories[traj_idx, 0] = base

        # Random permutation of parameter order
        order = rng.permutation(n_params)

        current = base.copy()
        for step, param_idx in enumerate(order):
            current = current.copy()
            # Perturb param_idx by +delta or -delta
            if current[param_idx] + delta <= 1.0:
                current[param_idx] += delta
            else:
                current[param_idx] -= delta
            trajectories[traj_idx, step + 1] = current

    return trajectories


def scale_trajectories(trajectories, perturbable):
    """Map [0,1]^k trajectories to full 430-dim parameter vectors.

    Args:
        trajectories: shape (r, k+1, k) in unit hypercube
        perturbable: list of perturbable param specs

    Returns:
        list of 430-dim parameter arrays, one per trajectory point
    """
    base_params, _ = build_params()
    r, steps, k = trajectories.shape

    param_sets = []
    for traj_idx in range(r):
        for step in range(steps):
            p = base_params.copy()
            x = trajectories[traj_idx, step]
            for i, spec in enumerate(perturbable):
                p[spec['index']] = spec['lo'] + x[i] * (spec['hi'] - spec['lo'])
            param_sets.append(p)

    return param_sets


def compute_elementary_effects(clinical_outputs, trajectories, perturbable):
    """Compute mean (mu*) and std (sigma) of elementary effects.

    Args:
        clinical_outputs: array shape (r*(k+1), n_clinical)
        trajectories: shape (r, k+1, k)
        perturbable: list of perturbable param specs

    Returns: dict per clinical output with ranked params.
        {output_name: [{param, mu_star, sigma, rank}, ...]}
    """
    r, steps, k = trajectories.shape
    n_clinical = clinical_outputs.shape[1]
    delta = trajectories[0, 1] - trajectories[0, 0]  # not exactly right, compute per step

    # Compute elementary effects for each param
    # EE_i = (f(x + delta_i) - f(x)) / delta
    effects = {i: [] for i in range(k)}  # param_idx -> list of EE vectors

    for traj_idx in range(r):
        base_idx_start = traj_idx * steps
        prev_point = trajectories[traj_idx, 0]

        for step in range(1, steps):
            curr_point = trajectories[traj_idx, step]
            diff = curr_point - prev_point
            changed_params = np.where(np.abs(diff) > 1e-10)[0]

            if len(changed_params) == 1:
                param_i = changed_params[0]
                d = diff[param_i]
                if abs(d) > 1e-10:
                    f_curr = clinical_outputs[base_idx_start + step]
                    f_prev = clinical_outputs[base_idx_start + step - 1]
                    ee = (f_curr - f_prev) / d
                    effects[param_i].append(ee)

            prev_point = curr_point

    # Compute statistics
    results = {}
    for c_idx, c_name in enumerate(CLINICAL_KEYS):
        param_stats = []
        for p_idx in range(k):
            ees = [e[c_idx] for e in effects[p_idx] if not np.isnan(e[c_idx])]
            if len(ees) > 0:
                mu_star = np.mean(np.abs(ees))
                sigma = np.std(ees)
            else:
                mu_star = 0.0
                sigma = 0.0
            param_stats.append({
                'param': perturbable[p_idx]['name'],
                'param_index': perturbable[p_idx]['index'],
                'mu_star': float(mu_star),
                'sigma': float(sigma),
            })

        # Rank by mu_star
        param_stats.sort(key=lambda x: x['mu_star'], reverse=True)
        for rank, ps in enumerate(param_stats):
            ps['rank'] = rank + 1

        results[c_name] = param_stats

    return results


def run_morris(r=20, t_hours=168, n_workers=None, levels=4, seed=42):
    """Full Morris analysis pipeline.

    1. Identify perturbable params (~350 of 430)
    2. Generate r trajectories -> r*(k+1) param vectors
    3. Run batch_integrate() on all param vectors (parallel)
    4. Extract 15 clinical outputs from final states
    5. Compute elementary effects
    6. Return ranked importance per output

    Returns:
        dict with keys: results, perturbable, metadata
    """
    import time
    t0 = time.time()

    # Step 1: Setup
    base_params, r_values = build_params()
    y0 = build_inits(r_values)
    perturbable = identify_perturbable_params(PARAM_MAP, base_params, r_values)
    k = len(perturbable)
    print(f"[Morris] {k} perturbable parameters (of {N_PARAM})")

    # Step 2: Generate trajectories
    trajectories = generate_morris_trajectories(k, r=r, levels=levels, seed=seed)
    n_sims = r * (k + 1)
    print(f"[Morris] {r} trajectories, {n_sims} total simulations")

    # Step 3: Scale to parameter vectors
    param_sets = scale_trajectories(trajectories, perturbable)

    # Step 4: Run simulations
    print(f"[Morris] Running {n_sims} simulations ({t_hours}h each)...")
    batch_results = batch_integrate(y0, param_sets, t_end=t_hours,
                                     dt_output=max(1.0, t_hours / 10),
                                     n_workers=n_workers)

    # Step 5: Extract clinical outputs
    n_clinical = len(CLINICAL_KEYS)
    clinical_outputs = np.full((n_sims, n_clinical), np.nan)

    n_ok = 0
    for idx, success, final_state in batch_results:
        if success and final_state is not None:
            try:
                clinical_outputs[idx] = extract_clinical_from_state(final_state)
                n_ok += 1
            except Exception:
                pass

    print(f"[Morris] {n_ok}/{n_sims} simulations completed successfully")

    # Step 6: Compute elementary effects
    results = compute_elementary_effects(clinical_outputs, trajectories, perturbable)

    elapsed = time.time() - t0

    # Find top params across all outputs (union of top-20 per output)
    top_set = set()
    for output_name, stats in results.items():
        for s in stats[:20]:
            top_set.add(s['param'])

    return {
        'results': results,
        'perturbable': perturbable,
        'top_params': sorted(top_set),
        'metadata': {
            'n_perturbable': k,
            'n_sims': n_sims,
            'n_ok': n_ok,
            'r': r,
            'levels': levels,
            't_hours': t_hours,
            'elapsed_seconds': round(elapsed, 1),
        }
    }
