"""
Generate training data for the Neural ODE surrogate.

Runs Latin Hypercube sampling over the knob space, integrates each
parameter set with the C solver, and extracts clinical trajectories.
Saves to HDF5 for training.

Usage:
    python -m neural_surrogate.generate_training_data [--n_samples 5000] [--t_hours 168]
"""

import numpy as np
import h5py
import time
import os
import sys
import argparse
import math
import signal
import multiprocessing as mp

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from neural_surrogate.model import (
    CLINICAL_KEYS, CRITICAL_STATE_INDICES, N_CLINICAL, N_CRITICAL
)

# Lazy imports — these load the C library, so only import in main/worker processes
_driver_loaded = False
def _ensure_driver():
    global _driver_loaded, build_params, build_inits, integrate
    global STATE_MAP, PARAM_MAP, N_STATE, N_PARAM
    if not _driver_loaded:
        from hallow_c_driver import (
            build_params, build_inits, integrate,
            STATE_MAP, PARAM_MAP, N_STATE, N_PARAM
        )
        _driver_loaded = True

# Import for top-level use (state map needed for extract functions)
from hallow_c_driver import STATE_MAP, PARAM_MAP, N_STATE, N_PARAM

# Knob ranges for Latin Hypercube Sampling
KNOB_RANGES = {
    'C_art_scale':  (0.3, 1.5),
    'TPR_mult':     (0.5, 2.5),
    'a1c':          (4.0, 12.0),
    'nephron_loss': (0.0, 0.6),
}
KNOB_NAMES = list(KNOB_RANGES.keys())
N_KNOBS = len(KNOB_NAMES)


def latin_hypercube(n_samples, n_dims, rng=None):
    """Simple LHS: stratified sampling in each dimension."""
    if rng is None:
        rng = np.random.default_rng(42)
    result = np.zeros((n_samples, n_dims))
    for d in range(n_dims):
        perm = rng.permutation(n_samples)
        for i in range(n_samples):
            result[perm[i], d] = (i + rng.random()) / n_samples
    return result


def knobs_to_params(base_params, r_values, knob_vec):
    """Apply a knob vector [C_art_scale, TPR_mult, a1c, nephron_loss] to params."""
    params = base_params.copy()
    c_art, tpr, a1c, nephron_loss = knob_vec

    idx = PARAM_MAP.get('C_art_scale')
    if idx is not None:
        params[idx] = c_art

    idx = PARAM_MAP.get('disease_effect_on_TPR_peripheral_resistance')
    if idx is not None:
        params[idx] = tpr

    glucose = max(3.0, 1.59 * a1c - 2.59)
    idx = PARAM_MAP.get('glucose_concentration')
    if idx is not None:
        params[idx] = glucose

    idx = PARAM_MAP.get('disease_effect_on_nephrons')
    if idx is not None:
        params[idx] = nephron_loss

    return params


def extract_clinical_from_state(y):
    """Extract clinical outputs from a raw 70-dim state vector.
    Mirrors server.py:extract_outputs but returns a flat numpy array."""
    SP = y[STATE_MAP['systolic_pressure']]
    DP = y[STATE_MAP['diastolic_pressure']]
    MAP = (SP / 3 + DP * 2 / 3) * 0.0075
    SBP = SP * 0.0075
    DBP = DP * 0.0075
    CO = y[STATE_MAP['CO_delayed']]
    BV = y[STATE_MAP['blood_volume_L']]
    Na = y[STATE_MAP['sodium_amount']] / BV if BV > 0 else 140.0
    EDV = y[STATE_MAP['LV_EDV']] * 1e6
    EDP = y[STATE_MAP['LV_EDP']] * 0.0075
    EDS = y[STATE_MAP['LV_EDS']]
    ASP = y[STATE_MAP['LV_active_stress_peak']]
    cmd = y[STATE_MAP['change_in_myocyte_diameter']]
    cml = y[STATE_MAP['change_in_myocyte_length']]
    sCr = y[STATE_MAP['serum_creatinine']] / BV if BV > 0 else 0.92

    # LV mass (same as server.py)
    Pi = 3.1416
    V_w_0 = 0.00012
    btmv = V_w_0 - V_w_0 * 0.02 - V_w_0 * 0.02 - V_w_0 * 0.22
    bsmv = btmv / 3.3e9
    bmd = 2 * math.sqrt(bsmv / (Pi * 0.000115))
    ml = 0.000115 + cml
    md = bmd + cmd
    smv = ml * Pi * (md ** 2) / 4
    tmv = smv * 3.3e9
    tnmv = V_w_0 * 0.02 + V_w_0 * 0.22 + V_w_0 * 0.02
    LV_mass = 1e6 * (tmv + tnmv) * 1.05

    # BNP
    BNP = math.exp(0.0008 * ((EDS + 1736) / 5.094) + 3.14)
    BNP = min(BNP, 5000.0)

    # Percentage changes
    pct_d = 100 * cmd / bmd if bmd > 0 else 0.0
    pct_l = 100 * cml / 0.000115 if 0.000115 > 0 else 0.0

    return np.array([
        MAP, SBP, DBP, CO,
        BV, Na, EDV, EDP,
        EDS, ASP, LV_mass, BNP,
        sCr, pct_d, pct_l,
    ], dtype=np.float64)


def extract_critical_state(y):
    """Extract the critical state subset needed for full-state reconstruction."""
    return y[CRITICAL_STATE_INDICES].copy()


def run_single(base_params, r_values, y0_base, knob_vec, t_hours, dt_output):
    """Run one simulation, return (times, clinical_traj, state_traj, full_state_traj, success)."""
    _ensure_driver()
    params = knobs_to_params(base_params, r_values, knob_vec)
    y0 = y0_base.copy()

    try:
        results, completed = integrate(y0, params, t_hours, dt_output=dt_output)
    except Exception as e:
        print(f"  Integration error: {e}")
        return None, None, None, None, False

    if not completed or len(results) < 2:
        return None, None, None, None, False

    n_pts = len(results)
    times = np.zeros(n_pts)
    clinical = np.zeros((n_pts, N_CLINICAL))
    crit_state = np.zeros((n_pts, N_CRITICAL))
    full_state = np.zeros((n_pts, N_STATE))

    for i, (t, y) in enumerate(results):
        times[i] = t
        if np.any(np.isnan(y)) or np.any(np.isinf(y)):
            return None, None, None, None, False
        clinical[i] = extract_clinical_from_state(y)
        crit_state[i] = extract_critical_state(y)
        full_state[i] = y

    if np.any(np.isnan(clinical)) or np.any(np.isinf(clinical)):
        return None, None, None, None, False

    return times, clinical, crit_state, full_state, True


def _subprocess_target(knob_vec, t_hours, dt_output, result_dict):
    """Target for subprocess — writes result into shared dict."""
    try:
        _ensure_driver()
        base_params, r_values = build_params()
        y0_base = build_inits(r_values)
        times, clinical, crit_state, full_state, ok = run_single(
            base_params, r_values, y0_base, knob_vec, t_hours, dt_output)
        if ok:
            result_dict['data'] = (knob_vec, times, clinical, crit_state,
                                   full_state[0], full_state[-1])
    except Exception:
        pass


def run_with_timeout(knob_vec, t_hours, dt_output, timeout=60):
    """Run a single simulation in a subprocess with hard timeout.
    Returns result tuple or None if failed/timed out."""
    manager = mp.Manager()
    result_dict = manager.dict()

    p = mp.Process(target=_subprocess_target,
                   args=(knob_vec, t_hours, dt_output, result_dict))
    p.start()
    p.join(timeout=timeout)

    if p.is_alive():
        p.kill()
        p.join(timeout=5)
        return None

    if p.exitcode != 0:
        return None

    return result_dict.get('data', None)


def generate_dataset(n_samples=5000, t_hours=168, dt_output=1.0, seed=42,
                     output_path=None, n_workers=None, timeout=120):
    """Generate the full training dataset.

    Each simulation runs in a subprocess pool to isolate C-level crashes.
    """
    if output_path is None:
        output_path = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), 'training_data.h5')

    if n_workers is None:
        n_workers = min(mp.cpu_count(), 4)

    rng = np.random.default_rng(seed)
    print(f"Generating {n_samples} samples, {t_hours}h each, dt={dt_output}h")
    print(f"Workers: {n_workers}, timeout: {timeout}s per sim")
    print(f"Output: {output_path}")

    # Latin Hypercube Sampling
    lhs = latin_hypercube(n_samples, N_KNOBS, rng)
    knob_samples = np.zeros((n_samples, N_KNOBS))
    for d, name in enumerate(KNOB_NAMES):
        lo, hi = KNOB_RANGES[name]
        knob_samples[:, d] = lo + lhs[:, d] * (hi - lo)

    n_timepoints = int(t_hours / dt_output) + 1

    # Pre-allocate HDF5
    with h5py.File(output_path, 'w') as f:
        f.create_dataset('knobs', shape=(0, N_KNOBS), maxshape=(None, N_KNOBS),
                         dtype='f4', chunks=True)
        f.create_dataset('times', shape=(0, n_timepoints), maxshape=(None, n_timepoints),
                         dtype='f4', chunks=True)
        f.create_dataset('clinical', shape=(0, n_timepoints, N_CLINICAL),
                         maxshape=(None, n_timepoints, N_CLINICAL), dtype='f4', chunks=True)
        f.create_dataset('critical_state', shape=(0, n_timepoints, N_CRITICAL),
                         maxshape=(None, n_timepoints, N_CRITICAL), dtype='f4', chunks=True)
        f.create_dataset('initial_state', shape=(0, N_STATE), maxshape=(None, N_STATE),
                         dtype='f4', chunks=True)
        f.create_dataset('final_state', shape=(0, N_STATE), maxshape=(None, N_STATE),
                         dtype='f4', chunks=True)

        f.attrs['n_state'] = N_STATE
        f.attrs['n_knobs'] = N_KNOBS
        f.attrs['n_clinical'] = N_CLINICAL
        f.attrs['n_critical'] = N_CRITICAL
        f.attrs['clinical_keys'] = [k.encode() for k in CLINICAL_KEYS]
        f.attrs['knob_names'] = [k.encode() for k in KNOB_NAMES]
        f.attrs['t_hours'] = t_hours
        f.attrs['dt_output'] = dt_output

    # Run simulations sequentially, each in its own subprocess with timeout
    n_success = 0
    n_failed = 0
    t_start = time.time()

    for i in range(n_samples):
        knob_vec = knob_samples[i]
        result = run_with_timeout(knob_vec, t_hours, dt_output, timeout=timeout)

        if result is not None:
            kv, times, clinical, crit_state, init_state, final_state = result
            if len(times) == n_timepoints:
                with h5py.File(output_path, 'a') as f:
                    idx = f['knobs'].shape[0]
                    for ds_name in ['knobs', 'times', 'clinical', 'critical_state',
                                    'initial_state', 'final_state']:
                        f[ds_name].resize(idx + 1, axis=0)
                    f['knobs'][idx] = kv
                    f['times'][idx] = times
                    f['clinical'][idx] = clinical
                    f['critical_state'][idx] = crit_state
                    f['initial_state'][idx] = init_state
                    f['final_state'][idx] = final_state
                n_success += 1
            else:
                n_failed += 1
        else:
            n_failed += 1

        if (i + 1) % 50 == 0 or (i + 1) == n_samples:
            elapsed = time.time() - t_start
            rate = (i + 1) / elapsed
            eta = (n_samples - i - 1) / rate if rate > 0 else 0
            print(f"  [{i+1}/{n_samples}] {n_success} ok, {n_failed} failed, "
                  f"{elapsed:.0f}s elapsed, ETA {eta:.0f}s")

    elapsed = time.time() - t_start
    print(f"\nComplete: {n_success}/{n_samples} successful in {elapsed:.0f}s")
    print(f"Dataset saved to {output_path}")

    with h5py.File(output_path, 'r') as f:
        print(f"\nDataset summary:")
        print(f"  Samples: {f['knobs'].shape[0]}")
        print(f"  Time points per segment: {n_timepoints}")
        print(f"  Clinical outputs: {N_CLINICAL}")
        print(f"  Critical state vars: {N_CRITICAL}")

    return output_path


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--n_samples', type=int, default=5000)
    parser.add_argument('--t_hours', type=float, default=168)
    parser.add_argument('--dt_output', type=float, default=1.0)
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--output', type=str, default=None)
    args = parser.parse_args()

    generate_dataset(
        n_samples=args.n_samples,
        t_hours=args.t_hours,
        dt_output=args.dt_output,
        seed=args.seed,
        output_path=args.output,
    )
