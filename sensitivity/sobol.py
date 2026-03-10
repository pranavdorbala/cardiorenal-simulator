"""
Sobol variance decomposition for the Hallow cardiorenal model.

Computes first-order (S1) and total-order (ST) Sobol indices for the
top parameters identified by Morris screening. Uses Saltelli sampling.

Cost: n_samples * (2*k + 2) simulations.
For k=50 params, n=1024 -> ~104K simulations (HPC only).
"""

import numpy as np
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from hallow_c_driver import (
    build_params, build_inits, batch_integrate,
    PARAM_MAP, N_PARAM
)
from neural_surrogate.generate_training_data import extract_clinical_from_state
from sensitivity.morris import CLINICAL_KEYS


def generate_sobol_samples(top_params, base_params, n_samples=1024, seed=42):
    """Generate Saltelli sampling matrices for Sobol analysis.

    Uses scipy.stats.qmc.Sobol for quasi-random sequences.

    Args:
        top_params: list of dicts with keys: name, index, lo, hi
        base_params: 430-dim baseline parameter array
        n_samples: base sample size (power of 2 recommended)
        seed: random seed

    Returns:
        (param_sets, sample_info) where:
        - param_sets: list of 430-dim parameter arrays
        - sample_info: dict with n_A, n_B, matrix structure for index computation
    """
    from scipy.stats.qmc import Sobol as SobolQMC

    k = len(top_params)
    # Generate 2 independent Sobol sequences (A and B matrices)
    sampler = SobolQMC(d=k, scramble=True, seed=seed)
    AB_raw = sampler.random(n_samples * 2)  # shape (2*n, k) in [0,1]^k
    A_raw = AB_raw[:n_samples]   # matrix A
    B_raw = AB_raw[n_samples:]   # matrix B

    def raw_to_params(raw_matrix):
        """Convert [0,1]^k matrix to list of 430-dim parameter arrays."""
        param_list = []
        for row in raw_matrix:
            p = base_params.copy()
            for i, spec in enumerate(top_params):
                p[spec['index']] = spec['lo'] + row[i] * (spec['hi'] - spec['lo'])
            param_list.append(p)
        return param_list

    # Collect all parameter sets: A, B, and AB_i matrices
    param_sets = []
    param_sets.extend(raw_to_params(A_raw))  # A: indices [0, n)
    param_sets.extend(raw_to_params(B_raw))  # B: indices [n, 2n)

    # AB_i matrices: A with column i replaced by B's column i
    ab_offsets = []
    for i in range(k):
        AB_i = A_raw.copy()
        AB_i[:, i] = B_raw[:, i]
        offset = len(param_sets)
        param_sets.extend(raw_to_params(AB_i))
        ab_offsets.append(offset)

    sample_info = {
        'n_samples': n_samples,
        'k': k,
        'a_offset': 0,
        'b_offset': n_samples,
        'ab_offsets': ab_offsets,
        'total_sims': len(param_sets),
    }

    return param_sets, sample_info


def compute_sobol_indices(clinical_outputs, sample_info):
    """Compute first-order (S1) and total-order (ST) Sobol indices.

    Uses Jansen estimator for ST and Saltelli estimator for S1.

    Args:
        clinical_outputs: array shape (total_sims, n_clinical)
        sample_info: from generate_sobol_samples

    Returns:
        dict per clinical output: {output_name: {S1: [...], ST: [...]}}
    """
    n = sample_info['n_samples']
    k = sample_info['k']
    a_off = sample_info['a_offset']
    b_off = sample_info['b_offset']
    ab_offsets = sample_info['ab_offsets']

    n_clinical = clinical_outputs.shape[1]

    # Extract A and B output matrices
    y_A = clinical_outputs[a_off:a_off + n]  # (n, n_clinical)
    y_B = clinical_outputs[b_off:b_off + n]  # (n, n_clinical)

    results = {}
    for c_idx, c_name in enumerate(CLINICAL_KEYS):
        fA = y_A[:, c_idx]
        fB = y_B[:, c_idx]

        # Total variance
        f_all = np.concatenate([fA, fB])
        f_all = f_all[~np.isnan(f_all)]
        if len(f_all) < 2:
            results[c_name] = {'S1': [], 'ST': []}
            continue
        var_total = np.var(f_all)
        if var_total < 1e-20:
            results[c_name] = {'S1': [0.0] * k, 'ST': [0.0] * k}
            continue

        S1_list = []
        ST_list = []

        for i in range(k):
            ab_off_i = ab_offsets[i]
            fABi = clinical_outputs[ab_off_i:ab_off_i + n, c_idx]

            # Mask NaN entries
            valid = ~(np.isnan(fA) | np.isnan(fB) | np.isnan(fABi))
            if np.sum(valid) < 10:
                S1_list.append(0.0)
                ST_list.append(0.0)
                continue

            fA_v = fA[valid]
            fB_v = fB[valid]
            fABi_v = fABi[valid]
            n_v = len(fA_v)

            # Saltelli S1 estimator
            S1 = np.mean(fB_v * (fABi_v - fA_v)) / var_total
            S1_list.append(float(np.clip(S1, -0.5, 1.5)))

            # Jansen ST estimator
            ST = np.mean((fA_v - fABi_v) ** 2) / (2 * var_total)
            ST_list.append(float(np.clip(ST, 0.0, 2.0)))

        results[c_name] = {
            'S1': S1_list,
            'ST': ST_list,
        }

    return results


def run_sobol(top_params, n_samples=1024, t_hours=168, n_workers=None, seed=42):
    """Full Sobol analysis on pre-selected parameters.

    Args:
        top_params: list of dicts from Morris analysis with name, index, lo, hi
        n_samples: base sample size
        t_hours: simulation duration
        n_workers: parallel workers
        seed: random seed

    Returns:
        dict with results, metadata
    """
    import time
    t0 = time.time()

    base_params, r_values = build_params()
    y0 = build_inits(r_values)
    k = len(top_params)

    print(f"[Sobol] {k} parameters, {n_samples} base samples")

    # Generate samples
    param_sets, sample_info = generate_sobol_samples(
        top_params, base_params, n_samples=n_samples, seed=seed)
    n_total = len(param_sets)
    print(f"[Sobol] {n_total} total simulations")

    # Run simulations
    print(f"[Sobol] Running {n_total} simulations ({t_hours}h each)...")
    batch_results = batch_integrate(y0, param_sets, t_end=t_hours,
                                     dt_output=max(1.0, t_hours / 10),
                                     n_workers=n_workers)

    # Extract clinical outputs
    n_clinical = len(CLINICAL_KEYS)
    clinical_outputs = np.full((n_total, n_clinical), np.nan)

    n_ok = 0
    for idx, success, final_state in batch_results:
        if success and final_state is not None:
            try:
                clinical_outputs[idx] = extract_clinical_from_state(final_state)
                n_ok += 1
            except Exception:
                pass

    print(f"[Sobol] {n_ok}/{n_total} simulations completed successfully")

    # Compute indices
    raw_results = compute_sobol_indices(clinical_outputs, sample_info)

    # Format results with param names
    results = {}
    for c_name, indices in raw_results.items():
        results[c_name] = []
        for i in range(k):
            results[c_name].append({
                'param': top_params[i]['name'],
                'param_index': top_params[i]['index'],
                'S1': indices['S1'][i] if i < len(indices['S1']) else 0.0,
                'ST': indices['ST'][i] if i < len(indices['ST']) else 0.0,
            })
        # Sort by ST descending
        results[c_name].sort(key=lambda x: x['ST'], reverse=True)

    elapsed = time.time() - t0

    return {
        'results': results,
        'metadata': {
            'k': k,
            'n_samples': n_samples,
            'n_total_sims': n_total,
            'n_ok': n_ok,
            't_hours': t_hours,
            'elapsed_seconds': round(elapsed, 1),
        }
    }
