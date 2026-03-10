"""
Slow-model integrator: drop-in replacement for hallow_c_driver.integrate().

Uses the QSS reduced model (48 slow ODEs instead of 70 stiff ODEs).
Allows much larger time steps because cardiac-cycle stiffness is eliminated.
"""

import numpy as np
from scipy.integrate import ode
import time
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from hallow_c_driver import N_STATE
from slow_model.qss_model import QSSModel, SLOW_INDICES, N_SLOW


def integrate_slow(y0_full, params, t_end, dt_output=1.0,
                   max_step=0.5, atol=1e-6, rtol=1e-4):
    """Integrate the slow model from initial conditions.

    Drop-in replacement for hallow_c_driver.integrate() using QSS model.
    Returns reconstructed full 70-dim states at each output point.

    Args:
        y0_full: full 70-dim initial state (from build_inits or committed segment)
        params: 430-dim parameter array
        t_end: duration in hours
        dt_output: output interval in hours (can be large, e.g. 1-24h)
        max_step: maximum integration step (0.5-1.0h recommended)
        atol: absolute tolerance
        rtol: relative tolerance

    Returns:
        (results, completed) in same format as hallow_c_driver.integrate()
        results = [(t, y_full_70dim), ...] where y_full is QSS-reconstructed
    """
    t0 = time.time()

    model = QSSModel(params)
    slow0 = model.extract_slow(y0_full)

    # Set up scipy LSODA integrator with relaxed settings
    solver = ode(model.rhs_scipy)
    solver.set_integrator('lsoda',
                          atol=atol, rtol=rtol,
                          max_step=max_step,
                          nsteps=100000)
    solver.set_initial_value(slow0, 0.0)

    # Initial full state
    full0 = model.slow_to_full(slow0)
    results = [(0.0, full0.copy())]

    t_points = np.arange(dt_output, t_end + dt_output / 2, dt_output)
    completed = True

    for t in t_points:
        solver.integrate(t)
        if not solver.successful():
            print(f"[slow_model] Integration failed at t={t:.4f}h")
            completed = False
            break

        slow_state = solver.y
        if np.any(np.isnan(slow_state)):
            nan_idx = np.where(np.isnan(slow_state))[0]
            print(f"[slow_model] NaN at t={t:.4f}h in slow indices: {nan_idx.tolist()}")
            completed = False
            break

        if np.any(np.isinf(slow_state)):
            inf_idx = np.where(np.isinf(slow_state))[0]
            print(f"[slow_model] Inf at t={t:.4f}h in slow indices: {inf_idx.tolist()}")
            completed = False
            break

        # Reconstruct full state for output
        full_y = model.slow_to_full(slow_state)
        results.append((t, full_y.copy()))

    elapsed = time.time() - t0
    status = 'completed' if completed else 'INCOMPLETE'
    print(f"[slow_model] {status}: {len(results)} points, "
          f"t=0..{results[-1][0]:.4f} of {t_end:.4f}h, {elapsed:.2f}s")

    return results, completed


def batch_integrate_slow(y0, param_sets, t_end=168.0, dt_output=1.0,
                         n_workers=None, progress=True):
    """Run many slow-model simulations in parallel.

    Same interface as hallow_c_driver.batch_integrate().

    Args:
        y0: initial state vector (shared across all sims)
        param_sets: list of parameter arrays
        t_end: simulation duration in hours
        dt_output: output interval
        n_workers: parallel workers
        progress: print progress

    Returns:
        list of (idx, success, final_state) sorted by index
    """
    import multiprocessing as mp

    if n_workers is None:
        n_workers = min(mp.cpu_count(), len(param_sets))

    def worker(args):
        idx, y0, params, t_end, dt_output = args
        try:
            results, completed = integrate_slow(y0, params, t_end, dt_output)
            if completed and len(results) > 0:
                return (idx, True, results[-1][1])
            return (idx, False, None)
        except Exception:
            return (idx, False, None)

    tasks = [(i, y0, p, t_end, dt_output) for i, p in enumerate(param_sets)]

    if progress:
        print(f"[slow_batch] {len(tasks)} sims, {n_workers} workers, "
              f"t_end={t_end}h")

    t0_batch = time.time()
    ctx = mp.get_context('fork')
    batch_results = []

    with ctx.Pool(n_workers) as pool:
        for result in pool.imap_unordered(worker, tasks, chunksize=1):
            batch_results.append(result)
            if progress and len(batch_results) % max(1, len(tasks) // 5) == 0:
                n_ok = sum(1 for r in batch_results if r[1])
                print(f"[slow_batch] {len(batch_results)}/{len(tasks)} "
                      f"({n_ok} ok) {time.time()-t0_batch:.1f}s")

    batch_results.sort(key=lambda x: x[0])
    elapsed = time.time() - t0_batch
    n_ok = sum(1 for r in batch_results if r[1])

    if progress:
        print(f"[slow_batch] Done: {n_ok}/{len(tasks)} ok in {elapsed:.1f}s")

    return batch_results
