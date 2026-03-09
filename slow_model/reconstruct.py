"""
Full-state reconstruction from slow-model state.

Given a slow-model state at time T, reconstructs the full 70-dim state
so the mechanistic C solver can resume without transient blow-up.

Two-phase approach:
  1. Algebraic expansion: QSS values for fast variables (instant)
  2. Settling integration: short full-model run to converge beat dynamics (~0.5h)
"""

import numpy as np
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from hallow_c_driver import integrate, N_STATE
from slow_model.qss_model import QSSModel


def algebraic_reconstruct(slow_state, params):
    """Phase 1: Algebraic reconstruction using QSS cardiac computation.

    Fast variables are set to their beat-averaged steady-state values.
    This produces a state that is close to the true manifold but may
    have small inconsistencies in the beat-level variables.

    Args:
        slow_state: N_SLOW-dim slow state vector
        params: 430-dim parameter array

    Returns:
        full 70-dim state vector
    """
    model = QSSModel(params)
    return model.slow_to_full(slow_state)


def full_reconstruct(slow_state, params, settle_hours=0.5, max_attempts=3):
    """Full reconstruction: algebraic + settling integration.

    Phase 1: Algebraic expansion (instant)
    Phase 2: Short full-model integration to let beat-level dynamics
             converge to their true limit cycle. The slow variables barely
             move in 0.5h, so they remain correct.

    Args:
        slow_state: N_SLOW-dim slow state vector
        params: 430-dim parameter array
        settle_hours: duration of settling integration (0.1-0.5h typical)
        max_attempts: number of attempts with increasing settle time

    Returns:
        (full_state_70dim, success)
    """
    y0 = algebraic_reconstruct(slow_state, params)

    for attempt in range(max_attempts):
        dt = settle_hours * (attempt + 1)
        try:
            results, completed = integrate(y0, params, dt, dt_output=dt)
        except Exception as e:
            print(f"[reconstruct] Settle attempt {attempt+1} failed: {e}")
            continue

        if completed and len(results) > 0:
            y_settled = results[-1][1]
            if not np.any(np.isnan(y_settled)) and not np.any(np.isinf(y_settled)):
                return y_settled, True
            print(f"[reconstruct] Settle attempt {attempt+1}: NaN/Inf in result")
        else:
            print(f"[reconstruct] Settle attempt {attempt+1}: integration incomplete")

    # Fall back to algebraic reconstruction if settling fails
    print("[reconstruct] All settle attempts failed, returning algebraic reconstruction")
    return y0, False


def reconstruct_from_full(y_full, params, settle_hours=0.5):
    """Reconstruct from a QSS-expanded full state (convenience wrapper).

    Use this when you already have a full 70-dim state from the slow model
    and want to settle it for the full C model.

    Args:
        y_full: 70-dim state (possibly from slow model output)
        params: 430-dim parameter array
        settle_hours: settling duration

    Returns:
        (settled_state, success)
    """
    model = QSSModel(params)
    slow = model.extract_slow(y_full)
    return full_reconstruct(slow, params, settle_hours)


def reconstruct_trajectory(slow_results, params, interval_hours=24,
                           settle_hours=0.1):
    """Reconstruct full states at regular intervals along a slow trajectory.

    For a 6-month slow simulation, reconstruct every 24h (180 points).
    Each reconstruction costs ~0.1-0.5h of C integration.

    Args:
        slow_results: list of (t, y_full) from integrate_slow()
        params: 430-dim parameter array
        interval_hours: interval between reconstructions
        settle_hours: settling time per reconstruction

    Returns:
        list of (t, settled_y_full, success) tuples
    """
    if not slow_results:
        return []

    model = QSSModel(params)
    reconstructed = []

    # Find points at regular intervals
    t_max = slow_results[-1][0]
    t_targets = np.arange(0, t_max + interval_hours / 2, interval_hours)

    result_idx = 0
    for t_target in t_targets:
        # Find closest result point
        while (result_idx < len(slow_results) - 1 and
               slow_results[result_idx + 1][0] <= t_target):
            result_idx += 1

        t, y_full = slow_results[result_idx]
        slow = model.extract_slow(y_full)
        y_settled, success = full_reconstruct(slow, params, settle_hours)
        reconstructed.append((t, y_settled, success))

    n_ok = sum(1 for _, _, s in reconstructed if s)
    print(f"[reconstruct] {n_ok}/{len(reconstructed)} points reconstructed "
          f"successfully")

    return reconstructed
