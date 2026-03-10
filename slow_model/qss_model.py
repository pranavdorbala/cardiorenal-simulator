"""
Quasi-steady-state reduced cardiorenal model.

Eliminates beat-level cardiac dynamics by replacing 22 fast variables
with algebraic QSS values. Keeps 48 slow variables as ODEs.

The key trick: calls the *same* C hallow_rhs() on a reconstructed full
state, then extracts only the slow derivatives. This ensures exact fidelity
with the original model for all renal equations.

Expected speedup: 50-100x for multi-day simulations, because max_step
can be 0.5-1.0h instead of 0.01h.
"""

import numpy as np
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from hallow_c_driver import HallowRHS, N_STATE, N_PARAM
from slow_model.cardiac_algebraic import compute_qss_fast_variables

# Fast variables: replaced by algebraic QSS at each RHS evaluation
# These are the beat-level cardiac volumes, flows, pressures, and their
# first-order delay filters (time constant C_cycle2 = 100, ~0.01h)
FAST_INDICES = [
    0,   # venous_volume (driven by blood volume conservation)
    1,   # LV_volume (oscillates within beat)
    2,   # arterial_volume
    3,   # peripheral_circulation_volume
    4,   # RV_volume
    5,   # pulmonary_arterial_volume
    6,   # pulmonary_venous_volume
    7,   # aortic_blood_flow_delayed (C_cycle2 filter)
    8,   # pulmonary_blood_flow_delayed (C_cycle2 filter)
    13,  # LV_sarcomere_length_delayed (C_cycle filter)
    14,  # RV_sarcomere_length_delayed (C_cycle filter)
    15,  # LV_EDV (C_cycle2 filter)
    16,  # LV_EDP (C_cycle2 filter)
    17,  # LV_EDS (C_cycle2 filter)
    18,  # arterial_pressure_delayed (C_cycle2 filter)
    19,  # arterial_pressure_bigger_delay (C_cycle2 filter)
    20,  # systolic_pressure (C_cycle2 filter)
    21,  # diastolic_pressure (C_cycle2 filter)
    22,  # venous_pressure_delayed (C_cycle2 filter)
    23,  # venous_pressure_bigger_delay (C_cycle2 filter)
    24,  # systolic_venous_pressure (C_cycle2 filter)
    25,  # diastolic_venous_pressure (C_cycle2 filter)
]
FAST_SET = set(FAST_INDICES)

# Slow variables: kept as ODEs
SLOW_INDICES = [i for i in range(N_STATE) if i not in FAST_SET]
N_SLOW = len(SLOW_INDICES)

# Reverse maps for fast conversion
_SLOW_POS = {full_idx: pos for pos, full_idx in enumerate(SLOW_INDICES)}
_FULL_FROM_SLOW = np.array(SLOW_INDICES, dtype=int)


class QSSModel:
    """Quasi-steady-state reduced cardiorenal model.

    Eliminates beat-level cardiac dynamics. All fast variables are
    computed algebraically from the slow state at each RHS evaluation.
    Allows dt = 0.5-1.0 hours (vs 0.01 hours for full model).
    """

    def __init__(self, params):
        self.params = params.copy()
        self._rhs = HallowRHS(self.params)
        self._full_y = np.zeros(N_STATE)

    def extract_slow(self, full_state):
        """Extract the N_SLOW-dim slow state from a full 70-dim state."""
        return full_state[_FULL_FROM_SLOW].copy()

    def slow_to_full(self, slow_state):
        """Expand slow state to full 70-dim state using algebraic QSS.

        1. Place slow variables at their indices
        2. Compute beat-average hemodynamics
        3. Fill in fast variables with QSS values

        Returns: full 70-dim state vector
        """
        y = self._full_y
        y[:] = 0.0

        # Place slow variables
        for pos, full_idx in enumerate(SLOW_INDICES):
            y[full_idx] = slow_state[pos]

        # Compute QSS fast variables
        qss = compute_qss_fast_variables(y, self.params)
        for idx, val in qss.items():
            y[idx] = val

        return y

    def rhs(self, t, slow_state):
        """Compute d(slow)/dt.

        1. Expand to full state via slow_to_full()
        2. Call C hallow_rhs() to get all 70 derivatives
        3. Extract only the N_SLOW slow derivatives

        Returns: N_SLOW-dim derivative vector
        """
        full_y = self.slow_to_full(slow_state)
        dydt_full = self._rhs(t, full_y)

        # Extract slow derivatives
        dslow = np.empty(N_SLOW)
        for pos, full_idx in enumerate(SLOW_INDICES):
            dslow[pos] = dydt_full[full_idx]

        return dslow

    def rhs_scipy(self, t, slow_state):
        """Scipy-compatible RHS (same signature as rhs but handles edge cases)."""
        dslow = self.rhs(t, slow_state)
        # Clamp any NaN/Inf derivatives to zero for solver stability
        dslow[~np.isfinite(dslow)] = 0.0
        return dslow
