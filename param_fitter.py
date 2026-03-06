"""
param_fitter.py — Disease parameter optimization for remodeling targets.

Uses scipy.optimize.differential_evolution to search the disease parameter
space and find trajectories that drive myocyte geometry through user-specified
target points.
"""

import numpy as np
import time
import threading
from scipy.optimize import differential_evolution

from hallow_c_driver import (
    build_params, build_inits, integrate,
    STATE_MAP, PARAM_MAP, N_STATE, N_PARAM
)


# Action space bounds
MECHANISMS = {
    'C_art_scale':  {'normal': 1.0, 'bounds': (0.3, 1.0)},
    'TPR_mult':     {'normal': 1.0, 'bounds': (0.8, 2.0)},
    'a1c':          {'normal': 5.7, 'bounds': (4.0, 12.0)},
    'nephron_loss': {'normal': 0.0, 'bounds': (0.0, 0.8)},
}
MECH_KEYS = list(MECHANISMS.keys())


def _extract_myocyte(y):
    """Extract pct_diameter and pct_length from a state vector."""
    cmd = y[STATE_MAP['change_in_myocyte_diameter']]
    cml = y[STATE_MAP['change_in_myocyte_length']]
    Pi = 3.1416
    V_w_0 = 0.00012
    btmv = V_w_0 - V_w_0*0.02 - V_w_0*0.02 - V_w_0*0.22
    bsmv = btmv / 3.3e9
    bmd = 2 * np.sqrt(bsmv / (Pi * 0.000115))
    pct_d = 100 * cmd / bmd if bmd > 0 else 0
    pct_l = 100 * cml / 0.000115
    return pct_d, pct_l


def _apply_knobs(params, base_params, knobs):
    """Apply knob values to parameter array (mirrors server.SimState.set_knobs)."""
    params[:] = base_params
    if 'C_art_scale' in knobs:
        idx = PARAM_MAP.get('C_art_scale')
        if idx is not None:
            params[idx] = float(knobs['C_art_scale'])
    if 'TPR_mult' in knobs:
        idx = PARAM_MAP.get('disease_effect_on_TPR_peripheral_resistance')
        if idx is not None:
            params[idx] = float(knobs['TPR_mult'])
    if 'a1c' in knobs:
        a1c = float(knobs['a1c'])
        glucose = max(3.0, 1.59 * a1c - 2.59)
        idx = PARAM_MAP.get('glucose_concentration')
        if idx is not None:
            params[idx] = glucose
    if 'nephron_loss' in knobs:
        idx = PARAM_MAP.get('disease_effect_on_nephrons')
        if idx is not None:
            params[idx] = float(knobs['nephron_loss'])


class ParamFitter:
    def __init__(self, targets, segment_days=30):
        """
        targets: list of 3 dicts [{day, pct_diameter, pct_length}, ...]
        """
        self.targets = sorted(targets, key=lambda t: t['day'])
        self.day_mid = self.targets[1]['day']
        self.day_final = self.targets[2]['day']
        self.segment_days = segment_days
        self.segment_hours = segment_days * 24
        self.n_segments = max(1, int(np.ceil(self.day_final / segment_days)))

        # Progress tracking
        self.progress = {
            'status': 'idle',
            'iteration': 0,
            'best_loss': float('inf'),
            'current_best_path': [],
            'params_used': None,
            'elapsed_seconds': 0,
        }
        self.lock = threading.Lock()
        self._stop = False

    def interpolate_params(self, x, day):
        """Piecewise linear interpolation: normal -> x_mid -> x_final."""
        result = {}
        for i, key in enumerate(MECH_KEYS):
            normal = MECHANISMS[key]['normal']
            v_mid = x[2 * i]
            v_final = x[2 * i + 1]
            if day <= 0:
                result[key] = normal
            elif day <= self.day_mid:
                frac = day / self.day_mid if self.day_mid > 0 else 1.0
                result[key] = normal + frac * (v_mid - normal)
            elif day <= self.day_final:
                span = self.day_final - self.day_mid
                frac = (day - self.day_mid) / span if span > 0 else 1.0
                result[key] = v_mid + frac * (v_final - v_mid)
            else:
                result[key] = v_final
        return result

    def evaluate(self, x):
        """Run model with parameter trajectory defined by x, return scalar loss."""
        if self._stop:
            return 1e6

        try:
            # Fresh simulation state
            params, r_values = build_params()
            base_params = params.copy()
            y = build_inits(r_values)
            t_hours = 0.0

            path = [{'x': 0.0, 'y': 0.0, 't_days': 0}]
            current_day = 0

            for seg_idx in range(self.n_segments):
                if self._stop:
                    return 1e6

                seg_end_day = min((seg_idx + 1) * self.segment_days, self.day_final)
                seg_mid_day = (current_day + seg_end_day) / 2
                knobs = self.interpolate_params(x, seg_mid_day)
                _apply_knobs(params, base_params, knobs)

                dt_hours = (seg_end_day - current_day) * 24
                if dt_hours <= 0:
                    break

                # Choose output spacing for speed
                dt_out = max(0.5, dt_hours / 5)
                results, completed = integrate(y, params, dt_hours, dt_output=dt_out)

                if len(results) < 2:
                    return 1e6

                # Advance state
                y = results[-1][1].copy()
                t_hours += results[-1][0]
                current_day = seg_end_day

                pct_d, pct_l = _extract_myocyte(y)
                path.append({
                    'x': round(pct_d, 3),
                    'y': round(pct_l, 3),
                    't_days': current_day,
                })

            # Compute loss
            loss = 0.0
            for target in self.targets[1:]:
                closest = min(path, key=lambda p: abs(p['t_days'] - target['day']))
                loss += (closest['x'] - target['pct_diameter']) ** 2
                loss += (closest['y'] - target['pct_length']) ** 2

            # Regularization
            for i, key in enumerate(MECH_KEYS):
                normal = MECHANISMS[key]['normal']
                loss += 0.001 * ((x[2*i] - normal)**2 + (x[2*i+1] - normal)**2)

            # Track best
            with self.lock:
                if loss < self.progress['best_loss']:
                    self.progress['best_loss'] = float(loss)
                    self.progress['current_best_path'] = path
                    self.progress['params_used'] = {
                        key: {'mid': float(x[2*i]), 'final': float(x[2*i+1])}
                        for i, key in enumerate(MECH_KEYS)
                    }

            return loss

        except Exception as e:
            print(f"[ParamFitter] evaluate error: {e}")
            return 1e6

    def fit(self, max_iter=12, popsize=5, callback=None):
        """Run differential_evolution. Returns result dict."""
        self._stop = False
        t0 = time.time()

        with self.lock:
            self.progress['status'] = 'running'
            self.progress['iteration'] = 0
            self.progress['best_loss'] = float('inf')

        bounds = []
        for key in MECH_KEYS:
            b = MECHANISMS[key]['bounds']
            bounds.append(b)  # mid
            bounds.append(b)  # final

        iter_count = [0]

        def de_callback(xk, convergence):
            iter_count[0] += 1
            with self.lock:
                self.progress['iteration'] = iter_count[0]
                self.progress['elapsed_seconds'] = round(time.time() - t0, 1)
            if callback:
                callback(self.progress)
            with self.lock:
                if self.progress['best_loss'] < 0.1:
                    return True
            if self._stop:
                return True
            return False

        result = differential_evolution(
            self.evaluate,
            bounds=bounds,
            maxiter=max_iter,
            popsize=popsize,
            tol=0.01,
            polish=False,
            callback=de_callback,
            seed=42,
        )

        elapsed = round(time.time() - t0, 1)

        with self.lock:
            self.progress['status'] = 'completed'
            self.progress['elapsed_seconds'] = elapsed
            self.progress['iteration'] = iter_count[0]

            return {
                'path': self.progress['current_best_path'],
                'params_used': self.progress['params_used'],
                'loss': float(result.fun),
                'iterations': iter_count[0],
                'elapsed_seconds': elapsed,
            }

    def stop(self):
        self._stop = True
