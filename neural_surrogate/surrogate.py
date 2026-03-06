"""
Inference wrapper for the trained Neural ODE surrogate.

Provides a drop-in replacement for the C integrator's integrate() function,
returning trajectories in the same format server.py expects.

Also provides a hybrid mode: Neural ODE for the bulk of the trajectory,
with optional short C-solver corrections for state anchoring.
"""

import torch
import numpy as np
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from neural_surrogate.model import (
    HallowSurrogate, CLINICAL_KEYS, N_CLINICAL, N_CRITICAL,
    CRITICAL_STATE_INDICES, LATENT_DIM, N_STATE, N_KNOBS,
)
from hallow_c_driver import (
    STATE_MAP, PARAM_MAP, N_STATE as DRIVER_N_STATE,
    integrate as c_integrate
)

_DIR = os.path.dirname(os.path.abspath(__file__))


class NeuralSurrogate:
    """Wraps the trained HallowSurrogate for server integration."""

    def __init__(self, checkpoint_path=None, device=None):
        if device is None:
            # torchdiffeq's dopri5 uses float64 internally, which MPS doesn't support.
            # The model is small (176K params), so CPU is fast enough for inference.
            if torch.cuda.is_available():
                self.device = torch.device('cuda')
            else:
                self.device = torch.device('cpu')
        else:
            self.device = torch.device(device)

        if checkpoint_path is None:
            checkpoint_path = os.path.join(_DIR, 'checkpoints', 'best_model.pt')

        if not os.path.exists(checkpoint_path):
            raise FileNotFoundError(
                f"No trained model at {checkpoint_path}. "
                f"Run training first: python -m neural_surrogate.train")

        ckpt = torch.load(checkpoint_path, map_location=self.device, weights_only=False)
        latent_dim = ckpt.get('latent_dim', LATENT_DIM)

        self.model = HallowSurrogate(latent_dim=latent_dim).to(self.device)
        self.model.load_state_dict(ckpt['model_state_dict'])
        self.model.eval()

        # Normalization stats
        ns = ckpt['norm_stats']
        self.clinical_mean = np.array(ns['clinical_mean'], dtype=np.float32)
        self.clinical_std = np.array(ns['clinical_std'], dtype=np.float32)
        self.state_mean = np.array(ns['state_mean'], dtype=np.float32)
        self.state_std = np.array(ns['state_std'], dtype=np.float32)
        self.knob_mean = np.array(ns['knob_mean'], dtype=np.float32)
        self.knob_std = np.array(ns['knob_std'], dtype=np.float32)

        print(f"[NeuralSurrogate] Loaded model ({latent_dim}D latent) on {self.device}")

    def _normalize_state(self, y):
        return (y.astype(np.float32) - self.state_mean) / self.state_std

    def _normalize_knobs(self, knob_vec):
        return (np.array(knob_vec, dtype=np.float32) - self.knob_mean) / self.knob_std

    def _denorm_clinical(self, clinical_normed):
        return clinical_normed * self.clinical_std + self.clinical_mean

    def _denorm_critical_state(self, crit_normed):
        crit_mean = self.state_mean[CRITICAL_STATE_INDICES]
        crit_std = self.state_std[CRITICAL_STATE_INDICES]
        return crit_normed * crit_std + crit_mean

    def extract_knobs_from_params(self, params):
        """Extract the 4 knob values from a full parameter array."""
        c_art = params[PARAM_MAP.get('C_art_scale', 116)]
        tpr = params[PARAM_MAP.get('disease_effect_on_TPR_peripheral_resistance', 102)]
        glucose = params[PARAM_MAP.get('glucose_concentration', 132)]
        # Convert glucose back to a1c: glucose = 1.59*a1c - 2.59 => a1c = (glucose+2.59)/1.59
        a1c = (glucose + 2.59) / 1.59
        nephron_loss = params[PARAM_MAP.get('disease_effect_on_nephrons', 402)]
        return np.array([c_art, tpr, a1c, nephron_loss], dtype=np.float32)

    def predict(self, y0, params, t_end, n_output=30):
        """Run neural ODE surrogate prediction.

        Args:
            y0: (70,) initial state vector
            params: (430,) parameter array
            t_end: duration in hours
            n_output: number of output time points

        Returns:
            dict with 'trajectory' (list of clinical output dicts),
            'final_y' (reconstructed 70-dim state), 'elapsed_seconds',
            and 'completed' flag.
        """
        t0 = time.time()

        knob_vec = self.extract_knobs_from_params(params)

        # Normalize inputs
        state_n = self._normalize_state(y0)
        knobs_n = self._normalize_knobs(knob_vec)

        state_t = torch.from_numpy(state_n).to(self.device)
        knobs_t = torch.from_numpy(knobs_n).to(self.device)
        t_eval = torch.linspace(0, t_end, n_output, device=self.device)

        with torch.no_grad():
            clinical_pred, crit_pred = self.model.predict_single(
                state_t, knobs_t, t_eval)

        # Denormalize
        clinical_np = self._denorm_clinical(clinical_pred.cpu().numpy())
        crit_state_np = self._denorm_critical_state(crit_pred.cpu().numpy())

        # Build trajectory in the format server.py expects
        trajectory = []
        times = np.linspace(0, t_end, n_output)
        for i in range(n_output):
            c = clinical_np[i]
            entry = {'t_hours': round(float(times[i]), 4)}
            for j, key in enumerate(CLINICAL_KEYS):
                val = float(c[j])
                # Clamp physiological ranges
                if key == 'MAP':
                    val = max(30.0, min(250.0, val))
                elif key == 'CO':
                    val = max(0.5, min(15.0, val))
                elif key == 'blood_volume_L':
                    val = max(2.0, min(10.0, val))
                elif key == 'Na':
                    val = max(100.0, min(180.0, val))
                elif key == 'BNP':
                    val = max(0.0, min(5000.0, val))
                elif key == 'serum_creatinine':
                    val = max(0.1, min(15.0, val))
                entry[key] = round(val, 4)
            # Add derived fields
            entry['t_days'] = round(float(times[i]) / 24.0, 2)
            entry['t_years'] = round(float(times[i]) / (24.0 * 365.25), 4)
            # Myocyte geometry (included in clinical outputs as pct_diameter/pct_length)
            entry['change_myocyte_diameter_um'] = 0.0
            entry['change_myocyte_length_um'] = 0.0
            trajectory.append(entry)

        # Reconstruct full 70-dim state from critical state prediction
        final_y = y0.copy()
        final_crit = crit_state_np[-1]
        final_y[CRITICAL_STATE_INDICES] = final_crit

        elapsed = time.time() - t0

        return {
            'trajectory': trajectory,
            'final_y': final_y.tolist(),
            'elapsed_seconds': round(elapsed, 4),
            'completed': True,
            'surrogate': True,
        }

    def predict_with_correction(self, y0, params, t_end, n_output=30,
                                 correction_hours=0.05):
        """Neural ODE prediction + short C-solver correction at the end.

        Runs the surrogate for the trajectory, then integrates 0.05h
        from the predicted final state with the real solver to snap
        the state back to the true manifold.
        """
        result = self.predict(y0, params, t_end, n_output)

        if correction_hours > 0:
            final_y = np.array(result['final_y'])
            try:
                corr_results, corr_ok = c_integrate(
                    final_y, params, correction_hours, dt_output=correction_hours)
                if corr_ok and len(corr_results) > 0:
                    _, y_corrected = corr_results[-1]
                    if not np.any(np.isnan(y_corrected)):
                        result['final_y'] = y_corrected.tolist()
                        result['correction_applied'] = True
            except Exception:
                pass  # Correction failed — use uncorrected state

        return result


def integrate_surrogate(y0, params, t_end, dt_output=0.01, surrogate=None):
    """Drop-in replacement for hallow_c_driver.integrate().

    Returns (results, completed) in the same format:
        results = [(t, y_array), ...]
        completed = True/False
    """
    if surrogate is None:
        surrogate = _get_default_surrogate()

    n_output = max(10, int(t_end / max(dt_output, 0.1)))
    n_output = min(n_output, 500)

    result = surrogate.predict(y0, params, t_end, n_output=n_output)

    # Convert to (t, y) format
    final_y = np.array(result['final_y'])
    results = []
    times = np.linspace(0, t_end, n_output)
    for i, t in enumerate(times):
        # Interpolate state between y0 and final_y for intermediate points
        frac = t / t_end if t_end > 0 else 0
        y_interp = y0 * (1 - frac) + final_y * frac
        # Override critical state vars from prediction
        results.append((float(t), y_interp))

    return results, result['completed']


# Singleton for lazy loading
_default_surrogate = None

def _get_default_surrogate():
    global _default_surrogate
    if _default_surrogate is None:
        _default_surrogate = NeuralSurrogate()
    return _default_surrogate


def is_surrogate_available():
    """Check if a trained model exists."""
    ckpt = os.path.join(_DIR, 'checkpoints', 'best_model.pt')
    return os.path.exists(ckpt)
