"""
Neural ODE surrogate for the Hallow cardiorenal model.

Architecture:
  Encoder:    (70 state vars + 4 knobs) -> latent z0 in R^d
  LatentODE:  dz/dt = f_theta(z, knobs), solved with adaptive solver at large dt
  Decoder:    z -> 15 clinical outputs

The latent ODE operates on the slow manifold, bypassing the stiff
cardiac-cycle dynamics that force the mechanistic solver to take tiny steps.
"""

import torch
import torch.nn as nn
from torchdiffeq import odeint

N_STATE = 70
N_KNOBS = 4
LATENT_DIM = 24

# Clinical outputs the surrogate must predict — order matters for loss weighting
CLINICAL_KEYS = [
    'MAP', 'SBP', 'DBP', 'CO',
    'blood_volume_L', 'Na', 'LV_EDV_mL', 'LV_EDP_mmHg',
    'LV_EDS', 'LV_active_stress_peak', 'LV_mass', 'BNP',
    'serum_creatinine', 'pct_diameter', 'pct_length',
]
N_CLINICAL = len(CLINICAL_KEYS)

# Per-output importance weights for loss function.
# Higher weight = tighter fit required. Hemodynamics and remodeling
# markers get the strongest weighting.
CLINICAL_WEIGHTS = {
    'MAP': 5.0, 'SBP': 3.0, 'DBP': 3.0, 'CO': 5.0,
    'blood_volume_L': 4.0, 'Na': 3.0, 'LV_EDV_mL': 3.0, 'LV_EDP_mmHg': 4.0,
    'LV_EDS': 2.0, 'LV_active_stress_peak': 2.0, 'LV_mass': 4.0, 'BNP': 3.0,
    'serum_creatinine': 3.0, 'pct_diameter': 4.0, 'pct_length': 4.0,
}

# State variable indices used for clinical output extraction (from state_map)
# These are the indices into the 70-dim state vector that the decoder
# needs to reconstruct to recover full state for continuation.
CRITICAL_STATE_INDICES = [
    0, 1, 2, 3, 4, 5, 6,       # volumes (venous, LV, arterial, peripheral, RV, pulm_art, pulm_ven)
    9, 10, 11,                   # myocyte length/diameter changes, active stress peak
    15, 16, 17,                  # LV_EDV, LV_EDP, LV_EDS
    20, 21,                      # systolic/diastolic pressure
    26, 27,                      # CO, CO_delayed
    28, 29, 30, 31, 32,         # RAAS: AngI, AngII, AT1, AT2, PRC
    33, 34,                      # blood_volume_L, interstitial_fluid_volume
    35, 36, 37,                  # sodium_amount, IF_sodium, stored_sodium
    38, 39,                      # TGF, aldosterone
    44,                          # vasopressin
    48,                          # renal_blood_flow_L_min_delayed
    60,                          # serum_creatinine
    69,                          # mitral_valve_leak
]
N_CRITICAL = len(CRITICAL_STATE_INDICES)


class Encoder(nn.Module):
    """Maps full state + knobs to latent initial condition."""

    def __init__(self, state_dim=N_STATE, knob_dim=N_KNOBS, latent_dim=LATENT_DIM):
        super().__init__()
        in_dim = state_dim + knob_dim
        self.net = nn.Sequential(
            nn.Linear(in_dim, 256),
            nn.LayerNorm(256),
            nn.GELU(),
            nn.Linear(256, 128),
            nn.LayerNorm(128),
            nn.GELU(),
            nn.Linear(128, latent_dim),
        )

    def forward(self, state, knobs):
        x = torch.cat([state, knobs], dim=-1)
        return self.net(x)


class LatentODEFunc(nn.Module):
    """Defines dz/dt = f(z, knobs). Knobs are concatenated at every eval."""

    def __init__(self, latent_dim=LATENT_DIM, knob_dim=N_KNOBS, hidden=192):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(latent_dim + knob_dim, hidden),
            nn.GELU(),
            nn.Linear(hidden, hidden),
            nn.GELU(),
            nn.Linear(hidden, hidden),
            nn.GELU(),
            nn.Linear(hidden, latent_dim),
        )
        # Store knobs for use during odeint (set before each solve)
        self._knobs = None

    def set_knobs(self, knobs):
        self._knobs = knobs

    def forward(self, t, z):
        # z: (batch, latent) or (latent,)
        if self._knobs is None:
            raise RuntimeError("Call set_knobs() before integrating")
        knobs = self._knobs
        if z.dim() == 1:
            inp = torch.cat([z, knobs])
        else:
            if knobs.dim() == 1:
                knobs = knobs.unsqueeze(0).expand(z.shape[0], -1)
            inp = torch.cat([z, knobs], dim=-1)
        return self.net(inp)


class ClinicalDecoder(nn.Module):
    """Maps latent z to clinical outputs + critical state variables."""

    def __init__(self, latent_dim=LATENT_DIM, n_clinical=N_CLINICAL,
                 n_critical_state=N_CRITICAL):
        super().__init__()
        self.n_clinical = n_clinical
        self.n_critical_state = n_critical_state
        out_dim = n_clinical + n_critical_state
        self.net = nn.Sequential(
            nn.Linear(latent_dim, 192),
            nn.GELU(),
            nn.Linear(192, 128),
            nn.GELU(),
            nn.Linear(128, out_dim),
        )

    def forward(self, z):
        out = self.net(z)
        clinical = out[..., :self.n_clinical]
        state = out[..., self.n_clinical:]
        return clinical, state


class HallowSurrogate(nn.Module):
    """Full surrogate: encode → ODE solve → decode."""

    def __init__(self, latent_dim=LATENT_DIM):
        super().__init__()
        self.encoder = Encoder(latent_dim=latent_dim)
        self.ode_func = LatentODEFunc(latent_dim=latent_dim)
        self.decoder = ClinicalDecoder(latent_dim=latent_dim)
        self.latent_dim = latent_dim

    def forward(self, state, knobs, t_eval, method='dopri5'):
        """
        Args:
            state: (batch, 70) or (70,) — full state vector (normalized)
            knobs: (batch, 4) or (4,) — knob values (normalized)
            t_eval: (T,) — output times (normalized hours)
            method: ODE solver method

        Returns:
            clinical: (T, batch, N_CLINICAL) or (T, N_CLINICAL)
            critical_state: (T, batch, N_CRITICAL) or (T, N_CRITICAL)
        """
        z0 = self.encoder(state, knobs)
        self.ode_func.set_knobs(knobs)
        z_traj = odeint(
            self.ode_func, z0, t_eval,
            method=method,
            rtol=1e-3, atol=1e-3,
        )
        clinical, critical_state = self.decoder(z_traj)
        return clinical, critical_state

    def predict_single(self, state, knobs, t_eval):
        """Convenience for single-sample inference (no batch dim)."""
        with torch.no_grad():
            clinical, crit = self.forward(state, knobs, t_eval)
        return clinical, crit
