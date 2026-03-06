"""
server.py — Interactive Cardiorenal HFpEF Simulator

Runs the exact Hallow et al. model (compiled to C) with interactive knobs.
The dashboard sends knob values, the server integrates the ODE system,
and returns trajectory data + mechanistic messages.

Segment-based interaction:
  - /api/preview: integrate forward without committing (dashed line on graph)
  - /api/commit:  integrate and store permanently (solid line on graph)
  - /api/reset:   return to initial conditions

Knobs:
  - C_art_scale:  arterial compliance (1.0=normal, <1=stiffer → aging/HTN)
  - TPR_mult:     peripheral resistance multiplier (1.0=normal, >1=HTN)
  - glucose:      plasma glucose in mmol/L (5.5=normal, set by A1C)
  - a1c:          HbA1c % (converts to glucose + nephron loss rate)

Start: python3 server.py
Open:  http://localhost:5010
"""

from flask import Flask, jsonify, request, send_from_directory, Response
import numpy as np
import json
import math
import time
import os
import threading

from hallow_c_driver import (
    build_params, build_inits, HallowRHS, integrate,
    STATE_MAP, PARAM_MAP, N_STATE, N_PARAM
)
from param_fitter import ParamFitter

app = Flask(__name__, static_folder='static')

def sanitize(obj):
    """Replace NaN/Inf with None so JSON serialization doesn't break."""
    if isinstance(obj, float) and (math.isnan(obj) or math.isinf(obj)):
        return None
    if isinstance(obj, dict):
        return {k: sanitize(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [sanitize(v) for v in obj]
    return obj

def safe_jsonify(obj):
    """jsonify with NaN/Inf replaced by null."""
    return Response(json.dumps(sanitize(obj)), mimetype='application/json')

# =========================================================================
# SIMULATOR STATE
# =========================================================================
class SimState:
    def __init__(self):
        self.lock = threading.Lock()
        self.reset()
    
    def reset(self):
        self.params, self.r_values = build_params()
        self.y = build_inits(self.r_values)
        self.t_hours = 0.0
        self.segments = []  # stored segments
        self.base_params = self.params.copy()  # for resetting knobs
    
    def set_knobs(self, knobs):
        """Apply knob values to the parameter array."""
        self.params = self.base_params.copy()
        
        if 'C_art_scale' in knobs and knobs['C_art_scale'] is not None:
            idx = PARAM_MAP.get('C_art_scale')
            if idx is not None:
                self.params[idx] = float(knobs['C_art_scale'])
        
        if 'TPR_mult' in knobs and knobs['TPR_mult'] is not None:
            idx = PARAM_MAP.get('disease_effect_on_TPR_peripheral_resistance')
            if idx is not None:
                self.params[idx] = float(knobs['TPR_mult'])
        
        if 'glucose' in knobs and knobs['glucose'] is not None:
            idx = PARAM_MAP.get('glucose_concentration')
            if idx is not None:
                self.params[idx] = float(knobs['glucose'])
        
        if 'a1c' in knobs and knobs['a1c'] is not None:
            a1c = float(knobs['a1c'])
            # ADAG: glucose_mmol = 1.59 * A1C - 2.59
            glucose = max(3.0, 1.59 * a1c - 2.59)
            idx = PARAM_MAP.get('glucose_concentration')
            if idx is not None:
                self.params[idx] = glucose

        if 'nephron_loss' in knobs and knobs['nephron_loss'] is not None:
            loss_frac = float(knobs['nephron_loss'])
            idx = PARAM_MAP.get('disease_effect_on_nephrons')
            if idx is not None:
                self.params[idx] = loss_frac
    
    def run_segment(self, dt_hours, n_output=30):
        """Integrate the model forward by dt_hours, return trajectory."""
        t0 = time.time()

        print(f"[DEBUG] run_segment: dt_hours={dt_hours}, t_hours={self.t_hours}")
        print(f"[DEBUG] y0 has {np.sum(np.isnan(self.y))} NaN, {np.sum(np.isinf(self.y))} Inf values")
        print(f"[DEBUG] params has {np.sum(np.isnan(self.params))} NaN, {np.sum(np.isinf(self.params))} Inf values")
        if np.any(np.isnan(self.y)):
            nan_idx = np.where(np.isnan(self.y))[0]
            rev_state = {v: k for k, v in STATE_MAP.items()}
            nan_names = [rev_state.get(i, f"idx_{i}") for i in nan_idx]
            print(f"[DEBUG] NaN in y0 at: {nan_names}")
        if np.any(np.isnan(self.params)):
            nan_idx = np.where(np.isnan(self.params))[0]
            rev_param = {v: k for k, v in PARAM_MAP.items()}
            nan_names = [rev_param.get(i, f"idx_{i}") for i in nan_idx]
            print(f"[DEBUG] NaN in params at: {nan_names}")

        # Choose output spacing: denser for short segments
        if dt_hours <= 1:
            dt_out = 0.01
        elif dt_hours <= 24:
            dt_out = 0.1
        else:
            dt_out = max(0.5, dt_hours / 200)

        results, completed = integrate(self.y, self.params, dt_hours, dt_output=dt_out)
        elapsed = time.time() - t0

        final_y = results[-1][1]
        n_nan = np.sum(np.isnan(final_y))
        if n_nan > 0:
            rev_state = {v: k for k, v in STATE_MAP.items()}
            nan_idx = np.where(np.isnan(final_y))[0]
            nan_names = [rev_state.get(i, f"idx_{i}") for i in nan_idx]
            print(f"[DEBUG] WARNING: final state has {n_nan} NaN values: {nan_names}")

        if len(results) < 2:
            return None

        # Sample n_output evenly spaced points
        indices = np.linspace(0, len(results) - 1, min(n_output, len(results)), dtype=int)
        indices = sorted(set(indices))

        trajectory = []
        for i in indices:
            t, y = results[i]
            trajectory.append(extract_outputs(y, self.t_hours + t))

        # Mechanistic messages
        h0 = trajectory[0]
        h1 = trajectory[-1]

        actual_hours = results[-1][0]
        return {
            'trajectory': trajectory,
            'kidney_to_heart': kidney_to_heart_msg(h0, h1),
            'heart_to_kidney': heart_to_kidney_msg(h0, h1),
            'summary': chain_summary(h0, h1, actual_hours),
            'elapsed_seconds': round(elapsed, 2),
            'n_ode_steps': len(results),
            'final_y': results[-1][1].tolist(),
            'completed': completed,
        }
    
    def commit(self, result):
        """Store a segment result and advance the simulation clock."""
        if result is None:
            return
        dt_hours = result['trajectory'][-1]['t_hours'] - self.t_hours
        self.y = np.array(result['final_y'])
        self.t_hours = result['trajectory'][-1]['t_hours']
        self.segments.append({
            'trajectory': result['trajectory'],
            'messages': {
                'kidney_to_heart': result['kidney_to_heart'],
                'heart_to_kidney': result['heart_to_kidney'],
                'summary': result['summary'],
            },
            'knobs': {
                'C_art_scale': float(self.params[PARAM_MAP.get('C_art_scale', 0)]),
                'TPR_mult': float(self.params[PARAM_MAP.get('disease_effect_on_TPR_peripheral_resistance', 0)]),
                'glucose': float(self.params[PARAM_MAP.get('glucose_concentration', 0)]),
            },
        })

sim = SimState()

# =========================================================================
# PARAMETER FITTER (background jobs)
# =========================================================================
fit_jobs = {}  # fit_id -> {fitter, thread, result}
_fit_counter = [0]

def _run_fit(fit_id):
    """Background thread entry point for parameter fitting."""
    job = fit_jobs[fit_id]
    fitter = job['fitter']
    try:
        result = fitter.fit(
            max_iter=job.get('max_iter', 12),
            popsize=job.get('popsize', 5),
        )
        job['result'] = result
    except Exception as e:
        job['result'] = {'error': str(e)}
        with fitter.lock:
            fitter.progress['status'] = 'error'

# =========================================================================
# OUTPUT EXTRACTION
# =========================================================================
def extract_outputs(y, t_hours):
    """Extract clinically meaningful outputs from state vector."""
    SP = y[STATE_MAP['systolic_pressure']]
    DP = y[STATE_MAP['diastolic_pressure']]
    MAP = (SP / 3 + DP * 2 / 3) * 0.0075
    SBP = SP * 0.0075
    DBP = DP * 0.0075
    CO = y[STATE_MAP['CO_delayed']]
    BV = y[STATE_MAP['blood_volume_L']]
    Na = y[STATE_MAP['sodium_amount']] / BV if BV > 0 else 140
    EDV = y[STATE_MAP['LV_EDV']] * 1e6
    EDP = y[STATE_MAP['LV_EDP']] * 0.0075
    EDS = y[STATE_MAP['LV_EDS']]
    ASP = y[STATE_MAP['LV_active_stress_peak']]
    cmd = y[STATE_MAP['change_in_myocyte_diameter']]
    cml = y[STATE_MAP['change_in_myocyte_length']]
    sCr = y[STATE_MAP['serum_creatinine']] / BV if BV > 0 else 0.92
    
    # LV mass from myocyte geometry (same formula as R model)
    Pi = 3.1416
    V_w_0 = 0.00012
    btmv = V_w_0 - V_w_0*0.02 - V_w_0*0.02 - V_w_0*0.22
    bsmv = btmv / 3.3e9
    bmd = 2 * np.sqrt(bsmv / (Pi * 0.000115))
    ml = 0.000115 + cml
    md = bmd + cmd
    smv = ml * Pi * (md**2) / 4
    tmv = smv * 3.3e9
    tnmv = V_w_0 * 0.02 + V_w_0 * 0.22 + V_w_0 * 0.02
    LV_mass = 1e6 * (tmv + tnmv) * 1.05
    
    # BNP
    BNP = np.exp(0.0008 * ((EDS + 1736) / 5.094) + 3.14)
    BNP = min(BNP, 5000)
    
    # Percentage changes
    pct_d = 100 * cmd / bmd if bmd > 0 else 0
    pct_l = 100 * cml / 0.000115 if 0.000115 > 0 else 0
    
    return {
        't_hours': round(t_hours, 4),
        't_days': round(t_hours / 24, 2),
        't_years': round(t_hours / (24 * 365.25), 4),
        'MAP': round(MAP, 1),
        'SBP': round(SBP, 1),
        'DBP': round(DBP, 1),
        'CO': round(CO, 2),
        'blood_volume_L': round(BV, 3),
        'Na': round(Na, 1),
        'LV_EDV_mL': round(EDV, 1),
        'LV_EDP_mmHg': round(EDP, 2),
        'LV_EDS': round(EDS, 1),
        'LV_active_stress_peak': round(ASP, 0),
        'LV_mass': round(LV_mass, 1),
        'BNP': round(BNP, 1),
        'serum_creatinine': round(sCr, 3),
        'pct_diameter': round(pct_d, 3),
        'pct_length': round(pct_l, 3),
        'change_myocyte_diameter_um': round(cmd * 1e6, 3),
        'change_myocyte_length_um': round(cml * 1e6, 3),
    }

# =========================================================================
# MECHANISTIC MESSAGES
# =========================================================================
def kidney_to_heart_msg(h0, h1):
    msgs = []
    dbv = h1['blood_volume_L'] - h0['blood_volume_L']
    if abs(dbv) > 0.005:
        d = "expanding" if dbv > 0 else "contracting"
        msgs.append(f"Blood volume {d} ({h0['blood_volume_L']:.2f} → {h1['blood_volume_L']:.2f} L) — {'increased' if dbv > 0 else 'decreased'} venous return")
    dna = h1['Na'] - h0['Na']
    if abs(dna) > 0.5:
        msgs.append(f"Sodium concentration {'rising' if dna > 0 else 'falling'} ({h0['Na']:.0f} → {h1['Na']:.0f} mEq/L)")
    if not msgs:
        msgs.append("Kidney function stable — no significant signals to heart")
    return msgs

def heart_to_kidney_msg(h0, h1):
    msgs = []
    dmap = h1['MAP'] - h0['MAP']
    if abs(dmap) > 0.5:
        msgs.append(f"MAP {'rose' if dmap > 0 else 'fell'} ({h0['MAP']:.0f} → {h1['MAP']:.0f} mmHg) — {'increased' if dmap > 0 else 'decreased'} renal perfusion pressure")
    dlvm = h1['LV_mass'] - h0['LV_mass']
    if abs(dlvm) > 0.5:
        rtype = "concentric" if h1['pct_diameter'] > h0['pct_diameter'] + 0.001 else "eccentric" if h1['pct_length'] > h0['pct_length'] + 0.001 else ""
        msgs.append(f"LV mass changed ({h0['LV_mass']:.0f} → {h1['LV_mass']:.0f} g){' — ' + rtype + ' remodeling' if rtype else ''}")
    dedp = h1['LV_EDP_mmHg'] - h0['LV_EDP_mmHg']
    if abs(dedp) > 0.1:
        msgs.append(f"Filling pressure {'rose' if dedp > 0 else 'fell'} (EDP: {h0['LV_EDP_mmHg']:.1f} → {h1['LV_EDP_mmHg']:.1f} mmHg)")
    if not msgs:
        msgs.append("Heart function stable — normal MAP and pressures")
    return msgs

def chain_summary(h0, h1, dt_hours):
    dt_days = dt_hours / 24
    if dt_days > 180:
        dt_str = f"{dt_days/365.25:.1f} years"
    elif dt_days > 1:
        dt_str = f"{dt_days:.0f} days"
    else:
        dt_str = f"{dt_hours:.1f} hours"
    
    lines = [f"Over {dt_str}:"]
    
    dd = h1['pct_diameter'] - h0['pct_diameter']
    dl = h1['pct_length'] - h0['pct_length']
    dlvm = h1['LV_mass'] - h0['LV_mass']
    dbv = h1['blood_volume_L'] - h0['blood_volume_L']
    
    if abs(dbv) > 0.005:
        lines.append(f"  Blood volume: {h0['blood_volume_L']:.2f} → {h1['blood_volume_L']:.2f} L")
    if abs(dd) > 0.001:
        lines.append(f"  Myocyte diameter: {'+' if dd > 0 else ''}{dd:.3f}% (concentric hypertrophy)")
    if abs(dl) > 0.001:
        lines.append(f"  Myocyte length: {'+' if dl > 0 else ''}{dl:.3f}% (eccentric remodeling)")
    if abs(dlvm) > 0.5:
        lines.append(f"  LV mass: {h0['LV_mass']:.0f} → {h1['LV_mass']:.0f} g")
    if h1['LV_EDP_mmHg'] > 15:
        lines.append(f"  ⚠ Elevated filling pressure (EDP={h1['LV_EDP_mmHg']:.1f} mmHg)")
    if h1['BNP'] > 100:
        lines.append(f"  ⚠ BNP elevated ({h1['BNP']:.0f} pg/mL)")
    
    if len(lines) == 1:
        lines.append("  System in homeostatic equilibrium")
    return "\n".join(lines)

# =========================================================================
# ROUTES
# =========================================================================
@app.route('/')
def index():
    return send_from_directory('static', 'index.html')

@app.route('/api/preview', methods=['POST'])
def preview():
    if not sim.lock.acquire(blocking=False):
        return jsonify({'error': 'Simulation already running — please wait'}), 429
    try:
        data = request.json or {}
        sim.set_knobs(data)
        dt_hours = float(data.get('dt_hours', 24))
        result = sim.run_segment(dt_hours)
        if result is None:
            return jsonify({'error': 'Integration failed — solver could not advance from current state'}), 500
        if not result['completed']:
            actual_t = result['trajectory'][-1]['t_hours'] - result['trajectory'][0]['t_hours']
            result['warning'] = f'Integration unstable — solver reached {actual_t:.1f}h of {dt_hours:.0f}h requested. These parameters may be too extreme for the model.'
            print(f"[DEBUG] Preview incomplete: reached {actual_t:.1f}h of {dt_hours:.0f}h")
        return safe_jsonify(result)
    finally:
        sim.lock.release()

@app.route('/api/commit', methods=['POST'])
def commit():
    if not sim.lock.acquire(blocking=False):
        return jsonify({'error': 'Simulation already running — please wait'}), 429
    try:
        data = request.json or {}
        sim.set_knobs(data)
        dt_hours = float(data.get('dt_hours', 24))
        result = sim.run_segment(dt_hours)
        if result is None:
            return jsonify({'error': 'Integration failed — solver could not advance from current state'}), 500
        if not result['completed']:
            actual_t = result['trajectory'][-1]['t_hours'] - result['trajectory'][0]['t_hours']
            print(f"[DEBUG] Commit blocked: integration only reached {actual_t:.1f}h of {dt_hours:.0f}h")
            return safe_jsonify({
                'error': f'Cannot commit — integration unstable at {actual_t:.1f}h of {dt_hours:.0f}h. Try shorter duration or less extreme parameters.',
                'trajectory': result['trajectory'],
                'warning': 'Partial result shown but not committed',
                'completed': False,
            }), 400
        sim.commit(result)
        return safe_jsonify({**result, 'segment_index': len(sim.segments) - 1})
    finally:
        sim.lock.release()

@app.route('/api/history')
def history():
    return jsonify({
        'segments': sim.segments,
        't_hours': sim.t_hours,
    })

@app.route('/api/reset', methods=['POST'])
def reset():
    sim.reset()
    return jsonify({'status': 'reset', 't_hours': 0})

@app.route('/api/state')
def state():
    out = extract_outputs(sim.y, sim.t_hours)
    return jsonify(out)

@app.route('/api/fit', methods=['POST'])
def start_fit():
    # Check if a fit is already running
    for fid, job in fit_jobs.items():
        if job['thread'].is_alive():
            return jsonify({'error': 'A fit is already running', 'fit_id': fid}), 429

    data = request.json or {}
    targets = data.get('targets', [])
    if len(targets) < 3:
        return jsonify({'error': 'Need at least 3 target points'}), 400

    max_iter = int(data.get('max_iterations', 12))
    segment_days = int(data.get('segment_days', 30))
    popsize = int(data.get('popsize', 5))

    fitter = ParamFitter(targets, segment_days=segment_days)
    _fit_counter[0] += 1
    fit_id = str(_fit_counter[0])

    job = {'fitter': fitter, 'thread': None, 'result': None,
           'max_iter': max_iter, 'popsize': popsize}
    fit_jobs[fit_id] = job

    t = threading.Thread(target=_run_fit, args=(fit_id,), daemon=True)
    job['thread'] = t
    t.start()

    return jsonify({'status': 'started', 'fit_id': fit_id})

@app.route('/api/fit/status')
def fit_status():
    # Find the most recent fit
    if not fit_jobs:
        return jsonify({'status': 'none'})

    fit_id = max(fit_jobs.keys(), key=int)
    job = fit_jobs[fit_id]
    fitter = job['fitter']

    with fitter.lock:
        progress = dict(fitter.progress)

    progress['fit_id'] = fit_id

    if job['result'] is not None:
        progress['status'] = 'completed'
        progress.update(job['result'])

    return safe_jsonify(progress)

@app.route('/api/fit/stop', methods=['POST'])
def stop_fit():
    for fid, job in fit_jobs.items():
        if job['thread'].is_alive():
            job['fitter'].stop()
            return jsonify({'status': 'stopping', 'fit_id': fid})
    return jsonify({'status': 'no_fit_running'})

if __name__ == '__main__':
    print("Cardiorenal HFpEF Simulator")
    print("=" * 40)
    print(f"Model: {N_STATE} state variables, {N_PARAM} parameters")
    print(f"Verifying initial state...")
    out = extract_outputs(sim.y, 0)
    print(f"  MAP={out['MAP']} CO={out['CO']} BV={out['blood_volume_L']} Na={out['Na']} EDV={out['LV_EDV_mL']} EDP={out['LV_EDP_mmHg']}")
    print()
    print("Starting server on http://localhost:5010")
    print("Press Ctrl+C to stop")
    app.run(host='0.0.0.0', port=5010, debug=False)
