"""
hallow_c_driver.py — Run the exact Hallow cardiorenal model via compiled C.

Two integration backends:
  1. C-level LSODA integrator (fast, no Python callback overhead)
  2. scipy LSODA fallback

Zero approximations. The C code is auto-generated from modelfile_commented.R.
"""

import ctypes
import numpy as np
from scipy.integrate import ode
import json
import time
import os

_DIR = os.path.dirname(__file__)

# =========================================================================
# Load compiled C libraries
# =========================================================================
LIB_PATH = os.path.join(_DIR, "hallow_rhs.so")
lib = ctypes.CDLL(LIB_PATH)

# void hallow_rhs(double t, const double *y, double *dydt, const double *p)
lib.hallow_rhs.restype = None
lib.hallow_rhs.argtypes = [
    ctypes.c_double,
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
]

# C-level LSODA integrator (optional — falls back to scipy if not compiled)
INTEGRATOR_PATH = os.path.join(_DIR, "hallow_integrator.so")
try:
    integ_lib = ctypes.CDLL(INTEGRATOR_PATH)
    integ_lib.integrate_hallow_lsoda.restype = ctypes.c_int
    integ_lib.integrate_hallow_lsoda.argtypes = [
        ctypes.POINTER(ctypes.c_double),  # y0
        ctypes.POINTER(ctypes.c_double),  # params
        ctypes.c_double,                  # t_end
        ctypes.c_double,                  # dt_output
        ctypes.POINTER(ctypes.c_double),  # out_t
        ctypes.POINTER(ctypes.c_double),  # out_y
        ctypes.c_int,                     # max_output
        ctypes.POINTER(ctypes.c_int),     # n_rhs_calls
        ctypes.c_double,                  # atol
        ctypes.c_double,                  # rtol
        ctypes.c_double,                  # max_step
    ]
    HAS_C_INTEGRATOR = True
    print("[hallow_c_driver] C LSODA integrator loaded")
except OSError:
    HAS_C_INTEGRATOR = False
    print("[hallow_c_driver] C integrator not found, using scipy LSODA fallback")

# Load index maps
STATE_MAP = json.load(open(os.path.join(_DIR, "hallow_rhs_state_map.json")))
PARAM_MAP = json.load(open(os.path.join(_DIR, "hallow_rhs_param_map.json")))
N_STATE = len(STATE_MAP)
N_PARAM = len(PARAM_MAP)

# =========================================================================
# Build parameter array from calcNomParams_timescale.R values
# =========================================================================
def build_params():
    """Build the parameter array matching the C index map.
    Values are transcribed directly from calcNomParams_timescale.R."""
    p = np.zeros(N_PARAM)
    
    # We need to load the R parameter file and extract values.
    # Parse calcNomParams_timescale.R to get name=value pairs.
    param_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "hallow_model", "calcNomParams_timescale.R")
    
    # Strategy: run a minimal R evaluation to get parameter values
    # Since we can't install packages, parse the R file directly
    import re
    
    r_values = {}
    lines = open(param_file).readlines()
    in_func = False
    brace = 0
    
    # Math functions for eval
    safe_ns = {
        '__builtins__': {},
        'log': np.log, 'exp': np.exp, 'sqrt': np.sqrt,
        'max': max, 'min': min, 'abs': abs, 'floor': np.floor,
        'ceiling': np.ceil, 'log2': np.log2, 'log10': np.log10,
    }
    
    for line in lines:
        s = line.strip()
        if 'calcNomParams' in s and '<-' in s or 'function' in s:
            in_func = True
            continue
        if not in_func:
            continue
        brace += s.count('{') - s.count('}')
        if 'return' in s and 'param' in s:
            break
        if s.startswith('#') or s == '' or s.startswith('for') or s.startswith('param'):
            continue
        if s.startswith('t=') and 'sort' in s:
            continue
        
        # Assignment
        m = re.match(r'^(\w+)\s*=\s*(.+?)(?:\s*#.*)?$', s)
        if m:
            name = m.group(1)
            expr = m.group(2).strip().rstrip(';')
            # Convert R syntax to Python: ^ → **
            expr = expr.replace('^', '**')
            try:
                val = eval(expr, {**safe_ns, **r_values})
                r_values[name] = float(val)
            except:
                pass
    
    # Map to parameter array
    for name, idx in PARAM_MAP.items():
        if name in r_values:
            p[idx] = r_values[name]
    
    return p, r_values

# =========================================================================
# Build initial conditions from getInits.R
# =========================================================================
def build_inits(r_values):
    """Build initial state vector matching the C state index map."""
    y0 = np.zeros(N_STATE)
    
    # Direct from getInits.R
    theta = r_values  # parameter values needed for inits
    
    y0[STATE_MAP['venous_volume']] = theta.get('venous_volume_0', 0.003278) - 0.0001
    y0[STATE_MAP['LV_volume']] = theta.get('LV_EDV_nom', 110e-6)
    y0[STATE_MAP['arterial_volume']] = theta.get('arterial_volume_0', 0.000454) + 0.0001
    y0[STATE_MAP['peripheral_circulation_volume']] = theta.get('peripheral_circulation_volume_0', 0.000424)
    y0[STATE_MAP['RV_volume']] = theta.get('RV_volume_0', 0.000075)
    y0[STATE_MAP['pulmonary_arterial_volume']] = theta.get('pulmonary_arterial_volume_0', 0.000041) + 0.00005
    y0[STATE_MAP['pulmonary_venous_volume']] = theta.get('pulmonary_venous_volume_0', 0.000258) - 0.00005
    y0[STATE_MAP['aortic_blood_flow_delayed']] = 0
    y0[STATE_MAP['pulmonary_blood_flow_delayed']] = 0
    y0[STATE_MAP['change_in_myocyte_length']] = 0
    y0[STATE_MAP['change_in_myocyte_diameter']] = 0
    y0[STATE_MAP['LV_active_stress_peak']] = 50000
    y0[STATE_MAP['sim_time']] = 0
    y0[STATE_MAP['LV_sarcomere_length_delayed']] = 2.1e-6
    y0[STATE_MAP['RV_sarcomere_length_delayed']] = 2e-6
    y0[STATE_MAP['LV_EDV']] = theta.get('LV_EDV_nom', 110e-6)
    y0[STATE_MAP['LV_EDP']] = 1000
    y0[STATE_MAP['LV_EDS']] = 3000
    
    map_sp = theta.get('nominal_map_setpoint', 85)
    y0[STATE_MAP['arterial_pressure_delayed']] = map_sp / 0.0075
    y0[STATE_MAP['arterial_pressure_bigger_delay']] = map_sp / 0.0075
    y0[STATE_MAP['systolic_pressure']] = (map_sp + 21) / 0.0075
    y0[STATE_MAP['diastolic_pressure']] = (map_sp - 10.5) / 0.0075
    y0[STATE_MAP['venous_pressure_delayed']] = 917
    y0[STATE_MAP['venous_pressure_bigger_delay']] = 917
    y0[STATE_MAP['systolic_venous_pressure']] = 1050
    y0[STATE_MAP['diastolic_venous_pressure']] = 850
    
    y0[STATE_MAP['CO']] = theta.get('CO_nom', 5.0)
    y0[STATE_MAP['CO_delayed']] = theta.get('CO_nom', 5.0)
    y0[STATE_MAP['AngI']] = 8.164
    y0[STATE_MAP['AngII']] = 5.17
    y0[STATE_MAP['AT1_bound_AngII']] = 16.6
    y0[STATE_MAP['AT2_bound_AngII']] = 5.5
    y0[STATE_MAP['plasma_renin_concentration']] = 17.845
    
    bv = theta.get('blood_volume_nom', 5.0)
    y0[STATE_MAP['blood_volume_L']] = bv
    y0[STATE_MAP['interstitial_fluid_volume']] = theta.get('IF_nom', 15.0)
    ref_na = theta.get('ref_Na_concentration', 140)
    y0[STATE_MAP['sodium_amount']] = bv * ref_na
    y0[STATE_MAP['IF_sodium_amount']] = theta.get('IF_nom', 15.0) * ref_na
    y0[STATE_MAP['stored_sodium']] = 0
    y0[STATE_MAP['tubulo_glomerular_feedback_effect']] = 1
    y0[STATE_MAP['normalized_aldosterone_level']] = 1
    y0[STATE_MAP['preafferent_pressure_autoreg_signal']] = 1
    y0[STATE_MAP['glomerular_pressure_autoreg_signal']] = 1
    y0[STATE_MAP['CO_error']] = 0
    y0[STATE_MAP['Na_concentration_error']] = 0
    y0[STATE_MAP['normalized_vasopressin_concentration_delayed']] = 1
    
    # nom_LoH_Na_outflow — need to compute from params
    nom_LoH = theta.get('nom_LoH_Na_outflow', 0.0)
    y0[STATE_MAP['F0_TGF']] = nom_LoH
    y0[STATE_MAP['P_bowmans']] = theta.get('Pc_pt_s1_mmHg', 20.2)
    y0[STATE_MAP['oncotic_pressure_difference']] = theta.get('nom_oncotic_pressure_difference', 28)
    y0[STATE_MAP['renal_blood_flow_L_min_delayed']] = theta.get('nom_renal_blood_flow_L_min', 1.0)
    y0[STATE_MAP['SN_macula_densa_Na_flow_delayed']] = nom_LoH / theta.get('baseline_nephrons', 2e6) if nom_LoH > 0 else 0
    y0[STATE_MAP['rsna_delayed']] = 1
    y0[STATE_MAP['disease_effects_increasing_Kf']] = 0
    y0[STATE_MAP['disease_effects_decreasing_CD_PN']] = 0
    y0[STATE_MAP['tubular_length_increase']] = 0
    y0[STATE_MAP['tubular_diameter_increase']] = 0
    y0[STATE_MAP['water_out_s1_delayed']] = 3e-8
    y0[STATE_MAP['water_out_s2_delayed']] = 1.9e-8
    y0[STATE_MAP['water_out_s3_delayed']] = 1.2e-8
    y0[STATE_MAP['reabsorbed_urea_cd_delayed']] = 0
    y0[STATE_MAP['UGE']] = 0
    y0[STATE_MAP['serum_creatinine']] = theta.get('equilibrium_serum_creatinine', 0.92) * bv
    y0[STATE_MAP['cumNaExcretion']] = 0
    y0[STATE_MAP['cumWaterExcretion']] = 0
    y0[STATE_MAP['cumCreatinineExcretion']] = 0
    y0[STATE_MAP['RTg_compensation']] = 0
    y0[STATE_MAP['SGLT2_inhibition_delayed']] = 1
    y0[STATE_MAP['RUGE_delayed']] = 0
    y0[STATE_MAP['postglomerular_pressure_delayed']] = theta.get('RIHP0', 9.32)
    y0[STATE_MAP['postglomerular_pressure_error']] = 0
    y0[STATE_MAP['mitral_valve_leak']] = 0
    
    return y0

# =========================================================================
# RHS wrapper for scipy — optimized to minimize Python↔C overhead
# =========================================================================
# Pre-build a C function handle with void* args to avoid per-call POINTER creation
_rhs_cfunc = ctypes.CFUNCTYPE(
    None, ctypes.c_double, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p
)(('hallow_rhs', lib))

class HallowRHS:
    """Fast RHS wrapper: pre-allocates dydt buffer and uses raw data pointers."""
    def __init__(self, params):
        self.params = params
        self._p_data = params.ctypes.data
        self._dydt = np.zeros(N_STATE)
        self._d_data = self._dydt.ctypes.data

    def __call__(self, t, y):
        _rhs_cfunc(t, y.ctypes.data, self._d_data, self._p_data)
        return self._dydt

# =========================================================================
# C-level integration (fast path)
# =========================================================================
def integrate_c(y0, params, t_end, dt_output=0.01, atol=1e-6, rtol=1e-4, max_step=0.01):
    """Integrate using the C-level LSODA integrator (no Python callbacks)."""
    max_output = int(t_end / dt_output) + 10
    out_t = np.zeros(max_output)
    out_y = np.zeros(max_output * N_STATE)
    n_rhs = ctypes.c_int(0)

    n_pts = integ_lib.integrate_hallow_lsoda(
        y0.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        params.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        t_end, dt_output,
        out_t.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        out_y.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        max_output,
        ctypes.byref(n_rhs),
        atol, rtol, max_step,
    )

    if n_pts < 0:
        error_names = {-1: "solver failure", -2: "NaN detected", -3: "max steps"}
        print(f"[DEBUG] C LSODA error: {error_names.get(n_pts, f'code {n_pts}')}")
        return [], False

    results = []
    completed = True
    for i in range(n_pts):
        yi = out_y[i*N_STATE:(i+1)*N_STATE].copy()
        results.append((out_t[i], yi))
        if np.any(np.isnan(yi)):
            print(f"[DEBUG] C LSODA: NaN at t={out_t[i]:.4f}")
            completed = False
            break

    if n_pts > 0 and results[-1][0] < t_end - 1e-6:
        completed = False

    print(f"[DEBUG] C LSODA {'completed' if completed else 'INCOMPLETE'}: {len(results)} points, "
          f"t=0..{results[-1][0]:.4f} of {t_end:.4f} hours")
    return results, completed

# =========================================================================
# scipy LSODA integration (fallback)
# =========================================================================
def integrate_scipy(y0, params, t_end, dt_output=0.01):
    """Integrate using scipy LSODA (slower but robust fallback)."""
    rhs = HallowRHS(params)

    solver = ode(rhs)
    solver.set_integrator('lsoda', atol=1e-6, rtol=1e-4, max_step=0.01,
                          nsteps=500000)
    solver.set_initial_value(y0, 0)

    results = [(0, y0.copy())]
    t_points = np.arange(dt_output, t_end + dt_output/2, dt_output)
    completed = True

    for t in t_points:
        solver.integrate(t)
        if not solver.successful():
            print(f"[DEBUG] Integration failed at t={t:.4f} hours")
            completed = False
            break
        if np.any(np.isnan(solver.y)):
            nan_indices = np.where(np.isnan(solver.y))[0]
            print(f"[DEBUG] NaN detected at t={t:.4f} hours in state indices: {nan_indices.tolist()}")
            completed = False
            break
        if np.any(np.isinf(solver.y)):
            inf_indices = np.where(np.isinf(solver.y))[0]
            print(f"[DEBUG] Inf detected at t={t:.4f} hours in state indices: {inf_indices.tolist()}")
            completed = False
            break
        results.append((t, solver.y.copy()))

    print(f"[DEBUG] scipy {'completed' if completed else 'INCOMPLETE'}: {len(results)} points, t=0..{results[-1][0]:.4f} of {t_end:.4f} hours")
    return results, completed

# =========================================================================
# Default integrator
# =========================================================================
def integrate(y0, params, t_end, dt_output=0.01):
    """Integrate the Hallow model from t=0 to t_end hours.

    Uses C-level LSODA when available (fastest), falling back to scipy LSODA.
    Returns (results, completed) where completed is True if integration
    reached t_end without solver failure or NaN/Inf.
    """
    if HAS_C_INTEGRATOR:
        results, completed = integrate_c(y0, params, t_end, dt_output)
        if completed:
            return results, completed
        print("[DEBUG] C LSODA failed, falling back to scipy")
    return integrate_scipy(y0, params, t_end, dt_output)

# =========================================================================
# Parallel batch runner for ML/optimization
# =========================================================================
def _batch_worker(args):
    """Worker function for parallel batch integration.
    Returns only (idx, success, final_state) to minimize pickle overhead."""
    idx, y0, params, t_end, dt_output = args
    try:
        results, completed = integrate_c(y0, params, t_end, dt_output)
        if completed and len(results) > 0:
            return (idx, True, results[-1][1])  # just final state
        return (idx, False, None)
    except Exception:
        return (idx, False, None)

def batch_integrate(y0, param_sets, t_end=24.0, dt_output=1.0,
                    n_workers=None, progress=True):
    """Run many simulations in parallel with different parameter sets.

    Args:
        y0: initial state vector (shared across all sims)
        param_sets: list of parameter arrays, one per simulation
        t_end: simulation duration in hours
        dt_output: output interval in hours
        n_workers: number of parallel workers (default: CPU count)
        progress: print progress updates

    Returns:
        list of (idx, success, final_state) tuples sorted by index.
        final_state is the state vector at t_end, or None if failed.
    """
    import multiprocessing as mp

    if n_workers is None:
        n_workers = min(mp.cpu_count(), len(param_sets))

    # Use fork context for fast subprocess creation (no reimport overhead)
    ctx = mp.get_context('fork')

    tasks = [(i, y0, p, t_end, dt_output) for i, p in enumerate(param_sets)]

    if progress:
        print(f"[batch] {len(tasks)} sims, {n_workers} workers, t_end={t_end}h")

    t0 = time.time()
    results = []

    with ctx.Pool(n_workers) as pool:
        for result in pool.imap_unordered(_batch_worker, tasks, chunksize=1):
            results.append(result)
            if progress and len(results) % max(1, len(tasks) // 5) == 0:
                n_ok = sum(1 for r in results if r[1])
                print(f"[batch] {len(results)}/{len(tasks)} "
                      f"({n_ok} ok) {time.time()-t0:.1f}s")

    results.sort(key=lambda x: x[0])
    elapsed = time.time() - t0
    n_ok = sum(1 for r in results if r[1])

    if progress:
        print(f"[batch] Done: {n_ok}/{len(tasks)} ok in {elapsed:.1f}s "
              f"({elapsed/len(tasks):.2f}s/sim, "
              f"{len(tasks)/elapsed:.1f} sims/s)")

    return results

def extract_batch_outputs(results):
    """Extract key clinical outputs from batch results.

    Returns dict of numpy arrays with NaN for failed sims.
    """
    n = len(results)
    outputs = {
        'MAP': np.full(n, np.nan),
        'CO': np.full(n, np.nan),
        'BV': np.full(n, np.nan),
        'Na': np.full(n, np.nan),
        'EDV': np.full(n, np.nan),
        'EDP': np.full(n, np.nan),
        'success': np.zeros(n, dtype=bool),
    }

    for idx, success, yf in results:
        if not success or yf is None:
            continue
        outputs['success'][idx] = True
        outputs['MAP'][idx] = (yf[STATE_MAP['systolic_pressure']] / 3 +
                                yf[STATE_MAP['diastolic_pressure']] * 2 / 3) * 0.0075
        outputs['CO'][idx] = yf[STATE_MAP['CO_delayed']]
        outputs['BV'][idx] = yf[STATE_MAP['blood_volume_L']]
        bv = yf[STATE_MAP['blood_volume_L']]
        outputs['Na'][idx] = yf[STATE_MAP['sodium_amount']] / bv if bv > 0 else np.nan
        outputs['EDV'][idx] = yf[STATE_MAP['LV_EDV']] * 1e6
        outputs['EDP'][idx] = yf[STATE_MAP['LV_EDP']] * 0.0075

    return outputs


# =========================================================================
# Test
# =========================================================================
if __name__ == '__main__':
    print("Building parameters...")
    params, r_values = build_params()
    print(f"  {N_PARAM} parameters loaded, {sum(params != 0)} non-zero")
    
    print("Building initial conditions...")
    y0 = build_inits(r_values)
    print(f"  {N_STATE} state variables")
    
    # Benchmark: single RHS call
    rhs = HallowRHS(params)
    dydt = rhs(0.0, y0)
    
    print("\nBenchmarking RHS evaluation...")
    t0 = time.time()
    N_calls = 10000
    for _ in range(N_calls):
        rhs(0.0, y0)
    elapsed = time.time() - t0
    print(f"  {N_calls} calls in {elapsed:.3f}s = {elapsed/N_calls*1e6:.1f} µs/call")
    print(f"  Estimated 1-hour sim (~7000 calls): {7000 * elapsed/N_calls:.3f}s")
    print(f"  Estimated 24-hour sim: {24 * 7000 * elapsed/N_calls:.3f}s")
    
    # Short integration test
    print("\nIntegrating 0.05 hours (~4 beats)...")
    t0 = time.time()
    results, _ = integrate(y0, params, 0.05, dt_output=0.01)
    elapsed = time.time() - t0
    print(f"  Done in {elapsed:.3f}s ({len(results)} output points)")

    # Extract outputs from final state
    _, y_final = results[-1]
    MAP = (y_final[STATE_MAP['systolic_pressure']] / 3 + 
           y_final[STATE_MAP['diastolic_pressure']] * 2 / 3) * 0.0075
    CO = y_final[STATE_MAP['CO_delayed']]
    BV = y_final[STATE_MAP['blood_volume_L']]
    Na = y_final[STATE_MAP['sodium_amount']] / BV if BV > 0 else 0
    EDV = y_final[STATE_MAP['LV_EDV']] * 1e6
    EDP = y_final[STATE_MAP['LV_EDP']] * 0.0075
    
    print(f"\n=== Results at t=0.05 hours ===")
    print(f"  MAP:  {MAP:.1f} mmHg  (target: ~85)")
    print(f"  CO:   {CO:.2f} L/min  (target: ~5.0)")
    print(f"  BV:   {BV:.3f} L     (target: ~5.0)")
    print(f"  Na:   {Na:.1f} mEq/L  (target: ~140)")
    print(f"  EDV:  {EDV:.1f} mL    (target: ~110)")
    print(f"  EDP:  {EDP:.1f} mmHg  (target: 7-16)")
    
    # Try longer integration
    print("\nIntegrating 1 hour (~70 beats)...")
    t0 = time.time()
    results_1h, _ = integrate(y0, params, 1.0, dt_output=0.1)
    elapsed_1h = time.time() - t0
    print(f"  Done in {elapsed_1h:.3f}s")
    _, y1h = results_1h[-1]
    MAP1h = (y1h[STATE_MAP['systolic_pressure']] / 3 + y1h[STATE_MAP['diastolic_pressure']] * 2 / 3) * 0.0075
    print(f"  MAP={MAP1h:.1f} CO={y1h[STATE_MAP['CO_delayed']]:.2f} BV={y1h[STATE_MAP['blood_volume_L']]:.3f}")
    
    print("\nIntegrating 24 hours (~1680 beats)...")
    t0 = time.time()
    results_24h, _ = integrate(y0, params, 24.0, dt_output=1.0)
    elapsed_24h = time.time() - t0
    print(f"  Done in {elapsed_24h:.3f}s")
    _, y24h = results_24h[-1]
    MAP24h = (y24h[STATE_MAP['systolic_pressure']] / 3 + y24h[STATE_MAP['diastolic_pressure']] * 2 / 3) * 0.0075
    print(f"  MAP={MAP24h:.1f} CO={y24h[STATE_MAP['CO_delayed']]:.2f} BV={y24h[STATE_MAP['blood_volume_L']]:.3f}")
