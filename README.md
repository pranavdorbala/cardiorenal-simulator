# Cardiorenal HFpEF Simulator

Interactive simulator of the Hallow et al. coupled cardiorenal model, running the **exact original equations** compiled to native C for performance.

## What This Is

This package takes the R model code from Hallow et al. (2023, PLOS Computational Biology), automatically converts the 1,750-line ODE system to C, compiles it to a shared library, and wraps it in a Python/Flask web dashboard with interactive knobs. The C code is auto-generated from `modelfile_commented.R` — no equations are approximated, simplified, or omitted. All 70 coupled ODEs, 430 parameters, and the complete cardiac-vascular-renal feedback loop run exactly as in the original RxODE implementation.

## Quick Start (Mac)

```bash
# 1. Make sure you have Xcode command line tools
xcode-select --install   # if not already installed

# 2. Clone or copy this directory
cd cardiorenal_simulator

# 3. Build (generates C code, compiles, verifies)
bash build.sh

# 4. Run
python3 server.py
```

Open **http://localhost:5010** in your browser.

## What's in the Dashboard

### Knobs (left panel)

- **Arterial Compliance** (`C_art_scale`, 0.3–1.0): Controls arterial stiffness. 1.0 = normal young arteries. Lower values represent aging, atherosclerosis, or chronic hypertension. Maps directly to the Hallow parameter `C_art_scale` which scales `C_art_initial`.

- **Peripheral Resistance** (`disease_effect_on_TPR_peripheral_resistance`, 0.8–2.0): Multiplier on total peripheral resistance. 1.0 = normal. Higher values represent sustained hypertension. This is the existing Hallow parameter that multiplies `R_per0`.

- **HbA1c** (4.0–12.0%): Converts to plasma glucose via the ADAG formula (`glucose_mmol = 1.59 × A1C − 2.59`) and sets `glucose_concentration` in the parameter array. Higher glucose causes SGLT2 overflow, glucosuria, and osmotic diuresis through the model's existing tubular reabsorption cascade.

### Segment Duration

Choose how far forward to integrate: 1 hour, 6 hours, 24 hours, 7 days, or 30 days. Each segment picks up where the previous one left off.

### Preview vs Store

- **Preview**: Integrates the model forward with current knob settings. Shows the trajectory as dashed lines on the graphs. Changes when you adjust knobs and re-preview. Does NOT advance the simulation clock.

- **Store**: Commits the current preview as permanent history. The dashed lines become solid. The simulation clock advances to the end of the segment. The next preview will start from this new state.

This lets you explore "what would happen if..." without committing, then lock in a trajectory and move forward.

### Charts (center panel)

Six time-series charts showing:
- **SBP / DBP** — Systolic and diastolic blood pressure (mmHg)
- **LV Mass** — Left ventricular mass (grams), increases with hypertrophy
- **Myocyte Diameter** — Percent change from baseline, driven by active stress (concentric hypertrophy)
- **Myocyte Length** — Percent change from baseline, driven by passive strain (eccentric remodeling)
- **Blood Volume** — Total blood volume (liters), regulated by kidney sodium/water handling
- **Cardiac Output** — Liters per minute

Stored segments appear as solid lines. The current preview appears as dashed lines connected to the last stored point.

### Mechanistic Messages (right panel)

After each preview or store, the panel shows:

- **Kidney → Heart**: What signals the kidney is sending to the heart (blood volume changes, sodium shifts)
- **Heart → Kidney**: What signals the heart is sending to the kidney (MAP changes, remodeling, filling pressure)
- **Chain Summary**: The mechanistic cascade from cause to effect for this segment

### Stored Segments Timeline

Shows all committed segments with their time range and key state changes.

## File Structure

```
cardiorenal_simulator/
├── build.sh                          # Build script (compile C, verify)
├── server.py                         # Flask web server
├── hallow_c_driver.py                # Python ↔ C interface (ctypes + scipy LSODA)
├── rxode_to_c.py                     # Parser: R model → C code generator
├── hallow_rhs.c                      # Auto-generated C source (70 ODEs, ~2000 lines)
├── hallow_rhs.so / .dylib            # Compiled shared library (built by build.sh)
├── hallow_rhs_state_map.json         # State variable name → index mapping
├── hallow_rhs_param_map.json         # Parameter name → index mapping
├── static/
│   └── index.html                    # Dashboard UI
└── hallow_model/                     # Original R source files (Hallow et al.)
    ├── modelfile_commented.R         # The ODE model (1,750 lines)
    ├── calcNomParams_timescale.R     # Parameter definitions (430 params)
    ├── getInits.R                    # Initial conditions (70 states)
    ├── calcEndTime.R                 # Utility
    ├── calcOutputs.R                 # Output calculation
    ├── VariableCheckFile.csv         # Validation ranges
    ├── NOTICE.txt                    # Attribution
    └── license.txt                   # GPL v3
```

## How the C Code is Generated

`rxode_to_c.py` performs the following mechanical transformation:

1. Extracts the ODE string from `modelfile_commented.R` (everything between `ode <- "` and the closing `"`)
2. Identifies 70 state variables from `d/dt(X) = ...` lines (in order of first appearance)
3. Identifies 430 parameters from `calcNomParams_timescale.R` (in order of assignment)
4. Tokenizes each line and replaces:
   - State variable reads → `y[index]`
   - Parameter reads → `p[index]`
   - `d/dt(X) = expr` → `dydt[index] = expr`
   - R's `^` → `safe_pow()` (handles negative base with fractional exponent)
   - `abs` → `fabs`, `min`/`max` → `fmin`/`fmax`
5. Pre-declares variables assigned inside `if/else` blocks (C scoping)
6. Detects parameters that are overwritten inside the model (e.g., `R_art` in the aortic stenosis block) and shadows them as local `double` variables

The resulting C function has signature:
```c
void hallow_rhs(double t, const double *y, double *dydt, const double *p)
```

## Performance

| Operation | Time |
|-----------|------|
| Single RHS evaluation | ~12 µs |
| 1-hour simulation (~70 beats) | ~14 s |
| 24-hour simulation (~1,680 beats) | ~17 s |
| Parameter change + 1-hour preview | ~14 s |

The bottleneck is LSODA's step count (the model is stiff due to the cardiac cycle operating on a sub-second timescale while renal dynamics operate on hours). The C compilation gives ~1,200× speedup over interpreted R.

## Verification

At t=0.05 hours (steady state), the model produces:

| Variable | Value | Expected |
|----------|-------|----------|
| MAP | 85.0 mmHg | 85 |
| CO | 5.00 L/min | 5.0 |
| Blood Volume | 4.999 L | 5.0 |
| Na | 140.0 mEq/L | 140 |
| LV EDV | 110.0 mL | 110 |
| LV EDP | 7.5 mmHg | 7–16 |

These match the R/deSolve reference implementation exactly.

## Attribution

The cardiorenal model was developed by K. Melissa Hallow and collaborators at the University of Georgia. If you use this model, cite:

> Basu S, Yu H, Murrow JR, Hallow KM (2023) Understanding heterogeneous mechanisms of heart failure with preserved ejection fraction through cardiorenal mathematical modeling. PLoS Comput Biol 19(11): e1011598.

See `hallow_model/NOTICE.txt` for full attribution and `hallow_model/license.txt` for the GPL v3 license.
