# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Interactive simulator of the Hallow et al. coupled cardiorenal HFpEF model. The R model code (1,750-line ODE system with 70 state variables and 430 parameters) is auto-converted to C, compiled to a shared library, and served via a Flask web dashboard.

**Citation**: Basu S, Yu H, Murrow JR, Hallow KM (2023) PLoS Comput Biol 19(11): e1011598. Licensed under GPL v3.

## Build & Run

```bash
# Build (generates C from R model, compiles shared library, verifies outputs)
bash build.sh

# Run the web server (serves on http://localhost:5010)
python3 server.py

# Run the C driver standalone (benchmarks + integration tests)
python3 hallow_c_driver.py

# Regenerate C code from R model (without compiling)
python3 rxode_to_c.py
```

Requirements: Python 3.8+, gcc/Xcode CLI tools, flask, scipy, numpy.

## Architecture

### Pipeline: R model → C → Python → Web

1. **`rxode_to_c.py`** — Parses `hallow_model/modelfile_commented.R` to extract the ODE string, tokenizes it, and emits `hallow_rhs.c`. Also parses `calcNomParams_timescale.R` for parameter names. Outputs JSON index maps (`hallow_rhs_state_map.json`, `hallow_rhs_param_map.json`).

2. **`hallow_rhs.c`** (auto-generated) — Contains a single function `void hallow_rhs(double t, const double *y, double *dydt, const double *p)`. State variables are accessed as `y[idx]`, parameters as `p[idx]`, derivatives written to `dydt[idx]`. Uses `safe_pow()` for R's `^` operator. Compiled to `.dylib` (macOS) or `.so` (Linux).

3. **`hallow_c_driver.py`** — Loads the compiled shared library via ctypes. Provides `build_params()` (parses R parameter file and eval's expressions), `build_inits()` (sets initial conditions from `getInits.R` logic), `HallowRHS` class (wraps the C function for scipy), and `integrate()` (uses scipy's LSODA solver). Exports `STATE_MAP`, `PARAM_MAP`, `N_STATE`, `N_PARAM`.

4. **`server.py`** — Flask app with `SimState` class managing simulation state. API endpoints: `/api/preview` (integrate without committing), `/api/commit` (integrate and store), `/api/reset`, `/api/state`, `/api/history`. Extracts clinical outputs (MAP, SBP/DBP, CO, LV mass, BNP, etc.) from raw state vectors and generates mechanistic messages (kidney→heart, heart→kidney feedback).

5. **`static/index.html`** — Single-file dashboard with knobs (arterial compliance, peripheral resistance, HbA1c), segment duration picker, preview/store workflow, six time-series charts, and mechanistic message panel.

### Key Concepts

- **Segment-based interaction**: Users preview a time segment (dashed lines), then optionally commit it (solid lines). Previews don't advance the simulation clock.
- **Knobs map to model parameters**: `C_art_scale` → arterial compliance, `disease_effect_on_TPR_peripheral_resistance` → peripheral resistance, `glucose_concentration` → derived from HbA1c via ADAG formula.
- **Unit conversions**: Internal pressure is in dyn/cm², displayed in mmHg (×0.0075). Volumes in m³ internally, displayed in mL (×1e6) or L.
- **The model is stiff**: Cardiac cycle runs sub-second while renal dynamics run on hours. LSODA handles this with `max_step=0.01`, `nsteps=500000`.

### R Source Files (in `hallow_model/`)

- `modelfile_commented.R` — The ODE model (source of truth for all equations)
- `calcNomParams_timescale.R` — All 430 parameter definitions with physiological values
- `getInits.R` — Initial conditions for the 70 state variables
- `calcOutputs.R`, `calcEndTime.R` — Utilities from original R implementation

## Key Conventions

- `hallow_rhs.c` is auto-generated — do not edit directly. Modify `rxode_to_c.py` or the R source instead.
- Parameter and state variable indices are defined by parse order in the R files and stored in the JSON map files.
- The `build_params()` function in `hallow_c_driver.py` uses Python `eval()` with a restricted namespace to evaluate R parameter expressions — be careful when modifying the R parameter file parser.
