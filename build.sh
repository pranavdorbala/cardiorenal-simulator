#!/bin/bash
# build.sh — Compile the Hallow cardiorenal model C code and set up the simulator
#
# Usage: cd cardiorenal_simulator && bash build.sh
#
# Requirements: Python 3.8+, pip, gcc (comes with Xcode command line tools on Mac)

set -e

echo "========================================"
echo "Cardiorenal HFpEF Simulator — Build"
echo "========================================"
echo

# Check Python
echo "Checking Python..."
python3 --version || { echo "ERROR: python3 not found. Install Python 3.8+"; exit 1; }

# Check gcc
echo "Checking C compiler..."
gcc --version 2>/dev/null | head -1 || cc --version 2>/dev/null | head -1 || {
    echo "ERROR: No C compiler found."
    echo "On Mac: xcode-select --install"
    echo "On Linux: sudo apt install gcc"
    exit 1
}

# Install Python dependencies
echo
echo "Installing Python dependencies..."
pip3 install flask scipy numpy --quiet

# Step 1: Generate C code from the R model (if not already done)
if [ ! -f hallow_rhs.c ]; then
    echo
    echo "Step 1: Generating C code from R model..."
    python3 rxode_to_c.py
else
    echo
    echo "Step 1: hallow_rhs.c already exists, skipping generation."
fi

# Step 2: Compile to shared libraries
echo
echo "Step 2: Compiling C code..."
UNAME=$(uname -s)
if [ "$UNAME" = "Darwin" ]; then
    # macOS
    echo "  Platform: macOS"
    gcc -O2 -dynamiclib -fPIC -o hallow_rhs.dylib hallow_rhs.c -lm
    ln -sf hallow_rhs.dylib hallow_rhs.so
    echo "  Output: hallow_rhs.dylib (+ .so symlink)"
    echo "  Compiling C integrator (LSODA)..."
    gcc -O2 -dynamiclib -fPIC -o hallow_integrator.dylib \
        hallow_integrator.c hallow_rhs.c \
        lsoda/lsoda.c lsoda/common.c lsoda/vmnorm.c lsoda/stoda.c \
        lsoda/prja.c lsoda/solsy.c lsoda/intdy.c \
        lsoda/cfode.c lsoda/cfode_static.c \
        lsoda/scaleh.c lsoda/methodswitch.c \
        lsoda/orderswitch.c lsoda/correction.c lsoda/corfailure.c \
        lsoda/fnorm.c lsoda/strdup_printf.c \
        lsoda/daxpy.c lsoda/ddot.c lsoda/dgefa.c lsoda/dgesl.c \
        lsoda/dscal.c lsoda/idamax.c \
        -Ilsoda -lm
    ln -sf hallow_integrator.dylib hallow_integrator.so
    echo "  Output: hallow_integrator.dylib (+ .so symlink)"
else
    # Linux
    echo "  Platform: Linux"
    gcc -O2 -shared -fPIC -o hallow_rhs.so hallow_rhs.c -lm
    echo "  Output: hallow_rhs.so"
    echo "  Compiling C integrator (LSODA)..."
    gcc -O2 -shared -fPIC -o hallow_integrator.so \
        hallow_integrator.c hallow_rhs.c \
        lsoda/lsoda.c lsoda/common.c lsoda/vmnorm.c lsoda/stoda.c \
        lsoda/prja.c lsoda/solsy.c lsoda/intdy.c \
        lsoda/cfode.c lsoda/cfode_static.c \
        lsoda/scaleh.c lsoda/methodswitch.c \
        lsoda/orderswitch.c lsoda/correction.c lsoda/corfailure.c \
        lsoda/fnorm.c lsoda/strdup_printf.c \
        lsoda/daxpy.c lsoda/ddot.c lsoda/dgefa.c lsoda/dgesl.c \
        lsoda/dscal.c lsoda/idamax.c \
        -Ilsoda -lm
    echo "  Output: hallow_integrator.so"
fi

# Step 3: Verify
echo
echo "Step 3: Verifying model..."
python3 -c "
from hallow_c_driver import build_params, build_inits, HallowRHS, integrate, STATE_MAP
import numpy as np
params, r_values = build_params()
y0 = build_inits(r_values)
results, _ = integrate(y0, params, 0.05, dt_output=0.01)
_, yf = results[-1]
MAP = (yf[STATE_MAP['systolic_pressure']]/3 + yf[STATE_MAP['diastolic_pressure']]*2/3) * 0.0075
CO = yf[STATE_MAP['CO_delayed']]
BV = yf[STATE_MAP['blood_volume_L']]
Na = yf[STATE_MAP['sodium_amount']] / BV
EDV = yf[STATE_MAP['LV_EDV']] * 1e6
EDP = yf[STATE_MAP['LV_EDP']] * 0.0075

print(f'  MAP = {MAP:.1f} mmHg  (expected: 85.0)')
print(f'  CO  = {CO:.2f} L/min  (expected: 5.00)')
print(f'  BV  = {BV:.3f} L     (expected: ~5.0)')
print(f'  Na  = {Na:.1f} mEq/L  (expected: 140)')
print(f'  EDV = {EDV:.1f} mL    (expected: 110)')
print(f'  EDP = {EDP:.1f} mmHg  (expected: 7-16)')

ok = abs(MAP - 85) < 1 and abs(CO - 5) < 0.1 and abs(BV - 5) < 0.1
if ok:
    print('  ✓ All outputs match expected values')
else:
    print('  ✗ OUTPUTS DO NOT MATCH — check compilation')
    exit(1)
"

echo
echo "========================================"
echo "Build successful!"
echo ""
echo "To start the simulator:"
echo "  python3 server.py"
echo ""
echo "Then open http://localhost:5010"
echo "========================================"
