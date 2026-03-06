"""
rxode_to_c.py — Parse the Hallow RxODE model and emit a correct C shared library.

Strategy:
  1. Extract the ODE string from modelfile_commented.R
  2. Extract state variable names (from d/dt lines) and parameter names (from calcNomParams)
  3. For each line, tokenize and replace variable references:
     - State variables → y[index] on read, dydt[index] on d/dt write
     - Parameters → p[index]
     - Everything else → local double variable
  4. Emit valid C with #include <math.h> and proper min/max macros
  5. Compile to shared library
"""

import re
import os
import subprocess
import sys

PROJECT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "hallow_model")

# =========================================================================
# Step 1: Extract ODE string from modelfile_commented.R
# =========================================================================
def extract_ode_string(model_file):
    lines = open(model_file).readlines()
    ode_lines = []
    in_ode = False
    for line in lines:
        if line.strip().startswith('ode <- "'):
            in_ode = True
            continue
        if in_ode and line.strip() == '"':
            in_ode = False
            continue
        if in_ode:
            ode_lines.append(line.rstrip('\n'))
    return ode_lines

# =========================================================================
# Step 2: Extract state variable names (in d/dt order) and parameter names
# =========================================================================
def extract_state_vars(ode_lines):
    """Extract state variable names in the order they first appear as d/dt targets."""
    state_vars = []
    for line in ode_lines:
        m = re.search(r'd/dt\((\w+)\)', line)
        if m:
            name = m.group(1)
            if name not in state_vars:
                state_vars.append(name)
    return state_vars

def extract_param_names(param_file):
    """Run the R calcNomParams function and extract parameter names."""
    # We parse the R file to find all variable assignments at the top level
    # This is simpler and more reliable than running R
    lines = open(param_file).readlines()
    params = []
    in_function = False
    brace_depth = 0
    for line in lines:
        stripped = line.strip()
        if stripped.startswith('calcNomParams'):
            in_function = True
            continue
        if in_function:
            brace_depth += stripped.count('{') - stripped.count('}')
            if brace_depth <= 0 and 'return' in stripped:
                break
            # Look for simple assignments: name = value or name = expr
            # Skip lines that are comments, control flow, or complex
            if stripped.startswith('#') or stripped == '' or stripped.startswith('if') or stripped.startswith('for'):
                continue
            if stripped.startswith('param') or stripped.startswith('t=') or stripped.startswith('return'):
                continue
            m = re.match(r'^(\w+)\s*=\s*', stripped)
            if m:
                name = m.group(1)
                # Skip temp variables used in calculations
                if name not in ['i', 't', 'param'] and name not in params:
                    params.append(name)
    return params

# =========================================================================
# Step 3: Tokenizer — split a line into tokens preserving structure
# =========================================================================
# R identifiers: letters, digits, underscores, dots (but dots are rare in this code)
IDENT_RE = re.compile(r'[A-Za-z_][A-Za-z0-9_]*')

def tokenize(line):
    """Split line into tokens: identifiers, numbers, operators, whitespace."""
    tokens = []
    i = 0
    while i < len(line):
        c = line[i]
        # Whitespace
        if c in ' \t':
            j = i
            while j < len(line) and line[j] in ' \t':
                j += 1
            tokens.append(('WS', line[i:j]))
            i = j
        # Identifier or keyword
        elif c.isalpha() or c == '_':
            m = IDENT_RE.match(line, i)
            tokens.append(('ID', m.group()))
            i = m.end()
        # Number (including scientific notation)
        elif c.isdigit() or (c == '.' and i+1 < len(line) and line[i+1].isdigit()):
            j = i
            while j < len(line) and (line[j].isdigit() or line[j] in '.eE+-'):
                if line[j] in 'eE' and j+1 < len(line) and line[j+1] in '+-':
                    j += 2
                elif line[j] in '+-' and j > i and line[j-1] not in 'eE':
                    break
                else:
                    j += 1
            tokens.append(('NUM', line[i:j]))
            i = j
        # Comment
        elif c == '#':
            tokens.append(('COMMENT', line[i:]))
            i = len(line)
        # Operators and punctuation
        elif c == '*' and i+1 < len(line) and line[i+1] == '*':
            tokens.append(('OP', '**'))
            i += 2
        elif c == '=' and i+1 < len(line) and line[i+1] == '=':
            tokens.append(('OP', '=='))
            i += 2
        elif c == '!' and i+1 < len(line) and line[i+1] == '=':
            tokens.append(('OP', '!='))
            i += 2
        elif c == '>' and i+1 < len(line) and line[i+1] == '=':
            tokens.append(('OP', '>='))
            i += 2
        elif c == '<' and i+1 < len(line) and line[i+1] == '=':
            tokens.append(('OP', '<='))
            i += 2
        elif c == '&' and i+1 < len(line) and line[i+1] == '&':
            tokens.append(('OP', '&&'))
            i += 2
        elif c == '|' and i+1 < len(line) and line[i+1] == '|':
            tokens.append(('OP', '||'))
            i += 2
        else:
            tokens.append(('OP', c))
            i += 1
    return tokens

# =========================================================================
# Step 4: Convert tokens to C
# =========================================================================
def convert_line_to_c(line, state_idx, param_idx, local_vars, in_ddt=False):
    """Convert a single line of RxODE code to C."""
    stripped = line.strip()
    
    # Empty line
    if not stripped:
        return ''
    
    # Comment
    if stripped.startswith('#'):
        return '  //' + stripped[1:]
    
    # Check if this is a d/dt line
    ddt_match = re.match(r'\s*d/dt\((\w+)\)\s*=\s*(.*)', line)
    if ddt_match:
        varname = ddt_match.group(1)
        expr = ddt_match.group(2)
        idx = state_idx[varname]
        c_expr = convert_expr(expr, state_idx, param_idx, local_vars)
        return f'  dydt[{idx}] = {c_expr}  /* d/dt({varname}) */'
    
    # Check if this is an if/else
    if stripped.startswith('if ') or stripped.startswith('if('):
        c_line = convert_expr(stripped, state_idx, param_idx, local_vars)
        return '  ' + c_line
    if stripped.startswith('} else'):
        c_line = convert_expr(stripped, state_idx, param_idx, local_vars)
        return '  ' + c_line
    if stripped == '}' or stripped == '{':
        return '  ' + stripped
    
    # Assignment: varname = expr  or  varname = expr;
    assign_match = re.match(r'\s*(\w+)\s*=\s*(.*)', line)
    if assign_match:
        lhs = assign_match.group(1)
        rhs_expr = assign_match.group(2)
        c_rhs = convert_expr(rhs_expr, state_idx, param_idx, local_vars)
        
        if lhs in state_idx:
            # Shouldn't happen outside d/dt, but handle it
            return f'  /* WARNING: direct state assignment */ y[{state_idx[lhs]}] = {c_rhs}'
        elif lhs in param_idx:
            # Parameter override (like R_art in disease effects)
            # Treat as local variable shadow
            if lhs not in local_vars:
                local_vars.add(lhs)
                return f'  double {lhs} = {c_rhs}'
            else:
                return f'  {lhs} = {c_rhs}'
        else:
            if lhs not in local_vars:
                local_vars.add(lhs)
                return f'  double {lhs} = {c_rhs}'
            else:
                return f'  {lhs} = {c_rhs}'
    
    # Fallback: convert as expression
    c_line = convert_expr(stripped, state_idx, param_idx, local_vars)
    return '  ' + c_line

def convert_expr(expr, state_idx, param_idx, local_vars):
    """Convert an R expression to C by replacing identifiers."""
    tokens = tokenize(expr)
    result = []
    
    i = 0
    while i < len(tokens):
        typ, val = tokens[i]
        
        if typ == 'COMMENT':
            # Convert R comment to C comment
            result.append('/* ' + val[1:].strip() + ' */')
            i += 1
            continue
        
        if typ == 'ID':
            # Check for R functions that need C equivalents
            if val == 'exp':
                result.append('exp')
            elif val == 'log':
                result.append('log')
            elif val == 'sqrt':
                result.append('sqrt')
            elif val == 'abs':
                result.append('fabs')
            elif val == 'floor':
                result.append('floor')
            elif val == 'max':
                result.append('fmax')
            elif val == 'min':
                result.append('fmin')
            elif val == 'log10':
                result.append('log10')
            elif val == 'sin':
                result.append('sin')
            elif val == 'cos':
                result.append('cos')
            elif val in state_idx and val not in local_vars:
                result.append(f'y[{state_idx[val]}]')
            elif val in param_idx and val not in local_vars:
                result.append(f'p[{param_idx[val]}]')
            else:
                # Local variable or C keyword
                result.append(val)
            i += 1
            continue
        
        if typ == 'OP' and val == '^':
            # R's ^ is C's pow(). Need to find the base (previous tokens) and exponent (next tokens).
            # This is tricky. Simple case: X^Y → pow(X, Y)
            # We'll handle the simple cases and flag complex ones.
            # For now, replace with pow wrapper
            # The base is everything before ^ up to the last operator or opening paren
            # This is hard to do perfectly with tokens. Let's use a simpler approach:
            # Replace A^B patterns in the string after token assembly.
            result.append('CARET_PLACEHOLDER')
            i += 1
            continue
        
        # R's ** is also pow
        if typ == 'OP' and val == '**':
            result.append('CARET_PLACEHOLDER')
            i += 1
            continue
        
        result.append(val)
        i += 1
    
    c_str = ''.join(result)
    
    # Post-process: handle ^ / ** → pow()
    # This needs regex on the assembled string
    c_str = fix_power_operator(c_str)
    
    # Ensure semicolons at end of statements (not for if/else/braces)
    stripped = c_str.strip()
    if stripped and not stripped.endswith('{') and not stripped.endswith('}') and \
       not stripped.startswith('if') and not stripped.startswith('} else') and \
       not stripped.startswith('/*') and not stripped.endswith('*/') and \
       not stripped.endswith(';'):
        c_str = c_str.rstrip() + ';'
    
    return c_str

def fix_power_operator(s):
    """Replace X CARET_PLACEHOLDER Y with pow(X, Y)."""
    # Handle cases like: (expr)CARET_PLACEHOLDER(expr) and simple var CARET_PLACEHOLDER number
    while 'CARET_PLACEHOLDER' in s:
        idx = s.index('CARET_PLACEHOLDER')
        
        # Find base (everything to the left that's part of this sub-expression)
        base_end = idx
        base_start = find_base_start(s, idx)
        base = s[base_start:base_end].strip()
        
        # Find exponent (everything to the right)
        exp_start = idx + len('CARET_PLACEHOLDER')
        exp_end = find_exp_end(s, exp_start)
        exponent = s[exp_start:exp_end].strip()
        
        # Replace
        s = s[:base_start] + f'safe_pow({base}, {exponent})' + s[exp_end:]
    
    return s

def find_base_start(s, caret_idx):
    """Find where the base expression starts, scanning left from caret."""
    i = caret_idx - 1
    # Skip whitespace
    while i >= 0 and s[i] == ' ':
        i -= 1
    if i < 0:
        return 0
    
    if s[i] == ')':
        # Matched parenthesized expression
        depth = 1
        i -= 1
        while i >= 0 and depth > 0:
            if s[i] == ')': depth += 1
            elif s[i] == '(': depth -= 1
            i -= 1
        # i is now one before the opening paren, check if there's a function name
        while i >= 0 and (s[i].isalnum() or s[i] == '_'):
            i -= 1
        return i + 1
    elif s[i].isalnum() or s[i] == '_' or s[i] == ']':
        # Simple identifier or array access
        if s[i] == ']':
            while i >= 0 and s[i] != '[':
                i -= 1
            while i > 0 and (s[i-1].isalnum() or s[i-1] == '_'):
                i -= 1
            return i
        else:
            while i >= 0 and (s[i].isalnum() or s[i] == '_' or s[i] == '.' or s[i] == '[' or s[i] == ']'):
                i -= 1
            return i + 1
    else:
        return i + 1

def find_exp_end(s, start):
    """Find where the exponent expression ends, scanning right from after caret."""
    i = start
    # Skip whitespace
    while i < len(s) and s[i] == ' ':
        i += 1
    if i >= len(s):
        return len(s)
    
    if s[i] == '(':
        # Parenthesized expression
        depth = 1
        i += 1
        while i < len(s) and depth > 0:
            if s[i] == '(': depth += 1
            elif s[i] == ')': depth -= 1
            i += 1
        return i
    elif s[i] == '-':
        # Negative sign + number or identifier
        i += 1
        while i < len(s) and (s[i].isalnum() or s[i] == '_' or s[i] == '.'):
            i += 1
        return i
    else:
        # Simple number or identifier
        while i < len(s) and (s[i].isalnum() or s[i] == '_' or s[i] == '.' or s[i] == '[' or s[i] == ']'):
            i += 1
        return i

# =========================================================================
# Step 5: Assemble the C file
# =========================================================================
def generate_c_file(ode_lines, state_vars, param_names, output_file):
    n_state = len(state_vars)
    n_param = len(param_names)
    
    state_idx = {name: i for i, name in enumerate(state_vars)}
    param_idx = {name: i for i, name in enumerate(param_names)}
    
    # Pre-process: join continuation lines (lines ending with operator or comma)
    joined_lines = []
    i = 0
    while i < len(ode_lines):
        line = ode_lines[i]
        # Strip trailing comment for continuation check
        code_part = line.split('#')[0] if '#' in line else line
        stripped_code = code_part.rstrip()
        
        # Check if line continues: ends with operator, comma, or the next line
        # starts with an operator
        while stripped_code and stripped_code[-1] in '+-*/,(&|' and i + 1 < len(ode_lines):
            i += 1
            next_line = ode_lines[i]
            line = line.rstrip() + '\n' + next_line
            code_part = next_line.split('#')[0] if '#' in next_line else next_line
            stripped_code = code_part.rstrip()
        
        # Also check if next line starts with an operator (continuation)
        if i + 1 < len(ode_lines):
            next_stripped = ode_lines[i+1].strip()
            if next_stripped and next_stripped[0] in '+-*/&|':
                # This line continues
                while i + 1 < len(ode_lines):
                    next_stripped = ode_lines[i+1].strip()
                    if not next_stripped or next_stripped[0] not in '+-*/&|':
                        break
                    i += 1
                    line = line.rstrip() + ' ' + ode_lines[i].strip()
        
        # Collapse multi-line to single line (remove internal newlines)
        line = ' '.join(line.split('\n'))
        joined_lines.append(line)
        i += 1
    
    # First pass: identify variables assigned inside braced blocks AND parameters
    # that are overwritten in the model (need to become local variables)
    if_else_vars = set()
    overwritten_params = set()
    brace_depth = 0
    for line in joined_lines:
        stripped = line.strip()
        brace_depth += stripped.count('{') - stripped.count('}')
        
        # Check for assignments
        m = re.match(r'\s*(\w+)\s*=\s*', stripped)
        if m and not stripped.startswith('d/dt') and not stripped.startswith('if') and not stripped.startswith('}'):
            varname = m.group(1)
            if varname in param_idx and varname not in state_idx:
                # Parameter being overwritten — needs local shadow
                overwritten_params.add(varname)
            elif varname not in state_idx and varname not in param_idx:
                if brace_depth > 0:
                    if_else_vars.add(varname)
    
    # Overwritten params also need pre-declaration (initialized from param array)
    all_predeclare = if_else_vars | overwritten_params
    
    local_vars = set()
    
    # Convert all lines
    c_body_lines = []
    
    # Pre-declare if/else variables and overwritten parameters
    if all_predeclare:
        c_body_lines.append('  /* Variables assigned in if/else blocks or shadowing parameters — pre-declared */')
        for var in sorted(all_predeclare):
            if var in overwritten_params:
                c_body_lines.append(f'  double {var} = p[{param_idx[var]}];  /* init from param, overwritten in model */')
            else:
                c_body_lines.append(f'  double {var} = 0.0;')
            local_vars.add(var)
        c_body_lines.append('')
    
    for line in joined_lines:
        c_line = convert_line_to_c(line, state_idx, param_idx, local_vars)
        if c_line is not None:
            c_body_lines.append(c_line)
    
    # Assemble
    header = f"""/*
 * Hallow et al. Cardiorenal Model — Auto-generated C code
 * 
 * Generated from modelfile_commented.R by rxode_to_c.py
 * {n_state} state variables, {n_param} parameters
 *
 * Function signature:
 *   void hallow_rhs(double t, const double *y, double *dydt, const double *p)
 *
 * Compile: gcc -O2 -shared -fPIC -o hallow_rhs.so hallow_rhs.c -lm
 */

#include <math.h>
#include <float.h>

/* Safe exp to prevent overflow */
static inline double safe_exp(double x) {{
    if (x > 500.0) return exp(500.0);
    if (x < -500.0) return 0.0;
    return exp(x);
}}

/* Safe pow: handles negative base with fractional exponent */
static inline double safe_pow(double base, double e) {{
    if (base < 0.0 && e != floor(e)) return -pow(-base, e);
    if (base == 0.0 && e <= 0.0) return 0.0;
    return pow(base, e);
}}

#define N_STATE {n_state}
#define N_PARAM {n_param}

void hallow_rhs(double t, const double *y, double *dydt, const double *p) {{
"""
    
    footer = """
}

/* Wrapper for scipy/ctypes: void rhs(int *n, double *t, double *y, double *dydt, double *p, int *ip) */
void rhs_desolve(int *n, double *t, double *y, double *dydt, double *rpar, int *ipar) {
    hallow_rhs(*t, y, dydt, rpar);
}
"""
    
    with open(output_file, 'w') as f:
        f.write(header)
        
        # Write state variable index comments
        f.write('  /* State variable indices:\n')
        for i, sv in enumerate(state_vars):
            f.write(f'   *   y[{i}] = {sv}\n')
        f.write('   */\n\n')
        
        # Write parameter index comments (first 20, then "...")
        f.write('  /* Parameter indices (see param_map.json for full list):\n')
        for i, pn in enumerate(param_names[:20]):
            f.write(f'   *   p[{i}] = {pn}\n')
        f.write(f'   *   ... ({n_param} total)\n')
        f.write('   */\n\n')
        
        # Write body
        for line in c_body_lines:
            f.write(line + '\n')
        
        f.write(footer)
    
    # Also write index maps as JSON for Python
    import json
    with open(output_file.replace('.c', '_state_map.json'), 'w') as f:
        json.dump({name: i for i, name in enumerate(state_vars)}, f, indent=2)
    with open(output_file.replace('.c', '_param_map.json'), 'w') as f:
        json.dump({name: i for i, name in enumerate(param_names)}, f, indent=2)
    
    return state_idx, param_idx

# =========================================================================
# Main
# =========================================================================
if __name__ == '__main__':
    print("Step 1: Extracting ODE string...")
    ode_lines = extract_ode_string(os.path.join(PROJECT_DIR, "modelfile_commented.R"))
    print(f"  {len(ode_lines)} lines")
    
    print("Step 2: Extracting state variables...")
    state_vars = extract_state_vars(ode_lines)
    print(f"  {len(state_vars)} state variables")
    
    print("Step 3: Extracting parameter names...")
    param_names = extract_param_names(os.path.join(PROJECT_DIR, "calcNomParams_timescale.R"))
    print(f"  {len(param_names)} parameters")
    
    print("Step 4: Generating C code...")
    output_c = "hallow_rhs.c"
    state_idx, param_idx = generate_c_file(ode_lines, state_vars, param_names, output_c)
    print(f"  Written to {output_c}")
    
    print("Step 5: Compiling...")
    compile_cmd = f"gcc -O2 -shared -fPIC -o hallow_rhs.so {output_c} -lm"
    result = subprocess.run(compile_cmd, shell=True, capture_output=True, text=True)
    if result.returncode == 0:
        print("  Compiled successfully to hallow_rhs.so")
    else:
        print("  Compilation failed:")
        print(result.stderr[:2000])
    
    print("\nDone. Files:")
    print(f"  {output_c} — C source")
    print(f"  hallow_rhs.so — Shared library")
    print(f"  hallow_rhs_state_map.json — State variable indices")
    print(f"  hallow_rhs_param_map.json — Parameter indices")
