/*
 * hallow_integrator.c — C-level LSODA integrator for the Hallow cardiorenal model.
 *
 * Uses liblsoda (sdwfrost/liblsoda, MIT license) — a C translation of the
 * Fortran LSODA algorithm. All RHS calls stay in C, eliminating the
 * Python↔C callback overhead that dominates scipy-based integration.
 */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "lsoda/lsoda.h"

#define NS 70   /* number of state variables */

/* Declared in hallow_rhs.c */
extern void hallow_rhs(double t, const double *y, double *dydt, const double *p);

/* ===================================================================
 * RHS adapter: bridges hallow_rhs to liblsoda's callback signature.
 *
 * liblsoda internally does y-- before calling the function with y+1,
 * so the callback receives a standard 0-indexed pointer.
 * =================================================================== */
static int lsoda_rhs_adapter(double t, double *y, double *ydot, void *data) {
    const double *params = (const double *)data;
    hallow_rhs(t, y, ydot, params);
    for (int i = 0; i < NS; i++) {
        if (ydot[i] != ydot[i]) return -1;  /* NaN check */
    }
    return 0;
}

/* ===================================================================
 * Main integration function.
 *
 * Returns: positive = number of output points stored (including t=0).
 *          Negative = error code (-1 solver failure, -2 NaN).
 * =================================================================== */
int integrate_hallow_lsoda(
    const double *y0,       /* initial state [NS], 0-indexed */
    const double *params,   /* parameters [430] */
    double t_end,           /* final time (hours) */
    double dt_output,       /* output interval (hours) */
    double *out_t,          /* output times buffer [max_output] */
    double *out_y,          /* output states buffer [max_output * NS] */
    int max_output,         /* size of output buffer */
    int *n_rhs_calls,       /* total RHS evaluations (for diagnostics) */
    double atol_val,        /* absolute tolerance */
    double rtol_val,        /* relative tolerance */
    double max_step         /* maximum step size */
) {
    /* liblsoda expects 0-indexed y (it does y-- internally) */
    double *y = (double *)malloc(NS * sizeof(double));
    if (!y) return -1;
    memcpy(y, y0, NS * sizeof(double));

    /* Store initial state */
    out_t[0] = 0.0;
    memcpy(out_y, y0, NS * sizeof(double));
    int out_idx = 1;

    /*
     * Per-component tolerances for liblsoda.
     * liblsoda accesses rtol/atol as 1-indexed via (opt->rtol - 1)[i],
     * so opt->rtol[0] maps to index 1, opt->rtol[NS-1] maps to index NS.
     * We allocate NS elements (0-indexed), and liblsoda's -1 offset
     * makes rtol_arr[0] = rtol[1], etc.
     */
    double *rtol_arr = (double *)malloc(NS * sizeof(double));
    double *atol_arr = (double *)malloc(NS * sizeof(double));
    if (!rtol_arr || !atol_arr) {
        free(y); free(rtol_arr); free(atol_arr);
        return -1;
    }
    for (int i = 0; i < NS; i++) {
        rtol_arr[i] = rtol_val;
        atol_arr[i] = atol_val;
    }

    /* Set up liblsoda context */
    struct lsoda_context_t ctx;
    memset(&ctx, 0, sizeof(ctx));
    ctx.function = lsoda_rhs_adapter;
    ctx.data = (void *)params;
    ctx.neq = NS;
    ctx.state = 1;

    struct lsoda_opt_t opt;
    memset(&opt, 0, sizeof(opt));
    opt.ixpr = 0;
    opt.mxstep = 500000;
    opt.mxhnil = 10;
    opt.hmxi = (max_step > 0.0) ? (1.0 / max_step) : 0.0;
    opt.itask = 1;
    opt.rtol = rtol_arr;
    opt.atol = atol_arr;

    if (lsoda_prepare(&ctx, &opt) != 1) {
        fprintf(stderr, "[C LSODA] Failed to prepare context\n");
        free(y); free(rtol_arr); free(atol_arr);
        return -1;
    }

    /* Integration loop */
    double t = 0.0;
    double tout = dt_output;
    int ret = 0;

    while (tout <= t_end + dt_output * 0.5 && out_idx < max_output) {
        double target = (tout > t_end) ? t_end : tout;

        lsoda(&ctx, y, &t, target);

        if (ctx.state <= 0) {
            fprintf(stderr, "[C LSODA] Solver failure at t=%.6f: state=%d",
                    t, ctx.state);
            if (ctx.error) fprintf(stderr, " (%s)", ctx.error);
            fprintf(stderr, "\n");
            ret = -1;
            break;
        }

        /* Check for NaN */
        int has_nan = 0;
        for (int i = 0; i < NS; i++) {
            if (y[i] != y[i]) { has_nan = 1; break; }
        }
        if (has_nan) {
            fprintf(stderr, "[C LSODA] NaN detected at t=%.6f\n", t);
            ret = -2;
            break;
        }

        /* Store output point */
        out_t[out_idx] = t;
        memcpy(out_y + out_idx * NS, y, NS * sizeof(double));
        out_idx++;
        tout += dt_output;
    }

    *n_rhs_calls = 0;

    lsoda_free(&ctx);
    free(y);
    free(rtol_arr);
    free(atol_arr);

    return (ret < 0) ? ret : out_idx;
}
