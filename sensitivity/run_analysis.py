"""
CLI entry point for sensitivity analysis.

Usage:
    # Morris only (fast screening)
    python -m sensitivity.run_analysis --method morris --r 5 --t_hours 24

    # Full analysis: Morris then Sobol on top params
    python -m sensitivity.run_analysis --method full --r 20 --t_hours 168

    # Sobol only (requires prior Morris results)
    python -m sensitivity.run_analysis --method sobol --morris_results sensitivity/results/morris.json
"""

import argparse
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def run_morris_analysis(args):
    from sensitivity.morris import run_morris
    print(f"Running Morris screening: r={args.r}, t={args.t_hours}h")
    result = run_morris(
        r=args.r,
        t_hours=args.t_hours,
        n_workers=args.n_workers,
        levels=args.levels,
        seed=args.seed,
    )
    return result


def run_sobol_analysis(args, morris_result=None, morris_path=None):
    from sensitivity.sobol import run_sobol

    # Load top params from Morris results
    if morris_result is not None:
        top_param_names = set(morris_result['top_params'][:args.top_k])
        perturbable = morris_result['perturbable']
    elif morris_path is not None:
        with open(morris_path) as f:
            morris_data = json.load(f)
        top_param_names = set(morris_data['top_params'][:args.top_k])
        perturbable = morris_data['perturbable']
    else:
        raise ValueError("Need either morris_result or morris_path for Sobol")

    top_params = [p for p in perturbable if p['name'] in top_param_names]
    print(f"Running Sobol on {len(top_params)} top parameters")

    result = run_sobol(
        top_params=top_params,
        n_samples=args.sobol_n,
        t_hours=args.t_hours,
        n_workers=args.n_workers,
        seed=args.seed,
    )
    return result


def main():
    parser = argparse.ArgumentParser(description='Sensitivity analysis for Hallow model')
    parser.add_argument('--method', choices=['morris', 'sobol', 'full'],
                        default='morris', help='Analysis method')
    parser.add_argument('--r', type=int, default=20,
                        help='Morris replicates (default: 20)')
    parser.add_argument('--levels', type=int, default=4,
                        help='Morris grid levels (default: 4)')
    parser.add_argument('--t_hours', type=float, default=168,
                        help='Simulation duration in hours (default: 168)')
    parser.add_argument('--n_workers', type=int, default=None,
                        help='Parallel workers (default: CPU count)')
    parser.add_argument('--seed', type=int, default=42,
                        help='Random seed')
    parser.add_argument('--top_k', type=int, default=50,
                        help='Top-k params for Sobol (default: 50)')
    parser.add_argument('--sobol_n', type=int, default=1024,
                        help='Sobol base sample size (default: 1024)')
    parser.add_argument('--morris_results', type=str, default=None,
                        help='Path to Morris results JSON (for sobol-only mode)')
    parser.add_argument('--output_dir', type=str, default=None,
                        help='Output directory')
    args = parser.parse_args()

    if args.output_dir is None:
        args.output_dir = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), 'results')
    os.makedirs(args.output_dir, exist_ok=True)

    t0 = time.time()
    combined = {}

    # Morris
    if args.method in ('morris', 'full'):
        morris_result = run_morris_analysis(args)
        morris_path = os.path.join(args.output_dir, 'morris.json')
        with open(morris_path, 'w') as f:
            json.dump(morris_result, f, indent=2)
        print(f"\nMorris results saved to {morris_path}")
        combined['morris'] = morris_result['results']
        combined['top_params'] = morris_result['top_params']
        combined['morris_metadata'] = morris_result['metadata']

    # Sobol
    if args.method in ('sobol', 'full'):
        if args.method == 'sobol':
            sobol_result = run_sobol_analysis(
                args, morris_path=args.morris_results)
        else:
            sobol_result = run_sobol_analysis(args, morris_result=morris_result)
        sobol_path = os.path.join(args.output_dir, 'sobol.json')
        with open(sobol_path, 'w') as f:
            json.dump(sobol_result, f, indent=2)
        print(f"\nSobol results saved to {sobol_path}")
        combined['sobol'] = sobol_result['results']
        combined['sobol_metadata'] = sobol_result['metadata']

    # Combined results
    elapsed = time.time() - t0
    combined['total_elapsed_seconds'] = round(elapsed, 1)
    combined_path = os.path.join(args.output_dir, 'sensitivity_results.json')
    with open(combined_path, 'w') as f:
        json.dump(combined, f, indent=2)
    print(f"\nCombined results saved to {combined_path}")
    print(f"Total elapsed: {elapsed:.1f}s")

    # Print summary
    if 'morris' in combined:
        print("\n=== Top 10 Parameters by Morris mu* (MAP) ===")
        for s in combined['morris'].get('MAP', [])[:10]:
            print(f"  {s['rank']:3d}. {s['param']:45s}  mu*={s['mu_star']:.4f}  sigma={s['sigma']:.4f}")


if __name__ == '__main__':
    main()
