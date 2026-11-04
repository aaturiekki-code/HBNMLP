#!/usr/bin/env python3
"""
merge_all_to_one.py

Recursively scan a parent directory (default: ../runs_interlayer) for per-task
NPZ/CSV files and merge EVERYTHING into a single combined NPZ and a single
combined CSV.

- Finds NPZs by glob: <parent>/**/task_*/features_labels_task_*.npz (recursive)
- Merges all NPZs into one big NPZ at --out-npz (default under parent dir)
- Similarly, merges all CSVs into one big CSV at --out-csv
- Validates SOAP parameters, aggregate, stress_unit, and feature dimension across all inputs

Usage (from bashscripts/):
  mamba run -n coeff-env python3 ../pyscripts/merge_all_to_one.py \
    --parent-dir ../runs_interlayer \
    --out-npz ../runs_interlayer/combined_features_labels_all.npz \
    --out-csv ../runs_interlayer/combined_features_labels_all.csv

Options:
  --include / --exclude: globs applied to the sweep folder basename (two levels up from task_*).
  --overwrite: allow overwriting existing outputs
  --dry-run: print what would be merged without writing

No SLURM env variables required.
"""

import os
import glob
import json
import csv
import argparse
import fnmatch
import numpy as np


def find_all_task_npz(parent_dir):
    # Recursive search for per-task NPZs
    pattern = os.path.join(parent_dir, "**", "task_*", "features_labels_task_*.npz")
    return sorted(glob.glob(pattern, recursive=True))


def find_all_task_csv(parent_dir):
    pattern = os.path.join(parent_dir, "**", "task_*", "features_labels_task_*.csv")
    return sorted(glob.glob(pattern, recursive=True))


def sweep_name_for_npz(npz_path):
    # Sweep dir is two levels up: .../<sweep>/task_*/file
    task_dir = os.path.dirname(npz_path)
    sweep_dir = os.path.dirname(task_dir)
    return os.path.basename(sweep_dir)


def passes_filters(name, include_patterns, exclude_patterns):
    if include_patterns:
        if not any(fnmatch.fnmatch(name, pat) for pat in include_patterns):
            return False
    if exclude_patterns:
        if any(fnmatch.fnmatch(name, pat) for pat in exclude_patterns):
            return False
    return True


def load_task_npz(npz_path):
    with np.load(npz_path, allow_pickle=True) as z:
        X = z["X"]
        energies = z["energies"]
        stresses = z["stresses"]
        dielectrics = z["dielectrics"]
        ids = z["ids"]
        soap_params = json.loads(z["soap_params"][0])
        aggregate = str(z["aggregate"])
        stress_unit = str(z["stress_unit"])
    return X, energies, stresses, dielectrics, ids, soap_params, aggregate, stress_unit


def header_from_first_csv(csv_path):
    with open(csv_path, newline="") as f:
        r = csv.reader(f)
        header = next(r)
    return header


def merge_all_npz(npz_paths, out_npz, dry_run=False):
    if not npz_paths:
        raise FileNotFoundError("No per-task NPZ files found under the parent directory.")

    Xs, Es, Ss, Ds, IDs = [], [], [], [], []
    soap_params_ref = None
    aggregate_ref = None
    stress_unit_ref = None
    feat_dim_ref = None

    for p in npz_paths:
        X, E, S, D, I, sp, agg, su = load_task_npz(p)
        if soap_params_ref is None:
            soap_params_ref = sp
            aggregate_ref = agg
            stress_unit_ref = su
            feat_dim_ref = X.shape[1]
        else:
            # Validate consistency
            if (sp.get("r_cut") != soap_params_ref.get("r_cut") or
                sp.get("n_max") != soap_params_ref.get("n_max") or
                sp.get("l_max") != soap_params_ref.get("l_max") or
                abs(sp.get("sigma", 0.0) - soap_params_ref.get("sigma", 0.0)) > 1e-12 or
                bool(sp.get("periodic")) != bool(soap_params_ref.get("periodic"))):
                raise ValueError(f"SOAP parameters mismatch (e.g., {p}).")
            if agg != aggregate_ref:
                raise ValueError(f"Aggregate mismatch: {agg} vs {aggregate_ref} in {p}")
            if su != stress_unit_ref:
                raise ValueError(f"Stress unit mismatch: {su} vs {stress_unit_ref} in {p}")
            if X.shape[1] != feat_dim_ref:
                raise ValueError(f"Feature dimension mismatch in {p}: {X.shape[1]} vs {feat_dim_ref}")

        Xs.append(X); Es.append(E); Ss.append(S); Ds.append(D); IDs.append(I)

    X = np.vstack(Xs)
    energies = np.concatenate(Es)
    stresses = np.vstack(Ss).reshape((-1, 3, 3))
    dielectrics = np.vstack(Ds).reshape((-1, 3, 3))
    ids = np.concatenate(IDs)

    if dry_run:
        return len(ids), X.shape, soap_params_ref, aggregate_ref, stress_unit_ref

    os.makedirs(os.path.dirname(out_npz) or ".", exist_ok=True)
    np.savez_compressed(
        out_npz,
        X=X,
        energies=energies,
        stresses=stresses,
        dielectrics=dielectrics,
        ids=ids,
        soap_params=np.array([json.dumps(soap_params_ref)], dtype=object),
        aggregate=aggregate_ref,
        stress_unit=stress_unit_ref,
    )
    return len(ids), X.shape, soap_params_ref, aggregate_ref, stress_unit_ref


def merge_all_csv(csv_paths, out_csv, dry_run=False):
    if not csv_paths:
        raise FileNotFoundError("No per-task CSV files found under the parent directory.")
    header = header_from_first_csv(csv_paths[0])
    if dry_run:
        return header
    os.makedirs(os.path.dirname(out_csv) or ".", exist_ok=True)
    with open(out_csv, "w", newline="") as fout:
        w = csv.writer(fout)
        w.writerow(header)
        for p in csv_paths:
            with open(p, newline="") as fin:
                r = csv.reader(fin)
                hdr = next(r, None)
                if hdr != header:
                    raise ValueError(f"CSV header mismatch in {p}")
                for row in r:
                    w.writerow(row)
    return header


def main():
    ap = argparse.ArgumentParser(description="Merge ALL per-task NPZ/CSV into ONE combined NPZ and ONE combined CSV.")
    ap.add_argument("--parent-dir", default="../runs_interlayer", help="Parent directory to scan recursively.")
    ap.add_argument("--out-npz", default=None, help="Combined NPZ path (default: <parent-dir>/combined_features_labels_all.npz)")
    ap.add_argument("--out-csv", default=None, help="Combined CSV path (default: <parent-dir>/combined_features_labels_all.csv)")
    ap.add_argument("--include", nargs="*", default=None, help="Include sweeps: globs on sweep folder basename.")
    ap.add_argument("--exclude", nargs="*", default=None, help="Exclude sweeps: globs on sweep folder basename.")
    ap.add_argument("--overwrite", action="store_true", help="Overwrite existing outputs.")
    ap.add_argument("--dry-run", action="store_true", help="Only report what would be merged; do not write.")
    args = ap.parse_args()

    parent = os.path.abspath(args.parent_dir)
    if not os.path.isdir(parent):
        raise FileNotFoundError(f"Parent directory not found: {parent}")

    out_npz = args.out_npz or os.path.join(parent, "combined_features_labels_all.npz")
    out_csv = args.out_csv or os.path.join(parent, "combined_features_labels_all.csv")

    if not args.overwrite and (os.path.isfile(out_npz) or os.path.isfile(out_csv)):
        raise FileExistsError(f"Output file(s) already exist. Use --overwrite to overwrite:\n  {out_npz}\n  {out_csv}")

    # Discover NPZs and optionally filter by sweep name
    all_npz = find_all_task_npz(parent)
    if not all_npz:
        print(f"No per-task NPZs found under: {parent}")
        return

    if args.include or args.exclude:
        filtered = []
        for p in all_npz:
            name = sweep_name_for_npz(p)
            if passes_filters(name, args.include, args.exclude):
                filtered.append(p)
        all_npz = filtered
        if not all_npz:
            print(f"No per-task NPZs matched filters under: {parent}")
            return

    # Discover CSVs similarly (for the global CSV)
    all_csv = find_all_task_csv(parent)
    if args.include or args.exclude:
        filtered_csv = []
        for p in all_csv:
            name = sweep_name_for_npz(p)  # same function works for csv path
            if passes_filters(name, args.include, args.exclude):
                filtered_csv.append(p)
        all_csv = filtered_csv

    print("Merging across ALL sweeps found under:", parent)
    print("  NPZ files:", len(all_npz))
    print("  CSV files:", len(all_csv))
    print("  Output NPZ:", out_npz)
    print("  Output CSV:", out_csv)
    if args.dry_run:
        print("Dry-run mode: no files will be written.")

    # Merge NPZs (global)
    n, xshape, sp, agg, su = merge_all_npz(all_npz, out_npz, dry_run=args.dry_run)
    # Merge CSVs (global)
    header = merge_all_csv(all_csv, out_csv, dry_run=args.dry_run) if all_csv else None

    print("==== Merge ALL summary ====")
    print(f"Samples total: {n}")
    print(f"Feature dim:   {xshape[1]}")
    print(f"Aggregate:     {agg}")
    print(f"Stress unit:   {su}")
    if header:
        print(f"CSV header columns: {len(header)}")
    print("Done.")

if __name__ == "__main__":
    main()
