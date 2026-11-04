#!/usr/bin/env python3
import argparse
import pathlib
import subprocess
import sys
import numpy as np

def parse_range(range_args, default_val=None):
    if range_args is None:
        return [default_val] if default_val is not None else [None]
    if len(range_args) != 3:
        raise ValueError("Each range must contain exactly 3 arguments: [start] [end] [points]")
    start, end, count = range_args[0], range_args[1], int(range_args[2])
    if count <= 0:
        raise ValueError("The amount of points must be >= 1")
    if count == 1:
        return [start]
    return np.linspace(start, end, count).tolist()

def filter_overridden_args(arg_list, flags_to_remove):
    filtered = []
    skip_count = 0
    for i, item in enumerate(arg_list):
        if skip_count > 0:
            skip_count -= 1
            continue
        if item in flags_to_remove:
            skip_count = 2 if item == '--slide' else 1
            continue
        filtered.append(item)
    return filtered

def write_master_submit_script(parent_dir: pathlib.Path):
    submit_script_path = parent_dir / "submit_all.sh"
    parent_name = parent_dir.name
    
    text = f"""#!/bin/bash
#SBATCH --job-name={parent_name}_submit_all
#SBATCH --time=00:30:00
#SBATCH --mem=1G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --output=slurm-%x-%j.out

if [ -n "$SLURM_SUBMIT_DIR" ]; then
    PARENT_DIR="$SLURM_SUBMIT_DIR"
else
    PARENT_DIR="$( cd "$( dirname "${{BASH_SOURCE[0]}}" )" && pwd )"
fi

cd "$PARENT_DIR"
shopt -s nullglob

for job_dir in */; do
    job_dir="${{job_dir%/}}"
    if [ -f "$PARENT_DIR/$job_dir/run.sh" ]; then
        echo "Submitting: $job_dir ..."
        (cd "$PARENT_DIR/$job_dir" && sbatch run.sh)
    fi
done
"""
    submit_script_path.write_text(text, encoding="utf-8")
    try:
        submit_script_path.chmod(submit_script_path.stat().st_mode | 0o111)
    except Exception:
        pass

def main():
    ap = argparse.ArgumentParser(description="Sweeps structures over ranges of configurations.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--parent-dir", required=True, type=pathlib.Path)
    ap.add_argument("--distance-range", nargs=3, type=float, default=None)
    ap.add_argument("--scale-range", nargs=3, type=float, default=None)
    ap.add_argument("--slide-x-range", nargs=3, type=float, default=None)
    ap.add_argument("--slide-y-range", nargs=3, type=float, default=None)

    args, forwarded_args = ap.parse_known_args()

    script_dir = pathlib.Path(__file__).resolve().parent
    createscf_path = script_dir / "createscf.py"
    if not createscf_path.exists():
        print(f"Error: createscf.py not found in {script_dir}", file=sys.stderr)
        sys.exit(1)

    forwarded_args = filter_overridden_args(forwarded_args, ['--distance', '--scale-xy', '--slide', '-o', '--output', '--name', '--slurm-name'])

    try:
        distances = parse_range(args.distance_range, default_val=None)
        scales = parse_range(args.scale_range, default_val=1.44)
        slides_x = parse_range(args.slide_x_range, default_val=None)
        slides_y = parse_range(args.slide_y_range, default_val=None)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(2)

    slide_pairs = []
    if slides_x[0] is not None and slides_y[0] is not None:
        for sx in slides_x:
            for sy in slides_y:
                slide_pairs.append((sx, sy))
    elif slides_x[0] is not None:
        for sx in slides_x:
            slide_pairs.append((sx, 0.0))
    elif slides_y[0] is not None:
        for sy in slides_y:
            slide_pairs.append((0.0, sy))
    else:
        slide_pairs.append((None, None))

    args.parent_dir.mkdir(parents=True, exist_ok=True)
    parent_folder_name = args.parent_dir.name

    for dist in distances:
        for scale in scales:
            for sx, sy in slide_pairs:
                folder_parts = []
                job_name_parts = [parent_folder_name]

                if dist is not None:
                    folder_parts.append(f"dist_{dist:.3f}")
                    job_name_parts.append(f"d{dist:.2f}")
                folder_parts.append(f"scale_{scale:.3f}")
                job_name_parts.append(f"s{scale:.2f}")
                if sx is not None and sy is not None:
                    folder_parts.append(f"slide_{sx:.3f}_{sy:.3f}")
                    job_name_parts.append(f"sl{sx:.2f}_{sy:.2f}")

                subfolder_name = "_".join(folder_parts)
                unique_job_name = "_".join(job_name_parts)
                job_dir = args.parent_dir / subfolder_name

                cmd = [sys.executable, str(createscf_path), "-o", str(job_dir), "--scale-xy", f"{scale:.6f}", "--slurm-name", unique_job_name]
                if dist is not None:
                    cmd += ["--distance", f"{dist:.6f}"]
                if sx is not None and sy is not None:
                    cmd += ["--slide", f"{sx:.6f}", f"{sy:.6f}"]
                cmd += forwarded_args

                try:
                    subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
                except subprocess.CalledProcessError as e:
                    print(f"FAILED generation for {subfolder_name}!\n{e.stderr}", file=sys.stderr)
                    sys.exit(e.returncode)

    write_master_submit_script(args.parent_dir)
    print(f"Sweep layouts complete. Launch with 'sbatch submit_all.sh' inside: {args.parent_dir.resolve()}")

if __name__ == "__main__":
    main()
