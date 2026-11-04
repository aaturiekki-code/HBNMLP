#!/usr/bin/env python3
import argparse
import json
import os
import pathlib
import re
import subprocess
import sys

def parse_pw_out(filepath):
    energy = None
    if not os.path.exists(filepath):
        return energy
    with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
        content = f.read()
    energy_match = re.search(r'!\s+total energy\s+=\s+([\d\.\-]+)\s+Ry', content)
    if energy_match:
        energy = float(energy_match.group(1))
    return energy

def parse_dielectric_tensor(content, marker_regex):
    match = re.search(marker_regex, content, re.IGNORECASE)
    if not match:
        return None
    start_idx = match.end()
    lines = content[start_idx:].split('\n')
    tensor = []
    float_or_asterisk_pattern = re.compile(r'([-+]?\d+\.\d+(?:[eE][-+]?\d+)?|[-+]?\d+|\*+)')
    for line in lines:
        if len(tensor) == 3:
            break
        parts = float_or_asterisk_pattern.findall(line)
        if len(parts) >= 3:
            row = []
            for p in parts[:3]:
                if '*' in p:
                    row.append(float('nan'))
                else:
                    try:
                        row.append(float(p))
                    except ValueError:
                        row.append(float('nan'))
            tensor.append(row)
    if len(tensor) == 3:
        return tensor
    return None

def parse_dynmat_out(filepath):
    elec_dielectric, total_dielectric = None, None
    if not os.path.exists(filepath):
        return elec_dielectric, total_dielectric
    with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
        content = f.read()
    elec_dielectric = parse_dielectric_tensor(content, r'Electronic dielectric permittivity tensor')
    total_dielectric = parse_dielectric_tensor(content, r'with zone-center polar mode contributions')
    return elec_dielectric, total_dielectric

def extract_vacuum_from_pw_in(pw_in_path):
    if not pw_in_path.exists():
        return None
    content = pw_in_path.read_text(encoding='utf-8', errors='ignore')
    match = re.search(
        r'CELL_PARAMETERS\s*\(angstrom\)\s*\n\s*[\d\.\-]+\s+[\d\.\-]+\s+[\d\.\-]+\s*\n\s*[\d\.\-]+\s+[\d\.\-]+\s+[\d\.\-]+\s*\n\s*[\d\.\-]+\s+[\d\.\-]+\s+([\d\.\-]+)', 
        content, re.IGNORECASE
    )
    return float(match.group(1)) if match else None

def get_or_build_job_name(job_dir, parent_dir, meta):
    slurm_name = meta.get('slurm_name')
    prefix = meta.get('prefix', '13')
    
    # Clean out generic values to force-build the unique naming scheme
    if slurm_name and slurm_name != '13' and slurm_name != prefix:
        return slurm_name
    
    job_name_parts = [parent_dir.name]
    dist_match = re.search(r'dist_([\d\.]+)', job_dir.name)
    if dist_match:
        job_name_parts.append(f"d{float(dist_match.group(1)):.2f}")
    scale_match = re.search(r'scale_([\d\.]+)', job_dir.name)
    if scale_match:
        job_name_parts.append(f"s{float(scale_match.group(1)):.2f}")
    slide_match = re.search(r'slide_([\d\.]+)_([\d\.]+)', job_dir.name)
    if slide_match:
        job_name_parts.append(f"sl{float(slide_match.group(1)):.2f}_{float(slide_match.group(2)):.2f}")
    
    if len(job_name_parts) == 1:
        job_name_parts.append(job_dir.name)
    return "_".join(job_name_parts)

def rebuild_createscf_args(meta):
    cmd_args = []
    for key, val in meta.items():
        if val is None or key in ['output', 'name']:
            continue
        flag = f"--{key.replace('_', '-')}"
        if key == 'overwrite':
            cmd_args.append('--overwrite' if val else '--no-overwrite')
            continue
        if isinstance(val, bool):
            cmd_args.append(flag)
            cmd_args.append(str(val))
            continue
        if isinstance(val, list):
            cmd_args.append(flag)
            cmd_args.extend(map(str, val))
            continue
        cmd_args.append(flag)
        cmd_args.append(str(val))
    return cmd_args

def modify_run_sh_time(run_sh_path, new_time):
    if not run_sh_path.exists():
        return False
    content = run_sh_path.read_text(encoding="utf-8")
    modified_content, count = re.subn(r'(#SBATCH\s+--time\s*=\s*)([^\s]+)', f'\\g<1>{new_time}', content)
    if count > 0:
        run_sh_path.write_text(modified_content, encoding="utf-8")
        return True
    return False

def main():
    ap = argparse.ArgumentParser(description="Find unfinished sweep jobs and regenerate/resubmit them.")
    ap.add_argument("-d", "--dir", default=".", type=pathlib.Path)
    ap.add_argument("--run", action="store_true")
    ap.add_argument("--time", type=str, default=None)
    ap.add_argument("--dev-xy", type=float, default=None)
    ap.add_argument("--dev-z", type=float, default=None)

    # Legacies parameters reconstruction fallbacks
    ap.add_argument("--fallback-mode", choices=["twist", "bilayer"], default="bilayer")
    ap.add_argument("--fallback-m0", type=int, default=None)
    ap.add_argument("--fallback-n1", type=int, default=1)
    ap.add_argument("--fallback-n2", type=int, default=1)
    ap.add_argument("--fallback-mirror", type=int, default=0)
    ap.add_argument("--fallback-ab", type=int, default=0)
    ap.add_argument("--fallback-bilayer-stack", choices=["AA", "AB"], default=None)

    args = ap.parse_args()
    parent_dir = args.dir.resolve()

    script_dir = pathlib.Path(__file__).resolve().parent
    createscf_path = script_dir / "createscf.py"

    subdirs = sorted([d for d in parent_dir.iterdir() if d.is_dir() and (d / "pw.in").exists()])
    if not subdirs:
        print("No calculation folders found.")
        return

    unfinished_jobs = []
    for subdir in subdirs:
        pw_file = subdir / "pw.out"
        dyn_file = subdir / "dynmat.out"
        if not (subdir / "run.sh").exists():
            continue

        energy = parse_pw_out(pw_file)
        elec_d, total_d = parse_dynmat_out(dyn_file)

        reasons = []
        if energy is None: reasons.append("SCF failed")
        if elec_d is None: reasons.append("Electronic tensor missing")
        if total_d is None: reasons.append("Polar total tensor missing")

        if reasons:
            unfinished_jobs.append({'path': subdir, 'name': subdir.name, 'reasons': reasons})

    if not unfinished_jobs:
        print("All jobs finished successfully!")
        return

    print(f"Found {len(unfinished_jobs)} incomplete jobs.")
    if not args.run:
        print("Dry run completed. Run with --run to apply changes and submit.")
        return

    for job in unfinished_jobs:
        job_dir = job['path']
        run_sh = job_dir / "run.sh"
        pw_in = job_dir / "pw.in"
        meta_file = job_dir / "createscf_meta.json"
        
        meta = None
        if meta_file.exists():
            try:
                meta = json.loads(meta_file.read_text(encoding="utf-8"))
            except Exception:
                pass

        if meta is None:
            # Fallback legacy reconstruction
            meta = {
                'mode': args.fallback_mode, 'm0': args.fallback_m0, 'n1': args.fallback_n1, 'n2': args.fallback_n2,
                'mirror': args.fallback_mirror, 'ab': args.fallback_ab, 'bilayer_stack': args.fallback_bilayer_stack,
                'dev_xy': 0.0, 'dev_z': 0.0
            }
            v_val = extract_vacuum_from_pw_in(pw_in)
            meta['vacuum'] = v_val if v_val is not None else 30.0

            dist_match = re.search(r'dist_([\d\.]+)', job_dir.name)
            if dist_match: meta['distance'] = float(dist_match.group(1))
            scale_match = re.search(r'scale_([\d\.]+)', job_dir.name)
            if scale_match: meta['scale_xy'] = float(scale_match.group(1))
            slide_match = re.search(r'slide_([\d\.]+)_([\d\.]+)', job_dir.name)
            if slide_match: meta['slide'] = [float(slide_match.group(1)), float(slide_match.group(2))]

        # Secure Slurm Name preservation
        meta['slurm_name'] = get_or_build_job_name(job_dir, parent_dir, meta)

        if args.dev_xy is not None: meta['dev_xy'] = args.dev_xy
        if args.dev_z is not None: meta['dev_z'] = args.dev_z
        if args.time: meta['slurm_time'] = args.time

        if meta.get('seed') is not None:
            meta['seed'] = int(meta['seed']) + 1

        try:
            rebuilt_args = rebuild_createscf_args(meta)
            cmd = [sys.executable, str(createscf_path), "-o", str(job_dir)] + rebuilt_args
            subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
        except Exception as e:
            if args.time: modify_run_sh_time(run_sh, args.time)

        # Clear old output logs
        for out_name in ["pw.out", "ph.out", "dynmat.out"]:
            out_file = job_dir / out_name
            if out_file.exists():
                try: out_file.unlink()
                except Exception: pass

        # Submit run script
        try:
            subprocess.run(["sbatch", "run.sh"], cwd=str(job_dir), stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
            print(f" -> Successfully queued: {job['name']}")
        except Exception as e:
            print(f"Failed to submit: {job['name']}", file=sys.stderr)

if __name__ == "__main__":
    main()
