#!/usr/bin/env python3
import argparse
import importlib
import inspect
import pathlib
import sys
import json
import numpy as np

# Make local geometries.py importable even if called from elsewhere
SCRIPT_DIR = pathlib.Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

def fmt_f(x, prec=9):
    return f"{float(x):.{prec}f}"

def fmt_vec(v, prec=9):
    return "   " + "  ".join(fmt_f(v[i], prec) for i in range(3))

def call_twisted_multilayer(mod, m0, r, rot, shift, rotate_center, ab, mirror):
    fn = getattr(mod, "twisted_multilayer", None)
    if fn is None:
        raise RuntimeError("geometries.twisted_multilayer not found")
    params = inspect.signature(fn).parameters
    kw = {"m0": m0}
    if "r" in params: kw["r"] = r
    if "rot" in params: kw["rot"] = rot
    if "shift" in params: kw["shift"] = shift
    if "rotate_center" in params: kw["rotate_center"] = rotate_center
    if "ab" in params: kw["ab"] = ab
    if "mirror" in params: kw["mirror"] = mirror
    try:
        return fn(**kw)
    except TypeError:
        return fn(m0=m0)

def call_bilayer(mod, l, mirror):
    fn = getattr(mod, "bilayer", None)
    if fn is None:
        raise RuntimeError("geometries.bilayer not found")
    try:
        return fn(l=l, mirror=mirror)
    except TypeError:
        return fn(l, mirror)

def call_supercell2d(mod, g, n1, n2):
    fn = getattr(mod, "supercell2d", None)
    if fn is None:
        raise RuntimeError("geometries.supercell2d not found")
    try:
        return fn(g, n1=n1, n2=n2)
    except TypeError:
        return fn(g, n1, n2)

def build_geometry(mode: str, m0: int | None, vacuum: float, scale_xy: float, ab_flag: int, mirror: int, 
                   n1: int, n2: int, dev_xy: float, dev_z: float, seed: int | None, 
                   target_distance: float | None, slide: list[float] | None, 
                   bilayer_stack: str):
    try:
        mod = importlib.import_module("geometries")
    except Exception as e:
        print(f"Error: could not import local module 'geometries.py'. {e}", file=sys.stderr)
        sys.exit(10)

    if mode == "twist":
        print(f"Generating twisted multilayer with m0={m0}, mirror={mirror}...")
        r, rot, shift = 1, [0, 1], [[0, 0], [0, 0]]
        gs = call_twisted_multilayer(mod, m0=m0, r=r, rot=rot, shift=shift, rotate_center="atom", ab=ab_flag, mirror=mirror)
    elif mode == "bilayer":
        print(f"Generating non-twisted bilayer with stacking '{bilayer_stack}', mirror={mirror}...")
        gs = call_bilayer(mod, l=bilayer_stack, mirror=mirror)
        
        # Specially aligned mirror constraints for bilateral structures
        if hasattr(gs, "sublattice"):
            if bilayer_stack == "AA":
                if mirror == 0:
                    print("AA Stacking (mirror=0): Globally inverting species (B <-> N)...")
                    gs.sublattice = gs.sublattice * -1
            else:
                if mirror == 1:
                    print(f"{bilayer_stack} Stacking (mirror=1): Globally inverting species (B <-> N)...")
                    gs.sublattice = gs.sublattice * -1
    else:
        raise ValueError(f"Unknown geometry mode: {mode}")

    # Process supercell expansion if requested
    if n1 > 1 or n2 > 1:
        print(f"Applying supercell2d with expansion factors n1={n1}, n2={n2}...")
        gs = call_supercell2d(mod, gs, n1=n1, n2=n2)

    if hasattr(gs, "set_origin"):
        gs.set_origin([0.5, 0, 0])

    if hasattr(gs, "get_fractional"):
        gs.get_fractional()
    elif hasattr(gs, "real2fractional"):
        gs.real2fractional()

    v1 = np.array(getattr(gs, "a1"), float) * float(scale_xy)
    v2 = np.array(getattr(gs, "a2"), float) * float(scale_xy)
    v3 = np.array([0.0, 0.0, float(vacuum)], float)

    frac = np.array(getattr(gs, "frac_r"), float)
    xfrac = frac[:, 0] - np.min(frac[:, 0])
    yfrac = frac[:, 1] - np.min(frac[:, 1])
    zfrac = frac[:, 2] * (scale_xy / float(vacuum)) + 0.5

    nat = len(frac)

    if hasattr(gs, "sublattice"):
        subl = np.array(gs.sublattice)
        if len(subl) < nat:
            reps = (nat + len(subl) - 1) // len(subl)
            subl = np.tile(subl, reps)[:nat]
    else:
        subl = np.array([1 if i % 2 == 0 else -1 for i in range(nat)])

    atoms = ['B' if int(s) == -1 else 'N' for s in subl]

    # 1. Apply interlayer sliding relative to primitive monolayer
    if slide is not None:
        s1, s2 = slide
        scale_factor = np.sqrt(nat / 4.0)
        s1_scaled = s1 / scale_factor
        s2_scaled = s2 / scale_factor
        print(f"Applying slide [{s1:.4f}, {s2:.4f}] -> scaled: [{s1_scaled:.6f}, {s2_scaled:.6f}]...")
        upper_mask = zfrac >= 0.5
        xfrac[upper_mask] += s1_scaled
        yfrac[upper_mask] += s2_scaled

    # 2. Apply Gaussian random perturbations
    if dev_xy > 0.0 or dev_z > 0.0:
        if seed is not None:
            np.random.seed(seed)
        print(f"Applying Gaussian deviations (dev_xy={dev_xy} Å, dev_z={dev_z} Å)...")

        if dev_xy > 0.0:
            dx_cart = np.random.normal(0.0, dev_xy, nat)
            dy_cart = np.random.normal(0.0, dev_xy, nat)
            M = np.array([[v1[0], v2[0]], [v1[1], v2[1]]])
            M_inv = np.linalg.inv(M)
            delta_frac = M_inv @ np.vstack([dx_cart, dy_cart])
            xfrac += delta_frac[0, :]
            yfrac += delta_frac[1, :]

        if dev_z > 0.0:
            dz_cart = np.random.normal(0.0, dev_z, nat)
            zfrac += dz_cart / float(vacuum)

    # 3. Apply interlayer distance correction
    if target_distance is not None:
        z_lower_init = zfrac[zfrac < 0.5]
        z_upper_init = zfrac[zfrac >= 0.5]
        if len(z_lower_init) > 0 and len(z_upper_init) > 0:
            d_init = (np.mean(z_upper_init) - np.mean(z_lower_init)) * float(vacuum)
            delta_d = target_distance - d_init
            delta_z = delta_d / float(vacuum)
            zfrac[zfrac < 0.5] -= delta_z / 2.0
            zfrac[zfrac >= 0.5] += delta_z / 2.0

    # 4. Eliminate negative coordinates
    min_x, min_y = np.min(xfrac), np.min(yfrac)
    if min_x < 0.0:
        xfrac -= min_x
    if min_y < 0.0:
        yfrac -= min_y

    B_lower = [i for i in range(nat) if atoms[i] == 'B' and zfrac[i] < 0.5]
    B_upper = [i for i in range(nat) if atoms[i] == 'B' and zfrac[i] >= 0.5]
    N_lower = [i for i in range(nat) if atoms[i] == 'N' and zfrac[i] < 0.5]
    N_upper = [i for i in range(nat) if atoms[i] == 'N' and zfrac[i] >= 0.5]
    ordered_indices = B_lower + B_upper + N_lower + N_upper

    cell_block = "CELL_PARAMETERS (angstrom)\n" + \
                 fmt_vec(v1, 9) + "\n" + fmt_vec(v2, 9) + "\n" + fmt_vec(v3, 9) + "\n"

    lines = ["ATOMIC_POSITIONS (crystal)"]
    for i in ordered_indices:
        sp = atoms[i]
        lines.append(f"{sp:<2}            {fmt_f(xfrac[i], 10)}        {fmt_f(yfrac[i], 10)}        {fmt_f(zfrac[i], 10)}")
    pos_block = "\n".join(lines) + "\n"
    return cell_block, pos_block, nat, v1, v2

def write_pw(prefix: str, outdir: str, pseudo_dir: str, nat: int, cell_block: str, pos_block: str,
             ecutwfc, ecutrho, vdw_corr, input_dft, assume_isolated, nosym: bool, noinv: bool,
             tprnfor: bool, tstress: bool, conv_thr, mixing_beta, startingpot,
             pp_N, mass_N, pp_B, mass_B, kgrid, kpoints_mode: str, filename: pathlib.Path):
    control = f"&CONTROL\n  calculation = 'scf',\n  prefix      = '{prefix}',\n  outdir      = '{outdir}',\n  pseudo_dir  = '{pseudo_dir}',\n"
    if tprnfor: control += "  tprnfor     = .true.,\n"
    if tstress: control += "  tstress     = .true.,\n"
    control += "/\n\n"

    system = f"&SYSTEM\n  ibrav     = 0,\n  nat  = {nat},\n  ntyp = 2,\n  ecutwfc = {ecutwfc},\n  ecutrho = {ecutrho},\n  vdw_corr = '{vdw_corr}',\n  input_dft = '{input_dft}',\n  assume_isolated = '{assume_isolated}',\n"
    if nosym: system += "  nosym = .true.,\n"
    if noinv: system += "  noinv = .true.,\n"
    system += "/\n\n"

    electrons = f"&ELECTRONS\n  conv_thr = {conv_thr},\n  mixing_beta = {mixing_beta},\n  startingpot = '{startingpot}'\n/\n\n"
    atomic_species = f"ATOMIC_SPECIES\n    N {mass_N} {pp_N}\n    B {mass_B}  {pp_B}\n\n"

    if kpoints_mode.lower() == "gamma":
        kpts = "K_POINTS gamma\n"
    else:
        kpts = f"K_POINTS automatic\n  {kgrid[0]} {kgrid[1]} {kgrid[2]} {kgrid[3]} {kgrid[4]} {kgrid[5]}\n"

    text = control + system + electrons + atomic_species + cell_block.rstrip() + "\n\n" + pos_block.rstrip() + "\n\n" + kpts
    filename.write_text(text, encoding="utf-8")

def write_ph(prefix: str, outdir: str, filename: pathlib.Path):
    text = f"&INPUTPH\n  prefix = '{prefix}',\n  outdir = '{outdir}',\n  tr2_ph = 1.0d-12,\n  fildyn = '{prefix}.dyn',\n  epsil = .true.,\n  asr = .true.,\n/\n0.0 0.0 0.0\n"
    filename.write_text(text, encoding="utf-8")

def write_dynmat(prefix: str, filename: pathlib.Path):
    text = f"&INPUT\n  fildyn = '{prefix}.dyn'\n  asr = 'crystal'\n  lperm = .true.\n/\n0.0 0.0 0.0\n"
    filename.write_text(text, encoding="utf-8")

def write_bash_runner(job_name: str, filename: pathlib.Path, cores: int, time_limit: str, mem_limit: str):
    text = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --time={time_limit}
#SBATCH --mem={mem_limit}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node={cores}
#SBATCH --cpus-per-task=1
#SBATCH --output=slurm-%x-%j.out
set -e

module restore qemod

echo "=== [1/3] Starting SCF calculations (pw.x) ==="
srun pw.x -np {cores} -inp pw.in > pw.out
echo "-> Done with PW!"

echo "=== [2/3] Starting phonon calculations (ph.x) ==="
srun ph.x -np {cores} -inp ph.in > ph.out
echo "-> Done with PH!"

echo "=== [3/3] Starting dynamical matrix calculations (dynmat.x) ==="
srun dynmat.x -in dynmat.in > dynmat.out
echo "-> Done with DYNMAT!"

echo "=== All steps finished successfully! ==="
"""
    filename.write_text(text, encoding="utf-8")
    try:
        filename.chmod(filename.stat().st_mode | 0o111)
    except Exception:
        pass

def main():
    ap = argparse.ArgumentParser(description="Produce a QE calculation folder (pw -> ph -> dynmat) with Slurm script.")
    ap.add_argument("--mode", choices=["twist", "bilayer"], default="twist")
    ap.add_argument("--m0", type=int)
    ap.add_argument("--ab", type=int, choices=[0, 1], default=0)
    ap.add_argument("--vacuum", type=float, required=True)
    ap.add_argument("--scale-xy", type=float, default=1.44)
    ap.add_argument("--mirror", type=int, default=0)
    ap.add_argument("--bilayer-stack", choices=["AA", "AB"], default=None)
    ap.add_argument("--distance", type=float, default=None)
    ap.add_argument("--slide", nargs=2, type=float, default=None)
    ap.add_argument("--n1", type=int, default=1)
    ap.add_argument("--n2", type=int, default=1)
    ap.add_argument("--dev-xy", type=float, default=0.0)
    ap.add_argument("--dev-z", type=float, default=0.0)
    ap.add_argument("--seed", type=int, default=None)

    # Slurm Settings
    ap.add_argument("--slurm-name", default=None)
    ap.add_argument("--slurm-cores", type=int, default=16)
    ap.add_argument("--slurm-time", default="00:30:00")
    ap.add_argument("--slurm-mem", default="10G")

    # QE settings
    ap.add_argument("--prefix", default="13")
    ap.add_argument("--outdir", default="/tmp")
    ap.add_argument("--pseudo-dir", default="../../pseudos/")
    ap.add_argument("--input-dft", default="PBEsol")
    ap.add_argument("--assume-isolated", default="2D")
    ap.add_argument("--ecutwfc", type=float, default=70)
    ap.add_argument("--ecutrho", type=float, default=280)
    ap.add_argument("--vdw-corr", default="grimme-d2")
    ap.add_argument("--tprnfor", type=bool, default=True)
    ap.add_argument("--tstress", type=bool, default=True)
    ap.add_argument("--nosym", type=bool, default=True)
    ap.add_argument("--noinv", type=bool, default=True)
    ap.add_argument("--conv-thr", default="1.D-10")
    ap.add_argument("--mixing-beta", type=float, default=0.1)
    ap.add_argument("--startingpot", default="file")
    ap.add_argument("--pp-N", default="N.pbe-n-kjpaw_psl.1.0.0.UPF")
    ap.add_argument("--mass-N", type=float, default=14.007)
    ap.add_argument("--pp-B", default="B.pbe-n-kjpaw_psl.0.1.UPF")
    ap.add_argument("--mass-B", type=float, default=10.81)

    # Kpoints settings
    ap.add_argument("--kpoints", choices=["automatic", "gamma"], default="automatic")
    ap.add_argument("--k-product", type=float, default=20.0)
    ap.add_argument("--kgrid", nargs=6, type=int, default=None)

    # Output Controls
    ap.add_argument("-o", "--output", type=pathlib.Path)
    ap.add_argument("--name")
    ow = ap.add_mutually_exclusive_group()
    ow.add_argument("--overwrite", dest="overwrite", action="store_true")
    ow.add_argument("--no-overwrite", dest="overwrite", action="store_false")
    ap.set_defaults(overwrite=True)

    args = ap.parse_args()

    if args.mode == "twist" and args.m0 is None:
        ap.error("--m0 is required when --mode is 'twist'.")

    bi_stack = args.bilayer_stack if args.bilayer_stack is not None else ("AB" if args.ab == 1 else "AA")

    cell_block, pos_block, nat, v1, v2 = build_geometry(
        mode=args.mode, m0=args.m0, vacuum=args.vacuum, scale_xy=args.scale_xy, ab_flag=args.ab,
        mirror=args.mirror, n1=args.n1, n2=args.n2, dev_xy=args.dev_xy, dev_z=args.dev_z, seed=args.seed,
        target_distance=args.distance, slide=args.slide, bilayer_stack=bi_stack
    )

    if args.kpoints.lower() == "automatic":
        if args.kgrid is not None:
            kgrid = tuple(args.kgrid)
        else:
            k1 = max(1, int(round(args.k_product / np.linalg.norm(v1))))
            k2 = max(1, int(round(args.k_product / np.linalg.norm(v2))))
            kgrid = (k1, k2, 1, 0, 0, 0)
    else:
        kgrid = (1, 1, 1, 0, 0, 0)

    out_dir = args.output if args.output else (pathlib.Path(args.name) if args.name else pathlib.Path(f"job_{args.prefix}"))
    if out_dir.exists() and not args.overwrite:
        print(f"Error: {out_dir} exists.", file=sys.stderr)
        sys.exit(2)

    out_dir.mkdir(parents=True, exist_ok=True)
    pw_filepath = out_dir / "pw.in"

    write_pw(
        prefix=args.prefix, outdir=args.outdir, pseudo_dir=args.pseudo_dir, nat=nat,
        cell_block=cell_block, pos_block=pos_block, ecutwfc=args.ecutwfc, ecutrho=args.ecutrho,
        vdw_corr=args.vdw_corr, input_dft=args.input_dft, assume_isolated=args.assume_isolated,
        nosym=args.nosym, noinv=args.noinv, tprnfor=args.tprnfor, tstress=args.tstress,
        conv_thr=args.conv_thr, mixing_beta=args.mixing_beta, startingpot=args.startingpot,
        pp_N=args.pp_N, mass_N=args.mass_N, pp_B=args.pp_B, mass_B=args.mass_B,
        kgrid=kgrid, kpoints_mode=args.kpoints, filename=pw_filepath
    )

    try:
        helpers = importlib.import_module("helpers")
        avg_dist = helpers.interlayer_distance(pw_filepath, nlayer=2, c=args.vacuum)
        pw_filepath.write_text(f"! Average Interlayer Distance: {avg_dist:.6f} Å\n" + pw_filepath.read_text(encoding="utf-8"), encoding="utf-8")
    except Exception:
        pass

    write_ph(prefix=args.prefix, outdir=args.outdir, filename=out_dir / "ph.in")
    write_dynmat(prefix=args.prefix, filename=out_dir / "dynmat.in")

    resolved_job_name = args.slurm_name if args.slurm_name is not None else args.prefix
    write_bash_runner(
        job_name=resolved_job_name, filename=out_dir / "run.sh",
        cores=args.slurm_cores, time_limit=args.slurm_time, mem_limit=args.slurm_mem
    )

    # Dump metadata configuration file
    try:
        meta_dict = {k: str(v) if isinstance(v, pathlib.Path) else v for k, v in vars(args).items()}
        (out_dir / "createscf_meta.json").write_text(json.dumps(meta_dict, indent=2), encoding="utf-8")
    except Exception:
        pass

    print(f"Generated workflow inside: {out_dir.resolve()}")

if __name__ == "__main__":
    main()
