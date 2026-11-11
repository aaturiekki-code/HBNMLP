#!/usr/bin/env python3
"""Submit/reconfigure QE calculations, optionally regenerating geometry."""

import argparse
import datetime
import fcntl
import getpass
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import uuid

from compile_dataset import read_structure, read_forcefield_results
from qe_response import write_workflow_runner
from qe_response import (
    actual_pw_identifier,
    read_folder_tensor,
    read_meta,
    write_ph,
)


ROOT = Path(__file__).resolve().parent.parent
WORKFLOWS = ("scf-only", "electronic-only", "full-response")


def now():
    return datetime.datetime.now(datetime.timezone.utc)


def active_jobs():
    result = subprocess.run(
        [
            "squeue", "--noheader", "--user", getpass.getuser(),
            "--format=%i|%Z",
        ],
        capture_output=True, text=True, check=True,
    )
    ids, folders = set(), set()

    for line in result.stdout.splitlines():
        if "|" not in line:
            continue
        job_id, folder = line.split("|", 1)
        ids.add(job_id.strip())
        folder = folder.strip()
        if folder and folder not in {"N/A", "(null)"}:
            folders.add(Path(folder).resolve())

    return ids, folders


def is_active(folder, ids, folders):
    state = {}
    path = folder / "submission_state.json"

    if path.is_file():
        try:
            state = json.loads(path.read_text())
        except (ValueError, OSError):
            pass

    return folder in folders or str(state.get("job_id", "")) in ids


def saved_workflow(meta):
    value = meta.get("workflow")
    if value is None:
        value = (
            "electronic-only"
            if meta.get("electronic_only", False)
            else "full-response"
        )
    if value not in WORKFLOWS:
        raise ValueError(f"Invalid saved workflow: {value}")
    return value


def desired_workflow(meta, requested):
    if requested == "keep":
        return saved_workflow(meta)
    if requested == "auto":
        if meta.get("mode") == "bilayer":
            return "electronic-only"
        if meta.get("mode") == "twist":
            if meta.get("m0") == 1:
                return "electronic-only"
            if meta.get("m0") == 2:
                return "scf-only"
        raise ValueError("Auto supports untwisted and m0=1,2")
    return requested


def done(path):
    return path.is_file() and "JOB DONE." in path.read_text(errors="replace")


def classify(folder, args):
    try:
        meta = read_meta(folder)
        target = desired_workflow(meta, args.workflow)
        atoms = read_structure(folder / "pw.in")
        if not (folder / "run.sh").is_file():
            raise ValueError("Missing run.sh")
    except Exception as error:
        return "preparation-error", None, [str(error)]

    errors = []
    try:
        read_forcefield_results(folder / "pw.out", len(atoms))
    except Exception as error:
        errors.append(f"SCF: {error}")

    if target != "scf-only":
        try:
            read_folder_tensor(
                folder,
                args.dielectric_symmetry_tol,
                args.max_dielectric_component,
            )
            if target == "full-response":
                if not done(folder / "ph.out"):
                    raise ValueError("PH workflow incomplete")
                if not done(folder / "dynmat.out"):
                    raise ValueError("dynmat workflow incomplete")
        except Exception as error:
            errors.append(f"Response: {error}")

    if not errors:
        return "finished", target, []

    expected = ["pw.out"]
    if target != "scf-only":
        expected.append("ph.out")
    if target == "full-response":
        expected.append("dynmat.out")

    if not any((folder / name).exists() for name in expected):
        return "new", target, errors

    state = (
        "invalid" if all(done(folder / name) for name in expected)
        else "unfinished"
    )
    return state, target, errors


def eligible(state, args):
    if state == "preparation-error":
        return False
    if args.force:
        return True
    return state in {"new", "unfinished"} or (
        state == "invalid" and args.retry_invalid
    )


def archive(folder, names, move=False):
    paths = [
        folder / name for name in names
        if (folder / name).is_file()
    ]
    if not paths:
        return

    category = "restart_archive" if move else "workflow_archive"
    stamp = now().strftime("%Y%m%dT%H%M%SZ")
    destination = (
        folder / category / f"{stamp}_{uuid.uuid4().hex[:8]}"
    )
    destination.mkdir(parents=True)

    for path in paths:
        if move:
            shutil.move(str(path), str(destination / path.name))
        else:
            shutil.copy2(path, destination / path.name)


def directive(text, name, fallback):
    match = re.search(
        rf"(?m)^\s*#SBATCH\s+--{re.escape(name)}\s*=\s*(\S+)",
        text,
    )
    return match.group(1) if match else fallback


def prepare_workflow(folder, args):
    """Rewrite response workflow only; preserve geometry."""
    meta = read_meta(folder)
    target = desired_workflow(meta, args.workflow)
    runner = (folder / "run.sh").read_text()

    if target == saved_workflow(meta) and args.time is None:
        required = ["pw.in", "run.sh"]
        if target != "scf-only":
            required.append("ph.in")
        if target == "full-response":
            required.append("dynmat.in")

        if all((folder / name).is_file() for name in required):
            return target

    prefix = actual_pw_identifier(folder, "prefix")
    outdir = "/tmp"

    time = args.time or directive("".join(runner), "time", meta.get("slurm_time"))
    memory = directive(runner, "mem", meta.get("slurm_mem"))
    cores = directive(runner, "ntasks-per-node", meta.get("slurm_cores"))
    name = directive(runner, "job-name", meta.get("slurm_name") or folder.name)

    if any(value is None for value in (time, memory, cores)):
        raise ValueError("Missing Slurm resource settings")

    archive(
        folder,
        ("ph.in", "dynmat.in", "run.sh", "createscf_meta.json"),
    )

    if target != "scf-only":
        write_ph(
            prefix, outdir, folder / "ph.in",
            electronic_only=(target == "electronic-only"),
        )

    if target == "full-response":
        (folder / "dynmat.in").write_text(
            f"&INPUT\n"
            f"  fildyn = '{prefix}.dyn',\n"
            f"  asr = 'crystal',\n"
            f"  lperm = .true.,\n"
            f"/\n"
        )

    write_workflow_runner(
        name, folder / "run.sh", int(cores), time, memory, target
    )

    meta["workflow"] = target
    meta["electronic_only"] = target == "electronic-only"
    meta["slurm_time"] = time
    (folder / "createscf_meta.json").write_text(
        json.dumps(meta, indent=2, allow_nan=False) + "\n"
    )
    return target


def regenerate(folder, args):
    """Generate in a staging directory before replacing existing inputs."""
    meta = read_meta(folder)
    target = desired_workflow(meta, args.workflow)

    for key, value in (
        ("dev_xy", args.dev_xy),
        ("dev_z", args.dev_z),
        ("seed", args.seed),
        ("slurm_time", args.time),
    ):
        if value is not None:
            meta[key] = value

    # New saved-data location for the new geometry.
    scratch = folder / "qe_scratch" / uuid.uuid4().hex[:12]
    scratch.mkdir(parents=True)
    meta["outdir"] = "/tmp"
    meta["startingpot"] = "atomic"
    meta["slurm_name"] = meta.get("slurm_name") or folder.name

    allowed = {
        "mode", "m0", "ab", "vacuum", "scale_xy", "mirror",
        "bilayer_stack", "distance", "slide", "n1", "n2",
        "dev_xy", "dev_z", "seed",
        "slurm_name", "slurm_cores", "slurm_time", "slurm_mem",
        "prefix", "outdir", "pseudo_dir", "input_dft",
        "assume_isolated", "ecutwfc", "ecutrho", "vdw_corr",
        "tprnfor", "tstress", "nosym", "noinv",
        "conv_thr", "mixing_beta", "startingpot",
        "pp_N", "mass_N", "pp_B", "mass_B",
        "kpoints", "k_product", "kgrid",
    }

    generator = Path(__file__).resolve().parent / "createscf.py"

    with tempfile.TemporaryDirectory(
        prefix=".qe_regenerate_", dir=folder.parent
    ) as temporary:
        staging = Path(temporary) / "generated"

        command = [
            sys.executable, str(generator),
            "-o", str(staging),
            "--no-overwrite",
            f"--{target}",
        ]

        for key in sorted(allowed):
            value = meta.get(key)
            if value is None:
                continue

            flag = "--" + key.replace("_", "-")

            if isinstance(value, bool):
                if not value:
                    raise ValueError(
                        f"{key}=False cannot safely be reconstructed with "
                        "the original generator's type=bool parser."
                    )
                command += [flag, "True"]
            elif isinstance(value, list):
                command += [flag, *map(str, value)]
            else:
                command += [flag, str(value)]

        result = subprocess.run(
            command, capture_output=True, text=True
        )
        if result.returncode:
            raise RuntimeError(
                f"Generation failed:\n{result.stdout}\n{result.stderr}"
            )

        read_structure(staging / "pw.in")
        if not (staging / "run.sh").is_file():
            raise ValueError("Generated run.sh missing")

        generated_meta = json.loads(
            (staging / "createscf_meta.json").read_text()
        )
        generated_meta["output"] = str(folder)
        (staging / "createscf_meta.json").write_text(
            json.dumps(generated_meta, indent=2) + "\n"
        )

        archive(
            folder,
            ("pw.in", "ph.in", "dynmat.in", "run.sh", "createscf_meta.json"),
        )

        for filename in (
            "pw.in", "ph.in", "dynmat.in",
            "run.sh", "createscf_meta.json",
        ):
            source = staging / filename
            destination = folder / filename

            if source.is_file():
                shutil.copy2(source, destination)
            elif destination.is_file() and filename in {"ph.in", "dynmat.in"}:
                destination.unlink()

    return target


def submit(folder, workflow):
    result = subprocess.run(
        ["sbatch", "--parsable", "run.sh"],
        cwd=str(folder), capture_output=True, text=True, check=True,
    )

    job_id = result.stdout.strip().split(";", 1)[0]
    if not job_id.isdigit():
        raise RuntimeError(f"Unexpected sbatch response: {result.stdout}")

    state = {
        "job_id": job_id,
        "workflow": workflow,
        "submitted_at_utc": now().isoformat(),
        "working_directory": str(folder),
    }
    temporary = folder / "submission_state.json.tmp"
    temporary.write_text(json.dumps(state, indent=2) + "\n")
    temporary.replace(folder / "submission_state.json")
    return job_id


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-d", "--dir", default="QEinputs")
    parser.add_argument("--run", action="store_true")
    parser.add_argument(
        "--workflow",
        choices=["keep", "auto", *WORKFLOWS],
        default="keep",
    )
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--regenerate", action="store_true")
    parser.add_argument("--retry-invalid", action="store_true")
    parser.add_argument("--dev-xy", type=float)
    parser.add_argument("--dev-z", type=float)
    parser.add_argument("--seed", type=int)
    parser.add_argument("--time")
    parser.add_argument("--dielectric-symmetry-tol", type=float, default=1e-3)
    parser.add_argument("--max-dielectric-component", type=float)
    args = parser.parse_args()

    if any(
        value is not None
        for value in (args.dev_xy, args.dev_z, args.seed)
    ) and not args.regenerate:
        parser.error("Deviation/seed overrides require --regenerate")

    for value in (args.dev_xy, args.dev_z):
        if value is not None and not 0 <= value < float("inf"):
            parser.error("Deviations must be finite and nonnegative")

    if args.seed is not None and not 0 <= args.seed < 2**32:
        parser.error("Seed must be between 0 and 2**32-1")

    if not 0 <= args.dielectric_symmetry_tol < float("inf"):
        parser.error("Invalid symmetry tolerance")

    if args.max_dielectric_component is not None:
        if not 0 < args.max_dielectric_component < float("inf"):
            parser.error("Invalid dielectric component limit")

    if args.time and not re.fullmatch(
        r"(?:\d+-)?\d{1,3}:[0-5]\d:[0-5]\d", args.time
    ):
        parser.error("Time must use HH:MM:SS or D-HH:MM:SS")

    root = Path(args.dir).expanduser()
    if not root.is_absolute():
        root = ROOT / root
    root = root.resolve()

    if not root.is_dir():
        parser.error(f"Missing directory: {root}")

    folders = sorted({
        path.parent.resolve()
        for path in root.rglob("pw.in")
        if not any(
            part in {"restart_archive", "workflow_archive"}
            or part.startswith(".qe_regenerate_")
            for part in path.relative_to(root).parts
        )
    })

    try:
        ids, active_folders = active_jobs()
    except Exception as error:
        parser.error(f"Cannot query Slurm; refusing submission: {error}")

    candidates = submitted = failed = 0

    for folder in folders:
        if is_active(folder, ids, active_folders):
            print(f"[ACTIVE: SKIP] {folder}")
            continue

        state, target, errors = classify(folder, args)
        if state == "finished" and not args.force:
            continue

        print(f"[{state.upper()}] {folder}")
        print(f"    Target: {target}")
        for error in errors:
            print(f"    {error}")

        if not eligible(state, args):
            continue

        candidates += 1
        if not args.run:
            continue

        with (folder / ".submission.lock").open("a") as lock:
            try:
                fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError:
                print("    Submission lock held; skipped")
                continue

            try:
                ids, active_folders = active_jobs()
                if is_active(folder, ids, active_folders):
                    continue

                state, _, _ = classify(folder, args)
                if not eligible(state, args):
                    continue

                workflow = (
                    regenerate(folder, args)
                    if args.regenerate
                    else prepare_workflow(folder, args)
                )

                archive(folder, ("submission_state.json",))
                archive(
                    folder,
                    ("pw.out", "ph.out", "dynmat.out"),
                    move=True,
                )

                job_id = submit(folder, workflow)
                submitted += 1
                print(f"    Submitted {job_id}: {workflow}")

            except Exception as error:
                failed += 1
                print(f"    FAILED: {error}")

    print(
        f"\nCandidates={candidates}, submitted={submitted}, errors={failed}"
    )
    if not args.run:
        print("Dry run only. Add --run to apply changes and submit.")

    if failed:
        sys.exit(1)


if __name__ == "__main__":
    main()
