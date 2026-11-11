#!/usr/bin/env python3
"""Shared QE electronic-response parsing and workflow utilities."""

import json
from pathlib import Path
import re

import numpy as np


NUMBER = r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?"
NUMBER_RE = re.compile(NUMBER)

PH_HEADER = re.compile(
    r"Dielectric\s+constant\s+in\s+cartesian\s+ax(?:is|es)",
    re.IGNORECASE,
)

DYNMAT_HEADER = re.compile(
    r"Electronic\s+dielectric\s+permittivity\s+tensor",
    re.IGNORECASE,
)


def qe_float(value):
    """Convert a QE number, including Fortran D exponents."""
    return float(value.replace("D", "E").replace("d", "e"))


def read_meta(folder):
    """Read existing structure-generation metadata."""
    path = Path(folder) / "createscf_meta.json"
    if not path.is_file():
        raise FileNotFoundError(f"Missing metadata: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def actual_pw_identifier(folder, name):
    """Read prefix or outdir from pw.in."""
    path = Path(folder) / "pw.in"
    text = path.read_text(encoding="utf-8", errors="replace")

    text = "\n".join(
        line.split("!", 1)[0] for line in text.splitlines()
    )

    match = re.search(
        rf"\b{re.escape(name)}\s*=\s*['\"]([^'\"]+)['\"]",
        text,
        flags=re.IGNORECASE,
    )

    if not match:
        raise ValueError(f"Cannot read {name} from {path}")

    return match.group(1)


def read_tensor_rows(lines, header_index, parenthesized):
    """Read exactly three numeric rows after a tensor header."""
    rows = []

    for line in lines[header_index + 1:]:
        if not line.strip():
            continue

        if "*" in line:
            raise ValueError(
                "Electronic dielectric tensor contains overflow markers (*)"
            )

        if parenthesized:
            match = re.fullmatch(r"\s*\(\s*(.*?)\s*\)\s*", line)
            if match is None:
                raise ValueError(
                    f"Malformed dielectric tensor row: {line.strip()}"
                )
            line = match.group(1)

        tokens = line.split()

        if (
            len(tokens) != 3
            or not all(NUMBER_RE.fullmatch(token) for token in tokens)
        ):
            raise ValueError(
                f"Invalid dielectric tensor row: {line.strip()}"
            )

        rows.append([qe_float(token) for token in tokens])

        if len(rows) == 3:
            break

    if len(rows) != 3:
        raise ValueError("Incomplete electronic dielectric tensor")

    return np.asarray(rows, dtype=np.float64)


def validate_tensor(tensor, symmetry_tol=1e-3, max_component=None):
    if not np.isfinite(tensor).all():
        raise ValueError("Nonfinite electronic dielectric tensor")

    relative_asymmetry = (
        np.linalg.norm(tensor - tensor.T)
        / max(np.linalg.norm(tensor), 1.0)
    )

    if relative_asymmetry > symmetry_tol:
        raise ValueError(
            f"Relative tensor asymmetry {relative_asymmetry:.3e} "
            f"exceeds tolerance {symmetry_tol:.3e}"
        )

    if (
        max_component is not None
        and np.max(np.abs(tensor)) > max_component
    ):
        raise ValueError(
            f"Dielectric component exceeds configured limit {max_component}"
        )

    # Preserve printed values; no silent symmetrization.
    return tensor


def read_ph_tensor(
    path,
    symmetry_tol=1e-3,
    max_component=None,
    require_job_done=False,
):
    """Read electronic permittivity from ph.out.

    With require_job_done=False, a completed electronic-response block
    can be retained even if later legacy phonon calculations stopped.
    """
    path = Path(path)
    text = path.read_text(encoding="utf-8", errors="replace")
    lines = text.splitlines()

    headers = [
        index
        for index, line in enumerate(lines)
        if PH_HEADER.search(line)
    ]

    if not headers:
        raise ValueError("Electronic dielectric tensor missing from ph.out")

    header_index = headers[-1]

    field_starts = [
        index
        for index, line in enumerate(lines[:header_index])
        if re.fullmatch(
            r"\s*Electric\s+Fields\s+Calculation\s*",
            line,
            flags=re.IGNORECASE,
        )
    ]

    if not field_starts:
        raise ValueError("Electric-field calculation start marker missing")

    section = "\n".join(
        lines[field_starts[-1]:header_index]
    )

    if not re.search(
        r"End\s+of\s+electric\s+fields\s+calculation",
        section,
        flags=re.IGNORECASE,
    ):
        raise ValueError("Electric-field calculation end marker missing")

    if re.search(
        r"convergence\s+NOT\s+achieved",
        section,
        flags=re.IGNORECASE,
    ):
        raise ValueError("Electric-field response did not converge")

    if require_job_done and "JOB DONE." not in text:
        raise ValueError("ph.x normal-completion marker missing")

    tensor = read_tensor_rows(
        lines, header_index, parenthesized=True
    )

    return validate_tensor(tensor, symmetry_tol, max_component)


def read_dynmat_tensor(path, symmetry_tol=1e-3, max_component=None):
    """Read the electronic-only tensor from legacy dynmat.out."""
    path = Path(path)
    lines = path.read_text(
        encoding="utf-8", errors="replace"
    ).splitlines()

    headers = [
        index
        for index, line in enumerate(lines)
        if DYNMAT_HEADER.search(line)
    ]

    if not headers:
        raise ValueError(
            "Electronic dielectric tensor missing from dynmat.out"
        )

    tensor = read_tensor_rows(
        lines, headers[-1], parenthesized=False
    )

    return validate_tensor(tensor, symmetry_tol, max_component)


def read_folder_tensor(folder, symmetry_tol=1e-3, max_component=None):
    """Return (electronic tensor, source path).

    Prefer ph.out. Fall back to dynmat.out only if ph.out has no tensor
    header. A present but invalid ph.out tensor is not hidden by fallback.
    """
    folder = Path(folder)
    ph_path = folder / "ph.out"

    if ph_path.is_file():
        text = ph_path.read_text(
            encoding="utf-8", errors="replace"
        )
        if PH_HEADER.search(text):
            tensor = read_ph_tensor(
                ph_path,
                symmetry_tol,
                max_component,
            )
            return tensor, ph_path

    dynmat_path = folder / "dynmat.out"

    if dynmat_path.is_file():
        tensor = read_dynmat_tensor(
            dynmat_path,
            symmetry_tol,
            max_component,
        )
        return tensor, dynmat_path

    raise ValueError(
        "No electronic tensor found in ph.out or dynmat.out"
    )


def write_ph(prefix, outdir, filename, electronic_only=False):
    """Write electronic-only or full phonon-response input."""
    lines = [
        "&INPUTPH",
        f"  prefix = '{prefix}',",
        f"  outdir = '{outdir}',",
        "  tr2_ph = 1.0d-12,",
        "  epsil = .true.,",
        f"  trans = {'.false.' if electronic_only else '.true.'},",
    ]

    if not electronic_only:
        lines.append(f"  fildyn = '{prefix}.dyn',")

    lines.extend([
        "/",
        "0.0 0.0 0.0",
        "",
    ])

    Path(filename).write_text(
        "\n".join(lines),
        encoding="utf-8",
    )


def write_runner(
    job_name,
    filename,
    cores,
    time_limit,
    mem_limit,
    electronic_only=False,
):
    """Compatibility writer for earlier generator versions."""
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

echo "Starting SCF"
srun pw.x -np {cores} -inp pw.in > pw.out

echo "Starting response calculation"
srun ph.x -np {cores} -inp ph.in > ph.out
"""

    if not electronic_only:
        text += """
echo "Starting dynmat"
srun --ntasks=1 dynmat.x -in dynmat.in > dynmat.out
"""

    text += '\necho "Workflow completed."\n'

    path = Path(filename)
    path.write_text(text, encoding="utf-8")
    path.chmod(path.stat().st_mode | 0o111)

def write_workflow_runner(
    job_name,
    filename,
    cores,
    time_limit,
    mem_limit,
    workflow,
):
    """Write a runner for SCF-only, electronic-only, or full response."""
    if workflow not in {
        "scf-only",
        "electronic-only",
        "full-response",
    }:
        raise ValueError(f"Unknown workflow: {workflow}")

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

echo "Starting SCF"
srun pw.x -np {cores} -inp pw.in > pw.out
"""

    if workflow != "scf-only":
        text += f"""
echo "Starting response calculation"
srun ph.x -np {cores} -inp ph.in > ph.out
"""

    if workflow == "full-response":
        text += """
echo "Starting dynmat"
srun --ntasks=1 dynmat.x -in dynmat.in > dynmat.out
"""

    text += '\necho "Workflow completed."\n'

    path = Path(filename)
    path.write_text(text, encoding="utf-8")
    path.chmod(path.stat().st_mode | 0o111)
