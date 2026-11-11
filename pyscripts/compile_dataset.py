#!/usr/bin/env python3
"""Compile independent Allegro and SchNetPack datasets.

Script location:
    project/pyscripts/compile_dataset.py

Each calculation directory:
    pw.in
    createscf_meta.json
    pw.out
    ph.out       optional for SCF-only calculations
    dynmat.out   optional legacy dielectric source

Relative paths are resolved from project/.

Allegro:
    Untwisted, m0=1, m0=2.
    Parent-grouped train/validation/test within each angle family.

Dielectric:
    Untwisted train/validation.
    m0=1 test.
    m0=2 excluded.

Requires qe_response.py in the same directory.
"""

import argparse
import csv
import glob
import hashlib
import json
from pathlib import Path
import re

import numpy as np
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import read as ase_read
from ase.io import write as ase_write
from ase.units import Bohr, Ry

from qe_response import read_folder_tensor


PROJECT_ROOT = Path(__file__).resolve().parent.parent
SUBSETS = ("train", "val", "test")

NUMBER = r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?"
NUMBER_RE = re.compile(NUMBER)

ENERGY_RE = re.compile(
    rf"!\s+total\s+energy\s*=\s*({NUMBER})\s+Ry",
    re.IGNORECASE,
)
FORCE_RE = re.compile(
    rf"atom\s+(\d+)\s+type\s+\d+\s+force\s*=\s*"
    rf"({NUMBER})\s+({NUMBER})\s+({NUMBER})",
    re.IGNORECASE,
)
STRESS_RE = re.compile(
    r"total\s+stress\s*\(\s*Ry/bohr",
    re.IGNORECASE,
)


def qe_float(value):
    return float(value.replace("D", "E").replace("d", "e"))


def project_path(value):
    path = Path(value).expanduser()
    return path if path.is_absolute() else PROJECT_ROOT / path


def file_hash(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def object_hash(value):
    text = json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False
    )
    return hashlib.sha256(text.encode()).hexdigest()


def collect_inputs(patterns):
    files = set()

    for pattern in patterns:
        expanded = str(project_path(pattern))
        for match in glob.glob(expanded) or [expanded]:
            directory = Path(match)
            if not directory.is_dir():
                print(f"[WARNING] Not a directory: {directory}")
                continue

            for path in directory.rglob("pw.in"):
                parts = path.relative_to(directory).parts
                if any(
                    part in {"restart_archive", "workflow_archive"}
                    or part.startswith(".qe_regenerate_")
                    for part in parts
                ):
                    continue
                files.add(path.resolve())

    return sorted(files)


def read_structure(path):
    """Read fixed input geometry; exported for the restart script."""
    path = Path(path)
    text = path.read_text(errors="replace")
    uncommented = "\n".join(
        line.split("!", 1)[0] for line in text.splitlines()
    )

    match = re.search(
        r"\bcalculation\s*=\s*['\"]([^'\"]+)['\"]",
        uncommented,
        re.IGNORECASE,
    )
    calculation = match.group(1).strip().lower() if match else "scf"

    if calculation != "scf":
        raise ValueError(
            f"Expected calculation='scf', found '{calculation}'"
        )

    atoms = ase_read(str(path), format="espresso-in")

    if not len(atoms):
        raise ValueError("Empty structure")
    if not np.isfinite(atoms.positions).all():
        raise ValueError("Nonfinite positions")
    if not np.isfinite(atoms.cell.array).all():
        raise ValueError("Nonfinite cell")
    if abs(np.linalg.det(atoms.cell.array)) < 1e-12:
        raise ValueError("Singular cell")

    return atoms


def read_forcefield_results(path, natoms):
    """Parse QE results directly; exported for the restart script."""
    lines = Path(path).read_text(errors="replace").splitlines()
    text = "\n".join(lines)

    if "JOB DONE." not in text:
        raise ValueError("QE completion marker missing")
    if re.search(r"convergence\s+NOT\s+achieved", text, re.I):
        raise ValueError("SCF convergence not achieved")
    if not re.search(r"convergence\s+has\s+been\s+achieved", text, re.I):
        raise ValueError("SCF convergence confirmation missing")

    energies = [
        qe_float(match.group(1))
        for line in lines
        if (match := ENERGY_RE.search(line))
    ]
    if not energies:
        raise ValueError("Total energy missing")

    force_headers = [
        index for index, line in enumerate(lines)
        if "forces acting on atoms" in line.lower()
    ]
    if not force_headers:
        raise ValueError("Force block missing")

    force_rows = {}
    for line in lines[force_headers[-1] + 1:]:
        match = FORCE_RE.search(line)
        if match:
            index = int(match.group(1)) - 1
            if index in force_rows:
                raise ValueError("Duplicate atom in force block")
            force_rows[index] = [
                qe_float(match.group(i)) for i in (2, 3, 4)
            ]
        elif force_rows and line.strip():
            break

    if set(force_rows) != set(range(natoms)):
        raise ValueError("Force block atom count/index mismatch")

    forces = np.asarray(
        [force_rows[index] for index in range(natoms)],
        dtype=np.float64,
    ) * Ry / Bohr

    stress_headers = [
        index for index, line in enumerate(lines)
        if STRESS_RE.search(line)
    ]
    if not stress_headers:
        raise ValueError("Stress block missing")

    rows = []
    for line in lines[stress_headers[-1] + 1:]:
        if not line.strip():
            continue
        tokens = line.split()
        if (
            "*" in line
            or len(tokens) not in (3, 6)
            or not all(NUMBER_RE.fullmatch(token) for token in tokens)
        ):
            raise ValueError("Malformed or overflowed stress block")

        rows.append([qe_float(token) for token in tokens[:3]])
        if len(rows) == 3:
            break

    if len(rows) != 3:
        raise ValueError("Incomplete stress block")

    stress = np.asarray(rows, dtype=np.float64)
    stress = -0.5 * (stress + stress.T) * Ry / Bohr**3

    # ASE Voigt order: xx, yy, zz, yz, xz, xy.
    stress = stress.reshape(-1)[[0, 4, 8, 5, 2, 1]]

    results = {
        "energy": energies[-1] * Ry,
        "forces": forces,
        "stress": stress,
    }
    if not all(np.isfinite(value).all() for value in results.values()):
        raise ValueError("Nonfinite force-field labels")

    return results


def read_metadata(folder):
    raw = json.loads(
        (folder / "createscf_meta.json").read_text()
    )
    mode = raw["mode"]

    if mode == "bilayer":
        m0 = None
        angle = 0.0
        angle_family = "untwisted"
        stacking = raw.get("bilayer_stack")
        if stacking is None:
            stacking = "AB" if raw.get("ab", 0) == 1 else "AA"
        ab = None

    elif mode == "twist":
        m0 = int(raw["m0"])
        if m0 not in (1, 2):
            raise ValueError("This dataset policy supports m0=1 and m0=2")

        # r=1 is explicitly fixed in your supplied generator.
        denominator = 3 * m0**2 + 3 * m0 + 1
        numerator = 3 * m0**2 + 3 * m0 + 0.5
        angle = float(np.degrees(np.arccos(numerator / denominator)))
        angle_family = f"m0_{m0}"
        stacking = None
        ab = int(raw.get("ab", 0))

    else:
        raise ValueError(f"Unknown mode: {mode}")

    slide = raw.get("slide")
    slide = [0.0, 0.0] if slide is None else list(map(float, slide))

    parent = {
        "mode": mode,
        "m0": m0,
        "r": 1 if mode == "twist" else None,
        "stacking": stacking,
        "ab": ab,
        "mirror": int(raw.get("mirror", 0)),
        "distance": raw.get("distance"),
        "scale_xy": raw.get("scale_xy"),
        "slide": slide,
        "vacuum": raw.get("vacuum"),
    }

    # Seed, deviations, and supercell expansion are deliberately excluded.
    return {
        "raw": raw,
        "parent_family_id": object_hash(parent),
        "angle_family": angle_family,
        "twist_angle_deg": angle,
        "m0": m0,
    }


def assign_groups(groups, val_fraction, test_fraction, rng):
    groups = sorted(set(groups))
    minimum = 3 if test_fraction > 0 else 2

    if len(groups) < minimum:
        raise ValueError(
            f"Need at least {minimum} parent groups; found {len(groups)}"
        )

    groups = rng.permutation(groups).tolist()
    n_val = max(1, round(len(groups) * val_fraction))
    n_test = (
        max(1, round(len(groups) * test_fraction))
        if test_fraction > 0 else 0
    )

    while n_val + n_test >= len(groups):
        if n_val > 1:
            n_val -= 1
        elif n_test > 1:
            n_test -= 1
        else:
            raise ValueError("Cannot make nonempty grouped subsets")

    assignments = {}
    for index, group in enumerate(groups):
        if index < n_test:
            subset = "test"
        elif index < n_test + n_val:
            subset = "val"
        else:
            subset = "train"
        assignments[group] = subset

    return assignments


def write_csv(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dirs", nargs="+", required=True)
    parser.add_argument("-o", "--output-dir", default="compiled_dataset")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--allegro-val-fraction", type=float, default=0.1)
    parser.add_argument("--allegro-test-fraction", type=float, default=0.1)
    parser.add_argument("--dielectric-val-fraction", type=float, default=0.2)
    parser.add_argument("--clean-dielectric-test-only", action="store_true")
    parser.add_argument("--dielectric-symmetry-tol", type=float, default=1e-3)
    parser.add_argument("--max-dielectric-component", type=float)
    args = parser.parse_args()

    fractions = (
        args.allegro_val_fraction,
        args.allegro_test_fraction,
        args.dielectric_val_fraction,
    )
    if not all(0 < value < 1 for value in fractions):
        parser.error("Split fractions must lie between 0 and 1")
    if args.allegro_val_fraction + args.allegro_test_fraction >= 1:
        parser.error("Allegro validation + test fractions must be below 1")
    if not 0 <= args.dielectric_symmetry_tol < float("inf"):
        parser.error("Invalid symmetry tolerance")
    if args.max_dielectric_component is not None:
        if not 0 < args.max_dielectric_component < float("inf"):
            parser.error("Invalid dielectric magnitude limit")

    output = project_path(args.output_dir).resolve()
    names = [
        "dataset.xyz", "train.xyz", "val.xyz", "test.xyz",
        "dielectric.db",
        "allegro_split_indices.npz", "dielectric_split.npz",
        "sample_manifest.csv", "compilation_report.csv",
        "dataset_metadata.json",
    ]
    existing = [
        output / name for name in names if (output / name).exists()
    ]
    if existing and not args.overwrite:
        parser.error("Outputs exist; use --overwrite or a new directory")

    paths = collect_inputs(args.dirs)
    if not paths:
        parser.error("No pw.in files found")

    samples = []
    for path in paths:
        try:
            metadata = read_metadata(path.parent)
        except Exception as error:
            parser.error(f"Metadata error in {path.parent}: {error}")
        samples.append({"path": path, "meta": metadata})

    angle_families = {s["meta"]["angle_family"] for s in samples}
    required_families = {"untwisted", "m0_1", "m0_2"}
    if not required_families.issubset(angle_families):
        parser.error(
            "Missing required angle families: "
            + ", ".join(sorted(required_families - angle_families))
        )

    allegro_rng = np.random.default_rng(args.seed)
    dielectric_rng = np.random.default_rng(args.seed + 1)
    allegro_assignments = {}

    try:
        for angle_family in sorted(angle_families):
            groups = [
                s["meta"]["parent_family_id"]
                for s in samples
                if s["meta"]["angle_family"] == angle_family
            ]
            allegro_assignments.update(assign_groups(
                groups,
                args.allegro_val_fraction,
                args.allegro_test_fraction,
                allegro_rng,
            ))

        dielectric_assignments = assign_groups(
            [
                s["meta"]["parent_family_id"]
                for s in samples
                if s["meta"]["angle_family"] == "untwisted"
            ],
            args.dielectric_val_fraction,
            0.0,
            dielectric_rng,
        )
    except ValueError as error:
        parser.error(str(error))

    frames, db_atoms, db_properties = [], [], []
    allegro_indices = {subset: [] for subset in SUBSETS}
    dielectric_indices = {subset: [] for subset in SUBSETS}
    manifest, report = [], []

    for sample in samples:
        path = sample["path"]
        folder = path.parent
        meta = sample["meta"]
        raw = meta["raw"]
        group = meta["parent_family_id"]
        allegro_split = allegro_assignments[group]

        if meta["angle_family"] == "untwisted":
            dielectric_split = dielectric_assignments[group]
        elif meta["m0"] == 1:
            perturbed = (
                float(raw.get("dev_xy", 0)) != 0
                or float(raw.get("dev_z", 0)) != 0
            )
            dielectric_split = (
                "excluded"
                if args.clean_dielectric_test_only and perturbed
                else "test"
            )
        else:
            dielectric_split = "excluded"

        try:
            folder_identifier = str(folder.relative_to(PROJECT_ROOT))
        except ValueError:
            folder_identifier = str(folder)

        sample_id = object_hash(folder_identifier)

        row = {
            "sample_id": sample_id,
            "calculation_dir": folder_identifier,
            "parent_family_id": group,
            "angle_family": meta["angle_family"],
            "twist_angle_deg": meta["twist_angle_deg"],
            "allegro_split": allegro_split,
            "dielectric_split": dielectric_split,
            "xyz_index": "",
            "db_index": "",
            "tensor_source": "",
            "geometry_sha256": file_hash(path),
            "metadata_json": json.dumps(raw, sort_keys=True),
        }
        status = {
            "sample_id": sample_id,
            "calculation_dir": folder_identifier,
            "structure_status": "",
            "allegro_status": "not processed",
            "dielectric_status": (
                "excluded" if dielectric_split == "excluded"
                else "not processed"
            ),
        }

        try:
            atoms = read_structure(path)
            atoms.info.update({
                "sample_id": sample_id,
                "parent_family_id": group,
                "twist_angle_deg": meta["twist_angle_deg"],
                "source_file": str(path),
            })
            status["structure_status"] = "accepted"
        except Exception as error:
            status["structure_status"] = f"skipped: {error}"
            manifest.append(row)
            report.append(status)
            continue

        try:
            results = read_forcefield_results(folder / "pw.out", len(atoms))
            frame = atoms.copy()
            frame.info["split"] = allegro_split
            frame.calc = SinglePointCalculator(frame, **results)

            index = len(frames)
            frames.append(frame)
            allegro_indices[allegro_split].append(index)
            row["xyz_index"] = index
            status["allegro_status"] = "accepted"
        except Exception as error:
            status["allegro_status"] = f"skipped: {error}"

        if dielectric_split != "excluded":
            try:
                tensor, source = read_folder_tensor(
                    folder,
                    args.dielectric_symmetry_tol,
                    args.max_dielectric_component,
                )
                frame = atoms.copy()
                frame.calc = None
                frame.info["split"] = dielectric_split
                frame.info["tensor_source"] = str(source)

                index = len(db_atoms)
                db_atoms.append(frame)
                db_properties.append({
                    "dielectric_tensor": tensor[None, :, :]
                })
                dielectric_indices[dielectric_split].append(index)

                row["db_index"] = index
                row["tensor_source"] = str(source)
                status["dielectric_status"] = "accepted"
            except Exception as error:
                status["dielectric_status"] = f"skipped: {error}"

        manifest.append(row)
        report.append(status)

    output.mkdir(parents=True, exist_ok=True)

    empty = [
        f"{model}/{subset}"
        for model, indices in (
            ("Allegro", allegro_indices),
            ("dielectric", dielectric_indices),
        )
        for subset in SUBSETS
        if not indices[subset]
    ]
    if empty:
        write_csv(output / "compilation_failed_report.csv", report)
        parser.error(
            "Empty subsets after filtering: " + ", ".join(empty)
            + ". Existing datasets were not replaced."
        )

    for indices, count in (
        (allegro_indices, len(frames)),
        (dielectric_indices, len(db_atoms)),
    ):
        flat = [i for subset in SUBSETS for i in indices[subset]]
        if sorted(flat) != list(range(count)):
            raise RuntimeError("Split indices overlap or omit samples")

    # Import before replacing outputs, so missing dependencies fail safely.
    from schnetpack.data import ASEAtomsData

    if args.overwrite:
        for path in existing:
            path.unlink()

    ase_write(
        str(output / "dataset.xyz"), frames,
        format="extxyz", write_results=True,
    )
    for subset in SUBSETS:
        ase_write(
            str(output / f"{subset}.xyz"),
            [frames[i] for i in allegro_indices[subset]],
            format="extxyz", write_results=True,
        )

    database = ASEAtomsData.create(
        str(output / "dielectric.db"),
        distance_unit="Ang",
        property_unit_dict={"dielectric_tensor": ""},
    )
    database.add_systems(
        property_list=db_properties,
        atoms_list=db_atoms,
    )

    for filename, indices in (
        ("allegro_split_indices.npz", allegro_indices),
        ("dielectric_split.npz", dielectric_indices),
    ):
        np.savez(
            output / filename,
            **{
                f"{subset}_idx": np.asarray(indices[subset], dtype=np.int64)
                for subset in SUBSETS
            },
        )

    write_csv(output / "sample_manifest.csv", manifest)
    write_csv(output / "compilation_report.csv", report)

    summary = {
        "seed": args.seed,
        "split_policy": {
            "allegro": (
                "Parent-grouped within untwisted, m0=1, and m0=2"
            ),
            "dielectric": (
                "Untwisted train/validation; m0=1 test; m0=2 excluded"
            ),
        },
        "split_fractions": {
            "allegro_validation": args.allegro_val_fraction,
            "allegro_test": args.allegro_test_fraction,
            "dielectric_validation": args.dielectric_val_fraction,
        },
        "clean_dielectric_test_only": args.clean_dielectric_test_only,
        "filters": {
            "symmetry_tolerance": args.dielectric_symmetry_tol,
            "max_dielectric_component": args.max_dielectric_component,
        },
        "units": {
            "energy": "eV",
            "forces": "eV/Angstrom",
            "stress": "eV/Angstrom^3, ASE convention",
            "dielectric_tensor": "dimensionless relative permittivity",
        },
        "dielectric_sample_shape": [1, 3, 3],
        "counts": {
            "allegro": {
                subset: len(allegro_indices[subset]) for subset in SUBSETS
            },
            "dielectric": {
                subset: len(dielectric_indices[subset]) for subset in SUBSETS
            },
        },
        "output_sha256": {
            name: file_hash(output / name)
            for name in names
            if name != "dataset_metadata.json"
        },
    }
    (output / "dataset_metadata.json").write_text(
        json.dumps(summary, indent=2, allow_nan=False) + "\n"
    )

    print(json.dumps(summary["counts"], indent=2))
    print(f"Outputs: {output}")


if __name__ == "__main__":
    main()
