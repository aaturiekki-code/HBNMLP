#!/usr/bin/env python3
import argparse
import glob
import os
import pathlib
import numpy as np
from ase.io import read as ase_read, write as ase_write

def main():
    parser = argparse.ArgumentParser(
        description="Compile Quantum ESPRESSO output files into a clean dataset for training."
    )
    parser.add_argument(
        "--dirs",
        nargs="+",
        required=True,
        help="List of directories containing the QE sweeps (supports wildcards, e.g., QEinputs/*_sweep*)"
    )
    parser.add_argument(
        "-o", "--output-dir",
        default="compiled_dataset",
        type=pathlib.Path,
        help="Directory where compiled database files will be stored"
    )
    # Parse but keep these arguments as dummy placeholders to keep workflow commands backward-compatible
    parser.add_argument("--l-max", type=int, default=2, help="Placeholder for descriptors")
    parser.add_argument("--n-max", type=int, default=8, help="Placeholder for descriptors")
    parser.add_argument("--r-cut", type=float, default=6.0, help="Placeholder for descriptors")

    args = parser.parse_args()

    print("=========================================================================")
    print(" Starting Quantum ESPRESSO Dataset Compiler")
    print("=========================================================================")

    # 1. Resolve and expand wildcard directory paths
    target_dirs = []
    for d in args.dirs:
        expanded = glob.glob(d)
        if expanded:
            target_dirs.extend(expanded)
        else:
            target_dirs.append(d)

    print(f"Scanning {len(target_dirs)} parent directories for completed 'pw.out' files...")

    # 2. Find all pw.out files recursively
    pw_out_files = []
    for parent in target_dirs:
        p_path = pathlib.Path(parent)
        if not p_path.exists():
            print(f"  -> Warning: Path does not exist, skipping: {parent}")
            continue
        # Search recursively for pw.out
        found = list(p_path.rglob("pw.out"))
        pw_out_files.extend(found)

    print(f"Found {len(pw_out_files)} total 'pw.out' files to process.")
    if not pw_out_files:
        print("Error: No 'pw.out' files found. Are your calculations still running?")
        return

    # Create the destination directory
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # 3. Parse files using ASE
    valid_frames = []
    skipped_count = 0

    print("\nParsing calculations...")
    for idx, f_path in enumerate(pw_out_files):
        try:
            # ASE reads the last image (fully completed step) from the QE output
            atoms = ase_read(str(f_path), format="espresso-out")
            
            # Verify required calculations exist (Energy and Forces are mandatory)
            _ = atoms.get_potential_energy()
            _ = atoms.get_forces()
            
            # Store structural metadata to prevent future parsing mismatches
            atoms.info["source_file"] = str(f_path.relative_to(f_path.parents[2]))
            valid_frames.append(atoms)
            
            if (len(valid_frames)) % 50 == 0:
                print(f"  -> Successfully parsed {len(valid_frames)} files...")
                
        except Exception as e:
            # Silently skip unconverged, crashed, or currently running calculations
            skipped_count += 1
            continue

    print(f"\nParsing Summary:")
    print(f"  -> Successfully parsed: {len(valid_frames)} frames")
    print(f"  -> Skipped/Unfinished:  {skipped_count} frames")

    if not valid_frames:
        print("Error: No valid finished calculations were found.")
        return

    # 4. Save to Extended XYZ Format (Standard for Allegro/NequIP)
    output_xyz = args.output_dir / "dataset.xyz"
    print(f"\nWriting clean consolidated dataset to: {output_xyz}")
    
    # Write as standard ExtXYZ with energy, forces, and stress automatically formatted
    ase_write(str(output_xyz), valid_frames, format="extxyz")

    # 5. Backward Compatibility Export (Numpy Arrays)
    # Note: Since atom counts vary, these are saved as object arrays if lengths differ
    print(f"Saving backup numpy files to: {args.output_dir}/")
    try:
        energies = np.array([atoms.get_potential_energy() for atoms in valid_frames])
        np.save(args.output_dir / "E.npy", energies)

        # Check if all systems have the same number of atoms
        num_atoms_list = [len(atoms) for atoms in valid_frames]
        all_same_size = len(set(num_atoms_list)) == 1

        if all_same_size:
            # Uniform dataset sizes: Save standard rectangular numpy arrays
            positions = np.array([atoms.get_positions() for atoms in valid_frames])
            forces = np.array([atoms.get_forces() for atoms in valid_frames])
            boxes = np.array([np.array(atoms.get_cell()) for atoms in valid_frames])
            atomic_numbers = valid_frames[0].get_atomic_numbers()

            np.save(args.output_dir / "R.npy", positions)
            np.save(args.output_dir / "F.npy", forces)
            np.save(args.output_dir / "box.npy", boxes)
            np.save(args.output_dir / "z.npy", atomic_numbers)
        else:
            # Variable dataset sizes: Save as object arrays to prevent shape crashes
            positions = np.array([atoms.get_positions() for atoms in valid_frames], dtype=object)
            forces = np.array([atoms.get_forces() for atoms in valid_frames], dtype=object)
            boxes = np.array([np.array(atoms.get_cell()) for atoms in valid_frames], dtype=object)
            atomic_numbers = np.array([atoms.get_atomic_numbers() for atoms in valid_frames], dtype=object)

            np.save(args.output_dir / "R.npy", positions)
            np.save(args.output_dir / "F.npy", forces)
            np.save(args.output_dir / "box.npy", boxes)
            np.save(args.output_dir / "z.npy", atomic_numbers)
            
        print("  -> Backup Numpy arrays successfully saved!")
    except Exception as e:
        print(f"  -> Warning during numpy backup export: {e}")

    print("\n=========================================================================")
    print(" Process Completed Successfully!")
    print(f" Ready to Train! Target File: {output_xyz}")
    print("=========================================================================")

if __name__ == "__main__":
    main()
