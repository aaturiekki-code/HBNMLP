#!/usr/bin/env python3
import argparse
import pathlib
import numpy as np
from ase import Atoms
from ase.io import read, write

def main():
    parser = argparse.ArgumentParser(
        description="Convert a LAMMPS .data file directly to a cleanly formatted, sorted VASP POSCAR."
    )
    parser.add_argument(
        "input_file",
        type=pathlib.Path,
        help="Path to the input LAMMPS data file (e.g., relaxed_structure.data)"
    )
    parser.add_argument(
        "-o", "--output",
        type=pathlib.Path,
        default=None,
        help="Optional: Path to the output POSCAR file (default: POSCAR_<input_filename_base>)"
    )
    
    args = parser.parse_args()

    # 1. Verify input file existence
    data_file = args.input_file
    if not data_file.exists():
        print(f"Error: Could not find the file '{data_file}'")
        return

    # 2. Determine output filename
    if args.output is None:
        # If input is 'relaxed_structure.data', base is 'relaxed_structure'
        # Output becomes 'POSCAR_relaxed_structure'
        output_file = pathlib.Path(f"POSCAR_{data_file.stem}")
    else:
        output_file = args.output

    print(f"Reading LAMMPS data file: {data_file}")

    # 3. Read the LAMMPS data file
    raw_atoms = read(str(data_file), format="lammps-data", style="atomic")

    # 4. Map integer types to their chemical species:
    # Type 1 -> Boron (B, atomic number 5), Type 2 -> Nitrogen (N, atomic number 7)
    types = raw_atoms.arrays["type"]
    numbers = [5 if t == 1 else 7 for t in types]

    # 5. Create a clean Atoms object (stripping extra LAMMPS dynamics columns)
    clean_atoms = Atoms(
        numbers=numbers,
        positions=raw_atoms.get_positions(),
        cell=raw_atoms.get_cell(),
        pbc=True
    )

    # 6. Sort atoms by atomic number to group B and N cleanly for VASP/VESTA standard
    sorted_indices = np.argsort(clean_atoms.numbers)
    sorted_atoms = clean_atoms[sorted_indices]

    # 7. Export directly to VASP POSCAR format
    write(str(output_file), sorted_atoms, format="vasp", vasp5=True)
    
    print("====================================================")
    print("Success! Created POSCAR File:")
    print(f"  -> {output_file}")
    print("\nYou can open this file directly in VESTA to visualize the box!")
    print("====================================================")

if __name__ == "__main__":
    main()
