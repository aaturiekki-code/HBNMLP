#!/usr/bin/env python3
import argparse
import os
import pathlib
import re
import numpy as np
from ase import Atoms
from ase.io import read as ase_read, write as ase_write
from ase.data import chemical_symbols

BOHR_TO_ANG = 0.529177210903

def parse_qe_input_fallback(infile: str) -> Atoms:
    """
    Robust manual fallback parser for Quantum ESPRESSO input files.
    Handles crystal, angstrom, and bohr coordinates natively.
    """
    def extract_unit(token_line: str, default: str = "angstrom") -> str:
        tl = token_line.strip().lower()
        m = re.search(r'\b(cell_parameters|atomic_positions)\b\s*\(?\s*([a-z]+)\s*\)?', tl)
        return m.group(2) if m else default

    cell = None
    pos = []
    symbols = []
    
    with open(infile, "r", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()
        
    i = 0
    n = len(lines)
    apos_unit = "angstrom"
    
    while i < n:
        line = lines[i].strip().lower()
        if line.startswith("cell_parameters"):
            unit = extract_unit(lines[i])
            i += 1
            if unit == "alat":
                raise ValueError(f"CELL_PARAMETERS in 'alat' not supported in fallback parser for {infile}")
            vecs = []
            for _ in range(3):
                if i >= n:
                    raise ValueError(f"Unexpected EOF reading CELL_PARAMETERS in {infile}")
                parts = lines[i].split()
                if len(parts) < 3:
                    raise ValueError(f"Malformed CELL_PARAMETERS line in {infile}: {lines[i]}")
                vecs.append([float(parts[0]), float(parts[1]), float(parts[2])])
                i += 1
            cell = np.array(vecs, dtype=float)
            if unit == "bohr":
                cell *= BOHR_TO_ANG
            continue
            
        if line.startswith("atomic_positions"):
            unit = extract_unit(lines[i])
            i += 1
            apos_unit = unit
            if unit == "alat":
                raise ValueError(f"ATOMIC_POSITIONS in 'alat' not supported in fallback parser for {infile}")
            while i < n:
                ln = lines[i].strip()
                if not ln:
                    break
                parts = ln.split()
                if parts[0] not in chemical_symbols:
                    break
                if len(parts) < 4:
                    raise ValueError(f"Malformed ATOMIC_POSITIONS line in {infile}: {lines[i]}")
                symbols.append(parts[0])
                pos.append([float(parts[1]), float(parts[2]), float(parts[3])])
                i += 1
            continue
        i += 1
        
    if cell is None or not pos:
        raise ValueError(f"Failed to parse geometry from {infile} (cell or positions missing).")
        
    positions = np.array(pos, dtype=float)
    if apos_unit == "crystal":
        positions = positions @ cell
    elif apos_unit == "bohr":
        positions *= BOHR_TO_ANG
        
    return Atoms(symbols=symbols, positions=positions, cell=cell, pbc=True)

def main():
    ap = argparse.ArgumentParser(description="Convert a Quantum ESPRESSO input file (e.g., pw.in) to a sorted VASP POSCAR.")
    ap.add_argument("-i", "--input", default="pw.in", type=pathlib.Path,
                    help="Path to the Quantum ESPRESSO input file (default: pw.in)")
    ap.add_argument("-o", "--output", default="POSCAR", type=pathlib.Path,
                    help="Path to the output VASP POSCAR file (default: POSCAR)")
    args = ap.parse_args()

    if not args.input.exists():
        print(f"Error: Input file '{args.input}' does not exist.")
        return

    print(f"Reading Quantum ESPRESSO input: {args.input}")
    
    # Try using ASE's native espresso reader first
    try:
        atoms = ase_read(str(args.input), format="espresso-in")
        print("  -> Parsed successfully via ASE native reader.")
    except Exception as e:
        print(f"  -> ASE native reader failed: {e}. Attempting manual parsing...")
        try:
            atoms = parse_qe_input_fallback(str(args.input))
            print("  -> Parsed successfully via manual fallback parser.")
        except Exception as e_inner:
            print(f"Error: Complete parsing failure: {e_inner}")
            return

    # CRITICAL: Sort atoms by atomic number to group species cleanly for VASP/VESTA standards
    sorted_indices = np.argsort(atoms.numbers)
    sorted_atoms = atoms[sorted_indices]

    # Write POSCAR file
    print(f"Saving sorted structure to: {args.output}")
    ase_write(str(args.output), sorted_atoms, format="vasp", vasp5=True)
    
    print("\n====================================================")
    print("Success! Conversion completed.")
    print("Atomic species are grouped cleanly by element types.")
    print("You can open this POSCAR file directly in VESTA!")
    print("====================================================")

if __name__ == "__main__":
    main()
