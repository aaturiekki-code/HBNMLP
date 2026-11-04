#!/usr/bin/env python3
import numpy as np
import pathlib
from ase import Atoms
from ase.io import write as ase_write

def main():
    data_dir = pathlib.Path(".")
    r = np.load(data_dir / "r.npy")          # Coordinates (N_frames, N_atoms, 3)
    R = np.load(data_dir / "R.npy")          # Lattice vectors (N_frames, 3, 3)
    E = np.load(data_dir / "E.npy")          # Energies (N_frames,)
    F = np.load(data_dir / "F.npy")          # Forces (N_frames, N_atoms, 3)
    stress = np.load(data_dir / "X.npy")      # Stresses (N_frames, 3, 3)
    z = np.load(data_dir / "z.npy")          # True atomic numbers (N_atoms,)

    n_frames, n_atoms, _ = r.shape
    print(f"Loaded dataset: {n_frames} frames with {n_atoms} atoms.")
    print(f"True atomic species sequence: {z}")

    frames = []
    for i in range(n_frames):
        # Create ASE Atoms object
        atoms = Atoms(
            numbers=z,
            positions=r[i],
            cell=R[i],
            pbc=True
        )
        
        # Attach energies, forces, and stresses
        atoms.info["energy"] = E[i]
        atoms.arrays["forces"] = F[i]
        atoms.info["stress"] = stress[i].flatten()
        
        frames.append(atoms)

    # Save to standard Extended XYZ format
    output_path = data_dir / "dataset.xyz"
    ase_write(str(output_path), frames, format="extxyz")
    print(f"Successfully compiled Extended XYZ database saved to: {output_path.resolve()}")

if __name__ == "__main__":
    main()
