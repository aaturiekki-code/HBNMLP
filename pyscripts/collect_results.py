#!/usr/bin/env python3
import argparse
import csv
import math
import os
import pathlib
import re

def fmt_val(v):
    """Ensure consistent formatting for None and NaN values in the CSV output."""
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return 'NaN'
    return v

def parse_pw_out(filepath):
    energy = None
    total_force = None
    atom_forces = []
    stress_tensor = None
    
    if not os.path.exists(filepath):
        return energy, total_force, atom_forces, stress_tensor
        
    with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
        content = f.read()
        
    # 1. Total Energy (Ry)
    energy_match = re.search(r'!\s+total energy\s+=\s+([\d\.\-]+)\s+Ry', content)
    if energy_match:
        energy = float(energy_match.group(1))
        
    # 2. Total Force (Ry/au)
    force_match = re.search(r'Total force =\s+([\d\.\-]+)', content)
    if force_match:
        total_force = float(force_match.group(1))
        
    # 3. Individual Atom Forces (Ry/au) - Line-by-line parsing
    force_blocks = list(re.finditer(r'Forces acting on atoms \(cartesian axes, Ry/au\):', content))
    if force_blocks:
        last_block_start = force_blocks[-1].end()
        lines_after = content[last_block_start:].split('\n')
        
        for line in lines_after:
            line_strip = line.strip()
            if not line_strip:
                continue
                
            # Stop parsing if we reach any of the sub-contribution headers
            if any(term in line_strip for term in ["contrib.", "contribution", "term", "Dispersion"]):
                break
                
            # Attempt to parse a standard atom force line
            m = re.search(r'atom\s+\d+\s+type\s+\d+\s+force\s+=\s+([\d\.\-]+)\s+([\d\.\-]+)\s+([\d\.\-]+)', line_strip)
            if m:
                atom_forces.append((float(m.group(1)), float(m.group(2)), float(m.group(3))))
            else:
                # If we hit a non-empty line that doesn't match, we have left the force block
                if len(atom_forces) > 0:
                    break
            
    # 4. Stress Tensor (kbar)
    stress_blocks = list(re.finditer(r'total\s+stress\s+\(Ry/bohr\*\*3\)\s+\(kbar\)\s+P=', content))
    if stress_blocks:
        last_stress_start = stress_blocks[-1].start()
        lines_after = content[last_stress_start:].split('\n')
        tensor = []
        for line in lines_after[1:4]:
            parts = line.split()
            if len(parts) >= 6:
                try:
                    # kbar components sit at columns indices 3, 4, 5
                    tensor.append([float(parts[3]), float(parts[4]), float(parts[5])])
                except ValueError:
                    pass
        if len(tensor) == 3:
            stress_tensor = tensor
            
    return energy, total_force, atom_forces, stress_tensor

def parse_dielectric_tensor(content, marker_regex):
    """
    Locates a dielectric tensor header in the file content, then scans succeeding lines,
    skipping any descriptive text or empty lines, to build exactly a 3x3 matrix.
    """
    match = re.search(marker_regex, content, re.IGNORECASE)
    if not match:
        return None
        
    start_idx = match.end()
    lines = content[start_idx:].split('\n')
    tensor = []
    
    # Matches valid floats (with decimals, signs, exponents) OR asterisks sequence from overflows
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
    elec_dielectric = None
    total_dielectric = None
    
    if not os.path.exists(filepath):
        return elec_dielectric, total_dielectric
        
    with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
        content = f.read()
        
    # 1. Parse Electronic Dielectric permittivity Tensor
    elec_dielectric = parse_dielectric_tensor(
        content, 
        r'Electronic dielectric permittivity tensor'
    )
    
    # 2. Parse Total Dielectric Tensor (including polar phonon mode contribution)
    total_dielectric = parse_dielectric_tensor(
        content, 
        r'with zone-center polar mode contributions'
    )
    
    return elec_dielectric, total_dielectric

def main():
    ap = argparse.ArgumentParser(description="Collect QE workflow calculation results from subfolders into separate CSV files.")
    ap.add_argument("-d", "--dir", default=".", type=pathlib.Path,
                    help="Target parent directory containing the sweep subfolders (defaults to current directory).")
    args = ap.parse_args()

    parent_dir = args.dir.resolve()
    print("=========================================")
    print(f"Scanning subfolders inside: {parent_dir}")
    print("=========================================")

    all_data = []
    max_atoms = 0
    
    # Locate directories containing pw.in files to identify run folders
    subdirs = sorted([d for d in parent_dir.iterdir() if d.is_dir() and (d / "pw.in").exists()])
    
    if not subdirs:
        print("No calculation directories (containing 'pw.in') found.")
        return

    for subdir in subdirs:
        pw_file = subdir / "pw.out"
        dyn_file = subdir / "dynmat.out"
        
        energy, total_force, atom_forces, stress = parse_pw_out(pw_file)
        elec_d, total_d = parse_dynmat_out(dyn_file)
        
        max_atoms = max(max_atoms, len(atom_forces))
        
        all_data.append({
            'folder': subdir.name,
            'energy_ry': energy,
            'total_force_ry_au': total_force,
            'atom_forces': atom_forces,
            'stress': stress,
            'elec_dielectric': elec_d,
            'total_dielectric': total_d
        })

    # 1. Write energies.csv
    with open(parent_dir / 'energies.csv', 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['folder', 'energy_ry', 'energy_ev'])
        for d in all_data:
            val_ry = fmt_val(d['energy_ry'])
            val_ev = fmt_val(d['energy_ry'] * 13.6056980659) if d['energy_ry'] is not None else 'NaN'
            writer.writerow([d['folder'], val_ry, val_ev])
            
    # 2. Write forces.csv (Contains ONLY folder and $3 \times N$ component-wise Cartesian forces)
    force_headers = ['folder']
    for i in range(max_atoms):
        force_headers += [f'fx_{i}', f'fy_{i}', f'fz_{i}']
    with open(parent_dir / 'forces.csv', 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(force_headers)
        for d in all_data:
            row = [d['folder']]
            for i in range(max_atoms):
                if i < len(d['atom_forces']):
                    row += [fmt_val(d['atom_forces'][i][0]), fmt_val(d['atom_forces'][i][1]), fmt_val(d['atom_forces'][i][2])]
                else:
                    row += ['NaN', 'NaN', 'NaN']
            writer.writerow(row)

    # 3. Write stresses.csv
    stress_headers = [
        'folder', 
        's_xx_kbar', 's_xy_kbar', 's_xz_kbar', 
        's_yx_kbar', 's_yy_kbar', 's_yz_kbar', 
        's_zx_kbar', 's_zy_kbar', 's_zz_kbar'
    ]
    with open(parent_dir / 'stresses.csv', 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(stress_headers)
        for d in all_data:
            row = [d['folder']]
            if d['stress'] is not None:
                row += [
                    fmt_val(d['stress'][0][0]), fmt_val(d['stress'][0][1]), fmt_val(d['stress'][0][2]),
                    fmt_val(d['stress'][1][0]), fmt_val(d['stress'][1][1]), fmt_val(d['stress'][1][2]),
                    fmt_val(d['stress'][2][0]), fmt_val(d['stress'][2][1]), fmt_val(d['stress'][2][2])
                ]
            else:
                row += ['NaN'] * 9
            writer.writerow(row)

    # 4. Write dielectric_electronic.csv (high-frequency clamped)
    dielec_headers = [
        'folder', 
        'eps_xx', 'eps_xy', 'eps_xz', 
        'eps_yx', 'eps_yy', 'eps_yz', 
        'eps_zx', 'eps_zy', 'eps_zz'
    ]
    with open(parent_dir / 'dielectric_electronic.csv', 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(dielec_headers)
        for d in all_data:
            row = [d['folder']]
            if d['elec_dielectric'] is not None:
                row += [
                    fmt_val(d['elec_dielectric'][0][0]), fmt_val(d['elec_dielectric'][0][1]), fmt_val(d['elec_dielectric'][0][2]),
                    fmt_val(d['elec_dielectric'][1][0]), fmt_val(d['elec_dielectric'][1][1]), fmt_val(d['elec_dielectric'][1][2]),
                    fmt_val(d['elec_dielectric'][2][0]), fmt_val(d['elec_dielectric'][2][1]), fmt_val(d['elec_dielectric'][2][2])
                ]
            else:
                row += ['NaN'] * 9
            writer.writerow(row)

    # 5. Write dielectric_total.csv (static unclamped)
    with open(parent_dir / 'dielectric_total.csv', 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(dielec_headers)
        for d in all_data:
            row = [d['folder']]
            if d['total_dielectric'] is not None:
                row += [
                    fmt_val(d['total_dielectric'][0][0]), fmt_val(d['total_dielectric'][0][1]), fmt_val(d['total_dielectric'][0][2]),
                    fmt_val(d['total_dielectric'][1][0]), fmt_val(d['total_dielectric'][1][1]), fmt_val(d['total_dielectric'][1][2]),
                    fmt_val(d['total_dielectric'][2][0]), fmt_val(d['total_dielectric'][2][1]), fmt_val(d['total_dielectric'][2][2])
                ]
            else:
                row += ['NaN'] * 9
            writer.writerow(row)

    print("Parsed values cleanly written to parent folder:")
    print(f"  -> {parent_dir}/energies.csv")
    print(f"  -> {parent_dir}/forces.csv")
    print(f"  -> {parent_dir}/stresses.csv")
    print(f"  -> {parent_dir}/dielectric_electronic.csv")
    print(f"  -> {parent_dir}/dielectric_total.csv")
    print("\nExtraction complete!")

if __name__ == "__main__":
    main()
