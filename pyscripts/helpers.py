import re
import argparse
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import pearsonr

def extract_cell_from_input(qe_input_file):
    """
    Extract lattice vectors from Quantum ESPRESSO input file.

    Parameters:
        qe_input_file (str): Path to the QE input file.

    Returns:
        list: Lattice vectors as a list of lists, or None if not found.
    """
    lattice_vectors = []
    with open(qe_input_file, 'r') as f:
        lines = f.readlines()

    cell_parameters_start = None
    for i, line in enumerate(lines):
        if line.strip().startswith("CELL_PARAMETERS"):
            cell_parameters_start = i + 1
            break

    if cell_parameters_start is not None:
        for i in range(cell_parameters_start, cell_parameters_start + 3):
            parts = lines[i].split()
            lattice_vectors.append([float(parts[0]), float(parts[1]), float(parts[2])])
        return lattice_vectors

    return None

def interlayer_distance(qe_input_file, nlayer=2, c=30):
    with open(qe_input_file, 'r') as f:
        lines = f.readlines()
    startpos = 0
    for i, line in enumerate(lines):
        if "ATOMIC_POSITIONS (crystal)" in line:
            startpos = i+1
        elif startpos and (line.strip() == "" or line.startswith("!END")):
            endpos = i
            break
    if nlayer == 2: z = [[], []] 
    elif nlayer ==3: z = [[], [], []]
    for line in lines[startpos:endpos]:
        zi = float(line.split()[3])
        if nlayer == 2:
            if zi < 0.5: z[0].append(zi)
            else: z[1].append(zi)
        elif nlayer == 3:
            if zi < 0.45: z[0].append(zi)
            elif zi > 0.55: z[2].append(zi)
            else: z[1].append(zi)
    if nlayer == 2:
        return (np.sum(z[1])/len(z[1])-np.sum(z[0])/len(z[0]))*c
    elif nlayer == 3:
        return [(np.sum(z[1])/len(z[1])-np.sum(z[0])/len(z[0]))*c, (np.sum(z[2])/len(z[2])-np.sum(z[1])/len(z[1]))*c]     


def create_shift(qe_input_file=None, qe_shifted_file=None, strain=0, e=0, nlayer=2, ntmp = 0, c=30):
    """
    Get a new shifted file from another input file

    Parameters:
        qe_input_file (str): Path to the QE input file.
        output_file (str): Path to save the new shifted input file
    """
    lattice_vectors = extract_cell_from_input(qe_input_file)
    if lattice_vectors is None:
            raise ValueError("Could not find lattice parameters in both QE output and input files.")

    # Extract atomic positions from output file
    atomic_positions = []
    element_order = []

    with open(qe_input_file, 'r') as f:
        lines = f.readlines()

    start_pos = None
    end_pos = None
    epos = None

    for i, line in enumerate(lines):
        if "  outdir      = './tmp'," in line:
            tmppos = i
        elif  "  eamp = 0," in line:
            epos = i
        elif "ATOMIC_POSITIONS (crystal)" in line:
            start_pos = i + 1
        elif start_pos and (line.strip() == "" or line.startswith("!END")):
            end_pos = i
            break

    if start_pos is None or end_pos is None:
        raise ValueError("Could not find atomic positions in the QE input file.")

    # Write to output file
    with open(qe_shifted_file, 'w') as f:
        if tmppos != None:
            for line in lines[0:tmppos]:
                f.write(line)
            f.write("  outdir      = './tmp" + str(ntmp) + "'" + "\n")

        if epos != None:
            for line in lines[tmppos+1:epos]:
                f.write(line)
            f.write("\teamp = " + str(e)+ "\n")
            for line in lines[epos+1:start_pos-1]:
                f.write(line)
   
                
        else:
            for line in lines[tmppos+1:start_pos-1]:
                f.write(line)

        f.write("ATOMIC_POSITIONS (crystal) \n")
        if nlayer == 2:
            d = interlayer_distance(qe_input_file, nlayer, c)
            shift = strain*d/c
            for line in lines[start_pos:end_pos]:
                parts = line.split()
                if float(parts[3]) < 0.5:
                    parts[3] = str(float(parts[3])-shift/2)
                else:
                    parts[3] = str(float(parts[3])+shift/2)
                f.write("\t".join(parts) + "\n")
        elif nlayer == 3:
             for line in lines[start_pos:end_pos]:
                parts = line.split()
                if float(parts[3]) < 0.45:
                    parts[3] = str(float(parts[3])-shift/40)
                elif float(parts[3]) > 0.55:
                    parts[3] = str(float(parts[3])+shift/40)
                f.write("\t".join(parts) + "\n")


        for line in lines[end_pos:len(lines)]:
            f.write(line)
    print(f"New shifted input saved to {qe_shifted_file}")

def findQ(qe_input_file, s, P, c, nlayer):
    d = interlayer_distance(qe_input_file, nlayer, c)
    s = np.array(s)
    P = np.array(P)
    P = P*(c/(d+s*d))
    P = P**2
    def func(data, a):
        return a*data
    qopt, qcov = curve_fit(func, P, s)
    plt.plot(P, s, "bo", label="DFT data")
    return(qopt, qcov)

def create_bash():
    return None

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create a new QE input file from old with changed interlayer distance.")
    parser.add_argument("qe_input_file", help="Path to the initial QE input file")
    parser.add_argument("qe_shifted_file", help="Path to the new shifted QE input file")
    args = parser.parse_args()
    create_shift(args.qe_output_file, args.qe_shifted_file, args.shift)