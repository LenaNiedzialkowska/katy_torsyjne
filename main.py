import os
import numpy as np
from Bio.PDB import PDBParser, calc_dihedral
from Bio.PDB import Polypeptide
import matplotlib.pyplot as plt

# Wczytanie struktury PDB
pdb_id = "1EHZ"
pdb_file = "1EHZ.pdb"
structure = PDBParser().get_structure(pdb_id, pdb_file)

def get_atom_vector(residue, atom_name):
    """Funkcja zwracająca wektor współrzędnych danego atomu w reszcie."""
    return residue[atom_name].get_vector() if atom_name in residue else None

def calculate_dihedral(v1, v2, v3, v4):
    """Funkcja obliczająca kąt dihedryczny pomiędzy 4 wektorami."""
    return calc_dihedral(v1, v2, v3, v4) if all((v1, v2, v3, v4)) else "N/A"

# Nagłówki do pliku wyjściowego
angles_names = ["alpha", "beta", "gamma", "delta", "epsilon", "zeta", "chi"]

with open("Output.txt", "w") as f:
    f.write("|".join(angles_names) + "\n")

    # Iteracja po resztach w strukturze
    residues = list(structure.get_residues())
    for i, residue in enumerate(residues):
        residue_name = residue.get_resname()

        # Pominięcie reszty HOH (woda)
        if residue_name == "HOH":
            continue

        # Pobranie atomów z poprzedniej i następnej reszty
        O3P = get_atom_vector(residues[i - 1], "O3'") if i > 0 else None
        O5N = get_atom_vector(residues[i + 1], "O5'") if i < len(residues) - 1 else None
        PN = get_atom_vector(residues[i + 1], "P") if i < len(residues) - 1 else None

        # Pobranie wektorów współrzędnych atomów
        P = get_atom_vector(residue, "P")
        O5 = get_atom_vector(residue, "O5'")
        C5 = get_atom_vector(residue, "C5'")
        C4 = get_atom_vector(residue, "C4'")
        C3 = get_atom_vector(residue, "C3'")
        O3 = get_atom_vector(residue, "O3'")
        C1 = get_atom_vector(residue, "C1'")
        O4 = get_atom_vector(residue, "O4'")
        C2C = get_atom_vector(residue, "C2")
        C4C = get_atom_vector(residue, "C4")

        N_name = "N9" if residue_name in ["A", "G"] else "N1"
        N = get_atom_vector(residue, N_name)

        # Obliczenie kątów dihedrycznych
        alpha = calculate_dihedral(O3P, P, O5, C5)
        beta = calculate_dihedral(P, O5, C5, C4)
        gamma = calculate_dihedral(O5, C5, C4, C3)
        delta = calculate_dihedral(C5, C4, C3, O3)
        epsilon = calculate_dihedral(C4, C3, O3, PN)
        zeta = calculate_dihedral(C3, O3, PN, O5N)
        chi = calculate_dihedral(O4, C1, N, C4C) if residue_name in ["A", "G"] else calculate_dihedral(O4, C1, N, C2C)

        # Zapis kątów do pliku
        angles = [alpha, beta, gamma, delta, epsilon, zeta, chi]
        f.write(" ".join(map(str, angles)) + "\n")

