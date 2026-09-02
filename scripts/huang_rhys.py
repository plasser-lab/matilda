#!/usr/bin/python3

from __future__ import print_function

"""
version 1.0.0
author: Felix Plasser, Sayan Ghosh
usage: Compute Huang-Rhys parameters and related quantities.
"""

import os, sys
import math
import numpy
import argparse

from matilda import vib_molden
from matilda import struc_linalg
from matilda import units


def compute_huang_rhys(gs_file, es_file, vib_file, file_type="xyz"):
    # Read in structures
    # TODO: One could also read gs_struc from gs vib_mld
    gs_struc = struc_linalg.structure()
    gs_struc.read_file(file_path=gs_file, file_type=file_type)

    es_struc = struc_linalg.structure()
    es_struc.read_file(file_path=es_file, file_type=file_type)
    # TODO: Could do superposition here

    diff_vec = es_struc.ret_vector() - gs_struc.ret_vector()
    print("Distance: %.4f Ang" % numpy.sum(diff_vec * diff_vec) ** 0.5)

    M = gs_struc.ret_mass_vector(power=0.5, rep=3)  # M_i = sqrt(m_i)
    mdiff_vec = diff_vec * M
    print("MW Dist: %.5f amu**.5 Ang"%(numpy.sum(mdiff_vec*mdiff_vec))**.5)
    print("MW MDist   : %.5f Ang"%(numpy.sum(mdiff_vec*mdiff_vec)/numpy.sum(M*M))**.5)
    mdiff_vec *= units.mass['amu']**(-.5) / units.length['A']

    mc = struc_linalg.mol_calc(gs_file, file_type=file_type)
    print("RMSD /3**.5: %.5f Ang"%(mc.RMSD(gs_struc, es_struc) / 3**.5))

    # Vibrational info
    vmol = vib_molden.vib_molden()
    vmol.read_molden_file(vib_file)
    Kmat = vmol.ret_vib_matrix()

    freqs_cm = vmol.ret_freqs()             # Frequencies in cm^-1
    Om = freqs_cm / units.energy['rcm']

    # TODO: put this into vib_molden
    for mode in Kmat:
        mode *= M
        norm = numpy.sum(mode * mode)
        if norm > 1E-8:
            mode *= norm ** (-0.5)

    dQ = Om ** 0.5 * numpy.dot(mdiff_vec, Kmat.T)
    S_factors = 0.5 * dQ ** 2

    # Filter out zero/imaginary frequencies
    valid_indices = freqs_cm > 1E-8
    freqs_cm = freqs_cm[valid_indices]
    Om = Om[valid_indices]
    dQ = dQ[valid_indices]
    S_factors = S_factors[valid_indices]

    reorg_energy_modes = S_factors * Om * units.energy['rcm']
    reorg_energy_total = numpy.sum(reorg_energy_modes)

    print("\nMode   Frequency(cm^-1)   dQ     Huang-Rhys S     Reorg_Energy (cm^-1)")
    for i in range(len(freqs_cm)):
        print(f"{i+1:3d}    {freqs_cm[i]:10.2f}   {dQ[i]:6.4f}   {S_factors[i]:10.6f}     {reorg_energy_modes[i]:12.2f}")

    print(f"\nTotal Reorganization Energy: {reorg_energy_total:.2f} cm^-1")

    return freqs_cm, S_factors, reorg_energy_modes, mdiff_vec, Kmat, Om, gs_struc


def show_top_modes(title, indices, freqs_cm, S_factors, reorg_energy_modes, N):
    print(f"\nTop {N} modes with {title}:")
    print(f"{'Mode':>4} {'Freq (cm^-1)':>15} {'S_i':>10} {'λ_i (cm^-1)':>12}")
    for idx in indices[:N]:
        print(f"{idx + 1:4d} {freqs_cm[idx]:15.2f} {S_factors[idx]:10.6f} {reorg_energy_modes[idx]:12.2f}")
    print()


def plot_spectrum(freqs_cm, S_factors, outfile):
    import matplotlib.pyplot as plt

    plt.figure(figsize=(8, 5))
    plt.bar(freqs_cm, S_factors, width=5, align='center', color='mediumblue', edgecolor='black', linewidth=0.7)
    plt.xlabel('Vibrational frequency (cm$^{-1}$)')
    plt.ylabel('Huang-Rhys factor $S_i$')
    plt.title('Huang-Rhys Spectrum')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig(outfile, dpi=300)
    print(f"Spectrum saved as '{outfile}'")


def save_data(freqs_cm, S_factors, outfile):
    with open(outfile, "w") as f:
        f.write("Mode\tFrequency(cm^-1)\tHuang-Rhys Factor (S_i)\n")
        for i in range(len(freqs_cm)):
            f.write(f"{i + 1}\t{freqs_cm[i]:.2f}\t{S_factors[i]:.6f}\n")
    print(f"Huang-Rhys data saved as '{outfile}' (further used for fc_mini.py)")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="huang_rhys.py",
        description="Compute Huang-Rhys parameters and reorganization energies."
    )

    parser.add_argument("gs_file", help="Ground state geometry file (e.g., gs.xyz)")
    parser.add_argument("es_file", help="Excited-state geometry file (e.g., es.xyz)")
    vib_file_group = parser.add_mutually_exclusive_group(required=True)
    vib_file_group.add_argument("-g", "--vib_gs", metavar="GS_VIB", help="Ground-state vibrational Molden file (default)")
    vib_file_group.add_argument("-e", "--vib_es", metavar="ES_VIB", help="Excited-state vibrational Molden file")

    parser.add_argument("-t", "--filetype", default="xyz", help="Filetype of structure files (default: xyz)")

    parser.add_argument("-s", "--topS", type=int, metavar="N", help="Show top N modes with highest Huang-Rhys factors")
    parser.add_argument("-l", "--topLambda", type=int, metavar="N", help="Show top N modes with highest reorganization energies")
    parser.add_argument("-p", "--plot", metavar="PNGFILE", help="Plot Huang-Rhys spectrum and save to file")

    args = parser.parse_args()

    if args.vib_gs:
        vib_file = args.vib_gs
    else:
        vib_file = args.vib_es
        print("NOTE: You are using excited-state vibrational modes (--vib_es).")

    freqs_cm, S_factors, reorg_energy_modes, mdiff_vec, Kmat, Om, gs_struc = compute_huang_rhys(
        gs_file=args.gs_file,
        es_file=args.es_file,
        vib_file=vib_file,
        file_type=args.filetype
    )

    if args.topS:
        top_indices = numpy.argsort(S_factors)[-args.topS:][::-1]
        show_top_modes("highest Huang-Rhys factors", top_indices, freqs_cm, S_factors, reorg_energy_modes, args.topS)

    if args.topLambda:
        top_reorg_indices = numpy.argsort(reorg_energy_modes)[-args.topLambda:][::-1]
        show_top_modes("highest reorganization energies", top_reorg_indices, freqs_cm, S_factors, reorg_energy_modes, args.topLambda)

    if args.plot:
        plot_spectrum(freqs_cm, S_factors, args.plot)

    save_data(freqs_cm, S_factors, "hr.txt")
