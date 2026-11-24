#!/usr/bin/env python3

import numpy as np 
import argparse
from matilda import matilda_header,struc_linalg,units

matilda_header.print_header('Geometries of lanthanide complexes')

parser = argparse.ArgumentParser(prog = 'rot_const',
                                description = 'Computes rotational constants (in MHz) of the molecule from geometry')
parser.add_argument('-f', '--filename' ,nargs='+', default = ['geom.xyz'])
parser.add_argument('--filetype', default = 'xyz')
args = parser.parse_args()

filename = args.filename
filetype = args.filetype

def rot_constant(filename):
    struc = struc_linalg.structure()
    struc.read_file(file_path = filename, file_type = filetype)

    I = struc.ret_moment_of_inertia() / units.constants['Nl'] * 1.E-23 #conversion of I from Da.A2 to kg.m2


    eigenvalues,_ = np.linalg.eigh(I)

    rot_const_MHz = units.constants['h'] / (8 * np.pi**2 * eigenvalues * 1e6)
    print("Rotational constants in MHz :", " ".join(f"{i:12.6f}" for i in rot_const_MHz))
for eachfile in filename:
    rot_constant(eachfile)
