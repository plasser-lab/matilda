#!/usr/bin/python3

"""
version 1.0
author: Felix Plasser
usage: Print out the value of a specified internal coordinate.
"""
import argparse
import sys
from openbabel import pybel
from matilda import matilda_header, file_handler, struc_linalg

matilda_header.print_header('Evaluation of internal coordinates')

parser = argparse.ArgumentParser(
    prog='int_coor_multi.py',
    description='Outputs the value of the specified internal coordinate'
)
parser.add_argument('-f', '--filename', nargs = '+')
parser.add_argument('--filetype', default = 'xyz')
parser.add_argument('-d', '--dist', nargs = 2, type = int, action='append', help = 'The indices of atoms for distance')
parser.add_argument('-b', '--bend', nargs = 3, type = int, action='append', help = 'The indices of atoms for bend angle')
parser.add_argument('-t', '--tors', nargs = 4, type = int, action='append', help = 'The indices of atoms for dihedral angle')
parser.add_argument('-dg', '--digits',type = int,  default = 4, help = 'Number of decimal points')
args = parser.parse_args()

dist = args.dist
bend = args.bend
tors = args.tors
files = args.filename
file_type = args.filetype
digits = args.digits # output digits

coors = []
if args.dist:
    for d in args.dist:
        coors.append(['dist'] + d)
if args.bend:
    for b in args.bend:
        coors.append(['bend'] + b)
if args.tors:
    for t in args.tors:
        coors.append(['tors'] + t)

        
struc = struc_linalg.structure()

tm = file_handler.table_maker(col_widths=[20]+len(coors)*[5+digits])
tm.print_line(['File']+[coor[0] for coor in coors])
for file in files:
  for mol in pybel.readfile(file_type, file):
    struc.get_mol(mol=mol.OBMol,file_path=file, file_type=file_type)
    out_line = [file]
    for coor in coors:
        if coor[0] == 'dist':
            dist = struc.ret_bond_length(*coor[1:])
            out_line += ["%.*f"%(digits, dist)]
        elif coor[0] == 'bend':
            bend = struc.ret_bend(*coor[1:])
            out_line += ["%.*f"%(digits,bend)]
        elif coor[0] == 'tors':
            tors = struc.ret_tors(*coor[1:])
            out_line += ["%.*f"%(digits,tors)]
    tm.print_line(out_line)
