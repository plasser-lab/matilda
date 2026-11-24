#!/usr/bin/python3

"""
version 1.0.0
author: Felix Plasser
usage: Add a multiple of a normal mode to a structure.
"""

import sys
from matilda import matilda_header
from matilda import vib_molden # from nma
from matilda import struc_linalg

matilda_header.print_header('Normal mode displacement')

if len(sys.argv) < 6:
   print('At least 5 arguments required.')
   print('Syntax: add_normal_mode.py <in_struc> <out_struc> <vib_file> <nm_ind> <disp> [<type>]')
   sys.exit()
   
# read in data
in_file = sys.argv[1]
out_file = sys.argv[2]
vib_file = sys.argv[3]
nm_ind = float(sys.argv[4])
disp =   float(sys.argv[5]) # displacement
try:
    file_type = sys.argv[6]
except IndexError:
    file_type = 'xyz'

        
in_struc = struc_linalg.structure()
in_struc.read_file(file_path=in_file, file_type=file_type)

vmol = vib_molden.vib_molden()
vmol.read_molden_file(vib_file)

new_vec = in_struc.ret_vector() + disp * vmol.vibs[nm_ind - 1].ret_joined_vector()

out_struc = struc_linalg.structure()
out_struc.read_file_vector(def_file_path=in_file, file_type=file_type, vector=new_vec)
out_struc.make_coord_file(out_file, file_type=file_type)
