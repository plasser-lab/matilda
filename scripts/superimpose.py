#!/usr/bin/python3

"""
version 1.0.0
author: Felix Plasser, University of Vienna, Institute for Theoretical Chemistry
Waehringerstr. 17, 1090, Vienna, Austria
usage: Superimpose structures <struc1>, <struc2>, ... onto <ref_struc>.
    Syntax: python superimpose.py <ref_struc> <struc1> <struc2> ... -t<type>
"""

import sys
from matilda import struc_linalg, matilda_header

matilda_header.print_header('Superimpose structures')

if len(sys.argv) < 3:
   print('usage: Superimpose structures <struc1>, <struc2>, ... onto <ref_struc>.')
   print('At least two arguments required.\n')
   print('Syntax: python superimpose.py <ref_struc> <struc1> <struc2> ... [-t<type=xyz>] [-mwp<mass_weight_power=1>]')
   sys.exit()

# defaults
mass_wt_pw = 1
file_type = 'xyz'

ref_file = sys.argv[1]
super_files = [] # list with filenames to be superimposed onto ref_struc
for arg in sys.argv[2:]:
    if arg[:4] == '-mwp':
        mass_wt_pw = eval(arg[4:])
    elif arg[:2] == '-t':
        file_type = arg[2:]
    else:
        super_files += [arg]


# for calculations an instance of this class with a default molecule has to be defined
    # from it the numbering of the atoms and their masses are taken
mc = struc_linalg.mol_calc(def_file_path=ref_file, file_type=file_type)

# load the structures (DK and MK are random names used for identification)
ref = struc_linalg.structure(name=ref_file[:9])
ref.read_file(file_path=ref_file, file_type=file_type)

strucs = [] #initial structures and superimposed structures
for i,super_file in enumerate(super_files):
    strucs += [struc_linalg.structure(name=super_file[:9])]
    strucs[-1].read_file(file_path=super_file, file_type=file_type)
    # superimpose the structures
    strucs += [strucs[-1].ret_superimposed_structure(struc=ref, mass_wt_pw=mass_wt_pw, name='SI_%s'%super_file)]
    strucs[-1].make_coord_file('SI_' + super_file,file_type=file_type)

# to follow the effect of the superposition, a matrix with the distances can be printed out
print('The structures have been aligned. New files SI_*** were created')
print('RMSDs (Ang per atom) between structures')
print(mc.RMSD_table([ref]+strucs, mass_wt_pw=mass_wt_pw))
