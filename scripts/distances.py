#!/usr/bin/python3

"""
version 1.0.0
author: Felix Plasser
usage: Print out a table with distances between structures.
They are superimposed onto the first structure.
"""
import sys
from matilda import struc_linalg, matilda_header

matilda_header.print_header('Print out a table with distances between structures')

if len(sys.argv) < 2:
   print('At least one argument required.')
   print('Syntax: python distances.py <struc1> <struc2> ... [-t<type>] [-mwp<mass_weight_power>] [-fit<fitting>] [-dig<digits>]')
   sys.exit()
   
# defaults
file_type = 'xyz'
mass_wt_pw = 1
fitting = 1
digits = 4

files = []
for arg in sys.argv[1:]:
    if arg[:4] == '-mwp':
        mass_wt_pw = eval(arg[4:])
    elif arg[:2] == '-t':
        file_type = arg[2:]
    elif arg[:4] == '-fit':
        fitting = arg[4:]
    elif arg[:4] == '-dig':
        digits = int(arg[4:])
    else:
        files += [arg]
        

struc0 = struc_linalg.structure(name=files[0][:4])
struc0.read_file(file_path=files[0], file_type=file_type)
out_strucs = [struc0]

for i,file in enumerate(files[1:]):
    struc = struc_linalg.structure(name=file[:digits+3])
    struc.read_file(file_path=file, file_type=file_type)
    # superimpose the structures
    if fitting == 1:
        print('- Superimposing the structures')
        out_strucs += [struc.ret_superimposed_structure(struc0, mass_wt_pw=mass_wt_pw)]
    else:
        out_strucs += [struc]

mc = struc_linalg.mol_calc(def_file_path=files[0], file_type=file_type)

print('Distances between structures')
print(mc.distance_table(out_strucs, mass_wt_pw=mass_wt_pw, digits=digits))


