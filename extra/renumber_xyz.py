#!/usr/bin/python

"""
version 1.0.0
author: Felix Plasser
usage: renumber an xyz file.
    renumber_xyz.py <in_file> <renumber_file> <out_file>
<renumber_file> contains the numbers of the atoms in the new order.
"""

import sys

if not len(sys.argv) == 4:
    print "Expected syntax:"
    print "python renumber_xyz.py <in_file> <renumber_file> <out_file>"
    sys.exit()
    
in_file = sys.argv[1]
ren_file = sys.argv[2]
out_file = sys.argv[3]
    
ren_list = [eval(num) for num in open(ren_file, 'r').readlines()] # contains the new order of the atom
in_list = open(in_file, 'r').readlines() # content of input file

# add numbers to the end of ren_list that were not entered into the renumber file
for i in xrange(1, len(in_list)-1):
    if not i in ren_list:
        ren_list += [i]

out_list = in_list[0:2]

for num in ren_list:
    out_list += [in_list[num+1]]
    
w_file = open(out_file, 'w')
w_file.writelines(out_list)
w_file.close()
