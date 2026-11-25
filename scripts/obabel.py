#!/usr/bin/python3

"""
version 1.0.0
author: Felix Plasser, University of Vienna, Institute for Theoretical Chemistry
Waehringerstr. 17, 1090, Vienna, Austria
felix.plasser@univie.ac.at

usage: molecular structure conversion based on the openbabel package
     - increased support for Columbus (col) and Tinker (txyz2) formats

    Syntax: babel.py <intype> <infile> <outtype> <outfile>
"""

#import argparse - only available in version 2.7
import sys
from matilda import struc_linalg, matilda_header

matilda_header.print_header('Molecular structure conversion based on the openbabel package')

if len(sys.argv) < 4+1:
    print("""     Molecular structure conversion based on the openbabel package
     - increased support for Columbus (col) and Tinker (txyz2) formats

    Syntax: babel.py <intype> <infile> <outtype> <outfile> 
    
    Supported file types:
    col -- Columbus and Newton-X format
    txyz2 -- Tinker format with reading possibility (verify the atom types in the output)
  additionally all formats from openbabel are included
    type 'babel -H' for a complete list""")
    
    print('\nFour arguments required.')
    sys.exit()
    
(intype,infile,outtype,outfile) = sys.argv[1:]

struc = struc_linalg.structure()
struc.read_file(file_path=infile, file_type=intype)
struc.make_coord_file(file_path=outfile,file_type=outtype)

print("Finished. Converted %s -> %s"%(infile, outfile))
