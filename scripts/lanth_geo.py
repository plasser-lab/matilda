#!/usr/bin/python3

"""
This script calculates distances, and torsion angles to decide
if the given structure is SAP, borderline or tSAP.
"""
import argparse
from openbabel import openbabel
from matilda import matilda_header
from matilda import struc_linalg

def dist_print(struc, i, j):
    dist = struc.ret_bond_length(i, j)
    labi = "%s(%02i)"%(struc.ret_symbol(i), i)
    labj = "%s(%02i)"%(struc.ret_symbol(j), j)
    print("Distance %s--%s: % 4.3f"%(labi, labj, dist))

def tors_print(struc, a, b, c, d):
    laba = "%s (%02i)"%(struc.ret_symbol(a), a)
    labb = "%s (%02i)"%(struc.ret_symbol(b), b)
    labc = "%s (%02i)"%(struc.ret_symbol(c), c)
    labd = "%s (%02i)"%(struc.ret_symbol(d), d)
    tors = struc.ret_tors(a,b,c,d)
    print("Torsion %s, Aux, %s, %s: % 4.3f"%(laba, labc, labd, tors))
    return tors

if __name__=='__main__':
    matilda_header.print_header('Geometries of lanthanide complexes')

    parser = argparse.ArgumentParser(
                         prog='lanth_geo.py',
                         description='Determine SAP or tSAP geomtry for lanthanide complex')

    parser.add_argument('-f', '--filename', default='geom.xyz')
    parser.add_argument('--filetype', default='xyz')
    parser.add_argument('-l', '--lanth_ind', type = int, default=1, help='index defining the lanthanide')
    parser.add_argument('-b', '--bottom_inds', nargs = '+', type = int, default=[2, 3, 4, 5], help='indices defining the bottom square')
    parser.add_argument('-t', '--top_inds', nargs = '+', type = int, default=[6, 7, 8, 9], help='indices defining the top square')
    args = parser.parse_args()

    infile = args.filename
    intype = args.filetype
    lanth_ind = args.lanth_ind
    bottom_inds = args.bottom_inds
    top_inds = args.top_inds

    struc = struc_linalg.structure()
    struc.read_file(file_path=infile, file_type=intype)
    print(infile, bottom_inds, top_inds, lanth_ind)
    print()


    nb_inds = len(bottom_inds)
    nt_inds = len(top_inds)


    if nb_inds != nt_inds:
        raise ValueError("Number of top and bottom indices should be equal.")
    elif nb_inds == nt_inds not in [3,4]:
        raise ValueError("The number of top and bottom indices should be either 3 or 4 ")
    
    if nb_inds == 4:
        print("Analysing twisting between squares")
    elif nb_inds == 3:
        print("Analysing twisting between triangles")

    # Find the average position of the bottom square
    pos = []
    for i in range(nb_inds):
        pos.append(struc.ret_pos(bottom_inds[i]))
    pos_aux = sum(pos) / nb_inds

    obatom = openbabel.OBAtom()
    obatom.SetAtomicNum(1)
    obatom.SetVector(*pos_aux)
    struc.mol.AddAtom(obatom)

    aux_ind = struc.mol.NumAtoms()
    
    print()
    for i in range(nb_inds):
        dist_print(struc, lanth_ind, bottom_inds[i])
    for i in range(nt_inds):
        dist_print(struc, lanth_ind, top_inds[i])

    print()
    tsum = 0
    for i in range(nb_inds):
        tsum += tors_print(struc, bottom_inds[i], aux_ind, lanth_ind, top_inds[i])
    print()
    average_tors = tsum/nb_inds

    print("Average torsion: %4.3f"%(average_tors))

    print()
    if nb_inds == 4:
        if abs(average_tors) >= 35 :
            print("The complex is classified as SAP")
        elif 35 > abs(average_tors) > 25:
            print("The complex is borderline (SAP/tSAP)")
        elif abs(average_tors) <= 25 :
            print("The complex is lassified as tSAP")
        else:
            print("Geometry not defined")
        print()
