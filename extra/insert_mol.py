#!/usr/bin/python

print """
version 1.0.0
author: Felix Plasser (felix.plasser@unive.ac.at)
 University of Vienna, Institute for Theoretical Chemistry
 Waehringerstr. 17, 1090, Vienna, Austria

usage: insert a molecule into a larger template structure, replacing the corresponding atoms.
 The molecule inserted is translated and rotated but internal coordinates are kept
 without changing.
"""
try:
    import numpy
except:
    print "Error: third party package numpy not installed"
    print "  please install: e.g. 'yum install numpy'"
    raise

try:
    from matilda import file_handler, struc_linalg, superposition
except:
    print "Error in importing module of this distribution"
    print "check the PYTHONPATH environment variable"
    raise

if __name__=='__main__':

    ref_file = file_handler.def_input("Template structure file")
    ins_file = file_handler.def_input("File name of the molecule to insert")
    file_type = file_handler.def_input("File type", "xyz")

    ins = struc_linalg.structure(name=ins_file[:4])
    ins.read_file(file_path=ins_file, file_type=file_type)
    num_at = ins.ret_num_at()

    ref = struc_linalg.structure(name=ref_file[:4])
    ref.read_file(file_path=ref_file, file_type=file_type)

    print
    print "Number of atoms to insert: %i"%num_at
    print """
    Enter indices of the atoms in the template structure corresponding to the
    indices of the molecule inserted
    """

    at_list_ref = []
    for iat in xrange(1,num_at+1):
        try:
            at_list_ref.append(int(file_handler.def_input("Ind. %3i"%iat)))
        except:
            print "Error: Please enter one integer number!"
            raise

        symb_ref = ref.ret_symbol(at_list_ref[-1])
        symb_ins = ins.ret_symbol(iat)
        
        print "  Symbols: %s %s"%(symb_ref,symb_ins)
        if symb_ref != symb_ins:
            print "Warning: different atom types!"

    ins_mat = ins.ret_3xN_matrix()
    ref_mat = ref.ret_3xN_matrix(at_list = at_list_ref)
    mass_vect = ins.ret_mass_vector(power=1)

# Perform the superposition
    sup = superposition.superposition()
    sup.superimpose(ref_points=ref_mat, mv_points=ins_mat, weights=mass_vect)
    
    print
    sup.print_all_info()

    ins_mat = ins_mat - sup.ret_mv_av()
    ins_mat = numpy.dot(ins_mat, sup.ret_rotation_matrix().transpose())
    ins_mat = ins_mat + sup.ret_ref_av()

    print ins_mat

# update the coordinates in the reference structure and print out the new file
    ref.read_3xN_matrix(coor_mat=ins_mat, at_list=at_list_ref)
    out_file = ins_file.partition('.')[0] + '_INS_' + ref_file
    ref.make_coord_file(file_path=out_file,file_type=file_type)
    
    print "created "+out_file
    print "check results for example with:"
    print "  pymol %s %s"%(ref_file, out_file)