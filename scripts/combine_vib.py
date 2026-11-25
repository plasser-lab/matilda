#!/usr/bin/env python3

"""
version 1.0.0
author: Felix Plasser (felix.plasser@univie.ac.at)
 University of Vienna, Institute for Theoretical Chemistry
 Waehringerstr. 17, 1090, Vienna, Austria

usage: Insert one ore more structure files into a template and optionally combine the vibrations.
"""

# todo: more solid handling of sigma bond rotations
try:
    import numpy
except:
    print("Error: third party package numpy not installed")
    print("  please install: e.g. 'yum install numpy'")
    raise

try:
    from matilda import file_handler, struc_linalg, superposition, vib_molden, matilda_header
except:
    print("Error in importing module of this distribution")
    print("check the PYTHONPATH environment variable")
    raise

class template:
    """
    Class for handling the template where the structures are inserted
    """
    def __init__(self,file_name,file_type='xyz',vib=False):
        self.struc = struc_linalg.structure(name=file_name[:4])
        self.struc.read_file(file_path=file_name, file_type=file_type)
        self.num_at = self.struc.ret_num_at()
        self.file_name = file_name
        self.file_type = file_type
        
        self.vib = vib
        if self.vib:
            self.freqs = numpy.zeros([3*self.num_at])
            self.vib_mat = numpy.zeros([3*self.num_at,3*self.num_at])
            self.curr_vib = 0
        
        print()
        print("Template %s read in"%file_name)
        print("Number of atoms: %i"%self.num_at)
        
    def insert(self,ins,ky_file,ins_vib=None):
        """
        Insert a molecule into the structure.
        """
        print("""
        Enter indices of the atoms in the template structure corresponding to the
        indices of the molecule inserted.
        Enter '<index> [<weight>]', <weight> specifies a manual weighting factor. Default is the atomic mass.
        Enter '-' to include an atom which has no corresponding atom in the template.
        Enter 'x' to discard an atom of the inserted molecule.
        Enter '--' if all following atoms have no corresponding atom in the template.
        Enter 'xx' to discard all following atoms.
        """)
    
        at_list_ref = []
        at_list_ins = []
        not_list_ins = []
        mass_vect = []
        for iat in range(1,ins.num_at+1):
            
            ind_inp = ky_file.def_input("Ind. %3i"%iat)
            if ind_inp == '-':
                not_list_ins.append(iat)
            elif ind_inp == '--':
                for jat in range(iat, ins.num_at+1):
                    not_list_ins.append(jat)
                break
            elif ind_inp == '-':
                not_list_ins.append(iat)
            elif ind_inp == 'x':
                print('  - skipping atom')
            elif ind_inp == 'xx':
                print('  - skipping all following atoms')
                break
            else:
                inp_list = ind_inp.split()
                at_list_ref.append(int(inp_list[0]))
                
                at_list_ins.append(iat)
                
                symb_ref = self.struc.ret_symbol(at_list_ref[-1])
                symb_ins = ins.struc.ret_symbol(iat)                
                print("  Symbols: %s %s"%(symb_ref,symb_ins))
                if symb_ref != symb_ins:
                    print("Warning: different atom types!")
                # this has to be the inserted molecule !!
                #mass_vect.append(self.struc.mol.GetAtom(iat).GetAtomicMass())
                if len(inp_list) == 2:
                    weight = float(inp_list[1])
                    print("  - Using manual weighting factor % .3f"%weight)
                else:
                    weight = ins.struc.mol.GetAtom(iat).GetAtomicMass()
                mass_vect.append(weight)
    
        ins_mat = ins.struc.ret_3xN_matrix(at_list_ins)
        ref_mat = self.struc.ret_3xN_matrix(at_list = at_list_ref)
#        mass_vect = ins.struc.ret_mass_vector(power=1)
    
    # Perform the superposition
        print('mass_vect: ',mass_vect)
        sup = superposition.superposition()
        sup.superimpose(ref_points=ref_mat, mv_points=ins_mat, weights=mass_vect)
        
        print()
        sup.print_all_info(prt_array=False)
        print()
    
        ins_mat = ins_mat - sup.ret_mv_av()
        ins_mat = numpy.dot(ins_mat, sup.ret_rotation_matrix().transpose())
        ins_mat = ins_mat + sup.ret_ref_av()
        
        
        #print ins_mat
    
    # update the coordinates in the reference structure and print out the new file
        self.struc.read_3xN_matrix(coor_mat=ins_mat, at_list=at_list_ref)
        
    # add atoms which were not in the template, e.g. link hydrogens.
        if not_list_ins != []:
            not_mat = ins.struc.ret_3xN_matrix(not_list_ins)
        
            not_mat = not_mat - sup.ret_mv_av()
            not_mat = numpy.dot(not_mat, sup.ret_rotation_matrix().transpose())
            not_mat = not_mat + sup.ret_ref_av()
            
            for i,iat in enumerate(not_list_ins):
                AtNum = ins.struc.mol.GetAtom(iat).GetAtomicNum()
                a = self.struc.mol.NewAtom()
                a.SetAtomicNum(AtNum)
                a.SetVector(not_mat[i][0],not_mat[i][1],not_mat[i][2])
    
            not_list_ins = []
            
            if self.vib:
                print("Cannot combine vibrations and add atoms at the same time!")
                print("Please perform in two steps!")
                print("Only combining geometries now ...\n")
                self.vib = False
        
        if self.vib:
            print("Inserting vibrations")
            vib_m = vib_molden.vib_molden()
            vib_m.read_molden_file(file_path=ins_vib)
            
            for vib in vib_m.vibs:
                vec_array_rot = numpy.dot(vib.vector_list, sup.ret_rotation_matrix().transpose())
                
                for iat,iref in enumerate(at_list_ref):
                    for coor in xrange(3):
                        self.vib_mat[self.curr_vib,3*(iref-1)+coor]=vec_array_rot[iat,coor]
                        
                self.freqs[self.curr_vib] = vib.frequency
                self.curr_vib += 1
            
            
    def write_output(self,prefix='INS_'):
        """
        Write output files.
        """
        out_file = 'combined_struc.'+self.file_type #prefix+self.file_name
        self.struc.make_coord_file(file_path=out_file, file_type=self.file_type, lvprt=1)
        print("check results for example with:")
        print("  pymol %s %s"%(self.file_name, out_file))
        
        if self.vib:
            out_file = 'combined_vib.mld' #prefix+self.file_name+'_vib.mld'
            vib_molden.make_molden_file(struc=self.struc, freqs=self.freqs, vibs=self.vib_mat, out_file=out_file, title='combined vibrations', num_at=None)
            print()
            print("created "+out_file)
            print("check results for example with:")
            print("  molden %s"%(out_file))

class molecule:
    """
    Class for handling the separate molecules inserted.
    """
    def __init__(self,file_name,file_type='xyz'):
        self.struc = struc_linalg.structure(name=file_name[:4])
        self.struc.read_file(file_path=file_name, file_type=file_type)
        self.num_at = self.struc.ret_num_at()
        
        print()
        print("%s read in"%file_name)
        print("Number of atoms: %i"%self.num_at)


if __name__=='__main__':
    matilda_header.print_header('Insert structure files into a template')

    ky_file = file_handler.ky_file(file_name="combine_vib.in")
    tmpyn = ky_file.def_input("Include vibrations? (y/n)")
    if tmpyn == "y":
        vib = True
    else:
        vib = False
        
    ref_file = ky_file.def_input("Template structure file")
    file_type = ky_file.def_input("File type", "xyz")
    
    templ = template(file_name=ref_file, file_type=file_type,vib=vib)
    
    num_ins = int(ky_file.def_input("How many molecules do you want to insert?"))
    
    # loop
    for i in range(1,num_ins+1):
        ins_file = ky_file.def_input("Structure-file of the molecule to insert")
        ins_typ = ky_file.def_input("File type", "xyz")
        if vib:
            ins_vib = ky_file.def_input("Molden vibration file (same structure as the above file)")
        else:
            ins_vib = None
        ins = molecule(file_name=ins_file, file_type=ins_typ)
        templ.insert(ins=ins,ky_file=ky_file,ins_vib=ins_vib)
    
    templ.write_output()
    

    
    ky_file.close_write()
