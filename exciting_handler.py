import numpy as np
import solid_state_tools as sst
import xml.etree.ElementTree as ET

import periodictable as pt

p_table = {i:el.__repr__() for i,el in enumerate(pt.elements)}
p_table_rev = {el.__repr__():i for i,el in enumerate(pt.elements)}


class Handler:

    def parse_input_file(self,filename):
        tree = ET.parse(filename)
        root = tree.getroot()

        structure = root.find('structure')
        crystal = structure.find('crystal')
        species_list = structure.findall('species')


        cart_bool = structure.get('cartesian')
        if cart_bool is None:
            cart_bool = False
        else:
            cart_bool = bool(cart_bool)

        basevecs = crystal.findall('basevect')
        try:
            scale = float(crystal.get('scale'))
        except Exception:
            scale = 1
        crystal_base = np.zeros((3, 3))

        for i, basevec in enumerate(basevecs):
            coords = np.array(self.split_and_remove_whitespace(basevec.text))
            crystal_base[i, :] = coords*scale

        atomic_cord_list = []
        for species in species_list:
            species_name = species.get('speciesfile').split('.')[0]
            species_number = p_table_rev[species_name]
            for j1, atom in enumerate(species):
                pos_vec=np.zeros(4)
                rel_coords = np.array(self.split_and_remove_whitespace(atom.get('coord')))
                if cart_bool:
                    pos_vec[:3] = np.dot(np.linalg.inv(crystal_base.T),rel_coords)
                else:
                    pos_vec[:3] = rel_coords
                pos_vec[3] = species_number
                atomic_cord_list.append(pos_vec)

        atom_array = np.array(atomic_cord_list)
        print(atom_array)
        print(crystal_base)
        crystal_structure = sst.CrystalStructure(crystal_base,atom_array)
        return crystal_structure

    def __init__(self):
        self.default_extension = '.xml'


    def split_and_remove_whitespace(self,string):
        l1 = string.split()
        l2 = [float(x) for x in l1 if len(x)>0]
        return l2

if __name__ == '__main__':
    handler = Handler()
    filename = 'testfiles/exciting_test.xml'
    handler.parse_input_file('testfiles/exciting_test.xml')



