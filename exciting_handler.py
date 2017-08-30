import numpy as np
import solid_state_tools as sst
import xml.etree.ElementTree as ET
import xml
from xml.dom import minidom
import subprocess
import periodictable as pt
import subprocess
import os
import time
import matplotlib.pyplot as plt
import re
plt.ion()

p_table = {i:el.__repr__() for i,el in enumerate(pt.elements)}
p_table_rev = {el.__repr__():i for i,el in enumerate(pt.elements)}

hartree = 27.211


def convert_greek(input):
  result = []
  for el in input:
    if el.lower() == 'gamma':
      result.append(r'$\Gamma$')
    else:
      result.append(el)
  return result


class Handler:

    def __init__(self):
        self.default_extension = '.xml'
        self.engine_command = ["excitingser"]
        self.working_dirctory = '/exciting_files/'
        self.engine_process = None
        self.info_file = 'INFO.OUT'
        self.project_directory = None
        self.input_filename = 'input.xml'
        self.exciting_folder = self.find_exciting_folder()
        self.scf_options = {'do':'fromscratch','nempty':'15','gmaxvr':'15','rgkmax': '5.0', 'ngridk': '5 5 5'}
        self.general_options = {'title':'title'}
        self.bs_options = {'steps':'300'}

    def find_exciting_folder(self):
        p = subprocess.Popen(['which','excitingser'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        res,err = p.communicate()
        res = res.split('bin')[0]
        return res.strip()


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
        crystal_structure = sst.CrystalStructure(crystal_base,atom_array)
        return crystal_structure

    def make_tree(self):
        root = ET.Element("input")
        tree = ET.ElementTree(root)
        ET.SubElement(root, "title").text=self.general_options['title']
        return tree

    def add_scf_to_tree(self,tree,crystal_structure):
        root = tree.getroot()
        structure = ET.SubElement(root, "structure",speciespath=self.exciting_folder+'species')
        crystal = ET.SubElement(structure,"crystal")
        for i in range(3):
            lattice_vector = crystal_structure.lattice_vectors[i,:]
            ET.SubElement(crystal, "basevect").text = "{0:1.6f} {1:1.6f} {2:1.6f}".format(*lattice_vector)

        abs_coord_atoms = crystal_structure.atoms
        species = set(abs_coord_atoms[:,3].astype(np.int))
        n_species = len(species)
        n_atoms = abs_coord_atoms.shape[0]

        for specie in species:
            species_mask = abs_coord_atoms[:,3].astype(np.int) == specie
            sub_coords = abs_coord_atoms[species_mask,:]
            specie_xml_el = ET.SubElement(structure, "species",speciesfile=p_table[specie]+'.xml')
            n_atoms_specie = sub_coords.shape[0]
            for i in range(n_atoms_specie):
                ET.SubElement(specie_xml_el, "atom",coord = "{0:1.6f} {1:1.6f} {2:1.6f}".format(*sub_coords[i,:]))

        groundstate = ET.SubElement(root, "groundstate", **self.scf_options)


    def add_bs_to_tree(self,tree,points):
        root = tree.getroot()
        properties = ET.SubElement(root, "properties")
        bandstructure = ET.SubElement(properties, "bandstructure")
        plot1d = ET.SubElement(bandstructure, "plot1d")
        path = ET.SubElement(plot1d, "path",**self.bs_options)
        for point in points:
            cords = point[0]
            label = point[1]
            ET.SubElement(path, "point", coord = "{0:1.6f} {1:1.6f} {2:1.6f}".format(*cords), label=label)


    def write_input_file(self,tree):
        if not os.path.isdir(self.project_directory+self.working_dirctory):
            os.mkdir(self.project_directory+self.working_dirctory)
        xmlstr = minidom.parseString(ET.tostring(tree.getroot())).toprettyxml(indent="   ")
        with open(self.project_directory+self.working_dirctory+self.input_filename, "w") as f:
            f.write(xmlstr)

        # tree.write(self.project_directory+self.working_dirctory+self.input_filename)

    def start_engine(self):
        os.chdir(self.project_directory+self.working_dirctory)
        command = self.engine_command
        self.engine_process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        os.chdir(self.project_directory)

    def start_ground_state_calculation(self):
        self.start_engine()

    def read_engine_status(self):
        if not self.is_engine_running():
            output, error = self.engine_process.communicate()
            return output,error
        else:
            return None

    def is_engine_running(self):
        if self.engine_process is None:
            return False
        if self.engine_process.poll() is None:
            return True
        else:
            return False

    def split_and_remove_whitespace(self,string):
        l1 = string.split()
        l2 = [float(x) for x in l1 if len(x)>0]
        return l2

    def find_bandgap(self, bands):
        for i in range(len(bands)):
            valence_band = bands[i]
            cond_band = bands[i + 1]
            if (np.max(cond_band[:, 1]) > 0) and (np.min(cond_band[:, 1]) < 0):
                return None, None
            if any(cond_band[:, 1] > 0):
                break
        # Direct bandgap
        band_diff = cond_band[:, 1] - valence_band[:, 1]
        bandgap_index = np.argmin(band_diff)
        bandgap = np.min(band_diff)
        k_bandgap = valence_band[bandgap_index, 0]
        # Indirect bandgap
        if np.abs((np.min(cond_band[:, 1]) - np.max(valence_band[:, 1]))  - bandgap) > 0.01:
            bandgap = (np.min(cond_band[:, 1]) - np.max(valence_band[:, 1]))
            k_bandgap = None
        return bandgap, k_bandgap

    def read_scf_status(self):
        try:
            f = open(self.project_directory+self.working_dirctory + self.info_file,'r')
        except IOError:
            return None
        info_text = f.read()
        f.close()
        scf_list = []
        matches = re.findall(r"SCF iteration number[\s\t]*:[\s\t]*\d*", info_text)
        for match in matches:
            ms = match.split(':')
            scf_list.append(int(ms[1]))

        scf_energy_list = []
        matches = re.findall(r"Total energy [\s\t]*:[\s\t]*[-+]?\d*\.\d+", info_text)
        for match in matches:
            ms = match.split(':')
            scf_energy_list.append(float(ms[1]))

        res = np.array(zip(scf_list,scf_energy_list))
        if len(res)<2:
            return None
        return res

    def read_bandstructure(self):
        e = xml.etree.ElementTree.parse(self.project_directory+self.working_dirctory+'bandstructure.xml').getroot()
        titlestring = e.find('title').itertext().next()

        bands = []

        for atype in e.findall('band'):
            band = []
            for point in atype.findall('point'):
                band.append([point.get('distance'), point.get('eval')])
            bands.append(np.array(band, dtype=np.float))

        special_k_points = []
        special_k_points_label = []
        for vertex in e.findall('vertex'):
            special_k_points.append(float(vertex.get('distance')))
            special_k_points_label.append(vertex.get('label'))

        special_k_points_label = convert_greek(special_k_points_label)

        empirical_correction = 0

        for band in bands:
            band[:,1] = band[:,1]*hartree
            if empirical_correction != 0:
                if (band[:, 1] > 0).any():
                    band[:, 1] = band[:, 1] + empirical_correction  / 2
                else:
                    band[:, 1] = band[:, 1] - empirical_correction  / 2

        bandgap, k_bandgap = self.find_bandgap(bands)
        special_k_points_together = zip(special_k_points,special_k_points_label)
        return sst.BandStructure(bands,bandgap,k_bandgap,special_k_points=special_k_points_together)

if __name__ == '__main__':
    os.chdir("/home/jannick/OpenDFT_projects/test/exciting_files")
    handler = Handler()
    handler.project_directory = "/home/jannick/OpenDFT_projects/test"

    print(handler.exciting_folder)
    tree = handler.make_tree()

    atoms = np.array([[0,0,0,6],[0.25,0.25,0.25,6]])
    unit_cell = 6.719*np.array([[0.5,0.5,0],[0.5,0,0.5],[0,0.5,0.5]])

    crystal_structure = sst.CrystalStructure(unit_cell,atoms)
    handler.add_scf_to_tree(tree,crystal_structure)
    handler.add_bs_to_tree(tree,[[np.array([0,0,0]),"GAMMA"]])

    handler.write_input_file(tree)
    # plt.figure(2)
    #
    # handler.start_engine()
    # while handler.is_engine_running():
    #     time.sleep(1)
    #     res = handler.read_scf_status()
    #     if res is not None and len(res) >1:
    #         print('updating plot')
    #         plt.clf()
    #         plt.plot(res[:,0],res[:,1])
    #         plt.pause(0.01)
    #         plt.draw()
    # print('finished')
    # band_structure = handler.read_bandstructure()
    # plt.figure(1)
    # for band in band_structure.bands:
    #     plt.plot(band[:,0],band[:,1],color='b',linewidth=2)
    # for xc,xl in band_structure.special_k_points:
    #     plt.axvline(x=xc, color='k', linewidth=1.5)
    #
    # unzipped_k = zip(*band_structure.special_k_points)
    # plt.xticks(unzipped_k[0], unzipped_k[1], rotation='horizontal')
