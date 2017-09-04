from __future__ import division
from __future__ import generators

import os
os.environ['ETS_TOOLKIT'] = 'qt4'
from pyface.qt import QtGui, QtCore
from traits.api import HasTraits, Instance, on_trait_change, Range, Bool, Button
from traitsui.api import View, Item, Group
from mayavi.core.ui.api import MayaviScene, MlabSceneModel, \
    SceneEditor

import numpy as np
import matplotlib as mpl
mpl.use('Qt4Agg')
mpl.rcParams['backend.qt4']='PySide'

# mpl.rc('font',**{'size': 22, 'family':'serif','serif':['Palatino']})
# mpl.rc('text', usetex=True)

mpl.rc('font',**{'size': 22})


from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
import matplotlib.pyplot as plt
import random


bohr = 0.52917721

cov_radii = np.loadtxt('./data/cov_radii.dat')/bohr

colors = {1: (0.8, 0.8, 0.8),3:(0,0.75,0.75),11:(0,0.75,0.75),19:(0,0.75,0.75),37:(0,0.75,0.75),5:(0.78,0.329,0.1176),7:(0,0,1),
          6:(0.25,.25,.25), 8: (1, 0, 0),9:(0,1,0),17:(0,1,0),35:(0,1,0),16:(1,1,0),13:(0.68,0.229,0.1176),31:(0.58,0.15,0.07),15:(0,0,0.8),33:(48/255,139/255,229/255)}

for i in range(21,31):
    colors[i] = (0,1,1)

for i in range(39,49):
    colors[i] = (0,1,1)

for i in range(71,80):
    colors[i] = (0,1,1)



def KnuthMorrisPratt(text, pattern):

    '''Yields all starting positions of copies of the pattern in the text.
Calling conventions are similar to string.find, but its arguments can be
lists or iterators, not just strings, it returns all matches, not just
the first one, and it does not need the whole text in memory at once.
Whenever it yields, it will have read the text exactly up to and including
the match that caused the yield.'''

    # allow indexing into pattern and protect against change during yield
    pattern = list(pattern)

    # build table of shift amounts
    shifts = [1] * (len(pattern) + 1)
    shift = 1
    for pos in range(len(pattern)):
        while shift <= pos and pattern[pos] != pattern[pos-shift]:
            shift += shifts[pos-shift]
        shifts[pos+1] = shift

    # do the actual search
    startPos = 0
    matchLen = 0
    for c in text:
        while matchLen == len(pattern) or \
              matchLen >= 0 and pattern[matchLen] != c:
            startPos += shifts[matchLen]
            matchLen -= shifts[matchLen]
        matchLen += 1
        if matchLen == len(pattern):
            yield startPos

class StructureVisualization(HasTraits):
    n_x    = Range(1, 10, 1, mode='spinner')#)
    n_y  = Range(1, 10, 1,  mode='spinner')#mode='spinner')
    n_z  = Range(1, 10, 1,  mode='spinner')#mode='spinner')
    prop_but = Button(label='Properties')

    show_unitcell = Bool(True)
    show_bonds = Bool(True)

    view = View(Item('scene', editor=SceneEditor(scene_class=MayaviScene),
                     height=450, width=500, show_label=False),Group('_', 'n_x', 'n_y', 'n_z','show_unitcell','show_bonds',orientation='horizontal'),
                resizable=True, # We need this to resize with the parent widget
                )
    scene = Instance(MlabSceneModel, ())

    def __init__(self, crystal_structure):
        super(StructureVisualization,self).__init__()
        self.crystal_structure = crystal_structure

    def clear_plot(self):
        self.scene.mlab.clf()

    @on_trait_change('scene.activated,show_unitcell,show_bonds,n_x,n_y,n_z')
    def update_plot(self):
        if self.crystal_structure is None:
            return
        self.scene.anti_aliasing_frames = 20
        # We can do normal mlab calls on the embedded scene.
        self.scene.mlab.clf()
        repeat = [self.n_x,self.n_y,self.n_z]

        self.plot_atoms(repeat=repeat)

        if self.show_bonds:
            self.plot_bonds(repeat=repeat)

        if self.show_unitcell:
            self.plot_unit_cell(repeat=repeat)

    def check_if_line_exists(self,p1,p2,list_of_lines):
        for line in list_of_lines:
            x1 = line[0]
            x2 = line[1]
            norm1 = np.linalg.norm(p1-x1)+np.linalg.norm(p2-x2)
            norm2 = np.linalg.norm(p1-x2)+np.linalg.norm(p2-x1)
            if norm1 < 0.05 or norm2 < 0.05:
                return True
        return False

    def plot_unit_cell(self, repeat=[1, 1, 1]):
        cell = self.crystal_structure.lattice_vectors
        existing_lines = []

        for j1 in range(repeat[0]):
            for j2 in range(repeat[1]):
                for j3 in range(repeat[2]):
                    offset = j1 * cell[0, :] + j2 * cell[1, :] + j3 * cell[2, :]
                    self.plot_single_unit_cell(offset)
                    # for i1, a in enumerate(cell):
                    #     i2 = (i1 + 1) % 3
                    #     i3 = (i1 + 2) % 3
                    #     for b in [np.zeros(3), cell[i2]]:
                    #         for c in [np.zeros(3), cell[i3]]:
                    #             p1 = b + c + offset
                    #             p2 = p1 + a
                    #             if not self.check_if_line_exists(p1,p2,existing_lines):
                    #                 self.scene.mlab.plot3d([p1[0], p2[0]], [p1[1], p2[1]], [p1[2], p2[2]], tube_radius=0.03)
                    #                 existing_lines.append([p1,p2])

    def plot_single_unit_cell(self,offset):
        cell = self.crystal_structure.lattice_vectors
        a1 = cell[0,:]
        a2 = cell[1,:]
        a3 = cell[2,:]

        p1 = offset
        p2 = offset+a1
        p3 = offset+a2
        p4 = offset+a3
        p5 = offset+a1+a2
        p6 = offset+a1+a3
        p7 = offset+a2+a3
        p8 = offset+a1+a2+a3

        coords = np.array([p1,p2,p5,p3,p1,p4,p6,p8,p7,p4,p6,p2,p5,p8,p7,p3])
        self.scene.mlab.plot3d(coords[:,0],coords[:,1],coords[:,2], tube_radius=0.05)

    def plot_atoms(self, repeat=[1, 1, 1]):
        abs_coord_atoms = self.crystal_structure.calc_absolute_coordinates(repeat=repeat)
        species = set(abs_coord_atoms[:,3].astype(np.int))
        n_species = len(species)
        n_atoms = abs_coord_atoms.shape[0]



        for specie in species:
            species_mask = abs_coord_atoms[:,3].astype(np.int) == specie
            sub_coords = abs_coord_atoms[species_mask,:]

            cov_radius = cov_radii[specie]
            atom_size = 0.4*np.log(specie)+0.6
            try:
                atomic_color = colors[specie]
            except KeyError:
                atomic_color = (0.8,0.8,0.8)
            self.scene.mlab.points3d(sub_coords[:,0],sub_coords[:,1],sub_coords[:,2],
                                     scale_factor=atom_size,
                                     resolution=30,
                                     color=atomic_color)

    def plot_density(self,ks_density):
        cam, foc = self.scene.mlab.move()
        print(cam,foc)
        unit_cell = self.crystal_structure.lattice_vectors
        cp = self.scene.mlab.contour3d(ks_density.density, contours=10, transparent=True,
                            opacity=0.5, colormap='hot')
        # Do some tvtk magic in order to allow for non-orthogonal unit cells:
        polydata = cp.actor.actors[0].mapper.input
        pts = np.array(polydata.points) - 1
        # Transform the points to the unit cell:
        polydata.points = np.dot(pts, unit_cell / np.array(ks_density.density.shape)[:, np.newaxis])
        cam_new, foc_new = self.scene.mlab.move()
        print(cam_new,foc_new)
        dist = (foc - foc_new)
        self.scene.mlab.move()
        self.scene.mlab.move(-dist[2],dist[1],dist[0])
        cam, foc = self.scene.mlab.move()

        print(cam,foc)

        self.scene.mlab.show()


    # def bonds_to_paths(self,bonds):
    # """See here my failed attempt on graph theory. Seems there is a reason that there is a mathematical field to it. Who would have thought?"""
    #     def is_bond_in_path(path,bond):
    #         a = KnuthMorrisPratt(path,bond)
    #         b = KnuthMorrisPratt(path,bond[::-1])
    #         try:
    #             next(a)
    #             return True
    #         except StopIteration:
    #             pass
    #
    #         try:
    #             next(b)
    #             return True
    #         except StopIteration:
    #             pass
    #         return False
    #
    #     def find_all_connection(bonds, point):
    #         connections = []
    #         for bond in bonds:
    #             if point in bond:
    #                 if point == bond[0]:
    #                     connections.append(bond[1])
    #                 else:
    #                     connections.append(bond[0])
    #         return connections
    #
    #     paths = []
    #     for bond in bonds:
    #         path = bond
    #         p1 = bond[0]
    #         while True:
    #             p1_con = find_all_connection(bonds,p1)
    #             if len(p1_con) == 1:
    #                 break
    #             for connection in p1_con:
    #                 if not is_bond_in_path(path,[p1,connection]):
    #                     path.append[p1_con[0]]





    def plot_bonds(self,repeat=[1,1,1]):
        abs_coord_atoms = self.crystal_structure.calc_absolute_coordinates(repeat=repeat)
        bonds = self.crystal_structure.find_bonds(abs_coord_atoms)

        for bond in bonds:
            i1 = bond[0]
            i2 = bond[1]
            p1 = abs_coord_atoms[i1,:3]
            p2 = abs_coord_atoms[i2,:3]
            self.scene.mlab.plot3d([p1[0], p2[0]], [p1[1], p2[1]], [p1[2], p2[2]], tube_radius=0.125)

class OpticalSpectrumVisualization(QtGui.QWidget):
    def __init__(self, parent=None):
        super(OpticalSpectrumVisualization, self).__init__()
        self.first_plot_bool = True
        # a figure instance to plot on
        self.figure = plt.figure(1)
        plt.close(plt.figure(1))
        self.ax = None

        self.canvas = FigureCanvas(self.figure)

        self.toolbar = NavigationToolbar(self.canvas, self)

        color = self.palette().color(QtGui.QPalette.Base)
        self.figure.patch.set_facecolor([color.red()/255,color.green()/255,color.blue()/255])
        # self.figure.patch.set_alpha(1.0)
        # self.figure.patch.set_facecolor('blue')

        layout = QtGui.QVBoxLayout()
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)
        # layout.addWidget(self.button)
        self.setLayout(layout)
        self.show()

    def clear_plot(self):
        if not self.first_plot_bool:
            self.figure.clf()
            self.first_plot_bool = True
            self.canvas.draw()

    def plot(self,optical_spectrum):
        if self.first_plot_bool:
            self.ax = self.figure.add_subplot(111)
        self.ax.cla()

        self.ax.plot(optical_spectrum.energy,optical_spectrum.epsilon2,linewidth=2,color='k')

        if self.first_plot_bool:
            self.first_plot_bool = False
            self.figure.tight_layout()
        self.canvas.draw()


class BandStructureVisualization(QtGui.QWidget):
    def __init__(self, parent=None):
        super(BandStructureVisualization, self).__init__()
        self.first_plot_bool = True
        # a figure instance to plot on
        self.figure = plt.figure(1)
        plt.close(plt.figure(1))
        self.ax = None

        self.canvas = FigureCanvas(self.figure)

        self.toolbar = NavigationToolbar(self.canvas, self)

        color = self.palette().color(QtGui.QPalette.Base)
        self.figure.patch.set_facecolor([color.red()/255,color.green()/255,color.blue()/255])
        # self.figure.patch.set_alpha(1.0)
        # self.figure.patch.set_facecolor('blue')

        layout = QtGui.QVBoxLayout()
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)
        # layout.addWidget(self.button)
        self.setLayout(layout)
        self.show()

    def clear_plot(self):
        if not self.first_plot_bool:
            self.ax.cla()
            self.canvas.draw()

    def plot(self,band_structure):
        if self.first_plot_bool:
            self.ax = self.figure.add_subplot(111)
        self.ax.cla()
        for band in band_structure.bands:
            self.ax.plot(band[:,0],band[:,1],color='b',linewidth=2)
        for xc,xl in band_structure.special_k_points:
            self.ax.axvline(x=xc, color='k', linewidth=1.5)

        unzipped_k = zip(*band_structure.special_k_points)
        bandgap = band_structure.bandgap
        k_bandgap = band_structure.k_bandgap
        special_k_points = unzipped_k[0]
        special_k_points_label = unzipped_k[1]

        xlength = band[:,0].max()-band[:,0].min()
        self.ax.set_xlim(band[:,0].min()-xlength/800,band[:,0].max())
        self.ax.plot([band[:,0].min(),band[:,0].max()], [0, 0], 'k--')

        self.ax.set_xticks(special_k_points)
        self.ax.set_xticklabels(special_k_points_label,rotation='horizontal',horizontalalignment='center')

        if k_bandgap is None and bandgap is not None:
            title_bandgap = ' $E_g = %1.1f $ eV' % bandgap + ' (indirect)'
        elif bandgap is None:
            title_bandgap = ' (metallic)'
        elif k_bandgap in special_k_points:
            k_bandgap_label = np.array(special_k_points_label)[k_bandgap == special_k_points][0]
            title_bandgap = ' $E_g$ = %1.1f eV' % bandgap + ' at ' + k_bandgap_label
        else:
            title_bandgap = ' $E_g$ = %1.1f eV' % bandgap


        self.ax.set_title('KS bandstructure,' + title_bandgap, fontsize=25)
        if self.first_plot_bool:
            self.first_plot_bool = False
            self.figure.tight_layout()
        self.canvas.draw()


class ScfVisualization(QtGui.QWidget):
    def __init__(self, parent=None):
        super(ScfVisualization, self).__init__()

        # a figure instance to plot on
        self.figure = plt.figure(2)
        plt.close(plt.figure(2))
        self.ax = None
        self.first_plot_bool = True
        self.canvas = FigureCanvas(self.figure)

        color = self.palette().color(QtGui.QPalette.Base)
        self.figure.patch.set_facecolor([color.red()/256,color.green()/256,color.blue()/256])
        # self.figure.patch.set_alpha(0.0)

        self.toolbar = NavigationToolbar(self.canvas, self)

        # Just some button connected to `plot` method
        # self.button = QtGui.QPushButton('Plot')
        # self.button.clicked.connect(self.plot)

        # set the layout
        layout = QtGui.QVBoxLayout()
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)
        # layout.addWidget(self.button)
        self.setLayout(layout)
        self.show()

    def clear_plot(self):
        if not self.first_plot_bool:
            self.ax.cla()
            self.canvas.draw()

    def plot(self,scf_data):
        if self.first_plot_bool:
            self.ax = self.figure.add_subplot(111)
        self.ax.cla()
        self.ax.plot(scf_data[:,0], scf_data[:,1],'o-',linewidth=2)
        self.ax.set_xlim(1,scf_data[:,0].max())
        self.ax.set_xlabel('Scf iteration')
        self.ax.set_ylabel('Total Energy')
        if self.first_plot_bool:
            self.first_plot_bool = False
            self.figure.tight_layout()
        self.canvas.draw()


if __name__ == "__main__":
    vis = StructureVisualization(None)
