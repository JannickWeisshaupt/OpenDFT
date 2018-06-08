# -*- coding: utf-8 -*-
from __future__ import division
from __future__ import generators

import os

os.environ['ETS_TOOLKIT'] = 'qt4'

from pyface.qt import QtGui, QtCore
from traits.api import HasTraits, Instance, on_trait_change, Range, Bool, Button, Array
from traitsui.api import View, Item, Group, HGroup

from mayavi import mlab
from mayavi.core.ui.api import MayaviScene, MlabSceneModel, \
    SceneEditor
from mayavi.core.api import Engine, PipelineBase, Source

from tvtk.api import tvtk
from tvtk.pyface.scene import Scene
from tvtk.tools import visual

import solid_state_tools as sst
import copy
import numpy as np
import matplotlib as mpl

mpl.use('Qt4Agg')

# mpl.rcParams['backend.qt4']='PySide'
from little_helpers import find_data_file, find_possible_sums, find_grid_connections
from bisect import bisect

# mpl.rc('font',**{'size': 22, 'family':'serif','serif':['Palatino']})
# mpl.rc('text', usetex=True)

mpl.rc('font', **{'size': 22})

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
import matplotlib.pyplot as plt
import random

bohr = 0.52917721

cov_radii = np.loadtxt(find_data_file('/data/cov_radii.dat')) / bohr

colors = {1: (0.8, 0.8, 0.8), 3: (0, 0.75, 0.75), 11: (0, 0.75, 0.75), 19: (0, 0.75, 0.75), 37: (0, 0.75, 0.75),
          5: (0.78, 0.329, 0.1176), 7: (0, 0, 1),
          6: (0.25, .25, .25), 8: (1, 0, 0), 9: (0, 1, 0), 17: (0, 1, 0), 35: (0, 1, 0), 16: (1, 1, 0),
          13: (0.68, 0.229, 0.1176), 31: (0.58, 0.15, 0.07), 15: (0, 0, 0.8), 33: (48 / 255, 139 / 255, 229 / 255)}

for i in range(21, 31):
    colors[i] = (0, 1, 1)

for i in range(39, 49):
    colors[i] = (0, 1, 1)

for i in range(71, 80):
    colors[i] = (0, 1, 1)

lut = np.zeros((255, 4))  # colormap for atoms
lut[:, :] = np.array([100, 100, 100, 255])
for key, color in colors.items():
    color_255 = [255 * x for x in color]
    opacity = 255
    lut[key - 1, :] = np.array(color_255 + [opacity, ])

try:
    with open(find_data_file('/data/colormaps.dat')) as f:
        t = f.read()
        t = t.replace("'", '')
        s = t.split()
except IOError:
    s = ['hot', 'viridis', 'jet']

colormap_list = sorted(s, key=str.lower)


def convert_to_greek(input):
    result = []
    for el in input:
        if el.lower().strip() == 'gamma':
            result.append(r'$\Gamma$')
        else:
            result.append(el)
    return result


class BrillouinVisualization(HasTraits):
    view = View(Item('scene', editor=SceneEditor(scene_class=MayaviScene),
                     height=450, width=500, show_label=False), Group('_', orientation='horizontal'),
                resizable=True,  # We need this to resize with the parent widget
                )

    scene = Instance(MlabSceneModel, ())

    def __init__(self, parent):
        super(BrillouinVisualization, self).__init__(parent=parent)
        self.parent = parent
        self.crystal_structure = None
        self.k_path = None
        self.brillouin_edges = None
        self.path_plot = None
        self.plot_of_vertices = None
        self.text_plots = []
        self.glyph_points = None
        self.picker = None

    def clear_plot(self):
        self.scene.mlab.clf(figure=self.scene.mayavi_scene)

    def set_crystal_structure(self, crystal_structure):
        self.crystal_structure = crystal_structure
        self.w_points = sst.construct_brillouin_vertices(crystal_structure)
        self.brillouin_edges = sst.construct_convex_hull(self.w_points)

    def set_path(self, k_path):
        self.k_path = k_path

    @on_trait_change('scene.activated')
    def update_plot(self, *args, **kwargs):
        if self.crystal_structure is None:
            return

        self.scene.mlab.clf(figure=self.scene.mayavi_scene)

        if self.k_path is not None:
            self.plot_path()

        self.plot_brillouin_zone()
        self.picker = self.scene.mayavi_scene.on_mouse_pick(self.picker_callback)
        self.picker.tolerance = 0.01

    def plot_unit_vectors(self):
        pass
        # for i in range(3):
        #     Arrow_From_A_to_B(0,0,0,*self.crystal_structure.inv_lattice_vectors[i,:],figure=self.scene.mayavi_scene)

    def plot_path(self):
        if self.path_plot is not None:
            self.path_plot.remove()
            for text_plot in self.text_plots:
                text_plot.remove()

        if len(self.k_path) == 0:
            self.path_plot = None
            return

        n_path = len(self.k_path)
        k_path_array = np.zeros((n_path, 3))

        for i in range(n_path):
            k_path_array[i, :] = np.dot(self.crystal_structure.inv_lattice_vectors.T, self.k_path[i][0])

        self.path_plot = self.scene.mlab.plot3d(k_path_array[:, 0], k_path_array[:, 1], k_path_array[:, 2],
                                                color=(0, 1, 0), reset_zoom=False, tube_radius=0.02,
                                                figure=self.scene.mayavi_scene)
        # self.scene.mlab.points3d(k_path_array[[0,-1],0],k_path_array[[0,-1],1],k_path_array[[0,-1],2], scale_factor=.1,reset_zoom=False, figure=self.scene.mayavi_scene)

        labels = [point[1] for point in self.k_path]
        self.text_plots = [None] * n_path
        for i in range(n_path):
            text_plot = self.scene.mlab.text3d(k_path_array[i, 0], k_path_array[i, 1], k_path_array[i, 2], labels[i],
                                               scale=0.1, figure=self.scene.mayavi_scene)
            self.text_plots[i] = text_plot

    def plot_brillouin_zone(self, plot_connections=True):

        self.wpoints_plot = np.append(self.w_points, np.array([[0, 0, 0]]), axis=0)
        # TODO find general way to make the center of facets points
        # for i in [-1,1]:
        #     for j in range(3):
        #         self.wpoints_plot = np.append(self.wpoints_plot, np.array([i*0.5*self.crystal_structure.inv_lattice_vectors[j,:]]), axis=0)

        # self.scene.mlab.points3d(0.0, 0.0, 0.0, color=(0.7, 0.7, 0.7), scale_factor=.1, figure=self.scene.mayavi_scene)
        self.plot_of_vertices = self.scene.mlab.points3d(self.wpoints_plot[:, 0], self.wpoints_plot[:, 1],
                                                         self.wpoints_plot[:, 2], color=(0.7, 0.7, 0.7),
                                                         scale_factor=.1, figure=self.scene.mayavi_scene)
        self.glyph_points = self.plot_of_vertices.glyph.glyph_source.glyph_source.output.points.to_array()

        self.scene.mlab.triangular_mesh(self.w_points[:, 0], self.w_points[:, 1], self.w_points[:, 2],
                                        self.brillouin_edges, opacity=0.3, color=(0.5, 0.5, 0.5), tube_radius=2,
                                        figure=self.scene.mayavi_scene)

        # if plot_connections:
        #     for i, connection in enumerate(self.brillouin_edges):
        #         for con in connection:
        #             bond = [i, con]
        #             self.scene.mlab.plot3d(self.w_points[bond, 0], self.w_points[bond, 1], self.w_points[bond, 2],figure=self.scene.mayavi_scene,tube_radius=0.01)

        self.plot_unit_vectors()

        self.outline = self.scene.mlab.outline(line_width=3, figure=self.scene.mayavi_scene)
        self.outline.outline_mode = 'cornered'

        self.outline.bounds = (- 0.001, + 0.001, - 0.001, + 0.001, - 0.001, + 0.001)

    def picker_callback(self, picker):
        """ Picker callback: this get called when on pick events.
        """
        if picker.actor in self.plot_of_vertices.actor.actors:
            # Find which data point corresponds to the point picked:
            # we have to account for the fact that each data point is
            # represented by a glyph with several points
            point_id = int(picker.point_id / self.glyph_points.shape[0])
            # If the no points have been selected, we have '-1'
            if point_id != -1:
                # Retrieve the coordinnates coorresponding to that data
                # point
                x, y, z = self.wpoints_plot[point_id, :]
                # Move the outline to the data point.
                self.outline.bounds = (x - 0.03, x + 0.03,
                                       y - 0.03, y + 0.03,
                                       z - 0.03, z + 0.03)
                k_point = np.array([x, y, z])
                k_point_conv = np.dot(np.linalg.inv(self.crystal_structure.inv_lattice_vectors.T), k_point)
                self.k_path.append([k_point_conv, ''])
                # self.update_plot()
                self.plot_path()
                self.parent.set_path(self.k_path)


class StructureVisualization(HasTraits):
    n_x = Range(1, 20, 1, mode='spinner')  # )
    n_y = Range(1, 20, 1, mode='spinner')  # mode='spinner')
    n_z = Range(1, 20, 1, mode='spinner')  # mode='spinner')
    prop_but = Button(label='Properties')

    show_unitcell = Bool(True)
    show_bonds = Bool(True)
    show_atoms = Bool(True)
    edge_atoms = Bool(False)

    view = View(Item('scene', editor=SceneEditor(scene_class=MayaviScene),
                     height=450, width=500, show_label=False),
                Group('_', 'n_x', 'n_y', 'n_z', 'show_unitcell', 'show_bonds', 'show_atoms','edge_atoms', orientation='horizontal'),
                resizable=True,  # We need this to resize with the parent widget
                )

    scene = Instance(MlabSceneModel, ())

    def __init__(self, crystal_structure):
        super(StructureVisualization, self).__init__()
        self.crystal_structure = crystal_structure
        self.density_plotted = None
        self.cp = None
        self.mayavi_atom = None
        self.mayavi_bonds = None
        self.mayavi_unitcell = None

    def clear_plot(self):
        self.scene.mlab.clf(figure=self.scene.mayavi_scene)

    @on_trait_change('scene.activated,n_x,n_y,n_z,edge_atoms')
    def update_plot(self, *args, **kwargs):
        if 'keep_view' in kwargs.keys():
            keep_view = kwargs['keep_view']
        else:
            keep_view = False
        if self.crystal_structure is None:
            return
        self.scene.anti_aliasing_frames = 20
        # We can do normal mlab calls on the embedded scene.
        self.scene.mlab.clf(figure=self.scene.mayavi_scene)
        repeat = [self.n_x, self.n_y, self.n_z]

        if keep_view:
            cur_view = self.scene.mlab.view(figure=self.scene.mayavi_scene)
            cur_roll = self.scene.mlab.roll(figure=self.scene.mayavi_scene)

        if self.show_atoms:
            self.plot_atoms(repeat=repeat,edge_atoms=self.edge_atoms)

        if self.show_bonds:
            self.plot_bonds(repeat=repeat)

        if self.show_unitcell:
            self.plot_unit_cell(repeat=repeat)

        if keep_view:
            self.scene.mlab.view(azimuth=cur_view[0], elevation=cur_view[1], distance=cur_view[2],
                                 focalpoint=cur_view[3], figure=self.scene.mayavi_scene)
            self.scene.mlab.roll(cur_roll, figure=self.scene.mayavi_scene)

    @on_trait_change('show_atoms')
    def _show_atoms_event(self):
        if self.show_atoms:
            repeat = [self.n_x, self.n_y, self.n_z]
            self.plot_atoms(repeat=repeat,edge_atoms=self.edge_atoms)
        else:
            self.mayavi_atom.remove()
            self.mayavi_atom = None

    @on_trait_change('show_bonds')
    def _show_bonds_event(self):
        if self.show_bonds:
            repeat = [self.n_x, self.n_y, self.n_z]
            self.plot_bonds(repeat)
        else:
            self.remove_bonds()

    @on_trait_change('show_unitcell')
    def _show_unitcell_event(self):
        if self.show_unitcell:
            repeat = [self.n_x, self.n_y, self.n_z]
            self.plot_unit_cell(repeat)
        else:
            if self.mayavi_unitcell:
                self.mayavi_unitcell.remove()
                self.mayavi_unitcell = None

    def check_if_line_exists(self, p1, p2, list_of_lines):
        for line in list_of_lines:
            x1 = line[0]
            x2 = line[1]
            norm1 = np.linalg.norm(p1 - x1) + np.linalg.norm(p2 - x2)
            norm2 = np.linalg.norm(p1 - x2) + np.linalg.norm(p2 - x1)
            if norm1 < 0.05 or norm2 < 0.05:
                return True
        return False

    def plot_unit_cell(self, repeat=(1, 1, 1)):
        if type(self.crystal_structure) is sst.MolecularStructure:
            return

        cell = self.crystal_structure.lattice_vectors
        a1 = cell[0, :]
        a2 = cell[1, :]
        a3 = cell[2, :]

        possible_sums = find_possible_sums(repeat)
        connections = find_grid_connections(possible_sums)

        n_grid = len(possible_sums)
        grid_points = np.zeros((n_grid, 3))

        for i, vec_sum in enumerate(possible_sums):
            grid_points[i, :] = vec_sum[0] * a1 + vec_sum[1] * a2 + vec_sum[2] * a3

        mayavi_grid = self.scene.mlab.points3d(grid_points[:, 0], grid_points[:, 1], grid_points[:, 2],
                                               np.ones_like(grid_points[:, 0]), figure=self.scene.mayavi_scene,
                                               scale_factor=0.0, resolution=5)
        mayavi_grid.mlab_source.dataset.lines = np.array(list(connections))

        tube = mlab.pipeline.tube(mayavi_grid, tube_radius=0.05)
        tube.filter.radius_factor = 1.
        # tube.filter.vary_radius = 'vary_radius_by_scalar'
        self.mayavi_unitcell = self.scene.mlab.pipeline.surface(tube, color=(0.8, 0.8, 0.8))
        mayavi_grid.mlab_source.update()

    def plot_atoms(self, repeat=[1, 1, 1],edge_atoms=False):
        abs_coord_atoms = self.crystal_structure.calc_absolute_coordinates(repeat=repeat,edges=edge_atoms)
        # species = set(abs_coord_atoms[:,3].astype(np.int))

        atom_size = 0.4 * np.log(abs_coord_atoms[:, 3]) + 0.6

        pts = self.scene.mlab.points3d(abs_coord_atoms[:, 0], abs_coord_atoms[:, 1], abs_coord_atoms[:, 2],
                                       scale_factor=0.6, vmin=1, vmax=256, figure=self.scene.mayavi_scene,
                                       resolution=20)

        pts.module_manager.scalar_lut_manager.lut.table = lut

        pts.glyph.scale_mode = 'scale_by_vector'
        pts.mlab_source.dataset.point_data.vectors = np.tile(atom_size, (3, 1)).T

        pts.mlab_source.dataset.point_data.scalars = abs_coord_atoms[:, 3]
        self.mayavi_atom = pts

        # n_species = len(species)
        # n_atoms = abs_coord_atoms.shape[0]
        #
        # for specie in species:
        #     species_mask = abs_coord_atoms[:,3].astype(np.int) == specie
        #     sub_coords = abs_coord_atoms[species_mask,:]
        #
        #     cov_radius = cov_radii[specie]
        #     atom_size = 0.4*np.log(specie)+0.6
        #     try:
        #         atomic_color = colors[specie]
        #     except KeyError:
        #         atomic_color = (0.8,0.8,0.8)
        #     mayavi_atom = self.scene.mlab.points3d(sub_coords[:,0],sub_coords[:,1],sub_coords[:,2],
        #                              scale_factor=atom_size,resolution=50,
        #                              color=atomic_color,figure=self.scene.mayavi_scene)
        #     self.mayavi_atoms.append(mayavi_atom)

    def clear_density_plot(self):
        if self.cp is not None:
            pass

    def plot_density(self, ks_density, contours=10, transparent=True, colormap='hot', opacity=0.5):
        repeat = [self.n_x, self.n_y, self.n_z]

        cur_view = self.scene.mlab.view()
        cur_roll = self.scene.mlab.roll()

        if cur_view is None:
            cur_view = self.scene.mlab.view()
        if cur_roll is None:
            cur_roll = self.scene.mlab.roll()

        dens = ks_density.density
        dens_plot = np.tile(dens, repeat)

        if type(contours) == int:
            color = None
        elif len(contours) > 1:
            color = None
        else:
            color = (1.0, 1.0, 0.2)
        self.cp = self.scene.mlab.contour3d(dens_plot, contours=contours, transparent=transparent,
                                            opacity=opacity, colormap=colormap, color=color,
                                            figure=self.scene.mayavi_scene)
        # Do some tvtk magic in order to allow for non-orthogonal unit cells:

        polydata = self.cp.actor.actors[0].mapper.input
        pts = np.array(polydata.points) - 1

        if type(ks_density) is sst.KohnShamDensity:
            unit_cell = self.crystal_structure.lattice_vectors
            # Transform the points to the unit cell:
            larger_cell = np.zeros((3, 3))
            larger_cell[0, :] = unit_cell[0, :] * repeat[0]
            larger_cell[1, :] = unit_cell[1, :] * repeat[1]
            larger_cell[2, :] = unit_cell[2, :] * repeat[2]

            polydata.points = np.dot(pts, larger_cell / np.array(dens_plot.shape)[:, np.newaxis])
        elif type(ks_density) is sst.MolecularDensity:
            lattice_vecs = ks_density.grid_vectors
            origin = ks_density.origin
            polydata.points = np.dot(pts, lattice_vecs / np.array(dens_plot.shape)[:, np.newaxis]) + origin

        else:
            raise ValueError('Invalid type for density')

        # self.scene.mlab.view(distance='auto')
        self.scene.mlab.view(azimuth=cur_view[0], elevation=cur_view[1], distance=cur_view[2], focalpoint=cur_view[3],
                             figure=self.scene.mayavi_scene)
        self.scene.mlab.roll(cur_roll, figure=self.scene.mayavi_scene)
        self.density_plotted = ks_density

    def plot_bonds(self, repeat=[1, 1, 1]):
        if self.mayavi_atom is None:
            self.plot_atoms(repeat=repeat,edge_atoms=self.edge_atoms)

        abs_coord_atoms = self.crystal_structure.calc_absolute_coordinates(repeat=repeat,edges=self.edge_atoms)
        bonds = self.crystal_structure.find_bonds(abs_coord_atoms)

        self.mayavi_atom.mlab_source.dataset.lines = np.array(bonds)

        tube = mlab.pipeline.tube(self.mayavi_atom, tube_radius=0.15)
        tube.filter.radius_factor = 1.
        # tube.filter.vary_radius = 'vary_radius_by_scalar'
        self.mayavi_bonds = self.scene.mlab.pipeline.surface(tube, color=(0.8, 0.8, 0.8))
        self.mayavi_atom.mlab_source.update()

        if not self.show_atoms:
            self.mayavi_atom.remove()
            self.mayavi_atom = None

            # paths = sst.bonds_to_path(bonds)
            #
            # for path in paths:
            #     x = abs_coord_atoms[path,0]
            #     y = abs_coord_atoms[path,1]
            #     z = abs_coord_atoms[path,2]
            #
            #     mayavi_bond = self.scene.mlab.plot3d(x,y,z, tube_radius=0.125,tube_sides=18,figure=self.scene.mayavi_scene)
            #     self.mayavi_bonds.append(mayavi_bond)

    def remove_bonds(self):
        if self.mayavi_bonds:
            self.mayavi_bonds.remove()


class VolumeSlicer(HasTraits):
    data = Array()
    crystal_structure = None

    scene3d = Instance(MlabSceneModel, ())

    # The data source
    data_src3d = Instance(Source)

    ipw_3d = Instance(PipelineBase)

    _axis_names = dict(x=0, y=1, z=2)

    show_unitcell = Bool(True)
    show_bonds = Bool(True)
    show_atoms = Bool(True)

    # ---------------------------------------------------------------------------
    def __init__(self, **traits):
        super(VolumeSlicer, self).__init__(**traits)
        # Force the creation of the image_plane_widgets:
        self.ipw_3d
        self.cut_plane = None
        self.mayavi_atoms = []
        self.mayavi_unitcell = None
        self.mayavi_bonds = []

    # def _ipw_3d_default(self):
    #     return self.make_ipw_3d()

    def _data_src3d_default(self):
        return mlab.pipeline.scalar_field(self.data,
                                          figure=self.scene3d.mayavi_scene)

    def rescale_polydata_points(self, polydata):
        pts = np.array(polydata.points) - 1
        unit_cell = self.crystal_structure.lattice_vectors
        polydata.points = np.dot(pts, unit_cell / np.array(self.data.shape)[:, np.newaxis])

    def make_ipw_3d(self):
        cut_plane = mlab.pipeline.scalar_cut_plane(self.data_src3d,
                                                   figure=self.scene3d.mayavi_scene, colormap=self.colormap)
        # cut_plane.implicit_plane.origin = (1, 1, 1)
        cut_plane.implicit_plane.widget.enabled = False

        polydata = cut_plane.actor.actors[0].mapper.input
        self.rescale_polydata_points(polydata)

        return cut_plane

    def set_data(self, data, crystal_structure):
        self.data = data
        self.data_src3d.scalar_data = data
        self.data_src3d.update()
        self.crystal_structure = crystal_structure

    def display_scene3d(self, colormap=None):
        self.clear_scene()

        if len(self.data) == 0 or self.crystal_structure is None:
            return

        if colormap is not None:
            self.colormap = colormap

        if self.show_atoms:
            self.plot_atoms()
        if self.show_bonds:
            self.plot_bonds()
        if self.show_unitcell:
            self.plot_unitcell()

        self.cut_plane = self.make_ipw_3d()
        self.scene3d.mlab.view(40, 50)

        self.scene3d.scene.interactor.interactor_style = tvtk.InteractorStyleTerrain()

    @on_trait_change('show_atoms')
    def _show_atoms_event(self):
        if self.show_atoms:
            self.plot_atoms()
        else:
            for mayavi_atom in self.mayavi_atoms:
                mayavi_atom.remove()
            self.mayavi_atoms = []

    @on_trait_change('show_bonds')
    def _show_bonds_event(self):
        if self.show_bonds:
            self.plot_bonds()
        else:
            for mayavi_bond in self.mayavi_bonds:
                mayavi_bond.remove()
            self.mayavi_bonds = []

    @on_trait_change('show_unitcell')
    def _show_unitcell_event(self):
        if self.show_unitcell:
            self.plot_unitcell()
        else:
            if self.mayavi_unitcell is not None:
                self.mayavi_unitcell.remove()
                self.mayavi_unitcell = None

    def clear_scene(self):
        self.scene3d.mlab.clf(figure=self.scene3d.mayavi_scene)
        if self.cut_plane is not None:
            self.cut_plane.remove()
            self.cut_plane = None
        if self.mayavi_unitcell is not None:
            self.mayavi_unitcell.remove()
            self.mayavi_unitcell = None

    def plot_unitcell(self):
        self.mayavi_unitcell = mlab.pipeline.outline(self.data_src3d,
                                                     figure=self.scene3d.mayavi_scene, reset_zoom=False,
                                                     colormap=self.colormap
                                                     )

        polydata = self.mayavi_unitcell.actor.actors[0].mapper.input
        self.rescale_polydata_points(polydata)

    def plot_atoms(self, repeat=[1, 1, 1]):
        self.mayavi_atoms = []
        abs_coord_atoms = self.crystal_structure.calc_absolute_coordinates(repeat=repeat)
        n_atoms = abs_coord_atoms.shape[0]

        species = set(abs_coord_atoms[:, 3].astype(np.int))
        n_species = len(species)

        for specie in species:
            species_mask = abs_coord_atoms[:, 3].astype(np.int) == specie
            sub_coords = abs_coord_atoms[species_mask, :]

            cov_radius = cov_radii[specie]
            atom_size = 0.4 * np.log(specie) + 0.6
            try:
                atomic_color = colors[specie]
            except KeyError:
                atomic_color = (0.8, 0.8, 0.8)
            mayavi_atom = self.scene3d.mlab.points3d(sub_coords[:, 0], sub_coords[:, 1], sub_coords[:, 2],
                                                     scale_factor=atom_size, resolution=150,
                                                     color=atomic_color, figure=self.scene3d.mayavi_scene)
            self.mayavi_atoms.append(mayavi_atom)

    def plot_bonds(self, repeat=[1, 1, 1]):
        self.mayavi_bonds = []
        abs_coord_atoms = self.crystal_structure.calc_absolute_coordinates(repeat=repeat)
        bonds = self.crystal_structure.find_bonds(abs_coord_atoms)
        n_atoms = abs_coord_atoms.shape[0]

        paths = sst.bonds_to_path(bonds)

        for path in paths:
            x = abs_coord_atoms[path, 0]
            y = abs_coord_atoms[path, 1]
            z = abs_coord_atoms[path, 2]

            mayavi_bond = self.scene3d.mlab.plot3d(x, y, z, tube_radius=0.125, tube_sides=18,
                                                   figure=self.scene3d.mayavi_scene)
            self.mayavi_bonds.append(mayavi_bond)


    view = View(
        Group(
            Item('scene3d',
                 editor=SceneEditor(scene_class=MayaviScene),
                 height=250, width=300),
            Group('_', 'show_unitcell', 'show_bonds', 'show_atoms', orientation='horizontal'),
            show_labels=False,
        ),
        resizable=True,
        title='Volume Slicer',
    )


class OpticalSpectrumVisualization(QtGui.QWidget):
    def __init__(self, parent=None):
        super(OpticalSpectrumVisualization, self).__init__()
        self.first_plot_bool = True
        self.last_optical_spectrum = None
        # a figure instance to plot on
        self.figure = plt.figure(1)
        plt.close(plt.figure(1))
        self.ax = None

        self.canvas = FigureCanvas(self.figure)

        self.toolbar = NavigationToolbar(self.canvas, self)

        color = self.palette().color(QtGui.QPalette.Base)
        self.figure.patch.set_facecolor([color.red() / 255, color.green() / 255, color.blue() / 255])

        if sum([color.red(), color.blue(), color.green()]) / 3 < 100:
            self.dark_mode = True
            self.bg_color = [color.red() / 255, color.green() / 255, color.blue() / 255]
        else:
            self.dark_mode = False
            self.bg_color = None

        # self.figure.patch.set_facecolor([236 / 255, 236 / 255, 236 / 255])
        # self.figure.patch.set_alpha(1.0)
        # self.figure.patch.set_facecolor('blue')

        layout = QtGui.QVBoxLayout()
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)

        option_widget = QtGui.QWidget()
        option_widget.setFixedHeight(60)
        option_layout = QtGui.QHBoxLayout(option_widget)
        option_layout.setAlignment(QtCore.Qt.AlignLeft)
        layout.addWidget(option_widget)
        from main import EntryWithLabel

        self.select_epsilon_cb = QtGui.QComboBox(self)
        option_layout.addWidget(self.select_epsilon_cb)
        self.select_epsilon_cb.addItem('Angular mean')
        self.select_epsilon_cb.addItem(u"ε_11")
        self.select_epsilon_cb.addItem(u"ε_22")
        self.select_epsilon_cb.addItem(u"ε_33")

        self.select_epsilon_cb.setCurrentIndex(0)
        self.select_epsilon_cb.currentIndexChanged.connect(
            lambda: self.plot(self.last_optical_spectrum[0], name_list=self.last_optical_spectrum[1]))

        self.imaginary_checkbox = QtGui.QCheckBox('Imag', self)
        option_layout.addWidget(self.imaginary_checkbox)
        self.imaginary_checkbox.toggle()
        self.imaginary_checkbox.stateChanged.connect(
            lambda: self.plot(self.last_optical_spectrum[0], name_list=self.last_optical_spectrum[1]))

        self.real_checkbox = QtGui.QCheckBox('Real', self)
        option_layout.addWidget(self.real_checkbox)
        self.real_checkbox.stateChanged.connect(
            lambda: self.plot(self.last_optical_spectrum[0], name_list=self.last_optical_spectrum[1]))

        width_text = 70
        width_label = 40

        self.Emin_entry = EntryWithLabel(option_widget, 'Emin', width_text=width_text, width_label=width_label)
        self.Emin_entry.connect_editFinished(
            lambda: self.plot(self.last_optical_spectrum[0], name_list=self.last_optical_spectrum[1]))
        option_layout.addWidget(self.Emin_entry)

        self.Emax_entry = EntryWithLabel(option_widget, 'Emax', width_text=width_text, width_label=width_label)
        self.Emax_entry.connect_editFinished(
            lambda: self.plot(self.last_optical_spectrum[0], name_list=self.last_optical_spectrum[1]))
        option_layout.addWidget(self.Emax_entry)

        self.eps_min_entry = EntryWithLabel(option_widget, u"ε min", width_text=width_text, width_label=width_label)
        self.eps_min_entry.connect_editFinished(
            lambda: self.plot(self.last_optical_spectrum[0], name_list=self.last_optical_spectrum[1]))
        option_layout.addWidget(self.eps_min_entry)

        self.eps_max_entry = EntryWithLabel(option_widget, u"ε max", width_text=width_text, width_label=width_label)
        self.eps_max_entry.connect_editFinished(
            lambda: self.plot(self.last_optical_spectrum[0], name_list=self.last_optical_spectrum[1]))
        option_layout.addWidget(self.eps_max_entry)

        self.broadening_entry = EntryWithLabel(option_widget, u"Γ", width_text=width_text, width_label=width_label)
        self.broadening_entry.connect_editFinished(
            lambda: self.plot(self.last_optical_spectrum[0], name_list=self.last_optical_spectrum[1]))
        option_layout.addWidget(self.broadening_entry)

        self.broadening_mode_cb = QtGui.QComboBox(self)
        option_layout.addWidget(self.broadening_mode_cb)
        self.broadening_mode_cb.addItem('Lorentzian')
        self.broadening_mode_cb.addItem('Gaussian')
        self.broadening_mode_cb.addItem('None')

        self.broadening_mode_cb.setCurrentIndex(0)
        self.broadening_mode_cb.currentIndexChanged.connect(
            lambda: self.plot(self.last_optical_spectrum[0], name_list=self.last_optical_spectrum[1]))
        self.broadening_mode_cb.setMaximumWidth(150)

        self.setLayout(layout)
        self.show()

    def clear_plot(self):
        if not self.first_plot_bool:
            self.figure.clf()
            self.first_plot_bool = True
            self.canvas.draw()

    def read_entries(self):
        try:
            Emin = float(self.Emin_entry.get_text())
        except Exception:
            Emin = None

        try:
            Emax = float(self.Emax_entry.get_text())
        except Exception:
            Emax = None

        try:
            eps_max = float(self.eps_max_entry.get_text())
        except Exception:
            eps_max = None

        try:
            eps_min = float(self.eps_min_entry.get_text())
        except Exception:
            eps_min = None

        try:
            gamma = float(self.broadening_entry.get_text())
        except Exception:
            gamma = None

        broaden_mode = self.broadening_mode_cb.currentText().lower()

        return {'Emin': Emin, 'Emax': Emax, 'eps min': eps_min, 'eps max': eps_max, 'Gamma': gamma,
                'broaden mode': broaden_mode}

    def plot(self, optical_spectrum_list, *args, **kwargs):
        name_list = kwargs.pop('name_list', None)
        if optical_spectrum_list is None:
            return

        if type(optical_spectrum_list) is not list:
            optical_spectrum_list = [optical_spectrum_list]

        self.last_optical_spectrum = [optical_spectrum_list, name_list]
        if self.first_plot_bool:
            self.ax = self.figure.add_subplot(111)
            self.ax.format_coord = lambda x, y: u'E = {0:1.2f} eV, ε = {1:1.3f}'.format(x, y)
        self.ax.cla()

        if self.dark_mode:
            set_dark_mode_matplotlib(self.figure, self.ax, self.bg_color)

        entry_values = self.read_entries()
        Emin = entry_values['Emin']
        Emax = entry_values['Emax']
        eps_min = entry_values['eps min']
        eps_max = entry_values['eps max']
        gamma = entry_values['Gamma']
        broaden_mode = entry_values['broaden mode']

        if Emin is not None:
            self.ax.set_xlim(left=Emin)
        else:
            self.ax.set_xlim(left=optical_spectrum_list[0].energy.min())

        if Emax is not None:
            self.ax.set_xlim(right=Emax)
        else:
            self.ax.set_xlim(right=optical_spectrum_list[0].energy.max())

        if eps_min is not None:
            self.ax.set_ylim(bottom=eps_min)
        if eps_max is not None:
            self.ax.set_ylim(top=eps_max)

        if name_list is None:
            name_list = len(optical_spectrum_list) * ['']

        handles = []
        for optical_spectrum, name in zip(optical_spectrum_list, name_list):

            if self.imaginary_checkbox.checkState():
                E_plot = optical_spectrum.energy

                cur_index_eps = self.select_epsilon_cb.currentIndex()
                if cur_index_eps == 0:
                    epsilon = optical_spectrum.epsilon2
                elif cur_index_eps == 1:
                    epsilon = optical_spectrum.epsilon2_11
                elif cur_index_eps == 2:
                    epsilon = optical_spectrum.epsilon2_22
                elif cur_index_eps == 3:
                    epsilon = optical_spectrum.epsilon2_33

                if gamma is None:
                    epsilon_plot = epsilon
                else:
                    E_plot, epsilon_plot = self.broaden_spectrum(E_plot, epsilon, width=gamma, mode=broaden_mode)

                p, = self.ax.plot(E_plot, epsilon_plot, linewidth=2, label=name + '_imag')
                handles.append(p)

            if self.real_checkbox.checkState():
                E_plot = optical_spectrum.energy

                cur_index_eps = self.select_epsilon_cb.currentIndex()
                if cur_index_eps == 0:
                    epsilon = optical_spectrum.epsilon1
                elif cur_index_eps == 1:
                    epsilon = optical_spectrum.epsilon1_11
                elif cur_index_eps == 2:
                    epsilon = optical_spectrum.epsilon1_22
                elif cur_index_eps == 3:
                    epsilon = optical_spectrum.epsilon1_33

                if gamma is None:
                    epsilon_plot = epsilon
                else:
                    E_plot, epsilon_plot = self.broaden_spectrum(E_plot, epsilon, width=gamma, mode=broaden_mode)

                p, = self.ax.plot(E_plot, epsilon_plot, linewidth=2, label=name + '_real')
                handles.append(p)

        self.ax.set_xlabel('Energy [eV]')
        self.ax.set_ylabel(r'Dielectric function $\varepsilon(\omega)$')

        if name_list is not None and len(handles) > 1:
            legend = self.ax.legend(loc='best', fancybox=True, framealpha=0.9)
            legend_frame = legend.get_frame()
            legend_frame.set_facecolor([0.95, 0.95, 0.95])
            legend_frame.set_linewidth(0)

        if self.first_plot_bool:
            self.first_plot_bool = False
            self.figure.tight_layout()
        self.canvas.draw()

    def broaden_spectrum(self, energy, epsilon, width, mode='lorentzian'):

        if width is None or mode == 'none':
            return energy, epsilon
        if mode == 'lorentzian':
            def broaden_function(x, width):
                return width ** 2 / (x ** 2 + width ** 2)
        elif mode == 'gaussian':
            def broaden_function(x, width):
                return np.exp(-x ** 2 / (2 * width ** 2))

        E_range = energy.max() - energy.min()
        dx = E_range / len(energy)

        conv_range = 50 * width
        gx = np.arange(-conv_range / 2, conv_range / 2, dx)
        broadenarray = broaden_function(gx, width)
        broadenarray = broadenarray / np.sum(broadenarray)

        epsilon_out = np.convolve(epsilon, broadenarray, mode="full")
        energy_out = np.linspace(energy.min() - conv_range / 2, energy.max() + conv_range / 2, len(epsilon_out))

        return energy_out, epsilon_out

    def export(self, filename, spectrum, code=False):
        entry_values = self.read_entries()
        gamma = entry_values['Gamma']
        broaden_mode = entry_values['broaden mode']

        energy = spectrum.energy

        E_plot, epsilon_plot = self.broaden_spectrum(energy, spectrum.epsilon2, width=gamma, mode=broaden_mode)

        all_epsilons = spectrum.all_epsilons
        spectrum_exists = [True if x is not None else False for x in all_epsilons]
        nr_of_spectra = sum(spectrum_exists)
        data = np.zeros((len(E_plot), nr_of_spectra + 1))
        data[:, 0] = E_plot

        i = 1
        for epsilon in all_epsilons:
            if epsilon is None:
                continue
            else:
                E_plot, epsilon_plot = self.broaden_spectrum(energy, epsilon, width=gamma, mode=broaden_mode)
                data[:, i] = epsilon_plot
                i += 1

        full_header = ['epsilon1', 'epsilon1_11', 'epsilon1_22', 'epsilon1_33', 'epsilon2', 'epsilon2_11',
                       'epsilon2_22', 'epsilon2_33']

        from itertools import compress
        header = 'Energy [eV] ' + ' '.join(list(compress(full_header, spectrum_exists)))

        np.savetxt(filename, data, header=header)


class BandStructureVisualization(QtGui.QWidget):
    def __init__(self, parent=None):
        super(BandStructureVisualization, self).__init__()
        self.first_plot_bool = True
        # a figure instance to plot on
        self.figure = plt.figure(1)
        plt.close(plt.figure(1))
        self.ax = None
        self.last_bandstructure = None
        self.canvas = FigureCanvas(self.figure)

        self.toolbar = NavigationToolbar(self.canvas, self)

        color = self.palette().color(QtGui.QPalette.Base)
        self.figure.patch.set_facecolor([color.red() / 255, color.green() / 255, color.blue() / 255])

        if sum([color.red(), color.blue(), color.green()]) / 3 < 100:
            self.dark_mode = True
            self.bg_color = [color.red() / 255, color.green() / 255, color.blue() / 255]
        else:
            self.dark_mode = False
            self.bg_color = None

        layout = QtGui.QVBoxLayout()
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)

        option_widget = QtGui.QWidget()
        option_widget.setFixedHeight(60)
        option_layout = QtGui.QHBoxLayout(option_widget)
        option_layout.setAlignment(QtCore.Qt.AlignLeft)
        layout.addWidget(option_widget)
        from main import EntryWithLabel
        self.Emin_entry = EntryWithLabel(option_widget, 'Emin')
        self.Emin_entry.connect_editFinished(lambda: self.plot(self.last_bandstructure))
        option_layout.addWidget(self.Emin_entry)
        self.Emax_entry = EntryWithLabel(option_widget, 'Emax')
        self.Emax_entry.connect_editFinished(lambda: self.plot(self.last_bandstructure))
        option_layout.addWidget(self.Emax_entry)

        # layout.addWidget(self.button)
        self.setLayout(layout)
        self.show()

    def clear_plot(self):
        if not self.first_plot_bool:
            self.ax.cla()
            self.canvas.draw()

    def plot(self, band_structure, *args, **kwargs):
        if band_structure is None:
            return
        if type(band_structure) is list:
            band_structure = band_structure[0]
        self.last_bandstructure = band_structure
        if self.first_plot_bool:
            self.ax = self.figure.add_subplot(111)

        if type(band_structure) is sst.EnergyDiagram:
            self.plot_energy_diagram(band_structure)
        elif type(band_structure) is sst.BandStructure:
            self.plot_bandstructure(band_structure)

        if self.first_plot_bool:
            self.first_plot_bool = False
            self.figure.tight_layout()

        try:
            Emin = float(self.Emin_entry.get_text())
        except Exception:
            Emin = None

        try:
            Emax = float(self.Emax_entry.get_text())
        except Exception:
            Emax = None

        if Emin is not None:
            self.ax.set_ylim(bottom=Emin)
        if Emax is not None:
            self.ax.set_ylim(top=Emax)
        self.canvas.draw()

    def plot_energy_diagram(self, energy_diagram):

        self.ax.format_coord = lambda x, y: 'E = {0:1.2f} eV, E_gap = {1:1.2f} eV'.format(
            *self.make_interactive_text_energy_diagram(x, y, energy_diagram))

        self.ax.cla()
        if self.dark_mode:
            set_dark_mode_matplotlib(self.figure, self.ax, self.bg_color)

        energies = energy_diagram.energies
        labels = energy_diagram.labels
        gap = energy_diagram.homo_lumo_gap

        self.ax.set_xlim(-1, 1)
        self.ax.get_xaxis().set_ticks([])
        energy_range = max(energies) - min(energies)
        self.ax.set_ylim(min(energies) - energy_range * 0.05, max(energies) + energy_range * 0.05)

        self.ax.plot([-0.4, 0.4], [energy_diagram.E_fermi] * 2, 'k--')
        for i, info in enumerate(zip(energies, labels)):
            energy = info[0]
            label = info[1]
            self.ax.plot([-0.3, 0.3], [energy, energy])
            if i % 2 == 0:
                text_pos = 0.38
                alignment = 'right'
            else:
                text_pos = -0.38
                alignment = 'left'
            self.ax.text(text_pos, energy, label, verticalalignment='center', horizontalalignment=alignment)

        self.ax.set_ylabel("Energy eV")
        if self.dark_mode:
            title_color = 'white'
        else:
            title_color = 'black'
        self.ax.set_title('Energy diagram. Homo-Lumo gap = {0:1.2f} eV'.format(gap), color=title_color)

    def plot_bandstructure(self, band_structure):
        self.ax.format_coord = lambda x, y: 'k_d = {0:1.1f}, E = {1:1.2f} eV, Gap = {2:1.2f} eV'.format(
            *self.make_interactive_text(x, y, band_structure))

        self.ax.cla()
        if self.dark_mode:
            set_dark_mode_matplotlib(self.figure, self.ax, self.bg_color)
            title_color = 'white'
        else:
            title_color = 'black'

        for band in band_structure.bands:
            self.ax.plot(band[:, 0], band[:, 1], color='#1f77b4', linewidth=2)

        xlength = band[:, 0].max() - band[:, 0].min()
        self.ax.set_xlim(band[:, 0].min() - xlength / 800, band[:, 0].max())
        self.ax.plot([band[:, 0].min(), band[:, 0].max()], [0, 0], 'k--')

        if band_structure.special_k_points is not None and len(band_structure.special_k_points) > 0:
            for xc, xl in band_structure.special_k_points:
                self.ax.axvline(x=xc, color='k', linewidth=1.5)

            unzipped_k = list(zip(*band_structure.special_k_points))
            special_k_points = unzipped_k[0]
            special_k_points_label = convert_to_greek(unzipped_k[1])

            self.ax.set_xticks(special_k_points)
            self.ax.set_xticklabels(special_k_points_label, rotation='horizontal', horizontalalignment='center')
        else:
            special_k_points = []
            special_k_points_label = []

        if band_structure.bs_type == 'electronic':
            bandgap = band_structure.bandgap
            k_bandgap = band_structure.k_bandgap
            if k_bandgap is None and bandgap is not None:
                title_bandgap = ' $E_g = %1.1f $ eV' % bandgap + ' (indirect)'
            elif bandgap is None:
                title_bandgap = ' (metallic)'
            elif k_bandgap in special_k_points:
                k_bandgap_label = np.array(special_k_points_label)[k_bandgap == special_k_points][0]
                title_bandgap = ' $E_g$ = %1.1f eV' % bandgap + ' at ' + k_bandgap_label
            else:
                title_bandgap = ' $E_g$ = %1.1f eV' % bandgap + ' (direct)'

            self.ax.set_title('KS bandstructure,' + title_bandgap, fontsize=25, color=title_color)
        else:
            self.ax.set_ylim(bottom=0)
            self.ax.set_title('Phonon bandstructure', fontsize=25, color=title_color)

    def make_interactive_text(self, k_in, E, band_structure):
        bands = band_structure.bands
        k = band_structure.bands[0][:, 0]
        index_k = np.argmin(np.abs(k - k_in))
        E_values_at_k = sorted([band[index_k, 1] for band in bands])
        E_index = bisect(E_values_at_k, E)
        try:
            E_above = E_values_at_k[E_index]
            E_below = E_values_at_k[E_index - 1]
            gap = E_above - E_below
        except IndexError:
            gap = 0
        if gap < 0:
            gap = 0
        return [k_in, E, gap]

    def make_interactive_text_energy_diagram(self, x, E, energy_diagram):
        energies = energy_diagram.energies
        E_index = bisect(energies, E)
        try:
            E_above = energies[E_index]
            E_below = energies[E_index - 1]
            gap = E_above - E_below
        except IndexError:
            gap = 0
        if gap < 0:
            gap = 0
        return [E, gap]

    def export(self, filename, band_structure, code=False):
        bands = band_structure.bands
        data = np.zeros((bands[0].shape[0], len(bands) + 1))
        data[:, 0] = bands[0][:, 0]
        for i, band in enumerate(bands):
            data[:, i + 1] = band[:, 1]
        np.savetxt(filename, data)

        if code:
            # TODO do this!
            code_string = """
import numpy as np
import matplotlib.pyplot as plt

            
self.first_plot_bool = True
self.figure = plt.figure(1)
plt.close(plt.figure(1))
self.ax = None
self.last_bandstructure = None"""


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
        self.figure.patch.set_facecolor([color.red() / 255, color.green() / 255, color.blue() / 255])

        if sum([color.red(), color.blue(), color.green()]) / 3 < 100:
            self.dark_mode = True
            self.bg_color = [color.red() / 255, color.green() / 255, color.blue() / 255]
        else:
            self.dark_mode = False
            self.bg_color = None

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

    def plot(self, scf_data):
        if self.first_plot_bool:
            self.ax = self.figure.add_subplot(111)
        self.ax.cla()

        if self.dark_mode:
            set_dark_mode_matplotlib(self.figure, self.ax, self.bg_color)

        xmax = int(round(scf_data[:, 0].max())) + 1
        xmin = int(scf_data[:, 0].min())
        dist = xmax - xmin
        if dist < 10:
            step = 1
            ms = 12
        elif dist < 20:
            step = 2
            ms = 12
        elif dist < 50:
            step = 5
            ms = 9
        else:
            step = int(round((dist // 10) / 5) * 5)
            ms = 6
        xint = list(range(0, xmax, step))[1:]

        self.ax.plot(scf_data[:, 0], scf_data[:, 1], marker='o', color='#1f77b4', ms=ms, markeredgecolor='none',
                     linewidth=0)
        self.ax.plot(scf_data[:, 0], scf_data[:, 1], linestyle='--', color='#1f77b4', linewidth=3, alpha=0.5)
        self.ax.set_xlim(scf_data[:, 0].min() - 0.1, scf_data[:, 0].max() + 0.1)

        self.ax.set_xticks(xint)
        self.ax.set_xlabel('Scf iteration')
        self.ax.set_ylabel('Total Energy')
        if self.first_plot_bool:
            self.first_plot_bool = False
            self.figure.tight_layout()
        self.canvas.draw()


def set_dark_mode_matplotlib(f, ax, color):
    ax.tick_params(axis='x', colors='white')
    ax.tick_params(axis='y', colors='white')
    ax.spines['bottom'].set_color('white')
    ax.spines['top'].set_color('white')
    ax.spines['left'].set_color('white')
    ax.spines['right'].set_color('white')
    ax.xaxis.label.set_color('white')
    ax.yaxis.label.set_color('white')
    ax.set_facecolor(color)


if __name__ == "__main__":
    atoms = np.array([[0, 0, 0, 6], [0.25, 0.25, 0.25, 6]])
    unit_cell = 6.719 * np.array([[0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]])

    from solid_state_tools import CrystalStructure

    crystal_structure = CrystalStructure(unit_cell, atoms)
    vis = StructureVisualization(crystal_structure)
    vis.configure_traits()

    # x, y, z = np.ogrid[-5:5:64j, -5:5:64j, -5:5:64j]
    # data = np.sin(3 * x) / x + 0.05 * z ** 2 + np.cos(3 * y)
    #
    # atoms = np.array([[0, 0, 0, 6], [0.25, .25, 0.25, 6]])
    # unit_cell = 6.719 * np.array([[0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0.0]])
    # from solid_state_tools import CrystalStructure
    # crystal_structure = CrystalStructure(unit_cell, atoms)
    #
    # m = VolumeSlicer(data=data,crystal_structure=crystal_structure)
    # m.configure_traits()
