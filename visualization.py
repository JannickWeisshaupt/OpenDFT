from traits.api import HasTraits, Instance, on_trait_change, Range, Bool, Button
from traitsui.api import View, Item, Group
from mayavi.core.ui.api import MayaviScene, MlabSceneModel, \
    SceneEditor

import numpy as np

colors = {0: (0.8, 0.8, 0.8), 8: (1, 0, 0)}


class StructureVisualization(HasTraits):
    n_x    = Range(1, 10, 1, mode='spinner')#)
    n_y  = Range(1, 10, 1,  mode='spinner')#mode='spinner')
    n_z  = Range(1, 10, 1,  mode='spinner')#mode='spinner')
    prop_but = Button(label='Properties')

    show_unitcell = Bool(True)

    view = View(Item('scene', editor=SceneEditor(scene_class=MayaviScene),
                     height=450, width=500, show_label=False),Group('_', 'n_x', 'n_y', 'n_z','show_unitcell','prop_but',orientation='horizontal'),
                resizable=True, # We need this to resize with the parent widget
                )
    scene = Instance(MlabSceneModel, ())

    def __init__(self, crystal_structure):
        super().__init__()

        self.crystal_structure = crystal_structure


    @on_trait_change('scene.activated,n_x,n_y,n_z,show_unitcell')
    def update_plot(self):
        if self.crystal_structure is None:
            return
        self.scene.anti_aliasing_frames = 5
        # We can do normal mlab calls on the embedded scene.
        self.scene.mlab.clf()
        repeat = [self.n_x,self.n_y,self.n_z]
        self.plot_atoms(repeat=repeat)
        if self.show_unitcell:
            self.plot_unit_cell(repeat=repeat)

    def plot_unit_cell(self, repeat=[1, 1, 1]):
        cell = self.crystal_structure.lattice_vectors

        for j1 in range(repeat[0]):
            for j2 in range(repeat[1]):
                for j3 in range(repeat[2]):
                    offset = j1 * cell[0, :] + j2 * cell[1, :] + j3 * cell[2, :]
                    for i1, a in enumerate(cell):
                        i2 = (i1 + 1) % 3
                        i3 = (i1 + 2) % 3
                        for b in [np.zeros(3), cell[i2]]:
                            for c in [np.zeros(3), cell[i3]]:
                                p1 = b + c + offset
                                p2 = p1 + a
                                self.scene.mlab.plot3d([p1[0], p2[0]], [p1[1], p2[1]], [p1[2], p2[2]], tube_radius=0.1)
                                # self.scene.mlab.points3d([p1[0], p2[0]],
                                #                          [p1[1], p2[1]],
                                #                          [p1[2], p2[2]], scale_factor=.2)

    def plot_atoms(self, repeat=[1, 1, 1]):
        abs_coord_atoms = self.crystal_structure.calc_absolute_coordinates(repeat=repeat)

        # self.scene.mlab.plot3d(*abs_coord_atoms[:,:3].T, [1, 2, 1],
        #             tube_radius=0.4, colormap='Reds')

        n_atoms = abs_coord_atoms.shape[0]

        for i in range(n_atoms):
            self.scene.mlab.points3d(*abs_coord_atoms[i, :],
                                     scale_factor=3,
                                     resolution=10,
                                     color=colors[abs_coord_atoms[i, 3]],
                                     scale_mode='none')