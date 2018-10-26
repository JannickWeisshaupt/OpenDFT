
import numpy as np
from mayavi import mlab
fig = mlab.figure(1, bgcolor = (1,1,1), fgcolor = (0.5, 0.5, 0.5), size=(350, 350))
mlab.clf()

# The position of the atoms
atoms_x = np.array([5,5.03148657000,3.03812145000,6.9887797000,5.03144969000   ])
atoms_y = np.array([4.59995459000,5.36509848000,3.33298297000,3.24804443000,6.45391950000])
atoms_z = np.array([5.20005652000,7.3733690000,4.82797246000,4.80229239000,3.73234170000])

res = 200

C = mlab.points3d(atoms_x[0], atoms_y[0], atoms_z[0],
                  scale_factor=3.7,
                  resolution=res,
                  color=(217/255, 78/255, 19/255),
                  scale_mode='none')

H = mlab.points3d(atoms_x[1:], atoms_y[1:], atoms_z[1:],
                   scale_factor=2.7,
                   resolution=res,
                   color=(103/255, 181/255, 44/255),
                   scale_mode='none')

H.glyph.glyph_source.glyph_source.phi_resolution = res
H.glyph.glyph_source.glyph_source.theta_resolution = res
fig.scene.anti_aliasing_frames = 20
mlab.draw()

# mlab.savefig('logo_methan.png')
mlab.show()