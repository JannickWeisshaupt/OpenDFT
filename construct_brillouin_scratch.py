from __future__ import division

import numpy as np
import mayavi.mlab as mlab



phi = 90/180*np.pi
l1 = 10*np.array([0.0,0.5,0.5])
l2 = 10*np.array([0.5,0,0.5])
l3 = 10*np.array([0.5,0.5,0])

# l1 = 10*np.array([1,0,0])
# l2 = 10*np.array([0,1,0])
# l3 = 10*np.array([0,0.1,0.9])

origin = 0*l1


point_list = []

for i in range(-1,2):
    for j in range(-1,2):
        for k in range(-1,2):

            point_list.append(l1*i+l2*j+k*l3)

N_points = len(point_list)
wigner_points = []

for i in range(N_points):
    for j in range(N_points):
        for k in range(N_points):

            if i==j or i==k or j==k:
                continue
            x1,y1,z1 = origin
            x2,y2,z2 = point_list[i]
            x3, y3, z3 = point_list[j]
            x4, y4, z4 = point_list[k]

            A = np.array([[2 * x1 - x2 - x3 -x4, 2 * y1 - y2 - y3 - y4, 2 * z1 - z2 - z3 - z4],
                          [2 * x2 - x1 - x3 -x4, 2 * y2 - y1 - y3 - y4, 2 * z2 - z1 - z3 - z4],
                          [2 * x3 - x1 - x2 -x4, 2 * y3 - y1 - y2 - y4, 2 * z3 - z1 - z2 - z4]
                          ])
            r = np.linalg.matrix_rank(A)
            if r<3:
                continue
            B = np.array([[x1 ** 2 - 0.5 * (x2 ** 2 + x3 ** 2 + x4**2) + y1 ** 2 - 0.5 * (y2 ** 2 + y3 ** 2 + y4**2)+ z1 ** 2 - 0.5 * (z2 ** 2 + z3 ** 2 + z4**2)],
                          [x2 ** 2 - 0.5 * (x1 ** 2 + x3 ** 2 + x4**2) + y2 ** 2 - 0.5 * (y1 ** 2 + y3 ** 2 + y4**2)+ z2 ** 2 - 0.5 * (z1 ** 2 + z3 ** 2 + z4**2)],
                          [x3 ** 2 - 0.5 * (x1 ** 2 + x2 ** 2 + x4**2) + y3 ** 2 - 0.5 * (y1 ** 2 + y2 ** 2 + y4**2) + z3 ** 2 - 0.5 * (z1 ** 2 + z2 ** 2 + z4**2)]

                          ])
            xout = np.dot(np.linalg.inv(A),B).T
            xout = np.array([xout[0,0],xout[0,1],xout[0,2]])
            wigner_points.append(xout)

wigner_points_cleaned = []

for w_point in wigner_points:
    dist = []
    for point in point_list:
        dist.append(np.linalg.norm(w_point-point))
    dist = np.array(dist)
    if np.all(np.linalg.norm(w_point-origin) <= dist*1.01):
        wigner_points_cleaned.append(w_point)



mlab.points3d(origin[0],origin[1],origin[2],color=(0,0,1))


# for point in point_list:
coords = zip(*point_list)

mlab.points3d(coords[0],coords[1],coords[2],color=(0,0,0), scale_factor=1)


coords2 = zip(*wigner_points_cleaned)


mlab.points3d(coords2[0],coords2[1],coords2[2],color=(1,0,0))


mlab.show()
