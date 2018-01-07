import scipy
from mayavi import mlab
import time

X = 100 * scipy.rand(300, 3)
figure = mlab.figure('myfig')
st = time.time()
figure.scene.disable_render = True  # Super duper trick that seems to have no effect here. Don't know why.
mlab.points3d(X[:, 0], X[:, 1], X[:, 2], scale_factor=0.4)
for i, x in enumerate(X):
    mlab.text3d(x[0], x[1], x[2], str(i), scale=(2, 2, 2))
figure.scene.disable_render = False  # Super duper trick

print(time.time()-st)
mlab.show()