import numpy as np
import matplotlib.pyplot as plt

from shapely.geometry.polygon import Polygon
from shapely.geometry import MultiPolygon

import cell_tree2d

# create a rotated Cartesian grid
xc, yc = np.mgrid[1:10:15j, 1:20:18j]
yc = yc**1.2 + xc**1.5

def rot2d(x, y, ang):
    '''rotate vectors by geometric angle'''
    xr = x*np.cos(ang) - y*np.sin(ang)
    yr = x*np.sin(ang) + y*np.cos(ang)
    return xr, yr

x, y = rot2d(xc, yc, 0.2)
y /= 10.0

x -= x.mean()
y -= y.mean()


# Create nodes and faces from grid
nodes = np.ascontiguousarray(np.column_stack((x[:].reshape(-1),
                                              y[:].reshape(-1)))).astype(np.float64)
y_size = x.shape[0]
x_size = y.shape[1]
faces = np.array([np.array([[xi, xi + 1, xi + x_size + 1, xi + x_size]
                            for xi in range(0, x_size - 1, 1)]) + yi * x_size for yi in range(0, y_size - 1)])
faces = np.ascontiguousarray(faces.reshape(-1, 4).astype(np.int32))

squares = [nodes[face] for face in faces]

## Convert to a bunch of shapely Polygon objects, for some unknown use.
# mesh = MultiPolygon([Polygon(p) for p in squares])

## Extra functions for plotting the grid
# for square in squares:
#      x, y = square.T
#      plt.fill(x, y)
#
# plt.gca().set_aspect(1.0)
# plt.show()


# Create some trial points and locate them using cell_tree
xyi = np.random.randn(10, 2)

ct = cell_tree2d.CellTree(nodes, faces)
idx = ct.locate(xyi)

