import numpy as np
import cell_tree2d

from shapely.geometry.polygon import Polygon
from shapely.geometry.point import Point
from shapely.geometry import MultiPolygon
from shapely.geometry import MultiPoint

#### NOTES

# Cell interpolation is done using this algorithm: 
# https://www.particleincell.com/2012/quad-interpolation/
#
# First the transformation from a quad in space to a normal unit square.
# The affine transformation is done using the equations:
#     x = a1 + a2*l + a3*m + a4*l*m
#     y = b1 + b2*l + b3*m + b4*l*m
# constrained by
#     x1 -> l=0 ; y1 -> m=0    ( lower left point )
#     x2 -> l=1 ; y2 -> m=0    ( lower right point )
#     x3 -> l=1 ; y3 -> m=1    ( upper right point )
#     x4 -> l=0 ; y4 -> m=1    ( upper left point )
# Which leads to the matrix equation
#     
#     [x1    [[1., 0., 0., 0.]   [a1 
#      x2  =  [1., 0., 1., 0.]    a2 
#      x3     [1., 1., 1., 1.]    a3 
#      x4]    [1., 1., 0., 0.]]   a4]
#
# Note the transformation defined by the matrix (A) is defined
# by starting in the lower left and cycling positively (counterclockwise)
# around the quad.
_AI = np.linalg.inv( np.array([[1., 0., 0., 0.], 
                               [1., 0., 1., 0.], 
                               [1., 1., 1., 1.], 
                               [1., 1., 0., 0.]]) )

def compute_coeffs(squares_i):
    """
    Input:
    -----
    squares :   array,  dimension = (Nsquares, 4, 2)
                x, y coordinates of the squares in the grid.
                Order matters for the second dimension: br, tr, tl, bl
    
    Output:
    ------
    a, b :      array, array
                a, the coefficients in x  (what coefficients?)
                b, the coefficients in y  (what coefficients?)
    
    """
    a = np.dot(_AI, squares_i[:,:,0].T)   # no need for complex conjugate. px is always real.
    b = np.dot(_AI, squares_i[:,:,1].T)   #                           .....so is py
    return a, b


# In short, the generated weights for values at the four node points
# are the ratio of the area defined by area of the rectangle defined
# by the trial point and the oposite node normalized by the total
# cell area.

# Using np.einsum(...) is probaby going to be useful for
# constructing the interpolated values in the end.




# Used by compute coeffs


def locate_faces(points, grid):
    """
    Returns the node grid indices, one per point.

    Points that are not on the node grid will have an index of -1

    If a single point is passed in, a single index will be returned.
    If a sequence of points is passed in an array of indexes will be returned.

    :param points:  The points that you want to locate -- (lon, lat). If the shape of point
                    is 1D, function will return a scalar index. If it is 2D, it will return
                    a 1D array of indices.
    :type points: array-like containing one or more points: shape (2,) for one point, shape (N, 2)
                 for more than one point.

    :param grid: The grid on which you want to locate the points
    :type grid: Name of the grid ('node', 'center', 'edge1', 'edge2)

    This version utilizes the CellTree data structure.

    """

    points = np.asarray(points, dtype=np.float64)
    just_one = (points.ndim == 1)
    points = points.reshape(-1, 2)

    tree = build_celltree(grid)
    indices = tree.locate(points)
    lon, lat = self._get_grid_vars(grid)
    x = indices % (lat.shape[1] - 1)
    y = indices // (lat.shape[1] - 1)
    ind = np.column_stack((y, x))
    ind[ind[:, 0] == -1] = [-1, -1]
    if just_one:
        res = ind[0]
        return res
    else:
        res = np.ma.masked_less(ind, 0)
        if _memo:
            self._add_memo(points, res, grid, self._ind_memo_dict, _copy, _hash)
        return res


def get_alphas(x, y, a, b):
    """
    Params:
    a: x coefficients
    b: y coefficients
    x: x coordinate of points defining squares (possibly only query points)
    y: y coordinate of points defining squares (possibly only query points)

    Returns:
    (l,m) - coordinate in logical space to use for interpolation

    Eqns:
    m = (-bb +- sqrt(bb^2 - 4*aa*cc))/(2*aa)
    l = (l-a1 - a3*m)/(a2 + a4*m)
    """
    def lin_eqn(l, m, ind_arr, aa, bb, cc):
        """
        AB is parallel to CD...need to solve linear equation instead.
        m = -cc/bb
        bb = Ei*Fj - Ej*Fi + Hi*Gj - Hj*Gi
        k0 = Hi*Ej - Hj*Ei
        """
        m[ind_arr] = -cc / bb
        l[ind_arr] = (x[ind_arr] - a[0][ind_arr] - a[2][ind_arr]
                      * m[ind_arr]) / (a[1][ind_arr] + a[3][ind_arr] * m[ind_arr])

    def quad_eqn(l, m, ind_arr, aa, bb, cc):
        """

        """
        if len(aa) is 0:
            return
        k = bb * bb - 4 * aa * cc
        k = np.ma.array(k, mask=(k < 0))

        det = np.ma.sqrt(k)
        m1 = (-bb - det) / (2 * aa)
        l1 = (x[ind_arr] - a[0][ind_arr] - a[2][ind_arr] *
              m1) / (a[1][ind_arr] + a[3][ind_arr] * m1)

        m2 = (-bb + det) / (2 * aa)
        l2 = (x[ind_arr] - a[0][ind_arr] - a[2][ind_arr] *
              m2) / (a[1][ind_arr] + a[3][ind_arr] * m2)

        m[ind_arr] = m1
        l[ind_arr] = l1

        t1 = np.logical_or(l1 < 0, l1 > 1)
        t2 = np.logical_or(m1 < 0, m1 > 1)
        t3 = np.logical_or(t1, t2)

        l[ind_arr[t3]] = l2[t3]
        m[ind_arr[t3]] = m2[t3]

    aa = a[3] * b[2] - a[2] * b[3]
    bb = a[3] * b[0] - a[0] * b[3] + a[1] * \
        b[2] - a[2] * b[1] + x * b[3] - y * a[3]
    cc = a[1] * b[0] - a[0] * b[1] + x * b[1] - y * a[1]

    m = np.zeros(bb.shape)
    l = np.zeros(bb.shape)
    t = aa[:] == 0
    
    # These functions pass both l and m as pointers, modified in place.
    lin_eqn(l, m, np.where(t)[0], aa[t], bb[t], cc[t])
    quad_eqn(l, m, np.where(~t)[0], aa[~t], bb[~t], cc[~t])
    
    aa = 1 - l - m + l * m
    ab = m - l * m
    ac = l * m
    ad = l - l * m

    return np.array((aa, ab, ac, ad)).T, l, m


def array2grid(x, y):
    '''Return faces and nodes given x and y grid node positions. Order: br, tr, tl, bl
    
    Input
    -----
    x, y :  ndarray
    
    Output
    ------
    nodes : flattened xy pairs
    faces : indices of square nodes
    '''
    nodes = np.ascontiguousarray(np.column_stack((x[:].reshape(-1),
                                                  y[:].reshape(-1)))).astype(np.float64)
    j_len, i_len = x.shape
    # faces = np.array([np.array([[xi + i_len, xi + i_len + 1, xi + 1, xi]
    #                             for xi in range(0, i_len - 1, 1)]) + yi * i_len for yi in range(0, j_len - 1)])
    faces = np.array([np.array([[xi, xi+1, xi+i_len+1, xi+i_len]
                                for xi in range(0, i_len-1, 1)]) + yi * i_len for yi in range(0, j_len-1)])
    faces = np.ascontiguousarray(faces.reshape(-1, 4).astype(np.int32))
    
    return nodes, faces



# sample grid
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

x_nodes, y_nodes = x.flatten(), y.flatten()

x_centers = 0.25*(x[1:, 1:] + x[1:, :-1] + x[:-1, 1:] + x[:-1, :-1]).flatten()
y_centers = 0.25*(y[1:, 1:] + y[1:, :-1] + y[:-1, 1:] + y[:-1, :-1]).flatten()

# create nodes and faces
nodes, faces = array2grid(x, y)
squares = np.array([nodes[face] for face in faces])
ct = cell_tree2d.CellTree(nodes, faces)

# create a set of trial points
# xyi = np.random.randn(10, 2)                 # random points
xyi = np.vstack((x_centers, y_centers)).T       # centers of gridded data squares
# xyi = np.vstack((x_nodes, y_nodes)).T           # nodes of gridded data squares, ie., at data points

# want only squares that contain query points
gridpoint_indices = ct.locate(xyi)

# remove gridpoint indices with value -1, 
# which are outside the grid domain
inside = gridpoint_indices >= 0
gridpoint_indices = gridpoint_indices[inside]
xyi = xyi[inside]


squares_i = squares[gridpoint_indices]


# First a sanity check. Are the trial points really inside squares_i?

def test_cell_tree2d(squares, points):
    '''
    squares: an array of quads, full grid [Nsquares, 4, 2]
    points:  an array of trial points [Npoints, 2]
    ''' 
    # find squares that contain query points
    gridpoint_indices = ct.locate(points)

    # remove indices with value -1, outside the grid domain
    inside = gridpoint_indices >= 0    # good points, inside domain
    gridpoint_indices = gridpoint_indices[inside]
    points = points[inside]

    mesh = MultiPolygon([Polygon(p) for p in squares])
    trial = MultiPoint([Point(p) for p in xyi])
    
    contains = [m.contains(p) for m, p in zip(mesh, trial)]
    
    assert(np.alltrue(contains))

##### make up some trial points
# 
#
a, b = compute_coeffs(squares_i)

xi, yi = xyi.T
alphas, l, m = get_alphas(xi, yi, a, b)


def z(x, y):
    return x**2 * y**2 + np.sin(x)*np.cos(y)

z_verts = z(*nodes.T)[faces][gridpoint_indices]

zi = (z_verts * alphas).sum(axis=-1)

zi_true = z(xi, yi)

# assert( np.allclose(zi, zi_true) )


#######################
# edited to hear
######################

#
# reflats = points[:, 1]
# reflons = points[:, 0]
#
# l, m = x_to_l(reflons, reflats, a, b)
#
# aa = 1 - l - m + l * m
# ab = m - l * m
# ac = l * m
# ad = l - l * m
# alphas = np.array((aa, ab, ac, ad)).T
#
# if _memo:
#     self._add_memo(points, alphas, grid, self._alpha_memo_dict, _copy, _hash)
# return alphas
