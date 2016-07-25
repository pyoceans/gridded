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
                Order matters for the second dimension: bl, br, tr, tl
                I.e., starting in the lower left and cycling CCW (positiv]]]]
    
    Output:
    ------
    a, b :      array, array
                a, the coefficients in x  (what coefficients?)
                b, the coefficients in y  (what coefficients?)
    
    """
    a = np.dot(_AI, squares_i[:,:,0].T)   # no need for complex conjugate. x is always real.
    b = np.dot(_AI, squares_i[:,:,1].T)   #                          .....so is y
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
        return res


def get_lm(x, y, a, b):
    """
    Input:
    -----
    a: x affine transformation coefficients
    b: y affine transformation coefficients
    x: x coordinate of points defining squares (possibly only query points)
    y: y coordinate of points defining squares (possibly only query points)

    Returns:
    weights - weights of grid nodes to use in interpolation
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

    def quad_eqn(l, m, ind_arr, aa, bb, cc, eps=1e-5):
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

        # the addition of eps is needed for numerical precision when the points are
        # on the cell edge. 
        t1 = np.logical_or(l1 < -eps, l1 > 1.0+eps)
        t2 = np.logical_or(m1 < -eps, m1 > 1.0+eps)
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
    
    return l, m

def get_weights(l, m):
    w1 = 1 - l - m + l * m
    w2 = m - l * m
    w3 = l * m
    w4 = l - l * m

    return np.array((w1, w2, w3, w4)).T


def array2grid(x, y):
    '''Return faces and nodes given x and y grid node positions. Order: br, tr, tl, bl
    
    Input
    -----
    x, y :  ndarray (each [Ny, Nx])
            The nodes defining the 2D grid.
    
    Output
    ------
    nodes : ndarray [Nnodes, 2]
            flattened xy pairs.    
    faces : ndarray [Nsquares, 4]
            connectivity indices of the nodes defining the grid squares. 
    '''
    nodes = np.ascontiguousarray( np.column_stack((x[:].reshape(-1),
                                                   y[:].reshape(-1))) ).astype(np.float64)
    j_len, i_len = x.shape

    # create square starting in the lower left, and continuing CCW (+)
    faces = np.array([np.array([[xi, xi+1, xi+i_len+1, xi+i_len]
                                for xi in range(0, i_len-1, 1)]) + yi * i_len for yi in range(0, j_len-1)])
    faces = np.ascontiguousarray(faces.reshape(-1, 4).astype(np.int32))
    
    return nodes, faces


class CellTree_interpolator(object):
    
    def __init__(self, squares, points, faces, gridpoint_indices):
        
        self.squares = squares
        self.points = points
        self.faces = faces
        self.gridpoint_indices = gridpoint_indices
        
        self.a, self.b = compute_coeffs(self.squares)
        xi, yi = self.points.T
    
        self.l, self.m = get_lm(xi, yi, self.a, self.b)
        self.alphas = get_weights(self.l, self.m)
    
    def interpolate(self, grid_z):
        # Interpolated values based on z at the verticies and the
        # calculated alpha weights.
        
        z = np.ascontiguousarray(grid_z[:].reshape(-1)).astype(np.float64)
        
        z_verts = z[self.faces][self.gridpoint_indices]

        return (z_verts * self.alphas).sum(axis=-1)
        

class CellTree(object):
    
    def __init__(self, grid_x, grid_y):
        self.nodes, self.faces = array2grid(grid_x, grid_y)
        self.squares = np.array([nodes[face] for face in faces])
        self.ct = cell_tree2d.CellTree(nodes, faces)
    
    def locate(self, points):
        
        # find the grid squares that contain the points
        gridpoint_indices = ct.locate(points)

        # remove gridpoint indices with value -1, 
        # which are outside the grid domain
        inside = gridpoint_indices >= 0
        points = points[inside]

        # get the corresponding squares to the points
        # inside the domain
        gridpoint_indices = gridpoint_indices[inside]
        squares = self.squares[gridpoint_indices]
        
        return CellTree_interpolator(squares, points, self.faces, gridpoint_indices)


def test_cell_tree2d(squares, points, eps=0.01):
    '''
    squares: an array of quads, full grid [Nsquares, 4, 2]
    points:  an array of trial points [Npoints, 2]
    ''' 
    # find squares that contain query points
    gridpoint_indices = ct.locate(points)

    # remove indices with value -1, outside the grid domain
    inside = gridpoint_indices >= 0    # good points, inside domain
    gridpoint_indices = gridpoint_indices[inside]

    points_i = points[inside]
    squares_i = squares[gridpoint_indices]
    
    mesh = MultiPolygon([Polygon(p).buffer(eps) for p in squares_i])
    trial = MultiPoint([Point(p) for p in points_i])
    
    contains = [m.contains(p) for m, p in zip(mesh, trial)]
        
    assert(np.alltrue(contains))
    
    return np.asarray(contains)


def test_interpolation(squares, points, zfunc):
    
    ct = CellTree(x, y)
    loc = ct.locate(points)
    zgrid = zfunc(x, y)
    zi = loc.interpolate(zgrid)
    print(zi.shape)
    
    # use loc.points, as this contains only the points in the domain.
    zi_true = zfunc(*loc.points.T)
    print(zi_true.shape)
    
    assert( np.allclose(zi, zi_true, rtol=1e-3) )
    

if __name__ == '__main__':
    
    def make_sample_grid(Ny, Nx):
        'return sample grid of dimension [Ny, Nx]'
        yc, xc = np.mgrid[1:10:Ny*1j, 1:20:Nx*1J]

        def rot2d(x, y, ang):
            '''rotate vectors by geometric angle'''
            xr = x*np.cos(ang) - y*np.sin(ang)
            yr = x*np.sin(ang) + y*np.cos(ang)
            return xr, yr

        x, y = rot2d(xc, (5+yc)**1.2*(3+xc)**0.3, 0.2)

        y /= y.ptp()/10.
        x /= x.ptp()/10.

        x -= x.mean()
        y -= y.mean()

        return x, y

    x, y = make_sample_grid(70, 50)

    # Some sample grid locations
    x_nodes, y_nodes = x.flatten(), y.flatten()

    x_u = 0.5*(x[:, 1:] + x[:, :-1]).flatten()
    y_u = 0.5*(y[:, 1:] + y[:, :-1]).flatten()

    x_v = 0.5*(x[1:, :] + x[:-1, :]).flatten()
    y_v = 0.5*(y[1:, :] + y[:-1, :]).flatten()

    x_centers = 0.25*(x[1:, 1:] + x[1:, :-1] + x[:-1, 1:] + x[:-1, :-1]).flatten()
    y_centers = 0.25*(y[1:, 1:] + y[1:, :-1] + y[:-1, 1:] + y[:-1, :-1]).flatten()

    # create nodes and faces
    nodes, faces = array2grid(x, y)
    squares = np.array([nodes[face] for face in faces])
    ct = cell_tree2d.CellTree(nodes, faces)
    

    def zfunc(x, y):
        'Sample field for interpolation'
        return np.sin(x/10.) + np.cos(y/10.)

    # create a set of trial points
    points = 10*np.random.randn(10000, 2)                 # random points
    test_cell_tree2d(squares, points)
    test_interpolation(squares, points, zfunc)

    points = np.vstack((x_centers, y_centers)).T   # cell centers
    test_cell_tree2d(squares, points)
    test_interpolation(squares, points, zfunc)

    points = np.vstack((x_nodes, y_nodes)).T       # cell corners
    test_cell_tree2d(squares, points)
    test_interpolation(squares, points, zfunc)
    
    points = np.vstack((x_u, y_u)).T               # x-direction edges
    test_cell_tree2d(squares, points)
    test_interpolation(squares, points, zfunc)

    points = np.vstack((x_v, y_v)).T               # y-direction edges
    test_cell_tree2d(squares, points)
    test_interpolation(squares, points, zfunc)
