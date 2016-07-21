
import numpy as np


def rot2d(x, y, ang):
    '''rotate vectors by geometric angle'''
    xr = x * np.s(ang) - y * np.sin(ang)
    yr = x * np.sin(ang) + y * np.cos(ang)
    return xr, yr


def generate_sgrid(x=None, y=None, rot=None):
    '''
    Generate the faces and notes for an sgrid

    Input
    -----
    x, y: slice (paramteres into the numpy.mgrid function)
    rot:  float

    Output
    ------
    nodes : flattened xy pairs
    faces : indices of squares
    '''
    x = x or slice(1, 10, 15j)
    y = y or slice(1, 20, 18j)
    rot = rot or 0.2

    # create a rotated Cartesian grid
    xc, yc = np.mgrid[x, y]
    yc = yc**1.2 + xc**1.5
    x, y = rot2d(xc, yc, rot)
    y /= 10.0
    x -= x.mean()
    y -= y.mean()

    return array_to_grid(x, y)


def generate_rgrid(x=None, y=None):
    '''
    Generate the faces and notes for an rgrid

    Input
    -----
    x, y: slice (paramteres into the numpy.mgrid function)

    Output
    ------
    nodes : flattened xy pairs
    faces : indices of squares
    '''
    x = x or slice(1, 10, 15j)
    y = y or slice(1, 20, 18j)

    # create a Cartesian grid
    xc, yc = np.mgrid[x, y]
    return array_to_grid(x, y)


def array_to_grid(x, y):
    '''
    Return faces and nodes given x and y grid node positions.
    Order: br, tr, tl, bl (counter-clockwise / positive)

    Input
    -----
    x, y :  ndarray

    Output
    ------
    nodes : flattened xy pairs
    faces : indices of squares
    '''

    # Create nodes and faces from grid
    nodes = np.ascontiguousarray(
        np.column_stack(
            (
                x[:].reshape(-1),
                y[:].reshape(-1)
            )
        )
    ).astype(np.float64)

    j_len, i_len = x.shape
    faces = np.array(
        [
            np.array(
                [
                    # Counter clockwise points (positive)
                    [ xi + i_len, xi + i_len + 1, xi + 1, xi ]
                    for xi in range(0, i_len - 1, 1)
                ]
            ) + yi * i_len for yi in range(0, j_len - 1)
        ]
    )
    faces = np.ascontiguousarray(
        faces.reshape(-1, 4).astype(np.int32)
    )

    return faces, nodes
