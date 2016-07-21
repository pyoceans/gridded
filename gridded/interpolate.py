
# Used by compute coeffs
A = np.array(([1, 0, 0, 0], 
              [1, 0, 1, 0], 
              [1, 1, 1, 1], 
              [1, 1, 0, 0]))
AI = np.linalg.inv(A)

def compute_coeffs(px, py):
    """
    Input:
    -----
    px, py :    array
                x, y coordinates of the polygon. Order matters: br, tr, tl, bl
    
    Output:
    ------
    a, b :      array, array
                a, the alpha coefficients in x
                b, the alpha coefficients in y
    
    """
    a = np.dot(AI, px.T)   # no need for complex conjugate. px is always real.
    b = np.dot(AI, py.T)   # so is py
    return np.array(a), np.array(b)

def x_to_l(x, y, a, b):
    """
    Params:
    a: x coefficients
    b: y coefficients
    x: x coordinate of point
    y: y coordinate of point

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

    return (l, m)


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
    faces = np.array([np.array([[xi + i_len, xi + i_len + 1, xi + 1, xi]
                                for xi in range(0, i_len - 1, 1)]) + yi * i_len for yi in range(0, j_len - 1)])
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


# create nodes and faces
nodes, faces = array2grid(x, y)

squares = np.array([nodes[face] for face in faces])



#######################
# edited to hear
######################



a, b = compute_coeffs(polyx, polyy)

reflats = points[:, 1]
reflons = points[:, 0]

l, m = x_to_l(reflons, reflats, a, b)

aa = 1 - l - m + l * m
ab = m - l * m
ac = l * m
ad = l - l * m
alphas = np.array((aa, ab, ac, ad)).T

if _memo:
    self._add_memo(points, alphas, grid, self._alpha_memo_dict, _copy, _hash)
return alphas
