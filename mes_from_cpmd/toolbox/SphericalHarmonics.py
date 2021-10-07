#!/usr/bin/python3
import numpy as np
import gc

from mes_from_cpmd.toolbox import CubeFileTools

def ComplexSphericalHarmonics(r_au, order=1):
    # http://en.wikipedia.org/wiki/Table_of_spherical_harmonics
    x    = r_au[0,:,:,:]
    y    = r_au[1,:,:,:]
    z    = r_au[2,:,:,:]
    xmiy = x-y*1j
    xpiy = x+y*1j
    z2   = z**2
    r2   = x**2+y**2+z2
    r    = np.sqrt(r2)
    pi1  = np.pi
    pi2  = 2*pi1
    r_nd = r[:]
    tup  = np.unravel_index(np.argmin(r_nd), r_nd.shape)
    if r[tup] < 1E-5:
        r_nd[tup] = 1E10
    if order == 0:
        sh = +1/2*np.sqrt(1/pi1)*np.ones((1, r_au.shape[1], r_au.shape[2], r_au.shape[3]))
    elif order == 1:
        sh = np.zeros((3, r_au.shape[1], r_au.shape[2], r_au.shape[3]), dtype=np.complex128)
        sh[0] = +1/2*np.sqrt(3/pi2)*xmiy
        sh[1] = +1/2*np.sqrt(3/pi1)*z
        sh[2] = -1/2*np.sqrt(3/pi2)*xpiy
        sh /= r_nd
    elif order == 2:
        sh = np.zeros((5, r_au.shape[1], r_au.shape[2], r_au.shape[3]), dtype=np.complex128)
        sh[0] = +1/4*np.sqrt(15/pi2)*xmiy**2
        sh[1] = +1/2*np.sqrt(15/pi2)*xmiy*z
        sh[2] = +1/4*np.sqrt( 5/pi1)*(3*z2-r2)
        sh[3] = -1/2*np.sqrt(15/pi2)*xpiy*z
        sh[4] = +1/4*np.sqrt(15/pi2)*xpiy**2
        sh /= r_nd**2
    elif order == 3:
        sh = np.zeros((7, r_au.shape[1], r_au.shape[2], r_au.shape[3]), dtype=np.complex128)
        sh[0] = +1/8*np.sqrt( 35/pi1)*xmiy**3
        sh[1] = +1/4*np.sqrt(105/pi2)*xmiy**2*z
        sh[2] = +1/8*np.sqrt( 21/pi1)*xmiy*(5*z2-r2)
        sh[3] = +1/4*np.sqrt(  7/pi1)*(5*z2-3*r2)*z
        sh[4] = -1/8*np.sqrt( 21/pi1)*xpiy*(5*z2-r2)
        sh[5] = +1/4*np.sqrt(105/pi2)*xpiy**2*z
        sh[6] = -1/8*np.sqrt( 35/pi1)*xpiy**3
        sh /= r_nd**3
    elif order == 4:
        sh = np.zeros((9, r_au.shape[1], r_au.shape[2], r_au.shape[3]), dtype=np.complex128)
        sh[0] = +3/16*np.sqrt(35/pi2)*xmiy**4
        sh[1] = +3/ 8*np.sqrt(35/pi1)*xmiy**3*z
        sh[2] = +3/ 8*np.sqrt( 5/pi2)*xmiy**2*(7*z2-r2)
        sh[3] = +3/ 8*np.sqrt( 5/pi1)*xmiy*(7*z2-3*r2)*z
        sh[4] = +3/16*np.sqrt( 1/pi1)*(35*z2**2-30*z2*r2+3*r2**2)
        sh[5] = -3/ 8*np.sqrt( 5/pi1)*xpiy*(7*z2-3*r2)*z
        sh[6] = +3/ 8*np.sqrt( 5/pi2)*xpiy**2*(7*z2-r2)
        sh[7] = -3/ 8*np.sqrt(35/pi1)*xpiy**3*z
        sh[8] = +3/16*np.sqrt(35/pi2)*xpiy**4
        sh /= r_nd**4
    elif order == 5:
        sh = np.zeros((11, r_au.shape[1], r_au.shape[2], r_au.shape[3]), dtype=np.complex128)
        sh[ 0] = +3/32*np.sqrt(  77/pi1)*xmiy**5
        sh[ 1] = +3/16*np.sqrt( 385/pi2)*xmiy**4*z
        sh[ 2] = +1/32*np.sqrt( 385/pi1)*xmiy**3*(9*z2-r2)
        sh[ 3] = +1/ 8*np.sqrt(1155/pi2)*xmiy**2*(3*z2-r2)*z
        sh[ 4] = +1/16*np.sqrt( 165/pi2)*xmiy*(21*z2**2-14*r2*z2+r2**2)
        sh[ 5] = +1/16*np.sqrt(  11/pi1)*(63*z2**2-70*z2*r2+15*r2**2)*z
        sh[ 6] = -1/16*np.sqrt( 165/pi2)*xpiy*(21*z2**2-14*r2*z2+r2**2)
        sh[ 7] = +1/ 8*np.sqrt(1155/pi2)*xpiy**2*(3*z2-r2)*z
        sh[ 8] = -1/32*np.sqrt( 385/pi1)*xpiy**3*(9*z2-r2)
        sh[ 9] = +3/16*np.sqrt( 385/pi2)*xpiy**4*z
        sh[10] = -3/32*np.sqrt(  77/pi1)*xpiy**5
        sh /= r_nd**5
    elif order == 6:
        sh = np.zeros((13, r_au.shape[1], r_au.shape[2], r_au.shape[3]), dtype=np.complex128)
        sh[ 0] = +1/64*np.sqrt(3003/pi1)*xmiy**6
        sh[ 1] = +3/32*np.sqrt(1001/pi1)*xmiy**5*z
        sh[ 2] = +3/32*np.sqrt(  91/pi2)*xmiy**4*(11*z2-r2)
        sh[ 3] = +1/32*np.sqrt(1365/pi1)*xmiy**3*(11*z2-3*r2)*z
        sh[ 4] = +1/64*np.sqrt(1365/pi1)*xmiy**2*(33*z2**2-18*z2*r2+r2**2)
        sh[ 5] = +1/16*np.sqrt( 273/pi2)*xmiy**1*(33*z2**2-30*z2*r2+5*r2**2)*z
        sh[ 6] = +1/32*np.sqrt(  13/pi1)*(231*z2**3-315*z2**2*r2+105*z2*r2**2-5*r2**3)
        sh[ 7] = -1/16*np.sqrt( 273/pi2)*xpiy**1*(33*z2**2-30*z2*r2+5*r2**2)*z
        sh[ 8] = +1/64*np.sqrt(1365/pi1)*xpiy**2*(33*z2**2-18*z2*r2+r2**2)
        sh[ 9] = -1/32*np.sqrt(1365/pi1)*xpiy**3*(11*z2-3*r2)*z
        sh[10] = +3/32*np.sqrt(  91/pi2)*xpiy**4*(11*z2-r2)
        sh[11] = -3/32*np.sqrt(1001/pi1)*xpiy**5*z
        sh[12] = +1/64*np.sqrt(3003/pi1)*xpiy**6
        sh /= r_nd**6
    else:
        raise Exception('Order %d not supported!'%order)
    return sh


def ComplexRegularSolidHarmonics(r_au, order=1):
    sh  = ComplexSphericalHarmonics(r_au, order=order)
    r   = np.sqrt(np.sum(r_au**2, axis=0))
    fac = np.sqrt(4*np.pi/(2*order+1))
    return fac*r**order*sh


def ComplexIrregularSolidHarmonics(r_au, order=1):
    sh     = ComplexSphericalHarmonics(r_au, order=order)
    r      = np.sqrt(np.sum(r_au**2, axis=0))
    fac    = np.sqrt(4*np.pi/(2*order+1))
    tup    = np.unravel_index(np.argmin(r), r.shape)
    if r[tup] < 1E-5:
        r[tup] = 1E10
    ish    = fac*sh/r**(order+1)
    return ish


def RealSphericalHarmonicsTransformation(order):
    # http://en.wikipedia.org/wiki/Spherical_harmonics
    # http://en.wikipedia.org/wiki/Solid_harmonics
    if order == 0:
        return np.ones(1).reshape((1,1))
    dim = 3+(order-1)*2
    u_mat = np.zeros((dim, dim), dtype=np.complex128)
    u_mat[order, order] = 1
    for m in range(1, order+1):
        m_p = order+m
        m_m = order-m
        u_mat[m_p,m_p] = +np.sqrt(0.5)*(-1)**m
        u_mat[m_m,m_p] = -np.sqrt(0.5)*(-1)**m*1j
        u_mat[m_p,m_m] = +np.sqrt(0.5)
        u_mat[m_m,m_m] = +np.sqrt(0.5)*1j
    return u_mat


def RealSphericalHarmonics(r_au, order=1):
    c_sh    = ComplexSphericalHarmonics(r_au, order=order)
    u_trans = RealSphericalHarmonicsTransformation(order)
    return np.real(np.dot(c_sh.T, u_trans.T).T) # compare http://docs.scipy.org/doc/numpy/reference/generated/numpy.dot.html


def RealRegularSolidHarmonics(r_au, order=1):
    c_rsh   = ComplexRegularSolidHarmonics(r_au, order=order)
    u_trans = RealSphericalHarmonicsTransformation(order)
    return np.real(np.dot(c_rsh.T, u_trans.T).T) # compare http://docs.scipy.org/doc/numpy/reference/generated/numpy.dot.html


def RealIrregularSolidHarmonics(r_au, order=1):
    c_ish   = ComplexIrregularSolidHarmonics(r_au, order=order)
    u_trans = RealSphericalHarmonicsTransformation(order)
    return np.real(np.dot(c_ish.T, u_trans.T).T) # compare http://docs.scipy.org/doc/numpy/reference/generated/numpy.dot.html



def SphericalHarmonicsDebug(min_order, max_order):
    # setup mesh and gaussian for decay
    n_x = 200
    a_x = 10
    dx = a_x/n_x
    d3r_loc = dx**3
    cell_loc = np.diag([dx]*3)
    mesh_loc = (n_x,n_x,n_x)
    origin_loc = np.array([n_x//2]*3)*cell_loc[0,0]
    r_loc = CubeFileTools.CalcGridPositions(cell_loc, mesh_loc, origin_au=origin_loc)
    
    r = np.sqrt(np.sum(r_loc**2, axis=0))
    tup = np.unravel_index(np.argmin(r), r.shape)
    
    gaussian  = np.sqrt(np.exp(-0.5*r**2))
    gaussian /= np.sqrt(CubeFileTools.Overlap(gaussian, gaussian, d3r_loc))
    gaussian *= np.sqrt(4*np.pi)
    #print(gaussian[0,0,0]) this should be dacayed!

    # output utility ...
    def PrintMatrix(mat, print_format=' % 5.3F'):
        for i in range(mat.shape[0]):
            for j in range(mat.shape[1]):
                print(print_format%mat[i,j], end='')
            print()
        print()

    # memory optimized overlap calculation
    ovlps = list()
    for i_order in range(min_order,max_order+1):
        tmp1 = ComplexSphericalHarmonics(r_loc, order=i_order)
        for j_order in range(i_order,max_order+1):
            if j_order == i_order:
                tmp2 = tmp1
            else:
                tmp2 = ComplexSphericalHarmonics(r_loc, order=j_order)
            ovlp_matrix = CubeFileTools.OverlapMatrix(tmp1*gaussian, tmp2*gaussian, d3r_loc, complex_data=True)
            ovlps.append(ovlp_matrix)
            del tmp2
            gc.collect()
        del tmp1
        gc.collect()

    print('debug complex spherical harmonics')
    for ovl in ovlps:
        dim1, dim2 = ovl.shape
        if dim1 == dim2:
            PrintMatrix(np.abs(np.real(ovl)))

    for ovl in ovlps:
        dim1, dim2 = ovl.shape
        if dim1 == dim2:
            tmp_r = np.linalg.norm(np.diag(np.ones(dim1)) - np.real(ovl))
            tmp_i = np.linalg.norm(np.imag(ovl))
            print('Diag:   % 5.3E + i*% 5.3E'%(tmp_r, tmp_i))
        else:
            tmp_r = np.linalg.norm(np.real(ovl))
            tmp_i = np.linalg.norm(np.imag(ovl))
            print('OffD-r: % 5.3E + i*% 5.3E'%(tmp_r, tmp_i))

    print('debug real spherical harmonics')
    r_ovlps = list()
    for i_order in range(min_order,max_order+1):
        r_sh = RealSphericalHarmonics(r_loc, order=i_order)
        ovlp_matrix = CubeFileTools.OverlapMatrix(r_sh*gaussian, r_sh*gaussian, d3r_loc, complex_data=False)
        r_ovlps.append(ovlp_matrix)
        del r_sh
        gc.collect()

    for ovl in r_ovlps:
        PrintMatrix(np.abs(np.real(ovl)))

    print('debug real regular solid harmoncis')
    for i_order in range(min_order,min(4, max_order)+1):
        r_sh1 = RealRegularSolidHarmonics(r_loc, order=i_order)
        r_sh2 = RealRegularSolidHarmonics_LEGACY(r_loc, order=i_order)
        print('%02d: % 5.3E'%(i_order, np.linalg.norm(r_sh1-r_sh2)))


def RealRegularSolidHarmonics_LEGACY(r_au, order=1):
    x = r_au[0,:,:,:]
    y = r_au[1,:,:,:]
    z = r_au[2,:,:,:]
    if order == 0:
        sh = np.ones((1, r_au.shape[1], r_au.shape[2], r_au.shape[3]))
    elif order == 1:
        sh = np.zeros((3, r_au.shape[1], r_au.shape[2], r_au.shape[3]))
        sh[0] = y                                                        # 1 -1 = y
        sh[1] = z                                                        # 1  0 = z
        sh[2] = x                                                        # 1 +1 = x
    elif order == 2:
        sh = np.zeros((5, r_au.shape[1], r_au.shape[2], r_au.shape[3]))
        sh[0] = np.sqrt(3)*x*y                                           # 2 -2 = sqrt(3) xy
        sh[1] = np.sqrt(3)*y*z                                           # 2 -1 = sqrt(3) yz
        sh[2] = 0.5*(2*z**2-x**2-y**2)                                   # 2  0 = 0.5 (2z²-x²-y²)
        sh[3] = np.sqrt(3)*x*z                                           # 2  1 = sqrt(3) xz
        sh[4] = np.sqrt(3)*(x**2-y**2)/2                                 # 2  2 = sqrt(3) (x²-y²)/2
    elif order == 3:
        sh = np.zeros((7, r_au.shape[1], r_au.shape[2], r_au.shape[3]))
        sh[0] = np.sqrt(5/8)*(3*x**2-y**2)*y                             # 3 -3 = sqrt(5/8) (3x²-y²)y
        sh[1] = np.sqrt(15)*x*y*z                                        # 3 -2 = sqrt(15) xyz
        sh[2] = np.sqrt(3/8)*(4*z**2-x**2-y**2)*y                        # 3 -1 = sqrt(3/8) (4z²-x²-y²)y
        sh[3] = 0.5*(2*z**2-3*x**2-3*y**2)*z                             # 3  0 = 0.5 (2z²-3x²-3y²)z
        sh[4] = np.sqrt(3/8)*(4*z**2-x**2-y**2)*x                        # 3  1 = sqrt(3/8) (4z²-x²-y²)x
        sh[5] = np.sqrt(15/4)*(x**2-y**2)*z                              # 3  2 = sqrt(15/4) (x²-y²)z
        sh[6] = np.sqrt(5/8)*(x**2-3*y**2)*x                             # 3  3 = sqrt(5/8) (x²-3y²)x
    elif order == 4:
        r2 = x**2+y**2+z**2
        sh = np.zeros((9, r_au.shape[1], r_au.shape[2], r_au.shape[3]))
        sh[0] = 0.5*np.sqrt(35)*x*y*(x**2-y**2)                          # 4 -4 = 0.5 \sqrt{35} xy(x^2-y^2)\\
        sh[1] = 0.5*np.sqrt(35/2)*(3*x**2-y**2)*y*z                      # 4 -3 = 0.5 \sqrt{35/2} (3x^2-y^2)yz\\
        sh[2] = 0.5*np.sqrt(5)*x*y*(7*z**2-r2)                           # 4 -2 = 0.5 \sqrt{5} xy(7z^2-r^2)\\
        sh[3] = 0.5*np.sqrt(5/2)*y*z*(7*z**2-3*r2)                       # 4 -1 = 0.5 \sqrt{5/2} yz(7z^2-3r^2)\\
        sh[4] = 0.125*(35*z**4-30*z**2*r2+3*r2**2)                        # 4  0 = 0.125 (35z^4-30z^2r^2+3r^4)\\
        sh[5] = 0.5*np.sqrt(5/2)*x*z*(7*z**2-3*r2)                       # 4 +1 = 0.5 \sqrt{5/2} xz(7z^2-3r^2)\\
        sh[6] = 0.25*np.sqrt(5)*(x**2-y**2)*(7*z**2-r2)                  # 4 +2 = 0.25 \sqrt{5} (x^2-y^2)(7z^2-r^2)\\
        sh[7] = 0.5*np.sqrt(35/2)*(x**2-3*y**2)*x*z                      # 4 +3 = 0.5 \sqrt{35/2} (x^2-3y^2)xz\\
        sh[8] = 0.125*np.sqrt(35)*(x**2*(x**2-3*y**2)-y**2*(3*x**2-y**2)) # 4 +4 = 0.125 \sqrt{35} (x^2(x^2-3y^2)-y^2(3x^2-y^2))
    else:
        raise Exception('Order %d not supported!'%order)
    return sh


def RealRegularSolidHarmonics_EXTENDED(r_au, order=1):
    rh = RealRegularSolidHarmonics(r_au, order=order)
    r2 = r_au[0,:,:,:]**2+r_au[1,:,:,:]**2+r_au[2,:,:,:]**2
    if order == 1:
        return rh
    elif order == 2:
        rh = np.vstack((rh, np.array([r2])))
    elif order == 3:
        rh = np.vstack((rh, np.array([r_au[0,:,:,:]*r2,  r_au[1,:,:,:]*r2, r_au[2,:,:,:]*r2])))
    else:
        raise Exception('Order %d not supported!'%order)
    return rh


def NumberOfStates(order, traceless=False):
    if traceless:
        n_states_in_order = 2*order + 1
        n_states_to_order = order*(order+2)
    else:
        n_states_in_order = int((order+1)*(order+2)/2)
        n_states_to_order = int(order*(order+1)*(2*order+10)/12 + order) # simpified version of int(order*(order+1)*(2*order+1)/12 + 3*order*(order+1)/4 + order
    return n_states_in_order, n_states_to_order
