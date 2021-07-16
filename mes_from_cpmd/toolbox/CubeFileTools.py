#!/usr/bin/python3
import numpy as np

#from lib import constants
from mes_from_cpmd.toolbox import constants
from mes_from_cpmd.toolbox import cube


def Overlap(array1, array2, differential_volume_element, complex_data=False):
    if complex_data:
        return np.dot(array1.flat, np.conj(array2).flat)*differential_volume_element
    else:
        return np.dot(array1.flat, array2.flat)*differential_volume_element


def Integration(array, differential_volume_element):
    return np.sum(array)*differential_volume_element


def OverlapMatrix(array1, array2, differential_volume_element, complex_data=False):
    if complex_data:
        overlap_matrix = np.zeros((array1.shape[0], array2.shape[0]), dtype=np.complex128)
    else:
        overlap_matrix = np.zeros((array1.shape[0], array2.shape[0]))
    for i1, d1 in enumerate(array1):
        for i2, d2 in enumerate(array2):
            overlap_matrix[i1,i2] = Overlap(d1, d2, differential_volume_element, complex_data=complex_data)
    return overlap_matrix



def LoadCellData(fn_cube):
    #if fn_cube[-5:] == '.cube':
    data, numbers, coords_au, cell_au, comment1, comment2, origin_au = cube.ReadCubeFile(fn_cube)
    # elif fn_cube[-3:] == '.gz':
    #     data, numbers, coords_au, cell_au, comment1, comment2, origin_au = cube.ReadGZIPCubeFile(fn_cube)
    # else:
    #     raise Exception('Unknown fromat?!')
    tmp = abs(cell_au[0,1])+abs(cell_au[0,2])+abs(cell_au[1,0])+abs(cell_au[1,2])+abs(cell_au[2,0])+abs(cell_au[2,1])
    if tmp > 1E-10:
        raise Exception('Only orthorhombic or simple cubic boxes supported!')
    cell_data = dict()
    cell_data['numbers'  ] = numbers
    cell_data['coords_au'] = coords_au
    cell_data['cell_au'  ] = cell_au
    cell_data['symbols'  ] = [constants.symbols[int(e)-1] for e in cell_data['numbers']]
    cell_data['species'  ] = list(set(cell_data['symbols']))
    cell_data['mesh'     ] = data.shape
    cell_data['d3r_au'   ] = cell_au[0,0]*cell_au[1,1]*cell_au[2,2]
    cell_data['volume_au'] = cell_data['d3r_au']*data.shape[0]*data.shape[1]*data.shape[2]
    cell_data['r_au'     ] = CalcGridPositions(cell_au, data.shape, origin_au=origin_au)
    cell_data['data'     ] = data
    cell_data['comment1' ] = comment1
    cell_data['comment2' ] = comment2
    cell_data['origin_au'] = origin_au
    return cell_data


def CalcGridPositions(cell_au, mesh, origin_au=None):
    """returns r_au of shape (mesh1,mesh2,mesh3,3) in units of cell, i.e. default AU"""
    n_x, n_y, n_z = mesh
    a_x, a_y, a_z = cell_au[0][0], cell_au[1][1], cell_au[2][2]
    if origin_au.all() == None:
        origin_au = np.array([n_x*a_x, n_y*a_y, n_z*a_z])/2
    x_au = np.arange(0, n_x*a_x, a_x) - origin_au[0]
    y_au = np.arange(0, n_y*a_y, a_y) - origin_au[1]
    z_au = np.arange(0, n_z*a_z, a_z) - origin_au[2]
    r_au = np.zeros((3, n_x, n_y, n_z)) #NB! Changed to row-major-order!
    r_au[0,:,:,:] = x_au[:,np.newaxis,np.newaxis]
    r_au[1,:,:,:] = y_au[np.newaxis,:,np.newaxis]
    r_au[2,:,:,:] = z_au[np.newaxis,np.newaxis,:]
    return r_au


def GenerateGrid(dimensions, origin_au, cell_au):
    """cell_au should be a diagonal matrix with a_lat_au/n_grid on the diagonal"""
    grid  = np.indices(dimensions, dtype=np.float)
    grid -= (origin_au/cell_au.diagonal())[:,np.newaxis,np.newaxis,np.newaxis]
    return np.sum(cell_au[:,:,np.newaxis,np.newaxis,np.newaxis]*grid, axis=0)


def Gaussian3D(r_loc_au, sigma_sqared_au):
    return (1/(2*np.pi*sigma_sqared_au))**(3/2)*np.exp(-0.5*np.sum(r_loc_au**2, axis=0)/sigma_sqared_au)

# def PBCMapping(r_au, a_lat_au):
#     return r_au - np.around(r_au/a_lat_au)*a_lat_au

def PBCMapping(r_au, cell_au):
    tmp = (cell_au.diagonal()*np.array(r_au.shape[1:]))[:,np.newaxis,np.newaxis,np.newaxis]
    return r_au - np.around(r_au/tmp)*tmp

def RadialProfile(data, r_abs_au, bin_width_au):
    r_abs_au      = np.floor(r_abs_au/bin_width_au).astype(np.int)
    tbin          = np.bincount(r_abs_au.ravel(), data.ravel())
    nr            = np.bincount(r_abs_au.ravel())
    radialprofile = tbin / nr
    return radialprofile, tbin, nr

def CartesianMoments(r_au, order=1):
    x = r_au[0,:,:,:]
    y = r_au[1,:,:,:]
    z = r_au[2,:,:,:]
    if order == 1:
        sh = np.zeros((3, r_au.shape[1], r_au.shape[2], r_au.shape[3]))
        sh[0] = x
        sh[1] = y
        sh[2] = z
    elif order == 2:
        sh = np.zeros((6, r_au.shape[1], r_au.shape[2], r_au.shape[3]))
        sh[0] = x**2                     
        sh[1] = x*y
        sh[2] = x*z
        sh[3] = y**2
        sh[4] = y*z
        sh[5] = z**2
    elif order == 3:
        sh = np.zeros((10, r_au.shape[1], r_au.shape[2], r_au.shape[3]))
        sh[0] = x**3
        sh[1] = x**2*y 
        sh[2] = x**2*z 
        sh[3] = y**2*x 
        sh[4] = y**3
        sh[5] = y**2*z
        sh[6] = z**2*x
        sh[7] = z**2*y
        sh[8] = z**3
        sh[9] = x*y*z
    elif order == 4:
        sh = np.zeros((15, r_au.shape[1], r_au.shape[2], r_au.shape[3])) #NB! 16 or 15 components ?! http://onlinelibrary.wiley.com/doi/10.1002/3527602747.app4/pdf
        sh[ 0] = x**4
        sh[ 1] = x**3*y
        sh[ 2] = x**3*z
        sh[ 3] = y**3*x
        sh[ 4] = y**4
        sh[ 5] = y**3*z
        sh[ 6] = z**3*x
        sh[ 7] = z**3*y
        sh[ 8] = z**4
        sh[ 9] = x**2*y**2
        sh[10] = x**2*y*z
        sh[11] = x**2*z**2
        sh[12] = y**2*z**2
        sh[13] = x*y**2*z
        sh[14] = x*y*z**2
    elif order == 5: # NB! not yet tested ...
        sh = np.zeros((21, r_au.shape[1], r_au.shape[2], r_au.shape[3]))
        sh[ 0] = x**5
        sh[ 1] = x**4*y
        sh[ 2] = x**4*z
        sh[ 3] = x**3*y**2
        sh[ 4] = x**3*z**2
        sh[ 5] = x**3*y*z
        sh[ 6] = y**5
        sh[ 7] = y**4*x
        sh[ 8] = y**4*z
        sh[ 9] = y**3*x**2
        sh[10] = y**3*z**2
        sh[11] = y**3*x*z
        sh[12] = z**5
        sh[13] = z**4*x
        sh[14] = z**4*y
        sh[15] = z**3*x**2
        sh[16] = z**3*y**2
        sh[17] = z**3*x*y
        sh[18] = x**2*y**2*z
        sh[19] = x**2*y*z**2
        sh[20] = x*y**2*z**2
    else:
        raise Exception('Order %d not supported!'%order)
    return sh


def PotentialOfChargeDistribution(origin_au, cell_au, mesh, charge_dist, d3r_au, verbose=False): # rewrite this in cython!!
    r_au = CalcGridPositions(cell_au, mesh, origin_au=origin_au)
    V = np.zeros(mesh)
    for i1 in range(mesh[0]):
        for i2 in range(mesh[1]):
            if verbose:
                print('PotentialOfChargeDistribution: %05d %05d'%(i1,i2), end='\r')
            for i3 in range(mesh[2]):
                tmp = V[i1,i2,i3]
                r_diff_au = np.sqrt(((r_au[:,i1,i2,i3,np.newaxis,np.newaxis,np.newaxis] - r_au)**2).sum(axis=0))
                V += charge_dist[i1,i2,i3]/r_diff_au
                V[i1,i2,i3]=tmp
    return V*d3r_au


def Potential(data, cell_au):
    from numpy.fft import ifftn,fftn
    n1, n2, n3 = data.shape
    a1, a2, a3 = tuple(cell_au.diagonal())
    R,K        = GetCell(n1, n2, n3, a1*n1, a2*n2, a3*n3)
    V_R        = ifftn(Vk(K)*fftn(data)).real
    return R, V_R

def Vk(k):
    """Fourier transform of Coulomb potential $1/r$"""
    with np.errstate(divide='ignore'):
        return np.where(k==0.0, 0.0, np.divide(4.0*np.pi, k**2))

def GetCell(n1, n2, n3, a1, a2, a3):
    from numpy.fft import fftfreq
    r1 = np.arange(n1)*(a1/n1)-a1/2
    k1 = 2*np.pi*fftfreq(n1,a1/n1)
    # r1 = np.arange(n1)*(a1/n1)-a1/2
    # k1 = 2*np.pi*fftfreq(n1,a1/n1)
    # r2 = np.arange(n2)*(a2/n2)-a2/2
    # k2 = 2*np.pi*fftfreq(n2,a2/n2)
    # r3 = np.arange(n3)*(a3/n3)-a3/2
    # k3 = 2*np.pi*fftfreq(n3,a3/n3)
    ix, iy, iz = (slice(None), None, None), (None, slice(None), None), (None, None, slice(None))
    (X, Kx), (Y, Ky), (Z, Kz) = [(r1[_i], k1[_i]) for _i in [ix, iy, iz]]
    R = np.sqrt(X**2 + Y**2 + Z**2)
    K = np.sqrt(Kx**2 + Ky**2 + Kz**2)
    return R,K


def GaussianChargeDistribution(position_au, charge, sigma_au, cell_au, mesh):
    r_au = CalcGridPositions(cell_au, mesh, origin_au=position_au)
    return charge*np.exp(-0.5*(r_au**2).sum(axis=0)/sigma_au**2)/(np.sqrt(2*np.pi)*sigma_au)**3

def GaussianPotential(position_au, charge, sigma_au, cell_au, mesh):
    import scipy.special
    # http://en.wikipedia.org/wiki/Poisson%27s_equation
    r_au   = CalcGridPositions(cell_au, mesh, origin_au=position_au)
    r      = np.sqrt((r_au**2).sum(axis=0))
    tup    = np.unravel_index(np.argmin(r), r.shape)
    c      = np.sqrt(2)*sigma_au
    erf    = scipy.special.erf(r/c)
    r[tup] = 1E10
    g      = erf/r
    g[tup] = 2/np.sqrt(np.pi)/c
    return charge*g
