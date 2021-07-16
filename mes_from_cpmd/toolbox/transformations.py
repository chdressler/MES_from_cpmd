#!/usr/bin/python3
import numpy as np
import numpy.random
import scipy.misc
import copy

from mes_from_cpmd.toolbox import constants, SphericalHarmonics


def GetMasses(symbols):
    return [constants.species[s]['MASS'] for s in symbols]


def GetCharges(symbols):
    return [constants.species[s]['Z'] for s in symbols]


def GetValenceCharges(symbols):
    return [constants.species[s]['ZV'] for s in symbols]

def CentersOfMass(traj_pos_aa, masses):
    if len(traj_pos_aa.shape) == 2:
        return np.sum(traj_pos_aa*np.array(masses)[:,np.newaxis], axis=0)/sum(masses)
    elif len(traj_pos_aa.shape) == 3:
        return np.sum(traj_pos_aa*np.array(masses)[np.newaxis,:,np.newaxis], axis=1)/sum(masses)
    else:
        raise Exception('Shape of traj_pos_aa is not supported!')

def CentersOfCharge(trajectory, charges): # use with care for neutral systems ...
    centers_of_charge = np.zeros((trajectory.shape[0],3))
    for i, Z in enumerate(charges):
        centers_of_charge += abs(Z)*trajectory[:,i,:]
    return centers_of_charge/np.sum(np.abs(charges))


def LinearMomenta(velocities, masses):
    linmoms = np.zeros((velocities.shape[0], 3))
    for i, mass in enumerate(masses):
        linmoms += velocities[:,i]*mass
    return linmoms


def CentersOfMassMotion(velocities, masses):
    return LinearMomenta(velocities, masses)/sum(masses)


def CalcAngMom(coords, velocities, masses):
    n_atoms = masses.shape[0]
    angmom = np.zeros(3)
    for atom in range(n_atoms):
        angmom += np.cross(coords[atom], velocities[atom])*masses[atom]
    return angmom


def RotateCoordinates(RMAT, coords):
    c = copy.copy(coords)
    for i in range(c.shape[0]):
        c[i] = np.dot(RMAT, coords[i])
    return c

def Rotate3RTensor(RMAT, third_rank_tensor):
    tmp = np.zeros(third_rank_tensor.shape)
    for i in range(tmp.shape[0]):
        tmp[i] = np.dot(np.dot(RMAT, third_rank_tensor[i]), RMAT.T)
    return tmp


def NearestNeighborImage(vec_aa, cell_aa):
    check = False
    for i in range(3):
        if abs(vec_aa[i]) > cell_aa[i]/2:
            check = True
            if vec_aa[i] < 0:
                vec_aa[i] += cell_aa[i]
            elif vec_aa[i] > 0:
                vec_aa[i] -= cell_aa[i]
    if check:
        vec_aa = NearestNeighborImage(vec_aa, cell_aa)
    return vec_aa


def WrapMoleculePBC(coords_aa, cell_aa):
    """expects coords_aa.shape == (n_frames, n_atoms, 3)"""
    for i_frame in range(coords_aa.shape[0]):
        for i_atom in range(1, coords_aa.shape[1]):
            vec_aa = NearestNeighborImage(coords_aa[i_frame, i_atom] - coords_aa[i_frame, i_atom-1], cell_aa)
            coords_aa[i_frame, i_atom] += vec_aa - (coords_aa[i_frame, i_atom] - coords_aa[i_frame, i_atom-1])
    return coords_aa


def EulerAngles(x, y, z, convention='zyz', verbose=False):
    x /= np.linalg.norm(x)
    y /= np.linalg.norm(y)
    z /= np.linalg.norm(z)
    # http://en.wikipedia.org/wiki/Euler_angles#Rotation_matrix
    # http://en.wikipedia.org/wiki/Gimbal_lock
    if convention == 'zyz':
        c2 = z[2]
        beta = np.arccos(c2)
        if c2 > +(1-1E-10):
            if verbose:
                print('Gimbal Lock: z -> +Z')
            alpha = 0.0
            gamma = np.arctan2(x[1],x[0])-alpha
        elif c2 < -(1-1E-10):
            if verbose:
                print('Gimbal Lock: z -> -Z')
            alpha = 0.0
            gamma = alpha-np.arctan2(-x[1],-x[0])
        else:
            s2 = np.sqrt(1-c2**2)
            alpha = np.arctan2(z[1],+z[0])
            gamma = np.arctan2(y[2],-x[2])
    else:
        raise Exception('Convention %s not supported!'%convention)
    return np.array([alpha, beta, gamma])
        

def RotationMatrix(euler_angles, convention='zyz', verbose=False):
    # http://en.wikipedia.org/wiki/Euler_angles#Rotation_matrix
    s = np.sin(euler_angles)
    c = np.cos(euler_angles)
    R = np.zeros((3,3))
    if convention == 'zyz':
        R[0,0] = c[0]*c[1]*c[2] - s[0]*s[2]
        R[0,1] =-c[2]*s[0] - c[0]*c[1]*s[2]
        R[0,2] = c[0]*s[1]
        R[1,0] = c[0]*s[2] + c[1]*c[2]*s[0]
        R[1,1] = c[0]*c[2] - c[1]*s[0]*s[2]
        R[1,2] = s[0]*s[1]
        R[2,0] =-c[2]*s[1]
        R[2,1] = s[1]*s[2]
        R[2,2] = c[1]
    else:
        raise Exception('Convention %s not supported!'%convention)
    return R


def TestRotations(n_tests, test_gimbal_lock=False, verbose=False):
    tmp    = numpy.random.uniform(size=6*n_tests).reshape((2,n_tests,3))
    errors = list()
    for x,y in zip(tmp[0], tmp[1]):
        if test_gimbal_lock:
            x[2] = 0
            y[2] = 0
        x /= np.linalg.norm(x)
        y -= x*np.dot(x,y)
        y /= np.linalg.norm(y)
        z = np.cross(x,y)
        xyz = np.vstack((x,y,z)).T
        euler_angles = EulerAngles(x, y, z, verbose=verbose)
        R = RotationMatrix(euler_angles, verbose=verbose)
        errors.append(np.linalg.norm(R-xyz))
    return np.array(errors)


def EulerAnglesOfRotationMatrix(R):
    return EulerAngles(R.T[0], R.T[1], R.T[2])


def EulerAnglesOfInverseRotation(euler_angles):
    R = RotationMatrix(euler_angles)
    return EulerAngles(R[0], R[1], R[2])


def Wigner_d_Matrix(beta, l):
    """http://en.wikipedia.org/wiki/Wigner_D-matrix"""
    if l == 0:
        return np.ones(1).reshape((1,1))
    dim = 3+2*(l-1)
    sbh = np.sin(beta/2)
    cbh = np.cos(beta/2)
    fac = scipy.misc.factorial(np.arange(0,4*l+1))
    wigner_d_mat = np.zeros((dim, dim))
    for m1 in range(-l,l+1):
        for m2 in range(-l,l+1):
            pref = np.sqrt(fac[l+m1]*fac[l-m1]*fac[l+m2]*fac[l-m2])
            for s in range(0,2*l+1):
                if l+m1-s >= 0 and l-m2-s >= 0 and m2-m1+s >= 0:
                    #print(s,l+m1-s,l-m2-s,m2-m1+s)
                    sign  = (-1)**(m2-m1+s)
                    denom = fac[l+m1-s]*fac[l-m2-s]*fac[m2-m1+s]*fac[s]
                    nom   = cbh**(2*(l-s)+m1-m2)*sbh**(2*s-m1+m2)
                    wigner_d_mat[l+m2,l+m1] += sign*pref*nom/denom
    return wigner_d_mat


def Wigner_D_Matrix(euler_angles, l):
    """http://en.wikipedia.org/wiki/Wigner_D-matrix"""
    if l == 0:
        return np.ones(1).reshape((1,1))
    alpha, beta, gamma = tuple(euler_angles)
    wigner_d_mat = Wigner_d_Matrix(beta, l)
    wigner_D_mat = np.zeros(wigner_d_mat.shape, dtype=np.complex128)
    for m1 in range(-l,l+1):
        for m2 in range(-l,l+1):
            wigner_D_mat[l+m2,l+m1] = wigner_d_mat[l+m2,l+m1]*np.exp(1j*(m1*gamma+m2*alpha)) #NB! no minus in the exponent here?!
    return wigner_D_mat


def RealMomentRotationMatrix(euler_angles, order):
    D = Wigner_D_Matrix(euler_angles, order)
    u_sh = SphericalHarmonics.RealSphericalHarmonicsTransformation(order)
    return np.real(np.dot(u_sh, np.dot(D, np.conj(u_sh.T))))


def ComplexMomentRotationOrder(euler_angles, moments, order):
    D = Wigner_D_Matrix(euler_angles, order)
    return np.dot(D, moments)


def RealMomentRotationOrder(euler_angles, moments, order):
    D_prime = RealMomentRotationMatrix(euler_angles, order)
    return np.dot(D_prime, moments)


def ComplexMomentRotation(euler_angles, moments):
    rot_moms = list()
    for i_order in range(0,len(moments)):
        rot_moms.append(ComplexMomentRotationOrder(euler_angles, moments[i_order], i_order))
    return rot_moms


def RealMomentsToComplexMoments(real_moments, order):
    u_sh = SphericalHarmonics.RealSphericalHarmonicsTransformation(order)
    return np.dot(np.conj(u_sh.T), real_moments)


def ComplexMomentsToRealMoments(comp_moments, order):
    u_sh = SphericalHarmonics.RealSphericalHarmonicsTransformation(order)
    return np.real(np.dot(u_sh, comp_moments))


# statewise functions

def StatewiseRealComplexMomentRotationOrder_SLOW(euler_angles, state_moments, order):
    # its not performat to calculate the D matrix for each state ....
    rot_moms = list()
    for i_state in range(state_moments.shape[0]):
        rot_moms.append(ComplexMomentRotationOrder(euler_angles, state_moments[i_state], order))
    return np.array(rot_moms)


def StatewiseRealMomentsToComplexMoments_SLOW(state_moments, order):
    cmom = list()
    for i_state in range(state_moments.shape[0]):
        cmom.append(RealMomentsToComplexMoments(state_moments[i_state], order))
    return np.array(cmom)


def StatewiseComplexMomentsToRealMoments_SLOW(state_moments, order):
    rmom = list()
    for i_state in range(state_moments.shape[0]):
        rmom.append(ComplexMomentsToRealMoments(state_moments[i_state], order))
    return np.array(rmom)


# performant versions not finally debugged. ...

def StatewiseRealMomentsToComplexMoments(real_moments, order): # use with care - debugged only for l=1
    n_states, n_moments = real_moments.shape
    u_sh = SphericalHarmonics.RealSphericalHarmonicsTransformation(order)
    return np.dot(np.conj(u_sh.T)[np.newaxis,:,:], real_moments[:,:,np.newaxis]).T.reshape(n_states, n_moments)


def StatewiseComplexMomentsToRealMoments(comp_moments, order): # use with care - debugged only for l=1
    n_states, n_moments = comp_moments.shape
    u_sh = SphericalHarmonics.RealSphericalHarmonicsTransformation(order)
    return np.real(np.dot(u_sh[np.newaxis,:,:], comp_moments[:,:,np.newaxis]).T.reshape(n_states, n_moments))

def KabschAlgorithm(P, Q, output=False):
    """ The Kabsch algorithm
http://en.wikipedia.org/wiki/Kabsch_algorithm
The algorithm starts with two sets of paired points P and Q.
P and Q should already be centered on top of each other.
Each vector set is represented as an NxD matrix, where D is the
the dimension of the space.
The algorithm works in three steps:
- a translation of P and Q
- the computation of a covariance matrix C
- computation of the optimal rotation matrix U
The optimal rotation matrix U is then used to
rotate P unto Q so the RMSD can be caculated
from a straight forward fashion.
"""
    # Computation of the covariance matrix
    C = np.dot(np.transpose(P), Q)
    # Computation of the optimal rotation matrix
    # This can be done using singular value decomposition (SVD)
    # Getting the sign of the det(V)*(W) to decide
    # whether we need to correct our rotation matrix to ensure a
    # right-handed coordinate system.
    # And finally calculating the optimal rotation matrix U
    # see http://en.wikipedia.org/wiki/Kabsch_algorithm
    V, S, W = np.linalg.svd(C)
    d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0
    if d:
        S[-1] = -S[-1]
        V[:,-1] = -V[:,-1]
    # Create Rotation matrix U
    U = np.dot(V, W)
    return U


def AdjustGeometry(traj_pos_aa, ref_geo_aa, weights, traj_vel_au=None):
    """Expects (n_frames,n_atoms,3) (n_atoms,3) (n_atoms) (n_frames,n_atoms,3) (n_frames,3,3)"""
    n_frames = traj_pos_aa.shape[0]
    # remove centers of mass
    ref_geo_aa     -= np.sum(ref_geo_aa*np.array(weights)[:,np.newaxis], axis=0)/np.sum(weights)
    traj_pos_com_aa = np.sum(traj_pos_aa*np.array(weights)[np.newaxis,:,np.newaxis], axis=1)/np.sum(weights)
    traj_pos_aa    -= traj_pos_com_aa[:,np.newaxis,:]
    # massweight coordinates
    mw_ref_aa = np.sqrt(weights)[:,np.newaxis]*ref_geo_aa
    mw_pos_aa = np.sqrt(weights)[np.newaxis,:,np.newaxis]*traj_pos_aa

    R = np.zeros((n_frames, 3, 3))
    for i_frame in range(n_frames):
        R[i_frame] = KabschAlgorithm(mw_pos_aa[i_frame], mw_ref_aa)

    rot_traj_pos_aa = np.sum(traj_pos_aa[:,:,:,np.newaxis]*R[:,np.newaxis,:,:], axis=2) + traj_pos_com_aa[:,np.newaxis,:]
    if traj_vel_au != None:
        rot_traj_vel_au = np.sum(traj_vel_au[:,:,:,np.newaxis]*R[:,np.newaxis,:,:], axis=2)
        return rot_traj_pos_aa, R, traj_pos_com_aa, rot_traj_vel_au
    else:
        return rot_traj_pos_aa, R, traj_pos_com_aa
