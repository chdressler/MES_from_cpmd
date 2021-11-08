import mes_from_cpmd.ext_fortran.fortran_io as fio
from mes_from_cpmd.toolbox import CubeFileTools
from mes_from_cpmd.toolbox import transformations
import numpy as np
from mes_from_cpmd.toolbox import lib_dme as lime
import ipdb
from mes_from_cpmd.misc import git_control
import mes_from_cpmd.toolbox.cube as cube 
from scipy.signal import convolve

import numpy as np
import os
import sys
import argparse
import subprocess
import ipdb
from scipy.ndimage import gaussian_filter

def PrintMatrix(mat):
    for i in range(mat.shape[0]):
        print(('%12.6F'*mat.shape[1])%tuple(mat[i]))


def gauss_kern_3d(size):
    size = int(size)
    #if not sizey:
    #    sizey = size
    #else:
    #    sizey = int(sizey)
    #from mes_from_cpmd.toolbox import CubeFileTools

    cell_au = np.array([[1,0,0],[0,1,0],[0,0,1]])
    mesh = np.array([2*size+1,2*size+1,2*size+1])
    origin = np.array([size,size,size])
    r_au =CubeFileTools.CalcGridPositions(cell_au, mesh, origin)
    
    size = 0.1 *size
    #x, y = mgrid[-size:size+1, -size:size+1]
    #g = exp(-(x**2/float(size)+y**2/float(sizey)))
    g = np.exp(-(r_au[0,:,:,:]**2/float(size)+r_au[1,:,:,:]**2/float(size)+r_au[2,:,:,:]**2/float(size)))
    return g / g.sum()



def main():
        print("command line:")
        print(sys.argv)
        git_control.get_git_version()

        parser = argparse.ArgumentParser('the most important part of this module. calculates the moment expanded states and performs the actual linear algebra operations. requires the difference densities (n_resp) and the basis functions of the perturbing potential.')
        parser.add_argument("path_to_DENSITY", help="path to DENSITY.cube file of target molecule for MES calculation)", default = "../DENSITY.cube")
        parser.add_argument("path_to_difference_densities", help="path to  the difference densities (= ntilde states)")
        parser.add_argument("path_to_potentials", help="path to  the basis function of the  potential")
        parser.add_argument("n", help="number of basis funtion of perturbing potentials used for calculation",type=int)
        parser.add_argument('--smooth', type = int,  help='Size of Gaussian for convolution')
        args = parser.parse_args()
        path_pot = args.path_to_potentials
        path_diff = args.path_to_difference_densities
        n_states = args.n
        fn_cube = args.path_to_DENSITY
        smooth = args.smooth



def perform_direct_me(fn_cube, path_diff, path_pot, n_states, smooth):
        cell_data_mom = CubeFileTools.LoadCellData(fn_cube)
        d3r_au = cell_data_mom['d3r_au']
        fn2 = path_pot + '/cartesian-functions-%05d.wan'
        states_compare = lime.load_states(fn_cube, fn2, n_states, pure = 'True', pert = 'True')
        fn3 = path_diff + '/difference_density-%05d.wan'
        #load raw density responses (=states_tilde) from cpmd
        states_tilde = lime.load_states(fn_cube, fn3, n_states, pure =  True, pert = False)
        

        
        print("overlaps of difference densities (=ntilde states) with basis functions of perturbing potential  ")
        #mom_tilde_self = lime.create_overlap_mat(states_tilde*d3r_au *1000, states_compare/0.001)
        
        if smooth:  
            #ipdb.set_trace()
            gauss_kernel = gauss_kern_3d(smooth)
            for i in range(states_tilde.shape[0]):
            #states_tilde[i] = convolve(states_tilde[i] ,gauss_kernel , mode='valid')
                #states_tilde[i] = convolve(states_tilde[i] ,gauss_kernel , mode='same')
                states_tilde[i] = gaussian_filter(states_tilde[i] , sigma=smooth)
            #states_tilde[i] = convolve(states_tilde[i] ,gauss_kernel )
            #data = gaussian_filter(data, sigma=10) 
        print("again overlaps of difference densities (=ntilde states) with basis functions of perturbing potential  ")
        mom_tilde_self = lime.create_overlap_mat(states_tilde*d3r_au *1000, states_compare/0.001)
        print("native overlap matrix")
        PrintMatrix(mom_tilde_self[1:n_states, 1:n_states])
        avg_mom_tilde = (mom_tilde_self[1:n_states, 1:n_states] + mom_tilde_self[1:n_states, 1:n_states].T)/2
        #avg_mom_tilde = (mom_tilde_self[1:n_states, 1:n_states].T + mom_tilde_self[1:n_states, 1:n_states].T)/2
        #avg_mom_tilde = (mom_tilde_self[1:n_states, 1:n_states] + mom_tilde_self[1:n_states, 1:n_states])/2
        print("symmetrized overlap matrix")
        PrintMatrix(avg_mom_tilde)

        print("direct moment expansion: dme")
        print("cholseky decomposition of tilde moment matrix")
        print(" M = R * R^T")
        print("R")
        chol_tmp = np.linalg.cholesky(avg_mom_tilde*-1)
        PrintMatrix(chol_tmp.T)
        chol1 = chol_tmp.T
        inv_chol_right = np.linalg.inv(chol1)
        mom_states_new = np.zeros(states_tilde.shape)
        print(n_states -1)
        for i in range(n_states-1):
            for j in range(n_states-1):
                #mom_states_new[i] += inv_chol_right[j,i] * states_tilde[j+1]*1000*np.sqrt(d3r_au)
                mom_states_new[i] += inv_chol_right[j,i] * states_tilde[j+1]*1000*np.sqrt(d3r_au) *  -1000 /8
                #mom_states_new[i] += inv_chol_right[j,i] * states_tilde[j+1]*d3r_au  

        bn_states_out = 'dme-%05d.wan'       
        bn_cube_out  =  'dme-%05d.cube'

        n_x, n_y, n_z = cell_data_mom['mesh']
        for i_state in range(n_states-1):
           state_data = np.asfortranarray(mom_states_new[i_state]).astype(np.float64)
           fio.fortran_write_unformatted(bn_states_out%(i_state+1), state_data, n_x, n_y, n_z)
           #cube.WriteCubeFile(bn_cube_out%(i_state+1),'','', cell_data_mom['numbers'], cell_data_mom['coords_au'], cell_data_mom['cell_au'], state_data, origin=cell_data_mom['origin_au'])


















        #ref_eq = CubeFileTools.LoadCellData(cube_ref_eq) 
        #        
        #
        #j = args.n
        #states_ntilde = []
        #for i in range(j):
        #    cube_in = path_resp  + str(i) + "/DENSITY_fullmesh.cube" 
        #    tmp = CubeFileTools.LoadCellData(cube_in)
        #    states_ntilde.append(tmp['data'] - ref_eq['data'])
        #states_ntilde = np.array(states_ntilde)
        #
        #bn_states_out = 'difference_density-%05d.wan'       
        #bn_cube_out  =  'difference_density-%05d.cube'

        #n_x, n_y, n_z = ref_eq['mesh']
        #for i_state in range(j):
        #   state_data = np.asfortranarray(states_ntilde[i_state]).astype(np.float64)
        #   fio.fortran_write_unformatted(bn_states_out%(i_state), state_data, n_x, n_y, n_z)
        #   cube.WriteCubeFile(bn_cube_out%(i_state),'','', ref_eq['numbers'], ref_eq['coords_au'], ref_eq['cell_au'], state_data, origin=ref_eq['origin_au'])















