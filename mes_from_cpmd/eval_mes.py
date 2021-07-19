#import fortran_io as fio
import mes_from_cpmd.ext_fortran.fortran_io as fio
#from chitools import rpa_tools
from mes_from_cpmd.toolbox import CubeFileTools
from mes_from_cpmd.toolbox import transformations

import numpy as np
from mes_from_cpmd.toolbox import lib_dme as lime
#from reduced_eigen import calc_dens as calc
import ipdb
from mes_from_cpmd.misc import git_control
import mes_from_cpmd.toolbox.cube as cube 

import numpy as np
import os
import sys
import argparse
import subprocess


def PrintMatrix(mat):
    for i in range(mat.shape[0]):
        print(('%12.6F'*mat.shape[1])%tuple(mat[i]))



def main():
        print("command line:")
        print(sys.argv)
        git_control.get_git_version()

        parser = argparse.ArgumentParser()
        parser.add_argument("path_to_DENSITY", help="path to DENSITY.cube file of target molecule for MES calculation)", default = "../DENSITY.cube")
        #parser.add_argument("path_to_input_for_cpmd", help="path to files required for density calculation of the target molecule due to the basis function of the  potential (for  CPMD: prepare_wfo, calc_wfo, wfo.job, wfo.inp, col_dens, new_run_cpmd2cube )")
        #parser.add_argument("path_to_density_responses", help="path to  the overall distorted densities")
        #parser.add_argument("path_to_difference_densities", help="path to  the difference densities (= ntilde states)")
        parser.add_argument("path_to_dmes", help="path to  the direct moment expanded states")
        parser.add_argument("path_to_mes", help="path to  the reference moment expanded states")
        parser.add_argument("path_to_potentials", help="path to  the basis function of the  potential")
        parser.add_argument("n", help="number of basis funtion of perturbing potentials used for calculation",type=int)
        args = parser.parse_args()
        path_pot = args.path_to_potentials
        #path_diff = args.path_to_difference_densities
        path_dmes = args.path_to_dmes
        path_mes = args.path_to_mes
        n_states = args.n
        fn_cube = args.path_to_DENSITY

        cell_data_mom = CubeFileTools.LoadCellData(fn_cube)

        #r_au =CubeFileTools.CalcGridPositions(cell_data_mom['cell_au'], cell_data_mom['mesh'], origin)
        d3r_au = cell_data_mom['d3r_au']
        #dummy, n_x, n_y, n_z = r_au.shape
        #volume_au = d3r_au*n_x*n_y*n_z


        #fn2 = '/home/dressler/post_doc/chi/gneral_chi_testing/gen_pot/pot_meth_v1/cartesian-functions-%05d.wan'
        fn2 = path_pot + '/cartesian-functions-%05d.wan'
        states_compare = lime.load_states(fn_cube, fn2, n_states, pure = 'True', pert = 'True')

        #fn3 = '/home/dressler/post_doc/chi/gneral_chi_testing/gen_pot/pot_meth_v1/density/col_density/out/normal-states-eq-mes-%05d.wan'
        fn3 = path_dmes + '/dme-%05d.wan'
        states_dme = lime.load_states(fn_cube, fn3, n_states, pure =  True, pert = False)

        fn4 = path_mes + '/mes-%05d.wan'
        states_mes = lime.load_states(fn_cube, fn4, n_states, pure =  True, pert = False)
       


        print("moment matrix for reference moment expanded states ")
        #mom_tilde_self = lime.create_overlap_mat(states_tilde*d3r_au *1000, states_compare/0.001)
        #avg_mom_tilde = (mom_tilde_self[1:n_states, 1:n_states] + mom_tilde_self[1:n_states, 1:n_states].T)/2
        
        #ref_mat = lime.create_overlap_mat(states_mes, states_compare[1:]/0.001)
#
        ref_mat = lime.create_overlap_mat(states_mes, states_compare[1:])
        PrintMatrix(ref_mat)
        #print(ref_mat)
        #avg_mom_tilde = (mom_tilde_self[1:n_states, 1:n_states] + mom_tilde_self[1:n_states, 1:n_states].T)/2
        #PrintMatrix(avg_mom_tilde)
        print("moment matrix for direct moment expanded states ")
        #mom_tilde_self = lime.create_overlap_mat(states_tilde*d3r_au *1000, states_compare/0.001) 
        #avg_mom_tilde = (mom_tilde_self[1:n_states, 1:n_states] + mom_tilde_self[1:n_states, 1:n_states].T)/2

        direct_mat = lime.create_overlap_mat(states_dme, states_compare[1:])
        PrintMatrix(direct_mat)
        #achtung falls die funktion cacl_moments_real_batch benutzt wird muss noch mit *np.sqrt(d3r_au) multipliziert werden
#        print()
#        PrintMatrix(direct_mat*np.sqrt(d3r_au))















