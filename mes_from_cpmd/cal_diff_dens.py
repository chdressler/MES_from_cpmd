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

def main():
        print("command line:")
        print(sys.argv)
        git_control.get_git_version()

        parser = argparse.ArgumentParser()
        parser.add_argument("path_to_DENSITY", help="path to DENSITY.cube file of target molecule for MES calculation)", default = "../DENSITY.cube")
        #parser.add_argument("path_to_input_for_cpmd", help="path to files required for density calculation of the target molecule due to the basis function of the  potential (for  CPMD: prepare_wfo, calc_wfo, wfo.job, wfo.inp, col_dens, new_run_cpmd2cube )")
        parser.add_argument("path_to_density_responses", help="path to  the overall distorted densities")
        #parser.add_argument("path_to_potentials", help="path to  the basis function of the  potential")
        parser.add_argument("n", help="number of basis funtion of perturbing potentials used for calculation",type=int)
        args = parser.parse_args()
        #path_pot = args.path_to_potentials
        cube_ref_eq = args.path_to_DENSITY
        path_resp = args.path_to_density_responses
        #path_inp = args.path_to_input_for_cpmd
        

        #cube_ref_eq =  '/net/shared/dressler/chi/paper_no4_theorem/data/new_chis/methanal/1/chi1/DENSITY.cube'
        #cube_ref_eq =  args.path_to_DENSITY 
        #'/DENSITY.cube'
        ref_eq = CubeFileTools.LoadCellData(cube_ref_eq) 
                
        
        #j = 12
        j = args.n
        states_ntilde = []
        #print('calc chol')
        #for i in range(j):
        for k in range(j):
            i = k + 1
            cube_in = path_resp  +"/"+ str(i) + "/DENSITY_fullmesh.cube" 
            #cube_in = '/home/dressler/post_doc/chi/gneral_chi_testing/pre_publish/data_pot/dens/col_density/DENSITY'+ str(i+1) + ".cube"
            tmp = CubeFileTools.LoadCellData(cube_in)
            states_ntilde.append(tmp['data'] - ref_eq['data'])
        states_ntilde = np.array(states_ntilde)
        
        bn_states_out = 'difference_density-%05d.wan'       
        bn_cube_out  =  'difference_density-%05d.cube'

        #for i_state in range(n_states):
        n_x, n_y, n_z = ref_eq['mesh']
        for i_state in range(j):
           state_data = np.asfortranarray(states_ntilde[i_state]).astype(np.float64)
           #fio.fortran_write_unformatted(bn_states_out%(key,(i_state+1)), state_data, n_x, n_y, n_z)
           #cube.WriteCubeFile(bn_cube_out%(key,(i_state+1)),'','', cell_data['numbers'], cell_data['coords_au'], cell_data['cell_au'], state_data, origin=cell_data['origin_au'])
           fio.fortran_write_unformatted(bn_states_out%(i_state+1), state_data, n_x, n_y, n_z)
           #cube.WriteCubeFile(bn_cube_out%(i_state+1),'','', cell_data['numbers'], cell_data['coords_au'], cell_data['cell_au'], state_data, origin=cell_data['origin_au'])
           cube.WriteCubeFile(bn_cube_out%(i_state+1),'','', ref_eq['numbers'], ref_eq['coords_au'], ref_eq['cell_au'], state_data, origin=ref_eq['origin_au'])















