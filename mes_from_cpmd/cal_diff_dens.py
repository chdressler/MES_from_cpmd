import mes_from_cpmd.ext_fortran.fortran_io as fio
from mes_from_cpmd.toolbox import CubeFileTools
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

        parser = argparse.ArgumentParser('calculate the difference (n_resp = n_pot - n_eq) of the perturbed densities (n_pot) and the equilibrium densitiy  (n_eq) of the target molecule. requires the path to the perturbed densities (n_pot) and the equilibrium densitiy  (n_eq)')
        parser.add_argument("path_to_DENSITY", help="path to DENSITY.cube file of target molecule for MES calculation)", default = "../DENSITY.cube")
        parser.add_argument("path_to_density_responses", help="path to  the overall distorted densities")
        parser.add_argument("n", help="number of basis funtion of perturbing potentials used for calculation",type=int)
        args = parser.parse_args()
        cube_ref_eq = args.path_to_DENSITY
        path_resp = args.path_to_density_responses
        

        ref_eq = CubeFileTools.LoadCellData(cube_ref_eq) 
                
        
        j = args.n
        states_ntilde = []
        for k in range(j):
            i = k + 1
            cube_in = path_resp  +"/"+ str(i) + "/DENSITY_fullmesh.cube" 
            tmp = CubeFileTools.LoadCellData(cube_in)
            states_ntilde.append(tmp['data'] - ref_eq['data'])
        states_ntilde = np.array(states_ntilde)
        
        bn_states_out = 'difference_density-%05d.wan'       
        bn_cube_out  =  'difference_density-%05d.cube'

        n_x, n_y, n_z = ref_eq['mesh']
        for i_state in range(j):
           state_data = np.asfortranarray(states_ntilde[i_state]).astype(np.float64)
           fio.fortran_write_unformatted(bn_states_out%(i_state+1), state_data, n_x, n_y, n_z)
           cube.WriteCubeFile(bn_cube_out%(i_state+1),'','', ref_eq['numbers'], ref_eq['coords_au'], ref_eq['cell_au'], state_data, origin=ref_eq['origin_au'])















