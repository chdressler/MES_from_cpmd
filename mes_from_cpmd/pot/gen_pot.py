from mes_from_cpmd.toolbox import CubeFileTools
from mes_from_cpmd.toolbox import transformations
from mes_from_cpmd.toolbox import lib_dme as lime
import ipdb
from mes_from_cpmd.misc import git_control
import numpy as np
import os
import sys
import argparse


def main():
        print("command line:")
        print(sys.argv)
        git_control.get_git_version()

        parser = argparse.ArgumentParser('calculate (here monomial) basis functions of perturbing potential, requires  the DENSITY.cube file from  cpmd single point calcaultion of the target molecule for the extraction of the  grid parameters.')
        parser.add_argument("--path_to_DENSITY", help="path to DENSITY.cube file of target molecule for MES calculation)", default = "../DENSITY.cube")
        parser.add_argument("n", help="number of monomial basis funtion for perturbing potentials, supported number of monomials: 4,10,20,35 or 56.",type=int)
        parser.add_argument('--fac', type = float,  help='Factor for rescaling of potential')
        args = parser.parse_args()
        
        path_dens = args.path_to_DENSITY
        fn_cube_ref_dens = path_dens 
        cell_data_ref_dens = CubeFileTools.LoadCellData(fn_cube_ref_dens)
        exp_order_dict = {}
        exp_order_dict[4] =  1
        exp_order_dict[10] =  2
        exp_order_dict[20] =  3
        exp_order_dict[35] =  4
        exp_order_dict[56] =  5
        n_states = args.n
        exp_order = exp_order_dict[n_states]
        n_states -= 1
        print(n_states, exp_order) 
        print("atoms", cell_data_ref_dens['numbers'])
        
        moa = np.copy(np.array(cell_data_ref_dens['numbers']))
        moa[moa == 6] = 12
        moa[moa == 8] = 16
        print("masses", list(moa))
        origin = transformations.CentersOfMass(cell_data_ref_dens['coords_au'], moa)
        
        #create lattice coordinates for example r_au[3,108,108,108]
        r_au =CubeFileTools.CalcGridPositions(cell_data_ref_dens['cell_au'], cell_data_ref_dens['mesh'], origin)
        
        mompol   = dict()
        moments  = dict()
        
        
        
        for i_order in range(1, exp_order+1):
                    #evaluate polynoms on grid
                    mompol['%02d'%i_order] = CubeFileTools.CartesianMoments(r_au, order=i_order)
        #mom_pol_ar = np.zeros((n_states + 1 , 108, 108 , 108))
        mom_pol_ar = np.zeros((n_states + 1 , cell_data_ref_dens['mesh'][0], cell_data_ref_dens['mesh'][1] , cell_data_ref_dens['mesh'][2]))
        #mom_pol_ar[0] = np.ones((108, 108 , 108))
        mom_pol_ar[0] = np.ones((cell_data_ref_dens['mesh'][0], cell_data_ref_dens['mesh'][1] , cell_data_ref_dens['mesh'][2]))
        
        
        
        #convert dictionary to array for example mom_pol_ar[55, 3, 108, 108,108]
        for j in range(n_states):
                k = j +1
                if (j < 3) and (j > -1):
                    mom_pol_ar[k] = mompol['01'][j]
                elif (j < 9) and (j > 2):
                    mom_pol_ar[k] = mompol['02'][j-3]
                elif (j < 19) and (j > 8):
                    mom_pol_ar[k] = mompol['03'][j-9]
                elif (j < 34) and (j > 18):
                    mom_pol_ar[k] = mompol['04'][j-19]
                elif (j < 55) and (j > 33):
                    mom_pol_ar[k] = mompol['05'][j-34]
                else:
                    print('index error')
        
        
        if args.fac:
            lime.print_states_dict(mom_pol_ar*0.001*args.fac, cell_data_ref_dens, n_states + 1, 'cartesian-functions-%05d', pure = True, pert = True)
        else:
            lime.print_states_dict(mom_pol_ar*0.001, cell_data_ref_dens, n_states + 1, 'cartesian-functions-%05d', pure = True, pert = True)

