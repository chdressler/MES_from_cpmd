#import fortran_io as fio
#from chitools import rpa_tools
from mes_from_cpmd.toolbox import CubeFileTools
from mes_from_cpmd.toolbox import transformations

import numpy as np
from mes_from_cpmd.toolbox import lib_dme as lime
#from reduced_eigen import calc_dens as calc
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

        parser = argparse.ArgumentParser()
        parser.add_argument("--path_to_DENSITY", help="path to DENSITY.cube file of target molecule for MES calculation)", default = "../")
        args = parser.parse_args()
        
        path_dens = args.path_to_DENSITY
        fn_cube5 = path_dens + 'DENSITY.cube'
        cell_data5 = CubeFileTools.LoadCellData(fn_cube5)
        n_states  = 3
        exp_order = 1
        
        print("atoms", cell_data5['numbers'])
        
        moa = np.copy(np.array(cell_data5['numbers']))
        moa[moa == 6] = 12
        moa[moa == 8] = 16
        print("masses", list(moa))
        origin = transformations.CentersOfMass(cell_data5['coords_au'], moa)
        ipdb.set_trace()

        #create lattice coordinates for example r_au[3,108,108,108]
        r_au =CubeFileTools.CalcGridPositions(cell_data5['cell_au'], cell_data5['mesh'], origin)
        
        mompol   = dict()
        moments  = dict()
        
        
        
        for i_order in range(1, exp_order+1):
                    #evaluate polynoms on grid
                    mompol['%02d'%i_order] = CubeFileTools.CartesianMoments(r_au, order=i_order)
        mom_pol_ar = np.zeros((n_states + 1 , 108, 108 , 108))
        mom_pol_ar[0] = np.ones((108, 108 , 108))
        
        
        
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
        
        
        
        lime.print_states_dict(mom_pol_ar*0.001, cell_data5, n_states + 1, 'cartesian-functions-%05d', pure = True, pert = True)

