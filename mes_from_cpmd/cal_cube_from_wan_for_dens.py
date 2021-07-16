#import fortran_io as fio
#from chitools import rpa_tools
#from mes_from_cpmd.toolbox import CubeFileTools
#from mes_from_cpmd.toolbox import transformations

import numpy as np
#from mes_from_cpmd.toolbox import lib_dme as lime
#from reduced_eigen import calc_dens as calc
import ipdb
from mes_from_cpmd.misc import git_control


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
        parser.add_argument("path_to_input_for_cpmd", help="path to files required for density calculation of the target molecule due to the basis function of the  potential (for  CPMD: prepare_wfo, calc_wfo, wfo.job, wfo.inp, col_dens, new_run_cpmd2cube )")
        #parser.add_argument("path_to_potentials", help="path to  the basis function of the  potential")
        #parser.add_argument("path_to_potentials", help="path to  the basis function of the  potential")
        parser.add_argument("n", help="number of basis funtion of perturbing potentials used for calculation",type=int)
        args = parser.parse_args()
        #path_pot = args.path_to_potentials
        path_inp = args.path_to_input_for_cpmd
        wd = os.getcwd()
        for i in range(args.n):
            #list_files = subprocess.run(["mkdir", str(i)]) 
            os.chdir("./"+ str(i))
            #list_files = subprocess.run(["cp", path_inp +"/wfo.inp", "." ])
            list_files = subprocess.run(["cp", path_inp +"/run_cpmd2cube", "." ])
            #list_files = subprocess.run(["cp", path_inp +"/wfo.job", "." ])
            #list_files = subprocess.run(["cp", path_pot +"/cartesian-functions-0000"+ str(i+1) +".wan", "extpot.unfo.grid" ])
            #os.chdir(wd)            
       # for i in range(args.n):
       #     os.chdir("./"+ str(i))
        #    print("caluate density for potential number " + str(i))
            #list_files = subprocess.run(["bash", "wfo.job"])
            list_files = subprocess.run(["bash", "run_cpmd2cube"])
            #list_files = subprocess.run(["cpmd2cube_old.x", "-fullmesh", "-o", "DENSITY_fullmesh", "-dens", "DENSITY", " >", "out_fullmesh.cpmd2cube_old"])
            os.chdir(wd)
            #cpmd2cube_old.x -fullmesh -o DENSITY_fullmesh -dens DENSITY > out_fullmesh.cpmd2cube_old
