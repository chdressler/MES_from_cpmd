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

        parser = argparse.ArgumentParser('calculate densities for the applied basis functions of the perturbing potential  (n_pot), requires path to the  basis functions of perturbing potential')
        parser.add_argument("path_to_input_for_cpmd", help="path to files required for density calculation of the target molecule due to the basis function of the  potential (for  CPMD: prepare_wfo, calc_wfo, wfo.job, wfo.inp, col_dens, new_run_cpmd2cube )")
        parser.add_argument("path_to_potentials", help="path to  the basis function of the  potential")
        parser.add_argument("n", help="number of basis funtion of perturbing potentials used for calculation",type=int)
        parser.add_argument("--noslurm",  action="store_true", help="Calc jobs on current node.")
        args = parser.parse_args()
        path_pot = args.path_to_potentials
        path_inp = args.path_to_input_for_cpmd
        wd = os.getcwd()
        for j in range(args.n):
            i = j+1
            list_files = subprocess.run(["mkdir", str(i)]) 
            os.chdir("./"+ str(i))
            list_files = subprocess.run(["cp", path_inp +"/wfo.inp", "." ])
            #list_files = subprocess.run(["cp", path_inp +"/wfo.job", "." ])
            list_files = subprocess.run(["cp", path_inp +"/calc.job", "." ])
            #path_pot_final = 
            list_files = subprocess.run(["cp", path_pot +"/cartesian-functions-%05d.wan"%(i), "extpot.unfo.grid" ])
            os.chdir(wd)            
        for j in range(args.n):
            i = j+1
            os.chdir("./"+ str(i))
            print("caluate density for potential number " + str(i))
            if args.noslurm:
                list_files = subprocess.run(["bash", "wfo.job"])
            else:    
                list_files = subprocess.run(["sbatch", "calc.job"])
            os.chdir(wd)
