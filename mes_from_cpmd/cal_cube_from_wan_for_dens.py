import numpy as np
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
        parser.add_argument("n", help="number of basis funtion of perturbing potentials used for calculation",type=int)
        args = parser.parse_args()
        path_inp = args.path_to_input_for_cpmd
        wd = os.getcwd()
        for j in range(args.n):
            i = j+1
            os.chdir("./"+ str(i))
            list_files = subprocess.run(["cp", path_inp +"/run_cpmd2cube", "." ])
            list_files = subprocess.run(["bash", "run_cpmd2cube"])
            os.chdir(wd)
