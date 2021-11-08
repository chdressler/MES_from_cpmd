import ipdb
from mes_from_cpmd.misc import git_control
import numpy as np
import os
import sys
import argparse
import subprocess
from  mes_from_cpmd.pot.gen_pot import create_potentials
from mes_from_cpmd.cal_dens_for_pot import calc_dens_resp
from mes_from_cpmd.cal_diff_dens import calc_diff_densities
from mes_from_cpmd.direct_moment_expansion import perform_direct_me
def main():
        print("command line:")
        print(sys.argv)
        git_control.get_git_version()

        parser = argparse.ArgumentParser(' calculates the first N moment expanded states for the first N monomial basis functions of the perturbing potential, requires only input file for CPMD single point calculation and path to slurm job')
        parser.add_argument("path_to_input_for_cpmd", help="path to cpmd input file  for density calculation of the target molecule. (recommended name: wfo.inp)")
        #parser.add_argument("path_to_potentials", help="path to  the basis function of the  potential")
        parser.add_argument("path_to_slurm", help="path to  slurm job file (including paths to cpmd and cpmd_to_cube executable,recommended name: wfo.job)")
        parser.add_argument("n", help="number of basis funtion of perturbing potentials used for calculation of the moment expanded states",type=int)
        parser.add_argument("--noslurm",  action="store_true", help="Calc jobs on current node.")
        args = parser.parse_args()
        path_input = args.path_to_input_for_cpmd
        path_slurm = args.path_to_slurm
        wd = os.getcwd()
        print(wd)
        if not os.path.isdir("single_point"):
            list_files = subprocess.run(["mkdir", "single_point"]) 
            os.chdir("./"+ "single_point")
            list_files = subprocess.run(["cp", path_input, "." ])
            list_files = subprocess.run(["cp", path_slurm , "." ])
            if args.noslurm: 
                list_files = subprocess.run(["bash", "wfo.job"])
                #list_files = subprocess.run(["bash", "wfo.job"])
            else:    
                i#list_files = subprocess.run(["sbatch", "calc.job"]) 
                list_files = subprocess.run(["sbatch", "--wait" , "calc.job"]) 
            os.chdir(wd)
        #input("Press Enter to continue (after slurm jobs are finished)...")

        if not os.path.isdir("potentials"):
            list_files = subprocess.run(["mkdir", "potentials"])
            os.chdir("./"+ "potentials")
            create_potentials(wd + "/single_point/DENSITY_fullmesh.cube", args.n, 0.1)
            os.chdir(wd)

        os.chdir(wd)
        file1 = open(path_input, 'r')
        Lines = file1.readlines()
        index1 =Lines.index('&CPMD\n')
        Lines.insert(index1+1, ' EXTERNAL POTENTIAL\n')
        os.chdir("./"+ "potentials")
        file2 = open('wfo.inp', 'w')
        file2.writelines(Lines)
        file2.close()        
        os.chdir(wd)
        if not os.path.isdir("densities_ext_pot"):
            list_files = subprocess.run(["mkdir", "densities_ext_pot"])
            os.chdir("./"+ "densities_ext_pot")     
            calc_dens_resp(wd + "/potentials/wfo.inp", path_slurm, wd + "/potentials", args.n, args.noslurm)
        input("Press Enter to continue (after slurm jobs are finished)...")

        os.chdir(wd)
        if not os.path.isdir("difference_densities"):
            list_files = subprocess.run(["mkdir", "difference_densities"])
            os.chdir("./"+ "difference_densities")

            calc_diff_densities(wd + "/densities_ext_pot", wd + "/single_point/DENSITY_fullmesh.cube", args.n)  
        
                 
        os.chdir(wd)
        if not os.path.isdir("moment_expanded_states"):
            list_files = subprocess.run(["mkdir", "moment_expanded_states"])
            os.chdir("./"+ "moment_expanded_states")
            perform_direct_me(wd + "/single_point/DENSITY_fullmesh.cube", wd + "/difference_densities", wd + "/potentials",args.n, False )
            
#        for j in range(args.n):
#            i = j+1
#            list_files = subprocess.run(["mkdir", str(i)]) 
#            os.chdir("./"+ str(i))
#            list_files = subprocess.run(["cp", path_inp +"/wfo.inp", "." ])
#            #list_files = subprocess.run(["cp", path_inp +"/wfo.job", "." ])
#            list_files = subprocess.run(["cp", path_inp +"/calc.job", "." ])
#            #path_pot_final = 
#            list_files = subprocess.run(["cp", path_pot +"/cartesian-functions-%05d.wan"%(i), "extpot.unfo.grid" ])
#            os.chdir(wd)            
#        for j in range(args.n):
#            i = j+1
#            os.chdir("./"+ str(i))
#            print("caluate density for potential number " + str(i))
#            if args.noslurm:
#                list_files = subprocess.run(["bash", "wfo.job"])
#            else:    
#                list_files = subprocess.run(["sbatch", "calc.job"])
#            os.chdir(wd)
