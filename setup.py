from setuptools import setup
import subprocess


def get_commit_hash():
        command = "git log -n 1 --format=%H%n%s%n%ad"
        try:
            commit_hash, commit_message, commit_date = subprocess.check_output(command.split()).strip().split(b"\n")
        except subprocess.CalledProcessError:
            print("Command '{}' could not be executed successfully.".format(command), file=sys.stderr)
        return commit_hash.decode(), commit_message.decode(), commit_date.decode()




setup(name='mes_from_cpmd',
     version='19.1',
     #description='Markov models for lithium transfer',
     #author='Christian Dreßler',
                               entry_points={
                                    'console_scripts': [
                                        #"jump_mat_ana = ana_tools.jump_mat_analyzer:main",
                                        #"create_jump_mat_li = markov.neighbor:main",
                                        #"jumps_from_grid = markov.jumps_from_grid:main",
                                        #"jump_trajek_from_grid =  markov.jump_trajek_from_grid:main", 
                                        #"msd_from_markov = markov.msd:main",
                                        #"get_jump_arrows = markov.get_jump_arrows:main",
                                        #"recover_old_mat_exe = markov.neighbor:recover_old_mat",
                                        #"print_eigenfunctions = markov.print_eigenstates:main",
                                        #"create_nn_dist_list = ana_tools.create_nn_dist_list:main",
                                        #"ana_charge_ar = ascripts.ana_charge_ar:main",
                                        #"ana_charge_ar = markov.ana_charge_ar:main",
                                        #"lattice_ana_charge_ar = scripts.ana_lattice_charge:main",
                                        #"wrap_li_to_box = ascripts.wrap_li_to_box:main",
                                        #"reduce_trajec = ascripts.reduce_trajec:main",
                                        #"reduce_trajec = bscripts.reduce_trajec:main",
                                        #"prepare_trajec = ascripts.prepare_trajec:main"]
                                        "first_try = mes_from_cpmd.1:main"]
                                                               },
                                    #zip_safe=False)
include_package_data=True,
      zip_safe=False,
      #ext_modules=ext_modules,
      #cmdclass={'build_ext': build_ext}
                                    
)

# Write hash, message and date of current commit to file
with open("commit_stuff/version_hash.py", "w") as f:
        commit_hash, commit_message, commit_date = get_commit_hash()
        print("commit_hash = \"{}\"".format(commit_hash), file=f)
        print("commit_message = \"{}\"".format(commit_message), file=f)
        print("commit_date = \"{}\"".format(commit_date), file=f)
