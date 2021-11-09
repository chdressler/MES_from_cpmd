import os
from setuptools import setup
import subprocess


from numpy.distutils.core import setup, Extension

def get_commit_hash():
        command = "git log -n 1 --format=%H%n%s%n%ad"
        try:
            commit_hash, commit_message, commit_date = subprocess.check_output(command.split()).strip().split(b"\n")
        except subprocess.CalledProcessError:
            print("Command '{}' could not be executed successfully.".format(command), file=sys.stderr)
        return commit_hash.decode(), commit_message.decode(), commit_date.decode()


def find_packages():
    """Find all packages (i.e. folders with an __init__.py file inside)"""
    packages = []
    for root, _, files in os.walk("mes_from_cpmd"):
        for f in files:
            if f == "__init__.py":
                packages.append(root)
    return packages


def readme():
    with open('README.md', 'r') as f:
        return f.read()



print(find_packages())



ext1 = Extension(name='mes_from_cpmd.ext_fortran.fortran_io',
                 sources=['mes_from_cpmd/ext_fortran/fortran_io.F90'],
                 f2py_options=['--quiet'],
                )



setup(name='mes_from_cpmd',
     version='19.1',
     packages=find_packages() + ['commit_stuff'],
     author='Christian Dreßler',
     long_description=readme(),
                               entry_points={
                                    'console_scripts': [
                                        "automatic_create_mes = mes_from_cpmd.automatic_create_mes:main",
                                        "direct_moment_expansion = mes_from_cpmd.direct_moment_expansion:main",
                                        "eval_mes = mes_from_cpmd.eval_mes:main",
                                        "wan2cube_cpmd = mes_from_cpmd.wan2cube_cpmd:main",
                                        "gen_pot = mes_from_cpmd.pot.gen_pot:main",
                                        "cal_dens_for_pot = mes_from_cpmd.cal_dens_for_pot:main",
                                        "cal_cube_from_wan_for_dens = mes_from_cpmd.cal_cube_from_wan_for_dens:main",
                                        "cal_diff_dens = mes_from_cpmd.cal_diff_dens:main"]
                                                               },
include_package_data=True,
      zip_safe=False,
       ext_modules=[ext1],
                                    
)

# Write hash, message and date of current commit to file
with open("commit_stuff/version_hash.py", "w") as f:
        commit_hash, commit_message, commit_date = get_commit_hash()
        print("commit_hash = \"{}\"".format(commit_hash), file=f)
        print("commit_message = \"{}\"".format(commit_message), file=f)
        print("commit_date = \"{}\"".format(commit_date), file=f)
