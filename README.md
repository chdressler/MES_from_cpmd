Direct Moment Expansion
================


Description
----------------
Calculation of Moment Expanded states of the static linear density-density response function (LDDRF)  using a the direct moment expansion algorithm 

Requirements
----------------
CPMD exectuable

Installation
----------------

python setup.py install

How To
----------------
execution order:

gen_pot: calculate (here monomial) basis functions of perturbing potential, requires  the DENSITY.cube file from  cpmd single point calcaultion of the target molecule for the extraction of the  grid parameters
.
cal_dens_for_pot: calculate densities for the applied basis functions of the perturbing potential  (n_pot), requires path to the  basis functions of perturbing potential

cal_cube_from_wan_for_dens: converts CPMD DENSITY files into .cube files, Warning: this script has to be executed on an AMD or Intel node (because it requires an old precompiled fortran script)

cal_diff_dens: calculate the difference (n_resp = n_pot - n_eq) of the perturbed densities (n_pot) and the equilibrium densitiy  (n_eq) of the target molecule. requires the path to the perturbed densities (n_pot) and the equilibrium densitiy  (n_eq)

direct_moment_expansion: the most important part of this module. calculates the moment expanded states and performs the actual linear algebra operations. requires the difference densities (n_resp) and the basis functions of the perturbing potential.

eval_mes: evaluates the calculated moment expanded states, requires the moment expanded stated calculated by this module and the moment expanded states obtained from the AG Sebastiani development version of CPMD.

AG Sebastiani members
---------------
complete methanal example: /net/shared/dressler/chi/how_to_MES_from_cpmd

