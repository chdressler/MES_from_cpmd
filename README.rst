Installation
================

python seup.py install


create DENSITY.cube file
------------

calculate elcectron density file (DENSITY.cube) of target molecule (for example by  CPMD single point calculation and cpmdtocube.x --fullmesh DENSITY subsequently)


calculate basis function of perturbing potential on a grid
------------
execute:
gen_pot --path_to_DENSITY


