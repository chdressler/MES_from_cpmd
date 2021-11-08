Direct Moment Expansion
================


Description
----------------
Calculation of Moment Expanded states of the static linear density-density response function (LDDRF)  using a the direct moment expansion algorithm 

Requirements
----------------
CPMD and cpmd2cube exectuable

Installation
----------------

python setup.py install

How To
----------------
 
calculation of mes (moment expanded states) for methanal:

path\_to\_repo=path\_to\_this\_directory

path\_to\_example=$path\_to\_repo/MES\_from\_cpmd/examples/automatic\_generation\_of\_me

cd $path\_to\_example 

automatic\_create\_mes $path\_to\_example/input\_files/wfo.inp  $path\_to\_example/input\_files/calc.job 10

and for benchmarking execute:

eval\_mes single\_point/DENSITY\_fullmesh.cube  moment\_expanded\_states ../referencs\_mes potentials 9
