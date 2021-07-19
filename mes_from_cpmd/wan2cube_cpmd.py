import mes_from_cpmd.ext_fortran.fortran_io as fio
from mes_from_cpmd.toolbox import CubeFileTools
from mes_from_cpmd.toolbox import transformations

import numpy as np
from mes_from_cpmd.toolbox import lib_dme as lime
import ipdb
from mes_from_cpmd.misc import git_control
import mes_from_cpmd.toolbox.cube as cube

import numpy as np
import os
import sys
import argparse
import subprocess






def main():
    parser=argparse.ArgumentParser(description="convert wan file created by CPMD to cub file", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("fn_cube_in",  help="Cubefile filename for header")
    parser.add_argument("fn_wan_in",   help="dat filename with volume data")
    parser.add_argument("fn_cube_out", help="Cubefile filename for output")
    parser.add_argument("--verbose",   default=False, action='store_true', help="Verbose output")
    args = vars(parser.parse_args())
    if args['verbose']:
        print(args)

    cell_data  = CubeFileTools.LoadCellData(args['fn_cube_in'])
    state_data = np.asfortranarray(np.zeros(cell_data['mesh'], dtype=np.float64))
    n_x, n_y, n_z = cell_data['mesh']
    fio.fortran_read_unformatted(args['fn_wan_in'], state_data, n_x, n_y, n_z)
    cube.WriteCubeFile(args['fn_cube_out'],
                       cell_data['comment1' ],
                       cell_data['comment2' ],
                       cell_data['numbers'  ],
                       cell_data['coords_au'],
                       cell_data['cell_au'  ],
                       state_data,
                       cell_data['origin_au'])
