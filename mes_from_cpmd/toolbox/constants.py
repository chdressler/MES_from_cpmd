#!/usr/bin/python3
import numpy as np

# misc constants
pi = np.pi # pi
avog = 6.02214129E+23 # Avogadro constant [mol-1]

# S.I. constants
c_si = 2.99792458E+08 # speed of light [m/s]
h_si = 6.62606957E-34 # Planck's constant [Js]
hbar_si = h_si/(2*pi) # reduced planck constant [Js]
m_p_si = 1.67262177E-27 # mass protron [kg]
m_e_si = 9.10938291E-31 # mass electron [kg]
m_amu_si = 1.660538921E-27 # atomic mass unit [kg]
e_si = 1.60217657E-19 # unit charge [C]
eps0_si = 8.854187817E-12 # vacuum permittivity
a0_si = 0.529177210859E-10 # Bohr radius [m]
k_B_si = 1.3806488E-23 # Boltzmann constant [J/K]

# A.U. conversion factors
l_au = a0_si # Bohr radius [m]
E_au = 4.35974417E-18 # Hartree energy [J]
t_au = hbar_si/E_au # time in A.U. [s]
m_p_au = m_p_si/m_e_si # proton mass in A.U. [1] ?
m_amu_au = m_amu_si/m_e_si # atomic mass unit in A.U. [1] ?
c_au = c_si/l_au*t_au # speed of light in A.U. [1]
k_B_au = k_B_si/E_au # Boltzmann constant [E_h/K]

# cgs units
c_cgs = c_si*1E2 # speed of light [cm/s]
e_cgs = e_si*c_si*1E1 # 4.80320427E-10 [esu]
h_cgs = h_si*1E7 # Planck's constant [erg s]
hbar_cgs = h_cgs/(2*pi) # reduced Planck constant [erg s]
a0_cgs = a0_si*1E2 # Bohr radius [cm]
k_B_cgs = k_B_si*1E7 # Boltzmann constant [erg/K]

# misc
a_lat = 2*pi/l_au # lattice constant
finestr = 1/c_au # finestructure constant
l_au2aa = l_au*1E10 # convertion A.U. to Angstrom: x(in au)*l_au2aa = x(in aa)
l_aa2au = 1/l_au2aa # convertion Angstrom to A.U. : x(in aa)*l_aa2au = x(in au)
t_au2fs = t_au*1E15
t_fs2au = 1/t_au2fs
#old
#l_au2aa = 1/l_au*1E-10 # convertion A.U. to Angstrom
#l_aa2au = 1/l_au2aa # convertion Angstrom to A.U.


machine_precision = 1E-14
max_memory_usage = 0.8 # NB! not yet added to test!

# element species # NB! not yet added to test!
species = dict()
species['H']  = {'SYMBOL': 'H', 'MASS':  1.00797, 'Z':  1, 'ZV':  1}
species['D']  = {'SYMBOL': 'H', 'MASS':  2.01410, 'Z':  1, 'ZV':  1}
species['He'] = {'SYMBOL':'He', 'MASS':  4.00260, 'Z':  2, 'ZV':  2}
species['Li'] = {'SYMBOL':'Li', 'MASS':  6.93900, 'Z':  3, 'ZV':  3}
species['Be'] = {'SYMBOL':'Be', 'MASS':  9.01220, 'Z':  4, 'ZV':  4}
species['C']  = {'SYMBOL': 'C', 'MASS': 12.01115, 'Z':  6, 'ZV':  4}
species['N']  = {'SYMBOL': 'N', 'MASS': 14.00670, 'Z':  7, 'ZV':  5}
species['O']  = {'SYMBOL': 'O', 'MASS': 15.99940, 'Z':  8, 'ZV':  6}
species['F']  = {'SYMBOL': 'F', 'MASS': 18.99840, 'Z':  9, 'ZV':  7}
species['Ne'] = {'SYMBOL':'Ne', 'MASS': 20.18300, 'Z': 10, 'ZV':  0}
species['Na'] = {'SYMBOL':'Na', 'MASS': 22.98980, 'Z': 11, 'ZV':  0}
species['P']  = {'SYMBOL':'P',  'MASS': 30.97376, 'Z': 15, 'ZV':  5}
species['S']  = {'SYMBOL':'S',  'MASS': 32.06400, 'Z': 16, 'ZV':  6}
species['Cl'] = {'SYMBOL':'Cl', 'MASS': 35.45300, 'Z': 17, 'ZV':  0}
species['Ca'] = {'SYMBOL':'Ca', 'MASS': 39.96259, 'Z': 20, 'ZV':  0}

# element symbols
symbols = ['H', 'D', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F' , 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca']
