import  numpy as np
#import fortran_io as fio
import mes_from_cpmd.ext_fortran.fortran_io as fio
from mes_from_cpmd.toolbox import CubeFileTools
#from  chitools import transformations
from  mes_from_cpmd.toolbox import transformations
import scipy
from mes_from_cpmd.toolbox import cube
#from reduced_eigen import elww

def over(xx, yy):
    return np.sum(xx * yy)

#def err(xx, yy):
#    return abs(over(xx, yy)) ** 0.5/ np.sum(yy**2)**0.5

def err(xx, yy):
    return np.sum((xx - yy)**2) ** 0.5/ np.sum(yy**2)**0.5


def shape_err(xx, yy):
    xx1 = xx / over(xx, xx)**0.5
    yy1 = yy / over(yy, yy)**0.5
    
    return abs(over(xx1, yy1))   


def err_mult(xx, yy):
    return abs(np.sum((xx * yy))) ** 0.5/ np.sum(yy**2)**0.5
    #xx1 = xx / over(xx, xx)**0.5
    #yy1 = yy / over(yy, yy)**0.5
    #return abs(over(xx1, yy1))   
def over(xx, yy):
    return np.sum(xx * yy)

def create_overlap_mat(*states4 ):
    #if states2.all() ==  None:
    #import ipdb
    #ipdb.set_trace()
    if len(states4) == 1:
        states = states4[0]
        n_states = states.shape[0]
        over_mat = np.zeros((n_states,n_states))
        for i_state in range(n_states):
            for j_state in range(n_states):
                over_mat[i_state,j_state] =  over(states[i_state], states[j_state])
        return over_mat  
    else:
        #import ipdb
        #ipdb.set_trace()
        states = states4[0]
        states2 = states4[1]
        #states2  = np.array(states3)
        n_states = states.shape[0]
        m_states = states2.shape[0]
        over_mat = np.zeros((n_states,m_states))
        for i_state in range(n_states):
            for j_state in range(m_states):
                over_mat[i_state,j_state] =  over(states[i_state], states2[j_state])
        return over_mat  


def pot_func1(data,scale,origin,diff_len):
   for i in range(data.shape[0]):
   #   for j in range(data.shape[0]):
   #      for k in range(data.shape[0]):
             data[i ,:, : ] = (i*diff_len -  origin[2]) * scale 
   return data
#def origin

def pot_func_x(data,scale,origin_single,diff_len):
    for i in range(data.shape[0]):
        data[i ,:, : ] = (i*diff_len -  origin_single) * scale
    return data

def pot_func_y(data,scale,origin,diff_len):
    for i in range(data.shape[0]):
        data[: ,i, : ] = (i*diff_len -  origin[1]) * scale
    return data


def pot_func_z(data,scale,origin,diff_len):
    for i in range(data.shape[0]):
        data[: ,:, i ] = (i*diff_len -  origin[1]) * scale
    return data



def re_cube(data):
   data_new = np.zeros(data.shape)
   list1 = []
   for i in range(data.shape[0]):
       for j in range(data.shape[1]):
           for k in range(data.shape[2]):
               list1.append(data[i,j,k])
   
   x = -1
   y = -1
   z= -1
   for l in range(len(list1)):
       #if y % data.shape[1]  == 0:
       if l % (data.shape[0]**2)  == 0:
           y = -1
           z += 1

       x += 1
       if l % data.shape[0] == 0:
           x = 0
           y += 1
       #if l < 140067:
       # print('index')   
       # print(l,x,y,z)
       data_new[x,y,z] = list1[l]
   return data_new

       #data.shape[0]
def pot_inv(data):
   data_new = np.zeros(data.shape)
   for i in range(data.shape[0]):
       for j in range(data.shape[1]):
           for k in range(data.shape[2]):
               data_new[k,j,i] = data[i,j,k]
   return data_new


def norm_states(statesT):
    states = np.copy(statesT)
    n_states = states.shape[0]
    for i_state in range(n_states):
        states[i_state] = states[i_state] / over(states[i_state], states[i_state])**0.5 
    return states

#def basis_expansion(coeff, states):
#    if len(coeff.shape) ==  1:
#         print('coeff vektor')
#         n_states = coeff.shape[0]
#         print('!!!!!!!!!!!!!!!!!!!')
#         print('basis',n_states)
#         calc_state = np.zeros(states[0].shape)
#         for i_states in range(n_states):
#             print(coeff[i_states])
#             calc_state += coeff[i_states] * states[i_states]
#         #moglich ohne for schleif in einer zeile
#         return calc_state
#     elif len(coeff.shape) ==  2:
#         print('coeff matrix')
#         for j in range(coeff.shape[1]):
#             coeff
       
def basis_expansion(coeff, states):
    n_states = coeff.shape[0]
    print('!!!!!!!!!!!!!!!!!!!')
    print('basis',n_states)
    calc_state = np.zeros(states[0].shape)
    for i_states in range(n_states):
        print(coeff[i_states])
        calc_state += coeff[i_states] * states[i_states]
    #moglich ohne for schleif in einer zeile
    return calc_state

def eval_dens(data_calc, data_ref, string1 = None):
    if string1 != None:
        print('#########################')
        print(string1)

    print('zero')
    print(np.sum(data_ref))
    print(np.sum(data_calc))
    print('square')
    print(np.sum(data_ref**2))
    print(np.sum(data_calc**2))
    print('abs')
    print(np.sum(abs(data_ref)))
    print(np.sum(abs(data_calc)))
    print(np.sum(abs(data_ref))/np.sum(abs(data_calc)))
    print("err subtraktiv", err(data_calc, data_ref ))
    print("overlap of normed states", shape_err(data_calc, data_ref ))
    print("normed overlap error", err_mult(data_calc, data_ref))
  
#print(


#fn_cube = '/home/dressler/promotion/chi/H2O-geos-speilwiese/direct_calc-modify/q1_0.9923-q2_0.9923-q3_1.7884/test.cube'
#cell_data = CubeFileTools.LoadCellData(fn_cube)
#print(cell_data.keys())
#bn_states={}
#bn_states['eq'] = '/shared/dressler/chi/2019-neue_mol_dme/hf/chi1/moment_expanded_state_name-%05d.wan'
#states={}
#for key in ['eq']:
# n_x, n_y, n_z = cell_data['mesh']
# state_data    = np.asfortranarray(np.zeros((n_x, n_y, n_z), dtype=np.float64))
# n_states      = 9
# states[key]        = np.zeros((n_states, n_x, n_y, n_z))
# for i_state in range(n_states):
#     fio.fortran_read_unformatted(bn_states[key]%(i_state+1), state_data, n_x, n_y, n_z)
#     states[key][i_state] = state_data/np.sqrt(n_x*n_y*n_z)
#     print(key,i_state, np.sum(states[key][i_state]**2))
#     
##bn_states_out = '/shared/dressler/chi/2019-neue_mol_dme/hf/cube1/normal-states-%s-mes-%05d.wan'
##bn_cube_out = '/shared/dressler/chi/2019-neue_mol_dme/hf/cube1/normal-states-%s-mes-%05d.cube'
##for key in ['eq']:
## for i_state in range(n_states):
#bn_states={}
##bn_states['eq'] = '/shared/dressler/chi/2019-neue_mol_dme/dimer/chi_from_1/chi1/moment_expanded_state_name-%05d.wan'
##bn_states['eq'] = '/shared/dressler/chi/2019-neue_mol_dme/dimer/chi_from_1/chi3_old_school/mes-%05d.wan'
##bn_states['dme'] = '/shared/dressler/chi/2019-neue_mol_dme/dimer/chi_from_1/chi1/chidm_state_name-%05d.wan'
#bn_states['dme'] = '/shared/dressler/chi/2019-neue_mol_dme/dimer/chi_from_1/chi3_old_school/mes-%05d.wan'
#bn_states['eq'] = '/shared/dressler/chi/2019-neue_mol_dme/dimer/chi_from_1/chi1/chidm_state_name-%05d.wan'
#states={}
#for key in ['eq', 'dme']:
# n_x, n_y, n_z = cell_data['mesh']
# state_data    = np.asfortranarray(np.zeros((n_x, n_y, n_z), dtype=np.float64))
# n_states      = 34
# states[key]        = np.zeros((n_states, n_x, n_y, n_z))
# calc_dens        = np.zeros(( n_x, n_y, n_z))
# for i_state in range(n_states):
#     fio.fortran_read_unformatted(bn_states[key]%(i_state+1), state_data, n_x, n_y, n_z)
#     states[key][i_state] = state_data/np.sqrt(n_x*n_y*n_z)
#     #states[key][i_state] = state_data
#     #states[key][i_state] = state_data*con_vol
#     #states[key][i_state] = state_data*con_vol
#     print(key,i_state, np.sum(states[key][i_state]**2), np.sum(abs(states[key][i_state])))
#
##CentersOfMass
##ELPOT_unnormed.cube
##/shared/dressler/chi/paper-no2/basis_fit/linear-nachgesampelt-2015/q1_0.9923-q2_0.9923-q3_1.7884/cpmd-density/DENSITY.cube
##fn_cube = '/shared/dressler/chi/paper-no2/basis_fit/linear-nachgesampelt-2015/q1_0.9923-q2_0.9923-q3_1.7884/cpmd-density/DENSITY.cube'
##fn_cube = '/home/dressler/promotion/chi/H2O-geos-speilwiese/direct_calc-modify/q1_0.9923-q2_0.9923-q3_1.7884/test.cube'
#fn_cube = '/net/shared/dressler/chi/2019-neue_mol_dme/dimer/geo3/DENSITY.cube'
#fn_cube0 = '/shared/dressler/chi/2019-neue_mol_dme/dimer/chi_from_1/chi1/DENSITY.cube'
##fn_cube0 = '/shared/dressler/chi/2019-neue_mol_dme/dimer/cpmd_calc_response_for_linear_pot_x/2_real_coord_eq/DENSITY.cube'
##fn_cube1 = '/shared/dressler/chi/2019-neue_mol_dme/dimer/2-pot_from_2/post/DENSITY.cube'
##fn_cube1 = '/shared/dressler/chi/2019-neue_mol_dme/dimer/cpmd_calc_response_for_linear_pot_y/pert1/DENSITY.cube'
##fn_cube1 = '/shared/dressler/chi/2019-neue_mol_dme/dimer/cpmd_calc_response_for_linear_pot_x/pert4_scale_0.1/DENSITY.cube'
##fn_cube1 = '/shared/dressler/chi/2019-neue_mol_dme/dimer/cpmd_calc_response_for_linear_pot_x/pert5_scale_0.001/DENSITY.cube'
##fn_cube1 = '/shared/dressler/chi/2019-neue_mol_dme/dimer/cpmd_calc_response_for_linear_pot_x/real_coord_0.001/DENSITY.cube'
##fn_cube1 = '/shared/dressler/chi/2019-neue_mol_dme/dimer/cpmd_calc_response_for_linear_pot_x/2real_coord_0.001/DENSITY.cube'
#fn_cube1 = '/shared/dressler/chi/2019-neue_mol_dme/dimer/cpmd_calc_response_for_linear_pot_y/pert2_scale_0.01/DENSITY.cube'
##fn_cube1 = '/shared/dressler/chi/2019-neue_mol_dme/dimer/cpmd_calc_response_for_linear_pot_y/pert3_scale_0.0001/DENSITY.cube'
##fn_cube1 = '/shared/dressler/chi/2019-neue_mol_dme/dimer/cpmd_calc_response_for_linear_pot_z/pert1/DENSITY.cube'
#fn_cube2 = '/shared/dressler/chi/2019-neue_mol_dme/dimer/2-pot_from_2/post/cpmd-script/ELPOT.cube'
#fn_cube3 = '/shared/dressler/chi/2019-neue_mol_dme/dimer/2-pot_from_2/post/cpmd-script/ELPOT_unnormed.cube'
#cell_data = CubeFileTools.LoadCellData(fn_cube)
#cell_data0 = CubeFileTools.LoadCellData(fn_cube0)
#cell_data1 = CubeFileTools.LoadCellData(fn_cube1)
#pot = CubeFileTools.LoadCellData(fn_cube2)
#pot_un = CubeFileTools.LoadCellData(fn_cube3)
#print(np.sum(pot['data']**2))
#print(np.sum(pot_un['data']**2))
#diff_pot = pot['data']- pot_un['data']
#print(np.sum(diff_pot**2))



def load_states(fn_cube, bn_states, n_states, pure =  False, pert = False ):
    #bn_cube_out1 = bn_out + '.cube'
    #fn_cube = '/net/shared/dressler/chi/2019-neue_mol_dme/dimer/geo3/DENSITY.cube'
    cell_data = CubeFileTools.LoadCellData(fn_cube)
    n_x, n_y, n_z = cell_data['mesh']
    state_data    = np.asfortranarray(np.zeros((n_x, n_y, n_z), dtype=np.float64))
    #n_states      = 34
    states        = np.zeros((n_states, n_x, n_y, n_z))
    #calc_dens        = np.zeros(( n_x, n_y, n_z))
    for i_state in range(n_states):
       if pert:
           fio.fortran_read_pert(bn_states%(i_state+1), state_data, n_x, n_y, n_z)
       else:
           fio.fortran_read_unformatted(bn_states%(i_state+1), state_data, n_x, n_y, n_z)
       if pure :
           states[i_state] = state_data
       else:
           states[i_state] = state_data/np.sqrt(n_x*n_y*n_z)
       print(i_state, np.sum(states[i_state]**2), np.sum(abs(states[i_state])))
    return states

def load_states_single(fn_cube, bn_states, pure =  False, pert = False):
    #bn_cube_out1 = bn_out + '.cube'
    #fn_cube = '/net/shared/dressler/chi/2019-neue_mol_dme/dimer/geo3/DENSITY.cube'
    cell_data = CubeFileTools.LoadCellData(fn_cube)
    n_x, n_y, n_z = cell_data['mesh']
    state_data    = np.asfortranarray(np.zeros((n_x, n_y, n_z), dtype=np.float64))
    #n_states      = 34
    states        = np.zeros(( n_x, n_y, n_z))
    #calc_dens        = np.zeros(( n_x, n_y, n_z))
    #for i_state in range(n_states):
    if pert:
        fio.fortran_read_pert(bn_states, state_data, n_x, n_y, n_z)
    else:
        fio.fortran_read_unformatted(bn_states, state_data, n_x, n_y, n_z)
    if pure :
        states = state_data
        print('true aa')
    else:    
        states = state_data/np.sqrt(n_x*n_y*n_z)
    print(np.sum(states**2), np.sum(abs(states)))
    return states



def print_states(data, cell_data, bn_out = '/storage/state', pure = False, pert = True):
    bn_cube_out1 = bn_out + '.cube'
    bn_states_out1 = bn_out + '.wan'
    n_x, n_y, n_z = cell_data['mesh']
    #state_data1 = np.asfortranarray(data).astype(np.float64)
    if pure:
         state_data1 = np.asfortranarray(data).astype(np.float64)
         print('true aa')
    else:
         state_data1 = np.asfortranarray(data*np.sqrt(n_x*n_y*n_z)).astype(np.float64)
    if pert:
        fio.fortran_write_pert(bn_states_out1, state_data1, n_x, n_y, n_z)
    else:
        fio.fortran_write_unformatted(bn_states_out1, state_data1, n_x, n_y, n_z)
    cube.WriteCubeFile(bn_cube_out1,'','', cell_data['numbers'], cell_data['coords_au'], cell_data['cell_au'], data, origin=cell_data['origin_au'])

#def print_states(data, cell_data, bn_out = '/storage/state', pure = False):
#    bn_cube_out1 = bn_out + '.cube'
#    bn_states_out1 = bn_out + '.wan'
#    n_x, n_y, n_z = cell_data['mesh']
#    #state_data1 = np.asfortranarray(data).astype(np.float64)
#    if pure:
#         state_data1 = np.asfortranarray(data).astype(np.float64)
#         print('true aa')
#    else:
#         state_data1 = np.asfortranarray(data*np.sqrt(n_x*n_y*n_z)).astype(np.float64)
#    fio.fortran_write_unformatted(bn_states_out1, state_data1, n_x, n_y, n_z)
#   cube.WriteCubeFile(bn_cube_out1,'','', cell_data['numbers'], cell_data['coords_au'], cell_data['cell_au'], data, origin=cell_data['origin_au'])

def print_states_cube(data, cell_data, bn_out = '/storage/state'):
    bn_cube_out1 = bn_out + '.cube'
    #bn_states_out1 = bn_out + '.wan'
    #n_x, n_y, n_z = cell_data['mesh']
    #state_data1 = np.asfortranarray(data).astype(np.float64)
    #state_data1 = np.asfortranarray(data*np.sqrt(n_x*n_y*n_z)).astype(np.float64)
    #fio.fortran_write_unformatted(bn_states_out1, state_data1, n_x, n_y, n_z)
    cube.WriteCubeFile(bn_cube_out1,'','', cell_data['numbers'], cell_data['coords_au'], cell_data['cell_au'], data, origin=cell_data['origin_au'])


def print_states_dict(states, cell_data, n_states, bn_out, pure = False, pert = False):
    bn_cube_out = bn_out + '.cube'
    bn_states_out = bn_out + '.wan'
    n_x, n_y, n_z = cell_data['mesh']
    for i_state in range(n_states):
        if pure:
             state_data = np.asfortranarray(states[i_state]).astype(np.float64)
        else:
             state_data = np.asfortranarray(states[i_state]*np.sqrt(n_x*n_y*n_z)).astype(np.float64)
        #fio.fortran_write_unformatted(bn_states_out%((i_state+1)), state_data, n_x, n_y, n_z)
        cube.WriteCubeFile(bn_cube_out%((i_state+1)),'','', cell_data['numbers'], cell_data['coords_au'], cell_data['cell_au'], states[i_state], origin=cell_data['origin_au'])
   
        if pert:
            fio.fortran_write_pert(bn_states_out%((i_state+1)), state_data, n_x, n_y, n_z)
        else:
            fio.fortran_write_unformatted(bn_states_out%((i_state+1)), state_data, n_x, n_y, n_z)


#print1 = 0
##bn_out = 'storage_lin_z/z_pot_0.0001'
##print_states(bn_out, y_pot ,cell_data)
#
#
##bn_out = 'storage_act/density-diff'
##print_states(bn_out, diff_dens ,cell_data)
#
#
#
#
#bn_states_out = '/shared/dressler/chi/2019-neue_mol_dme/hf/cube1/reduced-states-%s-mes-%05d.wan'
#bn_cube_out = '/shared/dressler/chi/2019-neue_mol_dme/hf/cube1/reduced-states-%s-mes-%05d.cube'
##for key in ['eq']:
#for i_state in range(n_states):
#   #state_data = np.asfortranarray(states[key][i_state]*np.sqrt(n_x*n_y*n_z)).astype(np.float64)
#   state_data = np.asfortranarray(reduced_states[i_state]*np.sqrt(n_x*n_y*n_z)).astype(np.float64)
#   fio.fortran_write_unformatted(bn_states_out%(key,(i_state+1)), state_data, n_x, n_y, n_z)
#   cube.WriteCubeFile(bn_cube_out%(key,(i_state+1)),'','', cell_data['numbers'], cell_data['coords_au'], cell_data['cell_au'], state_data, origin=cell_data['origin_au'])
#
charge = [-0.834, 0.417, 0.417]



#def tip3p(data, coords, charge, diff_len):
def tip3p(data, coords, charge, diff_len):
     for i in range(data.shape[0]):
         for j in range(data.shape[1]):
             for k in range(data.shape[2]):
                 ref = np.array([i,j,k])*diff_len
                 for l in range(3):
                     #print(coords[l,:], ref)
                     dist1 =  np.linalg.norm(coords[l,:] - ref)   
                     #data[i,j,k] +=1.0/dist1 * charge[l] * scipy.special.erf(dist1/(2**0.5 * 0.5 * 0.5*0.529177210))
                     data[k,j,i] +=1.0/dist1 * charge[l] * scipy.special.erf(dist1/(2**0.5 * 0.5*0.529177210))
     return data


# data = leerer array mit vorgegebene dimensionen for storpotential
#coords: koordinaten array von storenden molcule
#com: com of coords
#charge: ladungen der atome in coords
#diff_len: gitterkonstante
#com_box: com of coord wird in den ursprung gesetzt, com_box gibt mittelpunkt der box an, in der daserzeigt  Potential gespeichtert wird
def tip3p_oktaeder(data, coords, com, charge, diff_len, com_box):
    #sorgt dafur, dass mittelpunkt der erzeugten potential box com_box entspricht
    off_set = np.array(com_box - np.array([data.shape[0]/2, data.shape[1]/2, data.shape[2]/2])*diff_len)
    #verschiebt com von molekul mit stoerdichte in [0,0,0]
    coords = coords - com
    for i in range(data.shape[0]):
         for j in range(data.shape[1]):
             for k in range(data.shape[2]):
                 ref = np.array([i,j,k])*diff_len + off_set
                 if i == 62 and j == 59 and k == 52:
                     print('distance at', i,j,k, 'ist' ,ref)
                 for l in range(3):
                     #print(coords[l,:], ref)
                     dist1 =  np.linalg.norm(coords[l,:] - ref)   
                     #data[k,j,i] +=1.0/dist1 * charge[l] * scipy.special.erf(dist1/(2**0.5 * 0.5*0.529177210))
                     data[k,j,i] +=1.0/dist1 * charge[l] * scipy.special.erf(dist1/(2**0.5 * 0.5*0.529177210))
    return data



# data = leerer array mit vorgegebene dimensionen for storpotential
#coords: koordinaten array von storenden molcule
#com: com of coords
#charge: ladungen der atome in coords
#diff_len: gitterkonstante
#com_box: com of coord wird in den ursprung gesetzt, com_box gibt mittelpunkt der box an, in der daserzeigt  Potential gespeichtert wird
def tip3p_oktaeder_point(data, coords, com, charge, diff_len, com_box):
    #sorgt dafur, dass mittelpunkt der erzeugten potential box com_box entspricht
    off_set = np.array(com_box - np.array([data.shape[0]/2, data.shape[1]/2, data.shape[2]/2])*diff_len)
    #verschiebt com von molekul mit stoerdichte in [0,0,0]
    coords = coords - com
    for i in range(data.shape[0]):
         for j in range(data.shape[1]):
             for k in range(data.shape[2]):
                 ref = np.array([i,j,k])*diff_len + off_set
                 if i == 62 and j == 59 and k == 52:
                     print('distance at', i,j,k, 'ist' ,ref)
                 
                 for l in range(len(charge)):
                     #print(coords[l,:], ref)
                     dist1 =  np.linalg.norm(coords[l,:] - ref)   
                     #data[k,j,i] +=1.0/dist1 * charge[l] * scipy.special.erf(dist1/(2**0.5 * 0.5*0.529177210))
                     #data[k,j,i] +=1.0/dist1 * charge[l] * scipy.special.erf(dist1/(2**0.5 * 0.5*0.529177210))
                     data[k,j,i]  +=1.0/dist1 * charge[l]* scipy.special.erf(dist1/(2**0.5 * 0.5*0.529177210)) 
    return data

#def coulomb_pot_ang(dist1, charge1):
#    epsilon = 8.8541878128 * 10**(-12)
##    e = 1.602176634 * 10**(-19)
#    return charge1  * e / ( epsilon * 4 *  np.pi * dist1 *10**(-10) )
