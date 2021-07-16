#!/usr/bin/python3
import numpy as np
import gzip

def ReadCubeFile(filename):
    f = open(filename, 'r')
    inbuffer = f.read()
    f.close()
    data, numbers, coords, cell, comment1, comment2, origin = ParseCubeFile(inbuffer)
    return data, numbers, coords, cell, comment1, comment2, origin

def WriteCubeFile(filename, comment1, comment2, numbers, coords, cell, data, origin=np.zeros(3)):
    outbuffer = AssembleCubeFile(comment1, comment2, numbers, coords, cell, data, origin=origin)
    f = open(filename, 'w')
    f.write(outbuffer)
    f.close()
    # f.write(comment1.rstrip('\n').replace('\n','')+'\n')
    # f.write(comment2.rstrip('\n').replace('\n','')+'\n')
    # dim = list(data.shape)
    # n_atoms = coords.shape[0]
    # f.write('   %2d  %10.6f  %10.6f  %10.6f\n'%(n_atoms, 0,0,0))
    # for i in range(3):
    #     f.write('   %2d  %10.6F  %10.6f  %10.6f\n'%(dim[i], cell[i][0], cell[i][1], cell[i][2]))
    # for atom in range(n_atoms):
    #     f.write('   %2d  %10.6f  %10.6f  %10.6f  %10.6f\n'%(numbers[atom], numbers[atom], coords[atom][0], coords[atom][1], coords[atom][2]))
    # for i_x in range(dim[0]):
    #     for i_y in range(dim[1]):
    #         for i_z in range(dim[2]):
    #             if i_z % 6 == 0 and i_z != 0:
    #                 f.write('\n')
    #             f.write('%13.5E'%data[i_x][i_y][i_z])
    #         f.write('\n')
    # f.close()

def ReadGZIPCubeFile(filename):
    f = gzip.open(filename, 'rb')
    inbuffer = f.read().decode('UTF-8')
    f.close()
    data, numbers, coords, cell, comment1, comment2, origin = ParseCubeFile(inbuffer)
    return data, numbers, coords, cell, comment1, comment2, origin

def WriteGZIPCubeFile(filename, comment1, comment2, numbers, coords, cell, data, origin=np.zeros(3)):
    outbuffer = AssembleCubeFile(comment1, comment2, numbers, coords, cell, data, origin=origin)
    f = gzip.open(filename, 'wb')
    f.write(bytes(outbuffer, 'UTF-8'))
    f.close()

def ParseCubeFile(inbuffer):
    inbuffer = inbuffer.split('\n')
    comment1 = inbuffer[0]
    comment2 = inbuffer[1]
    tmp = inbuffer[2].split()
    n_atoms = int(tmp[0])
    origin = np.array([float(e) for e in tmp[1:4]])
    cell_x = inbuffer[3].split()
    cell_y = inbuffer[4].split()
    cell_z = inbuffer[5].split()
    n_x, a_lat_x = int(cell_x[0]), np.array([float(e) for e in cell_x[1:4]])
    n_y, a_lat_y = int(cell_y[0]), np.array([float(e) for e in cell_y[1:4]])
    n_z, a_lat_z = int(cell_z[0]), np.array([float(e) for e in cell_z[1:4]])
    cell = np.vstack((a_lat_x, a_lat_y, a_lat_z))
    coords = np.zeros((n_atoms, 3))
    numbers = np.zeros(n_atoms)
    for atom in range(n_atoms):
        tmp = inbuffer[6+atom].split()
        numbers[atom] = int(tmp[0])
        coords[atom] = np.array([float(e) for e in tmp[2:5]])
    data = np.zeros((n_x, n_y, n_z))
    counter = 6+n_atoms
    for i_x in range(n_x):
        for i_y in range(n_y):
            tmp = ''
            if n_z % 6 == 0:
                z_count = n_z//6
            else:
                z_count = n_z//6 + 1
            for i in range(z_count):
                tmp += inbuffer[counter]
                counter += 1
            data[i_x][i_y] = np.array(tmp.strip().split())
    return data, numbers, coords, cell, comment1, comment2, origin

def AssembleCubeFile(comment1, comment2, numbers, coords, cell, data, origin=np.zeros(3)):
    obuffer = ''
    obuffer += comment1.rstrip('\n').replace('\n','')+'\n'
    obuffer += comment2.rstrip('\n').replace('\n','')+'\n'
    dim = list(data.shape)
    n_atoms = coords.shape[0]
    obuffer += '   %2d  %10.6f  %10.6f  %10.6f\n'%(n_atoms, origin[0], origin[1], origin[2])
    for i in range(3):
        obuffer += '   %2d  %10.6F  %10.6f  %10.6f\n'%(dim[i], cell[i][0], cell[i][1], cell[i][2])
    for atom in range(n_atoms):
        obuffer += '   %2d  %10.6f  %10.6f  %10.6f  %10.6f\n'%(numbers[atom], numbers[atom], coords[atom][0], coords[atom][1], coords[atom][2])
    for i_x in range(dim[0]):
        for i_y in range(dim[1]):
            for i_z in range(dim[2]):
                if i_z % 6 == 0 and i_z != 0:
                    obuffer += '\n'
                obuffer += '%13.5E'%data[i_x][i_y][i_z]
            obuffer += '\n'
    return obuffer
