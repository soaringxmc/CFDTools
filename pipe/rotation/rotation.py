# rotating pipe planes for animations

import scipy.io as sp
import numpy as np
import time

dt = 0.1

with open('bou.in', 'r') as file:
    # Read the second line
    lines = file.readlines()
    parameters = lines[1].split()
    m1 = int(parameters[0])
    m2 = int(parameters[1])
    m3 = int(parameters[2])
    m1m = m1 - 1
    m2m = m2 - 1
    m3m = m3 - 1
    parameters = lines[5].split()
    alx3d = float(parameters[0])
    parameters = lines[9].split()
    ros = float(parameters[1])

with open('radcor.out', 'r') as file:
    rc = []
    rm = []
    for line in file:
        columns = line.split()
        rc.append(float(columns[1]))
        rm.append(float(columns[2]))
    rc = np.array(rc)
    rm = np.array(rm)

th = 2.0*np.pi/m1m*np.arange(m1)        
z  = alx3d*np.arange(m3m)/m3m
# both start from 0 for consistency with the fortran code   

rr = np.empty([2*m2m],dtype=float)
for j in range(2*m2m):
    if j > m2m-1:
        rr[j] =  rm[j-m2m]
    else:
        rr[j] = -rm[m2m-1-j]

jplane = np.empty([4],dtype=np.int32)
jplane[0] = m2m-1
jplane[1] = m2m-16
for j in range(m2m):
    if rm[j] > 0.8:
        jplane[2] = j
        break
for j in range(m2m):
    if rm[j] > 0.5:
        jplane[3] = j
        break

nsta = 20941
nend = 21001
for n in range(nsta,nend):

    print(n)
    # dth = dt*ros*(n-nsta)
    dth = 2.0*np.pi - dt*ros*(nend-1-n)

    # # planetr
    # planetr_x = np.empty([m1,m2m,3],dtype=float) 
    # rm_expand = rm[np.newaxis, :]
    # th_expand = th[:, np.newaxis]
    # planetr_x[:, :, 0] = rm_expand * np.cos(dth+th_expand)
    # planetr_x[:, :, 1] = rm_expand * np.sin(dth+th_expand)
    # planetr_x[:, :, 2] = z[0]
    # with sp.FortranFile('planetr_'+str("%05d" % n)+'.x','w') as f:
        # f.write_record(np.array([m1,m2m,1],np.int32))
        # planetr_x = planetr_x.reshape(m1*m2m*3,order='F')
        # f.write_record(planetr_x)

    # # planerz
    # planerz_x = np.empty([2 * m2m, m3m, 3], dtype=float)
    # planerz_x[:, :, 0] = rr[:, np.newaxis]*np.cos(dth)
    # planerz_x[:, :, 1] = rr[:, np.newaxis]*np.sin(dth)
    # planerz_x[:, :, 2] = z
    # with sp.FortranFile('planerz_'+str("%05d" % n)+'.x','w') as f:
        # f.write_record(np.array([1,2*m2m,m3m],np.int32))
        # planerz_x = planerz_x.reshape(2*m2m*m3m*3,order='F')
        # f.write_record(planerz_x)

    # planetz
    for l in range(1,2):
        j = jplane[l]
        planetz_x = np.empty([m1,m3m,3],dtype=float)
        # planetz_x[:,:,0] = rm[j]*np.cos(dth+th[:,np.newaxis])
        # planetz_x[:,:,1] = rm[j]*np.sin(dth+th[:,np.newaxis])
        # planetz_x[:,:,2] = z
        # unwrap_pos1
        planetz_x[:,:,0] = 1.0*np.cos(dth+th[:,np.newaxis]) - 1.0
        planetz_x[:,:,1] = 1.0*np.sin(dth+th[:,np.newaxis])
        planetz_x[:,:,2] = z
        with sp.FortranFile('planetz_'+str(l+1)+'_unwrap_pos1_'+str("%05d" % n)+'.x','w') as f:
            f.write_record(np.array([m1,1,m3m],np.int32))
            planetz_x = planetz_x.reshape(m1*m3m*3,order='F')
            f.write_record(planetz_x)