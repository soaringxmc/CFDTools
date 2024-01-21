import sys

import scipy.io as sp
import numpy as np
import time

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

dth = 0.1*ros*(21500-20941)

jpl = int(input())
r = rm[m2m-1-jpl]
pi2 = 2.0*np.pi
i_values = np.arange(m1)
k_values = np.arange(m3m)
planetz_x = np.empty([m1,m3m,3],float)
planetz_x[:, :, 0] = r*np.cos(dth+pi2/m1m*i_values[:,np.newaxis])
planetz_x[:, :, 1] = r*np.sin(dth+pi2/m1m*i_values[:,np.newaxis])
planetz_x[:, :, 2] = z[np.newaxis,k_values]
with sp.FortranFile('planetz_2_m2mp'+str("%05d" % jpl)+'.x','w') as f:
    f.write_record(np.array([m1,1,m3m],np.int32))
    planetz_x = planetz_x.reshape(m1*m3m*3,order='F')
    f.write_record(planetz_x)