#!/usr/bin/python
# Program to unwarp a cylinder into a plane
import sys
# from scipy.io import FortranFile
import scipy.io as sp
import numpy as np
# input: planetr.bin, planerz.bin, planetz.bin
#        bou.in, flowprop.dat, radcor.out
# output: *.x, *.q in plot3d form

with open('bou.in','r') as f:
    f.readline()
    line  = f.readline()
    value = line.split()
    m1 = int(value[0])
    m2 = int(value[1])
    m3 = int(value[2])
    f.readline()
    f.readline()
    f.readline()
    line  = f.readline()
    value = line.split()
    alx3d = float(value[0])
    f.readline()
    f.readline()
    f.readline()
    line  = f.readline()
    value = line.split()
    reb = float(value[0])

    print(m1,m2,m3,alx3d,reb)

m1m = m1 - 1
m2m = m2 - 1
m3m = m3 - 1

print(m1m,m2m,m3m)

with open('radcor.out','r') as f:
    lines = f.readlines()
    rc = []
    rm = []
    for line in lines:
        rc.append(line.split()[1])
        rm.append(line.split()[2])
    tmp = []
    for num in rc:
        tmp.append(float(num))
    rc = tmp
    tmp = []
    for num in rm:
        tmp.append(float(num))
    rm = tmp


with open('flowprop.dat','r') as f:
    line = f.readline()
    value= line.split()
    utau = float(value[3])
    ttau = float(value[16])
    pw   = float(value[32])

    utau = 0.5   # to get uz/ub

nvar   = 5
nplane = 4

# geometry

th = 2.0*np.pi/m1m*np.arange(m1)
    # start from 0, which should be pi/m1m strictly
    # keep it for consistency with the fortran version
    # so the Uz contour shifts pi/m1m deg, not a big deal                                    
z  = alx3d*np.arange(m3m)/m3m
    # start from 0 and end at alx3d*(m3m-1)/m3m
    # keep it for consistency with the fortran version
rr = np.empty([2*m2m],dtype=float)
for j in range(2*m2m):
    if j > m2m-1:
        rr[j] =  rm[j-m2m]
    else:
        rr[j] = -rm[m2m-1-j]

jplane = np.empty([4],dtype=np.int32)
jplane[0] = m2m-1
jplane[1] = m2m-17
for j in range(m2m):
    if rm[j] > 0.8:
        jplane[2] = j
        break
for j in range(m2m):
    if rm[j] > 0.5:
        jplane[3] = j
        break

# planetr

planetr_x = np.empty([m1,m2m,3],dtype=float)

for i in range(m1):
    for j in range(m2m):
        planetr_x[i,j,0] = rm[j]*np.cos(th[i])
        planetr_x[i,j,1] = rm[j]*np.sin(th[i])
        planetr_x[i,j,2] = z[0]

with sp.FortranFile('planetr.x','w') as f:
    f.write_record(np.array([m1,m2m,1],np.int32))
    planetr_x = planetr_x.reshape(m1*m2m*3,order='F')
    f.write_record(planetr_x)

with sp.FortranFile('planetr.bin','r') as f:
    planetr = f.read_reals(dtype=float)
    planetr = planetr.reshape(m1m,m2m,nvar,order='F')
    planetr_q = np.empty([m1,m2m],dtype=float)
    planetr_q[:m1m,:] = planetr[:m1m,:,2]
    planetr_q[ m1m,:] = planetr[   0,:,2]
    planetr_q = planetr_q/utau

with sp.FortranFile('planetr.q','w') as f:
    f.write_record(np.array([m1,m2m,1,1],np.int32))
    planetr_q = planetr_q.reshape(m1*m2m,order='F')
    f.write_record(planetr_q)

# planerz

planerz_x = np.empty([2*m2m,m3m,3],dtype=float)
for j in range(2*m2m):
    for k in range(m3m):
        planerz_x[j,k,0] = rr[j]*np.cos(th[0])
        planerz_x[j,k,1] = rr[j]*np.sin(th[0])
        planerz_x[j,k,2] = z[k]

with sp.FortranFile('planerz.x','w') as f:
    f.write_record(np.array([1,2*m2m,m3m],np.int32))
    planerz_x = planerz_x.reshape(2*m2m*m3m*3,order='F')
    f.write_record(planerz_x)

planerz = np.fromfile('planerz.bin', dtype=float)
planerz = planerz.reshape(nvar,2*m2m,m3m,order='F')
planerz_q = np.empty([2*m2m,m3m],dtype=float)
planerz_q[:,:] = planerz[2,:,:]
planerz_q = planerz_q/utau

with sp.FortranFile('planerz.q','w') as f:
    f.write_record(np.array([1,2*m2m,m3m,1],np.int32))
    planerz_q = planerz_q.reshape(2*m2m*m3m,order='F')
    f.write_record(planerz_q)

# planetz

for l in range(nplane):
    j = jplane[l]
    planetz_x = np.empty([m1,m3m,3],dtype=float)
    for i in range(m1):
        for k in range(m3m):
            planetz_x[i,k,0] = rm[j]*np.cos(th[i])
            planetz_x[i,k,1] = rm[j]*np.sin(th[i])
            planetz_x[i,k,2] = z[k]
    with sp.FortranFile('planetz_'+str(l+1)+'.x','w') as f:
        f.write_record(np.array([m1,1,m3m],np.int32))
        planetz_x = planetz_x.reshape(m1*m3m*3,order='F')
        f.write_record(planetz_x)

planetz = np.fromfile('planetz.bin', dtype=float)
planetz = planetz.reshape(nvar,nplane,m1m,m3m,order='F')

for l in range(nplane):
    with sp.FortranFile('planetz_'+str(l+1)+'.q','w') as f:
        f.write_record(np.array([m1,1,m3m,1],np.int32))
        planetz_q = np.empty([m1,m3m],dtype=float)
        planetz_q[:m1m,:] = planetz[2,l,:,:]
        planetz_q[ m1m,:] = planetz[2,l,0,:]
        planetz_q = planetz_q/utau
        planetz_q = planetz_q.reshape(m1*m3m,order='F')
        f.write_record(planetz_q)

# unrolled plane
planetz_x = np.empty([m1,m3m,3],dtype=float)
for i in range(m1):
    for k in range(m3m):
        planetz_x[i,k,0] = th[i]
        planetz_x[i,k,1] = 0.0
        planetz_x[i,k,2] = z[k]
with sp.FortranFile('planetz_unrolled.x','w') as f:
    f.write_record(np.array([m1,1,m3m],np.int32))
    planetz_x = planetz_x.reshape(m1*m3m*3,order='F')
    f.write_record(planetz_x)
    

# process of unwrapping a shell into plane
# uncommentted to save time
for l in range(1,2):
    j = jplane[l]
    len_arc = 2.0*np.pi*rm[j]
    for theta in np.arange(2.0*np.pi,0.0,-0.2):
        r = len_arc/theta
        planetz_x = np.empty([m1,m3m,3],float)
        for i in range(m1):
            for k in range(m3m):
                # planetz_x[i,k,0] = r*np.cos(theta/m1m*i) - r
                # planetz_x[i,k,1] = r*np.sin(theta/m1m*i)
                # planetz_x[i,k,2] = z[k]

                planetz_x[i,k,0] = r*np.sin(theta/m1m*i) 
                planetz_x[i,k,1] = -(r*np.cos(theta/m1m*i) - r)
                planetz_x[i,k,2] = -z[k]
                
        num = int((628-int(theta*100))/20)
        print(num,round(theta,2))
        with sp.FortranFile('planetz_'+str(l+1)+'_'+str("%05d" % num)+'.x','w') as f:
            f.write_record(np.array([m1,1,m3m],np.int32))
            planetz_x = planetz_x.reshape(m1*m3m*3,order='F')
            f.write_record(planetz_x)
