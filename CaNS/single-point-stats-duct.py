import numpy as np
import sys
import os
import re
import glob
#
# compute mean pressure gradient
#
use_new_arrays = False
shuffle_data   = False
args = sys.argv[1:]
casedir = args[3] if len(args) > 3 else ''
datadir = casedir+''
print('Data directory: ' + datadir+'\n')
resultsdir = casedir+'results/'
os.makedirs(resultsdir,exist_ok=True)
#
# import parameters
#
sys.path.append(casedir) # put input.py in the case directory
from input import *
visc = visci**(-1)
reb = 2*h*ub/visc
tbeg   = float(args[0]) if len(args) > 0 else tbeg
tend   = float(args[1]) if len(args) > 1 else tend
fldstp = int(  args[2]) if len(args) > 2 else fldstp
fname_dpdx = datadir+'forcing.out'
data = np.loadtxt(fname_dpdx)
tarr, ind = np.unique(data[:,0], return_index=True)
if(shuffle_data):
    print('Note: dataset will be randomly shuffled before fixing the time window\n')
    nlen = np.size(data[ind,1][np.where((tarr>tbeg) & (tarr<tend))])
    dpdx_arr = data[ind,1][np.where(tarr>tbeg)]
    np.random.shuffle(dpdx_arr[:])
else:
    dpdx_arr = data[ind,1][np.where((tarr>tbeg) & (tarr<tend))]
    nlen = np.size(dpdx_arr[:])
dpdx = -np.average(dpdx_arr[0:nlen])
#
utau = np.sqrt(dpdx*h/2.)
retau = utau*h/visc
dnu   = visc/utau/h
cf    = utau**2/(ub**2/2.)
np.savetxt(resultsdir+'stats.txt',np.c_[retau,utau,dnu])
#
print("Pressure gradient = ", dpdx)
print("u_tau/u_bulk = ", utau)
print("dnu/h        = ", dnu)
print("Friction Reynolds number = ", retau)
#
# compute single point statistics
#
saves_log = np.loadtxt(datadir+'time.out')
step_log, ind = np.unique(saves_log[:,0], return_index=True)
time_log = saves_log[ind,2]
def find_closest_save(time_log,step_log,time,stp):
    index = np.searchsorted(time_log, time, side="right")
    if(index < np.size(step_log[:])):
        isave = int(step_log[index])-stp//2
        isave = isave - isave%stp
    else:
        isave = -1
    return isave
fldbeg = find_closest_save(time_log,step_log,tbeg,fldstp)
fldend = find_closest_save(time_log,step_log,tend,fldstp)
files = np.sort(glob.glob(datadir+'velstats_fld_???????.out'))
lastfile = files[-1]
fldmax = list(map(int,re.findall(r'([0-9]{7,})', lastfile)))[-1]
if(fldend == -1):
    fldend = fldmax
    index = np.searchsorted(step_log[:], fldmax, side="left")
    print("\nN.B.: Maximum time exceeds number of saves, overwriting end time...")
    tend = time_log[index]
firstfile= files[0]
fldmin = list(map(int,re.findall(r'([0-9]{7,})',firstfile)))[-1]
if(fldmin > fldbeg):
    fldbeg = fldmin
    index = np.searchsorted(step_log[:], fldmin, side="left")
    print("\nN.B.: Minimum larger than the first found save, overwriting  time...")
    tbeg = time_log[index]
if(fldbeg > fldend):
    print("\nError in the selected time interval. Aborting...")
    exit()
t_dns = 1.
tb    = (2.*h)/ub
te    = h/utau
print("Max, min, and interval averaging time in half bulk times (h/ub) = ", tbeg         , tend         , (tend-tbeg)         ,"\n")
print("Max, min, and interval averaging time in turnover times         = ", tbeg*t_dns/te, tend*t_dns/te, (tend-tbeg)*t_dns/te,"\n")
#
fname_beg = datadir+'velstats_fld_'
fname_ext = '.out'

nflds = (fldend-fldbeg)//fldstp + 1
flds = np.array(range(fldbeg,fldend+fldstp,fldstp))
if(shuffle_data):
    np.random.shuffle(flds)
flds = flds[0:nflds]
norm = 1/nflds
grid = np.loadtxt(datadir+'geometry.out')
nz   = int(grid[0,2])
#
# single point statistics
#
for ifld in flds[:]:
    #
    fname = fname_beg+str(ifld).zfill(7)+fname_ext
    print(fname)
    data = np.loadtxt(fname)
    if ifld == flds[0]:
        data_avg = data
    else:
        data_avg += data
data_avg = data_avg*norm
data_avg[:,0] = data[:,0]
data_avg[:,1] = data[:,1]
#
fname = resultsdir + 'stats-single-point-duct-' + casename + fname_ext
with open(fname, 'w') as file:
  file.write('ZONE I=128, J=128, DATAPACKING=POINT\n')
  np.savetxt(file, data_avg, fmt='%16.6e', delimiter='') 
  # (fortran, e16.7) equivalent to (python,16,6e)

# utau_s = ((u1[0]+uconv)/zc[0]*visc)**.5
# uc     = u1[nz//2-1]+ub # convective reference frame
# uu_max = np.max(u2[:])
# uw_max = np.max(uw[:])
# print("Friction Reynolds number obtained from the velocity gradient at the wall: ", utau_s*h/visc)
# print("Centreline velocity:                  ", uc)
# print("Maximum streamwise velocity variance: ", uu_max)
# print("Maximum Reynolds shear stresses:      ", uw_max)

# do average of four corners
# plot centerline profiles as in Pedro's paper