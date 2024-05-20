import numpy as np
import sys
import os
import glob

sys.path.append('') # input.py in the working directory, rather than the directory of this script

from input import *

nt = 0
pattern = 'channel.???????.stats.txt'
with open('forcing.out', 'w') as fileout:
  for filename in sorted(glob.glob(pattern)):
    # glob.glob has the same order as ls -U, unordered by default
    with open(filename, 'r') as filein:
      first_line = filein.readline().split()
      fileout.write(' '.join(first_line[1:]) + '\n')
      t = float(first_line[1])
      if t > tbeg and t < tend:
        nt = nt + 1
        if(nt == 1):
          prop = np.array(first_line[1:], dtype=float)
          stat = np.loadtxt(filename, skiprows=1)
        else:
          prop = prop + np.array(first_line[1:], dtype=float)
          stat = stat + np.loadtxt(filename, skiprows=1)

prop = prop / nt
stat = stat / nt

retau = prop[1]
utau  = prop[2]
nu    = prop[3]
delv  = nu/utau

os.makedirs('results/',exist_ok=True)
np.savetxt('results/stats.txt', np.column_stack((retau, utau, delv)), fmt='%15.8f %15.8f %15.8f')
np.savetxt('results/stats-single-point-chan-'+casename+'.out', stat, fmt='%15.8f', delimiter='')
# default delimiter is a blank space
# delimiter='' must be specified to have the same format as F15.8 in Fortran