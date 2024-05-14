import numpy as np
import sys
import os
import glob

sys.path.append('')
from input import *

nt = 0
pattern = 'channel.???????.stats.txt'
with open('forcing.out', 'w') as fileout:
  for filename in glob.glob(pattern):
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

np.savetxt('stats.txt', np.column_stack((retau, utau, delv)), fmt='%15.8f %15.8f %15.8f')