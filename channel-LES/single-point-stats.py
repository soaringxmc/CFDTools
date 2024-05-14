import numpy as np
import sys
import os

sys.path.append('')
from input import *
import glob


with open('input.nml', 'r') as file:
  lines = file.readlines()
  npoints = int(lines[1].split()[1])
prop = np.zeros(8)
stats = np.zeros([npoints, 10])

nt = 0
pattern = 'channel.???????.stats.txt'
for filename in glob.glob(pattern):
  print(filename)
  with open(filename, 'r') as file:
    data = file.readline().split()
    t = float(data[1])
    if t > tbeg and t < tend:
      nt = nt + 1
      prop = prop + np.array(data[1:], dtype=float)
      stats = stats + np.loadtxt(filename, skiprows=1)
prop = prop / nt
stats = stats / nt

retau = prop[1]
utau  = prop[2]
nu    = prop[3]
delv  = nu/utau

np.savetxt('stats.txt', [retau,utau,delv], fmt='%15.8f', delimiter='')