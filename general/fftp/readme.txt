n_order		12                   	! 2**12 moments
n_window        3								! 
n_startline	1										! start line number, at least is 1. end line number = start line number + 2**14 moments
n_points	48										! 48 points on one line
windowtype	1										! window type, default
filename	Probe.P								! filename
SoundSpeed      340.0						! dimensional
MachNumber      0.17						! dimensional
Density         1.225						! dimensional
timestep        0.0014					! dimensional
length          0.4572					! dimensional




# a long signal is divided into (n_window*2 - 1) section
# average only occurs on the time dimension, without on spanwise dimension
# output includs n_points	48 curves


# attention
# 3 sections + 2 overlaps
# total points < 2**(n_order-1)*(n_window*2 - 1) = 2**11 * 5


# how to define n_window
# if n_window is too large, the maximum period will be vert small, the minimum frequency will be very large
# n_window depend on practical condition and the time period which is already calculated


# the influence of windowtypecan be neglected 