#!/usr/bin/python
# values should be consistent with dns.in
h     = 0.
ub    = 0.
visci = 0.
#
uconv = 0. # if we solve on a convective reference frame; else = 0.
#
# parameters for averaging
#
tbeg   = 0.8
tend   = 1.8
fldstp = 100
#
# case name (e.g., the Retau)
#
casename = '180'.zfill(5)
