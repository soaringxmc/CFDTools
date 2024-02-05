#!/usr/bin/python
# values should be consistent with dns.in
h     = 1.
ub    = 1.
visci = 39600.
#
uconv = 0. # if we solve on a convective reference frame; else = 0.
#
# parameters for averaging
#
tbeg   = 15.
tend   = 1.e9
fldstp = 100
#
# case name (e.g., the Retau)
#
casename = '1000'.zfill(5)
