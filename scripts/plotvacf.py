#!/usr/bin/python

import sys
from numpy import genfromtxt
from matplotlib import pyplot as plt

if len(sys.argv) < 2: sys.exit(1)

if len(sys.argv) >= 3: 
    scale = float(sys.argv[2])
else:
    scale = 0.005

data = genfromtxt(sys.argv[1], unpack=True)
xv = data[0]
yv = data[1]
plt.plot([x * 0.005 for x in xv], yv)
plt.xlabel('Time [fs]')
plt.ylabel('VAF')
plt.savefig(sys.argv[1]+'.pdf')

