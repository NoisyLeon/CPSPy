import numpy as np

inArr = np.loadtxt('SDISPR.ASC')

T = inArr[:,2]
V = inArr[:,4]
T10 = np.interp(10., T, V, left=None, right=None, period=None)