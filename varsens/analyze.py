from earm.lopez_embedded import model
import numpy as np
from varsens             import *

allkeys = model.parameters.keys()
keys = [k for k in allkeys if not k.endswith("_0")]
#keys = allkeys
k = len(keys)
n = 10000

s = Sample(k, n, None)
o = Objective(k, n, s, None)
o.load("/Users/garbetsp/Projects/earm/varsens/samples/objective-", ".csv", 184, scaling=1e-1)

o.fM_1  = np.log(o.fM_1)
o.fM_2  = np.log(o.fM_2)
o.fN_j  = np.log(o.fN_j)
o.fN_nj = np.log(o.fN_nj)

v = Varsens(o, sample=s)

v.E_2

bid   = sorted([ (keys[i], v.sens[i][0]) for i in range(0,106) ], key=lambda x: x[1], reverse=True) 
aSmac = sorted([ (keys[i], v.sens[i][1]) for i in range(0,106) ], key=lambda x: x[1], reverse=True) 
cPARP = sorted([ (keys[i], v.sens[i][2]) for i in range(0,106) ], key=lambda x: x[1], reverse=True) 