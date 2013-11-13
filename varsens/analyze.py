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
o.load("/Users/garbetsp/Projects/earm/varsens/samples/objective-", ".csv", 184, scaling=1e5)

v = Varsens(o, sample=s)

bid   = sorted([ (keys[i], v.sens[i][0]) for i in range(0,106) ], key=lambda x: x[1], reverse=True) 
aSmac = sorted([ (keys[i], v.sens[i][1]) for i in range(0,106) ], key=lambda x: x[1], reverse=True) 
cPARP = sorted([ (keys[i], v.sens[i][2]) for i in range(0,106) ], key=lambda x: x[1], reverse=True) 