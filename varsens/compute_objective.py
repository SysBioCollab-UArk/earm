# Load libraries
import numpy as np
from   earm.lopez_embedded import model
from   pysb.integrate import odesolve
from   varsens import *
import time
import sys
import os
import resource
import gc


# Generate Reference
t = np.linspace(0, 20000, 10001)

# Why the first lsoda fails on ACCRE is a mystery
try:
    odesolve(model, t, integrator='lsoda', nsteps=5000)
except:
    pass

reference = odesolve(model, t, integrator='lsoda', nsteps=5000)
allkeys = model.parameters.keys()
keys = [k for k in allkeys if not k.endswith("_0")]

# Given a set of parameters, determine outcome
def objective(x):
    for i in range(len(x)):
        model.parameters[keys[i]].value = x[i]
    x = odesolve(model, t, integrator='lsoda', nsteps=5000)
    return [np.sum((reference['mBid' ] - x['mBid' ])**2/10000/1e9),
            np.sum((reference['aSmac'] - x['aSmac'])**2/10000/1e9),
            np.sum((reference['cPARP'] - x['cPARP'])**2/10000/1e9)]


# Grab the sample space file
samplefile            = sys.argv[1]
(sampledir,filename)  = os.path.split(samplefile)
outfilename           = os.path.join(sampledir, "objective-"+(filename.split("-")[-1]))

print("Loading samples")
samples = np.loadtxt(samplefile)
print("Computing results")
result = []
for sample in samples:
    result.append(objective(sample))
    print resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    # Manage memory since Python gets it wrong
    if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss > 750000:
        gc.collect()

#result = np.array([objective(sample) for sample in samples])
print("Saving results")
np.savetxt(outfilename, np.array(result))
print("Output Saved (%s)" % outfilename)
