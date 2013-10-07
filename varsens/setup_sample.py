# Load libraries
import numpy as np
from   earm.lopez_embedded import model
from   pysb.integrate import odesolve
from   varsens import *
import time

import os, errno

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

mkdir_p("/scratch/garbetsp/varsens/samples")

# Why the first lsoda fails on ACCRE is a mystery
t = np.linspace(0, 20000, 10001)
try:
    odesolve(model, t, integrator='lsoda', nsteps=5000)
except:
    pass

# Generate Reference
t0=time.time()
reference = odesolve(model, t, integrator='lsoda', nsteps=5000)
exec_t = time.time() - t0


keys = model.parameters.keys()
k = len(keys) # Dimensions (i.e. parameters)
n = 5000              # Number of samples

# Just how bad is it going to be!
print "Evaluations required is %d\n" % (n*(2*k+2))
print "Hours estimated is %f\n" % (n*(2.0*k+2)*exec_t / 60.0 / 60.0)

# Given a set of parameters, determine outcome
def objective(x):
    for i in range(len(x)):
        model.parameters[keys[i]].value = x[i]
    x = odesolve(model, t, integrator='lsoda', nsteps=5000)
    return [np.sum((reference['mBid' ] - x['mBid' ])**2/10000),
            np.sum((reference['aSmac'] - x['aSmac'])**2/10000),
            np.sum((reference['cPARP'] - x['cPARP'])**2/10000)]

# Magnitude scaling a [0..1] set of points into the +/- 3 orders of magnitude
def scaling(points): return scale.magnitude(points, np.array([model.parameters[k].value for k in keys]))

sample = Sample(k, n, scaling)
sample.export("/scratch/garbetsp/varsens/samples/earm-batch-", ".csv", 11636)

