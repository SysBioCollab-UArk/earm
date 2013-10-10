#from earm.lopez_embedded import model
import numpy
from varsens             import *

k = 127 #len(model.parameters.keys())
n = 5000

s = Sample(k, n, None)
o = Objective(k, n, s, None)
o.load("/Users/garbetsp/Projects/earm/varsens/samples/objective-", ".csv", 111, scaling=1e5)

v = Varsens(o, sample=s)