import operator
import random
import scipy.optimize
import numpy
import pylab
import pysb.integrate
import pysb.util
import numpy as np

import matplotlib.pyplot as plt


from pysb.util import load_params
from pysb.bng import run_ssa
from pysb.integrate import odesolve

t = np.linspace(0,20000,1000)
def run(model,name):

    x = odesolve(model,t)
    #print model.odes
    for i,j in enumerate(model.odes):
        print i,j
#     for j,obs in enumerate(model.observables):
#         #plt.figure(str(obs))
#         plt.subplot(len(model.observables),1,j)
#         plt.title(str(obs))
#         plt.plot(t[:len(x[str(obs)])],x[str(obs)],lw=3)#label=name)
#         #plt.plot(t[:len(y[str(obs)])],y[str(obs)],lw=2,label=str(obs)+name)
#         plt.legend(loc=0)
#         plt.xlim(xmax=t[-1])
#         plt.xlabel("Time")
#         plt.ylabel("Concentration")
#         plt.ticklabel_format(useOffset=False)
#         plt.ticklabel_format(axis='x', style='sci', scilimits=(0,0))

def run2(model,name):
    x = run_ssa(model,t)
    #x = odesolve(model,t)
    #plt.figure(name)
    for j,obs in enumerate(model.observables):
        plt.title(obs)
        plt.subplot(len(model.observables),1,j)
        plt.plot(t[:len(x[str(obs)])],x[str(obs)],lw=3,label=str(name))
        plt.xlim(xmax=t[-1])
        plt.ticklabel_format(useOffset=False)
        plt.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
    plt.xlabel("Time")
    plt.ylabel("Concentration")
    plt.legend(loc=0)


from earm.lopez_direct import model as model1
from earm.lopez_indirect import model as model2
from earm.lopez_embedded import model as model3

#run(model1,' direct')
#run(model2,' indirect')
#run(model3,' embedded')
#plt.show()
def sample(model):
    #for i in ('Bid_0','BclxL_0','Mcl1_0','Bcl2_0','Bad_0','Noxa_0','CytoC_0','Smac_0','Bax_0','Bak_0'):
    for i in ('Noxa_0','Mcl1_0'):
        for j in (1e4,1e5,1e6,1e10):#np.linspace(0,2,5):
            tmp = np.copy(model.parameters[i].value)
            model.parameters[i].value = j
            run2(model,j)
            model.parameters[i].value = tmp
        plt.show()

#x = run_ssa(model1,t)
#x = odesolve(model3,t)
#for i in range(len(model3.species)):
#    plt.title(model3.species[i])
#    plt.plot(t,x['__s'+str(i)])
#    plt.show()
for i in xrange(1):
    run(model3,'h')