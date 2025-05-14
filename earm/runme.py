#from earm.albeck_11e import model
from earm.howells import model
import numpy as np
import matplotlib.pyplot as plt
from pysb.integrate import odesolve

print(model.monomers)

from pysb import Observable
from pysb.util import alias_model_components
alias_model_components()
Observable('Bid_', Bid(state='T') + Bid(state='M'))
Observable('mSmac_', Smac(state='C') + Smac(state='A'))
Observable('cPARP_', PARP(state='C'))

print(model.observables)

t = np.linspace(0, 20000, 10000)
x = odesolve(model, t)
Smac_0 = model.parameters['Smac_0'].value
Bid_0 = model.parameters['Bid_0'].value
PARP_0 = model.parameters['PARP_0'].value
plt.plot(t, 1 - (x['Bid_']/Bid_0), label='t/mBid')
plt.plot(t, 1 - (x['mSmac_']/Smac_0), label='c/aSmac')
plt.plot(t, x['cPARP_']/PARP_0, label='cPARP')
plt.legend(loc='best')
plt.tight_layout()
plt.show()

