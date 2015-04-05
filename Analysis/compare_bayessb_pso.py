import pylab as plt
import numpy as np


pso = np.loadtxt('/home/pinojc/Projects/earm/Analysis/pso_all_moves.txt')
bay = np.loadtxt('/home/pinojc/git/moo_and_mcmc/earm/bayessb_all_moves.txt')
print np.shape(pso)
print np.shape(bay)
plt.semilogy(pso,label='PSO')
plt.semilogy(bay,label='MC')
plt.legend(loc=0)
plt.xlabel('Iteration Number')
plt.ylabel('Log(cost function)')
plt.savefig('pso_vs_MC.png',dpi=300)
plt.show()