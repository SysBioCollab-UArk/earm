import numpy as np
import matplotlib.pyplot as plt
import os
from earm.lopez_embedded import model
from pysb.simulator import ScipyOdeSimulator

data = np.genfromtxt('forfits/EC-RP_IMS-RP_IC-RP_data_for_models.csv', dtype=None, delimiter=',', names=True,
                     encoding="utf_8_sig")
print(data.dtype.names)

plt.figure(constrained_layout=True)
reporters = ['ICRP', 'IMSRP', 'ECRP']
labels = {'ICRP': 'IC substrate (IC-RP)', 'ECRP': 'EC substrate (EC-RP)', 'IMSRP': 'MOMP (IMS-RP)'}
for reporter in reporters:
    # apply min-max normalization to the data
    min_max_norm = (data[reporter] - min(data[reporter])) / (max(data[reporter]) - min(data[reporter]))
    if reporter == 'IMSRP':
        yerr = [0.01] * len(data['VAR'])
    else:
        yerr = np.sqrt(data['nrm_var_%s' % reporter] / 50)  # 50 samples per data point
        print('avg:', np.mean(yerr))
    plt.errorbar(data['Time'] / 3600, min_max_norm, yerr=yerr, fmt='o', ms=4, capsize=4, label=labels[reporter])

plt.xlabel('Time (hr)')
plt.ylabel('Signal')
plt.xlim(right=6.25)
plt.legend(loc='best')

#############

params_file = os.path.join('..', 'EARM_2_0_M1a_fitted_params.txt')
params = np.genfromtxt(params_file, dtype=None, delimiter=',', encoding="utf_8_sig")
print(params)

tspan = np.linspace(0, 6*60*60, 1001)
sim = ScipyOdeSimulator(model, tspan, verbose=True)
param_values = {p_name: p_val for p_name, p_val in params}
output = sim.run(param_values=param_values)

plt.figure(constrained_layout=True)
for obs in model.observables:
    plt.plot(tspan / 3600, output.observables[obs.name] / output.observables[obs.name][-1], lw=2, label=obs.name)
plt.xlabel('Time (hr)')
plt.ylabel('Signal')
plt.legend(loc='best')

plt.show()
