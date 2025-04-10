import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt('EC-RP_IMS-RP_IC-RP_data_for_models.csv', dtype=None, delimiter=',', names=True,
                     encoding="utf_8_sig")
print(data.dtype.names)

plt.figure(constrained_layout=True)
reporters = ['ICRP', 'ECRP', 'IMSRP']
labels = {'ICRP': 'IC substrate (IC-RP)', 'ECRP': 'EC substrate (EC-RP)', 'IMSRP': 'MOMP (IMS-RP)'}
for reporter in reporters:
    # apply min-max normalization to the data
    min_max_norm = (data[reporter] - min(data[reporter])) / (max(data[reporter]) - min(data[reporter]))
    if reporter == 'IMSRP':
        plt.plot(data['Time'] / 3600, min_max_norm, 'o', ms=6, label=labels[reporter])
    else:
        plt.errorbar(data['Time'] / 3600, min_max_norm, yerr=np.sqrt(data['nrm_var_%s' % reporter] / 50), fmt='o', ms=6, capsize=4,
                     label=labels[reporter])

plt.xlabel('Time (hr)')
plt.ylabel('Signal')
plt.xlim(right=6.25)
plt.legend(loc='best')

plt.show()
