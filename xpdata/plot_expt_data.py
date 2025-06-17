from earm.lopez_embedded import model
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib.lines import Line2D
from matplotlib.collections import LineCollection
from matplotlib.container import ErrorbarContainer

# Plot experimental data

data = np.genfromtxt('forfits/EC-RP_IMS-RP_IC-RP_data_for_models.csv', dtype=None, delimiter=',', names=True,
                     encoding="utf_8_sig")
print(data.dtype.names)

plt.figure(constrained_layout=True)
reporters = ['ICRP', 'IMSRP', 'ECRP']
# labels = {'ICRP': 'IC substrate (IC-RP)', 'ECRP': 'EC substrate (EC-RP)', 'IMSRP': 'MOMP (IMS-RP)'}
labels = {'ICRP': 'IC-RP', 'ECRP': 'EC-RP', 'IMSRP': 'IMS-RP'}
markers = ['o', 's', '^', 'v', 'D', '<', '>', 'p']
colors = []
for i, reporter in enumerate(reporters):
    # apply min-max normalization to the data
    min_max_norm = (data[reporter] - min(data[reporter])) / (max(data[reporter]) - min(data[reporter]))
    if reporter == 'IMSRP':
        yerr = [0.01] * len(data['VAR'])
    else:
        yerr = np.sqrt(data['nrm_var_%s' % reporter] / 50)  # 50 samples per data point
    marker = markers[i % len(markers)]
    p = plt.errorbar(data['Time'] / 3600, min_max_norm, yerr=yerr, fmt=marker, ms=5, mfc='none', capsize=4, alpha=0.4,
                 label=labels[reporter], zorder=1)
    colors.append(p[0].get_color())
plt.xlabel('Time (hr)')
plt.ylabel('Signal')
plt.xlim(right=6.25)
plt.legend(loc='best')

# Plot simulation time courses

params_file = os.path.join('..', 'EARM_2_0_M1a_fitted_params.txt')
params = np.genfromtxt(params_file, dtype=None, delimiter=',', encoding="utf_8_sig")

tspan = np.linspace(0, 6*60*60, 1001)
sim = ScipyOdeSimulator(model, tspan, verbose=True)
param_values = {p_name: p_val for p_name, p_val in params}
output = sim.run(param_values=param_values)

for obs, color in zip(model.observables, colors):
    plt.plot(tspan / 3600, output.observables[obs.name] / output.observables[obs.name][-1], lw=2, color=color,
             label=obs.name, zorder=2)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel('Time (hr)', fontsize=14)
plt.ylabel('Signal', fontsize=14)
plt.legend(loc='best')

# Combine entries in the legend
handles, labels = plt.gca().get_legend_handles_labels()
new_handles = []
for handle, eb in zip(handles, handles[len(reporters):]):
    # Clone marker (eb[0] is a Line2D)
    marker = Line2D([0], [0], marker=eb[0].get_marker(), ls='', mfc=eb[0].get_markerfacecolor(),
                    mec=eb[0].get_markeredgecolor(), ms=eb[0].get_markersize(), label=eb[0].get_label(), alpha=1.0)
    # Clone line (eb[1] is a tuple of Line2D)
    line = Line2D([0], [0], ls=eb[1][0].get_linestyle(), lw=eb[1][0].get_linewidth(),
                  color=eb[1][0].get_color(), alpha=1.0)
    # Clone error bars (eb[2][0] is a LineCollection)
    bar = LineCollection(eb[2][0].get_segments(), colors=eb[2][0].get_edgecolor(), linewidths=eb[2][0].get_linewidths(),
                         linestyles=eb[2][0].get_linestyle(), alpha=1.0)
    opaque_eb = ErrorbarContainer((marker, (line,), (bar,)), has_yerr=True)
    new_handles.append((handle, opaque_eb))
new_labels = ['%s: %s' % (labels[i], labels[i + len(reporters)]) for i in range(len(reporters))]
plt.legend(new_handles, new_labels, loc='lower right', fontsize=12)

plt.show()
