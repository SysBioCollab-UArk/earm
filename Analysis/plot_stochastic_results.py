import pylab as plt
import numpy as np
import os

from pysb.util import load_params
from pysb.bng2 import run_ssa
from earm.lopez_embedded import model
import scipy.interpolate

# param_dict = load_params("EARM_pso_fitted_params.txt")
# for param in model.parameters:
#     if param.name in param_dict:
#         param.value = param_dict[param.name]
param = np.loadtxt('56_embedded_single_overall_best.txt')
rate_params = model.parameters_rules()
param_values = np.array([p.value for p in model.parameters])
rate_mask = np.array([p in rate_params for p in model.parameters])
k_ids = [p.value for p in model.parameters_rules()]
Y=np.copy(param)
param_values[rate_mask] = 10 ** Y
tspan = np.linspace(0, 20000, 1000)
#  
# 
proteins_of_interest = ['Bid_0','Bax_0','Bak_0','Bcl2_0','BclxL_0','Mcl1_0','Bad_0']
key_list = []
initial_tmp = dict()
for cp, value_obj in  model.initial_conditions:
    for j in proteins_of_interest:
        if  value_obj.name == j:
            initial_tmp[value_obj.name] = value_obj.value

 
def runner(changes,n_runs,directory):
    if not os.path.exists(directory):
        os.mkdir(directory)
    run_ssa(model,tspan,param_values,initial_changes=changes,n_runs=n_runs,output_dir=directory,verbose=False)
#changes = dict()
#changes['L_0'] = 3000*5
#runner(changes,1000,'Fit_250ns_ml_TRAIL') 

#changes = dict()
#changes['L_0'] = 3000/25
#runner(changes,1000,'Fit_2ng_ml_TRAIL') 
#quit()
#changes = dict()
#changes['Mcl1_0'] = initial_tmp['Mcl1_0']*.5
#runner(changes,1000,'Fit_Mcl1_timesonehalf') 
#quit()
#changes = dict()
#changes['Bad_0'] = initial_tmp['Bad_0']*2
#runner(changes,500,'Fit_Bad_times2')
#quit()
#changes = dict()
#changes['Bcl2_0'] = initial_tmp['Bcl2_0']*2
#runner(changes,1000,'Fit_Bcl2')
#quit()
# changes = dict()
# changes['Bax_0'] = initial_tmp['Bax_0']*2
# runner(changes,500,'Fit_values_500_runs_2timesBax')
# quit()
# changes = dict()
# changes['Bak_0'] = initial_tmp['Bak_0']*2
# runner(changes,500,'Fit_values_500_runs_2timesBak')
# changes = dict()
# changes['BclxL_0'] = initial_tmp['BclxL_0']*2
# runner(changes,500,'Fit_values_500_runs_2timesBclXL')
#quit()
# vals = [.1,1.1]
# done = []
#vals = np.logspace(0,1,25)
#vals = np.hstack((np.linspace(0,.9,9),np.logspace(0,1,10)))
#print vals
#quit()
# done = []
# vals = [.5,2.]
# for i  in proteins_of_interest:
#     x = initial_tmp[i]
#     for j in proteins_of_interest:
#         if i == j:
#             continue
#         if i+j in done:
#             continue
#         done.append(j+i)
#         samples = []
#         y = initial_tmp[j]
#         for a,c in enumerate(vals):
#            for b,d in enumerate(vals):
#                 tmp = dict()
#                 tmp[i] = x*c
#                 tmp[j] = y*d
#                 runner(tmp,1000,"Fit_%s_%s_%s_%s" % (i,str(c),j,str(d)))

#quit()
tbid = None
directory = '/home/pinojc/Projects/earm/Analysis/Fit_2ng_ml_TRAIL/'
for i in range(0,1000):
    FILE = directory+str(i)+'.gdat'
    if os.path.exists(FILE):
        #print i
        tmp = np.loadtxt(FILE,skiprows=1)
        if tbid == None:
            time = tmp[:,0]
            tbid = tmp[:,1]
            smac = tmp[:,2]
            cparp = tmp[:,3]
        else:
            tbid = np.column_stack((tbid,tmp[:,1]))
            smac = np.column_stack((smac,tmp[:,2]))
            cparp = np.column_stack((cparp,tmp[:,3]))
plt.plot(time,tbid)
plt.title('tbid IC-RP')
#plt.errorbar(time,np.average(tbid, axis=1),xerr=np.std(tbid, axis=1))
plt.show()

plt.plot(time,cparp)
plt.title('cparp')
#plt.errorbar(time,np.average(cparp, axis=1),xerr=np.std(cparp, axis=1))
plt.show()

plt.plot(time,smac)
plt.title('smac')
#plt.errorbar(time,np.average(smac, axis=1),xerr=np.std(smac, axis=1))
plt.show()
# plt.plot(tmp[:,0],tbid[:,:50]/np.max(tbid[:,:10]),'r')
# plt.plot(tmp[:,0],cparp[:,:50]/np.max(cparp[:,10]),'g')
# plt.plot(tmp[:,0],smac[:,:50]/np.max(smac[:,10]),'orange')
# plt.show()
global nodeath
nodeath = 0.
def likelihood(ysim_momp):
    global nodeath
    if np.nanmax(ysim_momp) == 0:
        print 'No aSmac!'
        ysim_momp_norm = ysim_momp
        t10 = 0
        t90 = 0
        
        nodeath +=1.
    else:
        ysim_momp_norm = ysim_momp / np.nanmax(ysim_momp)
        st, sc, sk = scipy.interpolate.splrep(tspan, ysim_momp_norm)
        try:
            t10 = scipy.interpolate.sproot((st, sc-0.10, sk))[0]
            t90 = scipy.interpolate.sproot((st, sc-0.90, sk))[0]
        except IndexError:
            t10 = 0
            t90 = 0

    td = (t10 + t90) / 2
    ts = t90 - t10
    return td/3600,ts/60
timeOfDeath = []
timeSwithching = []
for i in range(np.shape(cparp)[1]):
    tmp_tod,tmp_ts = likelihood(cparp[:,i])
    timeOfDeath.append(tmp_tod)
    timeSwithching.append(tmp_ts)
plt.hist(np.asarray(timeOfDeath),10)
plt.xlim(0,8)
plt.show()
plt.hist(np.asarray(timeSwithching),10)
plt.show()
print 'fraction of cells not dead = %s ' %str(nodeath/1000.*100.)
print 'time of death = %.4f , std of time of death %.4f' %(np.mean(timeOfDeath)*60,np.std(timeOfDeath)*60)
print 'time of switch = %.4f , std of time of switch %.4f' %(np.mean(timeSwithching),np.std(timeSwithching))
