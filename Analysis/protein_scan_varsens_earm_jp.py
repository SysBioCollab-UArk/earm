# -*- coding: utf-8 -*-
"""
Created on Fri Jul 11 17:27:41 2014

@author: pinojc
"""

#!/usr/bin/env python

from varsens import *
from pysb.integrate import Solver, odesolve
from pysb.util import load_params
import os
import copy
import pickle
import numpy as np
import scipy.interpolate
import matplotlib.pyplot as plt
import datetime
from multiprocessing import Pool, Value, Queue
import multiprocessing
import multiprocessing as mp
import math
from pysb.util import load_params
from pysb.bng import run_ssa
from earm.lopez_embedded import model


#param_dict = load_params("EARM_pso_fitted_params.txt")
#for param in model.parameters:
#    if param.name in param_dict:
#        param.value = param_dict[param.name]


tspan = np.linspace(0, 20000, 1000)

solver = Solver(model, tspan, integrator='vode', rtol=1e-6, atol=1e-6,)
proteins_of_interest = ['Bid_0','Bax_0','Bak_0','Bcl2_0','BclxL_0','Mcl1_0','Bad_0']
key_list = []
initial_tmp = dict()
for cp, value_obj in  model.initial_conditions:
    for j in proteins_of_interest:
        if  value_obj.name == j:
            initial_tmp[value_obj.name] = value_obj.value
#initial_vals = np.asarray(initial_vals.values())
#initial_tmp = dict(initial_vals)

            
momp_obs_total = model.parameters['Smac_0'].value
momp_data = np.array([9810.0, 180.0, momp_obs_total])
momp_var = np.array([7245000.0, 3600.0, 1e4])


           
def likelihood1(initial_tmp):
    
    solver.run(initial_changes=initial_tmp)
    #y = run_ssa(model, tspan, initial_changes=initial_tmp)
    ysim_momp = solver.yobs['aSmac']
    if np.nanmax(ysim_momp) == 0:
        print 'No aSmac!'
        ysim_momp_norm = ysim_momp
        t10 = 0
        t90 = 0
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
    yfinal = ysim_momp[-1]
    momp_sim = [td, ts, yfinal]
    e3 = np.sum((momp_data - momp_sim) ** 2 / (2 * momp_var)) / 3
    return   td 
tod = likelihood1(initial_tmp)
#quit()
print tod
def likelihood(initial_tmp):
    #print initial_tmp
    #solver.run(initial_changes=initial_tmp)
    yobs = run_ssa(model, tspan, initial_changes=initial_tmp)
    #quit()
    #ysim_momp = solver.yobs['aSmac']
    ysim_momp = yobs['aSmac']
    if np.nanmax(ysim_momp) == 0:
        #print 'No aSmac!'
        ysim_momp_norm = ysim_momp
        t10 = 0
        t90 = 0
        return tod*3
    else:
        ysim_momp_norm = ysim_momp / np.nanmax(ysim_momp)
        st, sc, sk = scipy.interpolate.splrep(tspan, ysim_momp_norm)
        try:
            t10 = scipy.interpolate.sproot((st, sc-0.10, sk))[0]
            t90 = scipy.interpolate.sproot((st, sc-0.90, sk))[0]
        except IndexError:
            t10 = 0
            t90 = 0
    plt.plot(tspan,ysim_momp)
    td = (t10 + t90) / 2
    ts = t90 - t10
    yfinal = ysim_momp[-1]
    momp_sim = [td, ts, yfinal]
    e3 = np.sum((momp_data - momp_sim) ** 2 / (2 * momp_var)) / 3
    return  (td-tod)/tod 
    

   
def init(sample,dictionary):
    global Sample
    global Dictionary
    Sample,Dictionary= sample,dictionary

def OBJ(block):
    obj_values[block]=likelihood(sample[block])

if "__main__":# main():
    import time
    start = time.time()
    vals = [0,.1,.25,.5,.75,1.,1.25,]
    done = []
    #vals = np.logspace(0,1,25)
    #vals = np.hstack((np.linspace(0,.9,9),np.logspace(0,1,10)))
    #print vals
    #quit()
    for i  in proteins_of_interest:
        x = initial_tmp[i]
        for j in proteins_of_interest:
            if i == j:
                continue
            if i+j in done:
                continue
            done.append(j+i)
            samples = []
            y = initial_tmp[j]
            im = np.zeros((len(vals),len(vals)))
            for a,c in enumerate(vals):
                for b,d in enumerate(vals):
                    tmp = dict()
                    tmp[i] = x*c
                    tmp[j] = y*d
                    im[a,b] = likelihood(tmp)
            plt.show()
            #plt.imshow(im,interpolation='nearest',origin='lower',vmin=-2,vmax=2,cmap=plt.get_cmap('bwr'))
            #plt.xticks(range(len(vals)),np.round(vals,decimals=2),rotation='vertical')
            #plt.yticks(range(len(vals)),np.round(vals,decimals=2))
            #plt.xlabel(j)
            #plt.ylabel(i)
            #plt.colorbar()
            #plt.title("%s_%s" % (i,j))
            #plt.savefig("%s_lit_values_%s_%s.png" % (model.name,i,j))
            #plt.show()
            #plt.clf()  
            np.savetxt("stoch_%s_lit_values_%s_%s.txt" % (model.name,i,j),im)  


