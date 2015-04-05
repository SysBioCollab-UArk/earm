# -*- coding: utf-8 -*-
"""
Created on Fri Jul 11 17:27:41 2014

@author: pinojc
"""

#!/usr/bin/env python

from varsens import *
from earm.lopez_embedded import model
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
param_dict = load_params("EARM_pso_fitted_params.txt")
for param in model.parameters:
    if param.name in param_dict:
        param.value = param_dict[param.name]


tspan = np.linspace(0, 20000, 1000)

solver = Solver(model, tspan, integrator='vode', rtol=1e-6, atol=1e-6,)
proteins_of_interest = ['Bid_0','Bax_0','Bak_0','Bcl2_0','BclxL_0','Mcl1_0','Bad_0']
key_list = []
initial_vals = dict()
for cp, value_obj in  model.initial_conditions:
    for j in proteins_of_interest:
        if  value_obj.name == j:
            initial_vals[value_obj.name] = value_obj.value
initial_vals = np.asarray(initial_vals.values())
initial_tmp = dict()
counter = 0
for cp, value_obj in  model.initial_conditions:
    for j in proteins_of_interest:
        if  value_obj.name == j:
            initial_tmp[value_obj.name] = 0#position[position]
            counter += 1


def likelihood1(position):
    initial_tmp = dict()
    counter = 0
    for cp, value_obj in  model.initial_conditions:
        for j in proteins_of_interest:
            if  value_obj.name == j:
                initial_tmp[value_obj.name] = position[counter]
                counter += 1
    solver.run(initial_changes=initial_tmp)
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
    return td

time_of_death = likelihood1(initial_vals)

def likelihood(position):
    initial_tmp = dict()
    counter = 0
    for cp, value_obj in  model.initial_conditions:
        for j in proteins_of_interest:
            if  value_obj.name == j:
                initial_tmp[value_obj.name] = position[counter]
                counter += 1
    solver.run(initial_changes=initial_tmp)
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
    return time_of_death - td
# 
# print initial_tmp
# solver.run(initial_changes=initial_tmp)
# print likelihood(initial_vals)
# quit()
   
def init(sample,dictionary):
    global Sample
    global Dictionary
    Sample,Dictionary= sample,dictionary

def OBJ(block):
    obj_values[block]=likelihood(sample[block])

if "__main__":# main():
    import time
    start = time.time()
    n_samples = 100
    sample = Sample(len(proteins_of_interest), n_samples, lambda x: scale.linear(x, lower_bound=0.5*initial_vals, upper_bound=5*initial_vals), verbose=True)
    sample = sample.flat()
    m = mp.Manager()
    obj_values = m.dict()
    p = mp.Pool(8,initializer = init, initargs=(sample,obj_values))
    allblocks =range(len(sample))
    p.imap_unordered(OBJ,allblocks)
    p.close()
    p.join()
    #print obj_values
    obj_vals = np.asarray(obj_values).reshape(len(obj_values),1)
    objective = Objective(len(proteins_of_interest), n_samples, objective_vals=obj_vals)
    v = Varsens(objective,verbose=True)
    print 'time of '+str(np.shape(sample)[0])+' calculations '+ str(time.time() - start)
    np.savetxt('sens_protein_concentration_earm.txt',v.sens)
    np.savetxt('senst_protein_concentration_earm.txt',v.sens_t)
    np.savetxt('sens2_protein_concentration_earm.txt',v.sens_2n[:,0,:,0])
    plt.plot(v.sens)
    plt.xticks(range(len(proteins_of_interest)),proteins_of_interest)
    plt.savefig('sens.png')
    plt.show()
    plt.plot(v.sens_t)
    plt.xticks(range(len(proteins_of_interest)),proteins_of_interest)
    plt.savefig('senst.png')
    plt.show()
    plt.plot(v.sens)
    plt.plot(v.sens_t)
    plt.show()
    plt.imshow(v.sens_2n[:,0,:,0],origin='lower',interpolation='nearest')
    plt.colorbar()
    plt.xticks(range(len(proteins_of_interest)),proteins_of_interest)
    plt.yticks(range(len(proteins_of_interest)),proteins_of_interest)
    plt.savefig('sens2.png')
    plt.show()
   


