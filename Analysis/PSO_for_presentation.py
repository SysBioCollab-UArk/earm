# -*- coding: utf-8 -*-
"""
Created on Mon Mar 31 16:25:43 2014

@author: james
"""

import operator
import random
import scipy.optimize
import numpy
import pylab
from deap import base
from deap import creator
from deap import tools
import pysb.integrate
import pysb.util
import numpy as np
import scipy.optimize
import scipy.interpolate
import matplotlib.pyplot as plt
import os
import inspect
import multiprocessing
import multiprocessing as mp
from earm.lopez_embedded import model
#from earm.lopez_direct import model
#from earm.lopez_indirect import model
import sys
from pysb.util import load_params



obs_names = ['mBid', 'cPARP']
data_names = ['norm_ICRP', 'norm_ECRP']
var_names = ['nrm_var_ICRP', 'nrm_var_ECRP']
# Total starting amounts of proteins in obs_names, for normalizing simulations
obs_totals = [model.parameters['Bid_0'].value,
              model.parameters['PARP_0'].value]

earm_path = '/home/pinojc/Projects/earm'
data_path = os.path.join(earm_path, 'xpdata', 'forfits',
                         'EC-RP_IMS-RP_IC-RP_data_for_models.csv')
exp_data = np.genfromtxt(data_path, delimiter=',', names=True)

# Model observable corresponding to the IMS-RP reporter (MOMP timing)
momp_obs = 'aSmac'
# Mean and variance of Td (delay time) and Ts (switching time) of MOMP, and
# yfinal (the last value of the IMS-RP trajectory)
momp_obs_total = model.parameters['Smac_0'].value
momp_data = np.array([9810.0, 180.0, momp_obs_total])
momp_var = np.array([7245000.0, 3600.0, 1e4])



ntimes = len(exp_data['Time'])
tmul = 10
tspan = np.linspace(exp_data['Time'][0], exp_data['Time'][-1],
                    (ntimes-1) * tmul + 1)

#solver = pysb.integrate.Solver(model, tspan, integrator='lsoda', rtol=1e-6, atol=1e-6, nsteps=20000)
solver = pysb.integrate.Solver(model, tspan, integrator='vode', rtol=1e-6, atol=1e-6,)


rate_params = model.parameters_rules()
param_values = np.array([p.value for p in model.parameters])
rate_mask = np.array([p in rate_params for p in model.parameters])
k_ids = [p.value for p in model.parameters_rules()]



def likelihood(position):
    Y=np.copy(position)
    param_values[rate_mask] = 10 ** Y
    solver.run(param_values)

    
    for obs_name, data_name, var_name, obs_total in \
            zip(obs_names, data_names, var_names, obs_totals):
        ysim = solver.yobs[obs_name][::tmul]
        ysim_norm = ysim / obs_total
        ydata = exp_data[data_name]
        yvar = exp_data[var_name]
        if obs_name == 'mBid':
            e1 = np.sum((ydata - ysim_norm) ** 2 / (2 * yvar)) / len(ydata)
        else:
            e2 = np.sum((ydata - ysim_norm) ** 2 / (2 * yvar)) / len(ydata)

    ysim_momp = solver.yobs[momp_obs]
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

    error = e1 + e2 
    return [error, e3]

def display(position):


    exp_obs_norm = exp_data[data_names].view(float).reshape(len(exp_data), -1).T
    var_norm = exp_data[var_names].view(float).reshape(len(exp_data), -1).T
    std_norm = var_norm ** 0.5
    Y=np.copy(position)
    param_values[rate_mask] = 10 ** Y
    solver.run(param_values)
    obs_names_disp = obs_names + ['aSmac']
    sim_obs = solver.yobs[obs_names_disp].view(float).reshape(len(solver.yobs), -1)
    totals = obs_totals + [momp_obs_total]
    sim_obs_norm = (sim_obs / totals).T
    colors = ('r', 'b')
    for exp, exp_err, sim, c in zip(exp_obs_norm, std_norm, sim_obs_norm, colors):
        plt.plot(exp_data['Time'], exp, color=c, marker='.', linestyle=':')
        plt.errorbar(exp_data['Time'], exp, yerr=exp_err, ecolor=c,
                     elinewidth=0.5, capsize=0, fmt=None)
        plt.plot(solver.tspan, sim, color=c)
    plt.plot(solver.tspan, sim_obs_norm[2], color='g')
    plt.vlines(momp_data[0], -0.05, 1.05, color='g', linestyle=':')
    plt.show()



    
def generate(size, speedmin, speedmax):
    u1 = np.random.uniform(-1,1,size)
    tmp = np.log10(k_ids) + u1
    part = creator.Particle(tmp)
    part.speed = np.random.uniform(speedmin,speedmax,size)
    part.smin = speedmin
    part.smax = speedmax
    return part

def updateParticle(part, best, phi1, phi2):


    u1 = numpy.random.uniform(0, phi1, len(part))
    u2 = numpy.random.uniform(0, phi2, len(part))
    v_u1 = u1 * (part.best - part)
    v_u2 = u2 * (best - part)
    part.speed += v_u1 + v_u2
    
    for i, speed in enumerate(part.speed):
        if speed < part.smin:
            part.speed[i] = part.smin
        elif speed > part.smax:
            part.speed[i] =  part.smax
    part += part.speed
    for i, pos in enumerate(part):
        if pos < lb[i]:
            part[i] = lb[i]
        elif pos > ub[i]:
            part[i] =  ub[i]
def dominates(new, old):
    sum_old = np.sum(old)
    sum_new = np.sum(new)
    diff = sum_old - sum_new
    change1= old[0] - new[0]
    change2= old[1] - new[1] 
    #change3 = old[2] - new[2]
    if  change1 >= 0. and change2 >= 0. :#and change3 >=0.:
        return True
    if change1 >= 0. or change2 >=0. :#or change3 >=0.: 
        if diff >= sum_old*.1:
            return True
        else:
            return False
    else:
        return False                
toolbox = base.Toolbox()
nominal_values = np.array([p.value for p in model.parameters])
xnominal = np.log10(nominal_values[rate_mask])
bounds_radius = 3
lb = xnominal - bounds_radius
ub = xnominal + bounds_radius


#creator.create("FitnessMin", base.Fitness,weights=(-1.00,-1.00,-1.00))
creator.create("FitnessMin", base.Fitness,weights=(-1.00,-1.00))
creator.create("Particle", np.ndarray, fitness=creator.FitnessMin, \
    speed=list,smin=list, smax=list, best=None)
toolbox.register("particle", generate, size=np.shape(rate_params)[0],\
                 speedmin=-.2,speedmax=.2)
toolbox.register("population", tools.initRepeat, list, toolbox.particle)
toolbox.register("update", updateParticle, phi1=2, phi2=2)
toolbox.register("evaluate", likelihood)


stats = tools.Statistics(lambda ind: ind.fitness.values)
stats.register("avg", numpy.mean, axis=0)
stats.register("std", numpy.std, axis=0)
stats.register("min", numpy.min, axis=0)
stats.register("max", numpy.max, axis=0)
logbook = tools.Logbook()
logbook.header = ["gen", "evals"] + stats.fields

def init(sample,dictionary):
    global Sample
    global Dictionary
    Sample,Dictionary = sample,dictionary

def OBJ(block):
    obj_values[block]=likelihood(sample[block])

if "__main__":# main():
    GEN = 1000
    num_particles = 10
    
    best = creator.Particle(xnominal)
    best.fitness.values = likelihood(xnominal)
    pop = toolbox.population(n=num_particles)
    best_values =np.zeros((GEN,2))
    evals = np.zeros(GEN)
    for g in range(1,GEN+1):
        m = mp.Manager()
        obj_values = m.dict()
        sample = []
        for p in pop:
            sample.append(p)
        p = mp.Pool(4,initializer = init, initargs=(sample,obj_values))
        allblocks =range(len(pop))
        p.imap_unordered(OBJ,allblocks)
        p.close()
        p.join()
        count=0

        for part in pop:
            part.fitness.values = obj_values[count]
            
            count+=1
            if not g == 1:
                #if part.fitness.values < part.best.fitness.values:
                if dominates(part.fitness.values,part.best.fitness.values):
                    part.best = creator.Particle(part)
                    part.best.fitness.values = part.fitness.values
            else:
                part.best = creator.Particle(part)
                part.best.fitness.values = part.fitness.values
            if dominates(part.fitness.values,best.fitness.values):
            #if part.fitness.values < best.fitness.values:
                best = creator.Particle(part)
                best.fitness.values = part.fitness.values
        for part in pop:
            toolbox.update(part, best)
 
        logbook.record(gen=g, evals=len(pop), **stats.compile(pop))
        print(logbook.stream),best.fitness.values
        best_values[g-1,:] = best.fitness.values
        evals[g-1] = g*num_particles
    
    plt.plot(evals,best_values[:,0])
    plt.show()
    plt.loglog(best_values[:,0],best_values[:,1])
    plt.show()
    np.savetxt('multi_pso_all_moves.txt',best_values)
    #np.savetxt('bayessb_kept_moves.txt',mcmc.likelihoods[mcmc.accepts])
    #quit()
    #display(best)
    
        




