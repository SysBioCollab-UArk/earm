from reformat_EARM_data import output_file
from earm.lopez_embedded import model
from param_calibration import *
from pysb.simulator import ScipyOdeSimulator

exp_data_file = output_file
solver = ScipyOdeSimulator(model)
sim_protocol = SimulationProtocol(solver)

custom_priors = {}
no_sample = []
obs_labels = {'mBid': 'mBid (IC-RP)', 'aSmac': 'aSmac (IMS-RP)', 'cPARP': 'cPARP (EC-RP)'}

if __name__ == '__main__':

    calibrator = ParameterCalibration(model,
                                      exp_data_file,
                                      sim_protocol,
                                      priors=custom_priors,
                                      no_sample=no_sample)

    calibrator.run(niterations=50000, nchains=5, obs_labels=obs_labels, plot_results=True,
                   plot_tc_args={'separate_plots': False, 'save_sim_data': True})
