from reformat_EARM_data import output_file
from earm.lopez_embedded import model
from param_calibration import *
from pysb.simulator import ScipyOdeSimulator

exp_data_file = output_file

params_file = os.path.join('..', '..', 'EARM_2_0_M1a_fitted_params.txt')
params = np.genfromtxt(params_file, dtype=None, delimiter=',', encoding="utf_8_sig")
for p_name, p_val in params:
    model.parameters[p_name].value = p_val

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
