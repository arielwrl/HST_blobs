"""

ariel@oapd
15/06/2022

Creates a grid of models and fits them

"""

import bagpipes as pipes
import numpy as np
from astropy.table import Table
import itertools
from hst_pipeutils import load_data_models as load_data
import matplotlib.pyplot as plt

data_dir = '/Data/'
test_id = '5err_doublescreen_eta_vary_oldmass10'

filter_files = [data_dir + 'filters/HST_WFC3_UVIS2.F275W.dat',
                data_dir + 'filters/HST_WFC3_UVIS2.F336W.dat',
                data_dir + 'filters/HST_WFC3_UVIS2.F606W.dat',
                data_dir + 'filters/HST_WFC3_UVIS2.F680N.dat',
                data_dir + 'filters/HST_WFC3_UVIS2.F814W.dat']

plot_dir = '/BAGPIPES/plots/Models/'
table_rows = []
model_dict = {}
model_ages = [0.005, 0.015, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25]
model_Avs = [0.0, 0.1, 0.2, 0.4, 0.6]

combinations = [(model_age, model_Av) for model_age, model_Av in itertools.product(model_ages, model_Avs)]

for i in range(len(combinations)):

    model_age, model_Av = combinations[i]
    print(i+1, model_age, model_Av)

    exp_young = {}
    exp_young["age"] = model_age
    exp_young["tau"] = 2.5
    exp_young["massformed"] = 5
    exp_young["metallicity"] = 0.5

    exp_old = {}
    exp_old["age"] = 10
    exp_old["tau"] = 3
    exp_old["massformed"] = 10
    exp_old["metallicity"] = 0.5

    dust = {}
    dust["type"] = "Cardelli"
    dust["Av"] = model_Av
    dust["eta"] = 2.

    nebular = {}
    nebular["logU"] = -2.5

    model_components = {}
    model_components["redshift"] = 0.04
    model_components["delayed1"] = exp_young
    model_components["delayed2"] = exp_old
    model_components["dust"] = dust
    model_components["t_bc"] = 0.02
    model_components["veldisp"] = 200.
    model_components["nebular"] = nebular

    model = pipes.model_galaxy(model_components, filt_list=filter_files)

    model_dict[str(model_age) + '_' + str(model_Av) + '_' + test_id] = model

for i in range(len(combinations)):

    model_age, model_Av = combinations[i]
    print(i+1, model_age, model_Av)

    model_id = str(model_age) + '_' + str(model_Av) + '_' + test_id

    blob_to_fit = pipes.galaxy(model_id, catalog=model_dict, load_data=load_data,
                               filt_list=filter_files, phot_units=model.phot_units, spectrum_exists=False)

    exp_young_fit = {}
    exp_young_fit["age"] = (0.0, 1.0)
    exp_young_fit["tau"] = (0, 5)
    exp_young_fit["massformed"] = (0., 10.)
    exp_young_fit["metallicity"] = (0.005, 2.5)

    exp_old_fit = {}
    exp_old_fit["age"] = (4, 14)
    exp_old_fit["tau"] = (0, 15)
    exp_old_fit["massformed"] = (0, 12)
    exp_old_fit["metallicity"] = 0.5

    dust = {}
    dust["type"] = "Cardelli"
    dust["Av"] = (0, 3)

    nebular = {}
    nebular["logU"] = -2.5
    dust["eta"] = (1, 2)
    # dust["eta"] = 2

    fit_instructions = {}
    fit_instructions["redshift"] = 0.04
    fit_instructions["delayed1"] = exp_young_fit
    fit_instructions["delayed2"] = exp_old_fit
    fit_instructions["dust"] = dust
    fit_instructions["nebular"] = nebular

    fit = pipes.fit(blob_to_fit, fit_instructions)

    fit.fit(verbose=True)

    table_rows.append([model_age, model_Av, np.median(fit.posterior.samples['delayed1:age']),
                      np.median(fit.posterior.samples['dust:Av']),
                      np.median(fit.posterior.samples['delayed1:tau']),
                      np.median(fit.posterior.samples['dust:eta']),
                      np.median(fit.posterior.samples['delayed1:massformed']),
                      np.median(fit.posterior.samples['delayed1:metallicity'])])

    fig = fit.plot_spectrum_posterior(save=False, show=False)
    plt.savefig(plot_dir + '/spec_' + model_id + '.png', dpi=100)
    fig = fit.plot_sfh_posterior(save=False, show=False)
    plt.savefig(plot_dir + '/sfh_' + model_id + '.png', dpi=100)
    fig = fit.plot_corner(save=False, show=False)
    plt.savefig(plot_dir + '/corner_' + model_id + '.png', dpi=100)

    plt.close('all')


table = Table(rows=table_rows, names=['model_age', 'model_Av', 'age', 'Av', 'tau', 'eta', 'massformed', 'metallicity'])
table.write('/home/ariel/Workspace/GASP/HST/Data/models/' + test_id + '.fits', overwrite=True)

# exp_old_fit = {}
# exp_old_fit["age"] = (8.0, 14.)
# exp_old_fit["tau"] = (0.0, 30.)
# exp_old_fit["massformed"] = (0., 12.)
# exp_old_fit["metallicity"] = (0.005, 2.5)

# exp_old = {}
# exp_old["age"] = 10.
# exp_old["tau"] = 4.
# exp_old["massformed"] = 5
# exp_old["metallicity"] = 0.5
