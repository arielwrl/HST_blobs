"""

ariel@oapd
30/03/2022

Check the SEDs outside an extraplanar blob selected by Marco, we should see if a significant stellar mass is detected
to test if the low stellar masses detected for old components in double delayed tau models are physical

"""

import numpy as np
import matplotlib.pyplot as plt
import bagpipes as pipes
import os
from astropy.table import Table
from matplotlib.backends.backend_pdf import PdfPages
from starlight_toolkit.synphot import synflux
from toolbox.stat_tools import int_to_bool_list

plt.ioff()

data_dir = '/home/ariel/Workspace/GASP/HST/Data/'
plot_dir = '/home/ariel/Workspace/GASP/HST/BAGPIPES/plots/'

filter_files = [data_dir + 'filters/HST_WFC3_UVIS2.F275W.dat',
                data_dir + 'filters/HST_WFC3_UVIS2.F336W.dat',
                data_dir + 'filters/HST_WFC3_UVIS2.F606W.dat',
                data_dir + 'filters/HST_WFC3_UVIS2.F680N.dat',
                data_dir + 'filters/HST_WFC3_UVIS2.F814W.dat']

hst_filters = ['F275W', 'F336W', 'F606W', 'F680N', 'F814W']


def load_data(blob_id):

    if blob_id is 'mark1_realer':
        fluxes = np.array([1.85e-19, 2.28e-19, 3.14e-19, 4.97e-19, 2.13e-19])
    if blob_id is 'mark2_realer':
        fluxes = np.array([7.81e-19, 1.12e-18, 1.86e-18, 1.82e-18, 1.51e-18])
    if blob_id is 'mark3_realer':
        fluxes = np.array([3.85e-19, 8.69e-19, 1.33e-18, 1.41e-18, 1.10e-18])

    # errors = 0.15 * fluxes
    errors = [1e-19, 1e-19, 5e-20, 1e-19, 1e-20]

    return np.array([fluxes, errors]).transpose()


id_list = ['mark1_realer', 'mark2_realer', 'mark3_realer']
results_list = []

for blob_id in id_list:

    pipes_blob = pipes.galaxy(blob_id, load_data, filt_list=filter_files, spectrum_exists=False, phot_units='ergscma')

    exp = {}
    exp["age"] = (8.0, 14.)
    exp["tau"] = (0.0, 30.)
    exp["massformed"] = (0., 12.)
    exp["metallicity"] = (0.005, 2.5)

    exp_young = {}
    exp_young["age"] = (0.0, 1.0)
    exp_young["tau"] = (0.0, 10.)
    exp_young["massformed"] = (0., 10.)
    exp_young["metallicity"] = (0.005, 2.5)

    dust = {}
    dust["type"] = "Cardelli"

    dust["Av"] = (0, 1.25)
    dust["eta"] = (0.5, 3.0)

    nebular = {}
    nebular["logU"] = -2.5

    fit_instructions = {}
    fit_instructions["redshift"] = 0.04675
    fit_instructions["delayed1"] = exp
    fit_instructions["delayed2"] = exp_young
    fit_instructions["dust"] = dust
    fit_instructions["nebular"] = nebular
    fit_instructions["t_bc"] = 0.02

    fit = pipes.fit(pipes_blob, fit_instructions)

    fit.fit(verbose=False)

    fit.posterior.get_advanced_quantities()

    results_list.append(fit)

    fig = fit.plot_spectrum_posterior(save=False, show=False)
    fig[0].set_size_inches(10, 6)
    plt.savefig('plots/fit_'+blob_id+'.png', dpi=100)
    fig = fit.plot_sfh_posterior(save=False, show=False)
    fig[0].set_size_inches(10, 6)
    plt.savefig('plots//sfh_'+blob_id+'.png', dpi=100)
    fig = fit.plot_corner(save=False, show=False)
    plt.savefig('plots//par_'+blob_id+'.png', dpi=89)

    plt.close('all')


