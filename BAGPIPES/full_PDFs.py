""" 

Retrieves the full PDFs of the main variables and adds them to an astropy table

"""

import numpy as np
import matplotlib.pyplot as plt
import bagpipes as pipes
import os
from astropy.table import Table
from hst_pipeutils import load_data
from toolbox.stat_tools import gini, find_maximum, find_pdf_peaks
import pickle
from pathlib import Path

plt.ioff()

detection = 'f606w'
config = 'dexp_logprior_single'

pipes_dir = '/mnt/c/Users/ariel/Workspace/GASP/HST/BAGPIPES/'
data_dir = '/mnt/c/Users/ariel/Workspace/GASP/HST/Data/'

clump_catalog = Table.read(data_dir + detection + '_bagpipes_input.fits')
clump_catalog['fit_id'] = [clump_catalog['clump_id'][i] + '_' + config for i in range(len(clump_catalog))]

filter_files = [data_dir + 'filters/HST_WFC3_UVIS2.F275W.dat',
                data_dir + 'filters/HST_WFC3_UVIS2.F336W.dat',
                data_dir + 'filters/HST_WFC3_UVIS2.F606W.dat',
                data_dir + 'filters/HST_WFC3_UVIS2.F680N.dat',
                data_dir + 'filters/HST_WFC3_UVIS2.F814W.dat']


output_catalog = Table()
output_catalog['clump_id'] = clump_catalog['clump_id']
output_catalog['fit_id'] = clump_catalog['fit_id']
output_catalog['config'] = np.full_like(clump_catalog['fit_id'], config)
output_catalog['stellar_mass'] = np.zeros((len(clump_catalog), 500))
output_catalog['mwage'] = np.zeros((len(clump_catalog), 500))
output_catalog['Av'] = np.zeros((len(clump_catalog), 500))
output_catalog['sfr'] = np.zeros((len(clump_catalog), 500))


for clump_id in clump_catalog['fit_id']:

    clump_index = np.argwhere(clump_catalog['fit_id'] == clump_id)[0][0]

    print('>>> Loading', clump_id, '(', clump_index + 1, 'out of', len(clump_catalog), ')')

    pipes_clump = pipes.galaxy(clump_id, clump_catalog, load_data, filt_list=filter_files, spectrum_exists=False,
                               phot_units='ergscma', spec_wavs=np.arange(2400, 8100, 1))

    fit = pipes.fit(pipes_clump, fit_instructions={})

    fit.fit(verbose=False)

    fit.posterior.get_advanced_quantities()

    output_catalog['mwage'][clump_index] = fit.posterior.samples['mass_weighted_age']
    output_catalog['stellar_mass'][clump_index] = fit.posterior.samples['stellar_mass']
    output_catalog['sfr'][clump_index] = fit.posterior.samples['sfr']
    output_catalog['Av'][clump_index] = fit.posterior.samples['dust:Av'] * fit.posterior.samples['dust:eta']


output_catalog.write(data_dir + detection + '_' + config + '_full_PDFs.fits', overwrite=True)