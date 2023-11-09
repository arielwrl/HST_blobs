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

config = 'dexp_logprior_single'

pipes_dir = '/mnt/c/Users/ariel/Workspace/GASP/HST/BAGPIPES/'
data_dir = '/mnt/c/Users/ariel/Workspace/GASP/HST/Data/'

halpha_catalog = Table.read(data_dir + 'halpha_bagpipes_input.fits')
halpha_catalog['fit_id'] = [halpha_catalog['clump_id'][i] + '_' + config for i in range(len(halpha_catalog))]

f275w_catalog = Table.read(data_dir + 'f275w_bagpipes_input.fits')
f275w_catalog['fit_id'] = [f275w_catalog['clump_id'][i] + '_' + config for i in range(len(f275w_catalog))]

f606w_catalog = Table.read(data_dir + 'f606w_bagpipes_input.fits')
f606w_catalog['fit_id'] = [f606w_catalog['clump_id'][i] + '_' + config for i in range(len(f606w_catalog))]

filter_files = [data_dir + 'filters/HST_WFC3_UVIS2.F275W.dat',
                data_dir + 'filters/HST_WFC3_UVIS2.F336W.dat',
                data_dir + 'filters/HST_WFC3_UVIS2.F606W.dat',
                data_dir + 'filters/HST_WFC3_UVIS2.F680N.dat',
                data_dir + 'filters/HST_WFC3_UVIS2.F814W.dat']

halpha_id = 'JO201_A40_halpha_dexp_logprior_single'
f275w_id = 'JO201_A83_f275w_dexp_logprior_single'
f606w_id = 'JO201_1625_f606w_dexp_logprior_single'

halpha_clump = pipes.galaxy(halpha_id, halpha_catalog, load_data, filt_list=filter_files, spectrum_exists=False,
                            phot_units='ergscma', spec_wavs=np.arange(2400, 8100, 1))
halpha_fit = pipes.fit(halpha_clump, fit_instructions={})

f275w_clump = pipes.galaxy(f275w_id, f275w_catalog, load_data, filt_list=filter_files, spectrum_exists=False,
                            phot_units='ergscma', spec_wavs=np.arange(2400, 8100, 1))
f275w_fit = pipes.fit(f275w_clump, fit_instructions={})

f606w_clump = pipes.galaxy(f606w_id, f606w_catalog, load_data, filt_list=filter_files, spectrum_exists=False,
                            phot_units='ergscma', spec_wavs=np.arange(2400, 8100, 1))
f606w_fit = pipes.fit(f606w_clump, fit_instructions={})

example_pdfs = Table()
example_pdfs['clump_ids'] = [halpha_id, f275w_id, f606w_id]
example_pdfs['mwage'] = [halpha_fit.posterior.samples['mass_weighted_age'], f275w_fit.posterior.samples['mass_weighted_age'], 
                         f606w_fit.posterior.samples['mass_weighted_age']]
example_pdfs['stellar_mass'] = [halpha_fit.posterior.samples['stellar_mass'], f275w_fit.posterior.samples['stellar_mass'], 
                                f606w_fit.posterior.samples['stellar_mass']]
example_pdfs['sfr'] = [np.log10(halpha_fit.posterior.samples['sfr']), np.log10(f275w_fit.posterior.samples['sfr']), 
                                np.log10(f606w_fit.posterior.samples['sfr'])]
example_pdfs['Av'] = [halpha_fit.posterior.samples['dust:Av'] * halpha_fit.posterior.samples['dust:eta'], 
                      f275w_fit.posterior.samples['dust:Av'] * f275w_fit.posterior.samples['dust:eta'], 
                      f606w_fit.posterior.samples['dust:Av'] * f606w_fit.posterior.samples['dust:eta']]

example_pdfs.write(data_dir + 'example_pdfs.fits', overwrite=True)


# output_catalog = Table()
# output_catalog['clump_id'] = clump_catalog['clump_id']
# output_catalog['fit_id'] = clump_catalog['fit_id']
# output_catalog['config'] = np.full_like(clump_catalog['fit_id'], config)
# output_catalog['stellar_mass'] = np.zeros((len(clump_catalog), 500))
# output_catalog['mwage'] = np.zeros((len(clump_catalog), 500))
# output_catalog['Av'] = np.zeros((len(clump_catalog), 500))
# output_catalog['sfr'] = np.zeros((len(clump_catalog), 500))


# for clump_id in clump_catalog['fit_id']:

#     clump_index = np.argwhere(clump_catalog['fit_id'] == clump_id)[0][0]

#     print('>>> Loading', clump_id, '(', clump_index + 1, 'out of', len(clump_catalog), ')')

#     pipes_clump = pipes.galaxy(clump_id, clump_catalog, load_data, filt_list=filter_files, spectrum_exists=False,
#                                phot_units='ergscma', spec_wavs=np.arange(2400, 8100, 1))

#     fit = pipes.fit(pipes_clump, fit_instructions={})

#     fit.fit(verbose=False)

#     fit.posterior.get_advanced_quantities()

#     output_catalog['mwage'][clump_index] = fit.posterior.samples['mass_weighted_age']
#     output_catalog['stellar_mass'][clump_index] = fit.posterior.samples['stellar_mass']
#     output_catalog['sfr'][clump_index] = fit.posterior.samples['sfr']
#     output_catalog['Av'][clump_index] = fit.posterior.samples['dust:Av'] * fit.posterior.samples['dust:eta']


# output_catalog.write(data_dir + detection + '_' + config + '_full_PDFs.fits', overwrite=True)