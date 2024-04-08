"""

ariel@oapd
06/11/2023

Re-generates models to calculate the mass in young components

"""

import numpy as np
import matplotlib.pyplot as plt
import bagpipes as pipes
from astropy.table import Table
from hst_pipeutils import load_data

data_dir = '/mnt/c/Users/ariel/Workspace/GASP/HST/Data/'

filter_files = [data_dir + 'filters/HST_WFC3_UVIS2.F275W.dat',
                data_dir + 'filters/HST_WFC3_UVIS2.F336W.dat',
                data_dir + 'filters/HST_WFC3_UVIS2.F606W.dat',
                data_dir + 'filters/HST_WFC3_UVIS2.F680N.dat',
                data_dir + 'filters/HST_WFC3_UVIS2.F814W.dat']

input_table = Table.read(data_dir + 'disk_catalogs/disk_f275w_bagpipes_input.fits')
output_table = Table.read(data_dir + 'disk_catalogs/disk_f275w_dexp_logprior_bagpipes_results.fits')

flag = (input_table['leaf_flag'] == 1) & (output_table['mass_ratio'] < 2)

input_table = input_table[flag]
output_table = output_table[flag]

input_table['fit_id'] = [input_table['fit_id'][i] + '_dexp_logprior' for i in range(len(input_table))]
output_table['stellar_mass_young'] = np.zeros_like(output_table['stellar_mass'])
output_table['mwage_young'] = np.zeros_like(output_table['stellar_mass'])
output_table['mass_ratio'] = np.zeros_like(output_table['stellar_mass'])

for clump_index in range(len(output_table)):

    clump_fit_id = output_table['fit_id'][clump_index] + '_dexp_logprior'

    print('>>> Loading', clump_fit_id, '(', clump_index + 1, 'out of', len(output_table), ')')

    pipes_clump = pipes.galaxy(clump_fit_id, input_table, load_data, filt_list=filter_files, spectrum_exists=False,
                               phot_units='ergscma', spec_wavs=np.arange(2400, 9100, 1), )

    fit = pipes.fit(pipes_clump, fit_instructions={})

    fit.posterior.get_advanced_quantities()

    dexp = {}
    dexp["age"] = np.median(fit.posterior.samples['delayed2:age'])
    dexp["tau"] = np.median(fit.posterior.samples['delayed2:tau'])
    dexp["massformed"] = np.median(fit.posterior.samples['delayed2:massformed'])
    dexp["metallicity"] = np.median(fit.posterior.samples['delayed2:metallicity'])

    dust = {}
    dust["type"] = "Cardelli"

    dust["Av"] = np.median(fit.posterior.samples['dust:Av'])
    dust["eta"] = np.median(fit.posterior.samples['dust:eta'])

    nebular = {}
    nebular["logU"] = -2.5

    model_components = {}
    model_components["redshift"] = input_table['galaxy_redshift'][clump_index]

    model_components["delayed"] = dexp

    model_components["dust"] = dust
    model_components["nebular"] = nebular
    model_components["t_bc"] = 0.02

    model_young = pipes.model_galaxy(model_components, filt_list=filter_files, spec_wavs=np.arange(2400, 9100, 1))

    output_table['stellar_mass_young'][clump_index] = model_young.sfh.stellar_mass
    output_table['mwage_young'][clump_index] = model_young.sfh.mass_weighted_age
    output_table['mass_ratio'][clump_index] = output_table['formed_mass_old'][clump_index] - output_table['formed_mass_young'][clump_index]

output_table.write(data_dir + 'disk_catalogs/disk_f275w_dexp_logprior_bagpipes_results_flagged.fits', overwrite=True)
input_table.write(data_dir + 'disk_catalogs/disk_f275w_bagpipes_input_flagged.fits', overwrite=True)

