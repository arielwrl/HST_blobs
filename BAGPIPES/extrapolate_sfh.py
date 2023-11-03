"""

ariel@oapd
06/06/2026

Extrapolates SFHs to compute luminosity evolution

"""

import numpy as np
import matplotlib.pyplot as plt
import bagpipes as pipes
from astropy.table import Table
from hst_pipeutils import load_data
from toolbox.wololo import ergstoabmags, calc_dm

data_dir = '/home/ariel/Workspace/GASP/HST/Data/'

filter_files = [data_dir + 'filters/HST_WFC3_UVIS2.F275W.dat',
                data_dir + 'filters/HST_WFC3_UVIS2.F336W.dat',
                data_dir + 'filters/HST_WFC3_UVIS2.F606W.dat',
                data_dir + 'filters/HST_WFC3_UVIS2.F680N.dat',
                data_dir + 'filters/HST_WFC3_UVIS2.F814W.dat']

filter_files_extended = [data_dir + 'filters/HST_WFC3_UVIS2.F275W.dat',
                         data_dir + 'filters/HST_WFC3_UVIS2.F336W.dat',
                         data_dir + 'filters/HST_WFC3_UVIS2.F606W.dat',
                         data_dir + 'filters/HST_WFC3_UVIS2.F680N.dat',
                         data_dir + 'filters/HST_WFC3_UVIS2.F814W.dat',
                         data_dir + 'filters/SLOAN_SDSS.r.dat',
                         data_dir + 'filters/Generic_Bessell.B.dat',
                         data_dir + 'filters/Generic_Bessell.V.dat',
                         data_dir + 'filters/Generic_Bessell.R.dat',
                         data_dir + 'filters/SLOAN_SDSS.g.dat']

output_table = Table.read('/home/ariel/Workspace/GASP/HST/Data/f606w_dexp_logprior_single_bagpipes_results.fits')
input_table = Table.read('/home/ariel/Workspace/GASP/HST/Data/f606w_bagpipes_input.fits')
input_table['fit_id'] = [input_table['clump_id'][i] + '_' + 'dexp_logprior_single' for i in range(len(input_table))]

time_steps = np.arange(0, 10.3, 0.3)

complex_evolution = Table()
complex_evolution['complex_id'] = input_table['clump_id']
complex_evolution['time_steps'] = [time_steps for i in range(len(complex_evolution))]
complex_evolution['magnitudes'] = np.zeros_like(complex_evolution['time_steps'])
complex_evolution['magnitudes_r'] = np.zeros_like(complex_evolution['time_steps'])
complex_evolution['magnitudes_g'] = np.zeros_like(complex_evolution['time_steps'])
complex_evolution['magnitudes_B'] = np.zeros_like(complex_evolution['time_steps'])
complex_evolution['magnitudes_V'] = np.zeros_like(complex_evolution['time_steps'])
complex_evolution['magnitudes_R'] = np.zeros_like(complex_evolution['time_steps'])
complex_evolution['absolute_magnitudes'] = np.zeros_like(complex_evolution['time_steps'])
complex_evolution['absolute_magnitudes_r'] = np.zeros_like(complex_evolution['time_steps'])
complex_evolution['absolute_magnitudes_g'] = np.zeros_like(complex_evolution['time_steps'])
complex_evolution['absolute_magnitudes_B'] = np.zeros_like(complex_evolution['time_steps'])
complex_evolution['absolute_magnitudes_V'] = np.zeros_like(complex_evolution['time_steps'])
complex_evolution['absolute_magnitudes_R'] = np.zeros_like(complex_evolution['time_steps'])
complex_evolution['stellar_mass'] = np.zeros_like(complex_evolution['time_steps'])

for complex_index in range(len(output_table)):

    complex_fit_id = output_table['fit_id'][complex_index]

    print('>>> Loading', complex_fit_id, '(', complex_index + 1, 'out of', len(output_table), ')')

    pipes_complex = pipes.galaxy(complex_fit_id, input_table, load_data, filt_list=filter_files, spectrum_exists=False,
                                 phot_units='ergscma', spec_wavs=np.arange(2400, 8100, 1))

    fit = pipes.fit(pipes_complex, fit_instructions={})

    dexp = {}
    dexp["age"] = np.median(fit.posterior.samples['delayed:age'])
    dexp["tau"] = np.median(fit.posterior.samples['delayed:tau'])
    dexp["massformed"] = np.median(fit.posterior.samples['formed_mass'])
    dexp["metallicity"] = np.median(fit.posterior.samples['delayed:metallicity'])

    dust = {}
    dust["type"] = "Cardelli"

    dust["Av"] = np.median(fit.posterior.samples['dust:Av'])
    dust["eta"] = np.median(fit.posterior.samples['dust:eta'])

    nebular = {}
    nebular["logU"] = np.median(fit.posterior.samples['nebular:logU'])

    model_components = {}
    model_components["redshift"] = 0.0045

    model_components["delayed"] = dexp

    model_components["dust"] = dust
    model_components["nebular"] = nebular
    model_components["t_bc"] = 0.02

    model_original = pipes.model_galaxy(model_components, filt_list=filter_files, spec_wavs=np.arange(2400, 9100, 1))

    ages = model_original.sfh.ages
    age = model_original.model_comp['delayed']['age'] * 10 ** 9
    tau = model_original.model_comp['delayed']['tau'] * 10 ** 9
    custom_sfh = np.zeros_like(ages)
    t = age - ages[ages < age]
    custom_sfh[ages < age] = t * np.exp(-t/tau)
    sfh_norm = np.max(model_original.sfh.sfh)/np.max(custom_sfh)
    custom_sfh *= sfh_norm
    formed_mass = np.log10(np.sum(model_original.sfh.age_widths * custom_sfh))

    magnitudes = []
    magnitudes_r = []
    magnitudes_g = []
    magnitudes_B = []
    magnitudes_V = []
    magnitudes_R = []
    absolute_magnitudes = []
    absolute_magnitudes_r = []
    absolute_magnitudes_g = []
    absolute_magnitudes_B = []
    absolute_magnitudes_V = []
    absolute_magnitudes_R = []
    stellar_mass = []

    for time_step in time_steps:
        age_new = age + time_step * 10 ** 9
        custom_sfh_new = np.zeros_like(ages)
        t = age_new - ages[ages < age_new]
        custom_sfh_new[ages < age_new] = t * np.exp(-t/tau)
        custom_sfh_new *= sfh_norm
        formed_mass_new = np.log10(np.sum(model_original.sfh.age_widths * custom_sfh_new))

        dexp["age"] = age_new / 1e9
        dexp["massformed"] = formed_mass_new
        model_components["delayed"] = dexp

        model_new = pipes.model_galaxy(model_components, filt_list=filter_files_extended, spec_wavs=np.arange(2400, 9100, 1))

        dm = 31.51

        magnitudes.append(ergstoabmags(model_new.photometry[2], 5887.70))
        magnitudes_r.append(ergstoabmags(model_new.photometry[5], 6175.58))
        magnitudes_g.append(ergstoabmags(model_new.photometry[9], 4702.50))
        magnitudes_B.append(ergstoabmags(model_new.photometry[6], 4371.07))
        magnitudes_V.append(ergstoabmags(model_new.photometry[7], 5466.87))
        magnitudes_R.append(ergstoabmags(model_new.photometry[8], 6498.09))
        absolute_magnitudes.append(ergstoabmags(model_new.photometry[2], 5887.70) - dm)
        absolute_magnitudes_r.append(ergstoabmags(model_new.photometry[5], 6175.58) - dm)
        absolute_magnitudes_g.append(ergstoabmags(model_new.photometry[9], 4702.50) - dm)
        absolute_magnitudes_B.append(ergstoabmags(model_new.photometry[6], 4371.07) - dm)
        absolute_magnitudes_V.append(ergstoabmags(model_new.photometry[7], 5466.87) - dm)
        absolute_magnitudes_R.append(ergstoabmags(model_new.photometry[8], 6498.09) - dm)
        stellar_mass.append(model_new.sfh.stellar_mass)

    complex_evolution['magnitudes'][complex_index] = magnitudes
    complex_evolution['absolute_magnitudes'][complex_index] = absolute_magnitudes
    complex_evolution['magnitudes_r'][complex_index] = magnitudes_r
    complex_evolution['absolute_magnitudes_r'][complex_index] = absolute_magnitudes_r
    complex_evolution['magnitudes_g'][complex_index] = magnitudes_g
    complex_evolution['absolute_magnitudes_g'][complex_index] = absolute_magnitudes_g
    complex_evolution['magnitudes_B'][complex_index] = magnitudes_B
    complex_evolution['absolute_magnitudes_B'][complex_index] = absolute_magnitudes_B
    complex_evolution['magnitudes_V'][complex_index] = magnitudes_V
    complex_evolution['absolute_magnitudes_V'][complex_index] = absolute_magnitudes_V
    complex_evolution['magnitudes_R'][complex_index] = magnitudes_R
    complex_evolution['absolute_magnitudes_R'][complex_index] = absolute_magnitudes_R
    complex_evolution['stellar_mass'][complex_index] = stellar_mass

complex_evolution.write('/home/ariel/Workspace/GASP/HST/Data/complex_evolution_fornax_10.fits', overwrite=True)

