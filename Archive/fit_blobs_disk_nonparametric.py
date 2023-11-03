"""

ariel@oapd
22/04/2022

Implements the SINOPSIS way of fitting (12 age bins) in BAGPIPES and hopes for the best - Disk version

Streamlined on 16/05/2022

"""

import numpy as np
import matplotlib.pyplot as plt
import bagpipes as pipes
import os
from astropy.table import Table
from starlight_toolkit.synphot import synflux
from hst_pipeutils import load_data, nonparametric_plots
from pysinopsis.utils import calc_mwage

plt.ioff()

detection = 'halpha'
config = 'sfrprior_nonparametric'
run_flag = False

pipes_dir = '/BAGPIPES/'
data_dir = '/Data/'
plot_dir = pipes_dir + 'plots/'
sfh_dir = pipes_dir + 'sfh/'

if detection + '_' + config + '_disk' not in os.listdir(plot_dir):
    os.makedirs(plot_dir + detection + '_' + config + '_disk')
if detection + '_' + config + '_disk' not in os.listdir(sfh_dir):
    os.makedirs(sfh_dir + detection + '_' + config + '_disk')

blob_catalog = Table.read(data_dir + 'disk_' + detection + '_bagpipes_input.fits')

filter_files = [data_dir + 'filters/HST_WFC3_UVIS2.F275W.dat',
                data_dir + 'filters/HST_WFC3_UVIS2.F336W.dat',
                data_dir + 'filters/HST_WFC3_UVIS2.F606W.dat',
                data_dir + 'filters/HST_WFC3_UVIS2.F680N.dat',
                data_dir + 'filters/HST_WFC3_UVIS2.F814W.dat']

original_ids = blob_catalog['blob_id']
blob_catalog['blob_id'] = [blob_catalog['blob_id'][i] + '_' + config for i in range(len(blob_catalog))]

if run_flag is False:
    output_catalog = Table()
    output_catalog['blob_id'] = original_ids
    output_catalog['galaxy'] = np.zeros(len(blob_catalog), dtype='S32')
    output_catalog['metallicity'] = np.zeros((len(blob_catalog), 12))
    output_catalog['formed_mass'] = np.zeros((len(blob_catalog), 12))
    output_catalog['sfh'] = np.zeros((len(blob_catalog), 12))
    output_catalog['sfh25'] = np.zeros((len(blob_catalog), 12))
    output_catalog['sfh75'] = np.zeros((len(blob_catalog), 12))
    output_catalog['Av'] = np.zeros(len(blob_catalog))
    output_catalog['eta'] = np.zeros(len(blob_catalog))
    output_catalog['sfr'] = np.zeros(len(blob_catalog))
    output_catalog['stellar_mass'] = np.zeros(len(blob_catalog))
    output_catalog['mwage'] = np.zeros(len(blob_catalog))
    output_catalog['mwage_alt'] = np.zeros(len(blob_catalog))
    output_catalog['tform'] = np.zeros(len(blob_catalog))
    output_catalog['av_detected'] = np.ones(len(blob_catalog))
    output_catalog['photometry_syn'] = np.zeros((len(blob_catalog), 5))
    output_catalog['photometry_obs'] = np.zeros((len(blob_catalog), 5))
    output_catalog['photometry_nebular'] = np.zeros((len(blob_catalog), 5))
    output_catalog['chi2'] = np.zeros(len(blob_catalog))
    output_catalog['Ha'] = np.zeros(len(blob_catalog))
    output_catalog['Nii'] = np.zeros(len(blob_catalog))
    output_catalog['Oiii'] = np.zeros(len(blob_catalog))
    output_catalog['Hb'] = np.zeros(len(blob_catalog))

age_bins = [0, 2.0e6, 4.0e6, 7.0e6, 2.0e7, 5.5e7, 2.0e8, 5.5e8, 1.0e9, 3.0e9,
            5.75e9, 1.0e10, 1.4e10]
age_bins_len = [age_bins[i+1]-age_bins[i] for i in range(12)]
age_bins_mid = [(age_bins[i+1]+age_bins[i])/2 for i in range(12)]

# for blob_id in blob_catalog['blob_id']:
for blob_id in blob_catalog['blob_id'][0:20]:

    blob_index = np.argwhere(blob_catalog['blob_id'] == blob_id)[0][0]

    print('>>> Fitting', blob_id, '(', blob_index + 1, 'out of', len(blob_catalog), ')')
    file_list = os.listdir(pipes_dir + 'pipes/posterior/')
    if ((blob_id+'.h5' in file_list) or (blob_id+'_resume.dat' in file_list)) and run_flag is True:
        continue
    if (blob_id+'_resume.dat' in file_list) and run_flag is False:
        continue

    pipes_blob = pipes.galaxy(blob_id, blob_catalog, load_data, filt_list=filter_files, spectrum_exists=False,
                              phot_units='ergscma')

    square_bursts = []

    # for i in range(len(age_bins)-1):
    #     constant = {}
    #     constant["metallicity"] = (0.005, 2.5)
    #     constant["age_max"] = age_bins[i+1] / 1e9
    #     constant["age_min"] = age_bins[i] / 1e9
    #     constant["massformed"] = (0, 10)
    #     square_bursts.append(constant)

    for i in range(len(age_bins)-1):
        constant = {}
        constant["metallicity"] = (0.005, 2.5)
        constant["age_max"] = age_bins[i+1] / 1e9
        constant["age_min"] = age_bins[i] / 1e9
        constant["massformed"] = (0, np.log10(4000 * age_bins_len[i]))
        square_bursts.append(constant)

    dust = {}
    dust["type"] = "Cardelli"

    dust["Av"] = (0, 1.25)
    dust["eta"] = (0.5, 2.5)

    nebular = {}
    nebular["logU"] = -2.5

    fit_instructions = {}
    for i in range(len(age_bins) - 1):
        fit_instructions['constant'+str(i)] = square_bursts[i]
    fit_instructions["redshift"] = (blob_catalog['galaxy_redshift'][blob_index]-0.0001,
                                    blob_catalog['galaxy_redshift'][blob_index]+0.0001)
    fit_instructions["dust"] = dust
    fit_instructions["nebular"] = nebular
    fit_instructions["t_bc"] = 0.02

    fit = pipes.fit(pipes_blob, fit_instructions)

    fit.fit(verbose=False)

    fit.posterior.get_advanced_quantities()

    sfh = [np.median(10 ** fit.posterior.samples['constant' + str(i) + ':massformed']) /
            age_bins_len[i] for i in range(12)]
    sfh_25 = [np.percentile(10 ** fit.posterior.samples['constant' + str(i) + ':massformed'],
                            25) / age_bins_len[i] for i in range(12)]
    sfh_75 = [np.percentile(10 ** fit.posterior.samples['constant' + str(i) + ':massformed'],
                            75) / age_bins_len[i] for i in range(12)]

    sfh_full = np.array([(10 ** fit.posterior.samples['constant' + str(i) + ':massformed']) /
                          age_bins_len[i] for i in range(12)])

    sfh_minchi2 = fit.posterior.samples['sfh'][np.argmin(fit.posterior.samples['chisq_phot'])]

    masses_full = np.array([fit.posterior.samples['constant' + str(i) + ':massformed'] for i in range(12)])

    sfh_median_bp = np.median(fit.posterior.samples['sfh'], axis=0)
    sfh_25_bp = np.percentile(fit.posterior.samples['sfh'], axis=0, q=25)
    sfh_75_bp = np.percentile(fit.posterior.samples['sfh'], axis=0, q=75)
    ages_bp = fit.posterior.sfh.ages

    np.savetxt(sfh_dir + detection + '_' + config + '_disk/sfh_' + blob_id + '.txt',
               np.array([age_bins_mid, sfh, sfh_25, sfh_75]).transpose(),
               header='Age SFR SFR(25th percentile) SFR(75th percentile)')

    np.savetxt(sfh_dir + detection + '_' + config + '_disk/sfh_' + blob_id + '_full.txt', sfh_full)

    np.savetxt(sfh_dir + detection + '_' + config + '_disk/masses_' + blob_id + '_full.txt', masses_full)

    np.savetxt(sfh_dir + detection + '_' + config + '_disk/sfh_' + blob_id + '_bagpipes.txt',
               np.array([ages_bp, sfh_median_bp, sfh_25_bp, sfh_75_bp]).transpose(),
               header='Age SFR SFR(25th percentile) SFR(75th percentile)')

    np.savetxt(sfh_dir + detection + '_' + config + '_disk/sfh_' + blob_id + '_minchi2.txt',
               np.array([ages_bp, sfh_minchi2]).transpose(),
               header='Age SFR SFR(25th percentile) SFR(75th percentile)')

    nonparametric_plots(fit, blob_id, blob_index, blob_catalog, plot_dir + detection + '_' + config + '_disk')

    if run_flag is False:

        output_catalog['galaxy'][blob_index] = blob_id.split('_')[0]
        output_catalog['metallicity'][blob_index] = [
            np.median(fit.posterior.samples['constant' + str(i) + ':metallicity']) for i in range(12)]
        output_catalog['formed_mass'][blob_index] = [
            np.median(fit.posterior.samples['constant' + str(i) + ':massformed']) for i in range(12)]
        output_catalog['Av'][blob_index] = np.median(fit.posterior.samples['dust:Av'])
        output_catalog['eta'][blob_index] = np.median(fit.posterior.samples['dust:eta'])
        output_catalog['sfr'][blob_index] = np.median(fit.posterior.samples['sfr'])
        output_catalog['stellar_mass'][blob_index] = np.median(fit.posterior.samples['stellar_mass'])
        output_catalog['mwage'][blob_index] = np.median(fit.posterior.samples['mass_weighted_age'])
        output_catalog['mwage_alt'][blob_index] = calc_mwage(age_bins_mid, sfh, age_bins_len)
        output_catalog['tform'][blob_index] = np.median(fit.posterior.samples['tform'])
        output_catalog['photometry_syn'][blob_index] = np.median(fit.posterior.samples['photometry'], axis=0)
        output_catalog['photometry_obs'][blob_index] = fit.galaxy.photometry[:, 1]
        output_catalog['photometry_nebular'][blob_index] = [synflux(fit.posterior.model_galaxy.wavelengths,
                                                                    fit.posterior.model_galaxy.nebular_spectrum_full,
                                                                    filter_file) for filter_file in filter_files]

        print(output_catalog['photometry_nebular'][blob_index])
        print(output_catalog['photometry_syn'][blob_index])
        print(output_catalog['photometry_nebular'][blob_index] / output_catalog['photometry_syn'][blob_index])

        output_catalog['chi2'][blob_index] = np.min(fit.posterior.samples['chisq_phot'])
        output_catalog['Ha'][blob_index] = fit.posterior.model_galaxy.line_fluxes['H  1  6562.81A']
        output_catalog['Nii'][blob_index] = fit.posterior.model_galaxy.line_fluxes['N  2  6583.45A']
        output_catalog['Oiii'][blob_index] = fit.posterior.model_galaxy.line_fluxes['O  3  5006.84A']
        output_catalog['Hb'][blob_index] = fit.posterior.model_galaxy.line_fluxes['H  1  4861.33A']

        output_catalog['sfh'][blob_index] = sfh
        output_catalog['sfh25'][blob_index] = sfh_25
        output_catalog['sfh75'][blob_index] = sfh_75

# output_catalog.write(data_dir + 'disk_' + detection + '_' + config + '_bagpipes_results_test_sample.fits', overwrite=True)
