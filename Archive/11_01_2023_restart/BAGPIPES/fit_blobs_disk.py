"""

ariel@oapd
10/05/2022

Adapted from the newer version of fit_blobs_tail, but for disks

"""

import numpy as np
import matplotlib.pyplot as plt
import bagpipes as pipes
import os
from astropy.table import Table
from starlight_toolkit.synphot import synflux
from hst_pipeutils import load_data, parametric_plots
from toolbox.stat_tools import gini, find_maximum, find_pdf_peaks
import gc

plt.ioff()

detection = 'f275w'
config = 'dexp_logprior'
run_flag = False

pipes_dir = '/home/ariel/Workspace/GASP/HST/BAGPIPES/'
data_dir = '/home/ariel/Workspace/GASP/HST/Data/'
plot_dir = pipes_dir + 'plots/'
sfh_dir = pipes_dir + 'sfh/'

if detection + '_' + config + '_extraplanar' not in os.listdir(plot_dir):
    os.makedirs(plot_dir + detection + '_' + config + '_extraplanar')
if detection + '_' + config + '_extraplanar' not in os.listdir(sfh_dir):
    os.makedirs(sfh_dir + detection + '_' + config + '_extraplanar')

blob_catalog = Table.read(data_dir + 'extraplanar_' + detection + '_bagpipes_input.fits')
blob_catalog = blob_catalog[blob_catalog['sel_flag'] == 31]

filter_files = [data_dir + 'filters/HST_WFC3_UVIS2.F275W.dat',
                data_dir + 'filters/HST_WFC3_UVIS2.F336W.dat',
                data_dir + 'filters/HST_WFC3_UVIS2.F606W.dat',
                data_dir + 'filters/HST_WFC3_UVIS2.F680N.dat',
                data_dir + 'filters/HST_WFC3_UVIS2.F814W.dat']

original_ids = blob_catalog['fit_id']
blob_catalog['fit_id'] = [blob_catalog['fit_id'][i] + '_' + config for i in range(len(blob_catalog))]

if run_flag is False:

    # blacklist = []
    # gray_list = []
    # all_runs_rows = []

    output_catalog = Table()
    output_catalog['blob_id'] = blob_catalog['blob_id']
    output_catalog['fit_id'] = original_ids
    output_catalog['fit_id_fitconfig'] = blob_catalog['fit_id']
    output_catalog['config'] = np.full_like(original_ids, config)
    output_catalog['galaxy'] = np.zeros(len(blob_catalog), dtype='S32')
    output_catalog['age_young'] = np.zeros(len(blob_catalog))
    output_catalog['age_young_std'] = np.zeros(len(blob_catalog))
    output_catalog['age_young_iqr'] = np.zeros(len(blob_catalog))
    output_catalog['age_young_gini'] = np.zeros(len(blob_catalog))
    output_catalog['tau_young'] = np.zeros(len(blob_catalog))
    output_catalog['metallicity_young'] = np.zeros(len(blob_catalog))
    output_catalog['formed_mass_young'] = np.zeros(len(blob_catalog))
    output_catalog['age_old'] = np.zeros(len(blob_catalog))
    output_catalog['tau_old'] = np.zeros(len(blob_catalog))
    output_catalog['metallicity_old'] = np.zeros(len(blob_catalog))
    output_catalog['formed_mass_old'] = np.zeros(len(blob_catalog))
    output_catalog['Av'] = np.zeros(len(blob_catalog))
    output_catalog['eta'] = np.zeros(len(blob_catalog))
    output_catalog['sfr'] = np.zeros(len(blob_catalog))
    output_catalog['stellar_mass'] = np.zeros(len(blob_catalog))
    output_catalog['mwage'] = np.zeros(len(blob_catalog))
    output_catalog['tform'] = np.zeros(len(blob_catalog))
    output_catalog['photometry_syn'] = np.zeros((len(blob_catalog), 5))
    output_catalog['photometry_obs'] = np.zeros((len(blob_catalog), 5))
    output_catalog['photometry_nebular'] = np.zeros((len(blob_catalog), 5))
    output_catalog['chi2'] = np.zeros(len(blob_catalog))
    output_catalog['Ha'] = np.zeros(len(blob_catalog))
    output_catalog['Nii'] = np.zeros(len(blob_catalog))
    output_catalog['Oiii'] = np.zeros(len(blob_catalog))
    output_catalog['Hb'] = np.zeros(len(blob_catalog))
    output_catalog['age_young_peaks'] = np.zeros((len(blob_catalog), 15))
    output_catalog['Av_peaks'] = np.zeros((len(blob_catalog), 15))
    # output_catalog['mwage_peaks'] = np.zeros((len(blob_catalog), 15))
    output_catalog['n_age_young_peaks'] = np.zeros((len(blob_catalog)))
    output_catalog['n_Av_peaks'] = np.zeros((len(blob_catalog)))
    output_catalog['n_mwage_peaks'] = np.zeros((len(blob_catalog)))

skip = 0

for blob_id in blob_catalog['fit_id']:

    blob_index = np.argwhere(blob_catalog['fit_id'] == blob_id)[0][0]

    print(blob_catalog['tail_gal_flag'][blob_index])

    if blob_index < skip:
        continue

    print('>>> Fitting', blob_id, '(', blob_index + 1, 'out of', len(blob_catalog), ')')
    file_list = os.listdir(pipes_dir + 'pipes/posterior/')
    if ((blob_id+'.h5' in file_list) or (blob_id+'_resume.dat' in file_list)) and run_flag is True:
        continue

    pipes_blob = pipes.galaxy(blob_id, blob_catalog, load_data, filt_list=filter_files, spectrum_exists=False,
                              phot_units='ergscma')

    exp = {}
    exp["age"] = (8.0, 14.)
    exp["tau"] = (0.1, 10.)
    exp["massformed"] = (0., 12.)
    exp["metallicity"] = (0.005, 2.5)

    exp_young = {}
    exp_young["age"] = (0.0, 0.5)
    exp_young["tau"] = (0.0000000001, 0.5)
    exp_young["tau_prior"] = 'log_10'
    exp_young["massformed"] = (0., 10.)
    exp_young["metallicity"] = (0.005, 2.5)

    dust = {}
    dust["type"] = "Cardelli"

    dust["Av"] = (0, 1.25)
    dust["eta"] = (1, 2.5)

    nebular = {}
    nebular["logU"] = -2.5

    fit_instructions = {}
    fit_instructions["redshift"] = (blob_catalog['galaxy_redshift'][blob_index]-0.0001,
                                    blob_catalog['galaxy_redshift'][blob_index]+0.0001)
    fit_instructions["delayed1"] = exp
    fit_instructions["delayed2"] = exp_young
    fit_instructions["dust"] = dust
    fit_instructions["nebular"] = nebular
    fit_instructions["t_bc"] = 0.02

    fit = pipes.fit(pipes_blob, fit_instructions)

    fit.fit(verbose=False)

    if run_flag:
        sfh_median = np.median(fit.posterior.samples['sfh'], axis=0)
        sfh_25 = np.percentile(fit.posterior.samples['sfh'], axis=0, q=25)
        sfh_75 = np.percentile(fit.posterior.samples['sfh'], axis=0, q=75)
        ages = fit.posterior.sfh.ages

        np.savetxt(sfh_dir + detection + '_' + config + '_extraplanar/sfh_' + blob_id + '.txt',
                   np.array([ages, sfh_median, sfh_25, sfh_75]).transpose(),
                   header='Age SFR SFR(25th percentile) SFR(75th percentile)')


    else:
        # mass_ratio = np.log10(np.median(fit.posterior.samples['delayed1:massformed'])**10/
        #                       np.median(fit.posterior.samples['delayed2:massformed'])**10)
        #
        # if (mass_ratio > 0) & (mass_ratio < 2):
        #
        #     print(mass_ratio, 'GRAY LIST!')
        #
        #     # gray_list.append(blob_id)
        #
        #     parametric_plots(fit, blob_id, blob_index, blob_catalog, plot_dir + detection + '_' + config +
        #                      '_extraplanar_mass02')
        #
        #     gc.collect()
        #
        # elif mass_ratio > 2:
        #
        #     print(mass_ratio, 'BLACK LIST!')
        #
        #     # blacklist.append(blob_id)
        #
        #     parametric_plots(fit, blob_id, blob_index, blob_catalog, plot_dir + detection + '_' + config +
        #                      '_extraplanar_mass2')
        #
        #     gc.collect()
        #
        # else:
        #
        #     print(mass_ratio)
        #
        #     parametric_plots(fit, blob_id, blob_index, blob_catalog, plot_dir + detection + '_' + config +
        #                      '_extraplanar')
        #
        #     gc.collect()

        fit.posterior.get_advanced_quantities()
        output_catalog['galaxy'][blob_index] = blob_id.split('_')[0]
        output_catalog['age_old'][blob_index] = np.median(fit.posterior.samples['delayed1:age'])
        output_catalog['tau_old'][blob_index] = np.median(fit.posterior.samples['delayed1:tau'])
        output_catalog['metallicity_old'][blob_index] = np.median(fit.posterior.samples['delayed1:metallicity'])
        output_catalog['formed_mass_old'][blob_index] = np.median(fit.posterior.samples['delayed1:massformed'])
        output_catalog['age_young'][blob_index] = np.median(fit.posterior.samples['delayed2:age'])
        output_catalog['age_young_std'][blob_index] = np.std(fit.posterior.samples['delayed2:age'])
        output_catalog['age_young_iqr'][blob_index] = np.percentile(fit.posterior.samples['delayed2:age'], 75) - \
                                                      np.percentile(fit.posterior.samples['delayed2:age'], 25)
        output_catalog['age_young_gini'][blob_index] = gini(fit.posterior.samples['delayed2:age'])
        output_catalog['tau_young'][blob_index] = np.median(fit.posterior.samples['delayed2:tau'])
        output_catalog['metallicity_young'][blob_index] = np.median(fit.posterior.samples['delayed2:metallicity'])
        output_catalog['formed_mass_young'][blob_index] = np.median(fit.posterior.samples['delayed2:massformed'])
        output_catalog['Av'][blob_index] = np.median(fit.posterior.samples['dust:Av'])
        output_catalog['eta'][blob_index] = np.median(fit.posterior.samples['dust:eta'])
        output_catalog['sfr'][blob_index] = np.median(fit.posterior.samples['sfr'])
        output_catalog['stellar_mass'][blob_index] = np.median(fit.posterior.samples['stellar_mass'])
        output_catalog['mwage'][blob_index] = np.median(fit.posterior.samples['mass_weighted_age'])
        output_catalog['tform'][blob_index] = np.median(fit.posterior.samples['tform'])
        output_catalog['photometry_syn'][blob_index] = np.median(fit.posterior.samples['photometry'], axis=0)
        output_catalog['photometry_obs'][blob_index] = fit.galaxy.photometry[:, 1]
        output_catalog['photometry_nebular'][blob_index] = [synflux(fit.posterior.model_galaxy.wavelengths,
                                                                    fit.posterior.model_galaxy.nebular_spectrum_full,
                                                                    filter_file) for filter_file in filter_files]
        output_catalog['chi2'][blob_index] = np.min(fit.posterior.samples['chisq_phot'])
        output_catalog['Ha'][blob_index] = fit.posterior.model_galaxy.line_fluxes['H  1  6562.81A']
        output_catalog['Nii'][blob_index] = fit.posterior.model_galaxy.line_fluxes['N  2  6583.45A']
        output_catalog['Oiii'][blob_index] = fit.posterior.model_galaxy.line_fluxes['O  3  5006.84A']
        output_catalog['Hb'][blob_index] = fit.posterior.model_galaxy.line_fluxes['H  1  4861.33A']

        age_peaks = find_pdf_peaks(fit.posterior.samples['delayed2:age'], bandwidth=0.1)
        Av_peaks = find_pdf_peaks(fit.posterior.samples['dust:Av'], bandwidth=0.25)
        # mwage_peaks = find_pdf_peaks(fit.posterior.samples['mass_weighted_age'], bandwidth=0.1)

        output_catalog['age_young_peaks'][blob_index] = np.append(age_peaks, np.zeros(15-len(age_peaks)))
        output_catalog['Av_peaks'][blob_index] = np.append(Av_peaks, np.zeros(15-len(Av_peaks)))
        # output_catalog['mwage_peaks'][blob_index] = np.append(mwage_peaks, np.zeros(15-len(mwage_peaks)))
        output_catalog['n_age_young_peaks'][blob_index] = len(age_peaks)
        output_catalog['n_Av_peaks'][blob_index] = len(Av_peaks)
        # output_catalog['n_mwage_peaks'][blob_index] = len(mwage_peaks)

        # for i in range(len(fit.posterior.samples['delayed1:age'])):
        #     all_runs_rows.append([blob_id.split('_')[0], blob_id, i, blob_id+str(i),
        #                           len(fit.posterior.samples['delayed2:age']),
        #                           fit.posterior.samples['delayed2:age'][i],
        #                           fit.posterior.samples['delayed2:tau'][i],
        #                           fit.posterior.samples['delayed2:metallicity'][i],
        #                           fit.posterior.samples['delayed2:massformed'][i],
        #                           fit.posterior.samples['dust:Av'][i],
        #                           fit.posterior.samples['dust:eta'][i],
        #                           fit.posterior.samples['sfr'][i],
        #                           fit.posterior.samples['stellar_mass'][i],
        #                           fit.posterior.samples['mass_weighted_age'][i],
        #                           fit.posterior.samples['tform'][i],
        #                           fit.posterior.samples['chisq_phot'][i]])

if run_flag is False:

    # np.savetxt('graylist', np.array(gray_list).transpose())
    # np.savetxt('blacklist', np.array(blacklist).transpose())

    output_catalog.write(data_dir + 'extraplanar_' + detection + '_' + config + '_bagpipes_results.fits', overwrite=True)

    # all_runs_catalog = Table(rows=all_runs_rows, names=['galaxy', 'fit_id', 'run_index', 'run_id', 'n_samples', 'age',
    #                                                     'tau', 'metallicity', 'massformed', 'Av', 'eta', 'sfr',
    #                                                     'stellar_mass', 'mass_weighted_age', 'tform', 'chisq_phot'])
    # all_runs_catalog.write(data_dir + 'extraplanar_' + detection + '_' + config + '_bagpipes_results_samples.fits',
    #                        overwrite=True)

