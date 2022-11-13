"""

ariel@oapd
23/03/2022

Fit HST-defined blobs from Eric, this script should be used specifically for tail clumps. Was a general script but due
to different configurations between tail and disks it was simpler to divide it.

Note that the galaxy class in bagpipes was changed to allow a load_data function that takes as arguments an ID and a
catalog. This allows us to define the load_data function in a separate file, rendering this code more clean.

If run_flag is True, the code will run for all of the blobs that were not started before, so it is possible to run in
parallel. Then in the end of the run I have to delete all remaining inclompete runs using
hst_pipeutils.clear_incomplete_runs() and re run. Then set run_flag=False and fill the table.

"""

import numpy as np
import matplotlib.pyplot as plt
import bagpipes as pipes
import os
from astropy.table import Table
from starlight_toolkit.synphot import synflux
from hst_pipeutils import load_data, parametric_plots
from toolbox.stat_tools import gini, find_maximum, find_pdf_peaks

plt.ioff()

detection = 'f275w'
config = 'dexp_logprior'
run_flag = False

pipes_dir = '/home/ariel/Workspace/GASP/HST/BAGPIPES/'
data_dir = '/home/ariel/Workspace/GASP/HST/Data/'
plot_dir = pipes_dir + 'plots/'
sfh_dir = pipes_dir + 'sfh/'

if detection + '_' + config + '_tail' not in os.listdir(plot_dir):
    os.makedirs(plot_dir + detection + '_' + config + '_tail')
if detection + '_' + config + '_tail' not in os.listdir(sfh_dir):
    os.makedirs(sfh_dir + detection + '_' + config + '_tail')

blob_catalog = Table.read(data_dir+'tail_'+detection+'_bagpipes_input.fits')
blob_catalog = blob_catalog[blob_catalog['sel_flag'] == 31]

filter_files = [data_dir + 'filters/HST_WFC3_UVIS2.F275W.dat',
                data_dir + 'filters/HST_WFC3_UVIS2.F336W.dat',
                data_dir + 'filters/HST_WFC3_UVIS2.F606W.dat',
                data_dir + 'filters/HST_WFC3_UVIS2.F680N.dat',
                data_dir + 'filters/HST_WFC3_UVIS2.F814W.dat']

original_ids = blob_catalog['fit_id']
blob_catalog['fit_id'] = [blob_catalog['fit_id'][i] + '_' + config for i in range(len(blob_catalog))]

if run_flag is False:
    output_catalog = Table()
    output_catalog['blob_id'] = blob_catalog['blob_id']
    output_catalog['fit_id'] = original_ids
    output_catalog['fit_id_fitconfig'] = blob_catalog['fit_id']
    output_catalog['config'] = np.full_like(original_ids, config)
    output_catalog['galaxy'] = np.zeros(len(blob_catalog), dtype='S32')
    output_catalog['sel_flag'] = blob_catalog['sel_flag']
    output_catalog['age'] = np.zeros(len(blob_catalog))
    output_catalog['age_std'] = np.zeros(len(blob_catalog))
    output_catalog['age_iqr'] = np.zeros(len(blob_catalog))
    output_catalog['age_gini'] = np.zeros(len(blob_catalog))
    output_catalog['metallicity'] = np.zeros(len(blob_catalog))
    output_catalog['formed_mass'] = np.zeros(len(blob_catalog))
    output_catalog['Av'] = np.zeros(len(blob_catalog))
    output_catalog['sfr'] = np.zeros(len(blob_catalog))
    output_catalog['stellar_mass'] = np.zeros(len(blob_catalog))
    output_catalog['stellar_mass_std'] = np.zeros(len(blob_catalog))
    output_catalog['stellar_mass_iqr'] = np.zeros(len(blob_catalog))
    output_catalog['stellar_mass_gini'] = np.zeros(len(blob_catalog))
    output_catalog['mwage'] = np.zeros(len(blob_catalog))
    output_catalog['mwage_std'] = np.zeros(len(blob_catalog))
    output_catalog['mwage_iqr'] = np.zeros(len(blob_catalog))
    output_catalog['mwage_gini'] = np.zeros(len(blob_catalog))
    output_catalog['mwage_max'] = np.zeros(len(blob_catalog))
    output_catalog['tform'] = np.zeros(len(blob_catalog))
    output_catalog['photometry_syn'] = np.zeros((len(blob_catalog), 5))
    output_catalog['photometry_obs'] = np.zeros((len(blob_catalog), 5))
    output_catalog['photometry_nebular'] = np.zeros((len(blob_catalog), 5))
    output_catalog['chi2'] = np.zeros(len(blob_catalog))
    output_catalog['Ha'] = np.zeros(len(blob_catalog))
    output_catalog['Nii'] = np.zeros(len(blob_catalog))
    output_catalog['Oiii'] = np.zeros(len(blob_catalog))
    output_catalog['Hb'] = np.zeros(len(blob_catalog))
    output_catalog['age_peaks'] = np.zeros((len(blob_catalog), 15))
    output_catalog['Av_peaks'] = np.zeros((len(blob_catalog), 15))
    output_catalog['mwage_peaks'] = np.zeros((len(blob_catalog), 15))
    output_catalog['n_age_peaks'] = np.zeros((len(blob_catalog)))
    output_catalog['n_Av_peaks'] = np.zeros((len(blob_catalog)))
    output_catalog['n_mwage_peaks'] = np.zeros((len(blob_catalog)))

    if config in ['parametric', 'exp', 'exp_smalltau', 'exp_logprior', 'dexp_logprior']:
        output_catalog['tau'] = np.zeros(len(blob_catalog))
        output_catalog['eta'] = np.zeros(len(blob_catalog))
    if config == 'constant':
        output_catalog['age_min'] = np.zeros(len(blob_catalog))
        output_catalog['age_min_std'] = np.zeros(len(blob_catalog))
        output_catalog['age_min_iqr'] = np.zeros(len(blob_catalog))
        output_catalog['age_min_gini'] = np.zeros(len(blob_catalog))
        output_catalog['age_max'] = np.zeros(len(blob_catalog))
        output_catalog['age_max_std'] = np.zeros(len(blob_catalog))
        output_catalog['age_max_iqr'] = np.zeros(len(blob_catalog))
        output_catalog['age_max_gini'] = np.zeros(len(blob_catalog))
    # all_runs_rows = []

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

    if config in ['parametric', 'exp']:
        exp = {}
        exp["age"] = (0.0, 0.5)
        exp["tau"] = (0.0, 5.)
        exp["massformed"] = (0., 10.)
        exp["metallicity"] = (0.005, 2.5)

    elif config == 'exp_smalltau':
        exp = {}
        exp["age"] = (0.0, 0.5)
        exp["tau"] = (0.0, 0.05)
        exp["massformed"] = (0., 10.)
        exp["metallicity"] = (0.005, 2.5)

    elif config in ['exp_logprior', 'dexp_logprior']:
        exp = {}
        exp["age"] = (0.0, 0.5)
        exp["tau"] = (0.0000000000000000001, 0.5)
        exp["tau_prior"] = 'log_10'
        exp["massformed"] = (0., 10.)
        exp["metallicity"] = (0.005, 2.5)

    elif config == 'ssp':
        burst = {}
        burst['age'] = (0, 0.3)
        burst["massformed"] = (0., 10.)
        burst["metallicity"] = (0.005, 2.5)

    elif config == 'constant':
        constant = {}
        constant["metallicity"] = (0.005, 2.5)
        constant["age_max"] = (0, 0.5)
        constant["age_min"] = (0.0, 0.5)
        constant["massformed"] = (0, 10)

    else:
        print('SFH config not recognized')

    dust = {}
    dust["type"] = "Cardelli"

    dust["Av"] = (0, 1.25)

    if config in ['parametric', 'constant', 'exp', 'exp_logprior', 'dexp_logprior']:
        dust["eta"] = (1, 2.5)

    nebular = {}
    nebular["logU"] = -2.5

    fit_instructions = {}
    fit_instructions["redshift"] = (blob_catalog['galaxy_redshift'][blob_index]-0.0001,
                                    blob_catalog['galaxy_redshift'][blob_index]+0.0001)
    if config in ['parametric', 'dexp_logprior']:
        fit_instructions["delayed"] = exp
    elif config in ['exp', 'exp_smalltau', 'exp_logprior']:
        fit_instructions["exponential"] = exp
    elif config == 'ssp':
        fit_instructions["burst"] = burst
    elif config == 'constant':
        fit_instructions["constant"] = constant
    else:
        print('Ahhmmm.... this is weird')
    fit_instructions["dust"] = dust
    fit_instructions["nebular"] = nebular
    fit_instructions["t_bc"] = 0.02

    fit = pipes.fit(pipes_blob, fit_instructions)

    fit.fit(verbose=False)

    if run_flag:
        print('NOT PLOTTING THIS TIME')

    else:
        # sfh_median = np.median(fit.posterior.samples['sfh'], axis=0)
        # sfh_25 = np.percentile(fit.posterior.samples['sfh'], axis=0, q=25)
        # sfh_75 = np.percentile(fit.posterior.samples['sfh'], axis=0, q=75)
        # ages = fit.posterior.sfh.ages
        #
        # np.savetxt(sfh_dir + detection + '_' + config + '_tail/sfh_' + blob_id + '.txt',
        #            np.array([ages, sfh_median, sfh_25, sfh_75]).transpose(),
        #            header='Age SFR SFR(25th percentile) SFR(75th percentile)')
        #
        parametric_plots(fit, blob_id, blob_index, blob_catalog, plot_dir + detection + '_' + config + '_tail')

#         fit.posterior.get_advanced_quantities()
#         output_catalog['galaxy'][blob_index] = blob_id.split('_')[0]
#
#         if config in ['parametric', 'dexp_logprior', 'exp', 'exp_smalltau', 'exp_logprior']:
#             output_catalog['age'][blob_index] = np.median(fit.posterior.samples['delayed:age'])
#             output_catalog['age_std'][blob_index] = np.std(fit.posterior.samples['delayed:age'])
#             output_catalog['age_iqr'][blob_index] = np.percentile(fit.posterior.samples['delayed:age'], 75) \
#                                                     - np.percentile(fit.posterior.samples['delayed:age'], 25)
#             output_catalog['age_gini'][blob_index] = gini(fit.posterior.samples['delayed:age'])
#             output_catalog['tau'][blob_index] = np.median(fit.posterior.samples['delayed:tau'])
#             output_catalog['metallicity'][blob_index] = np.median(fit.posterior.samples['delayed:metallicity'])
#             output_catalog['formed_mass'][blob_index] = np.median(fit.posterior.samples['delayed:massformed'])
#             output_catalog['eta'][blob_index] = np.median(fit.posterior.samples['dust:eta'])
#
#         elif config == 'ssp':
#             output_catalog['age'][blob_index] = np.median(fit.posterior.samples['burst:age'])
#             output_catalog['age_std'][blob_index] = np.std(fit.posterior.samples['burst:age'])
#             output_catalog['age_iqr'][blob_index] = np.percentile(fit.posterior.samples['burst:age'], 75) - np.percentile(fit.posterior.samples['burst:age'], 25)
#             output_catalog['age_gini'][blob_index] = gini(fit.posterior.samples['burst:age'])
#             output_catalog['metallicity'][blob_index] = np.median(fit.posterior.samples['burst:metallicity'])
#             output_catalog['formed_mass'][blob_index] = np.median(fit.posterior.samples['burst:massformed'])
#
#         elif config == 'constant':
#             output_catalog['age_min'][blob_index] = np.median(fit.posterior.samples['constant:age_min'])
#             output_catalog['age_min_std'][blob_index] = np.std(fit.posterior.samples['constant:age_min'])
#             output_catalog['age_min_iqr'][blob_index] = np.percentile(fit.posterior.samples['constant:age_min'], 75) \
#                                                         - np.percentile(fit.posterior.samples['constant:age_min'], 25)
#             output_catalog['age_min_gini'][blob_index] = gini(fit.posterior.samples['constant:age_min'])
#             output_catalog['age_max'][blob_index] = np.median(fit.posterior.samples['constant:age_max'])
#             output_catalog['age_max_std'][blob_index] = np.std(fit.posterior.samples['constant:age_max'])
#             output_catalog['age_max_iqr'][blob_index] = np.percentile(fit.posterior.samples['constant:age_max'], 75) \
#                                                         - np.percentile(fit.posterior.samples['constant:age_max'], 25)
#             output_catalog['age_max_gini'][blob_index] = gini(fit.posterior.samples['constant:age_max'])
#             output_catalog['metallicity'][blob_index] = np.median(fit.posterior.samples['constant:metallicity'])
#             output_catalog['formed_mass'][blob_index] = np.median(fit.posterior.samples['constant:massformed'])
#
#         output_catalog['Av'][blob_index] = np.median(fit.posterior.samples['dust:Av'])
#         output_catalog['sfr'][blob_index] = np.median(fit.posterior.samples['sfr'])
#         output_catalog['stellar_mass'][blob_index] = np.median(fit.posterior.samples['stellar_mass'])
#         output_catalog['stellar_mass_std'][blob_index] = np.std(fit.posterior.samples['stellar_mass'])
#         output_catalog['stellar_mass_iqr'][blob_index] = np.percentile(fit.posterior.samples['stellar_mass'], 75) - np.percentile(fit.posterior.samples['stellar_mass'], 25)
#         output_catalog['stellar_mass_gini'][blob_index] = gini(fit.posterior.samples['stellar_mass'])
#         output_catalog['mwage'][blob_index] = np.median(fit.posterior.samples['mass_weighted_age'])
#         output_catalog['mwage_std'][blob_index] = np.std(fit.posterior.samples['mass_weighted_age'])
#         output_catalog['mwage_iqr'][blob_index] = np.percentile(fit.posterior.samples['mass_weighted_age'], 75) - np.percentile(fit.posterior.samples['mass_weighted_age'], 25)
#         output_catalog['mwage_gini'][blob_index] = gini(fit.posterior.samples['mass_weighted_age'])
#         output_catalog['mwage_max'][blob_index] = find_maximum(fit.posterior.samples['mass_weighted_age'])
#         output_catalog['tform'][blob_index] = np.median(fit.posterior.samples['tform'])
#         output_catalog['photometry_syn'][blob_index] = np.median(fit.posterior.samples['photometry'], axis=0)
#         output_catalog['photometry_obs'][blob_index] = fit.galaxy.photometry[:, 1]
#         output_catalog['photometry_nebular'][blob_index] = [synflux(fit.posterior.model_galaxy.wavelengths,
#                                                                     fit.posterior.model_galaxy.nebular_spectrum_full,
#                                                                     filter_file) for filter_file in filter_files]
#         output_catalog['chi2'][blob_index] = np.min(fit.posterior.samples['chisq_phot'])
#         output_catalog['Ha'][blob_index] = fit.posterior.model_galaxy.line_fluxes['H  1  6562.81A']
#         output_catalog['Nii'][blob_index] = fit.posterior.model_galaxy.line_fluxes['N  2  6583.45A']
#         output_catalog['Oiii'][blob_index] = fit.posterior.model_galaxy.line_fluxes['O  3  5006.84A']
#         output_catalog['Hb'][blob_index] = fit.posterior.model_galaxy.line_fluxes['H  1  4861.33A']
#
#         if config in ['parametric', 'dexp_logprior']:
#             age_peaks = find_pdf_peaks(fit.posterior.samples['delayed:age'], bandwidth=0.1)
#         elif config == 'ssp':
#             age_peaks = find_pdf_peaks(fit.posterior.samples['burst:age'], bandwidth=0.1)
#         elif config in ['exp', 'exp_smalltau', 'exp_logprior']:
#             age_peaks = find_pdf_peaks(fit.posterior.samples['exponential:age'], bandwidth=0.1)
#         elif config == 'constant':
#             age_peaks = find_pdf_peaks(fit.posterior.samples['constant:age_max'], bandwidth=0.1)
#         Av_peaks = find_pdf_peaks(fit.posterior.samples['dust:Av'], bandwidth=0.25)
#         mwage_peaks = find_pdf_peaks(fit.posterior.samples['mass_weighted_age'], bandwidth=0.1)
#
#         output_catalog['age_peaks'][blob_index] = np.append(age_peaks, np.zeros(15-len(age_peaks)))
#         output_catalog['Av_peaks'][blob_index] = np.append(Av_peaks, np.zeros(15-len(Av_peaks)))
#         output_catalog['mwage_peaks'][blob_index] = np.append(mwage_peaks, np.zeros(15-len(mwage_peaks)))
#         output_catalog['n_age_peaks'][blob_index] = len(age_peaks)
#         output_catalog['n_Av_peaks'][blob_index] = len(Av_peaks)
#         output_catalog['n_mwage_peaks'][blob_index] = len(mwage_peaks)
#
#         # for i in range(len(fit.posterior.samples['delayed:age'])):
#         #     all_runs_rows.append([blob_id.split('_')[0], blob_id, i, blob_id+str(i),
#         #                           len(fit.posterior.samples['delayed:age']),
#         #                           fit.posterior.samples['delayed:age'][i],
#         #                           fit.posterior.samples['delayed:tau'][i],
#         #                           fit.posterior.samples['delayed:metallicity'][i],
#         #                           fit.posterior.samples['delayed:massformed'][i],
#         #                           fit.posterior.samples['dust:Av'][i],
#         #                           fit.posterior.samples['dust:eta'][i],
#         #                           fit.posterior.samples['sfr'][i],
#         #                           fit.posterior.samples['stellar_mass'][i],
#         #                           fit.posterior.samples['mass_weighted_age'][i],
#         #                           fit.posterior.samples['tform'][i],
#         #                           fit.posterior.samples['chisq_phot'][i]])
#
#
# if run_flag is False:
#
#     if detection == 'optical_only':
#         output_catalog['match_id'] = [blob_catalog['blob_id'][i].split('_')[0] + '_' +
#                                       blob_catalog['blob_id'][i].split('_')[1] + '_f606w' for i in
#                                       range(len(blob_catalog))]
#
#     output_catalog.write(data_dir + 'tail_' + detection + '_' + config + '_bagpipes_results.fits', overwrite=True)
#
#     # all_runs_catalog = Table(rows=all_runs_rows, names=['galaxy', 'fit_id', 'run_index', 'run_id', 'n_samples', 'age',
#     #                                                     'tau', 'metallicity', 'massformed', 'Av', 'eta', 'sfr',
#     #                                                     'stellar_mass', 'mass_weighted_age', 'tform', 'chisq_phot'])
#     # all_runs_catalog.write(data_dir + 'tail_' + detection + '_' + config + '_bagpipes_results_samples.fits',
#     #                        overwrite=True)



















# if blob_catalog['out_of_muse'][blob_index]:
#     print(blob_id + ' is out of MUSE field, using generic dust prior')
#     dust["Av"] = (0, 1.25)
# elif blob_catalog['av_muse'][blob_index] == 0:
#     print(blob_id + ' has no dust detected on MUSE, setting to unconstrained')
#     dust["Av"] = (0, 1.25)
# elif (blob_catalog['av_muse'][blob_index] != 0) & (blob_catalog['av_muse_std'][blob_index] == 0):
#     print(blob_id + ' has 0 Av std, using generic std')
#     dust["Av"] = (blob_catalog['av_muse'][blob_index]/2 - 0.5, blob_catalog['av_muse'][blob_index]/2 + 0.5)
# else:
#     dust["Av"] = (blob_catalog['av_muse'][blob_index]/2 - 3 * blob_catalog['av_muse_std'][blob_index],
#                   blob_catalog['av_muse'][blob_index]/2 + 3 * blob_catalog['av_muse_std'][blob_index])
#     dust["Av_prior"] = "Gaussian"
#     dust["Av_prior_mu"] = blob_catalog['av_muse'][blob_index]/2
#     dust["Av_prior_sigma"] = blob_catalog['av_muse_std'][blob_index]
