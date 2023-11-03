"""

ariel@oapd
29/03/2022

Implements the SINOPSIS way of fitting (12 age bins) in BAGPIPES and hopes for the best

"""

import numpy as np
import matplotlib.pyplot as plt
import bagpipes as pipes
import os
from astropy.table import Table
from starlight_toolkit.synphot import synflux
from toolbox.stat_tools import int_to_bool_list

plt.ioff()

galaxy = 'tail_all31'
detection = 'halpha'
test_id = 'default_nonparametric_hacked'

data_dir = '/Data/'
plot_dir = '/BAGPIPES/plots/'
sfh_dir = '/BAGPIPES/sfh/'

if galaxy + test_id + detection not in os.listdir(plot_dir):
    print(plot_dir + galaxy + test_id + detection, 'not found, creating directory...')
    os.makedirs(plot_dir + galaxy + test_id + detection)
if galaxy + '_' + test_id + '_' + detection not in os.listdir(sfh_dir):
    print(sfh_dir + galaxy + '_' + test_id + '_' + detection, 'not found, creating directory...')
    os.makedirs(sfh_dir + galaxy + '_' + test_id + '_' + detection)

blob_catalog = Table.read(data_dir + galaxy + '_' + detection + '_bagpipes_input.fits')

filter_files = [data_dir + 'filters/HST_WFC3_UVIS2.F275W.dat',
                data_dir + 'filters/HST_WFC3_UVIS2.F336W.dat',
                data_dir + 'filters/HST_WFC3_UVIS2.F606W.dat',
                data_dir + 'filters/HST_WFC3_UVIS2.F680N.dat',
                data_dir + 'filters/HST_WFC3_UVIS2.F814W.dat']

hst_filters = ['F275W', 'F336W', 'F606W', 'F680N', 'F814W']

original_ids = blob_catalog['blob_id']
blob_catalog['blob_id'] = [blob_catalog['blob_id'][i] + '_' + test_id for i in range(len(blob_catalog))]


def load_data(blob_id):

    blob = blob_catalog[blob_catalog['blob_id'] == blob_id]

    fluxes = [blob[hst_filter] for hst_filter in hst_filters]
    errors = [blob['err'+hst_filter] for hst_filter in hst_filters]

    flag = int_to_bool_list(blob['sel_flag'][0])
    for i in range(len(errors)):
        if flag[i] is False:
            errors[i] = errors[i] * 1000

    return np.array([fluxes, errors]).transpose()[0]


output_catalog = Table()
output_catalog['blob_id'] = original_ids
output_catalog['galaxy'] = np.zeros(len(blob_catalog), dtype='S32')
output_catalog['metallicity'] = np.zeros((len(blob_catalog), 12))
output_catalog['formed_mass'] = np.zeros((len(blob_catalog), 12))
output_catalog['sfh'] = np.zeros((len(blob_catalog), 12))
output_catalog['sfh_full'] = np.zeros((len(blob_catalog), 12))
output_catalog['sfh25'] = np.zeros((len(blob_catalog), 12))
output_catalog['sfh75'] = np.zeros((len(blob_catalog), 12))
output_catalog['Av'] = np.zeros(len(blob_catalog))
output_catalog['eta'] = np.zeros(len(blob_catalog))
output_catalog['sfr'] = np.zeros(len(blob_catalog))
output_catalog['stellar_mass'] = np.zeros(len(blob_catalog))
output_catalog['mwage'] = np.zeros(len(blob_catalog))
output_catalog['mwage_alt'] = np.zeros(len(blob_catalog))
output_catalog['mwage_alt_pipes'] = np.zeros(len(blob_catalog))
output_catalog['tform'] = np.zeros(len(blob_catalog))
output_catalog['av_detected'] = np.ones(len(blob_catalog))
output_catalog['logU'] = np.zeros(len(blob_catalog))
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

skip = 0

for blob_id in blob_catalog['blob_id']:

    blob_index = np.argwhere(blob_catalog['blob_id'] == blob_id)[0][0]

    if blob_index < skip:
        continue

    print('>>> Fitting', blob_id, '(', blob_index + 1, 'out of', len(blob_catalog), ')')
    if blob_id+'.h5' in os.listdir('/home/ariel/Workspace/GASP/HST/BAGPIPES/pipes/posterior/'):
        continue
    elif blob_id+'_resume.dat' in os.listdir('/home/ariel/Workspace/GASP/HST/BAGPIPES/pipes/posterior/'):
        continue

    pipes_blob = pipes.galaxy(blob_id, load_data, filt_list=filter_files, spectrum_exists=False, phot_units='ergscma')

    if type(blob_catalog['av_muse'][blob_index]) is not np.float64:
        blob_catalog['av_detected'][blob_index] = False
        continue

    square_bursts = []

    for i in range(len(age_bins)-1):
        constant = {}
        constant["metallicity"] = (0.005, 2.5)
        constant["age_max"] = age_bins[i+1] / 1e9
        constant["age_min"] = age_bins[i] / 1e9
        constant["massformed"] = (0, 10)
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

    output_catalog['galaxy'][blob_index] = blob_id.split('_')[0]
    output_catalog['metallicity'][blob_index] = [np.median(fit.posterior.samples['constant'+str(i)+':metallicity'])
                                                 for i in range(12)]
    output_catalog['formed_mass'][blob_index] = [np.median(fit.posterior.samples['constant'+str(i)+':massformed'])
                                                 for i in range(12)]
    output_catalog['Av'][blob_index] = np.median(fit.posterior.samples['dust:Av'])
    output_catalog['eta'][blob_index] = np.median(fit.posterior.samples['dust:eta'])
    output_catalog['sfr'][blob_index] = np.median(fit.posterior.samples['sfr'])
    output_catalog['stellar_mass'][blob_index] = np.median(fit.posterior.samples['stellar_mass'])
    output_catalog['mwage'][blob_index] = np.median(fit.posterior.samples['mass_weighted_age'])
    output_catalog['tform'][blob_index] = np.median(fit.posterior.samples['tform'])
    # output_catalog['logU'][blob_index] = fit.posterior.fitted_model.model_components['nebular']['logU']
    output_catalog['photometry_syn'][blob_index] = np.median(fit.posterior.samples['photometry'], axis=0)
    output_catalog['photometry_obs'][blob_index] = fit.galaxy.photometry[:, 1]
    output_catalog['photometry_nebular'][blob_index] = [synflux(fit.posterior.model_galaxy.wavelengths,
                                                                fit.posterior.model_galaxy.nebular_spectrum_full,
                                                                filter_file) for filter_file in filter_files]

    print(output_catalog['photometry_nebular'][blob_index])
    print(output_catalog['photometry_syn'][blob_index])
    print(output_catalog['photometry_nebular'][blob_index]/output_catalog['photometry_syn'][blob_index])

    output_catalog['chi2'][blob_index] = np.min(fit.posterior.samples['chisq_phot'])
    output_catalog['Ha'][blob_index] = fit.posterior.model_galaxy.line_fluxes['H  1  6562.81A']
    output_catalog['Nii'][blob_index] = fit.posterior.model_galaxy.line_fluxes['N  2  6583.45A']
    output_catalog['Oiii'][blob_index] = fit.posterior.model_galaxy.line_fluxes['O  3  5006.84A']
    output_catalog['Hb'][blob_index] = fit.posterior.model_galaxy.line_fluxes['H  1  4861.33A']

    output_catalog['sfh'][blob_index] = [np.median(10 ** fit.posterior.samples['constant'+str(i)+':massformed']) /
                                         age_bins_len[i] for i in range(12)]
    output_catalog['sfh25'][blob_index] = [np.percentile(10 ** fit.posterior.samples['constant'+str(i)+':massformed'],
                                                         25) / age_bins_len[i] for i in range(12)]
    output_catalog['sfh75'][blob_index] = [np.percentile(10 ** fit.posterior.samples['constant'+str(i)+':massformed'],
                                                         75) / age_bins_len[i] for i in range(12)]

    mass_formed_lin = np.array([np.median(10 ** fit.posterior.samples['constant' + str(i) + ':massformed']) for i in
                                range(12)])
    mwage = np.sum(age_bins_mid * mass_formed_lin) / np.sum(mass_formed_lin)

    output_catalog['mwage_alt'][blob_index] = mwage

    # Saving SFH
    ages = age_bins_mid
    sfh_median = output_catalog['sfh'][blob_index]
    sfh_25 = output_catalog['sfh25'][blob_index]
    sfh_75 = output_catalog['sfh75'][blob_index]

    np.savetxt(sfh_dir + galaxy + '_' + test_id + '_' + detection + '/sfh_' + blob_id + '.txt',
               np.array([ages, sfh_median, sfh_25, sfh_75]).transpose(),
               header='Age SFR SFR(25th percentile) SFR(75th percentile)')

    # Saving sfh BAGPIPES style (pointwise):
    sfh_median_bp = np.median(fit.posterior.samples['sfh'], axis=0)
    sfh_25_bp = np.percentile(fit.posterior.samples['sfh'], axis=0, q=25)
    sfh_75_bp = np.percentile(fit.posterior.samples['sfh'], axis=0, q=75)
    ages_bp = fit.posterior.sfh.ages

    np.savetxt(sfh_dir + galaxy + '_' + test_id + '_' + detection + '/sfh_' + blob_id + '_bagpipes.txt',
               np.array([ages_bp, sfh_median_bp, sfh_25_bp, sfh_75_bp]).transpose(),
               header='Age SFR SFR(25th percentile) SFR(75th percentile)')

    # Plotting
    split_id = blob_id.split('_')
    plot_id = split_id[0] + '\_' + split_id[1] + '\_' + split_id[2]

    fig = fit.plot_spectrum_posterior(save=False, show=False)
    plt.gca().set_title(plot_id + ', x= %0.2f , y= %0.2f'
                        % (blob_catalog['xc_pix(HST)'][blob_index],
                           blob_catalog['yc_pix(HST)'][blob_index]))
    fig[0].set_size_inches(10, 6)
    plt.savefig('plots/' + galaxy + test_id + detection + '/fit_'+blob_id+'.png', dpi=100)
    try:
        fig = fit.plot_sfh_posterior(save=False, show=False)
        plt.gca().set_title(plot_id + ', x= %0.2f , y= %0.2f'
                            % (blob_catalog['xc_pix(HST)'][blob_index],
                               blob_catalog['yc_pix(HST)'][blob_index]))
        fig[0].set_size_inches(10, 6)
        plt.gca().set_ylim(0, 1.2 * np.max(output_catalog['photometry_syn'][blob_index]))
        plt.savefig('plots/' + galaxy + test_id + detection + '/sfh_'+blob_id+'.png', dpi=100)
    except Exception:
        pass
    # fig = fit.plot_corner(save=False, show=False)
    # fig.suptitle(plot_id + ', x= %0.2f , y= %0.2f'
    #                     % (blob_catalog['xc_pix(HST)'][blob_index],
    #                        blob_catalog['yc_pix(HST)'][blob_index]))
    # plt.savefig('plots/' + galaxy + test_id + detection + '/par_'+blob_id+'.png', dpi=89)

    plt.close('all')

output_catalog.write(data_dir + galaxy + '_' + test_id + '_' + detection + '_bagpipes_results.fits', overwrite=True)
