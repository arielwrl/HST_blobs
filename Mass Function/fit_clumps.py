"""

ariel@oapd
11/01/2023

Fits clumps, _single in the end of the config variable for single component, _double for double component.

"""

import numpy as np
import matplotlib.pyplot as plt
import bagpipes as pipes
import os
from astropy.table import Table
from hst_pipeutils import load_data
from toolbox.stat_tools import gini, find_maximum, find_pdf_peaks
from pydrive.auth import GoogleAuth 
from pydrive.drive import GoogleDrive
import glob

plt.ioff()

gauth = GoogleAuth()
gauth.LocalWebserverAuth()

drive = GoogleDrive(gauth)

detection = 'f275w'
config = 'dexp_logpriorteeeeeeTEEEEESTTTSTST_double'
run_flag = True
region = 'disk'

pipes_dir = '/home/ariel/Workspace/Mass Function/'
data_dir = '/home/ariel/Workspace/Mass Function/Data/'

clump_catalog = Table.read(data_dir + detection + '_' + region + '_bagpipes_input.fits')
clump_catalog['fit_id'] = [clump_catalog['clump_id'][i] + '_' + config for i in range(len(clump_catalog))]

filter_files = [data_dir + 'filters/HST_WFC3_UVIS2.F275W.dat',
                data_dir + 'filters/HST_WFC3_UVIS2.F336W.dat',
                data_dir + 'filters/HST_WFC3_UVIS2.F606W.dat',
                data_dir + 'filters/HST_WFC3_UVIS2.F680N.dat',
                data_dir + 'filters/HST_WFC3_UVIS2.F814W.dat']

if run_flag is False:
    output_catalog = Table()
    output_catalog['clump_id'] = clump_catalog['clump_id']
    output_catalog['fit_id'] = clump_catalog['fit_id']
    output_catalog['config'] = np.full_like(clump_catalog['fit_id'], config)
    output_catalog['galaxy'] = np.zeros(len(clump_catalog), dtype='S32')
    output_catalog['sel_flag'] = clump_catalog['sel_flag']
    output_catalog['age'] = np.zeros(len(clump_catalog))
    output_catalog['age_std'] = np.zeros(len(clump_catalog))
    output_catalog['age_iqr'] = np.zeros(len(clump_catalog))
    output_catalog['age_gini'] = np.zeros(len(clump_catalog))
    output_catalog['tau'] = np.zeros(len(clump_catalog))
    output_catalog['tau_std'] = np.zeros(len(clump_catalog))
    output_catalog['tau_iqr'] = np.zeros(len(clump_catalog))
    output_catalog['metallicity'] = np.zeros(len(clump_catalog))
    output_catalog['metallicity_iqr'] = np.zeros(len(clump_catalog))
    output_catalog['formed_mass'] = np.zeros(len(clump_catalog))
    output_catalog['formed_mass_iqr'] = np.zeros(len(clump_catalog))
    output_catalog['logU'] = np.zeros(len(clump_catalog))
    output_catalog['logU_iqr'] = np.zeros(len(clump_catalog))

    if config.split('_')[-1] == 'single':
        output_catalog['age_peaks'] = np.zeros((len(clump_catalog), 15))
        output_catalog['Av_peaks'] = np.zeros((len(clump_catalog), 15))
        output_catalog['mwage_peaks'] = np.zeros((len(clump_catalog), 15))
        output_catalog['n_age_peaks'] = np.zeros((len(clump_catalog)))
        output_catalog['n_Av_peaks'] = np.zeros((len(clump_catalog)))
        output_catalog['n_mwage_peaks'] = np.zeros((len(clump_catalog)))

    elif config.split('_')[-1] == 'double':
        output_catalog['age_old'] = np.zeros(len(clump_catalog))
        output_catalog['tau_old'] = np.zeros(len(clump_catalog))
        output_catalog['metallicity'] = np.zeros(len(clump_catalog))
        output_catalog['metallicity_old'] = np.zeros(len(clump_catalog))
        output_catalog['formed_mass_old'] = np.zeros(len(clump_catalog))
        output_catalog['mass_ratio'] = np.zeros(len(clump_catalog))

    else:
        print('Dude you need more coffee.')

    output_catalog['Av'] = np.zeros(len(clump_catalog))
    output_catalog['Av_iqr'] = np.zeros(len(clump_catalog))
    output_catalog['eta'] = np.zeros(len(clump_catalog))
    output_catalog['eta_std'] = np.zeros(len(clump_catalog))
    output_catalog['eta_iqr'] = np.zeros(len(clump_catalog))
    output_catalog['sfr'] = np.zeros(len(clump_catalog))
    output_catalog['isfr'] = np.zeros(len(clump_catalog))
    output_catalog['max_sfr'] = np.zeros(len(clump_catalog))
    output_catalog['t_max_sfr'] = np.zeros(len(clump_catalog))
    output_catalog['stellar_mass'] = np.zeros(len(clump_catalog))
    output_catalog['stellar_mass_std'] = np.zeros(len(clump_catalog))
    output_catalog['stellar_mass_iqr'] = np.zeros(len(clump_catalog))
    output_catalog['stellar_mass_gini'] = np.zeros(len(clump_catalog))
    output_catalog['mwage'] = np.zeros(len(clump_catalog))
    output_catalog['mwage_std'] = np.zeros(len(clump_catalog))
    output_catalog['mwage_iqr'] = np.zeros(len(clump_catalog))
    output_catalog['mwage_gini'] = np.zeros(len(clump_catalog))
    output_catalog['mwage_max'] = np.zeros(len(clump_catalog))
    output_catalog['tform'] = np.zeros(len(clump_catalog))
    output_catalog['photometry_syn'] = np.zeros((len(clump_catalog), 5))
    output_catalog['photometry_syn_iqr'] = np.zeros((len(clump_catalog), 5))
    output_catalog['photometry_obs'] = np.zeros((len(clump_catalog), 5))
    output_catalog['chi2'] = np.zeros(len(clump_catalog))
    output_catalog['Ha'] = np.zeros(len(clump_catalog))
    output_catalog['Nii'] = np.zeros(len(clump_catalog))
    output_catalog['Nii_6583'] = np.zeros(len(clump_catalog))
    output_catalog['Nii_6548'] = np.zeros(len(clump_catalog))
    output_catalog['Oiii'] = np.zeros(len(clump_catalog))
    output_catalog['Hb'] = np.zeros(len(clump_catalog))

skip = 0

for clump_id in clump_catalog['fit_id']:

    clump_index = np.argwhere(clump_catalog['fit_id'] == clump_id)[0][0]

    if clump_index < skip:
        continue

    print('>>> Fitting', clump_id, '(', clump_index + 1, 'out of', len(clump_catalog), ')')
    file_list = os.listdir(pipes_dir + 'pipes/posterior/')
    if ((clump_id+'.h5' in file_list) or (clump_id+'_resume.dat' in file_list)) and run_flag is True:
        continue

    double_runs = drive.CreateFile({'id': '1fqq4n3J99VvRvQbIMw-d9SfszXJqyaYa'})
    double_runs_string = double_runs.GetContentString()
    double_runs_list = double_runs_string.split("\n")
    if clump_id + '.h5' in double_runs_list:
    	print("FILE IN GDRIVE LIST")
    else:
        double_runs_string += clump_id + '.h5\n'
        file = drive.CreateFile({'title': 'double_runs', 'parents' : [{'id' : '1Z1-ILMs-Z94gqODvfLlvLWWPA_mh4Xd_'}]})  
    file.SetContentString(double_runs_string)
    file.Upload()	

    pipes_clump = pipes.galaxy(clump_id, clump_catalog, load_data, filt_list=filter_files, spectrum_exists=False,
                               phot_units='ergscma', spec_wavs=np.arange(2400, 8100, 1))

    dexp = {}
    dexp["age"] = (0.0, 0.5)
    dexp["tau"] = (0.0001, 0.5)
    dexp["tau_prior"] = 'log_10'
    dexp["massformed"] = (0., 10.)
    dexp["metallicity"] = (0.25, 1.75)
    dexp["metallicity_prior"] = 'Gaussian'
    dexp["metallicity_prior_mu"] = 1.0
    dexp["metallicity_prior_sigma"] = 0.25

    if config.split('_')[-1] == 'double':
        dexp_old = {}
        dexp_old["age"] = (5.0, 14.)
        dexp_old["tau"] = (0.1, 14.)
        dexp_old["massformed"] = (0., 12.)
        dexp_old["metallicity"] = (0.25, 1.75)
        dexp_old["metallicity_prior"] = 'Gaussian'
        dexp_old["metallicity_prior_mu"] = 1.0
        dexp_old["metallicity_prior_sigma"] = 0.25

    dust = {}
    dust["type"] = "Cardelli"

    dust["Av"] = (0, 1.25)
    dust["eta"] = (1, 2.5)

    nebular = {}
    nebular["logU"] = (-3.5, -2)

    fit_instructions = {}
    fit_instructions["redshift"] = clump_catalog['galaxy_redshift'][clump_index]

    if config.split('_')[-1] == 'single':
        fit_instructions["delayed"] = dexp
    if config.split('_')[-1] == 'double':
        fit_instructions["delayed1"] = dexp
        fit_instructions["delayed2"] = dexp_old

    fit_instructions["dust"] = dust
    fit_instructions["nebular"] = nebular
    fit_instructions["t_bc"] = 0.02

    fit = pipes.fit(pipes_clump, fit_instructions)

    fit.fit(verbose=False)

    if run_flag is False:

        fit.posterior.get_advanced_quantities()

        # wl = fit.posterior.model_galaxy.wavelengths

        # flux_med = np.median(fit.posterior.samples['spectrum_full'], axis=0) / 1e-18
        # flux_25 = np.percentile(fit.posterior.samples['spectrum_full'], 25, axis=0) / 1e-18
        # flux_75 = np.percentile(fit.posterior.samples['spectrum_full'], 75, axis=0) / 1e-18

        # np.savetxt('/home/ariel/Workspace/GASP/HST/Data/full_spectra/' + clump_id + '.dat',
        #            np.array([wl, flux_med, flux_25, flux_75]).transpose())

        output_catalog['galaxy'][clump_index] = clump_id.split('_')[0]

        if config.split('_')[-1] == 'single':
            output_catalog['age'][clump_index] = np.median(fit.posterior.samples['delayed:age'])
            output_catalog['age_std'][clump_index] = np.std(fit.posterior.samples['delayed:age'])
            output_catalog['age_iqr'][clump_index] = np.percentile(fit.posterior.samples['delayed:age'], 75) \
                                                     - np.percentile(fit.posterior.samples['delayed:age'], 25)
            output_catalog['age_gini'][clump_index] = gini(fit.posterior.samples['delayed:age'])
            output_catalog['tau'][clump_index] = np.median(fit.posterior.samples['delayed:tau'])
            output_catalog['tau_std'][clump_index] = np.std(fit.posterior.samples['delayed:tau'])
            output_catalog['tau_iqr'][clump_index] = np.percentile(fit.posterior.samples['delayed:tau'], 75) \
                                                     - np.percentile(fit.posterior.samples['delayed:tau'], 25)
            output_catalog['metallicity'][clump_index] = np.median(fit.posterior.samples['delayed:metallicity'])
            output_catalog['metallicity_iqr'][clump_index] = np.percentile(fit.posterior.samples['delayed:metallicity'], 75) \
                                                             - np.percentile(fit.posterior.samples['delayed:metallicity'], 25)
            output_catalog['formed_mass'][clump_index] = np.median(fit.posterior.samples['delayed:massformed'])
            output_catalog['formed_mass_iqr'][clump_index] = np.percentile(fit.posterior.samples['delayed:massformed'], 75) \
                                                             - np.percentile(fit.posterior.samples['delayed:massformed'], 25)
            output_catalog['logU'][clump_index] = np.median(fit.posterior.samples['nebular:logU'])
            output_catalog['logU_iqr'][clump_index] = np.percentile(fit.posterior.samples['nebular:logU'],
                                                                           75) \
                                                             - np.percentile(fit.posterior.samples['nebular:logU'], 25)

            age_peaks = find_pdf_peaks(fit.posterior.samples['delayed:age'], bandwidth=0.1)
            Av_peaks = find_pdf_peaks(fit.posterior.samples['dust:Av'], bandwidth=0.25)
            mwage_peaks = find_pdf_peaks(fit.posterior.samples['mass_weighted_age'], bandwidth=0.1)

            output_catalog['age_peaks'][clump_index] = np.append(age_peaks, np.zeros(15 - len(age_peaks)))
            output_catalog['Av_peaks'][clump_index] = np.append(Av_peaks, np.zeros(15 - len(Av_peaks)))
            output_catalog['mwage_peaks'][clump_index] = np.append(mwage_peaks, np.zeros(15 - len(mwage_peaks)))
            output_catalog['n_age_peaks'][clump_index] = len(age_peaks)
            output_catalog['n_Av_peaks'][clump_index] = len(Av_peaks)
            output_catalog['n_mwage_peaks'][clump_index] = len(mwage_peaks)

        elif config.split('_')[-1] == 'double':
            output_catalog['age'][clump_index] = np.median(fit.posterior.samples['delayed1:age'])
            output_catalog['age_std'][clump_index] = np.std(fit.posterior.samples['delayed1:age'])
            output_catalog['age_iqr'][clump_index] = np.percentile(fit.posterior.samples['delayed1:age'], 75) \
                                                     - np.percentile(fit.posterior.samples['delayed1:age'], 25)
            output_catalog['age_gini'][clump_index] = gini(fit.posterior.samples['delayed1:age'])
            output_catalog['age_old'][clump_index] = np.median(fit.posterior.samples['delayed2:age'])
            output_catalog['tau'][clump_index] = np.median(fit.posterior.samples['delayed1:tau'])
            output_catalog['tau_std'][clump_index] = np.std(fit.posterior.samples['delayed1:tau'])
            output_catalog['tau_iqr'][clump_index] = np.percentile(fit.posterior.samples['delayed1:tau'], 75) \
                                                     - np.percentile(fit.posterior.samples['delayed1:tau'], 25)
            output_catalog['tau_old'][clump_index] = np.median(fit.posterior.samples['delayed2:tau'])
            output_catalog['metallicity'][clump_index] = np.median(fit.posterior.samples['delayed1:metallicity'])
            output_catalog['metallicity_iqr'][clump_index] = np.percentile(fit.posterior.samples['delayed1:metallicity'],
                                                                           75) - np.percentile(fit.posterior.samples['delayed1:metallicity'], 25)
            output_catalog['formed_mass'][clump_index] = np.median(fit.posterior.samples['delayed1:massformed'])
            output_catalog['formed_mass_iqr'][clump_index] = np.percentile(fit.posterior.samples['delayed1:massformed'],
                                                                           75) - np.percentile(fit.posterior.samples['delayed1:massformed'], 25)
            output_catalog['metallicity_old'][clump_index] = np.median(fit.posterior.samples['delayed2:metallicity'])
            output_catalog['formed_mass_old'][clump_index] = np.median(fit.posterior.samples['delayed2:massformed'])
            output_catalog['mass_ratio'][clump_index] = np.log10((output_catalog['formed_mass_old'][clump_index]**10) /
                                                                 (output_catalog['formed_mass'][clump_index]**10))

        output_catalog['eta'][clump_index] = np.median(fit.posterior.samples['dust:eta'])
        output_catalog['eta_std'][clump_index] = np.std(fit.posterior.samples['dust:eta'])
        output_catalog['eta_iqr'][clump_index] = np.percentile(fit.posterior.samples['dust:eta'], 75) \
                                                - np.percentile(fit.posterior.samples['dust:eta'], 25)
        output_catalog['Av'][clump_index] = np.median(fit.posterior.samples['dust:Av'])
        output_catalog['Av_iqr'][clump_index] = np.percentile(fit.posterior.samples['dust:Av'], 75) \
                                                 - np.percentile(fit.posterior.samples['dust:Av'], 25)
        output_catalog['sfr'][clump_index] = np.median(fit.posterior.samples['sfr'])

        sfh = np.median(fit.posterior.samples['sfh'], axis=0)
        ages = fit.posterior.sfh.ages
        output_catalog['isfr'][clump_index] = sfh[0]
        output_catalog['max_sfr'][clump_index] = np.max(sfh)
        output_catalog['t_max_sfr'][clump_index] = ages[np.argmax(sfh)]

        output_catalog['stellar_mass'][clump_index] = np.median(fit.posterior.samples['stellar_mass'])
        output_catalog['stellar_mass_std'][clump_index] = np.std(fit.posterior.samples['stellar_mass'])
        output_catalog['stellar_mass_iqr'][clump_index] = np.percentile(fit.posterior.samples['stellar_mass'], 75) - np.percentile(fit.posterior.samples['stellar_mass'], 25)
        output_catalog['stellar_mass_gini'][clump_index] = gini(fit.posterior.samples['stellar_mass'])
        output_catalog['mwage'][clump_index] = np.median(fit.posterior.samples['mass_weighted_age'])
        output_catalog['mwage_std'][clump_index] = np.std(fit.posterior.samples['mass_weighted_age'])
        output_catalog['mwage_iqr'][clump_index] = np.percentile(fit.posterior.samples['mass_weighted_age'], 75) - np.percentile(fit.posterior.samples['mass_weighted_age'], 25)
        output_catalog['mwage_gini'][clump_index] = gini(fit.posterior.samples['mass_weighted_age'])
        output_catalog['mwage_max'][clump_index] = find_maximum(fit.posterior.samples['mass_weighted_age'])
        output_catalog['tform'][clump_index] = np.median(fit.posterior.samples['tform'])
        output_catalog['photometry_syn'][clump_index] = np.median(fit.posterior.samples['photometry'], axis=0)
        output_catalog['photometry_syn_iqr'][clump_index] = np.percentile(fit.posterior.samples['photometry'], 75, axis=0) - np.percentile(fit.posterior.samples['photometry'], 25, axis=0)
        output_catalog['photometry_obs'][clump_index] = fit.galaxy.photometry[:, 1]
        output_catalog['chi2'][clump_index] = np.min(fit.posterior.samples['chisq_phot'])
        output_catalog['Ha'][clump_index] = fit.posterior.model_galaxy.line_fluxes['H  1  6562.81A']
        output_catalog['Nii_6583'][clump_index] = fit.posterior.model_galaxy.line_fluxes['N  2  6583.45A']
        output_catalog['Nii_6548'][clump_index] = fit.posterior.model_galaxy.line_fluxes['N  2  6548.05A']
        output_catalog['Oiii'][clump_index] = fit.posterior.model_galaxy.line_fluxes['O  3  5006.84A']
        output_catalog['Hb'][clump_index] = fit.posterior.model_galaxy.line_fluxes['H  1  4861.33A']

if run_flag is False:
    output_catalog.write(data_dir + detection + '_' + region + '_' + config + '_bagpipes_results.fits', overwrite=True)

