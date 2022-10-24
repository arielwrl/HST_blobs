"""

ariel@oapd
05/05/2022

Dumping some functions here to clean the rest of the code

"""

import numpy as np
import matplotlib.pyplot as plt
import os
from toolbox.stat_tools import int_to_bool_list
import seaborn as sns
import bagpipes as pipes

hst_filters = ['F275W', 'F336W', 'F606W', 'F680N', 'F814W']

data_dir = '/home/ariel/Workspace/GASP/HST/Data/'

filter_files = [data_dir + 'filters/HST_WFC3_UVIS2.F275W.dat',
                data_dir + 'filters/HST_WFC3_UVIS2.F336W.dat',
                data_dir + 'filters/HST_WFC3_UVIS2.F606W.dat',
                data_dir + 'filters/HST_WFC3_UVIS2.F680N.dat',
                data_dir + 'filters/HST_WFC3_UVIS2.F814W.dat']


def smart_remove_file(file_name):
    try:
        os.remove(file_name)
    except FileNotFoundError:
        print('Did not find', file_name)


def clear_incomplete_runs(pipes_directory='./'):

    file_list = os.listdir(pipes_directory + 'pipes/posterior/')

    for file in file_list:
        if file.split('.')[1] in ['points', 'dat', 'txt']:
            smart_remove_file(pipes_directory + 'pipes/posterior/' + file)

    print('>>> Removed all files!')


def load_data(blob_id, blob_catalog):

    blob = blob_catalog[blob_catalog['fit_id'] == blob_id]

    fluxes = [blob[hst_filter] for hst_filter in hst_filters]
    errors = [blob['err'+hst_filter] for hst_filter in hst_filters]

    flag = int_to_bool_list(blob['sel_flag'][0])
    for i in range(len(errors)):
        if flag[i] is False:
            errors[i] = errors[i] * 1000

    return np.array([fluxes, errors]).transpose()[0]


def parametric_plots(fit, blob_id, blob_index, blob_catalog, plot_dir):

    split_id = blob_id.split('_')
    plot_id = split_id[0] + '\_' + split_id[1] + '\_' + split_id[2]

    fig = fit.plot_spectrum_posterior(save=False, show=False)
    plt.gca().set_title(plot_id + ', x= %0.2f , y= %0.2f'
                        % (blob_catalog['xc_pix(HST)'][blob_index],
                           blob_catalog['yc_pix(HST)'][blob_index]))
    fig[0].set_size_inches(10, 6)
    plt.savefig(plot_dir + '/fit_'+blob_id+'.png', dpi=100)
    try:
        fig = fit.plot_sfh_posterior(save=False, show=False)
        plt.gca().set_title(plot_id + ', x= %0.2f , y= %0.2f'
                            % (blob_catalog['xc_pix(HST)'][blob_index],
                               blob_catalog['yc_pix(HST)'][blob_index]))
        fig[0].set_size_inches(10, 6)
        plt.savefig(plot_dir + '/sfh_' + blob_id + '.png', dpi=100)
    except Exception:
        pass
    fig = fit.plot_corner(save=False, show=False)
    fig.suptitle(plot_id + ', x= %0.2f , y= %0.2f'
                        % (blob_catalog['xc_pix(HST)'][blob_index],
                           blob_catalog['yc_pix(HST)'][blob_index]))
    plt.savefig(plot_dir + '/par_' + blob_id + '.png', dpi=89)

    plt.close('all')


def nonparametric_plots(fit, blob_id, blob_index, blob_catalog, plot_dir):

    split_id = blob_id.split('_')
    plot_id = split_id[0] + '\_' + split_id[1] + '\_' + split_id[2]

    fig = fit.plot_spectrum_posterior(save=False, show=False)
    plt.gca().set_title(plot_id + ', x= %0.2f , y= %0.2f'
                        % (blob_catalog['xc_pix(HST)'][blob_index],
                           blob_catalog['yc_pix(HST)'][blob_index]))
    fig[0].set_size_inches(10, 6)
    plt.savefig(plot_dir + '/fit_'+blob_id+'.png', dpi=100)
    try:
        fig = fit.plot_sfh_posterior(save=False, show=False)
        plt.gca().set_title(plot_id + ', x= %0.2f , y= %0.2f'
                            % (blob_catalog['xc_pix(HST)'][blob_index],
                               blob_catalog['yc_pix(HST)'][blob_index]))
        fig[0].set_size_inches(10, 6)
        plt.savefig(plot_dir + '/sfh_' + blob_id + '.png', dpi=100)
    except Exception:
        pass

    plt.close('all')


def load_data_models(id, model_dict):
    model = model_dict[id]
    return np.array([model.photometry, 0.05 * model.photometry]).transpose()


def plot_example(blob_id, input_catalog, ax=None, plot_distributions=True, spec_color='#758BFD', phot_color='#FF8600'):

    if ax is None:
        ax = plt.gca()

    pipes_blob = pipes.galaxy(blob_id, input_catalog, load_data, filt_list=filter_files, spectrum_exists=False,
                              phot_units='ergscma', spec_wavs=np.arange(2400, 8100, 1))
    blob = pipes.fit(pipes_blob, {})
    blob.posterior.get_advanced_quantities()

    z = blob.posterior.fitted_model.model_components['redshift']

    photometry = blob.posterior.galaxy.photometry
    photometry[:, 0] /= (1 + z)
    photometry[:, 1] *= (1 + z)/1e-18

    wl = blob.posterior.model_galaxy.wavelengths
    wl_highres = np.arange(2450, 9100, 1)
    spec_range = (wl > 2450) & (wl < 9100)

    phot_med = np.median(blob.posterior.samples['photometry'] * (1 + z), axis=0)/1e-18
    phot_25 = np.percentile(blob.posterior.samples['photometry'] * (1 + z), 25, axis=0)/1e-18
    phot_75 = np.percentile(blob.posterior.samples['photometry'] * (1 + z), 75, axis=0)/1e-18

    flux_med = np.median(blob.posterior.samples['spectrum_full'], axis=0)/1e-18
    flux_25 = np.percentile(blob.posterior.samples['spectrum_full'], 25, axis=0)/1e-18
    flux_75 = np.percentile(blob.posterior.samples['spectrum_full'], 75, axis=0)/1e-18

    ages = blob.posterior.sfh.ages / 1e6

    sfh_med = np.median(blob.posterior.samples['sfh'], axis=0)
    sfh_25 = np.percentile(blob.posterior.samples['sfh'], 25, axis=0)
    sfh_75 = np.percentile(blob.posterior.samples['sfh'], 75, axis=0)

    ax.plot(wl, flux_med, color=spec_color)
    ax.fill_between(wl, flux_25, flux_75, alpha=0.25, color=spec_color)
    ax.set_xlim(2450, 8550)
    ax.set_ylim(0, flux_med[spec_range][0])

    ax.bar(x=photometry[:, 0], height=phot_75 - phot_25, bottom=phot_25, width=150, zorder=10, color=phot_color,
           **{'edgecolor': phot_color, 'alpha': 0.8})
    ax.errorbar(photometry[:, 0], photometry[:, 1], photometry[:, 2], barsabove=True, fmt='ok', elinewidth=1,
                markersize=5,
                ecolor='k', zorder=15)
    ax.set_ylabel(r'$F_\lambda\,\mathrm{[10^{-18}\,erg\,s^{-1}\,cm^{-2}\,\AA^{-1}]}$', fontsize=20)
    ax.set_xlabel(r'$\lambda\,\mathrm{[\AA]}$', fontsize=20)

    if plot_distributions:
        # inset = ax.inset_axes([0.75, 0.79, 0.245, 0.20])
        # inset_age = ax.inset_axes([0.75, 0.535, 0.115, 0.18])
        # inset_mass = ax.inset_axes([0.88, 0.535, 0.115, 0.18])

        inset_age = ax.inset_axes([0.75, 0.65, 0.115, 0.3])
        inset_mass = ax.inset_axes([0.88, 0.65, 0.115, 0.3])

        sns.despine(ax=inset_age, left=True)
        sns.despine(ax=inset_mass, left=True)

        inset_age.set_yticks([])
        inset_mass.set_yticks([])

        # inset.plot(ages, sfh_med, color='#333891')
        # inset.fill_between(ages, sfh_25, sfh_75, alpha=0.25, color='#333891')
        # inset.set_xlim(ages[0], 1.25 * np.median(blob.posterior.samples['delayed:age'] * 1e3))
        # inset.set_ylim(0, sfh_75[0])
        # inset.set_xlabel(r'$t\,\mathrm{[Myr]}$', fontsize=13)
        # inset.set_ylabel(r'SFR\,$\,\mathrm{[M_\odot/yr]}$', fontsize=13)

        sns.kdeplot(1e3 * blob.posterior.samples['delayed:age'], ax=inset_age, color=spec_color, fill=True, alpha=0.5)
        sns.kdeplot(blob.posterior.samples['stellar_mass'], ax=inset_mass, color=spec_color, fill=True, alpha=0.5)

        # inset_age.set_ylabel(r'$P$', fontsize=13)
        inset_mass.set_ylabel(r'')
        inset_age.set_xlabel(r'$t_\star\mathrm{[Myr]}$', fontsize=13)
        inset_mass.set_xlabel(r'$\log\,M/M_\odot$', fontsize=13)


def plot_example_disk(blob_id, input_catalog, ax=None, plot_distributions=True, spec_color='#758BFD', phot_color='#FF8600'):

    if ax is None:
        ax = plt.gca()

    pipes_blob = pipes.galaxy(blob_id, input_catalog, load_data, filt_list=filter_files, spectrum_exists=False,
                              phot_units='ergscma', spec_wavs=np.arange(2400, 8100, 1))
    blob = pipes.fit(pipes_blob, {})
    blob.posterior.get_advanced_quantities()

    z = blob.posterior.fitted_model.model_components['redshift']

    photometry = blob.posterior.galaxy.photometry
    photometry[:, 0] /= (1 + z)
    photometry[:, 1] *= (1 + z)/1e-18

    wl = blob.posterior.model_galaxy.wavelengths
    wl_highres = np.arange(2450, 9100, 1)
    spec_range = (wl > 2450) & (wl < 9100)

    phot_med = np.median(blob.posterior.samples['photometry'] * (1 + z), axis=0)/1e-18
    phot_25 = np.percentile(blob.posterior.samples['photometry'] * (1 + z), 25, axis=0)/1e-18
    phot_75 = np.percentile(blob.posterior.samples['photometry'] * (1 + z), 75, axis=0)/1e-18

    flux_med = np.median(blob.posterior.samples['spectrum_full'], axis=0)/1e-18
    flux_25 = np.percentile(blob.posterior.samples['spectrum_full'], 25, axis=0)/1e-18
    flux_75 = np.percentile(blob.posterior.samples['spectrum_full'], 75, axis=0)/1e-18

    ages = blob.posterior.sfh.ages / 1e6

    sfh_med = np.median(blob.posterior.samples['sfh'], axis=0)
    sfh_25 = np.percentile(blob.posterior.samples['sfh'], 25, axis=0)
    sfh_75 = np.percentile(blob.posterior.samples['sfh'], 75, axis=0)

    ax.plot(wl, flux_med, color=spec_color)
    ax.fill_between(wl, flux_25, flux_75, alpha=0.25, color=spec_color)
    ax.set_xlim(2450, 8550)
    ax.set_ylim(0, 2*flux_med[spec_range][700])

    ax.bar(x=photometry[:, 0], height=phot_75 - phot_25, bottom=phot_25, width=150, zorder=10, color=phot_color,
           **{'edgecolor': phot_color, 'alpha': 0.8})
    ax.errorbar(photometry[:, 0], photometry[:, 1], photometry[:, 2], barsabove=True, fmt='ok', elinewidth=1,
                markersize=5,
                ecolor='k', zorder=15)
    ax.set_ylabel(r'$F_\lambda\,\mathrm{[10^{-18}\,erg\,s^{-1}\,cm^{-2}\,\AA^{-1}]}$', fontsize=20)
    ax.set_xlabel(r'$\lambda\,\mathrm{[\AA]}$', fontsize=20)

    if plot_distributions:
        inset = ax.inset_axes([0.78, 0.79, 0.215, 0.20])
        inset_age = ax.inset_axes([0.78, 0.5, 0.1, 0.18])
        inset_mass = ax.inset_axes([0.895, 0.5, 0.1, 0.18])

        sns.despine(ax=inset_age, left=True)
        sns.despine(ax=inset_mass, left=True)
        sns.despine(ax=inset)

        inset_age.set_yticks([])
        inset_mass.set_yticks([])

        inset.plot(ages, sfh_med, color=spec_color)
        inset.fill_between(ages, sfh_25, sfh_75, alpha=0.25, color=spec_color)
        inset.set_xlim(ages[0], 3.25 * np.median(blob.posterior.samples['delayed2:age'] * 1e3))
        inset.set_ylim(0, sfh_75[0])
        inset.set_xlabel(r'$t\,\mathrm{[Myr]}$', fontsize=13)
        inset.set_ylabel(r'SFR\,$\,\mathrm{[M_\odot/yr]}$', fontsize=13)

        sns.kdeplot(1e3 * blob.posterior.samples['delayed2:age'], ax=inset_age, color=spec_color, fill=True, alpha=0.5)
        sns.kdeplot(blob.posterior.samples['delayed2:massformed'], ax=inset_mass, color=spec_color, fill=True, alpha=0.5)

        inset_mass.set_ylabel(r'')
        inset_age.set_xlabel(r'$t_\star\mathrm{[Myr]}$', fontsize=13)
        inset_mass.set_xlabel(r'$\log\,M/M_\odot$', fontsize=13)
