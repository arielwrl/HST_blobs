"""

ariel@oapd
23/03/2023

Some trends with distance

"""

import numpy as np
import matplotlib
# matplotlib.use("Qt5Agg")
import matplotlib.pyplot as plt
import pandas as pd
from astropy.table import Table
import seaborn as sns
from astropy.cosmology import FlatLambdaCDM
from toolbox import plot_tools
import scipy.stats

cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

sns.set_style('ticks')
halpha_palette = sns.light_palette('goldenrod', 4)
f275w_palette = sns.light_palette('mediumvioletred', 4)
f606w_palette = sns.light_palette('indigo', 4)

halpha_input = Table.read('/home/ariel/Workspace/GASP/HST/Data/halpha_bagpipes_input.fits')
f275w_input = Table.read('/home/ariel/Workspace/GASP/HST/Data/f275w_bagpipes_input.fits')
f606w_input = Table.read('/home/ariel/Workspace/GASP/HST/Data/f606w_bagpipes_input.fits')

output_halpha = Table.read('/home/ariel/Workspace/GASP/HST/Data/halpha_dexp_logprior_single_bagpipes_results.fits')
output_f275w = Table.read('/home/ariel/Workspace/GASP/HST/Data/f275w_dexp_logprior_single_bagpipes_results.fits')
output_f606w = Table.read('/home/ariel/Workspace/GASP/HST/Data/f606w_dexp_logprior_single_bagpipes_results.fits')

output_halpha = output_halpha[(~halpha_input['disk']) & ((halpha_input['level'] == 0) | (halpha_input['leaf_flag'] == 1))]
output_f275w = output_f275w[(~f275w_input['disk']) & ((f275w_input['level'] == 0) | (f275w_input['leaf_flag'] == 1))]

halpha_input = halpha_input[(~halpha_input['disk']) & ((halpha_input['level'] == 0) | (halpha_input['leaf_flag'] == 1))]
f275w_input = f275w_input[(~f275w_input['disk']) & ((f275w_input['level'] == 0) | (f275w_input['leaf_flag'] == 1))]

halpha_flag = f606w_input['halpha_match']

halpha_input = halpha_input[~output_halpha['bad_double_fit'] & ~output_halpha['bad_fit']]
output_halpha = output_halpha[~output_halpha['bad_double_fit'] & ~output_halpha['bad_fit']]

f275w_input = f275w_input[~output_f275w['bad_double_fit'] & ~output_f275w['bad_fit']]
output_f275w = output_f275w[~output_f275w['bad_double_fit'] & ~output_f275w['bad_fit']]

f606w_input = f606w_input[~output_f606w['bad_fit']]
output_f606w = output_f606w[~output_f606w['bad_fit']]

# halpha_frac_flag = output_halpha['mass_frac_20'] > 0.1
# f275w_frac_flag = output_f275w['mass_frac_20'] > 0.1
# f606w_frac_flag = output_f606w['mass_frac_20'] > 0.1

halpha_frac_flag = output_halpha['late_decliner']
f275w_frac_flag = output_f275w['late_decliner']
f606w_frac_flag = output_f606w['late_decliner']

output_halpha['Av'] = output_halpha['Av'] * output_halpha['eta']
output_f275w['Av'] = output_f275w['Av'] * output_f275w['eta']
output_f606w['Av'] = output_f606w['Av'] * output_f606w['eta']

output_halpha['mwage'] = np.log10(1e3 * output_halpha['mwage'])
output_f275w['mwage'] = np.log10(1e3 * output_f275w['mwage'])
output_f606w['mwage'] = np.log10(1e3 * output_f606w['mwage'])

# output_halpha['sfr'] = np.log10(output_halpha['sfr'])
# output_f275w['sfr'] = np.log10(output_f275w['sfr'])
# output_f606w['sfr'] = np.log10(output_f606w['sfr'])

# sf_flag_halpha = output_halpha['sfr'] > 0.70 * output_halpha['stellar_mass'] - 7.25
# sf_flag_f275w = output_f275w['sfr'] > 0.70 * output_f275w['stellar_mass'] - 7.25
# sf_flag_f606w = output_f606w['sfr'] > 0.70 * output_f606w['stellar_mass'] - 7.25

radius_dict = {'JO201': 7.12,
               'JO204': 5.16,
               'JO206': 9.39,
               'JW100': 6.9,
               'JW39': 8.33,
               'JO175': 3.61}

redshift_dict = {'JO201': 0.044631,
                 'JO204': 0.042372,
                 'JO206': 0.051089,
                 'JW100': 0.061891,
                 'JW39': 0.066319,
                 'JO175': 0.046750}

kpc_radius_dict = {}
halpha_input['galaxy_radius_kpc'] = np.zeros_like(halpha_input['dist_kpc'])
halpha_input['galaxy_radius_arcsec'] = np.zeros_like(halpha_input['dist_kpc'])
f275w_input['galaxy_radius_kpc'] = np.zeros_like(f275w_input['dist_kpc'])
f275w_input['galaxy_radius_arcsec'] = np.zeros_like(f275w_input['dist_kpc'])
f606w_input['galaxy_radius_kpc'] = np.zeros_like(f606w_input['dist_kpc'])
f606w_input['galaxy_radius_arcsec'] = np.zeros_like(f606w_input['dist_kpc'])

for key in redshift_dict.keys():

    d_A = cosmo.angular_diameter_distance(z=redshift_dict[key])

    theta = radius_dict[key]
    theta_radian = theta * np.pi / 180 / 3600
    kpc_radius_dict[key] = d_A.value * theta_radian * 1000

    halpha_input['galaxy_radius_kpc'][halpha_input['gal'] == key] = np.full((halpha_input['gal'] == key).sum(),
                                                                               kpc_radius_dict[key])
    halpha_input['galaxy_radius_arcsec'][halpha_input['gal'] == key] = np.full((halpha_input['gal'] == key).sum(),
                                                                                 radius_dict[key])
    f275w_input['galaxy_radius_kpc'][f275w_input['gal'] == key] = np.full((f275w_input['gal'] == key).sum(),
                                                                             kpc_radius_dict[key])
    f275w_input['galaxy_radius_arcsec'][f275w_input['gal'] == key] = np.full((f275w_input['gal'] == key).sum(),
                                                                                radius_dict[key])
    f606w_input['galaxy_radius_kpc'][f606w_input['gal'] == key] = np.full((f606w_input['gal'] == key).sum(),
                                                                             kpc_radius_dict[key])
    f606w_input['galaxy_radius_arcsec'][f606w_input['gal'] == key] = np.full((f606w_input['gal'] == key).sum(),
                                                                                radius_dict[key])


fig, ax = plt.subplots(4, 3, figsize=(12.5, 12.5), sharex='col')

tail = halpha_input['tail'] & (np.log10(output_halpha['sfr']) > -20)
extraplanar = halpha_input['extraplanar']

ax[0, 0].scatter(np.log10(halpha_input['dist_kpc']/halpha_input['galaxy_radius_kpc'])[extraplanar],
                 np.log10(output_halpha['mwage'])[extraplanar], color=halpha_palette[3], edgecolors='w',
                 label='Extraplanar', s=25)
ax[0, 0].scatter(np.log10(halpha_input['dist_kpc']/halpha_input['galaxy_radius_kpc'])[tail],
                 np.log10(output_halpha['mwage'])[tail], color=halpha_palette[1], edgecolors='w', label='Tail', s=25)
plot_tools.plot_median_in_bins(np.log10(halpha_input['dist_kpc']/halpha_input['galaxy_radius_kpc']),
                               np.log10(output_halpha['mwage']), nbins=6, percentile_style='errorbar',
                               color='k', percentiles_color='k', ax=ax[0, 0], point_size=50, point_edgewidths=1,
                               point_edgecolors='w')

ax[1, 0].scatter(np.log10(halpha_input['dist_kpc']/halpha_input['galaxy_radius_kpc'])[extraplanar],
                 output_halpha['stellar_mass'][extraplanar], color=halpha_palette[3], edgecolors='w',
                 label='Extraplanar', s=25)
ax[1, 0].scatter(np.log10(halpha_input['dist_kpc']/halpha_input['galaxy_radius_kpc'])[tail],
                 output_halpha['stellar_mass'][tail], color=halpha_palette[1], edgecolors='w', label='Tail', s=25)
plot_tools.plot_median_in_bins(np.log10(halpha_input['dist_kpc']/halpha_input['galaxy_radius_kpc']),
                               output_halpha['stellar_mass'], nbins=6, percentile_style='errorbar',
                               color='k', percentiles_color='k', ax=ax[1, 0], point_size=50, point_edgewidths=1,
                               point_edgecolors='w')

ax[2, 0].scatter(np.log10(halpha_input['dist_kpc']/halpha_input['galaxy_radius_kpc'])[extraplanar],
                 np.log10(output_halpha['sfr'])[extraplanar], color=halpha_palette[3], edgecolors='w',
                 label='Extraplanar', s=25)
ax[2, 0].scatter(np.log10(halpha_input['dist_kpc']/halpha_input['galaxy_radius_kpc'])[tail],
                 np.log10(output_halpha['sfr'])[tail], color=halpha_palette[1], edgecolors='w', label='Tail', s=25)
plot_tools.plot_median_in_bins(np.log10(halpha_input['dist_kpc']/halpha_input['galaxy_radius_kpc']),
                               np.log10(output_halpha['sfr']), nbins=6, percentile_style='errorbar',
                               color='k', percentiles_color='k', ax=ax[2, 0], point_size=50, point_edgewidths=1,
                               point_edgecolors='w')

ax[3, 0].scatter(np.log10(halpha_input['dist_kpc']/halpha_input['galaxy_radius_kpc'])[extraplanar & halpha_frac_flag],
                 output_halpha['Av'][extraplanar & halpha_frac_flag], color=halpha_palette[3], edgecolors='w',
                 label='Extraplanar', s=25)
ax[3, 0].scatter(np.log10(halpha_input['dist_kpc']/halpha_input['galaxy_radius_kpc'])[tail & halpha_frac_flag],
                 output_halpha['Av'][tail & halpha_frac_flag], color=halpha_palette[1], edgecolors='w', label='Tail', s=25)
plot_tools.plot_median_in_bins(np.log10(halpha_input['dist_kpc']/halpha_input['galaxy_radius_kpc'])[halpha_frac_flag],
                               output_halpha['Av'][halpha_frac_flag], nbins=6, percentile_style='errorbar',
                               color='k', percentiles_color='k', ax=ax[3, 0], point_size=50, point_edgewidths=1,
                               point_edgecolors='w')


tail = f275w_input['tail']
extraplanar = f275w_input['extraplanar']

ax[0, 1].scatter(np.log10(f275w_input['dist_kpc']/f275w_input['galaxy_radius_kpc'])[extraplanar],
                 np.log10(output_f275w['mwage'])[extraplanar], color=f275w_palette[3], edgecolors='w',
                 label='Extraplanar', s=25)
ax[0, 1].scatter(np.log10(f275w_input['dist_kpc']/f275w_input['galaxy_radius_kpc'])[tail],
                 np.log10(output_f275w['mwage'])[tail], color=f275w_palette[1], edgecolors='w', label='Tail', s=25)
plot_tools.plot_median_in_bins(np.log10(f275w_input['dist_kpc']/f275w_input['galaxy_radius_kpc']),
                               np.log10(output_f275w['mwage']), nbins=6, percentile_style='errorbar', color='k',
                               percentiles_color='k', ax=ax[0, 1], point_size=50, point_edgewidths=1,
                               point_edgecolors='w')

ax[1, 1].scatter(np.log10(f275w_input['dist_kpc']/f275w_input['galaxy_radius_kpc'])[extraplanar],
                 output_f275w['stellar_mass'][extraplanar], color=f275w_palette[3], edgecolors='w',
                 label='Extraplanar', s=25)
ax[1, 1].scatter(np.log10(f275w_input['dist_kpc']/f275w_input['galaxy_radius_kpc'])[tail],
                 output_f275w['stellar_mass'][tail], color=f275w_palette[1], edgecolors='w', label='Tail', s=25)
plot_tools.plot_median_in_bins(np.log10(f275w_input['dist_kpc']/f275w_input['galaxy_radius_kpc']),
                               output_f275w['stellar_mass'], nbins=6, percentile_style='errorbar', color='k',
                               percentiles_color='k', ax=ax[1, 1], point_size=50, point_edgewidths=1,
                               point_edgecolors='w')

sf_flag = np.log10(output_f275w['sfr']) > 0.70 * output_f275w['stellar_mass'] - 7.25
ax[2, 1].scatter(np.log10(f275w_input['dist_kpc']/f275w_input['galaxy_radius_kpc'])[extraplanar & sf_flag],
                 np.log10(output_f275w['sfr'][extraplanar & sf_flag]), color=f275w_palette[3], edgecolors='w',
                 label='Extraplanar', s=25)
ax[2, 1].scatter(np.log10(f275w_input['dist_kpc']/f275w_input['galaxy_radius_kpc'])[tail & sf_flag],
                 np.log10(output_f275w['sfr'][tail & sf_flag]), color=f275w_palette[1], edgecolors='w', label='Tail', s=25)
plot_tools.plot_median_in_bins(np.log10(f275w_input['dist_kpc']/f275w_input['galaxy_radius_kpc'])[sf_flag],
                               np.log10(output_f275w['sfr'])[sf_flag], nbins=6, percentile_style='errorbar', color='k',
                               percentiles_color='k', ax=ax[2, 1], point_size=50, point_edgewidths=1,
                               point_edgecolors='w')

inset_ax = ax[2, 1].inset_axes([0, 0, 1, 0.3])
inset_ax.set_xticklabels([])
inset_ax.set_yticklabels([])
inset_ax.spines["top"].set_visible(False)
inset_ax.spines["right"].set_visible(False)
inset_ax.patch.set_alpha(0)
inset_ax.tick_params(axis='both', length=0)

inset_ax.hist(np.log10(f275w_input['dist_kpc']/f275w_input['galaxy_radius_kpc'])[~sf_flag],
              bins=14, range=[-0.2, 1.2], histtype='step', color='k', ls='dashed', density=True)
inset_ax.hist(np.log10(f275w_input['dist_kpc']/f275w_input['galaxy_radius_kpc'])[sf_flag],
              bins=14, range=[-0.2, 1.2], histtype='step', color='k', density=True)

ax[3, 1].scatter(np.log10(f275w_input['dist_kpc']/f275w_input['galaxy_radius_kpc'])[extraplanar & f275w_frac_flag],
                 output_f275w['Av'][extraplanar & f275w_frac_flag], color=f275w_palette[3], edgecolors='w',
                 label='Extraplanar', s=25)
ax[3, 1].scatter(np.log10(f275w_input['dist_kpc']/f275w_input['galaxy_radius_kpc'])[tail & f275w_frac_flag],
                 output_f275w['Av'][tail & f275w_frac_flag], color=f275w_palette[1], edgecolors='w', label='Tail', s=25)
plot_tools.plot_median_in_bins(np.log10(f275w_input['dist_kpc']/f275w_input['galaxy_radius_kpc'])[f275w_frac_flag],
                               output_f275w['Av'][f275w_frac_flag], nbins=6, percentile_style='errorbar', color='k',
                               percentiles_color='k', ax=ax[3, 1], point_size=50, point_edgewidths=1,
                               point_edgecolors='w')

ax[0, 2].scatter(np.log10(f606w_input['dist_kpc']/f606w_input['galaxy_radius_kpc']),
                 np.log10(output_f606w['mwage']), color=f606w_palette[2], edgecolors='w', s=25,
                 label='Tail')
# ax[0, 2].scatter(np.log10(f606w_input['dist_kpc']/f606w_input['galaxy_radius_kpc'])[halpha_flag],
#                  np.log10(output_f606w['mwage'][halpha_flag]), color=f606w_palette[2], edgecolors='goldenrod',
#                  linewidths=2, s=25, label=r'With $H\alpha$ Emission')
plot_tools.plot_median_in_bins(np.log10(f606w_input['dist_kpc']/f606w_input['galaxy_radius_kpc']),
                               np.log10(output_f606w['mwage']), nbins=5, percentile_style='errorbar', color='k',
                               percentiles_color='k', ax=ax[0, 2], point_size=50, point_edgewidths=1,
                               point_edgecolors='w')

ax[1, 2].scatter(np.log10(f606w_input['dist_kpc']/f606w_input['galaxy_radius_kpc']),
                 output_f606w['stellar_mass'], color=f606w_palette[2], edgecolors='w', s=25,
                 label='Tail')
# ax[1, 2].scatter(np.log10(f606w_input['dist_kpc']/f606w_input['galaxy_radius_kpc'])[halpha_flag],
#                  output_f606w['stellar_mass'][halpha_flag], color=f606w_palette[2], edgecolors='goldenrod',
#                  linewidths=2, s=25, label=r'With $H\alpha$ Emission')
plot_tools.plot_median_in_bins(np.log10(f606w_input['dist_kpc']/f606w_input['galaxy_radius_kpc']),
                               output_f606w['stellar_mass'], nbins=5, percentile_style='errorbar', color='k',
                               percentiles_color='k', ax=ax[1, 2], point_size=50, point_edgewidths=1,
                               point_edgecolors='w')

sf_flag = np.log10(output_f606w['sfr']) > 0.70 * output_f606w['stellar_mass'] - 7.25
ax[2, 2].scatter(np.log10(f606w_input['dist_kpc']/f606w_input['galaxy_radius_kpc'])[sf_flag],
                 np.log10(output_f606w['sfr'])[sf_flag], color=f606w_palette[2], edgecolors='w', s=25,
                 label='Tail')
plot_tools.plot_median_in_bins(np.log10(f606w_input['dist_kpc']/f606w_input['galaxy_radius_kpc'])[sf_flag],
                               np.log10(output_f606w['sfr'])[sf_flag], nbins=5, percentile_style='errorbar', color='k',
                               percentiles_color='k', ax=ax[2, 2], point_size=50, point_edgewidths=1,
                               point_edgecolors='w')

inset_ax = ax[2, 2].inset_axes([0, 0, 1, 0.3])
inset_ax.set_xticklabels([])
inset_ax.set_yticklabels([])
inset_ax.spines["top"].set_visible(False)
inset_ax.spines["right"].set_visible(False)
inset_ax.patch.set_alpha(0)
inset_ax.tick_params(axis='both', length=0)

inset_ax.hist(np.log10(f606w_input['dist_kpc']/f606w_input['galaxy_radius_kpc'])[~sf_flag],
              bins=12, range=[0, 1.2], histtype='step', color='k', ls='dashed', density=True)
inset_ax.hist(np.log10(f606w_input['dist_kpc']/f606w_input['galaxy_radius_kpc'])[sf_flag],
              bins=12, range=[0, 1.2], histtype='step', color='k', density=True)


ax[3, 2].scatter(np.log10(f606w_input['dist_kpc']/f606w_input['galaxy_radius_kpc'])[f606w_frac_flag],
                 output_f606w['Av'][f606w_frac_flag], color=f606w_palette[2], edgecolors='w', s=25,
                 label='Tail')
plot_tools.plot_median_in_bins(np.log10(f606w_input['dist_kpc']/f606w_input['galaxy_radius_kpc'])[f606w_frac_flag],
                               output_f606w['Av'][f606w_frac_flag], nbins=5, percentile_style='errorbar', color='k',
                               percentiles_color='k', ax=ax[3, 2], point_size=50, point_edgewidths=1,
                               point_edgecolors='w')

sns.despine()

ax[2, 1].set_ylim(-5.6, -0.5)
ax[2, 2].set_ylim(-5.1, -0.5)

ax[0, 0].legend(frameon=True, fontsize=14, loc=3)
ax[0, 1].legend(frameon=True, fontsize=14, loc=3)
ax[0, 2].legend(frameon=True, fontsize=14, loc=3)

# ax[0, 0].annotate(r'$H\alpha$', fontsize=20, xy=(0.02, 0.02), xycoords='axes fraction')
# ax[0, 1].annotate(r'F275W', fontsize=20, xy=(0.02, 0.02), xycoords='axes fraction')
# ax[0, 2].annotate(r'Star-Forming Complexes', fontsize=20, xy=(0.02, 0.02), xycoords='axes fraction')

ax[0, 0].set_title(r'$H\alpha$', fontsize=20)
ax[0, 1].set_title(r'F275W', fontsize=20)
ax[0, 2].set_title(r'Star-Forming Complexes', fontsize=20)

ax[3, 1].set_xlabel(r'$\log\,d/R_e$', fontsize=20)
ax[0, 0].set_ylabel(r'$\log\;\langle\,t_\star\,\rangle_M\,\mathrm{[Myr]}$', fontsize=20)
ax[1, 0].set_ylabel(r'$\log\,M_\star/M_\odot$', fontsize=20)
ax[2, 0].set_ylabel(r'$\log \, \mathrm{SFR} \,\mathrm{[M_\odot/yr]}$', fontsize=20)
ax[3, 0].set_ylabel(r'$A_V\,\mathrm{[mag]}$', fontsize=20)

for subplot in ax.ravel():
    subplot.tick_params(axis='both', labelsize=15)

fig.tight_layout()
fig.subplots_adjust(left=0.075, bottom=0.055)

print('mwage', np.corrcoef(np.log10(halpha_input['dist_kpc']/halpha_input['galaxy_radius_kpc']),
                               np.log10(output_halpha['mwage'])))
print('mass', np.corrcoef(np.log10(halpha_input['dist_kpc']/halpha_input['galaxy_radius_kpc']),
                               output_halpha['stellar_mass']))
print('mwage', np.corrcoef(np.log10(f275w_input['dist_kpc']/f275w_input['galaxy_radius_kpc']),
                               np.log10(output_f275w['mwage'])))
print('mass', np.corrcoef(np.log10(f275w_input['dist_kpc']/f275w_input['galaxy_radius_kpc']),
                               output_f275w['stellar_mass']))
print('mwage', np.corrcoef(np.log10(f606w_input['dist_kpc']/f606w_input['galaxy_radius_kpc']),
                               np.log10(output_f606w['mwage'])))
print('mass', np.corrcoef(np.log10(f606w_input['dist_kpc']/f606w_input['galaxy_radius_kpc']),
                               output_f606w['stellar_mass']))

plt.savefig('/home/ariel/Workspace/GASP/HST/Analysis/Plots/Paper/distance_trends.png', dpi=300)
plt.close()

fig, ax = plt.subplots(1, 3, figsize=(12.5, 5))

tail = halpha_input['tail'] & (np.log10(output_halpha['sfr']) > -20)
extraplanar = halpha_input['extraplanar']

ax[0].scatter(np.log10(halpha_input['dist_kpc']/halpha_input['galaxy_radius_kpc'])[extraplanar],
                 np.log10(output_halpha['mwage'])[extraplanar], color=halpha_palette[3], edgecolors='w',
                 label='Extraplanar', s=25)
ax[0].scatter(np.log10(halpha_input['dist_kpc']/halpha_input['galaxy_radius_kpc'])[tail],
                 np.log10(output_halpha['mwage'])[tail], color=halpha_palette[1], edgecolors='w', label='Tail', s=25)
plot_tools.plot_median_in_bins(np.log10(halpha_input['dist_kpc']/halpha_input['galaxy_radius_kpc']),
                               np.log10(output_halpha['mwage']), nbins=6, percentile_style='errorbar',
                               color='k', percentiles_color='k', ax=ax[0], point_size=50, point_edgewidths=1,
                               point_edgecolors='w')

tail = f275w_input['tail']
extraplanar = f275w_input['extraplanar']

ax[1].scatter(np.log10(f275w_input['dist_kpc']/f275w_input['galaxy_radius_kpc'])[extraplanar],
                 np.log10(output_f275w['mwage'])[extraplanar], color=f275w_palette[3], edgecolors='w',
                 label='Extraplanar', s=25)
ax[1].scatter(np.log10(f275w_input['dist_kpc']/f275w_input['galaxy_radius_kpc'])[tail],
                 np.log10(output_f275w['mwage'])[tail], color=f275w_palette[1], edgecolors='w', label='Tail', s=25)
plot_tools.plot_median_in_bins(np.log10(f275w_input['dist_kpc']/f275w_input['galaxy_radius_kpc']),
                               np.log10(output_f275w['mwage']), nbins=6, percentile_style='errorbar', color='k',
                               percentiles_color='k', ax=ax[1], point_size=50, point_edgewidths=1,
                               point_edgecolors='w')

ax[2].scatter(np.log10(f606w_input['dist_kpc']/f606w_input['galaxy_radius_kpc']),
                 np.log10(output_f606w['mwage']), color=f606w_palette[2], edgecolors='w', s=25,
                 label='Tail')
plot_tools.plot_median_in_bins(np.log10(f606w_input['dist_kpc']/f606w_input['galaxy_radius_kpc']),
                               np.log10(output_f606w['mwage']), nbins=5, percentile_style='errorbar', color='k',
                               percentiles_color='k', ax=ax[2], point_size=50, point_edgewidths=1,
                               point_edgecolors='w')

sns.despine()

ax[1].set_xlabel(r'$\log\,d/R_e$', fontsize=20)

ax[0].set_title(r'$H\alpha$', fontsize=20)
ax[1].set_title(r'F275W', fontsize=20)
ax[2].set_title(r'Star-Forming Complexes', fontsize=20)

ax[0].set_ylabel(r'$\log\;\langle\,t_\star\,\rangle_M\,\mathrm{[Myr]}$', fontsize=20)

ax[0].legend(frameon=True, fontsize=14, loc=3)
ax[1].legend(frameon=True, fontsize=14, loc=3)
ax[2].legend(frameon=True, fontsize=14, loc=3)

for subplot in ax:
    subplot.tick_params(axis='both', labelsize=15)

fig.tight_layout()

ax[0].legend(frameon=True, fontsize=14, loc=3)
ax[1].legend(frameon=True, fontsize=14, loc=3)
ax[2].legend(frameon=True, fontsize=14, loc=3)

for subplot in ax:
    subplot.tick_params(axis='both', labelsize=15)

fig.tight_layout()

plt.savefig('/home/ariel/Workspace/GASP/HST/Analysis/Plots/Presentations/mwage_distance_trends.png', dpi=300)


fig, ax = plt.subplots(1, 3, figsize=(12.5, 5))

tail = halpha_input['tail'] & (np.log10(output_halpha['sfr']) > -20)
extraplanar = halpha_input['extraplanar']

ax[0].scatter(np.log10(halpha_input['dist_kpc']/halpha_input['galaxy_radius_kpc'])[extraplanar],
                 output_halpha['stellar_mass'][extraplanar], color=halpha_palette[3], edgecolors='w',
                 label='Extraplanar', s=25)
ax[0].scatter(np.log10(halpha_input['dist_kpc']/halpha_input['galaxy_radius_kpc'])[tail],
                 output_halpha['stellar_mass'][tail], color=halpha_palette[1], edgecolors='w', label='Tail', s=25)
plot_tools.plot_median_in_bins(np.log10(halpha_input['dist_kpc']/halpha_input['galaxy_radius_kpc']),
                               output_halpha['stellar_mass'], nbins=6, percentile_style='errorbar',
                               color='k', percentiles_color='k', ax=ax[0], point_size=50, point_edgewidths=1,
                               point_edgecolors='w')

tail = f275w_input['tail']
extraplanar = f275w_input['extraplanar']

ax[1].scatter(np.log10(f275w_input['dist_kpc']/f275w_input['galaxy_radius_kpc'])[extraplanar],
                 output_f275w['stellar_mass'][extraplanar], color=f275w_palette[3], edgecolors='w',
                 label='Extraplanar', s=25)
ax[1].scatter(np.log10(f275w_input['dist_kpc']/f275w_input['galaxy_radius_kpc'])[tail],
                 output_f275w['stellar_mass'][tail], color=f275w_palette[1], edgecolors='w', label='Tail', s=25)
plot_tools.plot_median_in_bins(np.log10(f275w_input['dist_kpc']/f275w_input['galaxy_radius_kpc']),
                               output_f275w['stellar_mass'], nbins=6, percentile_style='errorbar', color='k',
                               percentiles_color='k', ax=ax[1], point_size=50, point_edgewidths=1,
                               point_edgecolors='w')

ax[2].scatter(np.log10(f606w_input['dist_kpc']/f606w_input['galaxy_radius_kpc']),
                 output_f606w['stellar_mass'], color=f606w_palette[2], edgecolors='w', s=25,
                 label='Tail')
plot_tools.plot_median_in_bins(np.log10(f606w_input['dist_kpc']/f606w_input['galaxy_radius_kpc']),
                               output_f606w['stellar_mass'], nbins=5, percentile_style='errorbar', color='k',
                               percentiles_color='k', ax=ax[2], point_size=50, point_edgewidths=1,
                               point_edgecolors='w')

sns.despine()

ax[1].set_xlabel(r'$\log\,d/R_e$', fontsize=20)

ax[0].set_title(r'$H\alpha$', fontsize=20)
ax[1].set_title(r'F275W', fontsize=20)
ax[2].set_title(r'Star-Forming Complexes', fontsize=20)

ax[0].set_ylabel(r'$\log\,M_\star/M_\odot$', fontsize=20)

ax[0].legend(frameon=True, fontsize=14, loc=1)
ax[1].legend(frameon=True, fontsize=14, loc=1)
ax[2].legend(frameon=True, fontsize=14, loc=1)

for subplot in ax:
    subplot.tick_params(axis='both', labelsize=15)

fig.tight_layout()

plt.savefig('/home/ariel/Workspace/GASP/HST/Analysis/Plots/Presentations/mass_distance_trends.png', dpi=300)


fig, ax = plt.subplots(1, 3, figsize=(12.5, 5))

tail = halpha_input['tail'] & (np.log10(output_halpha['sfr']) > -20)
extraplanar = halpha_input['extraplanar']

ax[0].scatter(np.log10(halpha_input['dist_kpc']/halpha_input['galaxy_radius_kpc'])[extraplanar],
                 np.log10(output_halpha['sfr'])[extraplanar], color=halpha_palette[3], edgecolors='w',
                 label='Extraplanar', s=25)
ax[0].scatter(np.log10(halpha_input['dist_kpc']/halpha_input['galaxy_radius_kpc'])[tail],
                 np.log10(output_halpha['sfr'])[tail], color=halpha_palette[1], edgecolors='w', label='Tail', s=25)
plot_tools.plot_median_in_bins(np.log10(halpha_input['dist_kpc']/halpha_input['galaxy_radius_kpc']),
                               np.log10(output_halpha['sfr']), nbins=6, percentile_style='errorbar',
                               color='k', percentiles_color='k', ax=ax[0], point_size=50, point_edgewidths=1,
                               point_edgecolors='w')

tail = f275w_input['tail']
extraplanar = f275w_input['extraplanar']

sf_flag = np.log10(output_f275w['sfr']) > 0.70 * output_f275w['stellar_mass'] - 7.25
ax[1].scatter(np.log10(f275w_input['dist_kpc']/f275w_input['galaxy_radius_kpc'])[extraplanar & sf_flag],
                 np.log10(output_f275w['sfr'][extraplanar & sf_flag]), color=f275w_palette[3], edgecolors='w',
                 label='Extraplanar', s=25)
ax[1].scatter(np.log10(f275w_input['dist_kpc']/f275w_input['galaxy_radius_kpc'])[tail & sf_flag],
                 np.log10(output_f275w['sfr'][tail & sf_flag]), color=f275w_palette[1], edgecolors='w', label='Tail', s=25)
plot_tools.plot_median_in_bins(np.log10(f275w_input['dist_kpc']/f275w_input['galaxy_radius_kpc'])[sf_flag],
                               np.log10(output_f275w['sfr'])[sf_flag], nbins=6, percentile_style='errorbar', color='k',
                               percentiles_color='k', ax=ax[1], point_size=50, point_edgewidths=1,
                               point_edgecolors='w')

inset_ax = ax[1].inset_axes([0, 0, 1, 0.3])
inset_ax.set_xticklabels([])
inset_ax.set_yticklabels([])
inset_ax.spines["top"].set_visible(False)
inset_ax.spines["right"].set_visible(False)
inset_ax.patch.set_alpha(0)
inset_ax.tick_params(axis='both', length=0)

inset_ax.hist(np.log10(f275w_input['dist_kpc']/f275w_input['galaxy_radius_kpc'])[~sf_flag],
              bins=14, range=[-0.2, 1.2], histtype='step', color='k', ls='dashed', density=True)
inset_ax.hist(np.log10(f275w_input['dist_kpc']/f275w_input['galaxy_radius_kpc'])[sf_flag],
              bins=14, range=[-0.2, 1.2], histtype='step', color='k', density=True)

sf_flag = np.log10(output_f606w['sfr']) > 0.70 * output_f606w['stellar_mass'] - 7.25
ax[2].scatter(np.log10(f606w_input['dist_kpc']/f606w_input['galaxy_radius_kpc'])[sf_flag],
                 np.log10(output_f606w['sfr'])[sf_flag], color=f606w_palette[2], edgecolors='w', s=25,
                 label='Tail')
plot_tools.plot_median_in_bins(np.log10(f606w_input['dist_kpc']/f606w_input['galaxy_radius_kpc'])[sf_flag],
                               np.log10(output_f606w['sfr'])[sf_flag], nbins=5, percentile_style='errorbar', color='k',
                               percentiles_color='k', ax=ax[2], point_size=50, point_edgewidths=1,
                               point_edgecolors='w')

inset_ax = ax[2].inset_axes([0, 0, 1, 0.3])
inset_ax.set_xticklabels([])
inset_ax.set_yticklabels([])
inset_ax.spines["top"].set_visible(False)
inset_ax.spines["right"].set_visible(False)
inset_ax.patch.set_alpha(0)
inset_ax.tick_params(axis='both', length=0)

inset_ax.hist(np.log10(f606w_input['dist_kpc']/f606w_input['galaxy_radius_kpc'])[~sf_flag],
              bins=12, range=[0, 1.2], histtype='step', color='k', ls='dashed', density=True)
inset_ax.hist(np.log10(f606w_input['dist_kpc']/f606w_input['galaxy_radius_kpc'])[sf_flag],
              bins=12, range=[0, 1.2], histtype='step', color='k', density=True)

sns.despine()

ax[1].set_ylim(-5.6, -0.5)
ax[2].set_ylim(-5.1, -0.5)

ax[1].set_xlabel(r'$\log\,d/R_e$', fontsize=20)

ax[0].set_title(r'$H\alpha$', fontsize=20)
ax[1].set_title(r'F275W', fontsize=20)
ax[2].set_title(r'Star-Forming Complexes', fontsize=20)

ax[0].set_ylabel(r'$\log \, \mathrm{SFR} \,\mathrm{[M_\odot/yr]}$', fontsize=20)

ax[0].legend(frameon=True, fontsize=14, loc=1)
ax[1].legend(frameon=True, fontsize=14, loc=1)
ax[2].legend(frameon=True, fontsize=14, loc=1)

for subplot in ax:
    subplot.tick_params(axis='both', labelsize=15)

fig.tight_layout()

plt.savefig('/home/ariel/Workspace/GASP/HST/Analysis/Plots/Presentations/sfr_distance_trends.png', dpi=300)

fig, ax = plt.subplots(1, 3, figsize=(12.5, 5))

tail = halpha_input['tail'] & (np.log10(output_halpha['sfr']) > -20)
extraplanar = halpha_input['extraplanar']

ax[0].scatter(np.log10(halpha_input['dist_kpc']/halpha_input['galaxy_radius_kpc'])[extraplanar & halpha_frac_flag],
                 output_halpha['Av'][extraplanar & halpha_frac_flag], color=halpha_palette[3], edgecolors='w',
                 label='Extraplanar', s=25)
ax[0].scatter(np.log10(halpha_input['dist_kpc']/halpha_input['galaxy_radius_kpc'])[tail & halpha_frac_flag],
                 output_halpha['Av'][tail & halpha_frac_flag], color=halpha_palette[1], edgecolors='w', label='Tail', s=25)
plot_tools.plot_median_in_bins(np.log10(halpha_input['dist_kpc']/halpha_input['galaxy_radius_kpc'])[halpha_frac_flag],
                               output_halpha['Av'][halpha_frac_flag], nbins=6, percentile_style='errorbar',
                               color='k', percentiles_color='k', ax=ax[0], point_size=50, point_edgewidths=1,
                               point_edgecolors='w')

tail = f275w_input['tail']
extraplanar = f275w_input['extraplanar']

ax[1].scatter(np.log10(f275w_input['dist_kpc']/f275w_input['galaxy_radius_kpc'])[extraplanar & f275w_frac_flag],
                 output_f275w['Av'][extraplanar & f275w_frac_flag], color=f275w_palette[3], edgecolors='w',
                 label='Extraplanar', s=25)
ax[1].scatter(np.log10(f275w_input['dist_kpc']/f275w_input['galaxy_radius_kpc'])[tail & f275w_frac_flag],
                 output_f275w['Av'][tail & f275w_frac_flag], color=f275w_palette[1], edgecolors='w', label='Tail', s=25)
plot_tools.plot_median_in_bins(np.log10(f275w_input['dist_kpc']/f275w_input['galaxy_radius_kpc'])[f275w_frac_flag],
                               output_f275w['Av'][f275w_frac_flag], nbins=6, percentile_style='errorbar', color='k',
                               percentiles_color='k', ax=ax[1], point_size=50, point_edgewidths=1,
                               point_edgecolors='w')

ax[2].scatter(np.log10(f606w_input['dist_kpc']/f606w_input['galaxy_radius_kpc'])[f606w_frac_flag],
                 output_f606w['Av'][f606w_frac_flag], color=f606w_palette[2], edgecolors='w', s=25,
                 label='Tail')
plot_tools.plot_median_in_bins(np.log10(f606w_input['dist_kpc']/f606w_input['galaxy_radius_kpc'])[f606w_frac_flag],
                               output_f606w['Av'][f606w_frac_flag], nbins=5, percentile_style='errorbar', color='k',
                               percentiles_color='k', ax=ax[2], point_size=50, point_edgewidths=1,
                               point_edgecolors='w')
sns.despine()

ax[1].set_xlabel(r'$\log\,d/R_e$', fontsize=20)

ax[0].set_title(r'$H\alpha$', fontsize=20)
ax[1].set_title(r'F275W', fontsize=20)
ax[2].set_title(r'Star-Forming Complexes', fontsize=20)

ax[0].set_ylabel(r'$A_V\,\mathrm{[mag]}$', fontsize=20)

ax[0].legend(frameon=True, fontsize=14, loc=1)
ax[1].legend(frameon=True, fontsize=14, loc=1)
ax[2].legend(frameon=True, fontsize=14, loc=1)

for subplot in ax:
    subplot.tick_params(axis='both', labelsize=15)

fig.tight_layout()

plt.savefig('/home/ariel/Workspace/GASP/HST/Analysis/Plots/Presentations/Av_distance_trends.png', dpi=300)


