"""

ariel@oapd
16/08/2023

Galaxy-to-galaxy variations in all kinds of things

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

mass_dict = {'All Galaxies': 0.00001,
             'JO201': 44194800000,
             'JO204': 54968402000,
             'JW100': 292875993000,
             'JW39': 164373004000,
             'JO175': 33957900300,
             'JO206': 77743301000
}

output_halpha = output_halpha[(~halpha_input['disk']) & ((halpha_input['level'] == 0) | (halpha_input['leaf_flag'] == 1))]
output_f275w = output_f275w[(~f275w_input['disk']) & ((f275w_input['level'] == 0) | (f275w_input['leaf_flag'] == 1))]

halpha_input = halpha_input[(~halpha_input['disk']) & ((halpha_input['level'] == 0) | (halpha_input['leaf_flag'] == 1))]
f275w_input = f275w_input[(~f275w_input['disk']) & ((f275w_input['level'] == 0) | (f275w_input['leaf_flag'] == 1))]

f606w_input = f606w_input[output_f606w['stellar_mass'] > 3]
output_f606w = output_f606w[output_f606w['stellar_mass'] > 3]

halpha_flag = f606w_input['halpha_match']

halpha_input = halpha_input[~output_halpha['bad_double_fit'] & ~output_halpha['bad_fit']]
output_halpha = output_halpha[~output_halpha['bad_double_fit'] & ~output_halpha['bad_fit']]

f275w_input = f275w_input[~output_f275w['bad_double_fit'] & ~output_f275w['bad_fit']]
output_f275w = output_f275w[~output_f275w['bad_double_fit'] & ~output_f275w['bad_fit']]

f606w_input = f606w_input[~output_f606w['bad_fit']]
output_f606w = output_f606w[~output_f606w['bad_fit']]

output_halpha['Av'] = output_halpha['Av'] * output_halpha['eta']
output_f275w['Av'] = output_f275w['Av'] * output_f275w['eta']
output_f606w['Av'] = output_f606w['Av'] * output_f606w['eta']

output_halpha['mwage'] = np.log10(1e3 * output_halpha['mwage'])
output_f275w['mwage'] = np.log10(1e3 * output_f275w['mwage'])
output_f606w['mwage'] = np.log10(1e3 * output_f606w['mwage'])

output_halpha['sfr'] = np.log10(output_halpha['sfr'])
output_f275w['sfr'] = np.log10(output_f275w['sfr'])
output_f606w['sfr'] = np.log10(output_f606w['sfr'])

# halpha_frac_flag = np.log10(output_halpha['age']/output_halpha['tau']) < 0.5
# f275w_frac_flag = np.log10(output_f275w['age']/output_f275w['tau']) < 0.5
# f606w_frac_flag = np.log10(output_f606w['age']/output_f606w['tau']) < 0.5

# halpha_frac_flag = -2.5*np.log10(halpha_input['F680N']/halpha_input['F606W']) < 0
# f275w_frac_flag = -2.5*np.log10(f275w_input['F680N']/f275w_input['F606W']) < 0
# f606w_frac_flag = -2.5*np.log10(f606w_input['F680N']/f606w_input['F606W']) < 0

halpha_frac_flag = output_halpha['late_decliner']
f275w_frac_flag = output_f275w['late_decliner']
f606w_frac_flag = output_f606w['late_decliner']

sf_flag_halpha = output_halpha['sfr'] > 0.70 * output_halpha['stellar_mass'] - 7.25
sf_flag_f275w = output_f275w['sfr'] > 0.70 * output_f275w['stellar_mass'] - 7.25
sf_flag_f606w = (output_f606w['sfr'] > 0.70 * output_f606w['stellar_mass'] - 7.25) & (output_f606w['clump_id'] != 'JO201_9518_f606w')

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

distance_halpha = np.log10(halpha_input['dist_kpc']/halpha_input['galaxy_radius_kpc'])
distance_f275w = np.log10(f275w_input['dist_kpc']/f275w_input['galaxy_radius_kpc'])
distance_f606w = np.log10(f606w_input['dist_kpc']/f606w_input['galaxy_radius_kpc'])

variables = ['mwage', 'stellar_mass', 'sfr', 'Av']
labels = [r'$\langle\,t_\star\,\rangle_M\,\mathrm{[Myr]}$', r'$\log\,M_\star/M_\odot$',
          r'$\log \, \mathrm{SFR} \,\mathrm{[M_\odot/yr]}$', r'$A_V\,\mathrm{[mag]}$']

fig, ax = plt.subplots(12, 6, figsize=(9.25, 12), sharex=True, sharey='row')

galaxy_list = ['JO175', 'JO201', 'JO204', 'JO206', 'JW39', 'JW100']

for galaxy_index in range(6):

    galaxy = galaxy_list[galaxy_index]

    for variable_index in range(4):

        variable = variables[variable_index]

        if variable == 'Av':
            ax[variable_index, galaxy_index].scatter(x=distance_halpha[(halpha_input['gal'] == galaxy) & halpha_frac_flag],
                                                     y=output_halpha[variable][(output_halpha['galaxy'] == galaxy)
                                                                               & halpha_frac_flag],
                                                     color=halpha_palette[2], edgecolors='w', s=20)
            ax[variable_index+4, galaxy_index].scatter(x=distance_f275w[(f275w_input['gal'] == galaxy) & f275w_frac_flag],
                                                       y=output_f275w[variable][(output_f275w['galaxy'] == galaxy) &
                                                                                f275w_frac_flag],
                                                       color=f275w_palette[2], edgecolors='w', s=20)
            ax[variable_index+8, galaxy_index].scatter(x=distance_f606w[(f606w_input['gal'] == galaxy) & f606w_frac_flag],
                                                       y=output_f606w[variable][(output_f606w['galaxy'] == galaxy) &
                                                                                f606w_frac_flag],
                                                       color=f606w_palette[2], edgecolors='w', s=20)

            correlation_halpha = scipy.stats.linregress(x=distance_halpha[(halpha_input['gal'] == galaxy) &
                                                                          halpha_frac_flag],
                                                        y=output_halpha[variable][(output_halpha['galaxy'] == galaxy) &
                                                                                  halpha_frac_flag])
            correlation_f275w = scipy.stats.linregress(x=distance_f275w[(f275w_input['gal'] == galaxy) & f275w_frac_flag],
                                                        y=output_f275w[variable][(output_f275w['galaxy'] == galaxy)
                                                                                 & f275w_frac_flag])
            correlation_f606w = scipy.stats.linregress(x=distance_f606w[(f606w_input['gal'] == galaxy) & f606w_frac_flag],
                                                        y=output_f606w[variable][(output_f606w['galaxy'] == galaxy)
                                                                                 & f606w_frac_flag])

        elif variable == 'sfr':
            ax[variable_index, galaxy_index].scatter(x=distance_halpha[(halpha_input['gal'] == galaxy) & sf_flag_halpha],
                                                     y=output_halpha[variable][(output_halpha['galaxy'] == galaxy)
                                                                               & sf_flag_halpha],
                                                     color=halpha_palette[2], edgecolors='w', s=20)
            ax[variable_index+4, galaxy_index].scatter(x=distance_f275w[(f275w_input['gal'] == galaxy) & sf_flag_f275w],
                                                       y=output_f275w[variable][(output_f275w['galaxy'] == galaxy) &
                                                                                sf_flag_f275w],
                                                       color=f275w_palette[2], edgecolors='w', s=20)
            ax[variable_index+8, galaxy_index].scatter(x=distance_f606w[(f606w_input['gal'] == galaxy) & sf_flag_f606w],
                                                       y=output_f606w[variable][(output_f606w['galaxy'] == galaxy) &
                                                                                sf_flag_f606w],
                                                       color=f606w_palette[2], edgecolors='w', s=20)

            correlation_halpha = scipy.stats.linregress(x=distance_halpha[(halpha_input['gal'] == galaxy) &
                                                                          sf_flag_halpha],
                                                        y=output_halpha[variable][(output_halpha['galaxy'] == galaxy) &
                                                                                  sf_flag_halpha])
            correlation_f275w = scipy.stats.linregress(x=distance_f275w[(f275w_input['gal'] == galaxy) & sf_flag_f275w],
                                                        y=output_f275w[variable][(output_f275w['galaxy'] == galaxy)
                                                                                 & sf_flag_f275w])
            correlation_f606w = scipy.stats.linregress(x=distance_f606w[(f606w_input['gal'] == galaxy) & sf_flag_f606w],
                                                        y=output_f606w[variable][(output_f606w['galaxy'] == galaxy)
                                                                                 & sf_flag_f606w])

        else:
            ax[variable_index, galaxy_index].scatter(x=distance_halpha[halpha_input['gal'] == galaxy],
                                                     y=output_halpha[variable][output_halpha['galaxy'] == galaxy],
                                                     color=halpha_palette[2], edgecolors='w', s=20)
            ax[variable_index+4, galaxy_index].scatter(x=distance_f275w[f275w_input['gal'] == galaxy],
                                                       y=output_f275w[variable][output_f275w['galaxy'] == galaxy],
                                                       color=f275w_palette[2], edgecolors='w', s=20)
            ax[variable_index+8, galaxy_index].scatter(x=distance_f606w[f606w_input['gal'] == galaxy],
                                                       y=output_f606w[variable][output_f606w['galaxy'] == galaxy],
                                                       color=f606w_palette[2], edgecolors='w', s=20)

            correlation_halpha = scipy.stats.linregress(x=distance_halpha[halpha_input['gal'] == galaxy],
                                                        y=output_halpha[variable][output_halpha['galaxy'] == galaxy])
            correlation_f275w = scipy.stats.linregress(x=distance_f275w[f275w_input['gal'] == galaxy],
                                                        y=output_f275w[variable][output_f275w['galaxy'] == galaxy])
            correlation_f606w = scipy.stats.linregress(x=distance_f606w[f606w_input['gal'] == galaxy],
                                                        y=output_f606w[variable][output_f606w['galaxy'] == galaxy])

        x_lim_halpha = ax[variable_index, galaxy_index].set_xlim()
        x_halpha = np.linspace(x_lim_halpha[0], x_lim_halpha[1])
        ax[variable_index, galaxy_index].plot(x_halpha, correlation_halpha.slope * x_halpha +
                                              correlation_halpha.intercept, '--k', lw=0.8)

        x_lim_f275w = ax[variable_index+4, galaxy_index].set_xlim()
        x_f275w = np.linspace(x_lim_f275w[0], x_lim_f275w[1])
        ax[variable_index+4, galaxy_index].plot(x_f275w, correlation_f275w.slope * x_f275w +
                                                correlation_f275w.intercept, '--k', lw=0.8)

        x_lim_f606w = ax[variable_index+8, galaxy_index].set_xlim()
        x_f606w = np.linspace(x_lim_f606w[0], x_lim_f606w[1])
        ax[variable_index+8, galaxy_index].plot(x_f606w, correlation_f606w.slope * x_f606w +
                                                correlation_f606w.intercept, '--k', lw=0.8)

        for offset in [0, 4, 8]:
            ax[variable_index + offset, 0].set_ylabel(labels[variable_index], fontsize=10)

for i in range(len(galaxy_list)):
    ax[0, i].set_title(galaxy_list[i], fontsize=13)

fig.supxlabel(r'$\log\,d/R_e$', fontsize=13)

sns.despine()

fig.tight_layout()
fig.subplots_adjust(bottom=0.045, hspace=0.05, wspace=0.03, top=0.98)

plt.savefig('/home/ariel/Workspace/GASP/HST/Analysis/Plots/Paper/distance_galbygal.pdf')
plt.savefig('/home/ariel/Workspace/GASP/HST/Analysis/Plots/Paper/distance_galbygal.png', dpi=300)

# pivotted_halpha = df_halpha.pivot('Galaxy', 'Variable', 'Slope')
# pivotted_f275w = df_f275w.pivot('Galaxy', 'Variable', 'Slope')
# pivotted_f606w = df_f606w.pivot('Galaxy', 'Variable', 'Slope')
#
# fig, ax = plt.subplots(1, 3, figsize=(10, 5.5))
#
# sns.heatmap(pivotted_halpha, annot=True, cmap=sns.light_palette('goldenrod', as_cmap=True).reversed(),
#             vmin=-0.5, vmax=0, ax=ax[0], linewidths=0.5, linecolor='w', cbar=False)
# sns.heatmap(pivotted_f275w, annot=True, cmap=sns.light_palette('mediumvioletred', as_cmap=True).reversed(),
#             vmin=-0.5, vmax=0, ax=ax[1], linewidths=0.5, linecolor='w', cbar=False)
# sns.heatmap(pivotted_f606w, annot=True, cmap=sns.light_palette('indigo', as_cmap=True).reversed(),
#             vmin=-0.5, vmax=0, ax=ax[2], linewidths=0.5, linecolor='w', cbar=False)
#
# plt.show()

# evil: output_f606w[(output_f606w['Av']>1) & (output_f606w['galaxy']=='JO201') & sf_flag_f606w]

# correlation_list_halpha = []
# correlation_list_f275w = []
# correlation_list_f606w = []
#
# for variable in variables:
#
#     if variable == 'Av':
#         correlation_halpha = scipy.stats.linregress(x=np.log10(halpha_input['dist_kpc'] / halpha_input['galaxy_radius_kpc'])[halpha_frac_flag],
#                                                     y=output_halpha[variable][halpha_frac_flag])
#         correlation_f275w = scipy.stats.linregress(x=np.log10(f275w_input['dist_kpc'] / f275w_input['galaxy_radius_kpc'])[f275w_frac_flag],
#                                                    y=output_f275w[variable][f275w_frac_flag])
#         correlation_f606w = scipy.stats.linregress(x=np.log10(f606w_input['dist_kpc'] / f606w_input['galaxy_radius_kpc'])[f606w_frac_flag],
#                                                    y=output_f606w[variable][f606w_frac_flag])
#
#     elif variable == 'sfr':
#         correlation_halpha = scipy.stats.linregress(x=np.log10(halpha_input['dist_kpc'] / halpha_input['galaxy_radius_kpc'])[sf_flag_halpha],
#                                                     y=output_halpha[variable][sf_flag_halpha])
#         correlation_f275w = scipy.stats.linregress(x=np.log10(f275w_input['dist_kpc'] / f275w_input['galaxy_radius_kpc'])[sf_flag_f275w],
#                                                    y=output_f275w[variable][sf_flag_f275w])
#         correlation_f606w = scipy.stats.linregress(x=np.log10(f606w_input['dist_kpc'] / f606w_input['galaxy_radius_kpc'])[sf_flag_f606w],
#                                                    y=output_f606w[variable][sf_flag_f606w])
#
#     else:
#         correlation_halpha = scipy.stats.linregress(x=np.log10(halpha_input['dist_kpc'] / halpha_input['galaxy_radius_kpc']),
#                                                     y=output_halpha[variable])
#         correlation_f275w = scipy.stats.linregress(x=np.log10(f275w_input['dist_kpc'] / f275w_input['galaxy_radius_kpc']),
#                                                    y=output_f275w[variable])
#         correlation_f606w = scipy.stats.linregress(x=np.log10(f606w_input['dist_kpc'] / f606w_input['galaxy_radius_kpc']),
#                                                    y=output_f606w[variable])
#
#     correlation_list_halpha.append(['All Galaxies', variable, correlation_halpha])
#     correlation_list_f275w.append(['All Galaxies', variable, correlation_f275w])
#     correlation_list_f606w.append(['All Galaxies', variable, correlation_f606w])
#
#
# for galaxy in ['JO175', 'JO201', 'JO204', 'JO206', 'JW39', 'JW100']:
#
#     for variable in variables:
#
#         if variable == 'Av':
#             correlation_halpha = scipy.stats.linregress(x=np.log10(halpha_input['dist_kpc']/halpha_input['galaxy_radius_kpc'])[(halpha_input['gal'] == galaxy) & sf_flag_halpha],
#                                                         y=output_halpha[variable][(output_halpha['galaxy'] == galaxy) & sf_flag_halpha])
#             correlation_f275w = scipy.stats.linregress(x=np.log10(f275w_input['dist_kpc']/f275w_input['galaxy_radius_kpc'])[(f275w_input['gal'] == galaxy) & sf_flag_f275w],
#                                                         y=output_f275w[variable][(output_f275w['galaxy'] == galaxy) & sf_flag_f275w])
#             correlation_f606w = scipy.stats.linregress(x=np.log10(f606w_input['dist_kpc']/f606w_input['galaxy_radius_kpc'])[(f606w_input['gal'] == galaxy) & sf_flag_f606w],
#                                                         y=output_f606w[variable][(output_f606w['galaxy'] == galaxy) & sf_flag_f606w])
#
#         elif variable == 'sfr':
#             correlation_halpha = scipy.stats.linregress(x=np.log10(halpha_input['dist_kpc']/halpha_input['galaxy_radius_kpc'])[(halpha_input['gal'] == galaxy) & sf_flag_halpha],
#                                                         y=output_halpha[variable][(output_halpha['galaxy'] == galaxy) & sf_flag_halpha])
#             correlation_f275w = scipy.stats.linregress(x=np.log10(f275w_input['dist_kpc']/f275w_input['galaxy_radius_kpc'])[(f275w_input['gal'] == galaxy) & sf_flag_f275w],
#                                                         y=output_f275w[variable][(output_f275w['galaxy'] == galaxy) & sf_flag_f275w])
#             correlation_f606w = scipy.stats.linregress(x=np.log10(f606w_input['dist_kpc']/f606w_input['galaxy_radius_kpc'])[(f606w_input['gal'] == galaxy) & sf_flag_f606w],
#                                                         y=output_f606w[variable][(output_f606w['galaxy'] == galaxy) & sf_flag_f606w])
#
#         else:
#             correlation_halpha = scipy.stats.linregress(x=np.log10(halpha_input['dist_kpc'] / halpha_input['galaxy_radius_kpc'])[halpha_input['gal'] == galaxy],
#                                                         y=output_halpha[variable][output_halpha['galaxy'] == galaxy])
#             correlation_f275w = scipy.stats.linregress(x=np.log10(f275w_input['dist_kpc']/f275w_input['galaxy_radius_kpc'])[f275w_input['gal'] == galaxy],
#                                                         y=output_f275w[variable][output_f275w['galaxy'] == galaxy])
#             correlation_f606w = scipy.stats.linregress(x=np.log10(f606w_input['dist_kpc']/f606w_input['galaxy_radius_kpc'])[f606w_input['gal'] == galaxy],
#                                                         y=output_f606w[variable][output_f606w['galaxy'] == galaxy])
#
#         correlation_list_halpha.append([galaxy, variable, correlation_halpha])
#         correlation_list_f275w.append([galaxy, variable, correlation_f275w])
#         correlation_list_f606w.append([galaxy, variable, correlation_f606w])
#
#
# dict_halpha = {'Galaxy': [correlation_list_halpha[i][0] for i in range(len(correlation_list_halpha))],
#                'Variable': [correlation_list_halpha[i][1] for i in range(len(correlation_list_halpha))],
#                'Pearson_r': [correlation_list_halpha[i][2][2] for i in range(len(correlation_list_halpha))],
#                }
#
# dict_f275w = {'Galaxy': [correlation_list_f275w[i][0] for i in range(len(correlation_list_f275w))],
#               'Variable': [correlation_list_f275w[i][1] for i in range(len(correlation_list_f275w))],
#               'Pearson_r': [correlation_list_f275w[i][2][2] for i in range(len(correlation_list_f275w))],
#               }
#
# dict_f606w = {'Galaxy': [correlation_list_f606w[i][0] for i in range(len(correlation_list_f606w))],
#               'Variable': [correlation_list_f606w[i][1] for i in range(len(correlation_list_f606w))],
#               'Pearson_r': [correlation_list_f606w[i][2][2] for i in range(len(correlation_list_f606w))],
#               }
#
# df_halpha = pd.DataFrame(dict_halpha)
# df_f275w = pd.DataFrame(dict_f275w)
# df_f606w = pd.DataFrame(dict_f606w)
#
# pivotted_halpha = df_halpha.pivot('Galaxy', 'Variable', 'Pearson_r')
# pivotted_f275w = df_f275w.pivot('Galaxy', 'Variable', 'Pearson_r')
# pivotted_f606w = df_f606w.pivot('Galaxy', 'Variable', 'Pearson_r')
#
# pivotted_halpha = pivotted_halpha[variables]
# pivotted_f275w = pivotted_f275w[variables]
# pivotted_f606w = pivotted_f606w[variables]
#
# fig, ax = plt.subplots(1, 3, figsize=(9.25, 4.5), sharey=True)
#
# sns.heatmap(pivotted_halpha, annot=True, cmap=sns.light_palette('goldenrod', as_cmap=True).reversed(),
#             ax=ax[0], linewidths=0.5, linecolor='w', cbar=False)
# sns.heatmap(pivotted_f275w, annot=True, cmap=sns.light_palette('mediumvioletred', as_cmap=True).reversed(),
#             ax=ax[1], linewidths=0.5, linecolor='w', cbar=False)
# sns.heatmap(pivotted_f606w, annot=True, cmap=sns.light_palette('indigo', as_cmap=True).reversed(),
#             ax=ax[2], linewidths=0.5, linecolor='w', cbar=False)
#
# for axes in ax:
#     axes.set_ylabel(' ')
#     axes.set_xlabel(' ')
#
# sns.despine(left=True, bottom=True)
#
# fig.tight_layout()
# fig.subplots_adjust(wspace=0.05)
#
# plt.show()
# plt.close()
