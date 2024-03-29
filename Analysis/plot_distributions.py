"""

ariel@oapd
15/08/2022

Plots distribution of ages and masses for clumps in the tails

H-alpha: "goldenrod" (0.8549019607843137, 0.6470588235294118, 0.12549019607843137, 1.0)
UV: "mediumvioletred" (0.7803921568627451, 0.08235294117647059, 0.5215686274509804, 1.0)
Optical: "indigo" (0.29411764705882354, 0.0, 0.5098039215686274, 1.0)

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from astropy.table import Table
import seaborn as sns
from toolbox import plot_tools
import pandas as pd

sns.set_style('ticks')
palette = sns.diverging_palette(220, 20, n=7)
halpha_palette = sns.light_palette('goldenrod',  3)
f275w_palette = sns.light_palette('mediumvioletred', 3)
f606w_palette = sns.light_palette('indigo', 3)

mass_dict = {'JO201': 44194800000,
             'JO204': 54968402000,
             'JW100': 292875993000,
             'JW39': 164373004000,
             'JO175': 33957900300,
             'JO206': 77743301000
}

halpha_input = Table.read('/home/ariel/Workspace/GASP/HST/Data/halpha_bagpipes_input.fits')
f275w_input = Table.read('/home/ariel/Workspace/GASP/HST/Data/f275w_bagpipes_input.fits')
f606w_input = Table.read('/home/ariel/Workspace/GASP/HST/Data/f606w_bagpipes_input.fits')
optical_only_input = Table.read('/home/ariel/Workspace/GASP/HST/Data/optical_only_bagpipes_input.fits')

output_halpha = Table.read('/home/ariel/Workspace/GASP/HST/Data/halpha_dexp_logprior_single_bagpipes_results.fits')
output_f275w = Table.read('/home/ariel/Workspace/GASP/HST/Data/f275w_dexp_logprior_single_bagpipes_results.fits')
output_f606w = Table.read('/home/ariel/Workspace/GASP/HST/Data/f606w_dexp_logprior_single_bagpipes_results.fits')
output_optical_only = Table.read('/home/ariel/Workspace/GASP/HST/Data/optical_only_dexp_logprior_single_bagpipes_results.fits')

halpha_input = halpha_input[~output_halpha['bad_double_fit'] & ~output_halpha['bad_fit']]
output_halpha = output_halpha[~output_halpha['bad_double_fit'] & ~output_halpha['bad_fit']]

f275w_input = f275w_input[~output_f275w['bad_double_fit'] & ~output_f275w['bad_fit']]
output_f275w = output_f275w[~output_f275w['bad_double_fit'] & ~output_f275w['bad_fit']]

f606w_input = f606w_input[~output_f606w['bad_fit']]
output_f606w = output_f606w[~output_f606w['bad_fit']]

output_halpha['galaxy'] = output_halpha['galaxy'].astype(str)
output_f275w['galaxy'] = output_f275w['galaxy'].astype(str)
output_f606w['galaxy'] = output_f606w['galaxy'].astype(str)

output_halpha['Location'] = np.zeros_like(output_halpha['galaxy'])
output_halpha['Location'][halpha_input['tail_gal_flag'] == 0] = np.full((halpha_input['tail_gal_flag'] == 0).sum(),
                                                                        'Tail')
output_halpha['Location'][halpha_input['tail_gal_flag'] == 1] = np.full((halpha_input['tail_gal_flag'] == 1).sum(),
                                                                        'Extraplanar')

output_f275w['Location'] = np.zeros_like(output_f275w['galaxy'])
output_f275w['Location'][f275w_input['tail_gal_flag'] == 0] = np.full((f275w_input['tail_gal_flag'] == 0).sum(),
                                                                      'Tail')
output_f275w['Location'][f275w_input['tail_gal_flag'] == 1] = np.full((f275w_input['tail_gal_flag'] == 1).sum(),
                                                                      'Extraplanar')

output_f606w['Location'] = np.zeros_like(output_f606w['galaxy'])
output_f606w['Location'][f606w_input['tail_gal_flag'] == 0] = np.full((f606w_input['tail_gal_flag'] == 0).sum(),
                                                                      'Tail')
output_f606w['Location'][f606w_input['tail_gal_flag'] == 1] = np.full((f606w_input['tail_gal_flag'] == 1).sum(),
                                                                      'Extraplanar')

output_halpha['galaxy_mass'] = np.zeros_like(output_halpha['mwage'])
for i in range(len(output_halpha)):
    output_halpha['galaxy_mass'][i] = np.log10(mass_dict[output_halpha['galaxy'][i]])

output_f275w['galaxy_mass'] = np.zeros_like(output_f275w['mwage'])
for i in range(len(output_f275w)):
    output_f275w['galaxy_mass'][i] = np.log10(mass_dict[output_f275w['galaxy'][i]])

output_f606w['galaxy_mass'] = np.zeros_like(output_f606w['mwage'])
for i in range(len(output_f606w)):
    output_f606w['galaxy_mass'][i] = np.log10(mass_dict[output_f606w['galaxy'][i]])

output_halpha['mwage'] *= 1e3
output_f275w['mwage'] *= 1e3
output_f606w['mwage'] *= 1e3
output_optical_only['mwage'] *= 1e3

output_halpha['age'] *= 1e3
output_f275w['age'] *= 1e3
output_f606w['age'] *= 1e3
output_optical_only['age'] *= 1e3

# halpha_frac_flag = output_halpha['mass_frac_20'] > 0.1
# f275w_frac_flag = output_f275w['mass_frac_20'] > 0.1
# f606w_frac_flag = output_f606w['mass_frac_20'] > 0.1

halpha_frac_flag = output_halpha['late_decliner']
f275w_frac_flag = output_f275w['late_decliner']
f606w_frac_flag = output_f606w['late_decliner']

output_halpha['mass'] = output_halpha['stellar_mass']
output_f275w['mass'] = output_f275w['stellar_mass']
output_f606w['mass'] = output_f606w['stellar_mass']

output_halpha['sfr'] = np.log10(output_halpha['sfr'])
output_f275w['sfr'] = np.log10(output_f275w['sfr'])
output_f606w['sfr'] = np.log10(output_f606w['sfr'])

output_halpha['Av'] = output_halpha['Av'] * output_halpha['eta']
output_f275w['Av'] = output_f275w['Av'] * output_f275w['eta']
output_f606w['Av'] = output_f606w['Av'] * output_f606w['eta']

halpha_df = output_halpha['galaxy', 'clump_id', 'age', 'tau', 'mwage', 'mass', 'galaxy_mass', 'sfr', 'Av',
                          'metallicity', 'Location', 'mass_frac_20'].to_pandas()
f275w_df = output_f275w['galaxy', 'clump_id', 'age', 'tau', 'mwage', 'mass', 'galaxy_mass', 'sfr', 'Av',
                        'metallicity', 'Location', 'mass_frac_20'].to_pandas()
f606w_df = output_f606w['galaxy', 'clump_id', 'age', 'tau', 'mwage', 'mass', 'galaxy_mass', 'sfr', 'Av',
                        'metallicity', 'Location', 'mass_frac_20'].to_pandas()

halpha_df['detection'] = np.full(len(halpha_df), 'Halpha')
f275w_df['detection'] = np.full(len(f275w_df), 'F275W')
f606w_df['detection'] = np.full(len(f606w_df), 'F606W')

halpha_df.sort_values(by='galaxy_mass', inplace=True)
f275w_df.sort_values(by='galaxy_mass', inplace=True)
f606w_df.sort_values(by='galaxy_mass', inplace=True)

all_df = pd.concat([halpha_df, f275w_df, f606w_df])

# Now we plot!

variables = ['mwage', 'mass', 'sfr', 'Av']
labels = [r'$\langle\,t_\star\,\rangle_M\,\mathrm{[Myr]}$', r'$\log\,M_\star/M_\odot$',
          r'$\log \, \mathrm{SFR} \,\mathrm{[M_\odot/yr]}$', r'$A_V\,\mathrm{[mag]}$']

fig, ax = plt.subplots(4, 3, figsize=(9.25, 9), sharex='row')
fig_all, ax_all = plt.subplots(1, 4, figsize=(9.25, 5.5), sharey='row')

for i in range(4):

    distribution_variable = variables[i]
    label = labels[i]

    if distribution_variable == 'mwage':
        all_mwage = np.concatenate([output_halpha['mwage'].data, output_f275w['mwage'].data, output_f606w['mwage'].data])
        max_mwage = np.percentile(all_mwage, 95)
        print('MAX MWAGE', np.round(max_mwage))
        hist_range = [0, np.round(max_mwage)]
        ax[i, 0].annotate('$5\%$', xy=(0.85, 0.9), xycoords='axes fraction', fontsize=14, ha='left')
        ax[i, 0].arrow(x=0.85, y=0.83, dx=0.08, dy=0.0, color='k', head_length=0.03, head_width=0.05,
                       transform=ax[i, 0].transAxes)
        all_range = hist_range

    elif distribution_variable == 'mass':
        hist_range = [3, 8.4]

        inset_hist = ax[i, 0].inset_axes([0.53, 0.63, 0.47, 0.37])

        inset_hist.hist(output_halpha['sigma_m'][halpha_input['resolved_flag'].astype(bool)],
                        histtype='step', range=[5.45, 7.6],
                        bins=9, density=True, color='goldenrod')
        inset_hist.hist(output_f275w['sigma_m'][f275w_input['resolved_flag'].astype(bool)],
                        histtype='step', range=[5.45, 7.6],
                        bins=9, density=True, color='mediumvioletred')
        inset_hist.hist(output_f606w['sigma_m'][f606w_input['resolved_flag'].astype(bool)],
                        histtype='step', range=[5.45, 7.6],
                        bins=9, density=True, color='indigo')

        sns.despine(ax=inset_hist, left=True)
        inset_hist.set_yticks([])

        inset_hist.tick_params(labelsize=8)

        inset_hist.set_xlabel(r'$\log \Sigma_{M_\star}\,\mathrm{[M_\odot/kpc^2]}$', fontsize=10)

        all_range = [np.percentile(all_df[variables[i]], 1), np.percentile(all_df[variables[i]], 99)]

    elif distribution_variable == 'sfr':
        all_sfr = np.concatenate([output_halpha['sfr'], output_f275w['sfr'], output_f606w['sfr']])
        min_sfr = np.percentile(all_sfr, 15)
        print('MIN SFR', np.round(min_sfr, 1))
        hist_range = [np.round(min_sfr, 1), 0]
        ax[i, 0].annotate('$15\%$', xy=(0.17, 0.9), xycoords='axes fraction', fontsize=14, ha='right')
        ax[i, 0].arrow(x=0.17, y=0.83, dx=-0.1, dy=0.0, color='k', head_length=0.03, head_width=0.05,
                       transform=ax[i, 0].transAxes)
    elif distribution_variable == 'Av':
        output_halpha = output_halpha[halpha_frac_flag]
        output_f275w = output_f275w[f275w_frac_flag]
        output_f606w = output_f606w[f606w_frac_flag]
        halpha_input = halpha_input[halpha_frac_flag]
        f275w_input = f275w_input[f275w_frac_flag]
        f606w_input = f606w_input[f606w_frac_flag]
        all_Av = np.concatenate([output_halpha['Av'], output_f275w['Av'], output_f606w['Av']])
        max_Av = np.percentile(all_Av, 95)
        print('MAX AV', np.round(max_Av, 1))
        hist_range = [0, np.round(max_Av, 1)]
        ax[i, 0].annotate('$5\%$', xy=(0.83, 0.9), xycoords='axes fraction', fontsize=14, ha='left')
        ax[i, 0].arrow(x=0.85, y=0.82, dx=0.08, dy=0.0, color='k', head_length=0.03, head_width=0.05,
                       transform=ax[i, 0].transAxes)
        all_df = all_df[all_df['mass_frac_20'] > 0.1]
        all_range = [0, np.percentile(all_df[variables[i]], 98.5)]

    print('Halpha', labels[i], np.round(np.median(output_halpha[variables[i]])),
          np.round(np.median(output_halpha[variables[i]][output_halpha['Location'] == 'Tail'])),
          np.round(np.median(output_halpha[variables[i]][output_halpha['Location'] == 'Extraplanar'])))
    print('f275w', labels[i], np.round(np.median(output_f275w[variables[i]])),
          np.round(np.median(output_f275w[variables[i]][output_f275w['Location'] == 'Tail'])),
          np.round(np.median(output_f275w[variables[i]][output_f275w['Location'] == 'Extraplanar'])))
    print('f606w', labels[i], np.median(output_f606w[variables[i]]))

    print('5th percentile:')
    print('Halpha', labels[i], np.round(np.percentile(output_halpha[variables[i]], 5)),
          np.round(np.percentile(output_halpha[variables[i]][output_halpha['Location'] == 'Tail'], 5)),
          np.round(np.percentile(output_halpha[variables[i]][output_halpha['Location'] == 'Extraplanar'], 5)))
    print('f275w', labels[i], np.round(np.percentile(output_f275w[variables[i]], 5)),
          np.round(np.percentile(output_f275w[variables[i]][output_f275w['Location'] == 'Tail'], 5)),
          np.round(np.percentile(output_f275w[variables[i]][output_f275w['Location'] == 'Extraplanar'], 5)))
    print('f606w', labels[i], np.round(np.percentile(output_f606w[variables[i]], 5)))

    print('95th percentile:')
    print('Halpha', labels[i], np.round(np.percentile(output_halpha[variables[i]], 95)),
          np.round(np.percentile(output_halpha[variables[i]][output_halpha['Location'] == 'Tail'], 95)),
          np.round(np.percentile(output_halpha[variables[i]][output_halpha['Location'] == 'Extraplanar'], 95)))
    print('f275w', labels[i], np.round(np.percentile(output_f275w[variables[i]], 95)),
          np.round(np.percentile(output_f275w[variables[i]][output_f275w['Location'] == 'Tail'], 95)),
          np.round(np.percentile(output_f275w[variables[i]][output_f275w['Location'] == 'Extraplanar'], 95)))
    print('f606w', labels[i], np.round(np.percentile(output_f606w[variables[i]], 95)))

    if variables[i] == 'sfr':
        print('!!!!!!!!!!!!!!! \n FLAGGED SFR!!! \n !!!!!!!!!!!!!!!')
        f275w_sfr_flag = output_f275w['sfr'] > -4.5
        print(1-f275w_sfr_flag.sum()/len(f275w_sfr_flag))
        f606w_sfr_flag = output_f606w['sfr'] > -4.5
        print(1-f606w_sfr_flag.sum()/len(f606w_sfr_flag))
        print('f275w', labels[i], np.median(output_f275w[variables[i]][f275w_sfr_flag]),
              np.median(output_f275w[variables[i]][(output_f275w['Location'] == 'Tail') & f275w_sfr_flag]),
              np.median(output_f275w[variables[i]][(output_f275w['Location'] == 'Extraplanar') & f275w_sfr_flag]))
        print('f606w', labels[i], np.median(output_f606w[variables[i]][f606w_sfr_flag]))

    ax[i, 0].hist(output_f606w[distribution_variable], label='Star-Forming Complexes',
                  histtype='step',
                  bins=14, range=hist_range, density=True, color='indigo', lw=1.5)
    ax[i, 0].hist(output_f275w[distribution_variable], label='F275W Clumps', histtype='step',
                  bins=14, range=hist_range, density=True, color='mediumvioletred', lw=1.5)
    ax[i, 0].hist(output_halpha[distribution_variable], label=r'$H\alpha$ Clumps',
                  histtype='step',
                  bins=14, range=hist_range, density=True, color='goldenrod', lw=1.5)

    ax[i, 1].hist(output_f606w[distribution_variable][f606w_input['tail']], label='Star-Forming Complexes', histtype='step',
                  bins=14, range=hist_range, density=True, color='indigo', lw=1.5)
    ax[i, 1].hist(output_f275w[distribution_variable][f275w_input['tail']], label='F275W Clumps', histtype='step',
                  bins=14, range=hist_range, density=True, color='mediumvioletred', lw=1.5)
    ax[i, 1].hist(output_halpha[distribution_variable][halpha_input['tail']], label=r'$H\alpha$ Clumps', histtype='step',
                  bins=14, range=hist_range, density=True, color='goldenrod', lw=1.5)

    ax[i, 2].hist(output_f275w[distribution_variable][f275w_input['extraplanar']], label='F275W Clumps', histtype='step',
                  bins=14, range=hist_range, density=True, color='mediumvioletred', lw=1.5)
    ax[i, 2].hist(output_halpha[distribution_variable][halpha_input['extraplanar']], label=r'$H\alpha$ Clumps', histtype='step',
                  bins=14, range=hist_range, density=True, color='goldenrod', lw=1.5)

    ax[3, 1].legend(frameon=False, fontsize=15, loc='upper center', bbox_to_anchor=(0.5, -0.25), ncol=3,
                    shadow=True)

    ax[0, 0].set_title('Whole Sample', fontsize=13)
    ax[0, 1].set_title('Tails', fontsize=13)
    ax[0, 2].set_title('Extraplanar Regions', fontsize=13)

    # ax[i, 0].set_ylabel(r'Density', fontsize=15)
    ax[i, 0].set_xlabel(label, fontsize=13)
    ax[i, 1].set_xlabel(label, fontsize=13)
    ax[i, 2].set_xlabel(label, fontsize=13)

    ax[i, 0].set_xlim(hist_range[0], hist_range[1])

    if variables[i] == 'sfr':
        flag = all_df['sfr'] > -4.5
        box = sns.boxplot(x=variables[i], y='galaxy', data=all_df[flag], hue='detection', ax=ax_all[i], whis=[15, 85],
                          fliersize=0, palette=[halpha_palette[1], f275w_palette[1], f606w_palette[1]], saturation=1,
                          linewidth=0.85, width=0.6)
        ax_all[i].set_xlim(-4, np.percentile(all_df['sfr'][flag], 95))

        legend_handles, legend_labels = ax_all[i].get_legend_handles_labels()
        new_labels = [r'$H\alpha$ Clumps', 'F275W Clumps', 'Star-Forming Complexes']
        legend = box.legend(legend_handles, new_labels, frameon=False, loc='upper center', bbox_to_anchor=(-0.15, 1.075),
                            ncol=3, shadow=True, fontsize=13)

    else:
        sns.boxplot(x=variables[i], y='galaxy', data=all_df, hue='detection', ax=ax_all[i], whis=[15, 85], fliersize=0,
                    palette=[halpha_palette[1], f275w_palette[1], f606w_palette[1]], saturation=1, linewidth=0.85,
                    width=0.6)
        ax_all[i].legend_.remove()

        ax_all[i].set_xlim(all_range[0], all_range[1])

    ax_all[i].set_ylabel(' ')
    ax_all[i].set_xlabel(labels[i], fontsize=13)

    ax_all[i].tick_params(labelsize=15, axis='y')
    ax_all[i].tick_params(labelsize=12, axis='x')

fig.supylabel(r'Density', fontsize=15)

sns.despine(fig=fig)
sns.despine(fig=fig_all)

fig.tight_layout()
fig.subplots_adjust(left=0.075, bottom=0.095, wspace=0.1745, hspace=0.3, top=0.975, right=0.99)

fig_all.subplots_adjust(left=0.075, bottom=0.1, wspace=0.02, hspace=0.0, top=0.96, right=0.99)

fig.savefig('/home/ariel/Workspace/GASP/HST/Analysis/Plots/Paper/all_hist.png', dpi=350)
fig.savefig('/home/ariel/Workspace/GASP/HST/Analysis/Plots/Paper/all_hist.pdf')
fig_all.savefig('/home/ariel/Workspace/GASP/HST/Analysis/Plots/Paper/all_galbygal.pdf')
fig_all.savefig('/home/ariel/Workspace/GASP/HST/Analysis/Plots/Paper/all_galbygal.jpg', dpi=350)


# fig, ax = plt.subplots(12, 1, figsize=(9, 10.5), sharex=True)
#
# for i in range(4):
#
#     distribution_variable = variables[i]
#     label = labels[i]
#
#     p1 = sns.boxplot(x='galaxy', y=distribution_variable, data=halpha_df, whis=[15, 85], width=.6, palette=halpha_palette[1:3]
#                      , fliersize=0, ax=ax[0+3*i], hue='Location', orient='vertical')
#     p2 = sns.boxplot(x='galaxy', y=distribution_variable, data=f275w_df, whis=[15, 85], width=.6, palette=f275w_palette[1:3]
#                      , fliersize=0, ax=ax[1+3*i], hue='Location', orient='vertical')
#     p3 = sns.boxplot(x='galaxy', y=distribution_variable, data=f606w_df, whis=[15, 85], width=.6, palette=f606w_palette[1:3]
#                      , fliersize=0, ax=ax[2+3*i], hue='Location', orient='vertical')
#
#     if distribution_variable == 'mwage':
#         hist_range = [0, 200]
#     elif distribution_variable == 'mass':
#         hist_range = [3, 7.8]
#     elif distribution_variable == 'sfr':
#         hist_range = [-6, 0]
#     else:
#         hist_range = [0, 1]
#
#     p1.set_ylim(np.percentile(output_halpha[distribution_variable], 2.5),
#                 np.percentile(output_halpha[distribution_variable], 97))
#     p2.set_ylim(np.percentile(output_f275w[distribution_variable], 2.5),
#                 np.percentile(output_f275w[distribution_variable], 97))
#     p3.set_ylim(np.percentile(output_f606w[distribution_variable][~np.isnan(output_f606w[distribution_variable])], 2.5),
#                 np.percentile(output_f606w[distribution_variable][~np.isnan(output_f606w[distribution_variable])], 95))
#
#     p1.set_ylabel(' ')
#     p2.set_ylabel(label, fontsize=20)
#     p3.set_ylabel(' ')
#
#     p1.legend(frameon=False, fontsize=12, loc=1)
#     p2.legend(frameon=False, fontsize=12, loc=1)
#     p3.legend(frameon=False, fontsize=12, loc=1)
#
#     ax[2+3*i].tick_params(labelsize=20, axis='x')
#     # ax[0].tick_params(labelsize=14, axis='x')
#     # ax[1].tick_params(labelsize=14, axis='x')
#     # ax[2].tick_params(labelsize=14, axis='x')
#     # ax[1].tick_params(left=False, labelleft=False)
#     # ax[2].tick_params(left=False, labelleft=False)
#     sns.despine()
#
#
# ax[11].set_xlabel(' ')
# fig.subplots_adjust(top=0.98, bottom=0.05, left=0.085, right=0.9935, hspace=0.0, wspace=0.0)
#
# plt.savefig('/home/ariel/Workspace/GASP/HST/Analysis/Plots/Paper/all_galbygal.png', dpi=300)


# # Onset of SF for complexes
#
# fig, ax = plt.subplots(1, 1, figsize=(5.5, 7.5), sharex=True)
#
# f606w_palette = sns.light_palette('indigo', 4)
#
# box = sns.boxplot(y='galaxy', x='age', data=f606w_df, whis=[1, 99], width=.6, color=f606w_palette[1]
#                   , fliersize=0, ax=ax, orient='horizontal', saturation=1, linewidth=0.85)
#
# galaxy_list = ['JW100', 'JW39', 'JO206', 'JO204', 'JO201', 'JO175']
#
# lines = ax.get_lines()
# boxes = [c for c in ax.get_children() if type(c).__name__ == 'PathPatch']
# lines_per_box = int(len(lines) / len(boxes))
# for line in lines[1:len(lines):lines_per_box]:
#     print(line.get_xdata())
#     x = line.get_xdata()[1]
#     y = line.get_ydata()[0]
#     ax.annotate('$%i \,\mathrm{Myr}$' % np.round(x, 0), xy=(x+10, y), va='center', fontsize=13.5)
#
# ax.tick_params(labelsize=18, axis='y')
# ax.tick_params(labelsize=12, axis='x')
# ax.set_ylabel(' ')
# ax.set_xlabel('$t_0^\star \; \mathrm{[Myr]}$', fontsize=20)
#
# ax.set_xlim(0, 445)
#
# sns.despine()
#
# fig.subplots_adjust(top=0.981,
#                     bottom=0.091,
#                     left=0.193,
#                     right=0.962,
#                     hspace=0.2,
#                     wspace=0.2)
#
# plt.savefig('/home/ariel/Workspace/GASP/HST/Analysis/Plots/Paper/onset_age.png', dpi=350)
# plt.savefig('/home/ariel/Workspace/GASP/HST/Analysis/Plots/Paper/onset_age.pdf')




# grid = sns.jointplot(x=output_optical_only['stellar_mass'].tolist(), y=output_optical_only['age'].tolist())
# grid.ax_joint.set_xlabel(labels[1], fontsize=20)
# grid.ax_joint.set_ylabel(r'$t_0\,\mathrm{[Myr]}$', fontsize=20)
#
# grid.fig.subplots_adjust(top=0.99,
#                          bottom=0.102,
#                          left=0.116,
#                          right=0.975,
#                          hspace=0.2,
#                          wspace=0.2)
#
# plt.savefig('/home/ariel/Workspace/GASP/HST/Analysis/Plots/Paper/joint_age.png', dpi=350)
#
# grid = sns.jointplot(x=output_optical_only['stellar_mass'].tolist(), y=output_optical_only['mwage'].tolist())
# grid.ax_joint.set_xlabel(labels[1], fontsize=20)
# grid.ax_joint.set_ylabel(labels[0], fontsize=20)
#
# grid.fig.subplots_adjust(top=0.99,
#                          bottom=0.102,
#                          left=0.116,
#                          right=0.975,
#                          hspace=0.2,
#                          wspace=0.2)
#
# plt.savefig('/home/ariel/Workspace/GASP/HST/Analysis/Plots/Paper/joint_mwage.png', dpi=350)
#
# beta = np.log10(f275w_input['F336W']/f275w_input['F275W'])/np.log10(3360/2750)
#
# lowbeta = beta < -3.34
#
# fig, ax = plt.subplots(figsize=(6, 5), sharex=True)
#
# sns.kdeplot(output_f275w['stellar_mass'][~lowbeta].tolist(), (output_f275w['mwage'][~lowbeta]).tolist(),
#             color=palette[0], levels=6)
# sns.kdeplot(output_f275w['stellar_mass'][lowbeta].tolist(), (output_f275w['mwage'][lowbeta]).tolist(),
#             color=palette[-1], levels=6)
#
# ax.annotate(r'$\beta>-3.34$', xy=(0.75, 0.92), fontsize=20, color=palette[0], xycoords='axes fraction')
# ax.annotate(r'$\beta<-3.34$', xy=(0.75, 0.84), fontsize=20, color=palette[-1], xycoords='axes fraction')
#
# ax.set_xlabel(r'$\log\,M/M_\odot$', fontsize=20)
# ax.set_ylabel(r'$\langle\,t_\star\,\rangle_M\,\mathrm{[Myr]}$', fontsize=20)
#
# ax.set_ylim(-30, 250)
#
# fig.tight_layout()
#
# plt.savefig('/home/ariel/Workspace/GASP/HST/Analysis/Plots/Paper/bad_beta.png', dpi=300)

