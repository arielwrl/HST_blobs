"""

ariel@viapalermo
31/03/2023

Some trends with axial ratio

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from astropy.table import Table
import seaborn as sns
from toolbox import plot_tools

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

output_halpha_unresolved = output_halpha[(~halpha_input['disk']) & ((halpha_input['level'] == 0) | (halpha_input['leaf_flag'] == 1)) & (~halpha_input['resolved_flag'])]
output_f275w_unresolved = output_f275w[(~f275w_input['disk']) & ((f275w_input['level'] == 0) | (f275w_input['leaf_flag'] == 1)) & (~f275w_input['resolved_flag'])]
output_f606w_unresolved = output_f606w[f606w_input['resolved_flag'] == False]

output_halpha = output_halpha[(~halpha_input['disk']) & ((halpha_input['level'] == 0) | (halpha_input['leaf_flag'] == 1)) & (halpha_input['resolved_flag'])]
output_f275w = output_f275w[(~f275w_input['disk']) & ((f275w_input['level'] == 0) | (f275w_input['leaf_flag'] == 1)) & (f275w_input['resolved_flag'])]

flag_606 = (f606w_input['resolved_flag'] == True) & (output_f606w['stellar_mass'] > 3)
output_f606w = output_f606w[flag_606]

halpha_input_unresolved = halpha_input[(~halpha_input['disk']) & ((halpha_input['level'] == 0) | (halpha_input['leaf_flag'] == 1)) & (~halpha_input['resolved_flag'])]
f275w_input_unresolved = f275w_input[(~f275w_input['disk']) & ((f275w_input['level'] == 0) | (f275w_input['leaf_flag'] == 1)) & (~f275w_input['resolved_flag'])]
f606w_input_unresolved = f606w_input[f606w_input['resolved_flag'] == False]

halpha_input = halpha_input[(~halpha_input['disk']) & ((halpha_input['level'] == 0) | (halpha_input['leaf_flag'] == 1)) & (halpha_input['resolved_flag'])]
f275w_input = f275w_input[(~f275w_input['disk']) & ((f275w_input['level'] == 0) | (f275w_input['leaf_flag'] == 1)) & (f275w_input['resolved_flag'])]
f606w_input = f606w_input[flag_606]

# fig, ax = plt.subplots(2, 3, figsize=(12.5, 7.5), sharey='row', sharex='col')
#
# tail = halpha_input['tail']
# extraplanar = halpha_input['extraplanar']
#
# ax[0, 0].scatter(halpha_input['axial_ratio'][extraplanar], np.log10(output_halpha['mwage'])[extraplanar],
#                  color=halpha_palette[3], edgecolors='w', label='Extraplanar')
# ax[0, 0].scatter(halpha_input['axial_ratio'][tail], np.log10(output_halpha['mwage'])[tail], color=halpha_palette[1],
#                  edgecolors='w', label='Tail', s=60)
# plot_tools.plot_median_in_bins(halpha_input['axial_ratio'], np.log10(output_halpha['mwage']), nbins=6,
#                                percentile_style='errorbar', color='k', percentiles_color='k', ax=ax[0, 0],
#                                point_size=50, point_edgewidths=1, point_edgecolors='w')
# sns.kdeplot(halpha_input_unresolved['axial_ratio'].data.tolist(), np.log10(output_halpha_unresolved['mwage']).data.tolist(), ax=ax[0, 0],
#             fill=False, color='0.8', zorder=1, levels=5)
#
#
# ax[1, 0].scatter(halpha_input['axial_ratio'][extraplanar], output_halpha['stellar_mass'][extraplanar],
#                  color=halpha_palette[3], edgecolors='w', label='Extraplanar', s=60)
# ax[1, 0].scatter(halpha_input['axial_ratio'][tail], output_halpha['stellar_mass'][tail], color=halpha_palette[1],
#                  edgecolors='w', label='Tail', s=60)
# plot_tools.plot_median_in_bins(halpha_input['axial_ratio'], output_halpha['stellar_mass'], nbins=6,
#                                percentile_style='errorbar', color='k', percentiles_color='k', ax=ax[1, 0],
#                                point_size=50, point_edgewidths=1, point_edgecolors='w')
# sns.kdeplot(halpha_input_unresolved['axial_ratio'].data.tolist(), output_halpha_unresolved['stellar_mass'].data.tolist(), ax=ax[1, 0],
#             fill=False, color='0.8', zorder=1, levels=5)
#
# tail = f275w_input['tail']
# extraplanar = f275w_input['extraplanar']
#
# ax[0, 1].scatter(f275w_input['axial_ratio'][extraplanar], np.log10(output_f275w['mwage'])[extraplanar], s=60,
#                  color=f275w_palette[3], edgecolors='w', label='Extraplanar')
# ax[0, 1].scatter(f275w_input['axial_ratio'][tail], np.log10(output_f275w['mwage'])[tail], color=f275w_palette[1],
#                  edgecolors='w', label='Tail', s=60)
# plot_tools.plot_median_in_bins(f275w_input['axial_ratio'], np.log10(output_f275w['mwage']), nbins=6,
#                                percentile_style='errorbar', color='k', percentiles_color='k', ax=ax[0, 1], point_size=50, point_edgewidths=1, point_edgecolors='w')
# sns.kdeplot(f275w_input_unresolved['axial_ratio'].data.tolist(), np.log10(output_f275w_unresolved['mwage']).data.tolist(), ax=ax[0, 1],
#             fill=False, color='0.8', zorder=1, levels=5)
#
# ax[1, 1].scatter(f275w_input['axial_ratio'][extraplanar], output_f275w['stellar_mass'][extraplanar], color=f275w_palette[3], edgecolors='w', s=60,
#                  label='Extraplanar')
# ax[1, 1].scatter(f275w_input['axial_ratio'][tail], output_f275w['stellar_mass'][tail], color=f275w_palette[1], edgecolors='w', label='Tail', s=60)
# plot_tools.plot_median_in_bins(f275w_input['axial_ratio'], output_f275w['stellar_mass'], nbins=6,
#                                percentile_style='errorbar', color='k', percentiles_color='k', ax=ax[1, 1],
#                                point_size=50, point_edgewidths=1, point_edgecolors='w')
# sns.kdeplot(f275w_input_unresolved['axial_ratio'].data.tolist(), output_f275w_unresolved['stellar_mass'].data.tolist(), ax=ax[1, 1],
#             fill=False, color='0.8', zorder=1, levels=5)
#
# halpha_flag = f606w_input['halpha_match']
#
# ax[0, 2].scatter(f606w_input['axial_ratio'][~halpha_flag], np.log10(output_f606w['mwage'][~halpha_flag]),
#                  color=f606w_palette[2], edgecolors='w', s=60, label=r'No $H\alpha$ Emission')
# ax[0, 2].scatter(f606w_input['axial_ratio'][halpha_flag], np.log10(output_f606w['mwage'][halpha_flag]),
#                  color=f606w_palette[2], edgecolors='goldenrod', linewidths=2, s=60, label=r'With $H\alpha$ Emission')
# plot_tools.plot_median_in_bins(f606w_input['axial_ratio'], np.log10(output_f606w['mwage']), nbins=5,
#                                percentile_style='errorbar', color='k', percentiles_color='k', ax=ax[0, 2],
#                                point_size=50, point_edgewidths=1, point_edgecolors='w')
# sns.kdeplot(f606w_input_unresolved['axial_ratio'].data.tolist(), np.log10(output_f606w_unresolved['mwage']).data.tolist(), ax=ax[0, 2],
#             fill=False, color='0.8', zorder=1, levels=5)
#
# ax[1, 2].scatter(f606w_input['axial_ratio'][~halpha_flag], output_f606w['stellar_mass'][~halpha_flag],
#                  color=f606w_palette[2], edgecolors='w', s=60)
# ax[1, 2].scatter(f606w_input['axial_ratio'][halpha_flag], output_f606w['stellar_mass'][halpha_flag],
#                  color=f606w_palette[2], edgecolors='goldenrod', linewidths=2, s=60)
# plot_tools.plot_median_in_bins(f606w_input['axial_ratio'], output_f606w['stellar_mass'], nbins=5,
#                                percentile_style='errorbar', color='k', percentiles_color='k', ax=ax[1, 2],
#                                point_size=50, point_edgewidths=1, point_edgecolors='w')
# sns.kdeplot(f606w_input_unresolved['axial_ratio'].data.tolist(), output_f606w_unresolved['stellar_mass'].data.tolist(), ax=ax[1, 2],
#             fill=False, color='0.8', zorder=1, levels=5)
#
# sns.despine()
#
# ax[0, 0].legend(frameon=True, fontsize=14, loc=1)
# ax[0, 1].legend(frameon=True, fontsize=14, loc=1)
# ax[0, 2].legend(frameon=True, fontsize=14, loc=1)
#
# ax[1, 0].annotate(r'$H\alpha$', fontsize=20, xy=(0.02, 0.02), xycoords='axes fraction')
# ax[1, 1].annotate(r'F275W', fontsize=20, xy=(0.02, 0.02), xycoords='axes fraction')
# ax[1, 2].annotate(r'F606W', fontsize=20, xy=(0.02, 0.02), xycoords='axes fraction')
#
# ax[1, 1].set_xlabel(r'Axial Ratio', fontsize=20)
# ax[0, 0].set_ylabel(r'$\log\;\langle\,t_\star\,\rangle_M\,\mathrm{[Myr]}$', fontsize=20)
# ax[1, 0].set_ylabel(r'$\log\,M/M_\odot$', fontsize=20)
#
# fig.tight_layout()
#
# plt.savefig('/home/ariel/Workspace/GASP/HST/Analysis/Plots/Paper/morphology_trends.png', dpi=300)


fig = plt.figure(figsize=(6.5, 5.5))

scatter = plt.scatter(f606w_input['axial_ratio'], output_f606w['stellar_mass'], c=1e3 * output_f606w['mwage'],
                      edgecolors='w', vmin=0, vmax=150, cmap='viridis',
                      s=1200*f606w_input['r_core_corr']/np.max(f606w_input['r_core_corr']))
plot_tools.plot_median_in_bins(f606w_input['axial_ratio'], output_f606w['stellar_mass'], nbins=5,
                               percentile_style='errorbar', color='k', percentiles_color='k', point_edgewidths=1,
                               point_edgecolors='w', point_size=50)
ax = plt.gca()
inset_ax = ax.inset_axes([0, 0, 0.3, 1])
inset_ax.set_xticklabels([])
inset_ax.set_yticklabels([])
inset_ax.spines["top"].set_visible(False)
inset_ax.spines["right"].set_visible(False)
inset_ax.patch.set_alpha(0)
inset_ax.tick_params(axis='both', length=0)
inset_ax.set_ylim(3.5, 8)

inset_ax.hist(output_f606w_unresolved['stellar_mass'], orientation='horizontal',
              bins=np.arange(3.5, 8.5, 0.25), histtype='step', color='k', density=True)
# ax.set_ylim(-2.05, 0.32)


plt.xlabel(r'Axial Ratio', fontsize=20)
plt.ylabel(r'$\log\,M/M_\odot$', fontsize=20)

plt.ylim(3.5, 8.5)

cb_ax = plt.gca().inset_axes([0.6, 0.9, 0.33, 0.04])
cb = plt.colorbar(cax=cb_ax, orientation='horizontal', mappable=scatter)
cb.set_label(r'$\langle\,t_\star\,\rangle_M\,\mathrm{[Myr]}$', fontsize=14)

plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

sns.despine()

fig.tight_layout()

plt.savefig('/home/ariel/Workspace/GASP/HST/Analysis/Plots/Paper/mass_axratio_color.png', dpi=300)




