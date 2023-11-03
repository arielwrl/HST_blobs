"""

ariel@padova

Plot results from the spectral synthesis

"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from astropy.table import Table
from toolbox.wololo import arcsectokpc

sns.set_style('ticks')

testid = 'default'

results = Table.read('/home/ariel/Workspace/GASP/HST/Data/tail_all_' + testid + '_halpha_bagpipes_results.fits')
input = Table.read('/home/ariel/Workspace/GASP/HST/Data/tail_all_halpha_bagpipes_input.fits')

# results_uv = Table.read('/home/ariel/Workspace/GASP/HST/Data/tail_allblobs_' + testid + '_f275w_bagpipes_results.fits')
# input_uv = Table.read('/home/ariel/Workspace/GASP/HST/Data/tail_allblobs_f275w_bagpipes_input.fits')

# results_f606 = Table.read('/home/ariel/Workspace/GASP/HST/Data/tail_all_' + testid + '_f606w_bagpipes_results.fits')
# input_f606 = Table.read('/home/ariel/Workspace/GASP/HST/Data/tail_all_f606w_bagpipes_input.fits')

results['sfr'] = np.log10(results['sfr'])
# results_uv['sfr'] = np.log10(results_uv['sfr'])
# results_f606['sfr'] = np.log10(results_f606['sfr'])
results['galaxy'] = results['galaxy'].astype(str)
# results_uv['galaxy'] = results_uv['galaxy'].astype(str)
# results_f606['galaxy'] = results_f606['galaxy'].astype(str)
input = input[results['age'] < 0.5]
# input_uv = input_uv[results_uv['age'] < 0.5]
# input_f606 = input_f606[results_f606['age'] < 1.5]
results = results[results['age'] < 0.5]
# results_uv = results_uv[results_uv['age'] < 0.5]
# results_f606 = results_f606[results_f606['age'] < 1.5]
# results['mwage'] = results['mwage']*1000
# results_uv['mwage'] = results_uv['mwage']*1000
# results['age'] = results['age']*1000
# results_uv['age'] = results_uv['age']*1000
# results['Av'] = results['Av'] * results['eta']
# results_uv['Av'] = results_uv['Av'] * results_uv['eta']

df = results['blob_id', 'galaxy', 'age', 'tau', 'metallicity', 'Av', 'eta', 'sfr', 'stellar_mass',
             'formed_mass', 'mwage', 'tform'].to_pandas()
# df_uv = results_uv['blob_id', 'galaxy', 'age', 'tau', 'metallicity', 'Av', 'eta', 'sfr', 'stellar_mass',
#                    'formed_mass', 'mwage', 'tform'].to_pandas()
# df_f606 = results_f606['blob_id', 'galaxy', 'age', 'tau', 'metallicity', 'Av', 'eta', 'sfr', 'stellar_mass',
#                      'formed_mass', 'mwage', 'tform'].to_pandas()

fig = plt.figure()
ax = fig.add_subplot(111)
sns.boxplot(x='galaxy', y='mwage', data=df, whis=[2, 98], width=.6, palette='flare', fliersize=0, ax=ax)
ax.set_xlabel(' ')
ax.set_ylabel(r'$\langle t \rangle_M \,\mathrm{[Myr]}$', fontsize=20)
ax.tick_params(axis='x', labelsize=13, rotation=45, direction='out')
sns.despine()
fig.tight_layout()
# plt.ylim(0, 0.08)
fig.savefig('/home/ariel/Workspace/GASP/HST/BAGPIPES/plots/results/box_halpha_mwage.png')
plt.close()

# fig = plt.figure()
# ax = fig.add_subplot(111)
# sns.boxplot(x='galaxy', y='mwage', data=df_uv, whis=[2, 98], width=.6, palette='crest', fliersize=0, ax=ax)
# ax.set_xlabel(' ')
# ax.set_ylabel(r'$\langle t \rangle_M \,\mathrm{[Myr]}$', fontsize=20)
# ax.tick_params(axis='x', labelsize=12, rotation=45, direction='out')
# sns.despine()
# fig.tight_layout()
# # plt.ylim(0, 0.5)
# fig.savefig('/home/ariel/Workspace/GASP/HST/BAGPIPES/plots/results/box_uv_mwage.png')
# plt.close()

# fig = plt.figure()
# ax = fig.add_subplot(111)
# sns.boxplot(x='galaxy', y='mwage', data=df_f606, whis=[2, 98], width=.6, palette='crest', fliersize=0, ax=ax)
# ax.set_xlabel(' ')
# ax.set_ylabel(r'$\langle t \rangle_M \,\mathrm{[Myr]}$', fontsize=20)
# ax.tick_params(axis='x', labelsize=12, rotation=45, direction='out')
# sns.despine()
# fig.tight_layout()
# # plt.ylim(0, 0.5)
# fig.savefig('/home/ariel/Workspace/GASP/HST/BAGPIPES/plots/results/box_f606_mwage.png')
# plt.close()

fig = plt.figure()
ax = fig.add_subplot(111)
sns.boxplot(x='galaxy', y='age', data=df, whis=[2, 98], width=.6, palette='flare', fliersize=0, ax=ax)
ax.set_xlabel(' ')
ax.set_ylabel(r'$t_0 \,\mathrm{[Myr]}$', fontsize=20)
ax.tick_params(axis='x', labelsize=13, rotation=45, direction='out')
sns.despine()
fig.tight_layout()
# plt.ylim(0, 0.08)
fig.savefig('/home/ariel/Workspace/GASP/HST/BAGPIPES/plots/results/box_halpha_age.png')
plt.close()

# fig = plt.figure()
# ax = fig.add_subplot(111)
# sns.boxplot(x='galaxy', y='age', data=df_uv, whis=[2, 98], width=.6, palette='crest', fliersize=0, ax=ax)
# ax.set_xlabel(' ')
# ax.set_ylabel(r'$t_0 \,\mathrm{[Myr]}$', fontsize=20)
# ax.tick_params(axis='x', labelsize=12, rotation=45, direction='out')
# sns.despine()
# fig.tight_layout()
# # plt.ylim(0, 0.5)
# fig.savefig('/home/ariel/Workspace/GASP/HST/BAGPIPES/plots/results/box_uv_age.png')
# plt.close()

# fig = plt.figure()
# ax = fig.add_subplot(111)
# sns.boxplot(x='galaxy', y='age', data=df_f606, whis=[2, 98], width=.6, palette='crest', fliersize=0, ax=ax)
# ax.set_xlabel(' ')
# ax.set_ylabel(r'$t_0 \,\mathrm{[Myr]}$', fontsize=20)
# ax.tick_params(axis='x', labelsize=12, rotation=45, direction='out')
# sns.despine()
# fig.tight_layout()
# # plt.ylim(0, 0.5)
# fig.savefig('/home/ariel/Workspace/GASP/HST/BAGPIPES/plots/results/box_f606_age.png')
# plt.close()

fig = plt.figure()
ax = fig.add_subplot(111)
sns.boxplot(x='galaxy', y='stellar_mass', data=df, whis=[2, 98], width=.6, palette='flare', fliersize=0, ax=ax)
ax.set_xlabel(' ')
ax.set_ylabel(r'$\log\,M/M_\odot$', fontsize=20)
ax.tick_params(axis='x', labelsize=13, rotation=45, direction='out')
sns.despine()
fig.tight_layout()
# plt.ylim(0, 0.08)
fig.savefig('/home/ariel/Workspace/GASP/HST/BAGPIPES/plots/results/box_halpha_mass.png')
plt.close()

# fig = plt.figure()
# ax = fig.add_subplot(111)
# sns.boxplot(x='galaxy', y='stellar_mass', data=df_uv, whis=[2, 98], width=.6, palette='crest', fliersize=0, ax=ax)
# ax.set_xlabel(' ')
# ax.set_ylabel(r'$\log\,M/M_\odot$', fontsize=20)
# ax.tick_params(axis='x', labelsize=12, rotation=45, direction='out')
# sns.despine()
# fig.tight_layout()
# # plt.ylim(0, 0.5)
# fig.savefig('/home/ariel/Workspace/GASP/HST/BAGPIPES/plots/results/box_uv_mass.png')
# plt.close()

# fig = plt.figure()
# ax = fig.add_subplot(111)
# sns.boxplot(x='galaxy', y='stellar_mass', data=df_f606, whis=[2, 98], width=.6, palette='crest', fliersize=0, ax=ax)
# ax.set_xlabel(' ')
# ax.set_ylabel(r'$\log\,M/M_\odot$', fontsize=20)
# ax.tick_params(axis='x', labelsize=12, rotation=45, direction='out')
# sns.despine()
# fig.tight_layout()
# # plt.ylim(0, 0.5)
# fig.savefig('/home/ariel/Workspace/GASP/HST/BAGPIPES/plots/results/box_f606_mass.png')
# plt.close()

fig = plt.figure()
ax = fig.add_subplot(111)
sns.boxplot(x='galaxy', y='eta', data=df, whis=[2, 98], width=.6, palette='flare', fliersize=0, ax=ax)
ax.set_xlabel(' ')
ax.set_ylabel(r'$\eta$', fontsize=20)
ax.tick_params(axis='x', labelsize=13, rotation=45, direction='out')
sns.despine()
fig.tight_layout()
# plt.ylim(0, 0.08)
fig.savefig('/home/ariel/Workspace/GASP/HST/BAGPIPES/plots/results/box_halpha_eta.png')
plt.close()

# fig = plt.figure()
# ax = fig.add_subplot(111)
# sns.boxplot(x='galaxy', y='eta', data=df_uv, whis=[2, 98], width=.6, palette='crest', fliersize=0, ax=ax)
# ax.set_xlabel(' ')
# ax.set_ylabel(r'$\eta$', fontsize=20)
# ax.tick_params(axis='x', labelsize=12, rotation=45, direction='out')
# sns.despine()
# fig.tight_layout()
# # plt.ylim(0, 0.5)
# fig.savefig('/home/ariel/Workspace/GASP/HST/BAGPIPES/plots/results/box_uv_eta.png')
# plt.close()

# fig = plt.figure()
# ax = fig.add_subplot(111)
# sns.boxplot(x='galaxy', y='eta', data=df_f606, whis=[2, 98], width=.6, palette='crest', fliersize=0, ax=ax)
# ax.set_xlabel(' ')
# ax.set_ylabel(r'$\eta$', fontsize=20)
# ax.tick_params(axis='x', labelsize=12, rotation=45, direction='out')
# sns.despine()
# fig.tight_layout()
# # plt.ylim(0, 0.5)
# fig.savefig('/home/ariel/Workspace/GASP/HST/BAGPIPES/plots/results/box_f606_eta.png')
# plt.close()

fig = plt.figure()
ax = fig.add_subplot(111)
sns.boxplot(x='galaxy', y='Av', data=df, whis=[2, 98], width=.6, palette='flare', fliersize=0, ax=ax)
ax.set_xlabel(' ')
ax.set_ylabel(r'$A_V$', fontsize=20)
ax.tick_params(axis='x', labelsize=13, rotation=45, direction='out')
sns.despine()
fig.tight_layout()
# plt.ylim(0, 0.08)
fig.savefig('/home/ariel/Workspace/GASP/HST/BAGPIPES/plots/results/box_halpha_Av.png')
plt.close()

# fig = plt.figure()
# ax = fig.add_subplot(111)
# sns.boxplot(x='galaxy', y='Av', data=df_uv, whis=[2, 98], width=.6, palette='crest', fliersize=0, ax=ax)
# ax.set_xlabel(' ')
# ax.set_ylabel(r'$A_V$', fontsize=20)
# ax.tick_params(axis='x', labelsize=12, rotation=45, direction='out')
# sns.despine()
# fig.tight_layout()
# # plt.ylim(0, 0.5)
# fig.savefig('/home/ariel/Workspace/GASP/HST/BAGPIPES/plots/results/box_uv_Av.png')
# plt.close()

# fig = plt.figure()
# ax = fig.add_subplot(111)
# sns.boxplot(x='galaxy', y='Av', data=df_f606, whis=[2, 98], width=.6, palette='crest', fliersize=0, ax=ax)
# ax.set_xlabel(' ')
# ax.set_ylabel(r'$A_V$', fontsize=20)
# ax.tick_params(axis='x', labelsize=12, rotation=45, direction='out')
# sns.despine()
# fig.tight_layout()
# # plt.ylim(0, 0.5)
# fig.savefig('/home/ariel/Workspace/GASP/HST/BAGPIPES/plots/results/box_f606_Av.png')
# plt.close()

fig = plt.figure()
ax = fig.add_subplot(111)
sns.boxplot(x='galaxy', y='sfr', data=df, whis=[2, 98], width=.6, palette='flare', fliersize=0, ax=ax)
ax.set_xlabel(' ')
ax.set_ylabel(r'$\log \mathrm{SFR\,[M_\odot/yr]} $', fontsize=20)
ax.tick_params(axis='x', labelsize=13, rotation=45, direction='out')
sns.despine()
fig.tight_layout()
plt.ylim(-5, 0.0)
fig.savefig('/home/ariel/Workspace/GASP/HST/BAGPIPES/plots/results/box_halpha_sfr.png')
plt.close()

# fig = plt.figure()
# ax = fig.add_subplot(111)
# sns.boxplot(x='galaxy', y='sfr', data=df_uv, whis=[2, 98], width=.6, palette='crest', fliersize=0, ax=ax)
# ax.set_xlabel(' ')
# ax.set_ylabel(r'$\log \mathrm{SFR\,[M_\odot/yr]} $', fontsize=20)
# ax.tick_params(axis='x', labelsize=12, rotation=45, direction='out')
# sns.despine()
# fig.tight_layout()
# plt.ylim(-5, 0.0)
# fig.savefig('/home/ariel/Workspace/GASP/HST/BAGPIPES/plots/results/box_uv_sfr.png')
# plt.close()

# fig = plt.figure()
# ax = fig.add_subplot(111)
# sns.boxplot(x='galaxy', y='sfr', data=df_f606, whis=[2, 98], width=.6, palette='crest', fliersize=0, ax=ax)
# ax.set_xlabel(' ')
# ax.set_ylabel(r'$\log \mathrm{SFR\,[M_\odot/yr]} $', fontsize=20)
# ax.tick_params(axis='x', labelsize=12, rotation=45, direction='out')
# sns.despine()
# fig.tight_layout()
# plt.ylim(-5, 0.0)
# fig.savefig('/home/ariel/Workspace/GASP/HST/BAGPIPES/plots/results/box_f606_sfr.png')
# plt.close()






# plt.setp(plt.gca().get_xticklabels(), rotation=45, ha="right",
#          rotation_mode="anchor")

# galaxies = ['JO175', 'JO201', 'JO204', 'JO206', 'JW39', 'JW100']
#
# # results['sfr_average'] = results['formed_mass']/(results['age']*1e9)
#
# histogram_properties = ['age', 'tau', 'stellar_mass', 'Av', 'metallicity', 'sfr']
# histogram_labels = [r'$t^\star$', r'$\tau$', r'$\log\,M/M_\odot$', r'$A_V$', r'$Z/Z_\odot$', r'$SFR_{20\mathrm{Myr}}$']
#
# results['sfr'] = np.log10(results['sfr'])
# results_uv['sfr'] = np.log10(results_uv['sfr'])

# For H alpha selected blobs:

# for galaxy in galaxies:
#
#     galaxy_table = results[(results['galaxy'] == galaxy) & (results['age'] < 0.3)]
#
#     fig = plt.figure(figsize=(10, 7))
#
#     fig.suptitle(galaxy, fontsize=20)
#
#     for i in range(len(histogram_properties)):
#
#         ax = fig.add_subplot(2, 3, i+1)
#
#         hist_range = [np.min(galaxy_table[histogram_properties[i]].data),
#                       np.max(galaxy_table[histogram_properties[i]].data)]
#
#         sns.distplot(galaxy_table[histogram_properties[i]].data, bins=10, kde=False, hist_kws={'range': hist_range})
#
#         ax.set_xlabel(histogram_labels[i], fontsize=16)
#         ax.set_ylabel('Density', fontsize=16)
#
#         sns.despine()
#
#     fig.tight_layout()
#
#     plt.savefig('/home/ariel/Workspace/GASP/HST/BAGPIPES/plots/results/' + galaxy + testid +
#                 '_halpha_distributions.png', dpi=200)
#
# plt.close('all')
#
#
# # Now the same for UV:
#
# for galaxy in galaxies:
#
#     galaxy_table = results_uv[(results_uv['galaxy'] == galaxy) & (results_uv['age'] < 0.6)]
#
#     fig = plt.figure(figsize=(10, 7))
#
#     fig.suptitle(galaxy, fontsize=20)
#
#     for i in range(len(histogram_properties)):
#
#         ax = fig.add_subplot(2, 3, i+1)
#
#         hist_range = [np.min(galaxy_table[histogram_properties[i]].data),
#                       np.max(galaxy_table[histogram_properties[i]].data)]
#
#         sns.distplot(galaxy_table[histogram_properties[i]].data, bins=10, kde=False, hist_kws={'range': hist_range})
#
#         ax.set_xlabel(histogram_labels[i], fontsize=16)
#         ax.set_ylabel('Density', fontsize=16)
#
#         sns.despine()
#
#     fig.tight_layout()
#
#     plt.savefig('/home/ariel/Workspace/GASP/HST/BAGPIPES/plots/results/' + galaxy + testid +
#                 '_uv_distributions.png', dpi=200)
#
# plt.close('all')