
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack
import seaborn as sns
from toolbox import plot_tools
from astropy.io import fits
import pandas as pd

sns.set_style('ticks')
palette = sns.diverging_palette(220, 20, n=5)

input_uv = Table.read('/home/ariel/Workspace/GASP/HST/Data/f275w_bagpipes_input.fits')
output_uv = Table.read('/home/ariel/Workspace/GASP/HST/Data/f275w_dexp_logprior_single_bagpipes_results.fits')

input_halpha = Table.read('/home/ariel/Workspace/GASP/HST/Data/halpha_bagpipes_input.fits')
output_halpha = Table.read('/home/ariel/Workspace/GASP/HST/Data/halpha_dexp_logprior_single_bagpipes_results.fits')

input_f606w = Table.read('/home/ariel/Workspace/GASP/HST/Data/f606w_bagpipes_input.fits')
output_f606w = Table.read('/home/ariel/Workspace/GASP/HST/Data/f606w_dexp_logprior_single_bagpipes_results.fits')

input_halpha = input_halpha[(~input_halpha['disk']) & ((input_halpha['level'] == 0) | (input_halpha['leaf_flag'] == 1))]
input_uv = input_uv[(~input_uv['disk']) & ((input_uv['level'] == 0) | (input_uv['leaf_flag'] == 1))]

input_halpha = input_halpha[~output_halpha['bad_double_fit']]
output_halpha = output_halpha[~output_halpha['bad_double_fit']]

input_uv = input_uv[~output_uv['bad_double_fit']]
output_uv = output_uv[~output_uv['bad_double_fit']]

filter_list = ['F275W', 'F336W', 'F606W', 'F680N', 'F814W']

photometry_syn_uv = {}
photometry_syn_halpha = {}
photometry_syn_f606w = {}

photometry_syn_uv_iqr = {}
photometry_syn_halpha_iqr = {}
photometry_syn_f606w_iqr = {}

for i in range(5):
    photometry_syn_uv[filter_list[i]] = output_uv['photometry_syn'][:, i]
    photometry_syn_halpha[filter_list[i]] = output_halpha['photometry_syn'][:, i]
    photometry_syn_f606w[filter_list[i]] = output_f606w['photometry_syn'][:, i]
    photometry_syn_uv_iqr[filter_list[i]] = output_uv['photometry_syn_iqr'][:, i]
    photometry_syn_halpha_iqr[filter_list[i]] = output_halpha['photometry_syn_iqr'][:, i]
    photometry_syn_f606w_iqr[filter_list[i]] = output_f606w['photometry_syn_iqr'][:, i]

residuals_uv = []
residuals_halpha = []
residuals_f606w = []

for filter_name in filter_list:

    for i in range(len(input_uv)):

        residuals_uv.append([filter_name, (photometry_syn_uv[filter_name][i]-input_uv[filter_name][i]) /
                             (2 * input_uv['err' + filter_name][i] + photometry_syn_uv_iqr[filter_name][i])])

    for i in range(len(input_halpha)):

       residuals_halpha.append([filter_name, (photometry_syn_halpha[filter_name][i] - input_halpha[filter_name][i]) /
                                (2 * input_halpha['err' + filter_name][i] + photometry_syn_halpha_iqr[filter_name][i])])

    for i in range(len(input_f606w)):

       residuals_f606w.append([filter_name, (photometry_syn_f606w[filter_name][i] - input_f606w[filter_name][i]) /
                                (2 * input_f606w['err' + filter_name][i] + photometry_syn_f606w_iqr[filter_name][i])])


dataframe_uv = pd.DataFrame(residuals_uv, columns=['filters', 'residuals'])
dataframe_halpha = pd.DataFrame(residuals_halpha, columns=['filters', 'residuals'])
dataframe_f606w = pd.DataFrame(residuals_f606w, columns=['filters', 'residuals'])

fig, ax = plt.subplots(3, 1, figsize=(6, 8), sharex=True)

sns.boxplot(x='filters', y='residuals', data=dataframe_halpha, ax=ax[0], fliersize=0, whis=[15, 85], palette=palette)
sns.boxplot(x='filters', y='residuals', data=dataframe_uv, ax=ax[1], fliersize=0, whis=[15, 85], palette=palette)
sns.boxplot(x='filters', y='residuals', data=dataframe_f606w, ax=ax[2], fliersize=0, whis=[15, 85], palette=palette)

for i in range(3):
    ax[i].axhline(y=-1, color='k', linestyle='dotted', zorder=0)
    ax[i].axhline(y=1, color='k', linestyle='dotted', zorder=0)
    ax[i].axhline(y=0, color='k', linestyle='dashed', zorder=0)

ax[0].set_ylim(-1.95, 1.95)
ax[1].set_ylim(-1.95, 1.95)
ax[2].set_ylim(-1.95, 1.95)

ax[0].annotate(r'$H\alpha$ Selection', xy=(0.65, 0.9), xycoords='axes fraction', fontsize=16)
ax[1].annotate(r'F275W Selection', xy=(0.65, 0.9), xycoords='axes fraction', fontsize=16)
ax[2].annotate(r'F606W Selection', xy=(0.65, 0.9), xycoords='axes fraction', fontsize=16)

ax[1].set_ylabel(r'$\left(M_{l} - O_{l}\right)/(2\sigma_{l} + M_{IQR})$', fontsize=20)
ax[0].set_ylabel(' ')
ax[2].set_ylabel(' ')

ax[0].set_xlabel(' ')
ax[1].set_xlabel(' ')
ax[2].set_xlabel(' ')

ax[2].tick_params(axis='x', labelsize=20)

fig.tight_layout()
fig.subplots_adjust(hspace=0.0)

plt.savefig('/home/ariel/Workspace/GASP/HST/Analysis/Plots/Paper/uv_slope.png', dpi=200)

