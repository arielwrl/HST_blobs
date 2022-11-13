"""

ariel@oapd
November 8, 2022

Compare Halpha tail clumps fitted with ssps or dexp

"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack
import seaborn as sns
from toolbox import plot_tools
from astropy.io import fits

plt.ion()

sns.set_style('ticks')

input = Table.read('/home/ariel/Workspace/GASP/HST/Data/tail_f275w_bagpipes_input.fits')

output = Table.read('/home/ariel/Workspace/GASP/HST/Data/tail_f275w_parametric_bagpipes_results.fits')
output_ssp = Table.read('/home/ariel/Workspace/GASP/HST/Data/tail_f275w_ssp_bagpipes_results.fits')
output_constant = Table.read('/home/ariel/Workspace/GASP/HST/Data/tail_f275w_constant_bagpipes_results.fits')
output_smalltau = Table.read('/home/ariel/Workspace/GASP/HST/Data/tail_f275w_exp_smalltau_bagpipes_results.fits')
output_logprior = Table.read('/home/ariel/Workspace/GASP/HST/Data/tail_f275w_exp_logprior_bagpipes_results.fits')


flag = input['sel_flag'] == 31

input = input[flag]
output = output[flag]
output_ssp = output_ssp[flag]

phot_eq_obs = np.log10(output_ssp['photometry_obs'][:, 3]/output_ssp['photometry_obs'][:, 2])
phot_eq_def = np.log10(output['photometry_syn'][:, 3]/output['photometry_syn'][:, 2])
phot_eq_ssp = np.log10(output_ssp['photometry_syn'][:, 3]/output_ssp['photometry_syn'][:, 2])

f606n_diff = output['photometry_syn'][:, 3]/output['photometry_obs'][:, 3]
#
# fig = plt.figure(figsize=(8, 8))
#
# plt.scatter(x=output_logprior['mwage'].tolist(), y=output_ssp['mwage'].tolist())
# x = np.linspace(0, 0.25)
# y = x
# plt.plot(x, y, '--k')
# plt.xlim(0, 0.25)
# plt.ylim(0, 0.25)
#
# fig.subplots_adjust(left=0.1, bottom=0.1)


fig = plt.figure(figsize=(10, 7))

plt.hexbin(x=output_ssp['mwage'], y=output_logprior['mwage'], C=(output_logprior['chi2']/5-output_ssp['chi2']/5)
           , vmin=-0.2, vmax=0.2, gridsize=60, cmap='viridis', edgecolors='w')
x = np.linspace(0, 0.25)
y = x
plt.plot(x, y, '--k')

plt.ylabel(r'$\langle\,t_\star\,\rangle_M\,\mathrm{[Gyr]}$ (Exponential)', fontsize=20)
plt.xlabel(r'$\langle\,t_\star\,\rangle_M\,\mathrm{[Gyr]}$ (SSP)', fontsize=20)
plt.xlim(0, 0.25)
plt.ylim(0, 0.25)

cb = plt.colorbar()
cb.set_label(r'$\Delta\chi^2$', fontsize=20)

hist_palette = sns.color_palette('viridis', 5)
ax_inset = plt.gca().inset_axes([0.6, 0.1, 0.35, 0.3])
ax_inset.hist(output_logprior['chi2']/5-output_ssp['chi2']/5, range=[-0.5, 0.5], bins=21, color=hist_palette[2])
ax_inset.set_xlabel(r'$\Delta\chi^2$', fontsize=15)
ax_inset.set_xlim(-0.5, 0.5)

fig.subplots_adjust(right=0.99, left=0.08, top=0.97, bottom=0.1)

plt.savefig('exp_ssp_logprior.png')

#
# sns.scatterplot(x=np.log10(input['F606W']).tolist(),
#                 y=np.log10(input['F680N']/input['F606W']).tolist())
#
# sns.scatterplot(x=np.log10(output['photometry_syn'][:, 2]).tolist(),
#                 y=np.log10(output['photometry_syn'][:, 3]/output['photometry_syn'][:, 2]).tolist())
#
# sns.scatterplot(x=np.log10(output_ssp['photometry_syn'][:, 2]).tolist(),
#                 y=np.log10(output_ssp['photometry_syn'][:, 3]/output_ssp['photometry_syn'][:, 2]).tolist())
#
# plt.figure()
#
# phot_eq_obs = np.log10(output_ssp['photometry_obs'][:, 0]/output_ssp['photometry_obs'][:, 1])
# phot_eq_def = np.log10(output['photometry_syn'][:, 0]/output['photometry_syn'][:, 1])
# phot_eq_ssp = np.log10(output_ssp['photometry_syn'][:, 0]/output_ssp['photometry_syn'][:, 1])
#
# sns.scatterplot(x=(phot_eq_obs-phot_eq_ssp).tolist(), y=(phot_eq_obs-phot_eq_def).tolist())
# x = np.linspace(0, 0.6)
# y = x
# plt.plot(x, y, '--k')