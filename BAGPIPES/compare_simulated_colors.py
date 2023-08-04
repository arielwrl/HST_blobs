"""

ariel@viapalermo
10/02/2023

Apparently metallicity increases with distance because clumps get bluer in the UV, let's look into it using
simulations and stuff

- Make color-age for stellar-only, and for models with nebular component but no nebular continuum
- Make color-color with Ha flux, and look at nebular continuum somehow
- Check disk metallicities
- Compare parameters of fits that had Z>1.5 with metallicity limited

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from astropy.table import Table
import seaborn as sns
from toolbox import plot_tools
from matplotlib import cm

sns.set_style('ticks')

f275w_input = Table.read('/home/ariel/Workspace/GASP/HST/Data/f275w_bagpipes_input.fits')
output_f275w = Table.read('/home/ariel/Workspace/GASP/HST/Data/f275w_dexp_logprior_single_bagpipes_results.fits')
output_f275w = output_f275w[(~f275w_input['disk']) & ((f275w_input['level'] == 0) | (f275w_input['leaf_flag'] == 1))]
simulations = Table.read('/home/ariel/Workspace/GASP/HST/Data/simulated_spectra/simulations_table.fits')

# plt.scatter(-2.5*np.log10(output_f275w['photometry_obs'][:, 0]/output_f275w['photometry_obs'][:, 2]),
#             -2.5*np.log10(output_f275w['photometry_obs'][:, 3]/output_f275w['photometry_obs'][:, 2]),
#             edgecolors='w', color='darkviolet', alpha=0.8, s=10)
#
# plt.scatter(-2.5*np.log10(simulations['photometry_syn'][:, 0]/simulations['photometry_syn'][:, 2]),
#             -2.5*np.log10(simulations['photometry_syn'][:, 3]/simulations['photometry_syn'][:, 2]),
#             c=simulations['metallicity'])
# plt.colorbar()


plt.hexbin(-2.5*np.log10(output_f275w['photometry_syn'][:, 0]/output_f275w['photometry_syn'][:, 2]),
            -2.5*np.log10(output_f275w['photometry_syn'][:, 3]/output_f275w['photometry_syn'][:, 2]),
            C=output_f275w['metallicity'], gridsize=12, cmap='viridis_r', vmin=0.9)

plt.colorbar(label=r'$Z/Z_\odot$')

model_flag = (simulations['Av'] == 0.01) & (simulations['logU'] == -2.5)
plt.scatter(-2.5*np.log10(simulations['photometry_syn'][:, 0]/simulations['photometry_syn'][:, 2]),
            -2.5*np.log10(simulations['photometry_syn'][:, 3]/simulations['photometry_syn'][:, 2]),
            edgecolors='k', c=simulations['metallicity'], s=30, cmap='Greens')
plt.colorbar(label=r'$Z/Z_\odot$')

# plt.scatter(-2.5*np.log10(output_f275w['photometry_obs'][:, 0]/output_f275w['photometry_obs'][:, 2]),
#             -2.5*np.log10(output_f275w['photometry_obs'][:, 3]/output_f275w['photometry_obs'][:, 2]),
#             edgecolors='w', color='k', alpha=0.8, s=10)

plt.xlabel('F275W-F606W')
plt.ylabel('F680N-F606W')
#
# plt.figure()
# model_flag = (simulations['Av'] == 0.01) & (simulations['tau'] == 0.001) & (simulations['logU'] == -2.5)
# scatter = plt.scatter(simulations['age'][model_flag], -2.5*np.log10(simulations['photometry_syn'][:, 0]/simulations['photometry_syn'][:, 2])[model_flag],
#                       c=simulations['metallicity'][model_flag], cmap='Spectral_r')
# #
# metallicities = np.unique(simulations['metallicity'])
#
# colors = [cm.Spectral_r(x) for x in np.linspace(0, 1, len(metallicities))]
# colordict = dict(zip(metallicities, colors))
#
# for metallicity in metallicities:
#     model_flag = (simulations['logU'] == -2.5) & (simulations['Av'] == 0.01) & (simulations['tau'] == 0.001) &\
#                  (simulations['metallicity'] == metallicity)
#     plt.plot(simulations['age'][model_flag], -2.5*np.log10(simulations['photometry_syn'][:, 0]/simulations['photometry_syn'][:, 2])[model_flag],
#              color=colordict[metallicity])

# plt.colorbar(label=r'$Z/Z_\odot$')
# plt.xlabel('SSP Age')
# plt.ylabel('F275W-F606W')


#
# plt.figure()
#
# plt.scatter(-2.5*np.log10(output_f275w['photometry_syn'][:, 0]/output_f275w['photometry_syn'][:, 2]),
#             output_f275w['metallicity'], edgecolors='w', color='darkviolet', alpha=0.8, s=10)
#
# model_flag = (simulations['Av'] == 0.01) & (simulations['logU'] == -2.5)
# plt.scatter(-2.5*np.log10(simulations['photometry_syn'][:, 0]/simulations['photometry_syn'][:, 2])[model_flag],
#             simulations['metallicity'][model_flag], c=simulations['age'][model_flag], cmap='Spectral_r')
#
# plt.xlabel('F275W-F606W')
# plt.ylabel(r'$Z/Z_\odot$')
#
# plt.figure()
#
# plt.scatter(-2.5*np.log10(output_f275w['photometry_obs'][:, 0]/output_f275w['photometry_obs'][:, 2]),
#             output_f275w['metallicity'], edgecolors='w', color='darkviolet', alpha=0.8, s=10)
#
# model_flag = (simulations['Av'] == 0.01) & (simulations['logU'] == -2.5)
# plt.scatter(-2.5*np.log10(simulations['photometry_syn'][:, 0]/simulations['photometry_syn'][:, 2])[model_flag],
#             simulations['metallicity'][model_flag], c=simulations['age'][model_flag], cmap='Spectral_r')
#
# plt.xlabel('F275W-F606W')
# plt.ylabel(r'$Z/Z_\odot$')
#
