"""

ariel@oapd
07/03/2023

Compare results of vanilla detection with the ones with 1pix dilation

color: mediumvioletred

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from astropy.table import Table
import seaborn as sns
from toolbox import plot_tools

sns.set_style('ticks')

plt.ioff()

output = Table.read('/home/ariel/Workspace/GASP/HST/Data/bad_beta_match.fits')

f275w = output['photometry_obs_1'][:, 0]
f336w = output['photometry_obs_1'][:, 1]
beta = np.log10(f336w/f275w)/np.log10(3360/2750)

f275w_dilated = output['photometry_obs_2'][:, 0]
f336w_dilated = output['photometry_obs_2'][:, 1]
beta_dilated = np.log10(f336w_dilated/f275w_dilated)/np.log10(3360/2750)

bad_beta_flag = beta < -3.34

plt.figure()
plot_tools.plot_median_in_bins((beta-beta_dilated)[bad_beta_flag],
                               (100*(output['mwage_1']-output['mwage_2'])/output['mwage_1'])[bad_beta_flag],
                               color='red', percentile_style='errorbar',
                               percentiles_color='red')
plot_tools.plot_median_in_bins((beta-beta_dilated)[~bad_beta_flag],
                               (100*(output['mwage_1']-output['mwage_2'])/output['mwage_1'])[~bad_beta_flag],
                               color='blue', percentile_style='errorbar',
                               percentiles_color='blue')
plt.xlabel(r'$\Delta\beta$', fontsize=20)
plt.ylabel(r'$\Delta\;t_M$', fontsize=20)

plt.figure()
plot_tools.plot_median_in_bins(beta-beta_dilated, beta,
                               color='mediumvioletred', percentile_style='errorbar',
                               percentiles_color='mediumvioletred')
plt.xlabel(r'$\Delta\beta$', fontsize=20)
plt.ylabel(r'$\beta$', fontsize=20)


plt.figure()
plot_tools.plot_median_in_bins(1e3*output['mwage_1'],
                               (100*(output['mwage_1']-output['mwage_2'])/output['mwage_1']),
                               color='red', percentile_style='errorbar',
                               percentiles_color='red')
plt.xlabel(r'$t_M$', fontsize=20)
plt.ylabel(r'$\Delta\;t_M$', fontsize=20)


plt.figure()
plt.scatter(beta-beta_dilated, 100*(output['mwage_1']-output['mwage_2'])/output['mwage_1'], c=beta, vmax=-3.34, vmin=-7)
           # gridsize=20, mincnt=2)
plt.colorbar()


plt.figure()
plot_tools.plot_median_in_bins((beta-beta_dilated)[bad_beta_flag],
                               (output['stellar_mass_1']-output['stellar_mass_2'])[bad_beta_flag],
                               color='red', percentile_style='errorbar',
                               percentiles_color='red')
plot_tools.plot_median_in_bins((beta-beta_dilated)[~bad_beta_flag],
                               (output['stellar_mass_1']-output['stellar_mass_2'])[~bad_beta_flag],
                               color='blue', percentile_style='errorbar',
                               percentiles_color='blue')
plt.xlabel(r'$\Delta\beta$', fontsize=20)
plt.ylabel(r'$\Delta\;M$', fontsize=20)
