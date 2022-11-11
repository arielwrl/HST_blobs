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

flag = input['sel_flag'] == 31

input = input[flag]
output = output[flag]
output_ssp = output_ssp[flag]
output_constant = output_constant[flag]

phot_eq_obs = np.log10(output_ssp['photometry_obs'][:, 3]/output_ssp['photometry_obs'][:, 2])
phot_eq_def = np.log10(output['photometry_syn'][:, 3]/output['photometry_syn'][:, 2])
phot_eq_ssp = np.log10(output_ssp['photometry_syn'][:, 3]/output_ssp['photometry_syn'][:, 2])

f606n_diff = output['photometry_syn'][:, 3]/output['photometry_obs'][:, 3]

fig = plt.figure(figsize=(8, 8))

plt.scatter(x=output_smalltau['mwage'].tolist(), y=output_constant['mwage'].tolist())
x = np.linspace(0, 0.25)
y = x
plt.plot(x, y, '--k')
plt.xlim(0, 0.25)
plt.ylim(0, 0.25)

# plt.xlabel(r'$\langle\,t_\star\,\rangle_M\,\mathrm{[Gyr]}$ (Delayed)', fontsize=20)
# plt.ylabel(r'$\langle\,t_\star\,\rangle_M\,\mathrm{[Gyr]}$ (SSP)', fontsize=20)

fig.subplots_adjust(left=0.1, bottom=0.1)


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