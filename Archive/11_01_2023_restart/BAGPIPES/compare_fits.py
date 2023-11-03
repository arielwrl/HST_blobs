"""

ariel@oapd
06/04/22

Compare results for different bagpipes configurations

"""

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from astropy.table import Table

sns.set_style('ticks')
palette = sns.color_palette('plasma', 3)

input_606 = Table.read('/home/ariel/Workspace/GASP/HST/Data/tail_all_f606w_bagpipes_input.fits')
input_all = Table.read('/home/ariel/Workspace/GASP/HST/Data/tail_all_halpha_bagpipes_input.fits')
default = Table.read('/home/ariel/Workspace/GASP/HST/Data/tail_all31_default_halpha_bagpipes_results.fits')
default_all = Table.read('/home/ariel/Workspace/GASP/HST/Data/tail_all_default_halpha_bagpipes_results.fits')
default_hacked_all = Table.read('/home/ariel/Workspace/GASP/HST/Data/tail_all_default_nebcont_hacked_halpha_bagpipes_'
                                'results.fits')
default_hacked = default_hacked_all[(input_all['sel_flag'] == 31)]
nonparametric = Table.read('/home/ariel/Workspace/GASP/HST/Data/tail_all31_default_nonparametric_halpha_bagpipes'
                           '_results.fits')


# fig = plt.figure()
# plt.hist(np.absolute(default['chi2']/(5-7)), bins=25, color=palette[0], histtype='step', density=True, label='Default',
#          range=[0, 8])
# plt.hist(np.absolute(default_hacked['chi2']/(5-7)), bins=25, color=palette[1], histtype='step', density=True,
#          label='No Nebular Continuum', range=[0, 8])
# plt.hist(np.absolute(nonparametric['chi2']/(5-17)), bins=25, color=palette[2], histtype='step', density=True,
#          label='Nonparametric', range=[0, 8])
# plt.legend(frameon=False, fontsize=15)
# plt.xlabel(r'$\chi^2/(n-m)$', fontsize=20)
# plt.ylabel('Density', fontsize=20)
# sns.despine()
# fig.tight_layout()
# plt.show()
#
# fig = plt.figure()
# plt.hist(default['mwage'], bins=18, color=palette[0], histtype='step', density=True, label='Default', range=[0, 1.1])
# plt.hist(default_hacked['mwage'], bins=18, color=palette[1], histtype='step', density=True, label='No Nebular Continuum', range=[0, 1.1])
# plt.hist(nonparametric['mwage'], bins=18, color=palette[2], histtype='step', density=True, label='Nonparametric', range=[0, 1.1])
# plt.legend(frameon=False, fontsize=15)
# plt.xlabel(r'$\langle\;t^\star\rangle_M$', fontsize=20)
# plt.ylabel('Density', fontsize=20)
# sns.despine()
# fig.tight_layout()
# plt.show()

# fig = plt.figure()
# plt.hist(np.log10(input_606['F680N_line_flux']), bins=18, color=palette[0], histtype='step', density=True, label='606')
# plt.hist(np.log10(input_all['F680N_line_flux']), bins=18, color=palette[1], histtype='step', density=True, label='Halpha')
# plt.legend(frameon=False, fontsize=15)
# plt.xlabel(r'$H\alpha$', fontsize=20)
# plt.ylabel('Density', fontsize=20)
# sns.despine()
# fig.tight_layout()
# plt.show()
