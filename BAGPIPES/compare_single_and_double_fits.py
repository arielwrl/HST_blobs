"""

ariel@oapd
04/10/22

Compare results for different bagpipes configurations

"""

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from astropy.table import Table

sns.set_style('ticks')

results = Table.read('/home/ariel/Workspace/GASP/HST/Data/old/extraplanar_matched.fits')

fig, ax = plt.subplots(3, 1)

ax[0].hist((results['age'] - results['age_young'])/results['age'], alpha=0.5, range=[-1, 1], bins=15, density=True)
sns.kdeplot(((results['age'] - results['age_young'])/results['age']).tolist(), ax=ax[0])
# ax[0].set_xlabel(r'$(t^\star_\mathrm{single}-t^\star_\mathrm{double})/t^\star_\mathrm{single}$')
ax[0].set_xlim(-1, 1)

ax[1].hist((results['Av_1']*results['eta_1'] - results['Av_2']*results['eta_2'])/results['Av_1']*results['eta_1'],
           alpha=0.5, range=[-1, 1], bins=15, density=True)
sns.kdeplot(((results['Av_1']*results['eta_1'] - results['Av_2']*results['eta_2'])/results['Av_1']*results['eta_1']).tolist(), ax=ax[1])
# ax[1].set_xlabel(r'$(A_V_\mathrm{single}-A_V_\mathrm{double}/A_V_\mathrm{single}$')
ax[1].set_xlim(-1, 1)

ax[2].hist((results['formed_mass'] - results['formed_mass_young'])/results['formed_mass'], alpha=0.5,
           range=[-0.25, 0.25], bins=15, density=True)
sns.kdeplot(((results['formed_mass'] - results['formed_mass_young'])/results['formed_mass']).tolist(), ax=ax[2])
# ax[2].set_xlabel(r'$(\log\;M^\star/M_\odot_\mathrm{single}-\log\;M^\star/M_\odot_\mathrm{double}/\log\;M^\star'
#                  r'/M_\odot_\mathrm{single}$')
ax[2].set_xlim(-0.25, 0.25)

sns.despine()

plt.show()

