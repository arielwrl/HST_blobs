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

ax[0].hist((results['age'] - results['age_young'])/results['age'], alpha=0.5)
ax[0].set_xlabel(r'$(t^\star_\mathrm{single}-t^\star_\mathrm{double})/t^\star_\mathrm{single}$')

ax[1].hist(results['Av_1']*results['eta_1'] - results['Av_2']*results['eta_2'], alpha=0.5)
ax[1].set_xlabel(r'$(A_V_\mathrm{single}-A_V_\mathrm{double}/A_V_\mathrm{single}$')

ax[2].hist(results['formed_mass']- results['formed_mass_young'], alpha=0.5)
ax[2].set_xlabel(r'$(\log\;M^\star/M_\odot_\mathrm{single}-\log\;M^\star/M_\odot_\mathrm{double}/\log\;M^\star'
                 r'/M_\odot_\mathrm{single}$')

sns.despine()

plt.show()

