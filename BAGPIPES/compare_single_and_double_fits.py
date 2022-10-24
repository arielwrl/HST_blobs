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

results = Table.read('/home/ariel/Workspace/GASP/HST/Data/extraplanar_matched.fits')

fig, ax = plt.subplots(2, 2)

ax[0, 0].scatter(results['age'], results['age_young'], alpha=0.5)
ax[0, 0].set_xlim(0, 0.3)
ax[0, 0].set_ylim(0, 0.3)
ax[0, 0].set_xlabel(r'Age (Single)')
ax[0, 0].set_ylabel(r'Age (Double)')

ax[0, 1].scatter(results['Av_1']*results['eta_1'], results['Av_2']*results['eta_2'], alpha=0.5)
ax[0, 1].set_xlim(0, 1)
ax[0, 1].set_ylim(0, 1)
ax[0, 1].set_xlabel(r'$A_V$ (Single)')
ax[0, 1].set_ylabel(r'$A_V$ (Double)')

ax[1, 0].scatter(results['formed_mass'], results['formed_mass_young'], alpha=0.5)
ax[1, 0].set_xlim(4, 7.5)
ax[1, 0].set_ylim(4, 7.5)
ax[1, 0].set_xlabel(r'Mass Formed (Single)')
ax[1, 0].set_ylabel(r'Mass Formed (Double)')

ax[1, 1].scatter(results['tau'], results['tau_young'], alpha=0.5)
ax[1, 1].set_xlabel(r'$\tau$ (Single)')
ax[1, 1].set_ylabel(r'$\tau$ (Double)')


sns.despine()

plt.show()

