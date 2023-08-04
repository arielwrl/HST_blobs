"""

ariel@viapalermo
30/01/2023

Plot results of the clump matching

"""

import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_style('ticks')

data_dir = '/home/ariel/Workspace/GASP/HST/Data/'

detection = 'halpha'
config = 'dexp_logprior_single'

input_table = Table.read(data_dir + detection + '_bagpipes_input.fits')
output_table = Table.read(data_dir + detection + '_' + config + '_bagpipes_results.fits')

# tail_halpha = input_halpha['tail_gal_flag'] == 0
# extraplanar_halpha = input_halpha['tail_gal_flag'] == 1
# disk_halpha = input_halpha['tail_gal_flag'] == 2

output_table = output_table[output_table['match_flag_output']]

fig = plt.figure()

plt.errorbar(output_table['parent_mwage'], output_table['min_siblings_mwage'], xerr=output_table['parent_age_std'],
             yerr=output_table['min_siblings_age_std'], marker='o', elinewidth=0.25, fmt=' ', markeredgecolor='white')

x = np.linspace(0, 0.5)
y = x
plt.plot(x, y, '--k', zorder=10)

plt.xlim(0, 0.5)
plt.ylim(0, 0.5)

plt.xlabel('Parent MW Age', fontsize=20)
plt.ylabel('Minimum Siblings MW Age', fontsize=20)

fig.tight_layout()

plt.savefig('/home/ariel/Workspace/GASP/HST/Analysis/Plots/fatherandson' + '_' + config + '_' + detection + '.png')
plt.savefig('/home/ariel/Workspace/GASP/HST/Analysis/Plots/fatherandson' + '_' + config + '_' + detection + '.pdf')

fig = plt.figure()

plt.errorbar(output_table['parent_mwage'], output_table['parent_mwage'] - output_table['min_siblings_mwage'], xerr=output_table['parent_age_std'],
             yerr=output_table['min_siblings_age_std'], marker='o', elinewidth=0.25, fmt=' ', markeredgecolor='white')

x = np.linspace(0, 0.5)
y = x
plt.plot(x, y, '--k', zorder=10)

plt.xlim(0, 0.5)
plt.ylim(0, 0.5)

plt.xlabel('Parent MW Age', fontsize=14)
plt.ylabel('Parent MW Age - Minimum Siblings MW Age', fontsize=14)

fig.tight_layout()

plt.savefig('/home/ariel/Workspace/GASP/HST/Analysis/Plots/fatherandson' + '_' + config + '_' + detection + '_diff.png')
plt.savefig('/home/ariel/Workspace/GASP/HST/Analysis/Plots/fatherandson' + '_' + config + '_' + detection + '_diff.pdf')