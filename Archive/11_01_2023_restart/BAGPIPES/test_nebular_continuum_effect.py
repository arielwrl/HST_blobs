"""

ariel@sacrafamiglia
23/04/2022

Statistically testing the effect of the nebular continuum emission in the derived properties

"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from toolbox import plot_tools
from pysinopsis.utils import calc_mwage

data_dir = '/home/ariel/Workspace/GASP/HST/Data/'

parametric = Table.read(data_dir + 'tail_all31_default_halpha_bagpipes_results.fits')
nonparametric = Table.read(data_dir + 'tail_all31_default_nonparametric_halpha_bagpipes_results.fits')
nonparametric_hacked = Table.read(data_dir + 'tail_all31_default_nonparametric_halpha_bagpipes_results_nebcont_hacked'
                                             '.fits')

cases = (nonparametric['stellar_mass'] - parametric['stellar_mass'] > 1)
print(parametric['blob_id'][cases])

plt.figure()
plt.scatter(nonparametric['stellar_mass'], parametric['stellar_mass'])
plt.scatter(nonparametric['stellar_mass'][cases], parametric['stellar_mass'][cases])
x = np.linspace(np.min(nonparametric['stellar_mass']) - 1, np.max(nonparametric['stellar_mass']) + 1)
y = x
plt.plot(x, y, '--k')
plt.xlabel('Nonparametric Stellar Mass', fontsize=20)
plt.ylabel('Parametric Stellar Mass', fontsize=20)

age_bins = [0, 2.0e6, 4.0e6, 7.0e6, 2.0e7, 5.5e7, 2.0e8, 5.5e8, 1.0e9, 3.0e9,
            5.75e9, 1.0e10, 1.4e10]
age_bins_len = [age_bins[i+1]-age_bins[i] for i in range(12)]
age_bins_mid = [(age_bins[i+1]+age_bins[i])/2 for i in range(12)]

old_component = np.log10(np.array([np.sum(10**nonparametric['formed_mass'][i][9:])
                                   for i in range(len(nonparametric))]))

plt.figure()
plt.scatter(old_component, nonparametric['stellar_mass'] - parametric['stellar_mass'])
plt.scatter(old_component[cases], (nonparametric['stellar_mass'] - parametric['stellar_mass'])[cases])
plt.xlabel('Mass Formed in The Old Component', fontsize=20)
plt.ylabel('Current Stellar Mass Difference', fontsize=20)

# plt.scatter(nonparametric['mwage'], mwage_nonparametric)
# x = np.linspace(np.min(nonparametric['mwage']) - 1, np.max(nonparametric['mwage']) + 1)
# y = x
# plt.plot(x, y, '--k')

# plt.scatter(nonparametric['stellar_mass'], parametric['stellar_mass'])
# x = np.linspace(np.min(nonparametric['stellar_mass']) - 1, np.max(nonparametric['stellar_mass']) + 1)
# y = x
# plt.plot(x, y, '--k')

old_component = np.log10(np.array([np.sum(10**nonparametric['formed_mass'][i][9:])
                                   for i in range(len(nonparametric))]))
old_component_hacked = np.log10(np.array([np.sum(10**nonparametric_hacked['formed_mass'][i][9:])
                                          for i in range(len(nonparametric_hacked))]))

plt.scatter(old_component-old_component_hacked, nonparametric['stellar_mass']-nonparametric_hacked['stellar_mass'])