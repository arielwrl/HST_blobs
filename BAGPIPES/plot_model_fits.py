"""

ariel@oapd
20/06/2022

Plots results for model fits

"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from astropy.table import Table
from matplotlib import cm

sns.set_style('ticks')

data_dir = '/home/ariel/Workspace/GASP/HST/Data/models/'

tail = Table.read(data_dir + '5err_doublescreen_eta_vary.fits')
old_mass_25 = Table.read(data_dir + '5err_doublescreen_eta_vary_oldmass25.fits')
old_mass_5 = Table.read(data_dir + '5err_doublescreen_eta_vary_oldmass5.fits')
old_mass_10 = Table.read(data_dir + '5err_doublescreen_eta_vary_oldmass10.fits')

table_list = [tail, old_mass_25, old_mass_5, old_mass_10]
label_list = ['Tail', 'Old Mass = 2.5', 'Old Mass = 5', 'Old Mass = 10']

fig, ax_grid = plt.subplots(2, 2, figsize=(8, 8))

ax = ax_grid.ravel()

colors = [cm.rainbow(x) for x in np.linspace(0, 1, 5)]
colordict = dict(zip(np.unique(tail['model_Av']), colors))

for i in range(4):

    table = table_list[i]

    scatter = ax[i].scatter(table['model_age'], (table['model_age'] - table['age']) / table['model_age'],
                            c=table['model_Av'], cmap='rainbow')
    cb = plt.colorbar(scatter, ax=ax[i], label='Model Av')

    for Av in np.unique(table['model_Av']):
        flag = table['model_Av'] == Av
        ax[i].plot(table['model_age'][flag], (table['model_age'] - table['age'])[flag] / table['model_age'][flag],
                   c=colordict[Av])

    ax[i].set_ylabel('(Model Age - Output Age) / Model Age')
    ax[i].set_xlabel('Model Age')

    ax[i].set_title(label_list[i])

    sns.despine()

plt.show()