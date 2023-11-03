"""

ariel@berlim
27/05/2022

Plots the full sampling of the nonparametric SFHs

"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from toolbox import plot_tools
import seaborn as sns

sns.set_style('ticks')

detection = 'halpha'
config = 'sfrprior_nonparametric'

sfh_dir = '/home/ariel/Workspace/GASP/HST/BAGPIPES/sfh/'
blob_catalog = Table.read('/home/ariel/Workspace/GASP/HST/Data/disk_halpha_bagpipes_input.fits')

age_bins = [0, 2.0e6, 4.0e6, 7.0e6, 2.0e7, 5.5e7, 2.0e8, 5.5e8, 1.0e9, 3.0e9, 5.75e9, 1.0e10, 1.4e10]
age_bins_len = [age_bins[i+1]-age_bins[i] for i in range(12)]
age_bins_mid = [(age_bins[i+1]+age_bins[i])/2 for i in range(12)]

colors = sns.color_palette('rainbow', n_colors=len(age_bins))

for blob_id in blob_catalog['blob_id'][0:120]:
    print(blob_id)

    masses_full = np.genfromtxt(sfh_dir + detection + '_' + config + '_disk/masses_' + blob_id + '_' + config + '_full.txt')

    fig, ax = plt.subplots(4, 3, figsize=(10, 8))

    for i in range(len(age_bins_mid)):
        ax.ravel()[i].hist(masses_full[i], bins=25, color=colors[i])
        ax.ravel()[i].annotate(str(age_bins_mid[i]/1e9), xy=(0.7, 0.9), xycoords='axes fraction')

    plt.savefig('/home/ariel/Workspace/GASP/HST/BAGPIPES/plots/nonparametric_distributions/' + blob_id + '.png')

    plt.close()