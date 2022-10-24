"""

ariel@oapd
13/04/2022

Compare SFH of parametric and nonparametric fits with and without nebular continuum emission

"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from toolbox import plot_tools
import seaborn as sns

sns.set_style('ticks')

sfh_dir = '/home/ariel/Workspace/GASP/HST/BAGPIPES/sfh/'
blob_catalog = Table.read('/home/ariel/Workspace/GASP/HST/Data/disk_halpha_bagpipes_input.fits')

age_bins = [0, 2.0e6, 4.0e6, 7.0e6, 2.0e7, 5.5e7, 2.0e8, 5.5e8, 1.0e9, 3.0e9, 5.75e9, 1.0e10, 1.4e10]
age_bins_len = [age_bins[i+1]-age_bins[i] for i in range(12)]
age_bins_mid = [(age_bins[i+1]+age_bins[i])/2 for i in range(12)]

for i in range(200):

    print(i)

    blob_id = blob_catalog['blob_id'][i]
    print(blob_id)

    try:
        age, sfr, sfr25, sfr75 = np.genfromtxt(sfh_dir + 'halpha_parametric_disk/sfh_'
                                               + blob_id + '_parametric.txt').transpose()
        age_npar, sfr_npar, sfr25_npar, sfr75_npar = np.genfromtxt(sfh_dir + 'halpha_sfrprior_nonparametric_disk/sfh_'
                                                                   + blob_id +
                                                                   '_sfrprior_nonparametric_bagpipes.txt').transpose()
        age_npar_df, sfr_npar_df, sfr25_npar_df, sfr75_npar_df = np.genfromtxt(sfh_dir + 'halpha_default_nonparametric_disk/sfh_'
                                                                               + blob_id +
                                                                               '_default_nonparametric_bagpipes.txt').transpose()
        age_npar_man, sfr_npar_man, sfr25_npar_man, sfr75_npar_man = np.genfromtxt(sfh_dir + 'halpha_sfrprior_nonparametric_disk/sfh_'
                                                                                   + blob_id +
                                                                                   '_sfrprior_nonparametric.txt').transpose()

    except Exception:
        print('Skipped')
        continue

    fig = plt.figure()
    plt.plot(np.log10(age), sfr, color='steelblue', label='Parametric')
    plt.plot(np.log10(age_npar), sfr_npar, color='firebrick', label='Nonparametric')
    # plt.plot(np.log10(age_npar), sfr_npar_df, color='forestgreen', label='Nonparametric - Mass prior')
    plt.scatter(np.log10(age_npar_man), sfr_npar_man, color='k', label='Nonparametric - Calculated from masses formed')

    plt.legend()

    sns.despine()
    plt.xlabel(r'$\log\,t\,\mathrm{[yr]}$')
    plt.ylabel(r'$\mathrm{SFR [M_\odot/yr]}$')

    fig.tight_layout()

    plt.savefig('/home/ariel/Workspace/GASP/HST/BAGPIPES/plots/sfh_comp_disk/' + blob_id + '.png', dpi=200)

    plt.close()