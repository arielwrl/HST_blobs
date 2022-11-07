"""

ariel@padova
25/01/2021

Some code to hopefully help Eric with hst dust correction

"""

import numpy as np
from astropy.table import Table
from astropy.io import fits
import matplotlib.pyplot as plt


########################################################################################################################
# Functions
########################################################################################################################

def ccm(wl, r_v=3.1):
    """
    Calculates the Cardelli, Clayton & Mathis (CCM) extinction curve.
    Note that wl has to be provided as array.

    Parameters
    ----------
    wl : array like
        Wavelength in Angstroms.
    r_v : float, optional
        Ratio of total to selective extinction.

    Returns
    -------
    q : array
        Extinction A_lambda / A_V. Array of same length as wl.
    Reference: http://adsabs.harvard.edu/abs/1989ApJ...345..245C

    """

    a = np.zeros_like(wl)
    b = np.zeros_like(wl)
    f_a = np.zeros_like(wl)
    f_b = np.zeros_like(wl)
    x = np.zeros_like(wl)
    y = np.zeros_like(wl)
    q = np.zeros_like(wl)

    x = 10000. / wl
    y = 10000. / wl - 1.82

    # Far-Ultraviolet: 8 <= x <= 10 ; 1000 -> 1250 Angstrom
    i = np.bitwise_and(x >= 8, x <= 10)

    a[i] = -1.073 - 0.628 * (x[i] - 8.) + 0.137 * (x[i] - 8.)**2 - 0.070 * (x[i] - 8.)**3
    b[i] = 13.670 + 4.257 * (x[i] - 8.) - 0.420 * (x[i] - 8.)**2 + 0.374 * (x[i] - 8.)**3

    # Ultraviolet: 3.3 <= x <= 8 ; 1250 -> 3030 Angstrom
    i = np.bitwise_and(x >= 5.9, x < 8)
    f_a[i] = -0.04473 * (x[i] - 5.9)**2 - 0.009779 * (x[i] - 5.9)**3
    f_b[i] = 0.2130 * (x[i] - 5.9)**2 + 0.1207 * (x[i] - 5.9)**3

    i = np.bitwise_and(x >= 3.3, x < 8)

    a[i] = 1.752 - 0.316 * x[i] - 0.104 / ((x[i] - 4.67)**2 + 0.341) + f_a[i]
    b[i] = -3.090 + 1.825 * x[i] + 1.206 / ((x[i] - 4.62)**2 + 0.263) +f_b[i]

    # Optical/NIR: 1.1 <= x <= 3.3 ; 3030 -> 9091 Angstrom
    i = np.bitwise_and(x >= 1.1, x < 3.3)

    a[i] = 1.+ 0.17699 * y[i] - 0.50447 * y[i]**2 - 0.02427 * y[i]**3 + \
        0.72085 * y[i]**4 + 0.01979 * y[i]**5 - 0.77530 * y[i]**6 + 0.32999 * y[i]**7
    b[i] = 1.41338 * y[i] + 2.28305 * y[i]**2 + 1.07233 * y[i]**3 - \
        5.38434 * y[i]**4 - 0.62251 * y[i]**5 + 5.30260 * y[i]**6 - 2.09002 * y[i]**7

    # Infrared: 0.3 <= x <= 1.1 ; 9091 -> 33333 Angstrom
    i = np.bitwise_and(x >= 0.3, x < 1.1)

    a[i] = 0.574 * x[i]**1.61
    b[i] = -0.527 * x[i]**1.61

    q = a + b / r_v

    return q


def calc_extinction(wl, ebv, r_v=3.1):
    """

    Calculates galactic extinction at a given wavelength

    Parameters
    ----------
    wl: array like
        Wavelengths in angstroms
    r_v (optional): float
        Ratio of total to selective extinction, default is the average MW value of 3.1
    ebv: E(B-V) from dust map

    Returns
    -------
    A_lambda: array
        Galactic extinction at wl
    """

    a_v = r_v * ebv
    a_lambda = a_v * ccm(wl)

    return a_lambda


########################################################################################################################
# Script
########################################################################################################################

data_path = '/home/ariel/Workspace/GASP/HST/Data/'
galaxy_id = 'JO204'

hst_filters = ['f275w', 'f336w', 'f606w', 'f680n', 'f814w']

pivot_wavelengths = {'f275w': 2703.57323,
                     'f336w': 3354.83065,
                     'f606w': 5884.58579,
                     'f680n': 6877.42566,
                     'f814w': 8043.67929
                     }

dust_table = Table.read('dust_new.fits')
dust_table['id'] = np.array([dust_table['id'][i].strip() for i in range(len(dust_table))])
# The above will "clean" the blank spaces from the strings in the table so we can apply the next step
ebv_galaxy = dust_table['eBV_SF'][dust_table['id'] == galaxy_id].data

corrected_fluxes = {}
uncorrected_fluxes = {}

for hst_filter in hst_filters:

    hdu = fits.open(data_path + galaxy_id + '/' + galaxy_id + '_' + hst_filter + '_0.04px_north_drc.fits')

    uncorrected_flux = hdu[1].data
    uncorrected_fluxes[hst_filter] = uncorrected_flux

    extinction_at_pivot_wavelength = calc_extinction(wl=np.array(pivot_wavelengths[hst_filter]), ebv=ebv_galaxy)

    corrected_flux = uncorrected_flux * 10 ** (0.4 * extinction_at_pivot_wavelength)
    corrected_fluxes[hst_filter] = corrected_flux


# fig, ax = plt.subplots(1, 5, figsize=(25, 10))
#
# x = np.linspace(0, 600)
# y = x
#
# for i in range(5):
#
#     ax[i].scatter(uncorrected_fluxes[hst_filters[i]].ravel(), corrected_fluxes[hst_filters[i]].ravel(), alpha=0.5)
#     ax[i].plot(x, y, '-k')
#
#     ax[i].set_title(hst_filters[i], fontsize=20)
#     ax[i].set_xlabel('Uncorrected Flux', fontsize=20)
#     ax[i].set_ylabel('Corrected Flux', fontsize=20)
#
# fig.tight_layout()
# plt.savefig('dust_correction_test.png', dpi=150)

