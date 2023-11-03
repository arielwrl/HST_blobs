"""

ariel@oapd
06/01/2023

Pick up problems in double component fits and flag them out

"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.table import join
from astropy.io import fits

extraplanar_input_halpha = Table.read('/home/ariel/Workspace/GASP/HST/Data/extraplanar_halpha_bagpipes_input.fits')
extraplanar_double_halpha = Table.read('/home/ariel/Workspace/GASP/HST/Data/extraplanar_halpha_dexp_logprior_'
                                       'bagpipes_results.fits')
extraplanar_single_halpha = Table.read('/home/ariel/Workspace/GASP/HST/Data/tail_halpha_dexp_logprior_'
                                       'bagpipes_results.fits')
extraplanar_single_input_halpha = Table.read('/home/ariel/Workspace/GASP/HST/Data/tail_halpha_bagpipes_input.fits')
extraplanar_single_input_halpha = extraplanar_single_input_halpha[extraplanar_single_input_halpha['sel_flag'] == 31]
# extraplanar_double_f275w = Table.read('/home/ariel/Workspace/GASP/HST/Data/extraplanar_f275w_dexp_logprior_'
#                                        'bagpipes_results.fits')
# extraplanar_single_f275w = Table.read('/home/ariel/Workspace/GASP/HST/Data/tail_f275w_dexp_logprior_'
#                                        'bagpipes_results.fits')

extraplanar_double_halpha['mass_ratio'] = np.log10((extraplanar_double_halpha['formed_mass_old']**10)/
                                                   (extraplanar_double_halpha['formed_mass_young']**10))

extraplanar_flag = extraplanar_single_input_halpha['tail_gal_flag'] == 1

unconstrained_double_fit = []
unconstrained_single_fit = []
high_mass_ratio = []

for i in range(len(extraplanar_double_halpha)):

    if extraplanar_double_halpha['age_young_iqr'][i] > 0.2:
        unconstrained_double_fit.append(True)
    else:
        unconstrained_double_fit.append(False)

    if extraplanar_single_halpha['age_iqr'][extraplanar_flag][i] > 0.2:
        unconstrained_single_fit.append(True)
    else:
        unconstrained_single_fit.append(False)

    if extraplanar_double_halpha['mass_ratio'][i] > 1:
        high_mass_ratio.append(True)
    else:
        high_mass_ratio.append(False)

extraplanar_single_halpha['unconstrained_double_fit_flag'] = np.zeros_like(extraplanar_single_halpha['age'], dtype=bool)
extraplanar_single_halpha['unconstrained_single_fit_flag'] = np.zeros_like(extraplanar_single_halpha['age'], dtype=bool)
extraplanar_single_halpha['high_mass_ratio'] = np.zeros_like(extraplanar_single_halpha['age'], dtype=bool)

extraplanar_single_halpha['unconstrained_double_fit_flag'][extraplanar_flag] = unconstrained_double_fit
extraplanar_single_halpha['unconstrained_single_fit_flag'][extraplanar_flag] = unconstrained_single_fit
extraplanar_single_halpha['high_mass_ratio'][extraplanar_flag] = high_mass_ratio

extraplanar_single_halpha.write('/home/ariel/Workspace/GASP/HST/Data/tail_halpha_dexp_logprior_bagpipes_results.fits',
                                overwrite=True)
# Age vs Age
#
# plt.figure()
#
# plt.scatter(extraplanar_double_halpha['age_young'], extraplanar_single_halpha['age'], c=extraplanar_double_halpha['mass_ratio'], vmin=-3, vmax=1)
# plt.colorbar(label='Log Mass Ratio (Old/Young)')
#
# plt.xlabel('Age (two component fit)')
# plt.ylabel('Age (single component fit)')



# JW100 example:
#
# f606_image_raw = fits.open('/home/ariel/Workspace/GASP/HST/Data/images/JO206_f606w_0.04px_drc_MWdc.fits')
# f606_image = np.ma.masked_array(f606_image_raw[1].data, mask=(f606_image_raw[1].data < 0.1e-21))
#
# plt.figure()
#
# lim_flag = extraplanar_double_halpha['mass_ratio'] < 6
#
# input = extraplanar_input_halpha[extraplanar_input_halpha['sel_flag'] == 31]
# galaxy_flag = input['gal'] == 'JO206'
# suspicious_flag = (input['leaf_flag'] == 1) & (input['level'] > 0)
# input_suspicious = input[galaxy_flag & lim_flag & suspicious_flag]
# input = input[galaxy_flag & lim_flag]
# output = extraplanar_double_halpha[galaxy_flag & lim_flag]
# mass_ratio = extraplanar_double_halpha['mass_ratio'][galaxy_flag & lim_flag]
# mass_ratio_suspicious_leaves = extraplanar_double_halpha['mass_ratio'][galaxy_flag & lim_flag & suspicious_flag]
#
# plt.imshow(np.log10(f606_image), cmap='Greys', origin='lower')
# plt.scatter(input['xc_pix(HST)'], input['yc_pix(HST)'], s=30,
#             c=mass_ratio, cmap='viridis', vmin=-3, vmax=1)
# plt.scatter(input_suspicious['xc_pix(HST)'], input_suspicious['yc_pix(HST)'], s=30,
#             facecolors=None, edgecolors='red')
#
#
# plt.colorbar(label='Log Mass Ratio (Old/Young)')


