"""

ariel@isa
18/01/2022

"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from toolbox import plot_tools
from astropy.io import fits

f606_image_raw = fits.open('/home/ariel/Workspace/GASP/HST/Data/JO206/JO206_f606w_0.04px_drc_MWdc.fits')
f606_image = np.ma.masked_array(f606_image_raw[1].data, mask=(f606_image_raw[1].data < 0))

input = Table.read('/home/ariel/Workspace/GASP/HST/Data/tail_all_f606w_bagpipes_input.fits')
output = Table.read('/home/ariel/Workspace/GASP/HST/Data/tail_all_default_f606w_bagpipes_results.fits')

galaxy_flag = input['gal'].astype(str) == 'JO206'

input = input[galaxy_flag]
output = output[galaxy_flag]

plt.imshow(np.log10(f606_image), cmap='Greys', origin='lower')
# age_map = plt.hexbin(input['xc_pix(HST)'], input['yc_pix(HST)'], C=output['age'], mincnt=1,
#                              gridsize=10)
age_map = plt.scatter(input['xc_pix(HST)'], input['yc_pix(HST)'], c=output['mwage'], s=50,
                      cmap='Spectral_r')
plt.colorbar(age_map, label='mwage')
