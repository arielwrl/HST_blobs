"""

ariel@isa
18/01/2022

"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from toolbox import plot_tools
from astropy.io import fits

f606_image_raw = fits.open('/home/ariel/Workspace/GASP/HST/Data/images/JO206_f606w_0.04px_drc_MWdc.fits')
f606_image = np.ma.masked_array(f606_image_raw[1].data, mask=(f606_image_raw[1].data < 0))

# input = Table.read('/home/ariel/Workspace/GASP/HST/Data/old/disk_halpha_bagpipes_input.fits')
output = Table.read('/home/ariel/Workspace/GASP/HST/Data/old/disk_halpha_parametric_bagpipes_results.fits')
matched = Table.read('/home/ariel/Workspace/GASP/HST/Data/old/disk_matched.fits')

galaxy_flag = matched['galaxy'].astype(str) == 'JO206'
matched = matched[galaxy_flag]

plt.imshow(np.log10(f606_image), cmap='Greys', origin='lower')
age_map = plt.scatter(matched['xc_pix(HST)'], matched['yc_pix(HST)'], c=matched['formed_mass_old'], s=50,
                      cmap='Spectral_r')
plt.colorbar(age_map, label='oldmass')

plt.figure()
