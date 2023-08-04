"""

ariel@sardegna
30/08/2022


"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from toolbox import plot_tools
from astropy.io import fits

f606_image_raw = fits.open('/home/ariel/Workspace/GASP/HST/Data/images/JO206_f606w_0.04px_drc_MWdc.fits')
f606_image = np.ma.masked_array(f606_image_raw[1].data, mask=(f606_image_raw[1].data < 0))

input = Table.read('/home/ariel/Workspace/GASP/HST/Data/tail_f275w_bagpipes_input.fits')
output = Table.read('/home/ariel/Workspace/GASP/HST/Data/tail_f275w_dexp_logprior_bagpipes_results.fits')

input = input[input['sel_flag'] == 31]
output = output[input['sel_flag'] == 31]

galaxy_flag = (input['gal'] == 'JO206')

input = input[galaxy_flag]
output = output[galaxy_flag]

plt.imshow(np.log10(f606_image), cmap='Greys', origin='lower')
plt.scatter(input['xc_pix(HST)'], input['yc_pix(HST)'], s=50, c=output['mwage'],
            vmax=0.05,
            vmin=0,
            cmap='Spectral_r')
plt.colorbar()

plt.show()
