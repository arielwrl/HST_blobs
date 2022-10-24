"""

ariel@oapd
22/03/2022

Hack the nebular part of BAGPIPES, replacing all fluxes with zero

"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

original_grid = fits.open('/home/ariel/.local/lib/python3.8/site-packages/bagpipes/models/grids/bc03_miles_nebular_cont'
                          '_grids.fits')

zeros = [np.zeros_like(original_grid[i+1].data) for i in range(35)]

hdu_list = [original_grid[0]]

for i in range(35):

    hdu = fits.ImageHDU(zeros[i], header=original_grid[i+1].header)
    hdu_list.append(hdu)

zeros_grid = fits.HDUList(hdu_list)

zeros_grid.writeto('/home/ariel/.local/lib/python3.8/site-packages/bagpipes/models/grids/'
                   'hacked_no_cont_bc03_miles_nebular_cont_grids.fits')

