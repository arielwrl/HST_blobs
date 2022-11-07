"""

ariel@padova
01/01/2021

Generate a quasi-datacube from HST images

From Eric's email:

JO201: xmin=1250, xmax=4150, ymin=700, ymax=3600
JO204: xmin=2950, xmax=4600, ymin=1950, ymax=3600
JW100: xmin=2400, xmax=4500, ymin=2400, ymax=4500
JW39: xmin=950, xmax=3950, ymin=1550, ymax=4450

[1950:3600, 2950:4600]

"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from collections import OrderedDict
import pickle

galaxy = 'JO206'

filters = ['f275w', 'f336w', 'f606w', 'f680n', 'f814w']

data_dir = '/home/ariel/Workspace/GASP/HST/Data/' + galaxy + '/'

photo_cube = OrderedDict()

for hst_filter in filters:
    hdu = fits.open(data_dir + galaxy + '_' + hst_filter + '_0.04px_drc_subreg_MWdc.fits')

    image = np.ma.masked_array(hdu[1].data, mask=hdu[1].data < 0.0)

    photo_cube[hst_filter + '_flux'] = image
    photo_cube[hst_filter + '_header'] = hdu[1].header

pickle.dump(photo_cube, open(data_dir + galaxy + '_photocube.pkl', 'wb'))
