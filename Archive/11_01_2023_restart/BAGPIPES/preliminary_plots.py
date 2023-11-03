"""

ariel@padova
29/03/2021

Plot Eric's masks and find them in the Av map

"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from matplotlib.patches import Circle
from astropy.table import Table

galaxy = 'JO206'
data_dir = '/home/ariel/Workspace/GASP/HST/Data/'

f606_image = fits.open('/home/ariel/Workspace/GASP/HST/Data/' + galaxy + '/' + galaxy +
                       '_f606w_0.04px_drc_MWdc.fits')
muse_fov = fits.open('/home/ariel/Workspace/GASP/HST/Data/' + galaxy + '/' + galaxy + '_IMAGE_FOV_0001_v1.fits')

wcs_hst = WCS(f606_image[1].header)
wcs_muse = WCS(muse_fov[1].header)

emission_cube = fits.open('/home/ariel/Workspace/GASP/HST/Data/' + galaxy + '/' + galaxy + '_eo_dc.fits')
av_map = np.ma.masked_array(emission_cube[0].data[129, :, :], mask=emission_cube[0].data[130, :, :] != 0)

blob_catalog = Table.read(data_dir + galaxy + '/' + galaxy + '_f275w_0.04px_reg_mask_2.5sigma_5delta_denoise_2sigma_'
                                                             '0delta_dendrogram_auto_cleaned_newcoords.csv'
                          , format='ascii.ecsv')
blob_catalog = blob_catalog[blob_catalog['tail_gal_flag'].astype(bool) & (blob_catalog['sel_flag'] == 31)]

fig = plt.figure(figsize=(20, 10))
ax1 = fig.add_subplot(121, projection=wcs_hst)
ax2 = fig.add_subplot(122, projection=wcs_hst, sharex=ax1, sharey=ax1)

image = ax1.imshow(f606_image[1].data / 1e-18, vmin=0, vmax=0.02, cmap='Blues')
av = ax2.imshow(av_map, cmap='Reds', vmax=0.8, transform=ax2.get_transform(wcs_muse))

cb1 = plt.colorbar(mappable=image, ax=ax1)
cb1.set_label('F680N Flux $[\mathrm{10^{-18} erg/s/cm^2/\AA}]$', fontsize=20)

cb2 = plt.colorbar(mappable=av, ax=ax2)
cb2.set_label('$A_V \, [\mathrm{mag}]$', fontsize=20)

for i in range(len(blob_catalog)):
    circle1 = Circle(xy=(blob_catalog['x_cen(pix)'][i], blob_catalog['y_cen(pix)'][i]),
                     radius=blob_catalog['iso_radius'][i], fill=False, linewidth=2.5, color='red')
    circle2 = Circle(xy=(blob_catalog['x_cen(pix)'][i], blob_catalog['y_cen(pix)'][i]),
                     radius=1.5 * blob_catalog['iso_radius'][i], fill=False, linewidth=2.5, color='blue')

    ax1.add_patch(circle1)
    ax2.add_patch(circle2)

    identifier = blob_catalog['id_catalog'][i] + str(blob_catalog['_idx'][i])

    ax1.annotate(identifier, xy=(blob_catalog['x_cen(pix)'][i], blob_catalog['y_cen(pix)'][i]
                 + blob_catalog['iso_radius'][i]), color='k', fontsize=14)
    ax2.annotate(identifier, xy=(blob_catalog['x_cen(pix)'][i], blob_catalog['y_cen(pix)'][i]
                 + blob_catalog['iso_radius'][i]), color='k', fontsize=14)

ax1.set_xlim(0, f606_image[1].data.shape[0])
ax1.set_ylim(0, f606_image[1].data.shape[1])

ax1.set_xlabel('RA', fontsize=18)
ax2.set_xlabel('RA', fontsize=18)

ax1.set_ylabel('Dec', fontsize=18)
ax2.set_ylabel('Dec', fontsize=18)

fig.tight_layout()

plt.show()
