"""

ariel@padova
18/01/2021

"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.io import fits
from astropy.wcs import WCS
import pickle
from matplotlib.patches import Circle
from pysinopsis.masking import create_mask
from astropy.wcs import utils
from astropy.wcs.utils import proj_plane_pixel_scales

blob_catalog = Table.read('/home/ariel/Workspace/GASP/HST/Data/all_blobs_infos_allinfo_new_v2_dist.cat', format='ascii')
JO204_blobs = blob_catalog[(blob_catalog['gasp_gal'] == 'JO204') & (blob_catalog['flag_out_sn'] == 1)]
JO204_blobs['ID'] = [JO204_blobs['gasp_gal'][i] + '_' + str(int(JO204_blobs['blob_num'][i]))
                     for i in range(len(JO204_blobs))]

photo_cube = pickle.load(open('/home/ariel/Workspace/GASP/HST/Data/JO204/JO204_photocube.pkl', 'rb'))
emission_cube = fits.open('/home/ariel/Workspace/GASP/HST/Data/JO204/JO204_eo_dc.fits')

av_map_unmasked = emission_cube[0].data[129, :, :]
av_flag = emission_cube[0].data[130, :, :]
av_map = np.ma.masked_array(av_map_unmasked, mask=av_flag != 0)

muse_fov = fits.open('/home/ariel/Workspace/GASP/HST/Data/JO204/JO204_IMAGE_FOV_0001_v1.fits')

wcs_hst = WCS(photo_cube['f606w_header'])
wcs_muse = WCS(muse_fov[1].header)

fig = plt.figure(figsize=(20, 10))
ax1 = fig.add_subplot(121, projection=wcs_muse)
ax2 = fig.add_subplot(122, projection=wcs_muse, sharex=ax1, sharey=ax1)

av = ax1.imshow(av_map, vmin=0, vmax=1, cmap='Reds')
image = ax2.imshow(photo_cube['f606w_flux'], vmax=1e-22, transform=ax2.get_transform(wcs_hst), cmap='Blues')

cb1 = plt.colorbar(mappable=av, ax=ax1)
cb1.set_label('$A_V$', fontsize=20)

cb2 = plt.colorbar(mappable=image, ax=ax2)
cb2.set_label('F606W Flux', fontsize=20)

for i in range(len(JO204_blobs)):

    circle1 = Circle(xy=(JO204_blobs['x'][i], JO204_blobs['y'][i]), radius=JO204_blobs['r_pix'][i], fill=False,
                     linewidth=2)
    circle2 = Circle(xy=(JO204_blobs['x'][i], JO204_blobs['y'][i]), radius=JO204_blobs['r_pix'][i], fill=False,
                     linewidth=2)

    ax1.add_patch(circle1)
    ax2.add_patch(circle2)

    ax1.annotate(int(JO204_blobs['blob_num'][i]), xy=(JO204_blobs['x'][i], JO204_blobs['y'][i]+JO204_blobs['r_pix'][i]),
                 color='k', fontsize=20)
    ax2.annotate(int(JO204_blobs['blob_num'][i]), xy=(JO204_blobs['x'][i], JO204_blobs['y'][i]+JO204_blobs['r_pix'][i]),
                 color='k', fontsize=20)

ax1.set_xlim(50, 280)
ax1.set_ylim(70, 330)

ax1.set_xlabel('RA', fontsize=18)
ax2.set_xlabel('RA', fontsize=18)

ax1.set_ylabel('Dec', fontsize=18)
ax2.set_ylabel('Dec', fontsize=18)

fig.tight_layout()

# plt.savefig('map_image.png', dpi=200)


########################################################################################################################
# Testing Blob
########################################################################################################################
#
# fig = plt.figure(figsize=(20, 10))
#
# plt.imshow(photo_cube['f606w_flux'], vmin=0.1, vmax=1, cmap='Blues')
#
# for blob_id in JO204_blobs['ID']:
#
#     print(blob_id)
#
#     blob = JO204_blobs[JO204_blobs['ID'] == blob_id]
#
#     sky = wcs_muse.pixel_to_world(blob['x'], blob['y'])
#     x_hst, y_hst = wcs_hst.world_to_pixel(sky)
#
#     # mask = create_mask(blob['x'], blob['y'], blob['r_pix'], av_map)
#     # av_blob = np.ma.masked_array(av_map, mask=(mask == 1))
#
#     # mask_hst = create_mask(x_hst, y_hst, blob['r_pix']*5, photo_cube['f606w_flux'])
#     # hst_blob = np.ma.masked_array(photo_cube['f606w_flux'], mask=(mask_hst == 1))
#
#     # ax[0].imshow(av_map, cmap='Reds', vmax=1)
#     # circle = Circle(xy=(blob['x'], blob['y']), radius=blob['r_pix'], fill=False,
#     #                 linewidth=2.5)
#     # ax[0].add_patch(circle)
#     # ax[0].imshow(mask, alpha=0.5)
#
#     circle = Circle(xy=(x_hst, y_hst), radius=blob['r_pix']*5, fill=False,
#                     linewidth=2.5)
#     plt.gca().add_patch(circle)
#     # plt.imshow(mask_hst, alpha=0.5)
#
# plt.xlim(3100, 4300)
# plt.ylim(2100, 3500)
#
# plt.savefig('test_correspondence.png', dpi=200)
#
# #
# # plt.figure()
# # plt.imshow(photo_cube['f606w_flux'], cmap='viridis')
# # plt.imshow(mask_hst, cmap='Reds')
# # plt.savefig('testmask.png')
# #
# #
# #





#from MUSE pix to RADEC using wcs_MUSE
#from RADEC to HST pix using wcs_HST

# x_muse=[100,200,300]
# y_muse=[100,200,300]
# ra,de=wcs_MUSE.all_pix2world(x_muse,y_muse,0)
# x_hst,y_hst=wcs_HST.all_world2pix(ra,de,0)
#
# print (x_hst,y_hst)