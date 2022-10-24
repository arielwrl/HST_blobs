"""

ariel@oapd
08/11/2021

Include Av from MUSE in Eric's blog catalogs

ariel@oapd
03/05/2022

Updated as a general input creator, lines that include Av from MUSE are commented out

ariel@oapd
29/09/2022

Updated to include extraplanar clumps


"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table
from pysinopsis.masking import create_mask

plt.ioff()

detection = 'optical_only'
selection = 'tail'
data_dir = '/home/ariel/Workspace/GASP/HST/Data/'

if detection == 'f606w':
    blob_catalog = Table.read(data_dir + '/HST_' + detection + '_0.04px_reg_mask_denoised_disk275_masked_3sigma_0delta_'
                                                               'only_trunks_match_all_gals.csv',
                              format='ascii.ecsv')
elif detection == 'optical_only':
    blob_catalog = Table.read(data_dir + 'HST_f606w_0.04px_reg_mask_denoised_disk275_masked_3sigma_0delta_only_trunks_'
                                         'match_optical_only.csv', format='ascii.ecsv')
else:
    blob_catalog = Table.read(data_dir + '/HST_' + detection + '_0.04px_reg_mask_clumps_2.5sigma_5delta_'
                                                               'denoise_2sigma_0delta_all_gals_z_BPT_sel.csv',
                              format='ascii.ecsv')


if selection == 'disk':
    blob_catalog = blob_catalog[blob_catalog['tail_gal_flag'] == 2]
if (selection == 'tail') & (detection != 'f606w'):
    blob_catalog = blob_catalog[(blob_catalog['tail_gal_flag'] == 0) | (blob_catalog['tail_gal_flag'] == 1)]
if (selection == 'tail') & (detection == 'f606w'):
    blob_catalog = blob_catalog[blob_catalog['tail_gal_flag'] == 0]


if (detection == 'f606w') | (detection == 'optical_only'):
    blob_catalog['fit_id'] = [blob_catalog['gal'][i] + '_' + str(blob_catalog['f606w_idx'][i]) + '_'
                              + detection + '_' + selection for i in range(len(blob_catalog))]
    blob_catalog['blob_id'] = [blob_catalog['gal'][i] + '_' + str(blob_catalog['f606w_idx'][i]) + '_'
                               + detection for i in range(len(blob_catalog))]
else:
    blob_catalog['fit_id'] = [blob_catalog['gal'][i] + '_' + blob_catalog['id_catalog'][i] +
                              str(blob_catalog['_idx'][i]) + '_' + detection + '_' + selection
                              for i in range(len(blob_catalog))]
    blob_catalog['blob_id'] = [blob_catalog['gal'][i] + '_' + blob_catalog['id_catalog'][i] +
                               str(blob_catalog['_idx'][i])
                               for i in range(len(blob_catalog))]

redshift_dict = {'JO201': 0.044631,
                 'JO204': 0.042372,
                 'JO206': 0.051089,
                 'JW100': 0.061891,
                 'JW39': 0.066319,
                 'JO175': 0.046750}

blob_catalog['galaxy_redshift'] = [redshift_dict[galaxy] for galaxy in blob_catalog['gal']]

print(blob_catalog)

if detection == 'optical_only':
    blob_catalog = blob_catalog[blob_catalog['F606W'] != 0]

blob_catalog.write(data_dir + selection + '_' + detection + '_bagpipes_input.fits', overwrite=True)

# wcs_hst_dict = {}
# wcs_muse_dict = {}
# av_map_dict = {}
#
# for galaxy in np.unique(blob_catalog['gal']):
#
#     print('Reading files for galaxy ' + galaxy)
#
#     f606_image = fits.open('/home/ariel/Workspace/GASP/HST/Data/' + galaxy + '/' + galaxy +
#                            '_f606w_0.04px_drc_MWdc.fits')
#     muse_fov = fits.open('/home/ariel/Workspace/GASP/HST/Data/' + galaxy + '/' + galaxy + '_IMAGE_FOV_0001_v1.fits')
#
#     wcs_hst = WCS(f606_image[1].header)
#     wcs_muse = WCS(muse_fov[1].header)
#
#     emission_cube = fits.open('/home/ariel/Workspace/GASP/HST/Data/' + galaxy + '/' + galaxy + '_eo_dc.fits')
#     av_map = np.ma.masked_array(emission_cube[0].data[129, :, :], mask=emission_cube[0].data[130, :, :] != 0)
#
#     wcs_hst_dict[galaxy] = wcs_hst
#     wcs_muse_dict[galaxy] = wcs_muse
#     av_map_dict[galaxy] = av_map
#
# avs = {}
# coordinates = {}
# out_of_muse = np.zeros(len(blob_catalog), dtype=bool)

# for blob_id in blob_catalog['blob_id']:
#     print('>>> Finding Av for ' + blob_id)
#
#     blob = blob_catalog[blob_catalog['blob_id'] == blob_id]
#
#     wcs_hst = wcs_hst_dict[galaxy]
#     wcs_muse = wcs_muse_dict[galaxy]
#     av_map = av_map_dict[galaxy]
#
#     sky = wcs_hst.pixel_to_world(blob['xc_pix(HST)'], blob['yc_pix(HST)'])
#     x_muse, y_muse = wcs_muse.world_to_pixel(sky)
#
#     if (x_muse > av_map.shape[1]) | (x_muse < 0) | (y_muse > av_map.shape[0]) | (y_muse < 0):
#         avs[blob_id] = (-999, -999)
#         coordinates[blob_id] = (x_muse, y_muse)
#         out_of_muse[np.argwhere(blob_catalog['blob_id'] == blob_id)] = True
#     else:
#         mask = create_mask(x_muse, y_muse, 1.5*blob['iso_radius_pix(HST)'], av_map)
#         av_blob = np.ma.masked_array(av_map, mask=(mask == 1))
#
#         if np.all(av_blob.mask):
#             avs[blob_id] = (0, 0)
#
#         else:
#             avs[blob_id] = (np.mean(av_blob), np.std(av_blob))
#             coordinates[blob_id] = (x_muse, y_muse)
#
#         plt.title(str(x_muse)+','+str(y_muse)+','+str((~av_blob.mask).sum()))
#         plt.savefig('/home/ariel/Workspace/GASP/HST/Data/mask_plots/' + galaxy + blob_id + '.png')
#         plt.clf()

# blob_catalog['av_muse'] = np.array([avs[blob_id][0] for blob_id in blob_catalog['blob_id']])
# blob_catalog['av_muse_std'] = np.array([avs[blob_id][1] for blob_id in blob_catalog['blob_id']])
# blob_catalog['out_of_muse'] = out_of_muse


