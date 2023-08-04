"""

ariel@oapd
11/01/2023

Converts Eric's catalogs into fits, adds galaxy redshift and clump ID, selects only 31 clumps


"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from toolbox import stat_tools
from toolbox.wololo import degtokpc
from scipy.spatial.distance import cdist

plt.ioff()

detection = 'f275w'
data_dir = '/home/ariel/Workspace/GASP/HST/Data/'
catalog_dir = '/home/ariel/Workspace/GASP/HST/Data/eric_catalogs/'

redshift_dict = {'JO201': 0.044631,
                 'JO204': 0.042372,
                 'JO206': 0.051089,
                 'JW100': 0.061891,
                 'JW39': 0.066319,
                 'JO175': 0.046750}

if detection == 'f606w':
    clump_catalog = Table.read(catalog_dir + 'HST_' + detection + '_0.04px_reg_mask_denoised_disk275_masked_3sigma_'
                               '0delta_only_trunks_matched_quantities.csv', format='ascii.ecsv')
    clump_catalog = clump_catalog[clump_catalog['tail_gal_flag'] == 0]
    # _all_gals_sel_morphology
elif detection == 'optical_only':
    clump_catalog = Table.read(catalog_dir + 'HST_f606w_0.04px_reg_mask_denoised_disk275_masked_3sigma_0delta_only_'
                               'trunks_match_optical_only_sel.csv', format='ascii.ecsv')
    clump_catalog = clump_catalog[clump_catalog['tail_gal_flag'] == 0]
elif detection == 'bad_beta':
    clump_catalog = Table.read(catalog_dir + '1pix_dilation_fluxes.csv')
else:
    clump_catalog = Table.read(catalog_dir + 'HST_' + detection + '_0.04px_reg_mask_SNR2_clumps_2.5sigma_5delta_'
                               'denoise_2sigma_0delta_all_gals_z_BPT_sel_morphology.csv', format='ascii.ecsv')
    resolved_catalog = Table.read(catalog_dir + 'resolved/HST_' + detection + '_0.04px_reg_mask_SNR2_clumps_2.5sigma_'
                                                                              '5delta_denoise_2sigma_0delta_all_gals_z_'
                                                                              'BPT_sel_resolved_morphology.csv',
                                  format='ascii.ecsv')
    parent_input = Table.read(catalog_dir + 'HST_f606w_0.04px_reg_mask_denoised_disk275_masked_3sigma_0delta_'
                                            'only_trunks_matched_quantities.csv', format='ascii.ecsv')
    parent_input_optical_only = Table.read(catalog_dir + 'HST_f606w_0.04px_reg_mask_denoised_disk275_masked_3sigma_'
                                                         '0delta_only_trunks_match_optical_only_sel.csv',
                                           format='ascii.ecsv')


match_catalog = Table.read(catalog_dir + 'HST_f606w_0.04px_reg_mask_denoised_disk275_masked_3sigma_'
                                         '0delta_only_trunks_matched_list.csv', format='ascii.ecsv')

if detection == 'halpha':
    f275w_match_catalog = Table.read(catalog_dir + 'HST_f275w_0.04px_reg_mask_SNR2_clumps_2.5sigma_5delta_denoise_'
                                                   '2sigma_0delta_z_sel_trunks_matched_list.csv', format='ascii.ecsv')

clump_catalog['clump_id'] = clump_catalog['id_clump']

if detection == 'optical_only':
    clump_catalog['clump_id'] = [clump_catalog['id_clump'][i] + 'optical_only' for i in range(len(clump_catalog))]

if detection == 'bad_beta':
    clump_catalog['gal'] = [clump_catalog['id_clump'][i].split('_')[0] for i in range(len(clump_catalog))]

clump_catalog['galaxy_redshift'] = [redshift_dict[galaxy] for galaxy in clump_catalog['gal']]

if detection != 'bad_beta':
    clump_catalog['n_det'] = [np.sum(stat_tools.int_to_bool_list(clump_catalog['sel_flag'][i]))
                              for i in range(len(clump_catalog))]
    clump_catalog['tail'] = clump_catalog['tail_gal_flag'] == 0
    clump_catalog['extraplanar'] = clump_catalog['tail_gal_flag'] == 1
    clump_catalog['disk'] = clump_catalog['tail_gal_flag'] == 2
else:
    clump_catalog['sel_flag'] = np.full(len(clump_catalog), 31)


if detection == 'f606w':

    clump_catalog['halpha_match'] = np.zeros(len(clump_catalog), dtype=bool)

    for i in range(len(clump_catalog)):

        complex_matches = match_catalog[match_catalog['id_complex'] == clump_catalog['clump_id'][i]]
        if np.any(complex_matches['filter'] == 'halpha'):
            clump_catalog['halpha_match'][i] = True


if (detection == 'halpha') | (detection == 'f275w'):

    clump_catalog['parent_id'] = np.full(len(clump_catalog), 'None', dtype='S16')
    clump_catalog['match_flag'] = np.zeros(len(clump_catalog), dtype=bool)
    clump_catalog['quantities_match_flag'] = np.zeros(len(clump_catalog), dtype=bool)
    clump_catalog['resolved_flag'] = np.zeros(len(clump_catalog), dtype=bool)
    clump_catalog['parent_axial_ratio'] = np.zeros(len(clump_catalog))
    clump_catalog['parent_distance'] = np.zeros(len(clump_catalog))
    clump_catalog['parent_radius'] = np.zeros(len(clump_catalog))
    clump_catalog['parent_iso_radius'] = np.zeros(len(clump_catalog))
    clump_catalog['parent_resolved_flag'] = np.zeros(len(clump_catalog))
    clump_catalog['parent_area'] = np.zeros(len(clump_catalog))
    clump_catalog['parent_displacement_ha'] = np.zeros(len(clump_catalog))
    clump_catalog['parent_displacement_uv'] = np.zeros(len(clump_catalog))
    clump_catalog['parent_f_uv'] = np.zeros(len(clump_catalog))
    clump_catalog['parent_f_ha'] = np.zeros(len(clump_catalog))
    clump_catalog['parent_f_opt'] = np.zeros(len(clump_catalog))
    clump_catalog['displacement_deg'] = np.zeros(len(clump_catalog))
    clump_catalog['displacement_kpc'] = np.zeros(len(clump_catalog))
    clump_catalog['displacement_deg_optical_only'] = np.zeros(len(clump_catalog))
    clump_catalog['displacement_kpc_optical_only'] = np.zeros(len(clump_catalog))

    for i in range(len(clump_catalog)):

        if clump_catalog['clump_id'][i] in resolved_catalog['id_clump']:
            clump_catalog['resolved_flag'][i] = True

        if clump_catalog['clump_id'][i] not in match_catalog['id_clump']:
            continue

        else:
            clump_catalog['match_flag'][i] = True

            match_index = np.argwhere(clump_catalog['clump_id'][i] == match_catalog['id_clump']).ravel()[0]
            clump_catalog['parent_id'][i] = match_catalog['id_complex'][match_index]

        if clump_catalog['parent_id'][i] not in parent_input['id_clump']:
            print(clump_catalog['parent_id'][i])
            continue

        else:
            clump_catalog['quantities_match_flag'][i] = True

            parent_catalog_index = np.argwhere(clump_catalog['parent_id'][i] == parent_input['id_clump']).ravel()[0]
            clump_catalog['parent_axial_ratio'][i] = parent_input['axial_ratio'][parent_catalog_index]
            clump_catalog['parent_distance'][i] = parent_input['dist_kpc'][parent_catalog_index]
            clump_catalog['parent_radius'][i] = parent_input['r_core_corr'][parent_catalog_index]
            clump_catalog['parent_iso_radius'][i] = parent_input['iso_radius(pc)'][parent_catalog_index]/1000
            clump_catalog['parent_area'][i] = parent_input['area_exact'][parent_catalog_index]
            clump_catalog['parent_f_uv'][i] = parent_input['dA_uv'][parent_catalog_index]
            clump_catalog['parent_f_ha'][i] = parent_input['dA_ha'][parent_catalog_index]
            clump_catalog['parent_resolved_flag'][i] = parent_input['resolved_flag'][parent_catalog_index]
            clump_catalog['parent_displacement_ha'][i] = parent_input['dist_ha'][parent_catalog_index]
            clump_catalog['parent_displacement_uv'][i] = parent_input['dist_uv'][parent_catalog_index]

            displacement_deg = np.sqrt((clump_catalog['x_cen'][i]-parent_input['x_cen'][parent_catalog_index]) ** 2 +
                                       (clump_catalog['y_cen'][i]-parent_input['y_cen'][parent_catalog_index]) ** 2)
            displacement_kpc = displacement_deg * degtokpc(clump_catalog['galaxy_redshift'][i])

            clump_catalog['displacement_deg'][i] = displacement_deg
            clump_catalog['displacement_kpc'][i] = displacement_kpc

            parent_catalog_index_optical_only = np.argwhere(clump_catalog['parent_id'][i] == parent_input_optical_only['id_clump']).ravel()[0]

            displacement_deg_optical_only = np.sqrt((clump_catalog['x_cen'][i]-parent_input_optical_only['x_cen'][parent_catalog_index_optical_only]) ** 2 +
                                                    (clump_catalog['y_cen'][i]-parent_input_optical_only['y_cen'][parent_catalog_index_optical_only]) ** 2)
            displacement_kpc_optical_only = displacement_deg_optical_only * degtokpc(clump_catalog['galaxy_redshift'][i])

            clump_catalog['displacement_deg_optical_only'][i] = displacement_deg_optical_only
            clump_catalog['displacement_kpc_optical_only'][i] = displacement_kpc_optical_only

            clump_catalog['parent_f_opt'][i] = parent_input_optical_only['area_exact'][parent_catalog_index_optical_only]\
                                               /parent_input['area_exact'][parent_catalog_index]


if detection == 'halpha':

    clump_catalog['f275w_parent_id'] = np.full(len(clump_catalog), 'None', dtype='S32')
    clump_catalog['f275w_match_flag'] = np.zeros(len(clump_catalog), dtype=bool)

    for i in range(len(clump_catalog)):

        if clump_catalog['clump_id'][i] not in f275w_match_catalog['id_clump']:
            continue

        else:
            clump_catalog['f275w_match_flag'][i] = True
            f275w_match_index = np.argwhere(clump_catalog['clump_id'][i] == f275w_match_catalog['id_clump']).ravel()[0]
            clump_catalog['f275w_parent_id'][i] = f275w_match_catalog['id_complex'][f275w_match_index]


clump_catalog['dist_to_disk_deg'] = np.zeros(len(clump_catalog))

for galaxy in ['JO201', 'JO204', 'JO206', 'JW100', 'JW39', 'JO175']:

    galaxy_flag = clump_catalog['gal'] == galaxy

    galaxy_contours = np.genfromtxt('/home/ariel/Workspace/GASP/HST/Data/contours/'+galaxy+'_f814w_contours_2sigma.dat')
    galaxy_coordinates = clump_catalog['x_cen', 'y_cen'][galaxy_flag].as_array()
    galaxy_coordinates_list = np.array([[galaxy_coordinates[i][0], galaxy_coordinates[i][1]] for i
                                        in range(len(galaxy_coordinates))])
    all_distances = cdist(galaxy_coordinates_list, galaxy_contours)
    min_distances = np.array([np.min(all_distances[i]) for i in range(len(all_distances))])

    clump_catalog['dist_to_disk_deg'][galaxy_flag] = min_distances

clump_catalog['dist_to_disk_kpc'] = np.array([clump_catalog['dist_to_disk_deg'][i] *
                                             degtokpc(clump_catalog['galaxy_redshift'][i]) for i in
                                             range(len(clump_catalog))])

if detection == 'optical_only':
    clump_catalog = clump_catalog[clump_catalog['F606W'] != 0]

print(clump_catalog)
clump_catalog.write(data_dir + 'full_catalogs/' + detection + '_bagpipes_input_all.fits', overwrite=True)

if detection not in ['f606w', 'optical_only']:
    clump_catalog = clump_catalog[(clump_catalog['sel_flag'] == 31) & (clump_catalog['id_catalog'] != 'M') &
                                  (~clump_catalog['disk']) &
                                  ((clump_catalog['level'] == 0) | (clump_catalog['leaf_flag'] == 1))]
else:
    clump_catalog = clump_catalog[clump_catalog['sel_flag'] == 31]

print(clump_catalog)
clump_catalog.write(data_dir + detection + '_bagpipes_input.fits', overwrite=True)



