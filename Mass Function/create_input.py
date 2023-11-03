"""

ariel@oapd
11/01/2023

Converts Eric's catalogs into fits, adds galaxy redshift and clump ID, selects only 31 clumps


"""

from astropy.table import Table

detection = 'halpha'
data_dir = '/home/ariel/Workspace/GASP/HST/Mass Function/Data/'
catalog_dir = '/home/ariel/Workspace/GASP/HST/Mass Function/Data/eric_catalogs/'

redshift_dict = {'JO201': 0.044631,
                 'JO204': 0.042372,
                 'JO206': 0.051089,
                 'JW100': 0.061891,
                 'JW39': 0.066319,
                 'JO175': 0.046750}

clump_catalog = Table.read(catalog_dir + 'HST_' + detection + '_0.04px_reg_mask_no_nan_tail_simulated_clumps_full'
                                                                  '_flux_fluxes_SNR2.csv', format='ascii.ecsv')

clump_catalog['clump_id'] = clump_catalog['id_clump_iter']

if detection == 'optical_only':
    clump_catalog['clump_id'] = [clump_catalog['id_clump'][i] + 'optical_only' for i in range(len(clump_catalog))]

clump_catalog['galaxy_redshift'] = [redshift_dict[galaxy] for galaxy in clump_catalog['gal']]

print(clump_catalog)
clump_catalog.write(data_dir + 'full_catalogs/' + detection + '_bagpipes_input_all.fits', overwrite=True)

if detection not in ['f606w', 'optical_only']:
    clump_catalog = clump_catalog[(clump_catalog['sel_flag'] == 31) &
                                  ((clump_catalog['level'] == 0) | (clump_catalog['leaf_flag'] == 1))]
else:
    clump_catalog = clump_catalog[clump_catalog['sel_flag'] == 31]

print(clump_catalog)
clump_catalog.write(data_dir + detection + '_bagpipes_input.fits', overwrite=True)



