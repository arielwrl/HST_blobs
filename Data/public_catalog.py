"""

ariel@oapd
11/12/2023

Organizes data in a public catalog

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from astropy.table import Table
import seaborn as sns
from toolbox.wololo import redshift2lumdistance, arcsectokpc
import astropy.units as u
from starlight_toolkit.dust import CAL
from matplotlib.patches import Ellipse

halpha_input = Table.read('C:/Users/ariel/Workspace/GASP/HST/Data/halpha_bagpipes_input.fits')
f275w_input = Table.read('C:/Users/ariel/Workspace/GASP/HST/Data/f275w_bagpipes_input.fits')
f606w_input = Table.read('C:/Users/ariel/Workspace/GASP/HST/Data/f606w_bagpipes_input.fits')

output_halpha = Table.read('C:/Users/ariel/Workspace/GASP/HST/Data/halpha_dexp_logprior_single_bagpipes_results.fits')
output_f275w = Table.read('C:/Users/ariel/Workspace/GASP/HST/Data/f275w_dexp_logprior_single_bagpipes_results.fits')
output_f606w = Table.read('C:/Users/ariel/Workspace/GASP/HST/Data/f606w_dexp_logprior_single_bagpipes_results.fits')

halpha_input = halpha_input[~output_halpha['bad_double_fit'] & ~output_halpha['bad_fit']]
output_halpha = output_halpha[~output_halpha['bad_double_fit'] & ~output_halpha['bad_fit']]

f275w_input = f275w_input[~output_f275w['bad_double_fit'] & ~output_f275w['bad_fit']]
output_f275w = output_f275w[~output_f275w['bad_double_fit'] & ~output_f275w['bad_fit']]

f606w_input = f606w_input[~output_f606w['bad_fit']]
output_f606w = output_f606w[~output_f606w['bad_fit']]

output_halpha['detection'] = np.full_like(output_halpha['clump_id'], 'halpha')
output_f275w['detection'] = np.full_like(output_f275w['clump_id'], 'f275w')
output_f606w['detection'] = np.full_like(output_f606w['clump_id'], 'f606w')

output_halpha['sfh_class'] = np.full_like(output_halpha['clump_id'], 'Intermediate')
output_f275w['sfh_class'] = np.full_like(output_f275w['clump_id'], 'Intermediate')
output_f606w['sfh_class'] = np.full_like(output_f606w['clump_id'], 'Intermediate')

output_halpha['sfh_class'][output_halpha['early_decliner']] = 'Early Decliner'
output_f275w['sfh_class'][output_f275w['early_decliner']] = 'Early Decliner'
output_f606w['sfh_class'][output_f606w['early_decliner']] = 'Early Decliner'

output_halpha['sfh_class'][output_halpha['late_decliner']] = 'Late Decliner'
output_f275w['sfh_class'][output_f275w['late_decliner']] = 'Late Decliner'
output_f606w['sfh_class'][output_f606w['late_decliner']] = 'Late Decliner'

output_halpha['galaxy'] = output_halpha['galaxy'].astype(str)
output_f275w['galaxy'] = output_f275w['galaxy'].astype(str)
output_f606w['galaxy'] = output_f606w['galaxy'].astype(str)

output_halpha['location'] = np.zeros_like(output_halpha['galaxy'])
output_halpha['location'][halpha_input['tail_gal_flag'] == 0] = np.full((halpha_input['tail_gal_flag'] == 0).sum(),
                                                                        'Tail')
output_halpha['location'][halpha_input['tail_gal_flag'] == 1] = np.full((halpha_input['tail_gal_flag'] == 1).sum(),
                                                                        'Extraplanar')

output_f275w['location'] = np.zeros_like(output_f275w['galaxy'])
output_f275w['location'][f275w_input['tail_gal_flag'] == 0] = np.full((f275w_input['tail_gal_flag'] == 0).sum(),
                                                                      'Tail')
output_f275w['location'][f275w_input['tail_gal_flag'] == 1] = np.full((f275w_input['tail_gal_flag'] == 1).sum(),
                                                                      'Extraplanar')

output_f606w['location'] = np.zeros_like(output_f606w['galaxy'])
output_f606w['location'][f606w_input['tail_gal_flag'] == 0] = np.full((f606w_input['tail_gal_flag'] == 0).sum(),
                                                                      'Tail')
output_f606w['location'][f606w_input['tail_gal_flag'] == 1] = np.full((f606w_input['tail_gal_flag'] == 1).sum(),
                                                                      'Extraplanar')

output_halpha['Av'] = output_halpha['Av'] * output_halpha['eta']
output_f275w['Av'] = output_f275w['Av'] * output_f275w['eta']
output_f606w['Av'] = output_f606w['Av'] * output_f606w['eta']



public_catalog = Table()
public_catalog['id'] = np.concatenate([output_halpha['clump_id'].astype(str), output_f275w['clump_id'].astype(str), output_f606w['clump_id'].astype(str)])
public_catalog['RA'] = np.concatenate([halpha_input['x_cen'], f275w_input['x_cen'], f606w_input['x_cen']])
public_catalog['Dec'] = np.concatenate([halpha_input['y_cen'], f275w_input['y_cen'], f606w_input['y_cen']])
public_catalog['galaxy'] = np.concatenate([output_halpha['galaxy'], output_f275w['galaxy'], output_f606w['galaxy']])
public_catalog['galaxy_redshift'] = np.concatenate([halpha_input['galaxy_redshift'], f275w_input['galaxy_redshift'], f606w_input['galaxy_redshift']])
public_catalog['detection'] = np.concatenate([output_halpha['detection'].astype(str), output_f275w['detection'].astype(str), output_f606w['detection'].astype(str)])
public_catalog['location'] = np.concatenate([output_halpha['location'], output_f275w['location'], output_f606w['location']])
public_catalog['resolved_flag'] = np.concatenate([halpha_input['resolved_flag'], f275w_input['resolved_flag'], f606w_input['resolved_flag']])
public_catalog['mwage'] = np.concatenate([1e3 * output_halpha['mwage'].round(decimals=4), 1e3 * output_f275w['mwage'].round(decimals=4), 1e3 * output_f606w['mwage'].round(decimals=4)])
public_catalog['stellar_mass'] = np.concatenate([output_halpha['stellar_mass'].round(decimals=4), output_f275w['stellar_mass'].round(decimals=4), output_f606w['stellar_mass'].round(decimals=4)])
public_catalog['sigma_m'] = np.concatenate([output_halpha['sigma_m'].round(decimals=4), output_f275w['sigma_m'].round(decimals=4), output_f606w['sigma_m'].round(decimals=4)])
public_catalog['sfr'] = np.concatenate([np.log10(output_halpha['sfr']).round(decimals=4), np.log10(output_f275w['sfr']).round(decimals=4), np.log10(output_f606w['sfr']).round(decimals=4)])
public_catalog['Av'] = np.concatenate([output_halpha['Av'].round(decimals=4), output_f275w['Av'].round(decimals=4), output_f606w['Av'].round(decimals=4)])
public_catalog['sfh_class'] = np.concatenate([output_halpha['sfh_class'].astype(str), output_f275w['sfh_class'].astype(str), output_f606w['sfh_class'].astype(str)])

public_catalog.write('C:/Users/ariel/Workspace/GASP/HST/Data/public_catalog.dat', format='ascii', overwrite=True)