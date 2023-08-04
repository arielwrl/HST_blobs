"""

ariel@oapd
16/01/2023

Includes properties of parent and sibling clumps in the output tables, flags out unfortunate objects

"""

import numpy as np
from astropy.table import Table

detection = 'halpha'
config = 'dexp_logprior_single'

data_dir = '/home/ariel/Workspace/GASP/HST/Data/'

input_table = Table.read(data_dir + detection + '_bagpipes_input.fits')
output_table = Table.read(data_dir + detection + '_' + config + '_bagpipes_results.fits')
complex_input = Table.read(data_dir + 'f606w_bagpipes_input.fits')
complex_output = Table.read(data_dir + 'f606w_dexp_logprior_single_bagpipes_results.fits')
optical_only_output = Table.read(data_dir + 'optical_only_dexp_logprior_single_bagpipes_results.fits')

if detection == 'halpha':
    f275w_output = Table.read(data_dir + 'f275w_dexp_logprior_single_bagpipes_results.fits')

output_table['F275W_F336W_color_syn'] = -2.5*np.log10(output_table['photometry_syn'][:, 0]/
                                                      output_table['photometry_syn'][:, 1])
output_table['F275W_F336W_color_obs'] = -2.5*np.log10(output_table['photometry_obs'][:, 0]/
                                                      output_table['photometry_obs'][:, 1])

output_table['F275W_F606W_color_syn'] = -2.5*np.log10(output_table['photometry_syn'][:, 0]/
                                                      output_table['photometry_syn'][:, 2])
output_table['F275W_F606W_color_obs'] = -2.5*np.log10(output_table['photometry_obs'][:, 0]/
                                                      output_table['photometry_obs'][:, 2])


filter_list = ['F275W', 'F336W', 'F606W', 'F680N', 'F814W']

output_table['photometry_errors'] = np.zeros_like(output_table['photometry_obs'])

for i in range(len(output_table)):
    output_table['photometry_errors'][i] = np.array([2 * input_table['err' + filter_name][i] for filter_name in filter_list])

output_table['residuals'] = (output_table['photometry_obs']-output_table['photometry_syn'])/output_table['photometry_errors']
output_table['residuals_within_iqr'] = (output_table['photometry_obs']-output_table['photometry_syn'])/(output_table['photometry_errors'] + output_table['photometry_syn_iqr'])
output_table['absolute_residuals'] = np.absolute(output_table['residuals'])
output_table['absolute_residuals_within_iqr'] = np.absolute(output_table['residuals_within_iqr'])

output_table['n_bad_points'] = np.zeros(len(output_table))
output_table['n_bad_points_within_iqr'] = np.zeros(len(output_table))

for i in range(len(output_table)):
    output_table['n_bad_points'][i] = (output_table['residuals'][i] > 1).sum()
    output_table['n_bad_points_within_iqr'][i] = (output_table['residuals_within_iqr'][i] > 1).sum()

output_table['outside_errorbar'] = output_table['absolute_residuals'] > 1
output_table['outside_errorbar_soft'] = output_table['absolute_residuals_within_iqr'] > 1

output_table['bad_fit_npoints'] = output_table['n_bad_points'] > 1
output_table['bad_fit_within_iqr'] = output_table['n_bad_points_within_iqr'] > 0
output_table['bad_fit'] = output_table['bad_fit_npoints'] | output_table['bad_fit_within_iqr']

if detection != 'f606w':
    output_table['parent_id'] = input_table['parent_id']
    output_table['match_flag_output'] = np.ones_like(output_table['mwage'], dtype=bool)
    output_table['match_flag_output_optical_only'] = np.ones_like(output_table['mwage'], dtype=bool)

    output_table['parent_mwage'] = np.zeros_like(output_table['mwage'])
    output_table['parent_stellar_mass'] = np.zeros_like(output_table['mwage'])
    output_table['parent_age'] = np.zeros_like(output_table['mwage'])
    output_table['parent_sfr'] = np.zeros_like(output_table['mwage'])
    output_table['parent_age_std'] = np.zeros_like(output_table['mwage'])
    output_table['parent_tau'] = np.zeros_like(output_table['mwage'])
    output_table['parent_Av'] = np.zeros_like(output_table['mwage'])
    output_table['parent_eta'] = np.zeros_like(output_table['mwage'])
    output_table['parent_Ha'] = np.zeros_like(output_table['mwage'])
    output_table['parent_mwage_optical_only'] = np.zeros_like(output_table['mwage'])
    output_table['parent_stellar_mass_optical_only'] = np.zeros_like(output_table['mwage'])
    output_table['parent_age_optical_only'] = np.zeros_like(output_table['mwage'])
    output_table['parent_age_std_optical_only'] = np.zeros_like(output_table['mwage'])
    output_table['parent_tau_optical_only'] = np.zeros_like(output_table['mwage'])
    output_table['parent_Av_optical_only'] = np.zeros_like(output_table['mwage'])
    output_table['parent_eta_optical_only'] = np.zeros_like(output_table['mwage'])
    output_table['parent_Ha_optical_only'] = np.zeros_like(output_table['mwage'])

    for i in range(len(input_table)):

        if input_table['parent_id'][i] not in complex_output['clump_id']:
            output_table['match_flag_output'][i] = False
            continue

        else:
            match_index = np.argwhere(complex_output['clump_id'] == input_table['parent_id'][i]).ravel()[0]

            output_table['parent_mwage'][i] = complex_output['mwage'][match_index]
            output_table['parent_stellar_mass'][i] = complex_output['stellar_mass'][match_index]
            output_table['parent_age'][i] = complex_output['age'][match_index]
            output_table['parent_sfr'][i] = complex_output['sfr'][match_index]
            output_table['parent_age_std'][i] = complex_output['age_std'][match_index]
            output_table['parent_tau'][i] = complex_output['tau'][match_index]
            output_table['parent_Av'][i] = complex_output['Av'][match_index]
            output_table['parent_eta'][i] = complex_output['eta'][match_index]
            output_table['parent_Ha'][i] = complex_output['Ha'][match_index]


    # Optical only match

    optical_only_transformed_ids = np.array([optical_only_output['clump_id'][i].split('_')[0] + '_' +
                                             optical_only_output['clump_id'][i].split('_')[1] + '_f606w' for i in
                                             range(len(optical_only_output))])

    for i in range(len(input_table)):

        if input_table['parent_id'][i] not in optical_only_transformed_ids:
            output_table['match_flag_output_optical_only'][i] = False
            continue

        else:
            match_index = np.argwhere(optical_only_transformed_ids == input_table['parent_id'][i]).ravel()[0]

            output_table['parent_mwage_optical_only'][i] = optical_only_output['mwage'][match_index]
            output_table['parent_stellar_mass_optical_only'][i] = optical_only_output['mwage'][match_index]
            output_table['parent_age_optical_only'][i] = optical_only_output['age'][match_index]
            output_table['parent_age_std_optical_only'][i] = optical_only_output['age_std'][match_index]
            output_table['parent_tau_optical_only'][i] = optical_only_output['tau'][match_index]
            output_table['parent_Av_optical_only'][i] = optical_only_output['Av'][match_index]
            output_table['parent_eta_optical_only'][i] = optical_only_output['eta'][match_index]
            output_table['parent_Ha_optical_only'][i] = optical_only_output['Ha'][match_index]

    # Matching siblings:

    output_table['n_siblings'] = np.zeros_like(output_table['mwage'])

    output_table['youngest_sibling'] = np.zeros_like(output_table['mwage'], dtype=bool)
    output_table['mean_siblings_mwage'] = np.zeros_like(output_table['mwage'])
    output_table['min_siblings_mwage'] = np.zeros_like(output_table['mwage'])
    output_table['mean_siblings_stellar_mass'] = np.zeros_like(output_table['mwage'])
    output_table['min_siblings_stellar_mass'] = np.zeros_like(output_table['mwage'])
    output_table['mean_siblings_age'] = np.zeros_like(output_table['age'])
    output_table['min_siblings_age'] = np.zeros_like(output_table['age'])
    output_table['min_siblings_age_std'] = np.zeros_like(output_table['age'])

    for i in range(len(complex_output)):

        if complex_output['clump_id'][i] not in output_table['parent_id']:
            continue

        complex_flag = output_table['parent_id'] == complex_output['clump_id'][i]
        output_table['n_siblings'][complex_flag] = np.full(complex_flag.sum(), complex_flag.sum())
        output_table['mean_siblings_mwage'][complex_flag] = np.full(complex_flag.sum(),
                                                                    np.mean(output_table['mwage'][complex_flag]))
        output_table['min_siblings_mwage'][complex_flag] = np.full(complex_flag.sum(),
                                                                   np.min(output_table['mwage'][complex_flag]))
        output_table['mean_siblings_stellar_mass'][complex_flag] = np.full(complex_flag.sum(),
                                                                    np.mean(output_table['stellar_mass'][complex_flag]))
        output_table['min_siblings_stellar_mass'][complex_flag] = np.full(complex_flag.sum(),
                                                                   np.min(output_table['stellar_mass'][complex_flag]))
        output_table['mean_siblings_age'][complex_flag] = np.full(complex_flag.sum(),
                                                                    np.mean(output_table['age'][complex_flag]))
        output_table['min_siblings_age'][complex_flag] = np.full(complex_flag.sum(),
                                                                   np.min(output_table['age'][complex_flag]))
        output_table['min_siblings_age_std'][complex_flag] = np.full(complex_flag.sum(),
                                                                     output_table['age_std'][complex_flag][np.argmin(output_table['age'][complex_flag]).ravel()[0]])

    for i in range(len(output_table)):
        if output_table['min_siblings_mwage'][i] == output_table['mwage'][i]:
            output_table['youngest_sibling'][i] = True



    # F275W match

    if detection == 'halpha':

        output_table['parent_id_f275w'] = input_table['f275w_parent_id']
        output_table['match_flag_output_f275w'] = np.ones_like(output_table['mwage'], dtype=bool)
        output_table['parent_mwage_f275w'] = np.zeros_like(output_table['mwage'])
        output_table['parent_age_f275w'] = np.zeros_like(output_table['mwage'])
        output_table['parent_age_std_f275w'] = np.zeros_like(output_table['mwage'])
        output_table['parent_tau_f275w'] = np.zeros_like(output_table['mwage'])
        output_table['parent_Av_f275w'] = np.zeros_like(output_table['mwage'])
        output_table['parent_eta_f275w'] = np.zeros_like(output_table['mwage'])

        for i in range(len(input_table)):

            if input_table['f275w_parent_id'][i] not in f275w_output['clump_id']:
                output_table['match_flag_output_f275w'][i] = False
                continue

            else:
                match_index = np.argwhere(f275w_output['clump_id'] == input_table['f275w_parent_id'][i]).ravel()[0]

                output_table['parent_mwage_f275w'][i] = f275w_output['mwage'][match_index]
                output_table['parent_age_f275w'][i] = f275w_output['age'][match_index]
                output_table['parent_age_std_f275w'][i] = f275w_output['age_std'][match_index]
                output_table['parent_tau_f275w'][i] = f275w_output['tau'][match_index]
                output_table['parent_Av_f275w'][i] = f275w_output['Av'][match_index]
                output_table['parent_eta_f275w'][i] = f275w_output['eta'][match_index]

        output_table['n_siblings_f275w'] = np.zeros_like(output_table['mwage'])
        output_table['mean_siblings_mwage_f275w'] = np.zeros_like(output_table['mwage'])
        output_table['min_siblings_mwage_f275w'] = np.zeros_like(output_table['mwage'])
        output_table['mean_siblings_age_f275w'] = np.zeros_like(output_table['age'])
        output_table['min_siblings_age_f275w'] = np.zeros_like(output_table['age'])
        output_table['min_siblings_age_std_f275w'] = np.zeros_like(output_table['age'])

        for i in range(len(f275w_output)):

            if f275w_output['clump_id'][i] not in output_table['parent_id_f275w']:
                continue

            f275w_flag = output_table['parent_id_f275w'] == f275w_output['clump_id'][i]
            output_table['n_siblings_f275w'][f275w_flag] = np.full(f275w_flag.sum(), f275w_flag.sum())
            output_table['mean_siblings_mwage_f275w'][f275w_flag] = np.full(f275w_flag.sum(),
                                                                            np.mean(output_table['mwage'][f275w_flag]))
            output_table['min_siblings_mwage_f275w'][f275w_flag] = np.full(f275w_flag.sum(),
                                                                           np.min(output_table['mwage'][f275w_flag]))
            output_table['mean_siblings_age_f275w'][f275w_flag] = np.full(f275w_flag.sum(),
                                                                          np.mean(output_table['age'][f275w_flag]))
            output_table['min_siblings_age_f275w'][f275w_flag] = np.full(f275w_flag.sum(),
                                                                          np.min(output_table['age'][f275w_flag]))
            output_table['min_siblings_age_std_f275w'][f275w_flag] = np.full(f275w_flag.sum(),
                                                                             output_table['age_std'][f275w_flag][np.argmin(output_table['age'][f275w_flag]).ravel()[0]])

    if config.split('_')[-1] == 'single':
        double_output = Table.read(data_dir + detection + '_dexp_logprior_double_bagpipes_results.fits')
        output_table['bad_double_fit'] = double_output['mass_ratio'] > 2

output_table.write(data_dir + detection + '_' + config + '_bagpipes_results.fits', overwrite=True)
