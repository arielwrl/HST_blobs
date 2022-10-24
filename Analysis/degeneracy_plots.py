"""

ariel@oapd
18/08/22

Checks possible degeneracies using PDF information from BAGPIPES

"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack
import seaborn as sns
from toolbox import plot_tools
from astropy.io import fits

sns.set_style('ticks')

f606_image_raw = fits.open('/home/ariel/Workspace/GASP/HST/Data/images/JO201_f606w_0.04px_drc_MWdc.fits')
f606_image = np.ma.masked_array(f606_image_raw[1].data, mask=(f606_image_raw[1].data < 0.1e-21))

input = Table.read('/home/ariel/Workspace/GASP/HST/Data/disk_halpha_bagpipes_input.fits')
output = Table.read('/home/ariel/Workspace/GASP/HST/Data/disk_halpha_parametric_bagpipes_results.fits')

matched_complexes = Table.read('/home/ariel/Workspace/GASP/HST/Data/matched_tail_f606w.fits')

sns.jointplot(x=matched_complexes['age_1'].tolist(), y=matched_complexes['age_2'].tolist())

plt.ylabel(r'Age', fontsize=20)
plt.xlabel(r'Age (Optical Only)', fontsize=20)

x = np.linspace(-0, 0.5, 100)
y = x

plt.plot(x, y, '--k')

plt.show()


plt.figure()

sns.jointplot(x=matched_complexes['age_1'].tolist(), y=matched_complexes['age_2'].tolist())

plt.ylabel(r'MW Age', fontsize=20)
plt.xlabel(r'MW Age (Optical Only)', fontsize=20)

x = np.linspace(-0, 0.5, 100)
y = x

plt.plot(x, y, '--k')

plt.show()




# flag = (input['sel_flag'] == 31)
#
# input = input[flag]
# output = output[flag]
#
# mass_ratio = np.log10(output['formed_mass_old']**10/output['formed_mass_young']**10)
#
# plt.scatter(x=mass_ratio.tolist(), y=output['age_young_iqr'].tolist(), c=output['age_young'].tolist(),
#             cmap='plasma_r', edgecolors='white')
#
# cb = plt.colorbar()
# cb.set_label(r'$t_\star$ (Young)', fontsize=20)
#
# plt.xlabel(r'$\log \; M_{\mathrm{Old}}/M_{\mathrm{Young}}$', fontsize=20)
# plt.ylabel(r'$\mathrm{IQR}_{t_\star}$ (Young)', fontsize=20)
#
# sns.despine()
#
# plt.show()
#
#
# plt.figure()
#
# galaxy_flag = input['gal'] == 'JO201'
# lim_flag = mass_ratio < 3
#
# input = input[galaxy_flag & lim_flag]
# output = output[galaxy_flag & lim_flag]
# mass_ratio = mass_ratio[galaxy_flag & lim_flag]
#
# plt.imshow(np.log10(f606_image), cmap='Greys', origin='lower')
# plt.scatter(input['xc_pix(HST)'], input['yc_pix(HST)'], s=30,
#             c=mass_ratio, cmap='viridis_r')
#
# plt.colorbar()
#
# plt.show()
#
# plt.figure()
#
# flag = input['sel_flag'] == 31
#
# plt.scatter(output['formed_mass_old'][flag], output['formed_mass_young'][flag])
#
# plt.figure()
#
# galaxy_flag = input['gal'] == 'JO206'
#
# input = input[galaxy_flag]
# output = output[galaxy_flag]
#
# density = np.log10(output['formed_mass_old']**10/input['area_exact'])
#
# min_density = np.percentile(density, 20)
# max_density = np.percentile(density, 80)
#
# plt.imshow(np.log10(f606_image), cmap='Greys', origin='lower')
# plt.scatter(input['xc_pix(HST)'], input['yc_pix(HST)'], s=30,
#             c=density, cmap='viridis_r', vmin=min_density, vmax=max_density)
#
# plt.colorbar()
#
# plt.show()

# halpha_input = Table.read('/home/ariel/Workspace/GASP/HST/Data/tail_halpha_bagpipes_input.fits')
# f275w_input = Table.read('/home/ariel/Workspace/GASP/HST/Data/tail_f275w_bagpipes_input.fits')
# f606w_input = Table.read('/home/ariel/Workspace/GASP/HST/Data/tail_f606w_bagpipes_input.fits')
#
# tail_halpha = Table.read('/home/ariel/Workspace/GASP/HST/Data/tail_halpha_parametric_bagpipes_results.fits')
# tail_f275w = Table.read('/home/ariel/Workspace/GASP/HST/Data/tail_f275w_parametric_bagpipes_results.fits')
# tail_f606w = Table.read('/home/ariel/Workspace/GASP/HST/Data/tail_f606w_parametric_bagpipes_results.fits')
#
# tail_halpha['full_detection'] = np.zeros_like(tail_halpha['sel_flag'], dtype=bool)
# tail_halpha['full_detection'][tail_halpha['sel_flag'] == 31] = np.ones((tail_halpha['sel_flag'] == 31).sum(), dtype=bool)
#
# tail_f275w['full_detection'] = np.zeros_like(tail_f275w['sel_flag'], dtype=bool)
# tail_f275w['full_detection'][tail_f275w['sel_flag'] == 31] = np.ones((tail_f275w['sel_flag'] == 31).sum(), dtype=bool)
#
# tail_f606w['full_detection'] = np.zeros_like(tail_f606w['sel_flag'], dtype=bool)
# tail_f606w['full_detection'][tail_f606w['sel_flag'] == 31] = np.ones((tail_f606w['sel_flag'] == 31).sum(), dtype=bool)
#
# tail_halpha['blob_class'] = np.full_like(tail_halpha['galaxy'], fill_value=r'$H\alpha$')
# tail_f275w['blob_class'] = np.full_like(tail_f275w['galaxy'], fill_value='UV')
# tail_f606w['blob_class'] = np.full_like(tail_f606w['galaxy'], fill_value='SF Complexes')
#
# del tail_f606w['age_peaks', 'mwage_peaks', 'Av_peaks']
# tail_all = vstack([tail_halpha, tail_f275w, tail_f606w])
#
# blob_type = [r'$H\alpha$', 'UV', 'SF Complexes', r'$H\alpha$', 'UV', 'SF Complexes']
# n_peaks_perc = [100 * (tail_halpha['n_mwage_peaks'] == 1).sum() / len(tail_halpha),
#                 100 * (tail_f275w['n_mwage_peaks'] == 1).sum() / len(tail_f275w),
#                 100 * (tail_f606w['n_mwage_peaks'] == 1).sum() / len(tail_f606w),
#                 100 * (tail_halpha['n_mwage_peaks'][tail_halpha['full_detection']] == 1).sum() /
#                 tail_halpha['full_detection'].sum(),
#                 100 * (tail_f275w['n_mwage_peaks'][tail_f275w['full_detection']] == 1).sum() /
#                 tail_f275w['full_detection'].sum(),
#                 100 * (tail_f606w['n_mwage_peaks'][tail_f606w['full_detection']] == 1).sum() /
#                 tail_f606w['full_detection'].sum()
#                 ]
# full_detection = ['Not detected in all filters', 'Not detected in all filters', 'Not detected in all filters',
#                   'Detected in all filters', 'Detected in all filters', 'Detected in all filters']
#
# sns.barplot(x=blob_type, y=n_peaks_perc, hue=full_detection, palette='crest')
#
# plt.ylabel(r'\% of single-peaked $\langle\,t_\star\,\rangle_M$ PDFs', fontsize=20)
# plt.legend(frameon=False, fontsize=20)
# plt.gca().tick_params(axis='x', labelsize=20)
# plt.gca().tick_params(axis='y', labelsize=15)
# plt.ylim(0, 120)
# sns.despine()

# halpha_df = tail_halpha['galaxy', 'blob_id', 'age', 'tau', 'mwage', 'stellar_mass', 'n_mwage_peaks', 'Av',
#                         'detection'].to_pandas()
# f275w_df = tail_f275w['galaxy', 'blob_id', 'age', 'tau', 'mwage', 'stellar_mass', 'n_mwage_peaks', 'Av',
#                       'detection'].to_pandas()
# f606w_df = tail_f606w['galaxy', 'blob_id', 'age', 'tau', 'mwage', 'stellar_mass', 'n_mwage_peaks', 'Av',
#                       'detection'].to_pandas()
# all_df = tail_all['galaxy', 'blob_id', 'age', 'tau', 'mwage', 'stellar_mass', 'n_mwage_peaks', 'Av',
#                       'detection'].to_pandas()



