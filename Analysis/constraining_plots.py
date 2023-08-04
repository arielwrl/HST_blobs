"""

ariel@oapd
30/01/2023

Plot histograms of the constraining parameter for different BAGPIPES measurements

Fit limits:

dexp = {}
dexp["age"] = (0.0, 0.5)
dexp["tau"] = (0.0001, 0.5)
dexp["massformed"] = (0., 10.)
dexp["metallicity"] = (0.005, 2.5)
dust["Av"] = (0, 1.25)
dust["eta"] = (1, 2.5)

"""

from astropy.table import Table
import matplotlib.pyplot as plt
import seaborn as sns
from collections import OrderedDict

sns.set_style('ticks')
palette = sns.color_palette("hls", 3)

input_halpha = Table.read('/home/ariel/Workspace/GASP/HST/Data/halpha_bagpipes_input.fits')
output_halpha_single = Table.read('/home/ariel/Workspace/GASP/HST/Data/halpha_dexp_logprior'
                                  '_single_bagpipes_results.fits')
output_halpha_double = Table.read('/home/ariel/Workspace/GASP/HST/Data/halpha_dexp_logprior'
                                  '_double_bagpipes_results.fits')

tail_halpha = input_halpha['tail_gal_flag'] == 0
extraplanar_halpha = input_halpha['tail_gal_flag'] == 1
disk_halpha = input_halpha['tail_gal_flag'] == 2

# Define constraining parameters

constraining_parameters_single = OrderedDict()
constraining_parameters_single['age'] = 2*output_halpha_single['age_iqr']/0.5
constraining_parameters_single['tau'] = 2*output_halpha_single['tau_iqr']/0.4999
constraining_parameters_single['formed_mass'] = 2*output_halpha_single['formed_mass_iqr']/10
constraining_parameters_single['metallicity'] = 2*output_halpha_single['metallicity_iqr']/2.495
constraining_parameters_single['Av'] = 2*output_halpha_single['Av_iqr']/1.25
constraining_parameters_single['eta'] = 2*output_halpha_single['eta_iqr']/1.5

constraining_parameters_double = OrderedDict()
constraining_parameters_double['age'] = 2*output_halpha_double['age_iqr']/0.5
constraining_parameters_double['tau'] = 2*output_halpha_double['tau_iqr']/0.4999
constraining_parameters_double['formed_mass'] = 2*output_halpha_double['formed_mass_iqr']/10
constraining_parameters_double['metallicity'] = 2*output_halpha_double['metallicity_iqr']/2.495
constraining_parameters_double['Av'] = 2*output_halpha_double['Av_iqr']/1.25
constraining_parameters_double['eta'] = 2*output_halpha_double['eta_iqr']/1.5


# # Histograms:
#
# fig, ax_list = plt.subplots(2, 3, figsize=(11.5, 8))
#
# for i in range(len(constraining_parameters_single.keys())):
#
#     ax = ax_list.ravel()[i]
#     key = list(constraining_parameters_single.keys())[i]
#
#     ax.hist(constraining_parameters_single[key][tail_halpha].tolist(), density=True, histtype='step', range=[0, 1.2],
#             color=palette[0], bins=15, label='Tail')
#     ax.hist(constraining_parameters_single[key][extraplanar_halpha].tolist(), density=True, histtype='step', range=[0, 1.2],
#             color=palette[1], bins=15, label='Extraplanar')
#     ax.hist(constraining_parameters_single[key][disk_halpha].tolist(), density=True, histtype='step', range=[0, 1.2],
#             color=palette[2], bins=15, label='Disk')
#
#     ax.set_xlabel(key, fontsize=20)
#
#     fig.tight_layout()
#
# ax_list.ravel()[0].legend()
#
# plt.savefig('/home/ariel/Workspace/GASP/HST/Analysis/Plots/constraining_halpha_single.png')
#
#
# fig, ax_list = plt.subplots(2, 3, figsize=(11.5, 8))
#
# for i in range(len(constraining_parameters_double.keys())):
#
#     ax = ax_list.ravel()[i]
#     key = list(constraining_parameters_double.keys())[i]
#
#     ax.hist(constraining_parameters_double[key][tail_halpha].tolist(), density=True, histtype='step', range=[0, 1.2],
#             color=palette[0], bins=15, label='Tail')
#     ax.hist(constraining_parameters_double[key][extraplanar_halpha].tolist(), density=True, histtype='step', range=[0, 1.2],
#             color=palette[1], bins=15, label='Extraplanar')
#     ax.hist(constraining_parameters_double[key][disk_halpha].tolist(), density=True, histtype='step', range=[0, 1.2],
#             color=palette[2], bins=15, label='Disk')
#
#     ax.set_xlabel(key, fontsize=20)
#
#     fig.tight_layout()
#
# ax_list.ravel()[0].legend()
#
# plt.savefig('/home/ariel/Workspace/GASP/HST/Analysis/Plots/constraining_halpha_double.png')
#

output_halpha_double = output_halpha_double[disk_halpha]

palette = sns.color_palette("hls", 3)

joint = sns.jointplot(output_halpha_double['mass_ratio'].tolist(), constraining_parameters_double['age'][disk_halpha].tolist(), palette=palette)

joint.fig.set_size_inches(11.5, 8)

joint.ax_joint.set_xlabel('Mass Ratio', fontsize=20)
joint.ax_joint.set_ylabel('Age constraining index', fontsize=20)
joint.ax_joint.set_xlim(-4.2, 6)
joint.ax_joint.set_ylim(0, 1.2)

joint.fig.tight_layout()

plt.savefig('/home/ariel/Workspace/GASP/HST/Analysis/Plots/mass_ratio_age_const.png')
