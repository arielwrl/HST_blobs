"""

ariel@oapd
November 14, 2022

Plots comparissons between models

"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import seaborn as sns
import astropy.visualization as vis
from astropy.wcs import WCS
import matplotlib.gridspec as gridspec
from scipy.ndimage.filters import gaussian_filter
from astropy.table import Table
from hst_pipeutils import plot_example
from toolbox import plot_tools
import gc

sns.set_style('ticks')

input_catalog = Table.read('/home/ariel/Workspace/GASP/HST/Data/tail_f606w_bagpipes_input.fits')

output_ssp = Table.read('/home/ariel/Workspace/GASP/HST/Data/tail_f606w_ssp_bagpipes_results.fits')
output_constant = Table.read('/home/ariel/Workspace/GASP/HST/Data/tail_f606w_constant_bagpipes_results.fits')
output_dlogprior = Table.read('/home/ariel/Workspace/GASP/HST/Data/tail_f606w_dexp_logprior_bagpipes_results.fits')

flag = input_catalog['sel_flag'] == 31

input_catalog = input_catalog[flag]

skip = 700

for i in range(len(input_catalog)):

    fig, ax_list = plt.subplots(2, 2, figsize=(15, 10), sharex=True, sharey=True)

    plot_example(output_ssp['fit_id_fitconfig'][i], input_catalog, ax=ax_list[0, 0], config='ssp',
                 show_xlabel=False, show_ylabel=True)
    plot_example(output_constant['fit_id_fitconfig'][i], input_catalog, ax=ax_list[0, 1], config='constant',
                 show_xlabel=False, show_ylabel=False)
    plot_example(output_dlogprior['fit_id_fitconfig'][i], input_catalog, ax=ax_list[1, 0], config='dexp_logprior',
                 show_xlabel=True, show_ylabel=True)

    ax_list[0, 0].annotate(r'SSP, $\chi^2/N=%0.3f$' % (output_ssp['chi2'][i]/5), fontsize=20, xy=(0.01, 0.03),
                           xycoords='axes fraction')
    ax_list[0, 1].annotate(r'Constant, $\chi^2/N=%0.3f$' % (output_constant['chi2'][i]/5), fontsize=20, xy=(0.01, 0.03),
                           xycoords='axes fraction')
    ax_list[1, 0].annotate(r'Delayed (log prior), $\chi^2/N=%0.3f$' % (output_dlogprior['chi2'][i]/5), fontsize=20,
                           xy=(0.01, 0.03), xycoords='axes fraction')

    fig.subplots_adjust(wspace=0, hspace=0, left=0.055, bottom=0.08, right=0.98, top=0.95)

    if input_catalog['tail_gal_flag'][i] == 0:
        id_label = input_catalog['blob_id'][i].split('_')
        plt.suptitle(id_label[0] + '\_' + id_label[1] + ', Tail', fontsize=20)
    if input_catalog['tail_gal_flag'][i] == 1:
        id_label = input_catalog['blob_id'][i].split('_')
        plt.suptitle(id_label[0] + '\_' + id_label[1] + ', Extraplanar', fontsize=20)

    plt.savefig('/home/ariel/Workspace/GASP/HST/BAGPIPES/plots/comparison_f606w/' + input_catalog['fit_id'][i] + '.png')

    plt.close('all')

    gc.collect()


    # if i < skip:
    #     continue
    #
    # print('Plotting clump ', i, 'of ', len(input_catalog))
    #
    # fig, ax = plt.subplots(1, 1, figsize=(8, 5), sharex=True, sharey=True)
    #
    # plot_example(output_dlogprior_hacked['fit_id_fitconfig'][i], input_catalog, ax=ax, config='dexp_logprior_hacked',
    #              show_xlabel=True, show_ylabel=True)
    # ax.annotate(r'Hacked, $\chi^2/N=%0.3f$' % (output_dlogprior_hacked['chi2'][i]/5), fontsize=20, xy=(0.01, 0.03),
    #                        xycoords='axes fraction')
    #
    # fig.subplots_adjust(wspace=0, hspace=0, left=0.08, bottom=0.13, right=0.98, top=0.925)
    #
    # if input_catalog['tail_gal_flag'][i] == 0:
    #     id_label = input_catalog['blob_id'][i].split('_')
    #     plt.suptitle(id_label[0] + '\_' + id_label[1] + ', Tail', fontsize=20)
    # if input_catalog['tail_gal_flag'][i] == 1:
    #     id_label = input_catalog['blob_id'][i].split('_')
    #     plt.suptitle(id_label[0] + '\_' + id_label[1] + ', Extraplanar', fontsize=20)
    #
    # plt.savefig('/home/ariel/Workspace/GASP/HST/BAGPIPES/plots/hacked_fits/' + input_catalog['fit_id'][i] + '.png')
    #
    # plt.close('all')
    #
    # gc.collect()
