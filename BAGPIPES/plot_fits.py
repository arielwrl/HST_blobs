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

input_catalog = Table.read('/home/ariel/Workspace/GASP/HST/Data/f275w_bagpipes_input.fits')
output_catalog = Table.read('/home/ariel/Workspace/GASP/HST/Data/f275w_dexp_logprior_single_bagpipes_results.fits')

output_catalog = output_catalog[(~input_catalog['disk']) & ((input_catalog['level'] == 0) |
                                                            (input_catalog['leaf_flag'] == 1))]

input_catalog = input_catalog[(~input_catalog['disk']) & ((input_catalog['level'] == 0) |
                                                          (input_catalog['leaf_flag'] == 1))]

input_catalog['fit_id'] = output_catalog['fit_id']

for i in range(len(input_catalog)):

    print('Plotting clump ', i, 'of ', len(input_catalog))

    fig, ax = plt.subplots(1, 1, figsize=(12.5, 7), sharex=True, sharey=True)

    plot_example(output_catalog['fit_id'][i], input_catalog, ax=ax, show_xlabel=True, show_ylabel=True,
                 spec_color='royalblue')
    ax.annotate(r'$\chi^2/N=%0.3f$' % (output_catalog['chi2'][i]/5), fontsize=20, xy=(0.01, 0.03),
                xycoords='axes fraction')

    fig.subplots_adjust(wspace=0, hspace=0, left=0.08, bottom=0.13, right=0.98, top=0.925)

    if input_catalog['tail_gal_flag'][i] == 0:
        id_label = input_catalog['clump_id'][i].split('_')
        plt.suptitle(id_label[0] + '\_' + id_label[1] + ', Tail', fontsize=20)
    if input_catalog['tail_gal_flag'][i] == 1:
        id_label = input_catalog['clump_id'][i].split('_')
        plt.suptitle(id_label[0] + '\_' + id_label[1] + ', Extraplanar', fontsize=20)

    plt.savefig('/home/ariel/Workspace/GASP/HST/BAGPIPES/plots/f275w_single_allZs/' + input_catalog['fit_id'][i]
                + '.png')

    plt.close('all')

    gc.collect()
