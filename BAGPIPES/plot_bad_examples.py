"""

ariel@cordenons
02/09/2023

"""

import matplotlib.pyplot as plt
import seaborn as sns
from astropy.table import Table
from hst_pipeutils import plot_example
from toolbox import plot_tools

sns.set_style('ticks')

input_catalog = Table.read('/home/ariel/Workspace/GASP/HST/Data/f275w_bagpipes_input.fits')
output_catalog = Table.read('/home/ariel/Workspace/GASP/HST/Data/f275w_dexp_logprior_single_bagpipes_results.fits')

input_catalog['fit_id'] = output_catalog['fit_id']

id_list = ['JO201_A259_f275w_dexp_logprior_single', 'JO201_B614_f275w_dexp_logprior_single',
           'JO201_B322_f275w_dexp_logprior_single']
note_list = ['Accepted Fit', 'Rejected by criterion a', 'Rejected by criterion b']

fig, ax = plt.subplots(3, 1, figsize=(9.25, 8.5), sharex=True)

for i in range(3):

    fit_id = id_list[i]

    print('Plotting clump ', fit_id)

    xlabel = False
    ylabel = False

    if i == 1:
        ylabel = True
    if i == 2:
        xlabel = True

    test = plot_example(fit_id, input_catalog, ax=ax[i], show_xlabel=xlabel, show_ylabel=ylabel,
                        spec_color='royalblue', plot_distributions=False)

    ax[i].annotate(note_list[i], xy=(0.98, 0.85), xycoords='axes fraction', ha='right', fontsize=16)

    ax[i].tick_params(axis='both', labelsize=12)

    sns.despine()

fig.subplots_adjust(wspace=0, hspace=0.065, left=0.1, bottom=0.08, right=0.98, top=0.98)

plt.savefig('/home/ariel/Workspace/GASP/HST/BAGPIPES/plots/good_bad_ugly.png')
plt.savefig('/home/ariel/Workspace/GASP/HST/BAGPIPES/plots/good_bad_ugly.pdf')

plt.close('all')
