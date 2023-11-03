import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from toolbox import plot_tools

sns.set_style('ticks')


def plot_filters(paths_to_filters, filter_names, colors = sns.diverging_palette(220, 20, n=5)):

    fig = plt.figure(figsize=(5.5, 4.25))

    for i in range(len(paths_to_filters)):
        filter = paths_to_filters[i]


        if i in [0, 2, 3]:

            wl, transmission = np.genfromtxt(filter).transpose()
            plt.plot(wl, transmission, color=colors[i], lw=1, label=filter_names[i])
            
        else:
            
            wl, transmission = np.genfromtxt(filter).transpose()
            plt.plot(wl, transmission, color=colors[i], lw=1, label=filter_names[i], ls='dashed')
            
        
            # plt.fill_between(wl, np.zeros_like(transmission), transmission,
            #                 color=colors[i], alpha=0.1)


    plt.legend(frameon=False, fontsize=11, loc=0)

    plt.ylim(0)

    plt.xlabel(r'$\lambda \mathrm{[\AA]}$', fontsize=16)
    plt.ylabel(r'$T_\lambda$', fontsize=16)

    plt.tick_params(axis='both', labelsize=11)

    sns.despine()

    fig.subplots_adjust(top=0.965,
                        bottom=0.136,
                        left=0.125,
                        right=0.968,
                        hspace=0.2,
                        wspace=0.2)


filters = ['F275W', 'F336W', 'F606W', 'F680N', 'F814W']
paths = ['/home/ariel/Workspace/GASP/HST/Data/filters/HST_WFC3_UVIS2.' + filter_name + '.dat' for filter_name in filters]


colors = ['mediumvioletred', '#0058b0', 'indigo', 'goldenrod', '#990147']

plot_filters(paths, filters, colors=colors)
plt.savefig('filters.jpg', dpi=300)
plt.savefig('filters.pdf')
# plt.show()
