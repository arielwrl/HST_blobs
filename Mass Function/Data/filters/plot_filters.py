import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from toolbox import plot_tools

sns.set_style('ticks')


def plot_filters(paths_to_filters, filter_names):
    colors = sns.diverging_palette(220, 20, n=5)

    plt.figure(figsize=(6, 5))

    for i in range(len(paths_to_filters)):
        filter = paths_to_filters[i]

        wl, transmission = np.genfromtxt(filter).transpose()
        plt.plot(wl, transmission, color='k', lw=0.5)
        plt.fill_between(wl, np.zeros_like(transmission), transmission,
                         color=colors[i], alpha=0.5, label=filter_names[i])

    plt.legend(frameon=False, fontsize=16)

    plt.ylim(0)

    plt.xlabel(r'$\lambda \mathrm{[\AA]}$', fontsize=20)
    plt.ylabel(r'$T_\lambda$', fontsize=20)


filters = ['F275W', 'F336W', 'F606W', 'F680N', 'F814W']
paths = ['/home/ariel/Workspace/GASP/HST/Data/filters/HST_WFC3_UVIS2.' + filter_name + '.dat' for filter_name in filters]

plot_filters(paths, filters)