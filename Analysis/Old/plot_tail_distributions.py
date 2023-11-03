"""

ariel@oapd
15/08/2022

Plots distribution of ages and masses for clumps in the tails


"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from astropy.table import Table
import seaborn as sns
from toolbox import plot_tools

sns.set_style('ticks')

distribution_variable = 'mwage'

mass_dict = {'JO201': 44194800000,
             'JO204': 54968402000,
             'JW100': 292875993000,
             'JW39': 164373004000,
             'JO175': 33957900300,
             'JO206': 77743301000
}

halpha_input = Table.read('/home/ariel/Workspace/GASP/HST/Data/tail_halpha_bagpipes_input.fits')
f275w_input = Table.read('/home/ariel/Workspace/GASP/HST/Data/tail_f275w_bagpipes_input.fits')
f606w_input = Table.read('/home/ariel/Workspace/GASP/HST/Data/tail_f606w_bagpipes_input.fits')

halpha_input = halpha_input[(halpha_input['sel_flag'] == 31)]
f275w_input = f275w_input[(f275w_input['sel_flag'] == 31)]
f606w_input = f606w_input[(f606w_input['sel_flag'] == 31)]

tail_halpha = Table.read('/home/ariel/Workspace/GASP/HST/Data/tail_halpha_dexp_logprior_bagpipes_results.fits')
tail_f275w = Table.read('/home/ariel/Workspace/GASP/HST/Data/tail_f275w_dexp_logprior_bagpipes_results.fits')
tail_f606w = Table.read('/home/ariel/Workspace/GASP/HST/Data/tail_f606w_dexp_logprior_bagpipes_results.fits')

tail_halpha = tail_halpha[(halpha_input['sel_flag'] == 31) & (halpha_input['tail_gal_flag'] == 0)]
tail_f275w = tail_f275w[(f275w_input['sel_flag'] == 31) & (f275w_input['tail_gal_flag'] == 0)]
tail_f606w = tail_f606w[(f606w_input['sel_flag'] == 31) & (f606w_input['tail_gal_flag'] == 0)]

tail_halpha['galaxy'] = tail_halpha['galaxy'].astype(str)
tail_f275w['galaxy'] = tail_f275w['galaxy'].astype(str)
tail_f606w['galaxy'] = tail_f606w['galaxy'].astype(str)

tail_halpha['galaxy_mass'] = np.zeros_like(tail_halpha['mwage'])
for i in range(len(tail_halpha)):
    tail_halpha['galaxy_mass'][i] = np.log10(mass_dict[tail_halpha['galaxy'][i]])

tail_f275w['galaxy_mass'] = np.zeros_like(tail_f275w['mwage'])
for i in range(len(tail_f275w)):
    tail_f275w['galaxy_mass'][i] = np.log10(mass_dict[tail_f275w['galaxy'][i]])

tail_f606w['galaxy_mass'] = np.zeros_like(tail_f606w['mwage'])
for i in range(len(tail_f606w)):
    tail_f606w['galaxy_mass'][i] = np.log10(mass_dict[tail_f606w['galaxy'][i]])

tail_halpha['mwage'] *= 1e3
tail_f275w['mwage'] *= 1e3
tail_f606w['mwage'] *= 1e3

tail_halpha['mass'] = tail_halpha['stellar_mass']
tail_f275w['mass'] = tail_f275w['stellar_mass']
tail_f606w['mass'] = tail_f606w['stellar_mass']

halpha_df = tail_halpha['galaxy', 'blob_id', 'age', 'tau', 'mwage', 'mass', 'galaxy_mass', 'sfr', 'Av'].to_pandas()
f275w_df = tail_f275w['galaxy', 'blob_id', 'age', 'tau', 'mwage', 'mass', 'galaxy_mass', 'sfr', 'Av'].to_pandas()
f606w_df = tail_f606w['galaxy', 'blob_id', 'age', 'tau', 'mwage', 'mass', 'galaxy_mass', 'sfr', 'Av'].to_pandas()

halpha_df.sort_values(by='galaxy_mass', inplace=True)
f275w_df.sort_values(by='galaxy_mass', inplace=True)
f606w_df.sort_values(by='galaxy_mass', inplace=True)


fig, ax = plt.subplots(1, 3, figsize=(10, 6.5), sharey=True)

box_halpha = ax[0]
box_f275w = ax[1]
box_f606w = ax[2]

hist_halpha = box_halpha.inset_axes([0, 0, 1, 0.66], sharex=box_halpha)
hist_f275w = box_f275w.inset_axes([0, 0, 1, 0.66], sharex=box_f275w)
hist_f606w = box_f606w.inset_axes([0, 0, 1, 0.66], sharex=box_f606w)

hist_halpha.patch.set_alpha(0.0)
hist_f275w.patch.set_alpha(0.0)
hist_f606w.patch.set_alpha(0.0)

sns.boxplot(y='galaxy', x=distribution_variable, data=halpha_df, whis=[15, 85], width=.6, palette='light:#DD2D4A'
            , fliersize=0, ax=box_halpha)
sns.boxplot(y='galaxy', x=distribution_variable, data=f275w_df, whis=[15, 85], width=.6, palette='light:#758BFD'
            , fliersize=0, ax=box_f275w)
sns.boxplot(y='galaxy', x=distribution_variable, data=f606w_df, whis=[15, 85], width=.6, palette='light:#8DB580'
            , fliersize=0, ax=box_f606w)

sns.kdeplot(x=distribution_variable, data=halpha_df, color='#DD2D4A', ax=hist_halpha, shade=True, alpha=0.25)
sns.kdeplot(x=distribution_variable, data=f275w_df, color='#758BFD', ax=hist_f275w, shade=True, alpha=0.25)
sns.kdeplot(x=distribution_variable, data=f606w_df, color='#8DB580', ax=hist_f606w, shade=True, alpha=0.25)

hist_halpha.set_xlim(np.percentile(tail_halpha[distribution_variable], 2.5),
                     np.percentile(tail_halpha[distribution_variable], 95))
hist_f275w.set_xlim(np.percentile(tail_f275w[distribution_variable], 2.5),
                    np.percentile(tail_f275w[distribution_variable], 95))
hist_f606w.set_xlim(np.percentile(tail_f606w[distribution_variable][~np.isnan(tail_f606w[distribution_variable])], 2.5),
                    np.percentile(tail_f606w[distribution_variable][~np.isnan(tail_f606w[distribution_variable])], 95))

sns.despine(ax=hist_halpha, left=True, bottom=True)
sns.despine(ax=hist_f275w, left=True, bottom=True)
sns.despine(ax=hist_f606w, left=True, bottom=True)

hist_halpha.tick_params(left=False, labelleft=False, bottom=False, labelbottom=False)
hist_f275w.tick_params(left=False, labelleft=False, bottom=False, labelbottom=False)
hist_f606w.tick_params(left=False, labelleft=False, bottom=False, labelbottom=False)

hist_halpha.set_ylabel('')
hist_f275w.set_ylabel('')
hist_f606w.set_ylabel('')

hist_halpha.set_xlabel('')
hist_f275w.set_xlabel('')
hist_f606w.set_xlabel('')

box_halpha.set_ylabel('')
box_f275w.set_ylabel('')
box_f606w.set_ylabel('')

box_halpha.set_xlabel('')
box_f606w.set_xlabel('')

if distribution_variable == 'mwage':
    box_f275w.set_xlabel(r'$\langle\,t_\star\,\rangle_M\,\mathrm{[Myr]}$', fontsize=20)
elif distribution_variable == 'mass':
    box_f275w.set_xlabel(r'$\log\,M/M_\odot$', fontsize=20)
elif distribution_variable == 'sfr':
    box_f275w.set_xlabel(r'SFR$\mathrm{[M_\odot/yr]}$', fontsize=20)

box_halpha.tick_params(labelsize=20, axis='y')
box_f275w.tick_params(left=False, labelleft=False)
box_f606w.tick_params(left=False, labelleft=False)

fig.subplots_adjust(top=0.98, bottom=0.1, left=0.09, right=0.995, hspace=0.0, wspace=0.0)

plt.savefig('/home/ariel/Workspace/GASP/HST/Analysis/Plots/'+distribution_variable+'.png', dpi=300)

plt.show()

