"""

ariel@sãojoão
11/04/2023

Looks at quenched and unquenched complexes

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from astropy.table import Table
import seaborn as sns
from toolbox import plot_tools
from toolbox.wololo import arcsectokpc

sns.set_style('ticks')
f606w_palette = sns.light_palette('indigo', 4)
quench_palette = sns.color_palette('Spectral', 10)

input_f606w = Table.read('/home/ariel/Workspace/GASP/HST/Data/f606w_bagpipes_input.fits')
output_f606w = Table.read('/home/ariel/Workspace/GASP/HST/Data/f606w_dexp_logprior_single_bagpipes_results.fits')

# Looking at the drop from peak SFR

drop = output_f606w['isfr']/output_f606w['max_sfr']
quenched = drop < 0.05
active = drop > 0.95
transition = (drop < 0.95) & (drop > 0.05)

output_f606w['quenching_class'] = np.zeros_like(output_f606w['clump_id'])
output_f606w['quenching_class'][quenched] = np.full(quenched.sum(), 'Quenched')
output_f606w['quenching_class'][active] = np.full(active.sum(), 'Star-Forming')
output_f606w['quenching_class'][transition] = np.full(transition.sum(), 'Transition')

fig = plt.figure()

plt.hist(drop[quenched], range=[0, 0.05], bins=2, color=quench_palette[0])
plt.hist(drop[active], range=[0.95, 1], bins=2, color=quench_palette[-1])
plt.hist(drop[transition], range=[0.05, 0.95], bins=36, color=quench_palette[7])

plt.axvline(x=0.05, linestyle='dashed', color='k')
plt.axvline(x=0.95, linestyle='dashed', color='k')

plt.annotate('Quenched', xy=(0.1, 1.01), xycoords='axes fraction', fontsize=15, ha='center')
plt.annotate('Star-forming', xy=(0.9, 1.01), xycoords='axes fraction', fontsize=15, ha='center')
# plt.annotate('Transition', xy=(0.5, 1.01), xycoords='axes fraction', fontsize=15, ha='center')

sns.despine()

plt.xlabel(r'Quenching Parameter', fontsize=23)
plt.ylabel(r'N', fontsize=20)

fig.tight_layout()

plt.savefig('/home/ariel/Workspace/GASP/HST/Analysis/Plots/quenching/quenching.pdf')


# Check the distributions of several parameters

conversion_factor = arcsectokpc(input_f606w['galaxy_redshift'])

parameter_dict = {}
parameter_dict['mwage'] = 1e3*output_f606w['mwage']
parameter_dict['age_t0'] = 1e3*output_f606w['age']
parameter_dict['log_tau'] = np.log10(output_f606w['tau'])
parameter_dict['t_max_sfr'] = output_f606w['t_max_sfr']/1e6
parameter_dict['stellar_mass'] = output_f606w['stellar_mass']
parameter_dict['area'] = input_f606w['area_exact'] * conversion_factor ** 2
parameter_dict['stellar_mass_density'] = np.log10((10**output_f606w['stellar_mass'])/(input_f606w['area_exact'] * conversion_factor ** 2))
parameter_dict['max_sfr_density'] = np.log10(output_f606w['max_sfr']/(input_f606w['area_exact'] * conversion_factor ** 2))
parameter_dict['d'] = input_f606w['dist_kpc']
parameter_dict['max_sfr'] = np.log10(output_f606w['max_sfr'])

for parameter in parameter_dict:

    fig = plt.figure()
    sns.violinplot(x=parameter_dict[parameter].tolist(), y=output_f606w['quenching_class'].tolist(),
                   fliersize=0, order=['Star-Forming', 'Quenched'])
    plt.xlabel(parameter, fontsize=20)

    plt.xlim(np.percentile(parameter_dict[parameter], 2), np.percentile(parameter_dict[parameter], 98))

    fig.tight_layout()

    plt.savefig('/home/ariel/Workspace/GASP/HST/Analysis/Plots/quenching/' + parameter + '.png')

plt.close('all')

parameter_dict['axial_ratio'] = input_f606w['axial_ratio']
parameter_dict['f_Ha'] = input_f606w['dA_ha']
parameter_dict['f_UV'] = input_f606w['dA_uv']

for parameter in ['f_Ha', 'f_UV', 'axial_ratio']:

    flag = ~parameter_dict[parameter].mask

    fig = plt.figure()
    sns.violinplot(x=parameter_dict[parameter][flag].tolist(), y=output_f606w['quenching_class'][flag].tolist(),
                   fliersize=0, order=['Star-Forming', 'Quenched'])
    plt.xlabel(parameter, fontsize=20)

    plt.xlim(np.percentile(parameter_dict[parameter][flag], 2), np.percentile(parameter_dict[parameter][flag], 98))

    fig.tight_layout()

    plt.savefig('/home/ariel/Workspace/GASP/HST/Analysis/Plots/quenching/' + parameter + '.png')

plt.close('all')