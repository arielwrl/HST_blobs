"""

ariel@isa
11/01/2022

Simulate the effect of different parameters on the HST photometry

"""

import numpy as np
import matplotlib.pyplot as plt
import bagpipes as pipes
import itertools
from starlight_toolkit.synphot import synmag
from starlight_toolkit.plotting import plot_filter
from astropy.table import Table
import seaborn as sns

sns.set_style('ticks')

data_dir = '/home/ariel/Workspace/GASP/HST/Data/'

filter_files = [data_dir + 'filters/HST_WFC3_UVIS2.F275W.dat',
                data_dir + 'filters/HST_WFC3_UVIS2.F336W.dat',
                data_dir + 'filters/HST_WFC3_UVIS2.F606W.dat',
                data_dir + 'filters/HST_WFC3_UVIS2.F680N.dat',
                data_dir + 'filters/HST_WFC3_UVIS2.F814W.dat']

logUs = [-2, -2.25, -2.5, -2.75, -3, -3.25, -3.5]
Avs = [0, 0.2, 0.4, 0.6, 0.8, 1, 2, 2.5]
metallicities = [0.005, 0.002, 0.02, 0.2, 0.4, 1., 2.5, 5.]

test_table = Table()

test_table = Table()
test_table['F275W'] = np.zeros(len(combinations))
test_table['F680N'] = np.zeros(len(combinations))
test_table['F336W'] = np.zeros(len(combinations))
test_table['age'] = np.zeros(len(combinations))
test_table['Av'] = np.zeros(len(combinations))


combinations = [(age, Av, metallicity) for age, Av, metallicity in itertools.product(ages, Avs, metallicities)]

for i in range(len(combinations)):
    age, Av = combinations[i]
    print(i, age, Av)

    exp = {}
    exp["age"] = age
    exp["tau"] = 0.001
    exp["massformed"] = 6
    exp['metallicity'] = metallicity

    dust = {}
    dust["type"] = "Cardelli"
    dust["Av"] = Av

    nebular = {}
    nebular["logU"] = logU

    dust["eta"] = 2

    model_components = {}
    model_components["redshift"] = 0.044631  # Observed redshift
    model_components["delayed"] = exp  # The star-formation history dictionary
    model_components["dust"] = dust  # The dust component
    model_components["veldisp"] = 75
    model_components["nebular"] = nebular  # The dust component

    # Wavelength array to plot spectra.
    wl = np.arange(1000, 10000, 1)

    # Creating a model galaxy
    model = pipes.model_galaxy(model_components, spec_wavs=wl)

    model.plot(show=False)
    plt.savefig('/home/ariel/Workspace/GASP/HST/Data/simulated_spectra/dust' + str(Av) + '_logU' +
               str(logU) + '_metal' + str(metallicity) + 'age+' + str(exp['age']) + '.png')
    plt.close()

    np.savetxt('/home/ariel/Workspace/GASP/HST/Data/simulated_spectra/dust' + str(Av) + '_logU' +
               str(logU) + '_metal' + str(metallicity) + 'age+' + str(exp['age']) + '.spec', model.spectrum)

    np.savetxt('/home/ariel/Workspace/GASP/HST/Data/simulated_spectra/dust' + str(Av) + '_logU' +
               str(logU) + '_metal' + str(metallicity) + '_weirdcase.spec', model.spectrum)

    test_table['logU'][i] = logU
    test_table['Av'][i] = Av
    test_table['gas_metallicity'][i] = metallicity
    test_table['stellar_metallicity'][i] = model.model_comp['delayed']['metallicity']
    test_table['F275W'][i] = synmag(model.spectrum[:, 0], model.spectrum[:, 1], filter_curve=filter_files[0])
    test_table['F680N'][i] = synmag(model.spectrum[:, 0], model.spectrum[:, 1], filter_curve=filter_files[3])
    test_table['F606W'][i] = synmag(model.spectrum[:, 0], model.spectrum[:, 1], filter_curve=filter_files[2])
    test_table['Ha'][i] = model.line_fluxes['H  1  6562.81A']
    test_table['Nii'][i] = model.line_fluxes['N  2  6583.45A']
    test_table['Oiii'][i] = model.line_fluxes['O  3  5006.84A']
    test_table['Hb'][i] = model.line_fluxes['H  1  4861.33A']


test_table.write('/home/ariel/Workspace/GASP/HST/Data/simulated_spectra/simulations_table.fits', overwrite=True)






# ages = [0.02, 0.04, 0.08, 0.1]
# Avs = [0.0, 0.1, 0.25, 0.5]
# spectra_dict = {}
#
# combinations = [(age, Av) for age, Av in itertools.product(ages, Avs)]
#
# test_table = Table()
# test_table['F275W'] = np.zeros(len(combinations))
# test_table['F680N'] = np.zeros(len(combinations))
# test_table['F336W'] = np.zeros(len(combinations))
# test_table['age'] = np.zeros(len(combinations))
# test_table['Av'] = np.zeros(len(combinations))
#
# for i in range(len(combinations)):
#
#     age, Av = combinations[i]
#     print(i, age, Av)
#
#     exp = {}
#     exp["age"] = age
#     exp["tau"] = 1
#     exp["massformed"] = 5.0
#     exp['metallicity'] = 1
#
#     dust = {}
#     dust["type"] = "Cardelli"
#     dust["Av"] = Av
#
#     nebular = {}
#     nebular["logU"] = -2.5
#
#     dust["eta"] = 1.8
#
#     model_components = {}
#     model_components["redshift"] = 0.044631  # Observed redshift
#     model_components["delayed"] = exp  # The star-formation history dictionary
#     model_components["dust"] = dust  # The dust component
#     model_components["veldisp"] = 75
#     model_components["nebular"] = nebular  # The dust component
#
#     # Wavelength array to plot spectra.
#     wl = np.arange(1000, 10000, 1)
#
#     # Creating a model galaxy
#     model = pipes.model_galaxy(model_components, spec_wavs=wl)
#
#     spectra_dict[str(age)+'_'+str(Av)] = model.spectrum[:, 1]
#
#     test_table['age'][i] = age
#     test_table['Av'][i] = Av
#     test_table['F275W'][i] = synmag(model.spectrum[:, 0], model.spectrum[:, 1], filter_curve=filter_files[0])
#     test_table['F336W'][i] = synmag(model.spectrum[:, 0], model.spectrum[:, 1], filter_curve=filter_files[1])
#     test_table['F680N'][i] = synmag(model.spectrum[:, 0], model.spectrum[:, 1], filter_curve=filter_files[3])
#
#
# # Without dust
# fig = plt.figure()
#
# palette = sns.color_palette('Spectral_r', 6)
#
# plt.plot(wl, spectra_dict['0.02_0.0'], color=palette[0], label=r'$20\,\mathrm{Myr}$')
# plt.plot(wl, spectra_dict['0.04_0.0'], color=palette[1], label=r'$40\,\mathrm{Myr}$')
# plt.plot(wl, spectra_dict['0.08_0.0'], color=palette[5], label=r'$80\,\mathrm{Myr}$')
#
# plt.legend(frameon=False, fontsize=12)
#
# for file in filter_files:
#     plot_filter(file, scale_factor=0.1e-17, filter_ls='-')
#
# plt.ylabel(r'$\mathrm{Flux \, [erg/s/cm^2/\AA]}$', fontsize=15)
# plt.xlabel(r'$\lambda \, \mathrm{[\AA]}$', fontsize=15)
#
# plt.ylim(0, 0.5e-17)
# plt.xlim(2100, 9500)
#
# fig.tight_layout()
#
# fig.savefig('models_ages.png', dpi=300)
#
# plt.show()
#
#
# # What dust can do:
#
# fig = plt.figure()
#
# palette = sns.color_palette('magma', 6)
#
# plt.plot(wl, spectra_dict['0.02_0.25'], color=palette[2], label=r'$20\,\mathrm{Myr}, A_V=0.25$')
# plt.plot(wl, spectra_dict['0.08_0.0'], color=palette[4], label=r'$80\,\mathrm{Myr}$')
#
# plt.legend(frameon=False, fontsize=12)
#
# # for file in filter_files:
# #     plot_filter(file, scale_factor=1e-17)
#
# plt.ylabel(r'$\mathrm{Flux \, [erg/s/cm^2/\AA]}$', fontsize=15)
# plt.xlabel(r'$\lambda \, \mathrm{[\AA]}$', fontsize=15)
#
# plt.ylim(0, 2.5e-17)
# plt.xlim(2200, 9500)
#
# fig.tight_layout()
#
# fig.savefig('models_ages_dust.png', dpi=300)
#
# plt.show()





