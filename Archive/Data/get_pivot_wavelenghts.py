"""

ariel@ariel's bed
25/01/2021

Get pivot wavelengths for HST filters

"""

from starlight_toolkit.synphot import pivot_wavelength

filter_files = ['../Data/filters/HST_WFC3_UVIS2.F275W.dat',
                '../Data/filters/HST_WFC3_UVIS2.F336W.dat',
                '../Data/filters/HST_WFC3_UVIS2.F606W.dat',
                '../Data/filters/HST_WFC3_UVIS2.F680N.dat',
                '../Data/filters/HST_WFC3_UVIS2.F814W.dat']

pivot_wavelengths = [pivot_wavelength(filter_file) for filter_file in filter_files]

print(pivot_wavelengths)