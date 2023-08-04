"""

ariel@oapd
06/12/2022

Simulates JWST blobservations for Eric

"""

import numpy as np
from astropy.table import Table
from starlight_toolkit.synphot import synflux

config = 'dexp_logprior'
detection = 'optical_only'
spec_dir = '/home/ariel/Workspace/GASP/HST/Data/full_spectra/'
filter_dir = '/home/ariel/Workspace/GASP/HST/JWST Simulation/filters/'

bagpipes_input = Table.read('/home/ariel/Workspace/GASP/HST/Data/tail_'+detection +
                            '_bagpipes_input.fits')
bagpipes_results = Table.read('/home/ariel/Workspace/GASP/HST/Data/tail_'+detection +
                              '_dexp_logprior_bagpipes_results.fits')

flag = bagpipes_input['sel_flag'] == 31
bagpipes_input = bagpipes_input[flag]

filter_list = ['F200W', 'F300M', 'F335M', 'F360M', 'F770W', 'F1000W', 'F1130W', 'F2100W']

photometry_table = Table()
photometry_table['id'] = bagpipes_input['blob_id']

for filter_name in filter_list:

    photometry_table[filter_name] = np.zeros_like(bagpipes_input['F275W'])

    for i in range(len(bagpipes_input)):

        wl, flux, flux_25, flux_75 = np.genfromtxt(spec_dir+bagpipes_input['fit_id'][i]+'_'+config+'.dat').transpose()

        z = bagpipes_input['galaxy_redshift'][i]
        wl *= (1 + z)
        flux /= (1 + z)

        photometry_table[filter_name][i] = synflux(wl, flux, filter_curve=filter_dir+filter_name+'.dat') * 1e-18

photometry_table.write('/home/ariel/Workspace/GASP/HST/JWST Simulation/faketometry' + detection + '.fits',
                       overwrite=True)