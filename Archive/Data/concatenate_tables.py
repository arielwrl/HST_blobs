"""

ariel@padova
12/04/2021

Puts together tables with results from all galaxies

"""

from astropy.table import Table, vstack


galaxy_list = ['JO201', 'JO204', 'JO206', 'JW100']

table = Table.read('/home/ariel/Workspace/GASP/HST/Data/' + galaxy_list[0] + '_results_burst.fits')

for i in range(1, len(galaxy_list)):
    table = vstack([table, Table.read('/home/ariel/Workspace/GASP/HST/Data/' + galaxy_list[i] + '_results_burst.fits')])

table.write('bagpipes_results_burst.fits', overwrite=True)
