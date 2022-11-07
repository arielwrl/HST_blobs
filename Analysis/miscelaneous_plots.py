"""

ariel@oapd
18/10/2022

The file name explains itself

"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack
import seaborn as sns
from toolbox import plot_tools
from astropy.io import fits

sns.set_style('ticks')

input = Table.read('/home/ariel/Workspace/GASP/HST/Data/tail_f275w_bagpipes_input.fits')
output = Table.read('/home/ariel/Workspace/GASP/HST/Data/tail_f275w_parametric_bagpipes_results.fits')

flag = input['sel_flag'] == 31

input = input[flag]
output = output[flag]

plt.figure()

sns.kdeplot(y=output['mwage'].tolist(), x=output['stellar_mass'].tolist(),
            hue=input['tail_gal_flag'].tolist())

plt.ylabel(r'mwage', fontsize=20)
plt.xlabel(r'$\log M/M_\odot$', fontsize=20)

plt.show()


plt.figure()

sns.jointplot(x=np.log10(output['Ha']).tolist(), y=np.log10(input['F680N_line_flux']).tolist(), hue=input['gal'].tolist())

plt.ylabel(r'$\log \; H\alpha$ Flux (Eric)', fontsize=20)
plt.xlabel(r'$\log \; H\alpha$ Flux (BAGPIPES)', fontsize=20)

x = np.linspace(-18, -12, 100)
y = x

plt.plot(x, y, '--k')

plt.show()


# Disks:

# input = Table.read('/home/ariel/Workspace/GASP/HST/Data/disk_halpha_bagpipes_input.fits')
# output = Table.read('/home/ariel/Workspace/GASP/HST/Data/disk_halpha_parametric_bagpipes_results.fits')

# flag = input['sel_flag'] == 31
#
# input = input[flag]
# output = output[flag]
