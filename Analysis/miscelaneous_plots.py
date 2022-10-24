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

input = Table.read('/home/ariel/Workspace/GASP/HST/Data/disk_halpha_bagpipes_input.fits')
output = Table.read('/home/ariel/Workspace/GASP/HST/Data/disk_halpha_parametric_bagpipes_results.fits')

flag = input['sel_flag'] == 31

input = input[flag]
output = output[flag]

sns.jointplot(x=np.log10(output['Ha']).tolist(), y=np.log10(0.65*input['F680N_line_flux']).tolist(), hue=input['gal'].tolist())

plt.ylabel(r'$\log \; 0.65 \times H\alpha$ Flux (Eric)', fontsize=20)
plt.xlabel(r'$\log \; H\alpha$ Flux (BAGPIPES)', fontsize=20)

x = np.linspace(-18, -12, 100)
y = x

plt.plot(x, y, '--k')

plt.show()