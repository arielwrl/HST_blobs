"""

ariel@oapd
17/02/2022

Simulate the effect of different parameters on the HST photometry

"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

bc03 = fits.open('/home/ariel/.local/lib/python3.8/site-packages/bagpipes/models/grids/bc03_miles_stellar_grids.fits')
nebular_cont = fits.open('/home/ariel/.local/lib/python3.8/site-packages/bagpipes/models/grids/bc03_miles_nebular_cont_'
                         'grids.fits')

wl = bc03[-1].data
ages = bc03[-2].data

fig = plt.figure()

metallicity = 0.005

for i in range(0, 16):
    print(i)
    plt.plot(wl, nebular_cont[22].data[i, :-1], label=str(ages[i]/1e9))

plt.legend()
plt.xlim(2000, 7000)

plt.show()