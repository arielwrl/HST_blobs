"""

ariel@oapd
13/08/2022

Goes through all filter files and convert wavelengths to angstrom for BAGPIPES

"""

import numpy as np
import os

filter_dir = './filters/'

filter_files = os.listdir(filter_dir)

for filter_file in filter_files:

    print(filter_file)

    wl, transmission = np.genfromtxt(filter_dir + filter_file).transpose()
    np.savetxt(filter_dir + filter_file.split('_')[0] + '.dat', np.array([10000*wl, transmission]).transpose())
