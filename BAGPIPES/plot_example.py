"""

ariel@oapd
22/08/2022

Makes a cute plot with an SED fitting example

"""

import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt
import seaborn as sns
from toolbox import plot_tools
from hst_pipeutils import plot_example_disk

sns.set_style('ticks')
plt.ioff()

halpha_input = Table.read('/home/ariel/Workspace/GASP/HST/Data/disk_halpha_bagpipes_input.fits')

blob_id = 'JO175_A99_halpha_disk_parametric'

halpha_input['blob_id'] = [halpha_input['blob_id'][i] + '_' + 'parametric' for i in range(len(halpha_input))]

fig, ax = plt.subplots(1, 1, figsize=(12.5, 5.5))

plot_example_disk(blob_id, halpha_input)

fig.tight_layout()


