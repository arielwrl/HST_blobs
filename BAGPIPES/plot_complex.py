"""

ariel@oapd
24/08/2022

Plots blob images and fits

"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import seaborn as sns
import astropy.visualization as vis
from astropy.wcs import WCS
import matplotlib.gridspec as gridspec
from scipy.ndimage.filters import gaussian_filter
from astropy.table import Table
from hst_pipeutils import plot_example
from toolbox import plot_tools

sns.set_style('ticks')

halpha_input = Table.read('/home/ariel/Workspace/GASP/HST/Data/tail_halpha_bagpipes_input.fits')
uv_input = Table.read('/home/ariel/Workspace/GASP/HST/Data/tail_f275w_bagpipes_input.fits')
complexes_input = Table.read('/home/ariel/Workspace/GASP/HST/Data/tail_f606w_bagpipes_input.fits')

halpha_id = 'JO201_A40_halpha_tail_parametric'
uv_id = 'JO201_A83_f275w_tail_parametric'
complex_id = 'JO201_1625_f606w_tail_parametric'

halpha_input['blob_id'] = [halpha_input['blob_id'][i] + '_' + 'parametric' for i in range(len(halpha_input))]
uv_input['blob_id'] = [uv_input['blob_id'][i] + '_' + 'parametric' for i in range(len(uv_input))]
complexes_input['blob_id'] = [complexes_input['blob_id'][i] + '_' + 'parametric' for i in range(len(complexes_input))]

mask_dir = '/home/ariel/Workspace/GASP/HST/Data/blob_masks/'
image_dir = '/home/ariel/Workspace/GASP/HST/Data/images/'

mask = fits.open(mask_dir+'JO201_complex46_masks.fits')
image = fits.open(image_dir+'JO201_f606w_0.04px_drc_subreg_MWdc.fits')

wcs = WCS(image[1].header)
image = image[1].data

fig = plt.figure(figsize=(16.25, 8))

gs = gridspec.GridSpec(3, 4)

ax_image = fig.add_subplot(gs[:, 0], projection=wcs)
ax_halpha = fig.add_subplot(gs[0, 1:4])
ax_f275w = fig.add_subplot(gs[1, 1:4])
ax_f606w = fig.add_subplot(gs[2, 1:4])

zscale = vis.ZScaleInterval()
z1, z2 = zscale.get_limits(image)
normalization = vis.ImageNormalize(vmin=0, vmax=z2*5, stretch=vis.AsinhStretch())

ax_image.imshow(image, cmap='binary', norm=normalization, origin='lower', interpolation='Gaussian')
ax_image.contour(mask[1].data, levels=0, cmap='Reds_r', lw=5)
ax_image.contour(mask[2].data, levels=0, cmap='Blues_r', lw=5)
ax_image.contour(mask[3].data, levels=0, cmap='Greens_r', lw=5)

ax_image.set_xlim(2070, 2130)
ax_image.set_ylim(2100, 2187)

ax_image.set_xlabel('$RA\,\mathrm{[\prime\prime]}$', fontsize=20)
ax_image.set_ylabel('$Dec\,\mathrm{[\prime\prime]}$', fontsize=20)

plot_example(halpha_id, halpha_input, ax=ax_halpha, plot_distributions=True, spec_color='#DD2D4A')
plot_example(uv_id, uv_input, ax=ax_f275w, plot_distributions=True)
plot_example(complex_id, complexes_input, ax=ax_f606w, plot_distributions=True, spec_color='#8DB580')

ax_halpha.set_xlabel('')
ax_halpha.set_ylabel('')
ax_f275w.set_xlabel('')
ax_f606w.set_ylabel('')

ax_halpha.annotate(r'$H\alpha$ Clump', xy=(0.01, 0.05), xycoords='axes fraction', fontsize=20)
ax_f275w.annotate(r'UV Clump', xy=(0.01, 0.05), xycoords='axes fraction', fontsize=20)
ax_f606w.annotate(r'SF Complex', xy=(0.01, 0.05), xycoords='axes fraction', fontsize=20)

ax_f275w.set_ylim(0, 6.9)

fig.subplots_adjust(top=0.975,
                    bottom=0.085,
                    left=0.075,
                    right=0.99,
                    hspace=0.0,
                    wspace=0.23)

plt.show()
