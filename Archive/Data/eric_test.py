from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from astropy.wcs import utils
import astropy.visualization as vis
from astropy.wcs.utils import proj_plane_pixel_scales


wcs = WCS(hdr)

#HST pixel size
pixel_size=0.04 #arcsec

#We read the MUSE regions and we convert the MUSE pix coord to HST pix coord
catM = Table.read('../MUSE/all_blobs_infos_allinfo_v3.cat', format='ascii')

#MUSE WCS
muse=fits.open('../MUSE/{0}/{0}_IMAGE_FOV_0001_v1.fits'.format(gal[g]))[1]
ima_muse,hdr_muse=muse.data,muse.header
wcs_muse=WCS(hdr_muse)

sx_muse, sy_muse = proj_plane_pixel_scales(wcs_muse) * 3600

#RESCALING TO HST
rmuse_to_hst=catM['r_pix']*sx_muse/pixel_size

skyreg=utils.pixel_to_skycoord(catM['x'], catM['y'], wcs=wcs_muse)
xmuse,ymuse=wcs_subhst.world_to_pixel(skyreg)