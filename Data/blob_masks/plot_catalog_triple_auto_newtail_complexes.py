import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
from matplotlib.transforms import Bbox

from astrodendro import Dendrogram, pp_catalog

from astropy import units as u
from astropy.io import fits
from astropy.table import Table
from astropy.visualization import ZScaleInterval
from astropy.wcs import WCS
import astropy.visualization as vis
from astropy.coordinates import SkyCoord
from astropy.cosmology import FlatLambdaCDM

cosmo = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Tcmb0=2.725 * u.K, Om0=0.3)

zscale = ZScaleInterval()


def read_fits(gal, filt, pixel_size, MWdc, masking=True, denoise=False, smooth=False):
    if masking is False and denoise is False and smooth is False:
        idir='{0}/north_flux'.format(gal[g])
        idata='{0}_{1}_{2}px_drc_subreg{3}'.format(gal,filt,pixel_size,MWdc)
    elif masking==True and denoise==False and smooth==False:
        idir_ima='{0}/north_flux'.format(gal)
        idata_ima='{0}_{1}_{2}px_drc_subreg{3}'.format(gal,filt,pixel_size,MWdc)
        idir='/media/eric/Samsung_T5/Programmi/HST/ds9_reg_masks'
        #idir_ima=idir
        idata='{0}_{1}_{2}px_reg_mask'.format(gal,filt,pixel_size)
        #idata_ima=idata
    elif masking==True and denoise==True and smooth==False:
        idir='../pysap_fits'
        idata='{0}_{1}_{2}px_reg_mask_denoised'.format(gal,filt,pixel_size)
    elif masking==True and denoise==False and smooth==True:
        idir='./smoothed'
        idata='{0}_{1}_{2}px_reg_mask_smoothed'.format(gal,filt,pixel_size)

    ima,hdr=fits.getdata('{0}/{1}.fits'.format(idir_ima,idata_ima),1,header=True)
    wcs=WCS(hdr)
    return ima,wcs,hdr,idata


def plot_clumps(ax,id_cat,dendros,cat,c):
    for i in range(len(id_cat)):#to select the catalog
        print(id_cat[i])
        for d in dendros[i]:
            if any(d.idx==j for j in cat['_idx'][cat['id_catalog']==id_cat[i]]):
                cond=(cat['_idx']==d.idx) & (cat['id_catalog']==id_cat[i])
                mask=d.get_mask()
                ax.contour(mask,levels=0, colors=c,zorder=int(d.level+1),linewidths=3)


def plot_lengthscale(bbox_pos,ax,scale_length,pixel_as_scale,pixel_kpc_scale):
    #print(bbox_pos)
    scalebar = AnchoredSizeBar(ax.transData, pixel_as_scale,
                               r'\textbf{1 kpc}',1,frameon=False,color='black',
                              #pad=0.3,borderpad=0.3,
                               size_vertical=1,sep=3,label_top=True,
                               fontproperties=mpl.font_manager.FontProperties(size=25),
                              bbox_to_anchor=Bbox.from_bounds(0,0,bbox_pos[0],bbox_pos[1]),
                               bbox_transform=ax.figure.transFigure)
    ax.add_artist(scalebar)



# %%
MUSE_CAT=False
bpt_sel=False
resolved=False

#"only trunks and leaves"
otl=1
otl_ped=['','_trunks_leaves']

masking=True
denoise=False
smooth=False

g=0

nsigma=2.5
delta=5
m=0

osigma=2
odelta=0

MWdc=['_MWdc','']
predriz=['','_corr']

gal=['JO201','JO204','JW100','JW39','JO175','JO206']
z=[0.0559,0.0451,0.0548,0.0634,0.0457,0.0486]
DL=cosmo.luminosity_distance(z) #Mpc
DL=DL.data

xmin,xmax,ymin,ymax=2020,2133,2005,2190

filt=['f275w','halpha']

pixel_size=0.04 #arcsec/pix


# %%
colors=['royalblue','red']
f_image='f606w'

print(gal[g])

# Setting pixel lengthscales
scale_length=1
pixel_kpc_scale=scale_length/(0.04*(u.arcsec).to(u.rad)*DL[g]/(1+z[g])**2.*1e3)
print(pixel_kpc_scale)
pixel_as_scale=scale_length/0.04


# Reading images and contours
ima,wcs,hdr,idata=read_fits(gal[g],f_image,pixel_size,MWdc[m])
#icont_mask='./RGB/Disk/contours/{0}_f814w_contours_1235sigma.fits'.format(gal[g])
#cont_mask,cont_hdr=fits.getdata(icont_mask,1,header=True)


################################# PLOTTING FIGURE ###############################################
# Setting the plot
gridspec = dict(wspace=0,width_ratios=[0.07, 1],height_ratios=[1,0.05],hspace=0)
fig,axs=plt.subplots(2,2,figsize=(7,9),subplot_kw=dict(projection=wcs),gridspec_kw=gridspec,
                    constrained_layout=True)
axs=axs.flatten()
ax=axs[1]
for i in [0,2,3]:
    axs[i].set_visible(False)


ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)
ax.coords[0].set_ticklabel(size=20,exclude_overlapping=True)
ax.coords[1].set_ticklabel(size=20,exclude_overlapping=True)
ax.coords[0].set_axislabel('',size=1)
ax.coords[1].set_axislabel('',size=1)

ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax)

# Plotting the F606W image
z1, z2 = zscale.get_limits(ima)
print(z1,z2)
norm = vis.ImageNormalize(vmin=0, vmax=z2*5, stretch=vis.AsinhStretch())
ax.set_aspect('equal')
ax.imshow(ima,origin='lower',norm=norm,cmap='binary')
#cont=ax.contour(cont_mask,transform=ax.get_transform(WCS(cont_hdr)),colors='dimgray',levels=[1,],
#                linestyles='--')

k=0

for f in [0,1]:
    # Loading dendrograms and catalogs.
    # Here I just need to get "idata"
    ima,wcs,hdr,idata=read_fits(gal[g],filt[f],pixel_size,MWdc[m])
    id_cat=['A','B','C']
    dendros=[]
    for i in id_cat[:2]:
        dendros.append(Dendrogram.load_from('{0}/dendro/{1}_{2}sigma_{3}delta_dendro{4}.fits'.format(gal[g],
                                                                                    idata,nsigma,delta,i)))
    dendros.append(Dendrogram.load_from('{0}/dendro/{1}_{2}sigma_{3}delta_denoise_{4}sigma_{5}delta_dendro{6}.fits'.format(gal[g],idata,nsigma,delta,osigma,odelta,id_cat[2])))

    if resolved==False:
        icat='HST_{0}_{1}px_reg_mask_SNR2_clumps_{2}sigma_{3}delta_denoise_{4}sigma_{5}delta_all_gals_z_sel'.format(filt[f],pixel_size,nsigma,delta,osigma,odelta)
    elif resolved==True:
        icat='HST_{0}_{1}px_reg_mask_SNR2_clumps_{2}sigma_{3}delta_denoise_{4}sigma_{5}delta_all_gals_z_sel_resolved'.format(filt[f],pixel_size,nsigma,delta,osigma,odelta)
    ifile='{0}.csv'.format(icat)
    cat=Table.read(ifile,format='ascii.ecsv')
    cat=cat[(cat['id_catalog']!='M') & (cat['gal']==gal[g])]
    if bpt_sel==True:
        cat=cat[(cat['BPT_flag']<2.5)]
        icat=icat.replace('z_sel','z_BPT_sel')


    if resolved==False:
        cat=cat[((cat['leaf_flag']==1) | (cat['level']==0))]

    xreg,yreg=wcs.world_to_pixel(SkyCoord(ra=cat['x_cen'].data, dec=cat['y_cen'].data, frame='icrs',
                                          unit=cat['x_cen'].unit))

    # Plotting clumps
    subima_cond=((xreg>xmin) & (xreg<xmax) & (yreg>ymin) & (yreg<ymax))
    print(np.sum(subima_cond))
    plot_clumps(ax,id_cat,dendros,cat[subima_cond],colors[f])


# Reading and plotting star-forming complexes
icompl_cat='HST_f606w_0.04px_reg_mask_denoised_disk275_masked_3sigma_0delta_only_trunks_match_all_gals.csv'
compl=Table.read(icompl_cat,format='ascii.ecsv')
compl=compl[(compl['gal']==gal[g]) & (compl['tail_gal_flag']==0)]

icont_compl='./{0}/dendro/{0}_f606w_0.04px_reg_mask_denoised_disk275_masked_3sigma_0delta_only_trunks_match.fits'.format(gal[g])
cont_compl=fits.getdata(icont_compl,2)
print(cont_compl)
cont_compl=np.array(cont_compl>0.1).astype(int)

ax.contour(cont_compl,colors='darkgreen',levels=[1,],linestyles='-',linewidths=3)

# Plotting lengthscale
plot_lengthscale([0.90,0.17],ax,scale_length,pixel_as_scale,pixel_kpc_scale)

oname=icat.replace('HST',gal[g])
plt.savefig('gal_plot/{0}/{0}_{1}_clumps_complexes.jpg'.format(gal[g],f_image),dpi=200)

plt.show()
plt.close()

