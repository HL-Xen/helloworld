from __future__ import division
from numpy import *
import matplotlib.pyplot as plt
from matplotlib.collections import EllipseCollection
from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy.io import fits
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.table import Table

#%% Load data

img_f560w = fits.open('jw01328-c1006_t014_miri_f560w_i2d.fits', memmap=True)[1].data
sources_f560w = Table.read('jw01328-c1006_t014_miri_f560w_i2d_SSCs_matched.csv', format='ascii.csv')
norm_f560w = ImageNormalize(vmin=0,vmax=200,stretch=SqrtStretch())

img_f770w = fits.open('jw01328-c1006_t014_miri_f770w_i2d.fits', memmap=True)[1].data
sources_f770w = Table.read('jw01328-c1006_t014_miri_f770w_i2d_SSCs_matched.csv', format='ascii.csv')
norm_f770w = ImageNormalize(vmin=0,vmax=200,stretch=SqrtStretch())

img_f1500w = fits.open('jw01328-c1006_t014_miri_f1500w_i2d.fits', memmap=True)[1].data
sources_f1500w = Table.read('jw01328-c1006_t014_miri_f1500w_i2d_SSCs_matched.csv', format='ascii.csv')
norm_f1500w = ImageNormalize(vmin=0,vmax=200,stretch=SqrtStretch())

#%% Plot images with detected sources
fig, (ax_f560w, ax_f770w, ax_f1500w) = plt.subplots(1,3,figsize=(15,10))

# plot the image and mark the location of detected sources
im_f560w = ax_f560w.imshow(img_f560w, interpolation='none', cmap='hot', origin='lower', norm=norm_f560w)
sc_f560w = ax_f560w.add_collection(EllipseCollection(full(size(sources_f560w),8),full(size(sources_f560w),8),zeros(size(sources_f560w)),
                         units='xy',facecolor='none',edgecolor='white', offsets=list(zip(sources_f560w['X'],sources_f560w['Y'])),
                         transOffset=ax_f560w.transData))

ax_f560w.set_title('F560W')
divider = make_axes_locatable(ax_f560w)
cax_f560w = divider.append_axes("right", size="5%", pad=0.1)
fig.colorbar(im_f560w, cax=cax_f560w)

im_f770w = ax_f770w.imshow(img_f770w, interpolation='none', cmap='hot', origin='lower', norm=norm_f770w)
sc_f770w = ax_f770w.add_collection(EllipseCollection(full(size(sources_f770w),8),full(size(sources_f770w),8),zeros(size(sources_f770w)),
                         units='xy',facecolor='none',edgecolor='white', offsets=list(zip(sources_f770w['X'],sources_f770w['Y'])),
                         transOffset=ax_f770w.transData))

ax_f770w.set_title('FF770W')
divider = make_axes_locatable(ax_f770w)
cax_f770w = divider.append_axes("right", size="5%", pad=0.1)
fig.colorbar(im_f770w, cax=cax_f770w)

im_f1500w = ax_f1500w.imshow(img_f1500w, interpolation='none', cmap='hot', origin='lower', norm=norm_f1500w)
sc_f1500w = ax_f1500w.add_collection(EllipseCollection(full(size(sources_f1500w),8),full(size(sources_f1500w),8),zeros(size(sources_f1500w)),
                         units='xy',facecolor='none',edgecolor='white', offsets=list(zip(sources_f1500w['X'],sources_f1500w['Y'])),
                         transOffset=ax_f1500w.transData))

ax_f1500w.set_title('F1500W')
divider = make_axes_locatable(ax_f1500w)
cax_f1500w = divider.append_axes("right", size="5%", pad=0.1)
fig.colorbar(im_f1500w, cax=cax_f1500w)

annot_f560w = ax_f560w.annotate("", xy=(0,0), xytext=(-20,-20),textcoords="offset points",
                    bbox=dict(boxstyle="round", fc="w"),
                    arrowprops=dict(arrowstyle="->",color='w'))
annot_f770w = ax_f770w.annotate("", xy=(0,0), xytext=(-20,-20),textcoords="offset points",
                    bbox=dict(boxstyle="round", fc="w"),
                    arrowprops=dict(arrowstyle="->",color='w'))
annot_f1500w = ax_f1500w.annotate("", xy=(0,0), xytext=(-20,-20),textcoords="offset points",
                    bbox=dict(boxstyle="round", fc="w"),
                    arrowprops=dict(arrowstyle="->",color='w'))

annot_f560w.set_visible(False)
annot_f770w.set_visible(False)
annot_f1500w.set_visible(False)

def update_annot(ind):
    pos_f560w = sc_f560w.get_offsets()[ind["ind"][0]]
    annot_f560w.xy = pos_f560w
    text_f560w = "{}".format(" ".join(list(map(str,ind["ind"]))))
    annot_f560w.set_text(text_f560w)
    
    pos_f770w = sc_f770w.get_offsets()[ind["ind"][0]]
    annot_f770w.xy = pos_f770w
    text_f770w = "{}".format(" ".join(list(map(str,ind["ind"]))))
    annot_f770w.set_text(text_f770w)
    
    pos_f1500w = sc_f1500w.get_offsets()[ind["ind"][0]]
    annot_f1500w.xy = pos_f1500w
    text_f1500w = "{}".format(" ".join(list(map(str,ind["ind"]))))
    annot_f1500w.set_text(text_f1500w)

def hover(event):
    vis_f560w = annot_f560w.get_visible()
    vis_f770w = annot_f770w.get_visible()
    vis_f1500w = annot_f1500w.get_visible()

    if event.inaxes in [ax_f560w,ax_f770w,ax_f1500w]:
        if event.inaxes == ax_f560w:
            cont, ind = sc_f560w.contains(event)
        elif event.inaxes == ax_f770w:
            cont, ind = sc_f770w.contains(event)
        elif event.inaxes == ax_f1500w:
            cont, ind = sc_f1500w.contains(event)
        if cont:
            update_annot(ind)
            annot_f560w.set_visible(True)
            annot_f770w.set_visible(True)
            annot_f1500w.set_visible(True)
            fig.canvas.draw_idle()                
        else:
            if vis_f560w:
                annot_f560w.set_visible(False)          
            if vis_f770w:
                annot_f770w.set_visible(False)
            if vis_f1500w:
                annot_f1500w.set_visible(False)
            fig.canvas.draw_idle()
  
plt.tight_layout()
fig.canvas.mpl_connect("motion_notify_event", hover)
plt.show()

