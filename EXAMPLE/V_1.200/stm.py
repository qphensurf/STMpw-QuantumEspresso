#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import scipy.ndimage as ndimage
from matplotlib import cm

data=np.loadtxt('Topography.gnu')

fs=16

X=data[:,0]
Y=data[:,1]
xdim = np.count_nonzero(X == 0)
ydim = np.count_nonzero(Y== 0)
Z=data[:,2].reshape(ydim,xdim)

Z2 = ndimage.gaussian_filter(Z, sigma=1.0, order=0)
Z2=Z2.T

ratio=xdim/ydim
fig, axs = plt.subplots(1,1, figsize=(8,8*ratio))

IM=axs.imshow(Z2, extent=[0,np.max(X)*0.529, 0, np.max(Y)*0.529],cmap='afmhot')

divider = make_axes_locatable(axs)
cax = divider.append_axes("right", size="5%", pad=0.05)
cb = plt.colorbar(IM, cax=cax)
cb.ax.tick_params(labelsize=fs) 

axs.tick_params(axis='both', labelsize = fs)
axs.set_xlabel(r'x ($\AA$)', size=fs)
axs.set_ylabel(r'y ($\AA$)', size=fs)
axs.text(0.90, np.max(Y)*0.529*0.98, '%0.2f V' %1.200, color='white', size=fs*2, ha='left', va='top')

fig.tight_layout()
fig.savefig("Topo-8A.png", dpi=400)

# plt.close(fig)
