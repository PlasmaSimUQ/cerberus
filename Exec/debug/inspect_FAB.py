# -*- coding: utf-8 -*-
"""
Created on Thu Sep 20 05:46:07 2018

@author: uqdbond1
"""

import numpy as np
import re
import pylab as plt
import os, sys

np.set_printoptions(linewidth=1000)

#==============================================================================
# 
#==============================================================================

def multi_slice_viewer(volume, ax, error_limit=0):

    ax.volume = volume
    data = np.ma.masked_invalid(volume)
    ax.index = volume.shape[-1] / 2
    if error_limit:
        vmax = error_limit
    else:
        vmax = data.max()
    im = ax.imshow(volume[:,:,ax.index].T, origin="lower", aspect='equal',
                   vmin=data.min(), vmax=vmax)
    im.cmap.set_bad('magenta')
    if error_limit:
        im.cmap.set_over('red')
    
    hi = volume.shape
    
    # Minor ticks
    ax.set_xticks(np.arange(-.5, hi[0]+0.5, 1), minor=True);
    ax.set_yticks(np.arange(-.5, hi[1]+0.5, 1), minor=True);
    
    # Gridlines based on minor ticks
    ax.grid(which='minor', color='w', linestyle='-', linewidth=0.5)
    
    ax.figure.canvas.mpl_connect('scroll_event', process_key)
    plt.colorbar(im, orientation="horizontal")
    
    return ax

def process_key(event):
    fig = event.canvas.figure
    ax = fig.axes[0]
    if event.button == 'up':
        previous_slice(ax)
    elif event.button == 'down':
        next_slice(ax)
    fig.canvas.draw()

def previous_slice(ax):
    """Go to the previous slice."""
    volume = ax.volume
    ax.index = (ax.index - 1) % volume.shape[-1]  # wrap around using %
    ax.images[0].set_array(volume[:,:,ax.index].T)
    print ax.index

def next_slice(ax):
    """Go to the next slice."""
    volume = ax.volume
    ax.index = (ax.index + 1) % volume.shape[-1]
    ax.images[0].set_array(volume[:,:,ax.index].T)
    print ax.index
    
def get_data(fname):
    f = open(fname, "r")
    
    lines = f.readlines()
    
    # get the size of the data
    dims = np.zeros((6), dtype=int)
    
    lo_str = lines[0].split(",")
    hi_str = lines[1].split(",")
    
    for i in range(3):
        dims[i] = int(lo_str[i])
        dims[i+3] = int(hi_str[i])

    lo = dims[0:3]
    size = dims[3:6] - dims[0:3] + 1
    
    data = np.zeros(size)*np.NaN
    
    for line in lines[2::]:
        num = line.split(',')
        
        i = int(num[0]) - lo[0]
        j = int(num[1]) - lo[1]
        k = int(num[2]) - lo[2]
        f = float(num[3])        
        
        data[i,j,k] = f
    
    data = np.ma.masked_invalid(data)
    
    return data, dims
    
def combine_data(folder, fname):

    fid = open(os.path.join(folder, fname + ".fab"))
    files = fid.readlines()

    for i,f in enumerate(files):
        files[i] = f.strip()
    
    data = []
    large = sys.maxint
    lo = np.ones((3),dtype=int)*large
    hi = -np.ones((3),dtype=int)*large
    
    check = np.empty((3,2), dtype=int)
    
    #idf.append("fab")
    
    for f in files:
#        print f
        a, b = get_data(f)
        dat = {"data":a, "limits":b}
        data.append(dat)

        check[:,0] = lo[:]
        check[:,1] = b[0:3]
        lo = np.min(check, axis=1)

        check[:,0] = hi
        check[:,1] = b[3:6]
        hi = np.max(check, axis=1)
    
    size = hi - lo + 1
    
    large = np.ones((size))*np.nan
    
    
    
    for dat in data:
        lim = dat["limits"]
        a = lim[0:3] - lo
        b = lim[3:6] - lo

        for i in range(a[0],b[0]+1):
            for j in range(a[1],b[1]+1):
                for k in range(a[2],b[2]+1):
                    d = dat["data"][i-a[0], j-a[1], k-a[2]]
                    if np.isfinite(d):
                        large[i,j,k] = d
        
        #large[a[0]:b[0]+1, a[1]:b[1]+1, a[2]:b[2]+1] = dat["data"]
        
        
    large = np.ma.masked_invalid(large)
        
        
#    print lo, hi
        
    
    return large, lo, hi
    
    
        

#==============================================================================
    
cutoff = 0.0
    
if 1:
    fig = plt.figure(figsize=(9,9))
    
    ax = fig.add_subplot(1,3,1)
    dat, lo, hi = combine_data("debug", "original_psi_E_c0")
    multi_slice_viewer(dat, ax)
    
    ax = fig.add_subplot(1,3,2)
    dat, lo, hi = combine_data("debug", "filled_psi_E_c0")
    multi_slice_viewer(dat, ax)
    
    ax = fig.add_subplot(1,3,3)
    dat, lo, hi = combine_data("debug", "dudt_psi_E_c0")
    multi_slice_viewer(dat, ax)
    
    plt.show()

