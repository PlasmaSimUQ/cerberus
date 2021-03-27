
import sys
cmd_folder = "../../../vis"
if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)
    
from get_boxlib import ReadBoxLib, get_files

import numpy as np
import pylab as plt
import matplotlib as mpl
from matplotlib.image import NonUniformImage
from multiprocessing import Pool
import h5py
import matplotlib.gridspec as gridspec
from scipy.interpolate import RectBivariateSpline
from mpl_toolkits.axes_grid1 import make_axes_locatable


#==============================================================================
# 
#==============================================================================

def check():

    plt_file = "Orszag-Tang"
    component = "mhd"
    window = [[0, 0], [1,1]]

    # get a list of all the files in this directory
    files = get_files('.', include=[plt_file], exclude=["input", "temp", ".png"], get_all=True)
    # get tracer particle data

    # get fluid variables

    rh5 = ReadBoxLib(files[-1], max_level=-1, limits=window)

    xn, z = rh5.get("p-%s"%component, grid='node')

    xn, yn = xn
    
    x = xn[0:-1]+(xn[1] + xn[0])/2.0
    y = yn[0:-1]+(yn[1] + yn[0])/2.0
    y, x = np.meshgrid(y, x)

    # =============================================================================
    # 
    # =============================================================================

    # read in reference values
    # HIGH-ORDER UPWIND SCHEMES FOR MULTIDIMENSIONAL MAGNETOHYDRODYNAMICS
    # THE ASTROPHYSICAL JOURNAL, 530 : 508-524, 2000 February 10

    Londrillo_1D_hi = h5py.File("y=0-4277_t=0-5.hdf5", "r")


    # =============================================================================
    # 
    # =============================================================================

    plt.rc("font", family="serif")
    plt.rc("font", size=8)
    plt.rc("mathtext", fontset="cm")
    # matplotlib.rc('text', usetex = True)
    params= {'text.latex.preamble' : [r'\usepackage{amsmath}']}
    plt.rcParams.update(params)

    axes = []

    fig = plt.figure(figsize=(5,2.5))

    gs = gridspec.GridSpec(nrows=1, ncols=3, width_ratios=[0.6, 1, 0.05], wspace=0.01, bottom=0.14, top=0.97, left=0.1, right=0.88)

    #====

    ax = fig.add_subplot(gs[0,1]); axes.append(ax)
    ax_cb = fig.add_subplot(gs[0,2])

    # plot the contour
    im = NonUniformImage(ax, interpolation='bilinear', extent=[xn[0], yn[0], xn[-1], yn[-1]],
                        cmap="viridis")
    im.set_data(x[:,0], y[0,:], z.T)
    ax.images.append(im)

    # divider = make_axes_locatable(ax)
    # ax_cb = divider.new_horizontal(size="5%", pad=0.05)
    plt.colorbar(im, cax=ax_cb, label=r"$p$")

    ax.plot([0,1],2*[0.4277], "r--", lw=1)

    ax.set_aspect(1)
    ax.set_xlim(0,1)
    ax.set_ylim(0,1)
    ax.set_xticks([])
    ax.set_yticks([])

    #====

    xx = x[:,0]
    yy = y[0,:]

    #====

    ax = fig.add_subplot(gs[0,0]); axes.append(ax)

    # plot reference points
    cx = []
    cy = []
    for key, value in Londrillo_1D_hi.items():
        if "Line" not in key:
            continue
        if "CX" not in value:
            continue

        cx.append(value["CX"][()])
        cy.append(value["CY"][()])

    cx = np.array(cx)
    cy = np.array(cy)

    I = np.argsort(cx)
    cx = cx[I]
    cy = cy[I]

    ax.plot(cx, cy, "-ok", lw=0.25, ms=2.0, mfc='none', mew=0.3)

    p = RectBivariateSpline(xx, yy, z)
    sy = p(cx, cx.size*[0.4277], grid=False)

    ax.plot(cx, sy, 'k-', lw=0.5)

    rms = np.sqrt(np.mean(np.square(cy - sy)))
    # print("RMS = ",rms)


    ax.set_xlim(0,1)
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$p$")

    #====

    fig.savefig(plt_file+".png", dpi=300)
    plt.close(fig)

    if rms > 0.05:
        return 1
    else:
        return 0

if __name__ == "__main__":
    sys.exit(check())