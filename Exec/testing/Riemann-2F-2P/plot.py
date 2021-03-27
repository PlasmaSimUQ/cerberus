
import sys
cmd_folder = "../../../vis"
if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)
    
from get_boxlib import ReadBoxLib, get_files

import numpy as np
import pylab as plt
from matplotlib.image import NonUniformImage
import matplotlib.ticker as ticker
import h5py

#==============================================================================
# 
#==============================================================================

plt.rc("font", family="serif")
plt.rc("font", size=8)
plt.rc("mathtext", fontset="cm")
params= {'text.latex.preamble' : [r'\usepackage{amsmath}']}
plt.rcParams.update(params)

# =============================================================================
# 
# =============================================================================

rho_ref = "Loverich-1-rho.hdf5"
ux_ref  = "Loverich-1-ux.hdf5"
uy_ref  = "Loverich-1-uy.hdf5"
uz_ref  = "Loverich-1-uz.hdf5"
p_ref   = "Loverich-1-p.hdf5"
By_ref  = "Loverich-1-By.hdf5"
Bz_ref  = "Loverich-1-Bz.hdf5"

def ref_data(ax, fname, label=False):
    h5 = h5py.File(fname,'r')
    
    for item in h5.items():
        if "Line" in item[0]:
            X = item[1]['X'][()]
            Y = item[1]['Y'][()]
            
            if X.size > 8:
                if label:
                    set_label = label
                else:
                     set_label = None
                     
                ax.plot(X, Y,'ko', lw=0.3 ,ms=2.0, mew=0.3, mfc='none', markevery=20, label = set_label)
                # ax.plot(X, Y,'r--', lw=0.5, label = set_label)

def sim_data(ax, x, y, label=None):
    ax.plot(x, y,'k-',lw=0.75, label=label)

def mhd_data(ax, x, y, label=None):
    ax.plot(x, y,'k-.',lw=0.5, label=label)

# =============================================================================
# 
# =============================================================================




# get a list of all the files in this directory
files = get_files('.', include=['.plt'], exclude=["temp"], times=[], tol=1e-4, get_all=True)

for f in files:

    rh5 = ReadBoxLib(f)
    t = rh5.time
        
    rh5 = ReadBoxLib(f, max_level=-1)

    # standard two fluid

    x, rho_e, ratio = rh5.get('rho-electrons', get_refinement=True)
    x, ux_e  = rh5.get("x_vel-electrons")
    x, p_e   =  rh5.get("p-electrons")

    x, rho_i = rh5.get("rho-ions")
    x, ux_i  = rh5.get("x_vel-ions")
    x, p_i   = rh5.get("p-ions")

    x, By    = rh5.get("y_B-field")

    # two pressure

    x, rho_e_2p, ratio = rh5.get('rho-electrons2p', get_refinement=True)
    x, ux_e_2p  = rh5.get("x_vel-electrons2p")
    x, p_e_2p   =  rh5.get("p-electrons2p")
    x, pp_e_2p   =  rh5.get("pp-electrons2p")

    x, rho_i_2p = rh5.get("rho-ions2p")
    x, ux_i_2p  = rh5.get("x_vel-ions2p")
    x, p_i_2p   = rh5.get("p-ions2p")
    x, pp_i_2p   = rh5.get("pp-ions2p")

    x, By_2p    = rh5.get("y_B-field2p")

    # MHD

    # x, rho_mhd = rh5.get("rho-mhd")
    # x, ux_mhd  = rh5.get("x_vel-mhd")
    # x, p_mhd   =  rh5.get("p-mhd")
    # x, By_mhd    = rh5.get("y_B-mhd")



    axes = []
    fig = plt.figure(figsize=(5,4))

    nr = 4
    nc = 1

    x = x[0]

    ax = fig.add_subplot(nr,nc,1); axes.append(ax)
    sim_data(ax, x, rho_e + rho_i, r"2F")
    ax.plot(x, rho_e_2p + rho_i_2p, "r-.", label="2P")
    # ref_data(ax, rho_ref, "Loverich")
    # mhd_data(ax, x, rho_mhd, "MHD")
    ax.set_ylabel(r"$\rho$")
    ax.legend()

    # ax = ax.twinx()
    # ax.plot(x, ratio, 'k-.', lw=0.25)
    # ax.set_yticks([0,2,4])
    # ax.set_ylabel(r"level")

    ax = fig.add_subplot(nr,nc,2); axes.append(ax)
    sim_data(ax, x, (rho_i*ux_i + rho_e*ux_e)/(rho_i + rho_e))
    ax.plot(x, (rho_i_2p*ux_i_2p + rho_e_2p*ux_e_2p)/(rho_i_2p + rho_e_2p), "r-")
    # ref_data(ax, ux_ref, "Loverich")
    # mhd_data(ax, x, ux_mhd)
    ax.set_ylabel(r"$u_x$")

    ax = fig.add_subplot(nr,nc,3); axes.append(ax)
    sim_data(ax, x, p_i + p_e)
    ax.plot(x, p_i + p_e, "r-.")
    ax.plot(x, pp_i_2p + pp_e_2p, "r--")
    # ref_data(ax, p_ref)
    # mhd_data(ax, x, p_mhd)
    ax.set_ylabel(r"$p$")

    ax = fig.add_subplot(nr,nc,4); axes.append(ax)
    sim_data(ax, x, By,)
    ax.plot(x, By_2p, "r-.")
    # ref_data(ax, By_ref)
    # mhd_data(ax, x, By_mhd)
    ax.set_ylabel(r"$B_y$")


    for ax in axes:
        ax.set_xlim(0.25, 0.75)

    for ax in axes[0:-1]:
        ax.set_xticklabels([])

    axes[-1].set_xlabel(r"$x$")
    fig.tight_layout()

    fig.savefig(f+".png", dpi=300)
    plt.close(fig)
    
    
print("DONE")
