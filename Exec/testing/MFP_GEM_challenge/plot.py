import sys

cmd_folder = "../../../vis"
if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)

import pdb

import matplotlib.ticker as ticker
import numpy as np
import pylab as plt
from get_boxlib import ReadBoxLib, get_files
from matplotlib.image import NonUniformImage

# ==============================================================================
#
# ==============================================================================


# get a list of all the files in this directory
files = get_files(".", include=["plt"])

print(files)

N = 5

for i in range(N):
    fig = plt.figure(figsize=(10, 5))

    # get data
    f = files[i * int(len(files) / float(N - 1))]

    print(f"Get data for output file:\t{f}") 

    data = ReadBoxLib(f)
    t = data.time

    data = ReadBoxLib(f, max_level=-1)

    # -->>>> get any other variables you are interested in from here <<<<--#
    """
['T-electrons', 'T-ions', 'charge-electrons', 'charge-ions', 'cost', 'cp-electrons', 'cp-ions', 'ep-field', 'gamma-electrons', 'gamma-ions', 'mass-electrons', 'mass-ions', 'mu-field', 'nrg-electrons', 'nrg-ions', 'p-electrons', 'p-ions', 'phi-field', 'psi-field', 'rho-electrons', 'rho-ions', 'x_B-field', 'x_D-field', 'x_mom-electrons', 'x_mom-ions', 'x_vel-electrons', 'x_vel-ions', 'y_B-field', 'y_D-field', 'y_mom-electrons', 'y_mom-ions', 'y_vel-electrons', 'y_vel-ions', 'z_B-field', 'z_D-field', 'z_mom-electrons', 'z_mom-ions', 'z_vel-electrons', 'z_vel-ions']
    """
    #data format ( list of vectors holding spatial coordinates, and array of requested value
    # xc, vf = data.get("vfrac-field", grid="cell")
    xn, xD = data.get("x_D-field", grid="node")

    xn, xB = data.get("x_B-field", grid="node")
    xn, yB = data.get("y_B-field", grid="node")
    xn, zB = data.get("z_B-field", grid="node")

    B_mag = np.sqrt( xB*xB + yB*yB + zB*zB ) 

    # set a trace (otherwise known as a breakpoint) in the code at this exact point
    # pdb.set_trace()

    yn, xn = np.meshgrid(xn[1], xn[0])
    #yc, xc = np.meshgrid(xc[1], xc[0])

    # plot stuff
    total_panels = 2 
    # subpanel index, panel_index
    panel_index = 1
    ax = fig.add_subplot(1, total_panels, panel_index)
    pc = ax.pcolormesh(xn, yn, xD)
    plt.colorbar(pc, orientation="horizontal")
    ax.set_aspect(1)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(r"$D_x$ @ $t=%g$" % t)

    panel_index = 2
    ax = fig.add_subplot(1, total_panels, panel_index)
    pc = ax.pcolormesh(xn, yn, B_mag)
    plt.colorbar(pc, orientation="horizontal")
    ax.set_aspect(1)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(r"$|B|$ @ $t=%g$" % t)

    fig.tight_layout()

    fig.savefig(f"plot_t-{t}.png", dpi=300)


print("DONE")
