import sys
cmd_folder = "../../../vis"  # nopep8
if cmd_folder not in sys.path:  # nopep8
    sys.path.insert(0, cmd_folder)
from tile_mov import tile_movie
from make_mov import make_all, get_particle_trajectories
import pylab as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gridspec

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica"]})

# ==============================================================================
# MAKE MOVIES
# ==============================================================================

def smooth_limits(vmin, vmax):

    from scipy.signal import savgol_filter

    vmin = savgol_filter(vmin, 11, 3)
    vmax = savgol_filter(vmax, 11, 3)

    return vmin, vmax

def get_number_density(ds, c):
    x, r = ds.get("rho-%s"%c["component"])
    x, m = ds.get("mass-%s"%c["component"], grid='node')

    return {"x":x[0], "y":x[1], "value":r/m}

def get_D_mag(ds, c):
    x, Dx = ds.get("x_D-field")
    x, Dy = ds.get("y_D-field")

    return {"x":x[0], "y":x[1], "value":np.sqrt(Dx**2 + Dy**2)}

def get_Bz(ds, c):
    x, Bz = ds.get("z_B-field")

    return {"x":x[0], "y":x[1], "value":Bz}

def plot(frame, data, output_name):

    xn = data["nd-ion"]["x"][()]
    yn = data["nd-ion"]["y"][()]

    ni = data["nd-ion"]["value"][()]
    ni_min = frame["nd-ion"]["min"]
    ni_max = frame["nd-ion"]["max"]

    ne = data["nd-electron"]["value"][()]
    ne_min = frame["nd-electron"]["min"]
    ne_max = frame["nd-electron"]["max"]

    D = data["D"]["value"][()]
    D_min = frame["D"]["min"]
    D_max = frame["D"]["max"]

    B = data["B"]["value"][()]
    B_min = frame["B"]["min"]
    B_max = frame["B"]["max"]


    x = np.concatenate((-xn[::-1][0:-1], xn))
    y = np.concatenate((-yn[::-1][0:-1], yn))

    y, x = np.meshgrid(y, x)

    axes = []


    # join the data
    nx = xn.size - 1
    ny = yn.size - 1

    fig = plt.figure(figsize=(3,3))
    gs = gridspec.GridSpec(ncols=1, nrows=1, hspace=0.01, wspace=0.01)
    ax = fig.add_subplot(gs[0,0]); axes.append(ax)

    # number densities
    J = np.zeros((2*nx, 2*ny))*np.nan

    # J[0:nx, 0:ny] = np.rot90(ne.T,2)
    J[0:nx, ny::] = np.rot90(ne)

    # J[nx::, 0:ny] = np.rot90(ni.T,3)
    J[nx::, ny::] = ni

    vmin = min(ne_min, ni_min)
    vmax = max(ne_max, ni_max)
    pcm = ax.pcolormesh(x, y, J, vmin=vmin, vmax=vmax)

    ax.text(0.025, 0.975, r'$n_e$', horizontalalignment='left', 
    verticalalignment='top', transform=ax.transAxes, fontsize=10)

    ax.text(0.975, 0.975, r'$n_i$', horizontalalignment='right', 
    verticalalignment='top', transform=ax.transAxes, fontsize=10)

    # fields
    J  = np.zeros((2*nx, 2*ny))*np.nan
    J[0:nx, 0:ny] = np.rot90(D.T,2)
    pcm = ax.pcolormesh(x, y, J, vmin=D_min, vmax=D_max)

    J  = np.zeros((2*nx, 2*ny))*np.nan
    J[nx::, 0:ny] = np.rot90(B.T,3)
    big = max(abs(B_max), abs(B_min))
    pcm = ax.pcolormesh(x, y, J, vmin=-big, vmax=big, cmap="bwr")

    ax.text(0.025, 0.025, r'$\left|\vec{D}\right|$', horizontalalignment='left', 
    verticalalignment='bottom', transform=ax.transAxes, fontsize=10)

    ax.text(0.975, 0.025, r'$B_z$', horizontalalignment='right', 
    verticalalignment='bottom', transform=ax.transAxes, fontsize=10)

    for ax in axes:
        ax.set_xlim(-2, 2)
        ax.set_ylim(-2, 2)

        ax.set_aspect(1)
        ax.axes.xaxis.set_visible(False)
        ax.axes.yaxis.set_visible(False)

    # fig.tight_layout()
    fig.savefig(output_name, dpi=300, bbox_inches="tight")
    plt.close(fig)

    return

if 1:

    Q = []

    q = {}
    q["files_dir"] = "."
    q["level"] = -1


    q["get"] = [
        {"func":get_number_density, "tag":"nd-ion", "component":"ion"},
        {"func":get_number_density, "tag":"nd-electron", "component":"electron"},
        {"func":get_D_mag, "tag":"D"},
        {"func":get_Bz, "tag":"B"}
    ]

    q["plot"] = plot
    q["name"] = "movie"
    
    dt = 0.005

    ##
    q["framerate"] = 20
    q["mov_save"] = q["files_dir"] + "/mov"
    q["offset"] = [0.0, 0.0]
    q["xy_limits"] = [[0,0], [4,4]]
    q["file_include"] = ["TRMI.plt"]
    q["file_exclude"] = []
    q["cores"] = 11
    q["time_span"] = [] #np.arange(1.95,2+dt, dt)
    q["force_data"] = False
    q["force_frames"] = True
    q["only_frames"] = False
    q["redo_streaks"] = False
    q["dpi"] = 300

    q["normalize"] = "none" #{"smooth":smooth_limits}

    Q.append(q)

    make_all(Q)

print("DONE")
