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

def get_number_density(ds, c):
    x, r = ds.get("rho-%s"%c["component"])
    x, m = ds.get("mass-%s"%c["component"], grid='node')

    return {"x":x[0], "y":x[1], "value":r/m}

def get_vfrac(ds, c):
    x, vf = ds.get("vfrac-%s"%c["component"])

    return {"x":x[0], "y":x[1], "value":vf}

def plot(frame, data, output_name):

    xn = data["nd-ion1"]["x"]
    yn = data["nd-ion1"]["y"]

    ni1 = data["nd-ion1"]["value"]
    ni1_min = frame["nd-ion1"]["min"]
    ni1_max = frame["nd-ion1"]["max"]

    ne1 = data["nd-electron1"]["value"]
    ne1_min = frame["nd-electron1"]["min"]
    ne1_max = frame["nd-electron1"]["max"]

    ni2 = data["nd-ion2"]["value"]
    ni2_min = frame["nd-ion2"]["min"]
    ni2_max = frame["nd-ion2"]["max"]

    ne2 = data["nd-electron2"]["value"]
    ne2_min = frame["nd-electron2"]["min"]
    ne2_max = frame["nd-electron2"]["max"]

    vf = data["vfrac-fluid"]["value"]
    xc = data["vfrac-fluid"]["x"]
    yc = data["vfrac-fluid"]["y"]

    limits = frame["q"]["xy_limits"]

    yn, xn = np.meshgrid(yn, xn)

    yc, xc = np.meshgrid(yc, xc)

    axes = []

    fig = plt.figure(figsize=(8.5,5))

    gs = gridspec.GridSpec(nrows=2, ncols=2, wspace=0.05, hspace=0.05)

    ax = fig.add_subplot(gs[0,0]); axes.append(ax)
    pcm = ax.pcolormesh(xn, yn, ni1, vmin=ni1_min, vmax=ni1_max)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='1.5%', pad=0.05)
    fig.colorbar(pcm, cax=cax, label=r"$n_i^{j=0}$")


    ax = fig.add_subplot(gs[0,1]); axes.append(ax)
    pcm = ax.pcolormesh(xn, yn, ne1, vmin=ne1_min, vmax=ne1_max)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='1.5%', pad=0.05)
    fig.colorbar(pcm, cax=cax, label=r"$n_e^{j=0}$")


    ax = fig.add_subplot(gs[1,0]); axes.append(ax)
    pcm = ax.pcolormesh(xn, yn, ni2, vmin=ni2_min, vmax=ni2_max)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='1.5%', pad=0.05)
    fig.colorbar(pcm, cax=cax, label=r"$n_i^{j\neq 0}$")


    ax = fig.add_subplot(gs[1,1]); axes.append(ax)
    pcm = ax.pcolormesh(xn, yn, ne2, vmin=ne2_min, vmax=ne2_max)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='1.5%', pad=0.05)
    fig.colorbar(pcm, cax=cax, label=r"$n_e^{j\neq 0}$")

    for ax in axes:

        cs = ax.contour(xc, yc, vf, levels=[0.5], linewidths=[0.5], colors=['k'])
        for line in cs.allsegs[0]:

            x = line[:,0].tolist()
            y = line[:,1].tolist()

            ax.fill(x, y, 'w', hatch='////')

        ax.set_xlim(limits[0][0], limits[1][0])
        ax.set_ylim(limits[0][1], limits[1][1])

        ax.set_aspect(1)
        ax.axes.xaxis.set_visible(False)
        ax.axes.yaxis.set_visible(False)

    fig.savefig(output_name, dpi=300, bbox_inches="tight")
    plt.close(fig)

    return

if 1:
    dt = 0.005

    Q = []

    q = {}
    q["files_dir"] = "."
    q["level"] = -1


    q["get"] = [
        {"func":get_number_density, "tag":"nd-ion1", "component":"ion1"},
        {"func":get_number_density, "tag":"nd-electron1", "component":"electron1"},

        {"func":get_number_density, "tag":"nd-ion2", "component":"ion2"},
        {"func":get_number_density, "tag":"nd-electron2", "component":"electron2"},

        {"func":get_vfrac, "tag":"vfrac-fluid", "component":"ion1"},
    ]

    q["plot"] = plot
    q["name"] = "movie"
    

    ##
    q["framerate"] = 60
    q["mov_save"] = q["files_dir"] + "/mov"
    q["offset"] = [0.0, 0.0]
    q["xy_limits"] = [[-5,-5], [10,5]]
    q["file_include"] = ["wire.plt"]
    q["file_exclude"] = []
    q["cores"] = 10
    q["time_span"] = [] #np.arange(0, 0.1, 0.005)
    q["force_data"] = False
    q["force_frames"] = False
    q["only_frames"] = False
    q["redo_streaks"] = False
    q["dpi"] = 200

    q["normalize"] = "all"

    Q.append(q)

    make_all(Q)

print("DONE")
