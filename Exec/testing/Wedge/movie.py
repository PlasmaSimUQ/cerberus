
import sys  # nopep8
cmd_folder = "../../../vis"  # nopep8
if cmd_folder not in sys.path:  # nopep8
    sys.path.insert(0, cmd_folder)

from get_boxlib import get_files, get_time
from tile_mov import tile_movie
from make_mov import make_all, get_particle_trajectories
import numpy as np
import pylab as plt
import matplotlib as mpl
from matplotlib.image import NonUniformImage
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle


# ==============================================================================
# MAKE MOVIES
# ==============================================================================

def get_rho(ds, c):
    x, v = ds.get("vfrac-fluid")
    x, r = ds.get("rho-fluid")
    r = np.ma.masked_where(v <= 1e-2, r)
    return {"x":x[0], "y":x[1], "value":r}

def get_T(ds, c):
    x, v = ds.get("vfrac-fluid")
    x, T = ds.get("T-fluid")
    T = np.ma.masked_where(v <= 1e-2, T)
    return {"x":x[0], "y":x[1], "value":T}

def get_vfrac(ds, c):
    x, v = ds.get("vfrac-fluid")
    return {"x":x[0], "y":x[1], "value":v}

def get_boxes(ds, c):
    boxes = ds.get_boxes()
    return {"boxes":boxes}


def plot(frame, data, output_name):

    xc = data["vfrac"]["x"]
    yc = data["vfrac"]["y"]
    v = data["vfrac"]["value"][()]
    rho = data["rho"]["value"][()]
    T = data["T"]["value"][()][()]
    boxes = data["boxes"]["boxes"][()]
    
    rho_vmin = frame["rho"]["min"]
    rho_vmax = 5.5 #frame["rho"]["max"]

    T_vmin = frame["T"]["min"]
    T_vmax = 6.0 #frame["T"]["max"]

    limits = frame["q"]["xy_limits"]


    fig = plt.figure()
    axes = []

    ###
    ax = fig.add_subplot(121); axes.append(ax)

    norm = mpl.colors.Normalize(vmin=rho_vmin, vmax=rho_vmax)
    im = NonUniformImage(ax, interpolation='bilinear', extent=(limits[0][0], limits[1][0], limits[0][1], limits[1][1]),
        norm=norm,
        cmap="viridis")

    im.set_data(xc, yc, rho.T)
    ax.images.append(im)

    cs = ax.contour(xc, yc, v.T, levels=[0.5], colors=['k'], linewidths=[0.25])
    line_x = cs.allsegs[0][0][:,0].tolist()
    line_y = cs.allsegs[0][0][:,1].tolist()
    plt.fill(line_x+[limits[1][0]], line_y+[limits[0][1]], color='k')

    ax.text(0.05, 0.05, r'$\rho$', horizontalalignment='left', 
        verticalalignment='bottom', transform=ax.transAxes, fontsize=18, color='w')

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax)

    ###
    ax = fig.add_subplot(122); axes.append(ax)

    norm = mpl.colors.Normalize(vmin=T_vmin, vmax=T_vmax)
    im = NonUniformImage(ax, interpolation='bilinear', extent=(limits[0][0], limits[1][0], limits[0][1], limits[1][1]),
        norm=norm,
        cmap="viridis")

    im.set_data(xc, yc, T.T)
    ax.images.append(im)

    cs = ax.contour(xc, yc, v.T, levels=[0.5], colors=['k'], linewidths=[0.25])
    line_x = cs.allsegs[0][0][:,0].tolist()
    line_y = cs.allsegs[0][0][:,1].tolist()
    plt.fill(line_x+[limits[1][0]], line_y+[limits[0][1]], color='k')

    # plot boxes
    grid = []
    for box in boxes:
        sz = box[1] - box[0]
        rect = Rectangle((box[0][0], box[0][1]), sz[0], sz[1])
        grid.append(rect)

    pc = PatchCollection(grid, facecolor='none', alpha=1.0, edgecolor='w', linewidth=0.25)
    ax.add_collection(pc)

    ax.text(0.05, 0.05, r'$T$', horizontalalignment='left', 
        verticalalignment='bottom', transform=ax.transAxes, fontsize=18, color='w')

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax)

    for ax in axes:
        ax.set_xlim(limits[0][0], limits[1][0])
        ax.set_ylim(limits[0][1], limits[1][1])
        ax.set_aspect(1)
        ax.axes.xaxis.set_visible(False)
        ax.axes.yaxis.set_visible(False)

    fig.tight_layout()
    fig.savefig(output_name, dpi=300, bbox_inches="tight")
    plt.close(fig)

    return


Q = []

q = {}
q["files_dir"] = "."
q["level"] = -1

# all the data we need to retrieve
q["get"] = [
    {"func":get_rho, "tag":"rho"},
    {"func":get_T, "tag":"T"},
    {"func":get_vfrac, "tag":"vfrac"},
    {"func":get_boxes, "tag":"boxes"},
]

# how to make a frame
q["plot"] = plot
q["name"] = "movie"



##
q["framerate"] = 30
q["mov_save"] = q["files_dir"] + "/mov"
q["offset"] = [0.0, 0.0]
q["xy_limits"] = [[-0.009, 0.0], [0.001, 0.01]]
q["file_include"] = ["plt"]
q["file_exclude"] = ["chk"]
q["cores"] = 11
q["force_data"] = False
q["force_frames"] = True
q["only_frames"] = False
q["redo_streaks"] = False
q["dpi"] = 300

q["normalize"] = "all"

# non-linear time sampling
files = get_files(q["files_dir"], include=q["file_include"], exclude=q["file_exclude"], get_all=True)
n_files = len(files)
times = []
for f in files:
    times.append(get_time(f))
times = np.array(times)

log_i = np.logspace(0, np.log10(n_files-1), num=int(n_files/2), dtype=int)
log_i = list(dict.fromkeys(log_i))
q["time_span"] = times[log_i].tolist()

Q.append(q)

make_all(Q)

print("DONE")
