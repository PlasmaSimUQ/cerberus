
import sys  # nopep8
cmd_folder = "../../../vis"  # nopep8
if cmd_folder not in sys.path:  # nopep8
    sys.path.insert(0, cmd_folder)

from tile_mov import tile_movie
from make_mov import make_all, get_particle_trajectories
import numpy as np
import pylab as plt
import matplotlib as mpl
from matplotlib.image import NonUniformImage
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle


# ==============================================================================
# MAKE MOVIES
# ==============================================================================

def get_velocity_magnitude(ds, c):
    x, u = ds.get("x_vel-air")
    x, v = ds.get("y_vel-air", grid='node')

    return {"x":x[0], "y":x[1], "value":np.sqrt(u**2 + v**2)}

def get_alpha(ds, c):
    x, a = ds.get("alpha-air", grid='node')
    return {"x":x[0], "y":x[1], "value":a}

def get_vfrac(ds, c):
    x, v = ds.get("vfrac-air")
    return {"x":x[0], "y":x[1], "value":v}

def get_particles(ds, c):
    idat, rdat =  ds.get_particles('air')
    return {"i":idat, "r":rdat}

def get_boxes(ds, c):
    boxes = ds.get_boxes()
    return {"boxes":boxes}

def plot(frame, data, output_name):

    xc = data["vfrac-air"]["x"]
    yc = data["vfrac-air"]["y"]
    v = data["vfrac-air"]["value"]
    # boxes = data["boxes"]["boxes"][()]
    
    # xn = data["alpha-air"]["x"]
    # yn = data["alpha-air"]["y"]
    a = data["alpha-air"]["value"]
    vmin = frame["alpha-air"]["min"]
    vmax = frame["alpha-air"]["max"]

    limits = frame["q"]["xy_limits"]

    # particles
    px, py, pi = get_particle_trajectories(data["particles-air"], limits)

    fig = plt.figure()
    ax = fig.add_subplot(111)

    # yn, xn = np.meshgrid(yn, xn)

    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    im = NonUniformImage(ax, interpolation='bilinear', extent=(limits[0][0], limits[1][0], limits[0][1], limits[1][1]),
        norm=norm)

    im.set_data(xc, yc, np.transpose(a))
    ax.images.append(im)

    # ax.pcolormesh(xn, yn, a, vmin=vmin, vmax=vmax)

    cs = ax.contour(xc, yc, np.transpose(v), levels=[0.5], colors=['k'])

    for line in cs.allsegs[0]:
        plt.fill(line[:,0], line[:,1])

    l = 10 #streak length
    ax.plot(px[-l::], py[-l::], "k-", lw=0.5, alpha=0.5)
    ax.plot(px[-1], py[-1], "ko", ms=0.5, alpha=0.5)

    # plot boxes
    # grid = []
    # for box in boxes:
    #     sz = box[1] - box[0]
    #     rect = Rectangle((box[0][0], box[0][1]), sz[0], sz[1])
    #     grid.append(rect)

    # pc = PatchCollection(grid, facecolor='none', alpha=1.0, edgecolor='w', linewidth=0.25)
    # ax.add_collection(pc)

    ax.set_xlim(limits[0][0], limits[1][0])
    ax.set_ylim(limits[0][1], limits[1][1])

    ax.set_aspect(1)
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)
    fig.tight_layout()
    fig.savefig(output_name, dpi=300, bbox_inches="tight")
    plt.close(fig)

    return



dt = 0.005

Q = []

q = {}
q["files_dir"] = "."
q["level"] = -1

# all the data we need to retrieve
q["get"] = [
    {"func":get_alpha, "tag":"alpha-air"},
    {"func":get_vfrac, "tag":"vfrac-air"},
    {"func":get_particles, "tag":"particles-air", "get_streak":True},
    # {"func":get_boxes, "tag":"boxes"},
]

# how to make a frame
q["plot"] = plot
q["name"] = "movie"

##
q["framerate"] = 20
q["mov_save"] = q["files_dir"] + "/mov"
q["offset"] = [0.0, 0.0]
q["xy_limits"] = [[0.0, 0.0], [1.0, 1.0]]
q["file_include"] = ["plt"]
q["file_exclude"] = ["chk"]
q["cores"] = 8
q["time_span"] = [] #np.arange(0,0.1,0.005).tolist()
q["force_data"] = True
q["force_frames"] = True
q["only_frames"] = False
q["redo_streaks"] = False
q["dpi"] = 300

q["normalize"] = "all"

Q.append(q)

make_all(Q)

print("DONE")
