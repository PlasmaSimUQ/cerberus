
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


# ==============================================================================
# MAKE MOVIES
# ==============================================================================

def get_Dx(ds, c):
    x, Dx = ds.get("x_D-field")
    return {"x":x[0], "y":x[1], "value":Dx}

def get_Dy(ds, c):
    x, Dy = ds.get("y_D-field")
    return {"x":x[0], "y":x[1], "value":Dy}

def get_vfrac(ds, c):
    x, v = ds.get("vfrac-field")
    return {"x":x[0], "y":x[1], "value":v}


def plot(frame, data, output_name):

    xc = data["vfrac"]["x"]
    yc = data["vfrac"]["y"]
    v = data["vfrac"]["value"][()]
    
    Dx = data["Dx"]["value"][()]
    Dy = data["Dy"]["value"][()]
    vmin = 0.0
    vmax = frame["Dx"]["max"]

    D = np.sqrt(Dx**2 + Dy**2)

    limits = frame["q"]["xy_limits"]

    mask = D < 0.01

    Dx = np.ma.masked_where(mask, Dx)
    Dy = np.ma.masked_where(mask, Dy)

    Dx /= D
    Dy /= D


    fig = plt.figure()
    ax = fig.add_subplot(111)

    # yn, xn = np.meshgrid(yn, xn)

    scale = 0.25

    norm = mpl.colors.Normalize(vmin=vmin, vmax=scale*vmax)
    im = NonUniformImage(ax, interpolation='bilinear', extent=(limits[0][0], limits[1][0], limits[0][1], limits[1][1]),
        norm=norm,
        cmap="viridis")

    im.set_data(xc, yc, D.T)
    ax.images.append(im)

    # ax.pcolormesh(xn, yn, a, vmin=vmin, vmax=vmax)

    cs = ax.contour(xc, yc, v.T, levels=[0.5], colors=['k'])

    skip = 5
    ax.quiver(xc[0::skip], yc[0::skip], Dx.T[0::skip,0::skip], Dy.T[0::skip,0::skip], 
        scale=60.0, scale_units='width', pivot='mid', width=0.002)#, headwidth=0.02)


    ax.text(0.95, 0.05, r'$\left| \vec{D} \right|$', horizontalalignment='right', 
        verticalalignment='bottom', transform=ax.transAxes, fontsize=18, color='w')

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
    {"func":get_Dx, "tag":"Dx"},
    {"func":get_Dy, "tag":"Dy"},
    {"func":get_vfrac, "tag":"vfrac"},
]

# how to make a frame
q["plot"] = plot
q["name"] = "movie"

##
q["framerate"] = 15
q["mov_save"] = q["files_dir"] + "/mov"
q["offset"] = [0.0, 0.0]
q["xy_limits"] = [[-8.0, -8.0], [8.0, 8.0]]
q["file_include"] = ["plt"]
q["file_exclude"] = ["chk"]
q["cores"] = 8
q["time_span"] = [] #np.arange(0,0.1,0.005).tolist()
q["force_data"] = False
q["force_frames"] = True
q["only_frames"] = False
q["redo_streaks"] = False
q["dpi"] = 300

q["normalize"] = "all"

Q.append(q)

make_all(Q)

print("DONE")
