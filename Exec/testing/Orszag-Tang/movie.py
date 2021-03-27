
import sys  # nopep8
cmd_folder = "../../../vis"  # nopep8
if cmd_folder not in sys.path:  # nopep8
    sys.path.insert(0, cmd_folder)

from tile_mov import tile_movie
from make_mov import make_all, get_particle_trajectories
import numpy as np
import pylab as plt


# ==============================================================================
# MAKE MOVIES
# ==============================================================================

def get_velocity_magnitude(ds, c):
    x, u = ds.get("x_vel-mhd")
    x, v = ds.get("y_vel-mhd", grid='node')

    return {"x":x[0], "y":x[1], "value":np.sqrt(u**2 + v**2)}


def get_particles(ds, c):
    idat, rdat =  ds.get_particles('mhd')
    return {"i":idat, "r":rdat}

def plot(frame, data, output_name):

    dat = data["velocity"]
    xn = dat["x"]
    yn = dat["y"]
    vel = dat["value"]
    vmin = frame["velocity"]["min"]
    vmax = frame["velocity"]["max"]

    limits = frame["q"]["xy_limits"]

    # particles
    px, py, pi = get_particle_trajectories(data["particles"], limits)

    fig = plt.figure()
    ax = fig.add_subplot(111)

    yn, xn = np.meshgrid(yn, xn)

    ax.pcolormesh(xn, yn, vel, vmin=vmin, vmax=vmax)

    l = 10 #streak length
    ax.plot(px[-l::], py[-l::], "k-", lw=0.5, alpha=0.5)
    ax.plot(px[-1], py[-1], "ko", ms=0.5, alpha=0.5)

    ax.text(0.05, 0.05, r'$\left| \vec{u} \right|$', horizontalalignment='left', 
        verticalalignment='bottom', transform=ax.transAxes, fontsize=18)

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
    {"func":get_velocity_magnitude, "tag":"velocity"},
    {"func":get_particles, "tag":"particles", "get_streak":True},
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
q["cores"] = 4
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
