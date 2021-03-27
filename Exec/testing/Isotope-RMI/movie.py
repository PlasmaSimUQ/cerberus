import sys
cmd_folder = "../../../vis"  # nopep8
if cmd_folder not in sys.path:  # nopep8
    sys.path.insert(0, cmd_folder)
from tile_mov import tile_movie
from make_mov import make_all, get_particle_trajectories
import pylab as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable

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

def get_particles(ds, c):
    idat, rdat =  ds.get_particles(c["component"])
    return {"i":idat, "r":rdat}

def plot(frame, data, output_name):

    xn = data["nd-ion"]["x"]
    yn = data["nd-ion"]["y"]

    ni = data["nd-ion"]["value"]
    ni_min = frame["nd-ion"]["min"]
    ni_max = frame["nd-ion"]["max"]

    ne = data["nd-electron"]["value"]
    ne_min = frame["nd-electron"]["min"]
    ne_max = frame["nd-electron"]["max"]

    limits = frame["q"]["xy_limits"]

    # particles
    i_px, i_py, i_pi = get_particle_trajectories(data["particles-ion"], limits)
    e_px, e_py, e_pi = get_particle_trajectories(data["particles-electron"], limits)


    yn, xn = np.meshgrid(yn, xn)

    axes = []

    l = 0 #streak length

    fig = plt.figure(figsize=(6,3))
    ax = fig.add_subplot(211); axes.append(ax)

    pcm = ax.pcolormesh(xn, yn, ni, vmin=ni_min, vmax=ni_max)
    ax.plot(i_px[-l::], i_py[-l::], "w-", lw=0.5, alpha=0.5)
    ax.plot(i_px[-1], i_py[-1], "wo", ms=0.5, alpha=0.5)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='1.5%', pad=0.05)
    fig.colorbar(pcm, cax=cax, label=r"$n_i$")
    

    ax = fig.add_subplot(212); axes.append(ax)

    pcm = ax.pcolormesh(xn, yn, ne, vmin=ne_min, vmax=ne_max)
    ax.plot(e_px[-l::], e_py[-l::], "w-", lw=0.5, alpha=0.5)
    ax.plot(e_px[-1], e_py[-1], "wo", ms=0.5, alpha=0.5)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='1.5%', pad=0.05)
    fig.colorbar(pcm, cax=cax, label=r"$n_e$")

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

if 1:
    dt = 0.005

    Q = []

    q = {}
    q["files_dir"] = "."
    q["level"] = -1


    q["get"] = [
        {"func":get_number_density, "tag":"nd-ion", "component":"ion"},
        {"func":get_number_density, "tag":"nd-electron", "component":"electron"},
        {"func":get_particles, "tag":"particles-ion", "get_streak":True, "component":"ion"},
        {"func":get_particles, "tag":"particles-electron", "get_streak":True, "component":"electron"},
    ]

    q["plot"] = plot
    q["name"] = "movie"
    

    ##
    q["framerate"] = 20
    q["mov_save"] = q["files_dir"] + "/mov"
    q["offset"] = [0.0, 0.0]
    q["xy_limits"] = [[-1.25,0], [1.25,1]]
    q["file_include"] = ["IRMI.plt"]
    q["file_exclude"] = []
    q["cores"] = 8
    q["time_span"] = []  # np.arange(0,1+dt, dt)
    q["force_data"] = False
    q["force_frames"] = True
    q["only_frames"] = False
    q["redo_streaks"] = False
    q["dpi"] = 200

    q["normalize"] = "all"

    Q.append(q)

    make_all(Q)

print("DONE")
