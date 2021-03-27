
import sys
cmd_folder = "../../../vis"
if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)
    
from get_boxlib import ReadBoxLib, get_files

import numpy as np
import pylab as plt
import matplotlib.ticker as ticker

def check():

    #==============================================================================
    # Simulation results
    #==============================================================================

    # get a list of all the files in this directory
    files = get_files('.', include=['plt'], exclude=["temp"], get_all=True)

    f = files[-1]

    rh5 = ReadBoxLib(f)
    t = rh5.time
    c = rh5.data["lightspeed"]
        
    rh5 = ReadBoxLib(f, max_level=-1)

    name = "field"

    xn, By = rh5.get("y_B-%s"%name, grid='node')
    xc, Bz = rh5.get("z_B-%s"%name)

    x = xc[0]

    sim_data = [By, Bz]

    # =============================================================================
    # analytical solution
    # =============================================================================

    By0 = 1.0
    Bz0 = 0.0

    By1 = np.cos(1.5)
    Bz1 = np.sin(1.5)

    mBy = 0.5*(By0 + By1)
    mBz = 0.5*(Bz0 + Bz1)

    r_x = t*np.array([-10*c, -c, -c, c, c, 10*c])
    r_By =np.array([By0, By0, mBy, mBy, By1, By1])
    r_Bz =np.array([Bz0, Bz0, mBz, mBz, Bz1, Bz1])

    ref_data = [r_By, r_Bz]

    # =============================================================================
    # check
    # =============================================================================

    
    offset = 0.1
    sample_x = 2*[t*np.array([-(1+offset)*c, -(1-offset)*c, (1-offset)*c, (1+offset)*c]),]

    sample_y = []
    sample_err  = []

    success = 0

    for i in range(len(sim_data)):
        sim = sim_data[i]
        sy = np.interp(sample_x[i], x, sim)
        sample_y.append(sy)

        ref = ref_data[i]
        ry = np.interp(sample_x[i], r_x, ref)

        err = ry - sy

        Linf = np.max(np.abs(err))
        sample_err.append(err)

        if (Linf > 0.05):
            success = 1

    # =============================================================================
    # plot
    # =============================================================================

    axes = []
    fig = plt.figure(figsize=(8,8))

    nr = 2
    nc = 1

    count = 1

    ax = fig.add_subplot(nr,nc,count); axes.append(ax); count += 1
    ax.plot(x, By,'.',ms=1)
    ax.plot(r_x, r_By,'-',lw=0.5)
    ax.set_ylabel(r"$B_y$")

    ax = fig.add_subplot(nr,nc,count); axes.append(ax); count += 1
    ax.plot(x, Bz,'.',ms=1)
    ax.plot(r_x, r_Bz,'-',lw=0.5)
    ax.set_ylabel(r"$B_z$")


    for i, ax in enumerate(axes):
        ax.set_xlim(x[0], x[-1])
        ax.set_xlabel(r"$x$")

        ax.plot(sample_x[i], sample_y[i], 'ro', mfc='none')

    fig.tight_layout()

    fig.savefig("plot.png", dpi=300)
    plt.close(fig)

    return success
        

if __name__ == "__main__":
    sys.exit(check())