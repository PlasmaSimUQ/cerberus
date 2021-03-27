
import sys
cmd_folder = "../../../vis"
if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)
    
from get_boxlib import ReadBoxLib, get_files
from riemann_solver import riemannSolver

import numpy as np
import pylab as plt
import matplotlib.ticker as ticker

def check():

    #==============================================================================
    # Simulation results
    #==============================================================================

    # get a list of all the files in this directory
    files = get_files('.', include=['.plt'], exclude=["temp"], get_all=True)

    f = files[-1]

    rh5 = ReadBoxLib(f)
    t = rh5.time
        
    rh5 = ReadBoxLib(f, max_level=-1)

    name = "neutral"


    xc, rho = rh5.get("rho-%s"%name)
    xn, u = rh5.get("x_vel-%s"%name, grid='node')
    _, v = rh5.get("y_vel-%s"%name)
    _, w = rh5.get("z_vel-%s"%name)
    _, p = rh5.get("p-%s"%name)
    _, T = rh5.get("T-%s"%name)

    _, gam = rh5.get("gamma-%s"%name)

    g = gam[0]

    x = xc[0]

    N = len(x)

    x0 = xn[0]
    x1 = xn[-1]
    L = x1 - x0

    sim_data = [rho, u, p, T]

    # =============================================================================
    # analytical solution
    # =============================================================================

    rs = riemannSolver(x, 0.0, g, 1.0, rho[0], u[0], p[0], rho[-1], u[-1], p[-1],1,1,1,1)
    rdata = rs.compute(t)

    r_x = rdata[:,0]
    r_rho = rdata[:,1]
    r_u = rdata[:,2]
    r_p = rdata[:,3]
    r_T = rdata[:,4]

    ref_data = [r_rho, r_u, r_p, r_T]

    # =============================================================================
    # check
    # =============================================================================

    sample_x = [
        [-0.6, -0.5, -0.3, -0.2, 0.1, 0.24, 0.62, 0.68 ],
        [-0.6, -0.5, -0.3, -0.2, 0.62, 0.68],
        [-0.6, -0.5, -0.3, -0.2, 0.6, 0.68],
        [-0.6, -0.5, -0.3, 0.12, 0.23, 0.6, 0.68],
    ]

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
    nc = 2

    ax = fig.add_subplot(nr,nc,1); axes.append(ax)
    ax.plot(x, rho,'.',ms=1, label=name)
    ax.plot(r_x, r_rho,'-',lw=0.5)
    ax.set_ylabel(r"$\rho$")

    ax = fig.add_subplot(nr,nc,2); axes.append(ax)
    ax.plot(x, u,'.',ms=1)
    ax.plot(r_x, r_u,'-',lw=0.5)
    ax.set_ylabel(r"$u$")

    ax = fig.add_subplot(nr,nc,3); axes.append(ax)
    ax.plot(x, p,'.',ms=1)
    ax.plot(r_x, r_p,'-',lw=0.5)
    ax.set_ylabel(r"$p$")

    ax = fig.add_subplot(nr,nc,4); axes.append(ax)
    ax.plot(x, T,'.',ms=1)
    ax.plot(r_x, r_T,'-',lw=0.5)
    ax.set_ylabel(r"$T$")


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